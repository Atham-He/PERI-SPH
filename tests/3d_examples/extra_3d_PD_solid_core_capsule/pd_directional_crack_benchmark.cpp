#include "pd_snsni_capsule_fixture.h"
#include <Eigen/Geometry>
#include <filesystem>
#include <iomanip>
using namespace SPH;using namespace SPH::solid_dynamics;
namespace SPH{
class DirectionalCrackSlab;
template<>class ParticleGenerator<BaseParticles,DirectionalCrackSlab>:public ParticleGenerator<BaseParticles>{
    Real dp,t;int nx,ny;
  public:ParticleGenerator(SPHBody&b,BaseParticles&p,Real d,Real h,int n,int m):ParticleGenerator<BaseParticles>(b,p),dp(d),t(h),nx(n),ny(m){}
    void prepareGeometricData()override{ThicknessQuadrature rule(3);for(size_t k=0;k<3;++k)for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)addPositionAndVolumetricMeasure(Vecd((i+.5-nx/2.)*dp,(j+.5-ny/2.)*dp,.5*t*rule.xi[k]),dp*dp*.5*t*rule.weight[k]);}
};}
int main(int argc,char**argv){try{
    if(argc!=2)throw std::invalid_argument("Output directory required");std::filesystem::path out(argv[1]);std::filesystem::create_directories(out);
    const Real E=2.151e9,nu=.3,ft=46.4e6,Gc=100,dp=25e-6,t=5e-6;const int nx=16,ny=10;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-.002),Vecd::Constant(.002)),dp,8);PDBody body(system,makeShared<DefaultShape>("DirectionalFrozenCrackLayer"));
    body.defineAdaptation<PDAdaptation>(1.5075,1);body.defineMatterMaterial<HughesWingetSolid>(1200,E,nu);body.generateParticles<BaseParticles,DirectionalCrackSlab>(dp,t,nx,ny);
    InnerRelation inner(body);SimpleDynamics<PDTimeStepInitialization> init(body,makeShared<Gravity>(Vecd::Zero()));SimpleDynamics<NosbPDFirstStep> first(body);InteractionWithUpdate<NosbPDShapeMatrix> shape(inner);
    auto&p=body.getBaseParticles();size_t count=p.TotalRealParticles(),columns_count=count/3;StdVec<StdVec<size_t>>columns(columns_count);StdVec<Vecd>directors(columns_count,Vecd::UnitZ());for(size_t c=0;c<columns_count;++c)for(size_t k=0;k<3;++k)columns[c].push_back(k*columns_count+c);
    ThicknessQuadrature rule(3);peridynamics::LayeredParticleGeometry geometry(p,columns,directors,rule.xi);PDLayeredQuadraticReference ref(inner,geometry);
    InteractionWithUpdate<NosbPDQuadraticSecondStep> second(inner,ref);InteractionDynamics<NosbPDQuadraticThirdStep> third(inner,ref);InteractionDynamics<ProjectedPairStabilization> old_hg(inner,.5,.1);InteractionDynamics<PairNumericalDampingforPD> old_damping(inner);
    PDEnergyCrackBand model(inner,ref,E,nu,ft,Gc,.5,body.getSPHAdaptation().getKernel()->CutOffRadius(),true);
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();ref.build();model.initialize();State s(body);
    auto*acc=p.getVariableDataByName<Vecd>("PDAcceleration");auto*force=p.getVariableDataByName<Vecd>("PDInternalForce");auto*hgforce=p.getVariableDataByName<Vecd>("HourglassForce");auto*stress=p.getVariableDataByName<Matd>("PDStress");
    auto evaluate=[&](bool advance=false){for(size_t i=0;i<count;++i)model.evaluate(i,advance);for(size_t i=0;i<count;++i){acc[i].setZero();model.forces(i);}return model.statistics();};
    auto set=[&](const Matd&F,Real d){for(size_t i=0;i<count;++i){s.x[i]=F*s.x0[i];model.history[i]={0,d,0};model.crack_normal[i]=Vecd::UnitX();model.crack_projection_ready[i]=0;}};
    std::ofstream csv(out/"opening_closure_sliding.csv");csv<<std::setprecision(17)<<"mode,normal_height,sliding,normal_traction_Pa,shear_traction_Pa,stored_J,reference_stored_J,force_balance,torque_balance,force_energy_error,timestep_bound_s\n";
    Real volume=0;for(size_t i=0;i<count;++i)volume+=s.measure[i];Real max_open=0,max_sliding=0,max_close=0,max_fd=0,max_balance=0,max_torque=0,max_objective=0;
    struct Mode{const char*name;Real h,slide;};
    // A single frozen regularized crack-layer specimen; successive boundary
    // modes include finite opening, closure, frictionless slip and reopening.
    StdVec<Mode>modes={{"open",1.02,0},{"open",2,0},{"open",10,0},{"open_slip",2,8},{"close",.98,0},{"closed_slip",.98,2},{"reopen",1.1,2}};
    for(const auto&mode:modes){Matd F=Matd::Identity();F(0,0)=mode.h;F(1,0)=mode.slide;set(F,1);auto stat=evaluate();
        Real tn=0,ts=0;Vecd sum=Vecd::Zero(),torque=Vecd::Zero();Real scale=0,power=0;Matd V;V<<.13,.02,-.01,.04,-.07,.03,.01,.02,.09;
        for(size_t i=0;i<count;++i){tn+=stress[i](0,0)/count;ts+=stress[i](1,0)/count;Vecd f=force[i]+hgforce[i];sum+=f;torque+=s.x[i].cross(f);scale+=f.norm();s.v[i]=V*s.x0[i];power+=f.dot(s.v[i]);}
        Real balance=sum.norm()/std::max(scale,E*volume/(nx*dp)),moment=torque.norm()/std::max(scale*nx*dp,E*volume);
        Real eps=2e-7;for(size_t i=0;i<count;++i)s.x[i]+=eps*s.v[i];Real ep=evaluate().stored;for(size_t i=0;i<count;++i)s.x[i]-=2*eps*s.v[i];Real em=evaluate().stored;
        Real fd=std::abs((ep-em)/(2*eps)+power)/(E*volume);max_fd=std::max(max_fd,fd);max_balance=std::max(max_balance,balance);max_torque=std::max(max_torque,moment);
        Matd no_slip=Matd::Identity();no_slip(0,0)=mode.h;auto intact=model.elastic.evaluate(no_slip);Real expected=mode.h<1?volume*(intact.positive+intact.negative):0;
        if(mode.h>=1)max_open=std::max(max_open,std::abs(stat.stored)/(E*volume));else max_close=std::max(max_close,std::abs(stat.stored-expected)/std::max(expected,Real(1e-30)));
        max_sliding=std::max(max_sliding,std::abs(ts)/E);
        csv<<mode.name<<','<<mode.h<<','<<mode.slide<<','<<tn<<','<<ts<<','<<stat.stored<<','<<expected<<','<<balance<<','<<moment<<','<<fd<<','<<model.finiteStrainTimeStep(ref.elastic_dt_bound)<<'\n';
        Matd Q=Eigen::AngleAxis<Real>(.73,Vecd(1,2,3).normalized()).toRotationMatrix();auto a=model.crack_energy.evaluate(F,Vecd::UnitX()),b=model.crack_energy.evaluate(Q*F,Vecd::UnitX());max_objective=std::max(max_objective,(b.P-Q*a.P).norm()/E);
    }
    // The condensed Heaviside residual must not leave an artificial spring
    // across its LOCAL representable jump mode; other residual modes stay stiff.
    // A globally sharp jump is not exactly represented by every point-centered
    // plane in this diffuse crack-band discretization; report that residual.
    Real local_mode_error=0;for(size_t i=0;i<count;++i){model.buildCrackProjection(i,Vecd::UnitX());const auto &row=model.projection[i];Real norm=0;PDEnergyCrackBand::V9 orth=PDEnergyCrackBand::V9::Zero();for(const auto&e:row){orth+=e.weight*e.basis*e.crack_mode;norm+=e.weight*e.crack_mode*e.crack_mode;}local_mode_error=std::max(local_mode_error,orth.norm()/std::max(std::sqrt(norm),Real(1)));}
    for(size_t i=0;i<count;++i){model.history[i]={0,1,0};model.crack_normal[i]=Vecd::UnitX();model.crack_projection_ready[i]=0;s.x[i]=s.x0[i]+Vecd(s.x0[i][0]>0?1e-7:-1e-7,0,0);}auto jump=evaluate();Real jump_hg=jump.stabilization/(E*volume);
    for(size_t i=0;i<count;++i)s.x[i]=s.x0[i]+Vecd(0,0,1e-8*std::sin(13*s.x0[i][0]/dp)*std::cos(11*s.x0[i][1]/dp));auto nonpolynomial=evaluate();Real surviving_hg=nonpolynomial.stabilization;
    // Partial damage: verify the exact transpose material + residual force,
    // including energy, linear momentum and angular momentum in one state.
    Matd F=Matd::Identity();F(0,0)=1.04;F(1,1)=.993;set(F,.4);
    for(size_t i=0;i<count;++i){s.x[i]+=Vecd(0,0,2e-9*std::sin(9*s.x0[i][0]/dp));s.v[i]=Vecd(.02*std::sin(3*s.x0[i][0]/dp),.03*std::cos(4*s.x0[i][1]/dp),.01);}
    evaluate();Real power=0,scale=0;Vecd sum=Vecd::Zero(),torque=Vecd::Zero();for(size_t i=0;i<count;++i){Vecd f=force[i]+hgforce[i];power+=f.dot(s.v[i]);scale+=f.norm();sum+=f;torque+=s.x[i].cross(f);}
    Real eps=1e-12;for(size_t i=0;i<count;++i)s.x[i]+=eps*s.v[i];Real ep=evaluate().stored;for(size_t i=0;i<count;++i)s.x[i]-=2*eps*s.v[i];Real em=evaluate().stored;
    Real partial_fd=std::abs((ep-em)/(2*eps)+power)/std::max(std::abs(power),Real(1e-20));
    max_balance=std::max(max_balance,sum.norm()/scale);max_torque=std::max(max_torque,torque.norm()/(scale*nx*dp));
    // Undamped explicit velocity-Verlet hold of the SAME closed/sliding crack
    // layer. Fixed perimeter does no work; free interior masses must remain
    // stable and conserve total energy with the computed physical time bound.
    F=Matd::Identity();F(0,0)=.98;F(1,0)=.25;set(F,1);
    StdVec<int> fixed(count);for(size_t i=0;i<count;++i){size_t c=i%columns_count,ix=c/ny,iy=c%ny;fixed[i]=ix==0||ix==nx-1||iy==0||iy==ny-1;s.v[i].setZero();if(!fixed[i])s.x[i][2]+=1e-10*std::sin(7*s.x0[i][0]/dp)*std::sin(5*s.x0[i][1]/dp);}
    Real e0=evaluate().stored,dt=.2*std::min(model.finiteStrainTimeStep(ref.elastic_dt_bound),model.stabilizer_dt),energy_drift=0,min_dt_bound=std::numeric_limits<Real>::infinity();
    for(int step=0;step<200;++step){for(size_t i=0;i<count;++i)if(!fixed[i]){s.v[i]+=.5*dt*acc[i];s.x[i]+=dt*s.v[i];}auto stat=evaluate();Real kinetic=0;for(size_t i=0;i<count;++i)if(!fixed[i]){s.v[i]+=.5*dt*acc[i];kinetic+=.5*s.mass[i]*s.v[i].squaredNorm();}energy_drift=std::max(energy_drift,std::abs((stat.stored+kinetic)/e0-1));min_dt_bound=std::min(min_dt_bound,model.finiteStrainTimeStep(ref.elastic_dt_bound));}
    // Integrate a monotone uniaxial separation with the unchanged prescribed ft
    // and Gc. The remaining crack-parallel elastic energy is subtracted from
    // external work, rather than falsely counting it as fracture dissipation.
    auto sample=[&](Real axial){Matd f=model.elastic.uniaxial(axial);auto tr=model.elastic.evaluate(f);Real retained=model.crack_energy.energy(f,Vecd::UnitX());LinearEnergyCrackBand::History h;model.fracture.advance(h,std::max(Real(0),tr.positive+tr.negative-retained),true);auto cr=model.crack_energy.evaluate(f,Vecd::UnitX());Matd pk=(1-h.damage)*(tr.positive_P+tr.negative_P)+h.damage*cr.P;return std::make_tuple(f,pk,(1-h.damage)*(tr.positive+tr.negative)+h.damage*retained,h);};
    Real end=1.04;while(std::get<3>(sample(end)).damage<1){end=1+2*(end-1);if(end>8)throw std::runtime_error("Gc separation did not complete");}
    Real work=0;auto previous=sample(1);const int steps=12000;Real peak=0;
    for(int k=1;k<=steps;++k){auto next=sample(1+(end-1)*k/steps);work+=.5*((std::get<1>(previous)+std::get<1>(next)).array()*(std::get<0>(next)-std::get<0>(previous)).array()).sum();peak=std::max(peak,std::get<1>(next)(0,0));previous=next;}
    Real recovered=(work-std::get<2>(previous))*model.fracture.width,work_error=std::abs(recovered/Gc-1),diss_error=std::abs(std::get<3>(previous).dissipation*model.fracture.width/Gc-1);
    bool passed=max_open<1e-12&&max_sliding<1e-10&&max_close<1e-9&&max_fd<1e-7&&partial_fd<1e-5&&max_balance<1e-10&&max_torque<1e-9&&max_objective<1e-10&&local_mode_error<1e-12&&surviving_hg>0&&work_error<2e-4&&diss_error<1e-12&&energy_drift<1e-3&&dt<=min_dt_bound;
    std::ofstream report(out/"acceptance.json");report<<std::setprecision(17)<<"{\"status\":\""<<(passed?"PASSED":"FAILED")<<"\",\"specimens\":1,\"particles\":"<<count<<",\"E_Pa\":"<<E<<",\"nu\":"<<nu<<",\"ft_Pa\":"<<ft<<",\"Gc_J_m2\":"<<Gc<<",\"opening_energy_over_E_volume\":"<<max_open<<",\"sliding_traction_over_E\":"<<max_sliding<<",\"closure_energy_relative_error\":"<<max_close<<",\"force_energy_error_over_E_volume\":"<<max_fd<<",\"partially_damaged_force_energy_relative_error\":"<<partial_fd<<",\"force_balance\":"<<max_balance<<",\"torque_balance\":"<<max_torque<<",\"objectivity_stress_error_over_E\":"<<max_objective<<",\"jump_stabilizer_energy_over_E_volume\":"<<jump_hg<<",\"local_jump_mode_projection_error\":"<<local_mode_error<<",\"noncrack_mode_stabilizer_J\":"<<surviving_hg<<",\"integrated_Gc_J_m2\":"<<recovered<<",\"Gc_external_work_relative_error\":"<<work_error<<",\"Gc_history_relative_error\":"<<diss_error<<",\"observed_uniaxial_peak_Pa\":"<<peak<<",\"explicit_undamped_steps\":200,\"explicit_dt_s\":"<<dt<<",\"explicit_energy_relative_drift\":"<<energy_drift<<",\"minimum_explicit_bound_s\":"<<min_dt_bound<<",\"diffuse_directional_contact_not_sharp_fragment_geometry\":true}\n";
    std::cout<<(passed?"DIRECTIONAL_CRACK_BENCHMARK_PASSED":"DIRECTIONAL_CRACK_BENCHMARK_FAILED")<<std::endl;return passed?0:2;
 }catch(const std::exception&e){std::cerr<<"DIRECTIONAL_BENCHMARK_ERROR "<<e.what()<<std::endl;return 1;}}

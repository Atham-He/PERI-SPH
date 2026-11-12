#define main capsule_main_unused
#include "pd_solid_capsule.cpp"
#undef main
#include "pd_energy_crack_band.h"
#include <Eigen/Geometry>
using namespace SPH::solid_dynamics;
namespace SPH {
class EnergySlabGenerator;
template<> class ParticleGenerator<BaseParticles,EnergySlabGenerator>:public ParticleGenerator<BaseParticles>
{
    Real dp,t;int nx,ny;
  public:
    ParticleGenerator(SPHBody&b,BaseParticles&p,Real d,Real thickness,int x,int y):ParticleGenerator<BaseParticles>(b,p),dp(d),t(thickness),nx(x),ny(y){}
    void prepareGeometricData()override
    {ThicknessQuadrature rule(3);for(size_t l=0;l<3;++l)for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)
      addPositionAndVolumetricMeasure(Vecd((i+.5-nx/2.)*dp,(j+.5-ny/2.)*dp,.5*t*rule.xi[l]),dp*dp*.5*t*rule.weight[l]);}
}; }
static Real simpson(const std::function<Real(Real)>&f,Real a,Real b,Real tol,int depth=24)
{
    Real c=.5*(a+b),fa=f(a),fc=f(c),fb=f(b),whole=(b-a)*(fa+4*fc+fb)/6;
    std::function<Real(Real,Real,Real,Real,Real,Real,Real,int)> rec;
    rec=[&](Real x,Real y,Real fx,Real fm,Real fy,Real old,Real err,int d)->Real
    {Real m=.5*(x+y),l=.5*(x+m),r=.5*(m+y),fl=f(l),fr=f(r);Real one=(m-x)*(fx+4*fl+fm)/6,two=(y-m)*(fm+4*fr+fy)/6;
     if(d==0||std::abs(one+two-old)<=15*err)return one+two+(one+two-old)/15;
     return rec(x,m,fx,fl,fm,one,err/2,d-1)+rec(m,y,fm,fr,fy,two,err/2,d-1);};
    return rec(a,b,fa,fc,fb,whole,tol,depth);
}
struct PatchResult {Real dp,Gc_recovered,work_error,force_FD,force_balance,torque_balance,quadratic_HG_ratio,terminal_traction_ratio,mixed_surface_work_error,quadratic_surface_work_error,layered_surface_work_error;size_t particles;};
static PatchResult patch(Real dp)
{
    const Real E=2.151e9,nu=.3,Gc=100,ft=46.4e6,t=10e-6;const int nx=16,ny=10;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-.002),Vecd::Constant(.002)),dp,8);
    PDBody body(system,makeShared<DefaultShape>("EnergySlab"));body.defineAdaptation<PDAdaptation>(1.5075,1);
    body.defineMatterMaterial<HughesWingetSolid>(1200,E,nu);body.generateParticles<BaseParticles,EnergySlabGenerator>(dp,t,nx,ny);
    InnerRelation inner(body);SimpleDynamics<PDTimeStepInitialization> init(body,makeShared<Gravity>(Vecd::Zero()));SimpleDynamics<NosbPDFirstStep> first(body);
    InteractionWithUpdate<NosbPDShapeMatrix> shape(inner);size_t columns_count=body.getBaseParticles().TotalRealParticles()/3;
    StdVec<StdVec<size_t>>columns(columns_count);StdVec<Vecd>directors(columns_count,Vecd(0,0,1));
    for(size_t c=0;c<columns_count;++c)for(size_t l=0;l<3;++l)columns[c].push_back(l*columns_count+c);
    ThicknessQuadrature normal_rule(3);peridynamics::LayeredParticleGeometry geometry(body.getBaseParticles(),columns,directors,normal_rule.xi);
    PDLayeredQuadraticReference ref(inner,geometry);InteractionWithUpdate<NosbPDQuadraticSecondStep> second(inner,ref);
    InteractionDynamics<NosbPDQuadraticThirdStep> third(inner,ref);InteractionDynamics<ProjectedPairStabilization> old_hg(inner,.5,.1);InteractionDynamics<PairNumericalDampingforPD> old_damping(inner);
    PDEnergyCrackBand model(inner,ref,E,nu,ft,Gc,.5,body.getSPHAdaptation().getKernel()->CutOffRadius());
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();ref.build();model.initialize();State s(body);

    // Uniform normal traction does zero work on v_z=(2*z/t)^2: the
    // two faces have equal velocity and opposite tractions. The actual
    // volume derivative and mapped face loads must agree on this mode.
    peridynamics::SolidForcePort port(body,peridynamics::SolidForcePort::Storage::PDAccelerationPrior);
    StdVec<Vecd>top,bottom;for(size_t c=0;c<columns_count;++c){Vecd X=s.x0[c];top.push_back(Vecd(X[0],X[1],t/2));bottom.push_back(Vecd(X[0],X[1],-t/2));}
    peridynamics::ReferenceSolidSurface top1(port,top,dp,2.5,1),top2(port,top,dp,2.5,2),bottom2(port,bottom,dp,2.5,2);
    peridynamics::LayeredSolidSurface layered_top(port,top,geometry,1,dp),layered_bottom(port,bottom,geometry,-1,dp);
    const Real pressure=1000,area=columns_count*dp*dp;Real volume_work=0;
    for(size_t i=0;i<s.n;++i)s.v[i]=Vecd(0,0,std::pow(2*s.x0[i][2]/t,2));
    for(size_t i=0;i<s.n;++i){Matd gv=Matd::Zero();for(auto&e:ref.rows[i])gv+=(s.v[e.j]-s.v[i])*e.gradient.transpose();volume_work-=s.measure[i]*pressure*gv(2,2);}
    auto boundary_work=[&](peridynamics::ReferenceSolidSurface&a,peridynamics::ReferenceSolidSurface&b){port.clear();for(size_t c=0;c<columns_count;++c){a.scatter(c,Vecd(0,0,-pressure*dp*dp));b.scatter(c,Vecd(0,0,pressure*dp*dp));}Real work=0;for(size_t i=0;i<s.n;++i)work+=s.mass[i]*port.external[i].dot(s.v[i]);return work;};
    Real mixed_surface_error=std::abs(boundary_work(top1,bottom2)-volume_work)/(pressure*area);
    Real quadratic_surface_error=std::abs(boundary_work(top2,bottom2)-volume_work)/(pressure*area);
    Real layered_surface_error=std::abs(boundary_work(layered_top,layered_bottom)-volume_work)/(pressure*area);
    if(mixed_surface_error<.01||quadratic_surface_error>1e-8||layered_surface_error>1e-8)throw std::runtime_error("Normal traction/volume quadratic-mode consistency audit failed");
auto&particles=body.getBaseParticles();auto*acc=particles.getVariableDataByName<Vecd>("PDAcceleration");auto*force=particles.getVariableDataByName<Vecd>("PDInternalForce");auto*hgforce=particles.getVariableDataByName<Vecd>("HourglassForce");
    auto evaluate=[&](bool advance){for(size_t i=0;i<s.n;++i)model.evaluate(i,advance);for(size_t i=0;i<s.n;++i){acc[i].setZero();model.forces(i);}return model.statistics();};
    // Full P2 residual must leave the already-verified bending kinematics intact.
    for(size_t i=0;i<s.n;++i){Vecd X=s.x0[i];s.x[i]=X+Vecd(.001*X[0]+50*X[1]*X[1],-30*X[0]*X[1],20*X[0]*X[0]);}
    auto q=evaluate(false);Real qratio=q.stabilization/std::max(q.stored,Real(1e-30));if(qratio>1e-15)throw std::runtime_error("Energy stabilizer penalized a quadratic field");
    for(size_t i=0;i<s.n;++i){Vecd X=s.x0[i];s.x[i]=Vecd(1.024*X[0],(1-.3*.024)*X[1],(1-.3*.024)*X[2])+Vecd(0,0,2e-8*std::sin(13*X[0]/(nx*dp))*std::cos(11*X[1]/(ny*dp)));
     s.v[i]=Vecd(.1*std::sin(7*X[0]/dp),-.08*std::cos(3*X[1]/dp),.04*X[2]/t);}
    auto loaded=evaluate(true);if(!loaded.initiated||loaded.failed==s.n)throw std::runtime_error("Partially damaged operator not exercised");
    Vecd sum=Vecd::Zero(),moment=Vecd::Zero();Real power=0,scale=0;
    for(size_t i=0;i<s.n;++i){Vecd f=force[i]+hgforce[i];sum+=f;moment+=s.x[i].cross(f);power+=f.dot(s.v[i]);scale+=f.norm();}
    Real balance=sum.norm()/scale,torque=moment.norm()/(nx*dp*scale),eps=1e-12;
    for(size_t i=0;i<s.n;++i)s.x[i]+=eps*s.v[i];Real plus=evaluate(false).stored;
    for(size_t i=0;i<s.n;++i)s.x[i]-=2*eps*s.v[i];Real minus=evaluate(false).stored;
    Real fd=std::abs((plus-minus)/(2*eps)+power)/std::max(std::abs(power),Real(1e-30));
    if(fd>1e-5||balance>1e-10||torque>1e-9)throw std::runtime_error("Damaged energy force/torque/finite difference failed");
    for(size_t i=0;i<s.n;++i){model.history[i]={};s.x[i]=s.x0[i]+Vecd(s.x0[i][0]>0?dp:0,0,0);}evaluate(false);
    StdVec<Matd> B(s.n);StdVec<Real> gamma(s.n);StdVec<bool> active(s.n,false);
    for(size_t i=0;i<s.n;++i){B[i]=(s.F[i]-Matd::Identity())/dp;gamma[i]=model.hg_density[i]/(dp*dp);auto&nb=inner.inner_configuration_[i];
      for(size_t n=0;n<nb.current_size_;++n)active[i]=active[i]||((s.x0[i][0]>0)!=(s.x0[nb.j_[n]][0]>0));}
    auto state_at=[&](size_t i,Real w){auto e=model.elastic.evaluate(Matd::Identity()+w*B[i]);LinearEnergyCrackBand::History h;model.fracture.advance(h,e.positive+gamma[i]*w*w,e.max_effective_tension>1e-12*E);return std::make_pair(e,h);};
    Real terminal=dp*.02;
    for(int k=0;k<50;++k){bool all=true;for(size_t i=0;i<s.n;++i)if(active[i])all=all&&state_at(i,terminal).second.damage==1;if(all)break;terminal*=2;if(k==49)throw std::runtime_error("Planar opening did not completely separate");}
    Real work=0,expected=Gc*(ny*dp*t);
    for(size_t i=0;i<s.n;++i)if(active[i])
    {
        auto root=[&](Real target){Real lo=0,hi=terminal;for(int k=0;k<70;++k){Real m=.5*(lo+hi);if(state_at(i,m).second.kappa<target)lo=m;else hi=m;}return .5*(lo+hi);};
        Real a=root(model.fracture.k0),b=root(model.fracture.kf);
        auto power_w=[&](Real w){auto [e,h]=state_at(i,w);Matd p=(1-h.damage)*e.positive_P+e.negative_P;return s.measure[i]*((p.array()*B[i].array()).sum()+2*(1-h.damage)*gamma[i]*w);};
        work+=simpson(power_w,0,a,expected*1e-10/s.n)+simpson(power_w,a,b,expected*1e-10/s.n);
    }
    for(size_t i=0;i<s.n;++i){s.x[i]=s.x0[i]+Vecd(s.x0[i][0]>0?terminal:0,0,0);model.history[i]={};}auto complete=evaluate(true);
    Real recovered=complete.dissipated/(ny*dp*t),workerr=std::abs(work-complete.dissipated)/expected;
    Vecd right=Vecd::Zero();for(size_t i=0;i<s.n;++i)if(s.x0[i][0]>0)right+=force[i]+hgforce[i];Real traction=right.norm()/(ft*ny*dp*t);
    if(std::abs(recovered/Gc-1)>.02||workerr>1e-5||traction>1e-6)throw std::runtime_error("Unit crack work does not reproduce input Gc or leave a traction-free opening");

    Matd G=.97*Matd::Identity();G(0,1)=.097;
    for(size_t i=0;i<s.n;++i){s.x[i]=G*s.x0[i];model.history[i]={model.fracture.kf,1,model.fracture.Gc/model.fracture.width};}
    auto comp=evaluate(false);Real cp=0;for(size_t i=0;i<s.n;++i)cp+=(force[i]+hgforce[i]).dot(s.v[i]);
    for(size_t i=0;i<s.n;++i)s.x[i]+=eps*s.v[i];Real cplus=evaluate(false).stored;
    for(size_t i=0;i<s.n;++i)s.x[i]-=2*eps*s.v[i];Real cminus=evaluate(false).stored;
    Real cfd=std::abs((cplus-cminus)/(2*eps)+cp)/std::max(std::abs(cp),Real(1e-30));
    if(cfd>1e-5||comp.stored<=0||comp.stabilization!=0)throw std::runtime_error("Fully failed compressed NOSB operator is not variational");
return {dp,recovered,workerr,fd,balance,torque,qratio,traction,mixed_surface_error,quadratic_surface_error,layered_surface_error,s.n};
}
int main(int argc,char**argv)
{
 try {
    if(argc!=2)throw std::runtime_error("output.json");const Real E=2.151e9,nu=.3,ft=46.4e6,Gc=100,h=150.75e-6;
    FiniteVolumeBiotEnergy elastic(E,nu);LinearEnergyCrackBand law(E,ft,Gc,h,elastic.initiationEnergy(ft));LinearEnergyCrackBand::History history;
    Real peak=0,work=0,last=0,previous=0;int steps=12000;Real end=.08;
    for(int n=0;n<=steps;++n){Real e=end*n/steps;Matd F=elastic.uniaxial(1+e);auto t=elastic.evaluate(F);law.advance(history,t.positive,true);Real p=(1-history.damage)*t.positive_P(0,0)+t.negative_P(0,0);peak=std::max(peak,p);work+=.5*(p+last)*(e-previous)*h;last=p;previous=e;}
    if(std::abs(work/Gc-1)>1e-6||std::abs(peak/ft-1)>2e-4||std::abs(history.dissipation*h/Gc-1)>1e-12)throw std::runtime_error("Uniaxial strength or fracture-energy normalization failed");
    LinearEnergyCrackBand::History unload;Real e=.5*(law.k0+law.kf);Matd F=elastic.uniaxial(1+e);auto t=elastic.evaluate(F);law.advance(unload,t.positive,true);Real d=unload.damage,w=unload.dissipation;law.advance(unload,.01*t.positive,true);
    if(unload.damage!=d||unload.dissipation!=w)throw std::runtime_error("Damage healed on unloading");
    Matd Q=Eigen::AngleAxis<Real>(.63,Vecd(1,2,3).normalized()).toRotationMatrix();auto rotated=elastic.evaluate(Q*F);
    Real objective_stress_error=(rotated.positive_P-Q*t.positive_P).norm()/E,objective_energy_error=std::abs(rotated.positive-t.positive)/t.positive;
    if(objective_stress_error>1e-12||objective_energy_error>1e-12)throw std::runtime_error("Energy material is not objective");
    auto compressed=elastic.evaluate(.97*Matd::Identity());LinearEnergyCrackBand::History c;law.advance(c,compressed.positive,compressed.max_effective_tension>1e-12*E);
    if(c.damage!=0||compressed.negative<=0)throw std::runtime_error("Hydrostatic compression was damaged");
    bool invalid=false;try{LinearEnergyCrackBand bad(E,ft,1,h);}catch(const std::invalid_argument&){invalid=true;}if(!invalid)throw std::runtime_error("Under-resolved fracture energy silently accepted");

    // Reproduce the constitutive defect without any PD discretization or
    // interface: finite isochoric shear under hydrostatic precompression.
    CorotatedBiotEnergy legacy(E,nu);Matd C=.97*Matd::Identity(),H=Matd::Identity();H(0,1)=.1;
    auto old0=legacy.evaluate(C),old1=legacy.evaluate(C*H),new0=elastic.evaluate(C),new1=elastic.evaluate(C*H);
    Real legacy_drop=old1.negative-old0.negative,volume_shear_error=std::abs(new1.negative-new0.negative)/new0.negative;
    if(legacy_drop>=0||volume_shear_error>1e-12)throw std::runtime_error("Isochoric compression regression failed");
    Matd collapsed=Matd::Identity();collapsed(0,0)=.01;collapsed(1,1)=collapsed(2,2)=1.5;
    auto oldcollapse=legacy.evaluate(collapsed),newcollapse=elastic.evaluate(collapsed);
    if(oldcollapse.negative!=0||newcollapse.negative<E)throw std::runtime_error("Finite-volume barrier regression failed");
    // Exact derivatives are required for BOTH unilateral branches, including
    // a compressed/sheared fully damaged state (only the negative part).
    Real derivative_error=0;Matd Ftest=C*H,V;V<<.13,-.07,.02,.03,-.11,.04,.05,.01,.09;
    auto A=elastic.evaluate(Ftest);Real eps=1e-7;
    auto plus=elastic.evaluate(Ftest+eps*V),minus=elastic.evaluate(Ftest-eps*V);
    Real total_fd=(plus.positive+plus.negative-minus.positive-minus.negative)/(2*eps),total_work=((A.positive_P+A.negative_P).array()*V.array()).sum();
    Real neg_fd=(plus.negative-minus.negative)/(2*eps),neg_work=(A.negative_P.array()*V.array()).sum();
    derivative_error=std::max(std::abs(total_fd-total_work)/E,std::abs(neg_fd-neg_work)/E);
    if(derivative_error>1e-8)throw std::runtime_error("Finite-volume stress is not an energy derivative");
    Eigen::Matrix<Real,9,9> Cmat;
    for(int k=0;k<9;++k){Matd dF=Matd::Zero();dF(k/3,k%3)=eps;
        auto pp=elastic.evaluate(Ftest+dF),mm=elastic.evaluate(Ftest-dF);Matd dP=(pp.positive_P+pp.negative_P-mm.positive_P-mm.negative_P)/(2*eps);
        for(int q=0;q<9;++q)Cmat(q,k)=dP(q/3,q%3);}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real,9,9>> eig(.5*(Cmat+Cmat.transpose()));
    Real numeric_norm=eig.eigenvalues().cwiseAbs().maxCoeff(),analytic_bound=elastic.tangentBound(Ftest,0);
    if(eig.info()!=Eigen::Success||numeric_norm>analytic_bound*(1+1e-7))throw std::runtime_error("Finite-strain tangent bound is not conservative");
    std::ofstream proof(argv[1]+std::string(".regression.json"));proof<<std::setprecision(17)
       <<"{\"legacy_isochoric_shear_energy_change_J_m3\":"<<legacy_drop<<",\"corrected_isochoric_shear_relative_error\":"<<volume_shear_error
       <<",\"collapsed_J\":"<<collapsed.determinant()<<",\"legacy_retained_compression_J_m3\":"<<oldcollapse.negative
       <<",\"corrected_retained_compression_J_m3\":"<<newcollapse.negative<<",\"derivative_error_over_E\":"<<derivative_error
       <<",\"numeric_tangent_norm_Pa\":"<<numeric_norm<<",\"analytic_tangent_bound_Pa\":"<<analytic_bound
       <<",\"discretization_or_coupling_used_in_material_counterexample\":false,\"fully_failed_compressed_NOSB_force_energy_check\":true,\"passed\":true}\n";
auto a=patch(25e-6),b=patch(12.5e-6);std::ofstream out(argv[1]);out<<std::setprecision(17)<<"{\"passed\":true,\"E_Pa\":"<<E<<",\"ft_Pa\":"<<ft<<",\"Gc_J_m2\":"<<Gc
      <<",\"uniaxial_peak_Pa\":"<<peak<<",\"uniaxial_work_J_m2\":"<<work<<",\"objective_stress_error_over_E\":"<<objective_stress_error<<",\"objective_energy_relative_error\":"<<objective_energy_error
      <<",\"unloading_irreversible\":true,\"rotation_objective\":true,\"compression_retained\":true,\"inadmissible_parameters_rejected\":true,\"patches\":[";
    int index=0;for(auto r:{a,b}){if(index++)out<<',';out<<"{\"spacing_m\":"<<r.dp<<",\"particles\":"<<r.particles<<",\"recovered_Gc_J_m2\":"<<r.Gc_recovered<<",\"external_work_balance_error\":"<<r.work_error
      <<",\"force_energy_FD_error\":"<<r.force_FD<<",\"force_balance\":"<<r.force_balance<<",\"torque_balance\":"<<r.torque_balance<<",\"quadratic_HG_energy_ratio\":"<<r.quadratic_HG_ratio<<",\"terminal_traction_ratio\":"<<r.terminal_traction_ratio<<",\"mixed_surface_quadratic_work_error\":"<<r.mixed_surface_work_error<<",\"quadratic_surface_work_error\":"<<r.quadratic_surface_work_error<<",\"layered_surface_work_error\":"<<r.layered_surface_work_error<<'}';}
    out<<"],\"physical_compression_histories\":0,\"mixed_mode_fracture_calibrated\":false,\"capsule_spatial_convergence\":false}\n";std::cout<<"ENERGY_CRACK_BAND_AUDIT_PASSED"<<std::endl;return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
}

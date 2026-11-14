#include "capsule_core_fixture.h"
#include "pd_solid_reference_coupling.h"
#include "solid_snsni_plastic.h"
#include <tbb/task_arena.h>
#include "recovery_checkpoint.h"
using namespace SPH::solid_dynamics;
using Field=StdVec<Vecd>;
Real dot(const Field&a,const Field&b){Real x=0;for(size_t j=0;j<a.size();++j)x+=a[j].dot(b[j]);return x;}
class RKSurface:public peridynamics::ReferenceSolidSurface
{
 public:
    RKSurface(peridynamics::SolidForcePort&p,const StdVec<Vecd>&q,const SnsniReference&r)
      :ReferenceSolidSurface(p,q,1,DeferredBuildTag{})
    {for(const auto&x:q){Row row;for(const auto&e:r.shape(x)){row.index.push_back(e.j);row.weight.push_back(e.w);}rows.push_back(std::move(row));}}
};
void writeFields(const fs::path&path,const SnsniReference&r,const SnsniPlasticSolid&body,const Field&d)
{
    std::ofstream f(path);f<<std::setprecision(17);size_t n=r.cells.size();
    f<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<n<<"\" NumberOfVerts=\""<<n<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t l=0;l<n;++l){Vecd x=r.cells[l].center+r.value(l,d);f<<x[0]<<' '<<x[1]<<' '<<x[2]<<'\n';}
    f<<"</DataArray></Points><Verts><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";for(size_t l=0;l<n;++l)f<<l<<' ';
    f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";for(size_t l=0;l<n;++l)f<<l+1<<' ';f<<"</DataArray></Verts><PointData>";
    auto scalar=[&](const char*name,auto value){f<<"<DataArray type=\"Float64\" Name=\""<<name<<"\" format=\"ascii\">";for(size_t l=0;l<n;++l)f<<value(l)<<' ';f<<"</DataArray>";};
    scalar("Pressure_Pa",[&](size_t l){return body.pending[l].material.p;});
    scalar("J",[&](size_t l){return body.F[l].determinant();});
    scalar("PlasticShear",[&](size_t l){return body.history[l].alpha;});
    scalar("RelativeDensity",[&](size_t l){return body.pending[l].material.relative_density;});
    f<<"</PointData></Piece></PolyData></VTKFile>\n";
}
int qualify(SnsniReference&r,SnsniPlasticSolid&body,Real rho,Real cfl,const fs::path&out)
{
    size_t n=r.X.size();Field d(n),linear,nonlinear;
    for(size_t j=0;j<n;++j){Vecd z=r.X[j]/.0005;d[j]=5e-11*Vecd(z[0]*z[2],-.3*z[1]*z[2],z[2]*z[2]-.2*z[0]*z[0]);}
    r.apply(d,linear);body.evaluate(d,rho,cfl);nonlinear=body.internal;body.discard();
    Field diff(n);for(size_t j=0;j<n;++j)diff[j]=nonlinear[j]-linear[j];Real limit=std::sqrt(dot(diff,diff)/dot(linear,linear));
    // Nonuniform finite deformation, including plastic correction, then rotate.
    for(auto&v:d)v*=50000;
    body.evaluate(d,rho,cfl);nonlinear=body.internal;body.discard();
    Matd Q=Eigen::AngleAxisd(.73,Vecd(1,2,3).normalized()).toRotationMatrix();Field rotated(n);
    for(size_t j=0;j<n;++j)rotated[j]=Q*(r.X[j]+d[j])-r.X[j];
    body.evaluate(rotated,rho,cfl);Vecd net=Vecd::Zero(),torque=Vecd::Zero();Real force_sum=0;
    for(size_t j=0;j<n;++j){diff[j]=body.internal[j]-Q*nonlinear[j];net+=body.internal[j];torque+=(r.X[j]+rotated[j]).cross(body.internal[j]);force_sum+=body.internal[j].norm();}
    Real objective=std::sqrt(dot(diff,diff)/dot(nonlinear,nonlinear));Real net_error=net.norm()/force_sum,torque_error=torque.norm()/(.0005*force_sum);body.discard();
    bool shear=false,cap=false;Real native_error=0;size_t commit_start=body.commits;
    for(int k=0;k<2;++k){Matd F=Matd::Identity();if(k==0)F(0,1)=.015;else F*=.97;
        for(size_t j=0;j<n;++j)d[j]=(F-Matd::Identity())*r.X[j];
        auto old=body.history;body.evaluate(d,rho,cfl);
        for(size_t l=0;l<r.cells.size();++l){auto exact=body.model.update(body.F[l],old[l]);native_error=std::max(native_error,(exact.sigma-body.pending[l].material.sigma).norm()/std::max(Real(1),exact.sigma.norm()));shear|=exact.surface==1;cap|=exact.surface==2;}
        body.commit();}
    bool double_rejected=false;try{body.commit();}catch(const std::logic_error&){double_rejected=true;}
    // One undamped linear dynamics check of the actual central-difference path.
    SnsniCentralDifference time(r.mass);for(size_t j=0;j<n;++j){Vecd z=r.X[j]/.0005;time.displacement[j]=1e-11*Vecd(z[0]*z[2],0,z[2]*z[2]);}
    Field f;r.apply(time.displacement,f);for(size_t j=0;j<n;++j)time.acceleration[j]=-f[j]/r.mass[j];
    Real dt=cfl/std::sqrt(body.linear_frequency_bound_squared);
    auto hamilton=[&](){r.apply(time.displacement,f);Real e=time.kinetic()+.5*dot(time.displacement,f);for(size_t j=0;j<n;++j)e-=.125*dt*dt*r.mass[j]*time.acceleration[j].squaredNorm();return e;};
    Real initial=hamilton(),energy_error=0;
    for(size_t s=0;s<128;++s){time.step(dt,0,[&](const Field&u,const Field&,Field&force){r.apply(u,force);for(auto&v:force)v=-v;});energy_error=std::max(energy_error,std::abs(hamilton()/initial-1));}
    bool pass=limit<2e-5&&objective<1e-7&&net_error<1e-10&&torque_error<1e-9&&native_error==0&&shear&&cap&&double_rejected&&body.commits-commit_start==2&&energy_error<1e-9;
    std::ofstream result(out/"qualification.json");result<<std::setprecision(17)<<"{\n\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\n\"elastic_limit_relative_force_error\":"<<limit<<",\n\"finite_rotation_force_objectivity_error\":"<<objective<<",\n\"net_internal_force_error\":"<<net_error<<",\n\"internal_torque_error\":"<<torque_error<<",\n\"native_DPC_stress_difference\":"<<native_error<<",\n\"shear_and_cap_exercised\":"<<(shear&&cap?"true":"false")<<",\n\"double_commit_rejected\":"<<(double_rejected?"true":"false")<<",\n\"central_difference_invariant_error\":"<<energy_error<<",\n\"linear_frequency_bound_squared\":"<<body.linear_frequency_bound_squared<<",\n\"linear_dt_s\":"<<dt<<",\n\"physical_compression_histories\":0\n}\n";result.close();
    std::cout<<"QUALIFICATION "<<(pass?"PASSED":"FAILED")<<" elastic="<<limit<<" objective="<<objective<<" energy="<<energy_error<<std::endl;return pass?0:2;
}
int run(const PDSolidCapsuleOptions&o,const std::string&geometry,bool gate)
{
    if(o.coated||o.material!="dpc")throw std::invalid_argument("This qualified entry is frozen-input bare DPC; capsule remains gated");
    fs::path out=fs::absolute(o.output);if(fs::exists(out)&&!fs::is_empty(out))throw std::runtime_error("Output must be new or empty");fs::create_directories(out/"vtp");
    Real R=o.diameter/2;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-2*R),Vecd::Constant(2*R)),o.spacing,o.threads);
    SolidBody core(system,makeShared<GeometricShapeBall>(Vecd::Zero(),R,"RKCore"));
    core.defineMatterMaterial<DruckerPragerCapSolid>(o.core_rho,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
    core.generateParticles<BaseParticles,CapsuleCoreGenerator>(R,o.spacing);
    auto&p=core.getBaseParticles();size_t n=p.TotalRealParticles();auto*x=p.getVariableDataByName<Vecd>("Position");auto*v=p.registerStateVariableData<Vecd>("Velocity");auto*m=p.getVariableDataByName<Real>("Mass");auto*volume=p.getVariableDataByName<Real>("VolumetricMeasure");
    Field X(x,x+n);SnsniReference reference(X,o.spacing,o.core_E/(2*(1+o.core_nu)),o.core_E*o.core_nu/((1+o.core_nu)*(1-2*o.core_nu)),o.core_rho,geometry);
    SnsniPlasticSolid body(reference,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
    if(gate)return qualify(reference,body,o.core_rho,o.cfl,out);
    for(size_t j=0;j<n;++j){m[j]=reference.mass[j];volume[j]=m[j]/o.core_rho;}
    peridynamics::SolidForcePort port(core,peridynamics::SolidForcePort::Storage::ForcePrior);
    StdVec<Vecd>points;StdVec<Real>areas;for(const auto&face:reference.faces)if(face.boundary)for(const auto&q:face.q){points.push_back(q.x);areas.push_back(q.w);}
    RKSurface surface(port,points,reference);Real K=o.contact_factor*o.constraint_E/o.spacing;
    peridynamics::SolidSurfacePlaneContact contact(surface,areas,K,R);
    SnsniCentralDifference time(reference.mass);Field total(n,Vecd::Zero());
    const Real ramp=o.ramp_steps*o.max_dt,max_duration=o.max_stage_steps*o.max_dt,damping=o.relaxation_damping/o.max_dt;
    const Real linear_dt=o.cfl/std::sqrt(body.linear_frequency_bound_squared);
    const Real constant_bound=1/std::sqrt(1/(linear_dt*linear_dt)+1/(contact.stable_dt*contact.stable_dt));
    std::ofstream eq(out/"equilibria.csv"),progress(out/"progress.csv"),history(out/"material_history.csv");
    eq<<std::setprecision(17)<<"load_step,relaxation_steps,compression_m,strain,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,min_J\n";
    progress<<std::setprecision(17)<<"step,stage,stage_time_s,dt_s,strain,upper_N,lower_N,kinetic_J,force_residual,min_J,commits\n";
    history<<std::setprecision(17)<<"load_step,strain,core_volume_ratio,min_relative_density,max_relative_density,max_plastic_alpha,max_kappa,plastic_stress_work_J,maximum_yield_residual\n";
    Real stage_time=0;size_t stage=1,stage_steps=0,checks=0;Real last_dt=std::min(o.max_dt,constant_bound);
    auto delta=[&](Real t){Real a=std::clamp(t/ramp,Real(0),Real(1));return o.diameter*o.strain*(Real(stage-1)+.5*(1-std::cos(Pi*a)))/o.load_steps;};
    auto evaluate=[&](const Field&d,const Field&velocity,Real t,Field&force){body.evaluate(d,o.core_rho,o.cfl);
        for(size_t j=0;j<n;++j){x[j]=X[j]+d[j];v[j]=velocity[j];}port.clear();Real u=delta(t);contact.apply(-R+u/2,R-u/2);
        force.resize(n);for(size_t j=0;j<n;++j)force[j]=port.external[j]-body.internal[j];body.commit();};
    evaluate(time.displacement,time.velocity,0,total);for(size_t j=0;j<n;++j)time.acceleration[j]=total[j]/m[j];writeFields(out/"vtp/core_0.vtp",reference,body,time.displacement);
    std::ofstream config(out/"numerics.json");config<<std::setprecision(17)<<"{\"linear_dt_bound_s\":"<<linear_dt<<",\"contact_dt_bound_s\":"<<contact.stable_dt<<",\"combined_bound_s\":"<<constant_bound<<",\"contact_K_Pa_per_m\":"<<K<<",\"ramp_duration_s\":"<<ramp<<",\"relaxation_rate_s_inverse\":"<<damping<<",\"coefficient_dofs\":"<<n<<",\"material_sites\":"<<body.history.size()<<",\"mass_kg\":"<<std::accumulate(m,m+n,Real(0))<<"}\n";config.close();
    auto transfer=[&](recovery::Archive&a){
        a.fields(stage_time,stage,stage_steps,checks,last_dt,total);
        recovery::integrator(a,time);recovery::solid(a,body);recovery::contact(a,contact);
        a.span(x,n);a.span(v,n);
        if(time.displacement.size()!=n||body.history.size()!=reference.cells.size())throw std::runtime_error("Checkpoint geometry size mismatch");
    };
    auto save=[&](){eq.flush();progress.flush();history.flush();recovery::checkpoint.save(out,{"equilibria.csv","progress.csv","material_history.csv"},time.steps,transfer);};
    if(!recovery::checkpoint.resume.empty()){
        eq.close();progress.close();history.close();
        recovery::checkpoint.load(out,{"equilibria.csv","progress.csv","material_history.csv"},transfer);
        eq.open(out/"equilibria.csv",std::ios::app);progress.open(out/"progress.csv",std::ios::app);history.open(out/"material_history.csv",std::ios::app);
        eq<<std::setprecision(17);progress<<std::setprecision(17);history<<std::setprecision(17);
    }
    save();
    while(stage<=o.load_steps){size_t prior_stage=stage;Real dt=std::min({o.max_dt,constant_bound,body.acoustic_dt});if(stage_time<ramp)dt=std::min(dt,ramp-stage_time);
        if(!std::isfinite(dt)||dt<1e-15||time.steps>=o.max_steps)throw std::runtime_error("Explicit time/count bound exhausted");
        time.step(dt,damping,[&](const Field&d,const Field&vel,Field&force){evaluate(d,vel,stage_time+dt,force);total=force;});stage_time+=dt;++stage_steps;last_dt=dt;
        if(time.steps%200==0){Real mean=.5*(contact.upper+contact.lower),residual=std::sqrt(dot(total,total))/std::max(mean,Real(1e-20)),imbalance=std::abs(contact.upper-contact.lower)/std::max(mean,Real(1e-20));
            Real kr=time.kinetic()/std::max(mean*delta(stage_time),Real(1e-30)),minJ=1,vol=0,V0=0,rmin=1,rmax=0,alpha=0,kappa=0,work=0,yr=0;
            for(size_t l=0;l<body.history.size();++l){Real J=body.F[l].determinant(),V=reference.cells[l].volume;minJ=std::min(minJ,J);vol+=V*J;V0+=V;rmin=std::min(rmin,body.pending[l].material.relative_density);rmax=std::max(rmax,body.pending[l].material.relative_density);alpha=std::max(alpha,body.history[l].alpha);kappa=std::max(kappa,body.history[l].kappa);work+=V*body.accumulated_plastic_work[l];yr=std::max(yr,body.pending[l].material.yield_residual);}
            progress<<time.steps<<','<<stage<<','<<stage_time<<','<<dt<<','<<delta(stage_time)/o.diameter<<','<<contact.upper<<','<<contact.lower<<','<<time.kinetic()<<','<<residual<<','<<minJ<<','<<body.commits<<'\n';progress.flush();
            if(time.steps%1000==0)std::cout<<"PROGRESS step="<<time.steps<<" stage="<<stage<<" strain="<<delta(stage_time)/o.diameter<<" force="<<mean<<" residual="<<residual<<" kinetic="<<kr<<std::endl;
            bool good=stage_time>=ramp&&mean>1e-10&&residual<=o.force_tolerance&&imbalance<=o.balance_tolerance&&kr<=o.kinetic_tolerance;checks=good?checks+1:0;
            if(checks>=3){eq<<stage<<','<<stage_steps<<','<<delta(stage_time)<<','<<delta(stage_time)/o.diameter<<','<<contact.upper<<','<<contact.lower<<','<<residual<<','<<imbalance<<','<<kr<<','<<minJ<<'\n';eq.flush();history<<stage<<','<<delta(stage_time)/o.diameter<<','<<vol/V0<<','<<rmin<<','<<rmax<<','<<alpha<<','<<kappa<<','<<work<<','<<yr<<'\n';history.flush();writeFields(out/"vtp"/("core_"+std::to_string(stage)+".vtp"),reference,body,time.displacement);std::cout<<"ACCEPTED stage="<<stage<<" force="<<mean<<std::endl;++stage;stage_time=0;stage_steps=0;checks=0;}
        }
        if(stage_time>max_duration)throw std::runtime_error("Stage did not equilibrate; no completed curve claimed");
        if(time.steps%200==0 && recovery::checkpoint.due(time.steps,stage!=prior_stage))save();
        if(recovery::checkpoint.stop_step && time.steps>=recovery::checkpoint.stop_step){save();return 86;}
    }
    std::ofstream summary(out/"summary.json");summary<<"{\"status\":\"COMPLETE_FROZEN_INPUT_BARE_SNSNI_CURVE\",\"accepted_points\":"<<o.load_steps<<",\"steps\":"<<time.steps<<",\"material_evaluations\":"<<body.evaluations<<",\"history_commits\":"<<body.commits<<",\"material_parameters_changed\":false}\n";return 0;
}
int main(int argc,char**argv)
{try{std::string cells;bool gate=false;std::vector<char*>args{argv[0]};for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--qualify-snsni")gate=true;else if(a.rfind("--snsni-cells=",0)==0)cells=a.substr(14);else if(recovery::checkpoint.take(a)){}else args.push_back(argv[i]);}
    PDSolidCapsuleOptions o;o.parse(int(args.size()),args.data());if(cells.empty())throw std::invalid_argument("--snsni-cells required");tbb::task_arena arena(int(o.threads));return arena.execute([&]{return run(o,cells,gate);});
 }catch(const std::exception&e){std::cerr<<"SNSNI_ERROR "<<e.what()<<std::endl;return 1;}}

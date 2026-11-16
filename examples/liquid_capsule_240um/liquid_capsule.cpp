// Isolated reproduction of Zhang et al., IJMS 205 (2021) 106604, Fig.3d.
// Solver units: micrometre, nanosecond, density 1000 kg/m^3.
// Therefore stress=GPa, force=mN, energy=nJ. No material mass scaling.
#include "sphinxsys.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <csignal>
#include <cstring>
#include <limits>
#include <tbb/scalable_allocator.h>
#include <malloc.h>
using namespace SPH;
namespace fs=std::filesystem;
static volatile std::sig_atomic_t stop_requested=0;
static void on_stop(int){stop_requested=1;}
static long rss_kb(){std::ifstream f("/proc/self/status");std::string s;while(std::getline(f,s))if(s.rfind("VmRSS:",0)==0)return std::stol(s.substr(6));return -1;}
static void compact(ParticleConfiguration &c){
    tbb::parallel_for(tbb::blocked_range<size_t>(0,c.size()),[&](const auto&r){for(size_t i=r.begin();i<r.end();++i){
        auto&n=c[i];n.j_.shrink_to_fit();n.W_ij_.shrink_to_fit();n.dW_ijV_j_.shrink_to_fit();n.r_ij_.shrink_to_fit();n.e_ij_.shrink_to_fit();n.bondLive_.shrink_to_fit();
    }});
}
constexpr double R=120, RI=118, T=2, ES=3.5, NU=.33, SY=.040, RS=1.5, RF=.83, KF=1.66;
constexpr double CS=1.4142135623730951, DS=2./3., DF=2.;

class BallRegion: public Shape {
    double a,b;
public:
    BallRegion(std::string name,double inner,double outer):Shape(name),a(inner),b(outer){}
    bool checkContain(const Vecd& p,bool included=true) override {
        double r=p.norm(); return included?r>=a && r<=b:r>a && r<b;
    }
    Vecd findClosestPoint(const Vecd& p) override {
        double r=p.norm(); Vecd n=Vecd::UnitX();if(r>1e-12)n=p/r;
        double q=a>0 && std::abs(r-a)<std::abs(r-b)?a:b; return q*n;
    }
protected: BoundingBox findBounds() override {return BoundingBox(Vecd::Constant(-b),Vecd::Constant(b));}
};
class ShellGenerator:public ParticleGenerator {
public:
    explicit ShellGenerator(SPHBody& b):ParticleGenerator(b){}
    void initializeGeometricVariables() override {
        // Three ordinary volumetric layers; exact spherical-sector volumes.
        // Particle centres lie half a spacing inside the actual free surfaces.
        // A cut Cartesian lattice can otherwise preload the original wall law.
        const int nt=static_cast<int>(std::ceil(Pi*(R-.5*T)/DS));
        for(int l=0;l<3;++l){double a=RI+l*DS,b=a+DS,r=.5*(a+b);
            for(int j=0;j<nt;++j){double ta=Pi*j/nt,tb=Pi*(j+1)/nt,t=.5*(ta+tb);
                int np=std::max(4,4*static_cast<int>(std::ceil(2*Pi*r*std::sin(t)/(4*DS))));
                double v=(b*b*b-a*a*a)/3*(std::cos(ta)-std::cos(tb))*2*Pi/np;
                for(int k=0;k<np;++k){double p=2*Pi*(k+.5)/np;
                    initializePositionAndVolumetricMeasure(Vecd(r*std::sin(t)*std::cos(p),r*std::cos(t),r*std::sin(t)*std::sin(p)),v);
                }
            }
        }
    }
};
class CoreGenerator:public ParticleGenerator {
public:
    explicit CoreGenerator(SPHBody& b):ParticleGenerator(b){}
    void initializeGeometricVariables() override {
        for(int l=0;l<59;++l){double a=l*DF,b=a+DF,r=.5*(a+b);
            int nt=std::max(2,static_cast<int>(std::ceil(Pi*r/DF)));
            for(int j=0;j<nt;++j){double ta=Pi*j/nt,tb=Pi*(j+1)/nt,t=.5*(ta+tb);
                int np=std::max(4,4*static_cast<int>(std::ceil(2*Pi*r*std::sin(t)/(4*DF))));
                double v=(b*b*b-a*a*a)/3*(std::cos(ta)-std::cos(tb))*2*Pi/np;
                for(int k=0;k<np;++k){double p=2*Pi*(k+.5)/np;
                    initializePositionAndVolumetricMeasure(Vecd(r*std::sin(t)*std::cos(p),r*std::cos(t),r*std::sin(t)*std::sin(p)),v);
                }
            }
        }
    }
};
class PlateGenerator:public ParticleGenerator {
    double y; int sign;
public:
    PlateGenerator(SPHBody& b,double y0,int s):ParticleGenerator(b),y(y0),sign(s){}
    void initializeGeometricVariables() override {
        const int m=150;
        for(int i=-m;i<m;++i)for(int k=-m;k<m;++k)for(int l=0;l<4;++l){
            Vecd x((i+.5)*DS,y+sign*(l+.5)*DS,(k+.5)*DS);
            initializePositionAndVolumetricMeasure(x,DS*DS*DS);
        }
    }
};
class IdealJ2:public J2PlasticityforPD {
public: IdealJ2(double rho,double E,double nu,double sy):J2PlasticityforPD(rho,E,nu,sy,0,0){
    // The upstream constructor hard-codes nonzero Voce saturation even at H=0.
    sigmaY_infinite_=sigmaY_0_;
}
};
class ConsistentContactDensity:public solid_dynamics::ContactDensitySummation {
public:
    explicit ConsistentContactDensity(SurfaceContactRelation&r):ContactDensitySummation(r){
        auto&body=r.getSPHBody();
        for(size_t k=0;k<r.contact_bodies_.size();++k){
            auto&other=*r.contact_bodies_[k];
            // Match NeighborBuilderSurfaceContact's kernel exactly. The old
            // offset used the source PD influence function and a different h.
            KernelWendlandC2 kernel(.5*(body.sph_adaptation_->ReferenceSmoothingLength()+other.sph_adaptation_->ReferenceSmoothingLength()));
            double distance=.5*(body.sph_adaptation_->ReferenceSpacing()+other.sph_adaptation_->ReferenceSpacing());
            offset_W_ij_[k]=kernel.W(distance,ZeroVecd);
        }
    }
};
class PlateRegion:public Shape {
    double y; int sign;
public:
    PlateRegion(std::string n,double a,int s):Shape(n),y(a),sign(s){}
    bool checkContain(const Vecd&p,bool=true) override {return std::abs(p[0])<=100&&std::abs(p[2])<=100&&sign*(p[1]-y)>=0&&sign*(p[1]-y)<=4*DS;}
    Vecd findClosestPoint(const Vecd&p) override {return Vecd(p[0],y,p[2]);}
protected: BoundingBox findBounds() override {return BoundingBox(Vecd(-100,y+(sign<0?-4*DS:0),-100),Vecd(100,y+(sign>0?4*DS:0),100));}
};
class CSV {
    std::ofstream f;
public:
    explicit CSV(const fs::path&p,const std::string&header,bool append):f(p,append?std::ios::app:std::ios::out){
        f<<std::setprecision(15); if(!append)f<<header<<'\n';
    }
    template<class...A>void line(A...a){bool first=true;((f<<(first?(first=false,""):",")<<a),...);f<<'\n';f.flush();}
};
template<class V>void dumpvec(std::ofstream&f,const V&v,size_t n){
    for(size_t i=0;i<n;++i)f.write(reinterpret_cast<const char*>(&v[i]),sizeof(v[i]));
}
template<class V>void readvec(std::ifstream&f,V&v,size_t n){
    for(size_t i=0;i<n;++i)f.read(reinterpret_cast<char*>(&v[i]),sizeof(v[i]));
}
struct History {uint64_t step=0;int stage=1;double time=0,start_time=0,hold_start=-1,top=R,previous_force=0;};

int main(int argc,char**argv){
 try {
    if(argc==3 && std::string(argv[1])=="--material-audit"){
        SPHSystem sys(BoundingBox(Vecd::Constant(-1),Vecd::Constant(1)),1.,1);
        PDBody b(sys,makeShared<BallRegion>("PDBodyAudit",0,1));
        b.defineParticlesAndMaterial<NosbPDPlasticParticles,IdealJ2>(RS,ES,NU,SY);
        b.generateParticles<ObserverParticleGenerator>(StdVec<Vecd>{Vecd::Zero()});
        auto&p=dynamic_cast<NosbPDPlasticParticles&>(b.getBaseParticles());
        auto&m=dynamic_cast<IdealJ2&>(p.plastic_solid_);Matd stress=Matd::Zero(),ep=Matd::Zero(),g=Matd::Zero();
        g(0,0)=.001;g(1,1)=g(2,2)=-.0005;double err=0,last=0;
        for(int k=0;k<200;++k){stress=m.PlasticConstitutiveRelation(g,stress,ep,0,1.);Matd dev=stress-stress.trace()/3.*Matd::Identity();last=std::sqrt(1.5*dev.squaredNorm());if(k>20)err=std::max(err,std::abs(last-SY));}
        bool pass=err<1e-10 && ep.norm()>0;
        std::ofstream f(argv[2]);f<<std::setprecision(16)<<"{\"passed\":"<<(pass?"true":"false")<<",\"max_postyield_equivalent_stress_error_GPa\":"<<err<<",\"final_equivalent_stress_GPa\":"<<last<<",\"target_GPa\":"<<SY<<",\"plastic_strain_norm\":"<<ep.norm()<<"}\n";
        return pass?0:2;
    }
    fs::path out="results"; uint64_t max_steps=0; bool resume=false;size_t threads=6;bool initialize_only=false;
    for(int i=1;i<argc;++i){std::string a=argv[i];
        if(a=="--output"&&i+1<argc)out=argv[++i];
        else if(a=="--max-steps"&&i+1<argc)max_steps=std::stoull(argv[++i]);
        else if(a=="--resume")resume=true;
        else if(a=="--initialize-only")initialize_only=true;
        else if(a=="--threads"&&i+1<argc)threads=std::stoull(argv[++i]);
        else throw std::runtime_error("Unknown argument "+a);
    }
    fs::create_directories(out);std::signal(SIGTERM,on_stop);std::signal(SIGINT,on_stop);
    const auto wall0=std::chrono::steady_clock::now();
    SPHSystem sys(BoundingBox(Vecd::Constant(-132),Vecd::Constant(132)),DF,threads);
    sys.setReloadParticles(false);sys.setRunParticleRelaxation(false);
    IOEnvironment io(sys);
    auto mark=[&](const char*s){auto m=mallinfo2();std::cout<<"MEM "<<s<<" rss_kB="<<rss_kb()<<" glibc_live="<<m.uordblks<<" glibc_free="<<m.fordblks<<std::endl;};
    FluidBody fluid(sys,makeShared<BallRegion>("WaterBody",0,RI));
    fluid.defineParticlesAndMaterial<FluidParticles,WeaklyCompressibleFluid>(RF,CS);
    fluid.generateParticles<CoreGenerator>();
    mark("fluid_generated");
    PDBody shell(sys,makeShared<BallRegion>("PDBodyShell",RI,R));
    shell.defineAdaptation<SPHAdaptation>(1.5075,DF/DS);
    shell.defineParticlesAndMaterial<NosbPDPlasticParticles,IdealJ2>(RS,ES,NU,SY);
    shell.generateParticles<ShellGenerator>();
    mark("shell_generated");
    SolidBody top(sys,makeShared<PlateRegion>("TopPlate",R,1));
    SolidBody bottom(sys,makeShared<PlateRegion>("BottomPlate",-R,-1));
    top.defineAdaptationRatios(1.15,DF/DS);bottom.defineAdaptationRatios(1.15,DF/DS);
    top.defineParticlesAndMaterial<SolidParticles,Solid>(RS,ES/RS);
    bottom.defineParticlesAndMaterial<SolidParticles,Solid>(RS,ES/RS);
    top.generateParticles<PlateGenerator>(R,1);bottom.generateParticles<PlateGenerator>(-R,-1);
    mark("plates_generated");
    auto &sp=dynamic_cast<NosbPDPlasticParticles&>(shell.getBaseParticles());
    auto &iso=*sp.getVariableByName<Real>("IsotropicHardeningParam");
    auto &kin=*sp.getVariableByName<Matd>("KinematicHardeningParam");
    auto &fp=dynamic_cast<FluidParticles&>(fluid.getBaseParticles());
    auto &tp=dynamic_cast<SolidParticles&>(top.getBaseParticles());
    auto &bp=dynamic_cast<SolidParticles&>(bottom.getBaseParticles());
    const size_t ns=sp.total_real_particles_,nf=fp.total_real_particles_;
    std::vector<Vecd> X(ns);for(size_t i=0;i<ns;++i){X[i]=sp.pos_[i];sp.n_[i]=-X[i].normalized();}
    for(size_t i=0;i<tp.total_real_particles_;++i)tp.n_[i]=Vecd(0,-1,0);
    for(size_t i=0;i<bp.total_real_particles_;++i)bp.n_[i]=Vecd(0,1,0);
    InnerRelation shell_inner(shell);
    ComplexRelation fluid_rel(fluid,RealBodyVector{&shell});
    ContactRelation shell_fluid(shell,{&fluid});
    mark("relations_before_plates");
    SurfaceContactRelation shell_top(shell,{&top}),shell_bottom(shell,{&bottom});
    mark("relations_after_plates");
    SimpleDynamics<TimeStepInitialization> fluid_init(fluid);
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWallforPD> fluid_first(fluid_rel);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWallforPD> fluid_second(fluid_rel);
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidforPD> fsi(shell_fluid);
    solid_dynamics::AverageVelocityAndAcceleration average(shell);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> pd_first(shell);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic> pd_second(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> pd_third(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> pd_fourth(shell);
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass(shell_inner,shell.sph_adaptation_->getKernel());
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> damping(shell_inner,shell.sph_adaptation_->getKernel());
    InteractionDynamics<ConsistentContactDensity> top_density(shell_top),bottom_density(shell_bottom);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> top_contact(shell_top),bottom_contact(shell_bottom);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> pd_dt(shell,.10);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_dt(fluid,.20);
    // Fixed particle masses + continuity equation. No free-surface density reset.
    sys.initializeSystemCellLinkedLists();sys.initializeSystemConfigurations();shape.parallel_exec();
    long rss_before=rss_kb();
    compact(shell_inner.inner_configuration_);compact(fluid_rel.inner_configuration_);
    compact(fluid_rel.contact_configuration_[0]);compact(shell_fluid.contact_configuration_[0]);
    compact(shell_top.contact_configuration_[0]);compact(shell_bottom.contact_configuration_[0]);
    scalable_allocation_command(TBBMALLOC_CLEAN_ALL_BUFFERS,nullptr);
    malloc_trim(0);
    {auto m=mallinfo2();std::ofstream f(out/"memory_compaction.json");f<<"{\"rss_before_kB\":"<<rss_before<<",\"rss_after_kB\":"<<rss_kb()<<",\"glibc_live_bytes\":"<<m.uordblks<<",\"glibc_free_bytes\":"<<m.fordblks<<",\"physical_or_discrete_parameters_changed\":false}\n";}
    uint64_t npd=0,nff=0,nfs=0;
    for(size_t i=0;i<ns;++i)npd+=shell_inner.inner_configuration_[i].current_size_;
    for(size_t i=0;i<nf;++i){nff+=fluid_rel.inner_configuration_[i].current_size_;nfs+=fluid_rel.contact_configuration_[0][i].current_size_;}
    {std::ofstream f(out/"initialization.json");f<<"{\"shell_particles\":"<<ns<<",\"fluid_particles\":"<<nf<<",\"plate_particles_each\":"<<tp.total_real_particles_<<",\"PD_directed_neighbors\":"<<npd<<",\"fluid_directed_neighbors\":"<<nff<<",\"fluid_shell_neighbors\":"<<nfs<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-wall0).count()<<"}\n";}
    std::cout<<"Initialized PD="<<ns<<" fluid="<<nf<<" neighbors="<<npd+nff+nfs<<std::endl;
    if(initialize_only)return 0;
    History h;double force_top=0,force_bottom=0;
    auto checkpoint=[&](){
        std::ofstream f(out/"checkpoint.tmp",std::ios::binary);
        f.write("PUF240A1",8);f.write((char*)&ns,sizeof ns);f.write((char*)&nf,sizeof nf);f.write((char*)&h,sizeof h);
        for(auto*v:{&sp.pos_,&sp.vel_,&sp.acc_,&sp.acc_prior_,&sp.acc_old_,&sp.vel_ave_,&sp.acc_ave_})dumpvec(f,*v,ns);
        for(auto*v:{&sp.F_,&sp.F_half_,&sp.stress_,&sp.plastic_strain_,&kin})dumpvec(f,*v,ns);
        dumpvec(f,iso,ns);dumpvec(f,sp.rho_,ns);
        for(auto*v:{&fp.pos_,&fp.vel_,&fp.acc_,&fp.acc_prior_})dumpvec(f,*v,nf);
        dumpvec(f,fp.rho_,nf);dumpvec(f,fp.p_,nf);dumpvec(f,fp.drho_dt_,nf);
        f.close();fs::rename(out/"checkpoint.tmp",out/"checkpoint.bin");
    };
    if(resume){
        std::ifstream f(out/"checkpoint.bin",std::ios::binary);char magic[8];size_t a,b;
        f.read(magic,8);f.read((char*)&a,sizeof a);f.read((char*)&b,sizeof b);
        if(std::memcmp(magic,"PUF240A1",8)||a!=ns||b!=nf)throw std::runtime_error("Checkpoint identity mismatch");
        f.read((char*)&h,sizeof h);
        for(auto*v:{&sp.pos_,&sp.vel_,&sp.acc_,&sp.acc_prior_,&sp.acc_old_,&sp.vel_ave_,&sp.acc_ave_})readvec(f,*v,ns);
        for(auto*v:{&sp.F_,&sp.F_half_,&sp.stress_,&sp.plastic_strain_,&kin})readvec(f,*v,ns);
        readvec(f,iso,ns);readvec(f,sp.rho_,ns);
        for(auto*v:{&fp.pos_,&fp.vel_,&fp.acc_,&fp.acc_prior_})readvec(f,*v,nf);
        readvec(f,fp.rho_,nf);readvec(f,fp.p_,nf);readvec(f,fp.drho_dt_,nf);
        if(!f)throw std::runtime_error("Truncated checkpoint");
        for(size_t i=0;i<tp.total_real_particles_;++i)tp.pos_[i][1]+=h.top-R;
    }
    const std::string header="step,stage,strain,displacement_um,physical_time_s,top_force_mN,bottom_force_mN,force_imbalance,KE_nJ,recoverable_energy_nJ,kinetic_ratio,mean_pressure_MPa,pressure_cv,SPH_volume_um3,shell_volume_um3,p_from_shell_volume_MPa,max_eq_plastic_strain,fluid_max_radius_um,dt_ns,wall_seconds,max_inner_surface_crossing_um,outside_surface_box_particles,pressure_volume_relative_mismatch";
    CSV trace(out/"diagnostics.csv",header,resume),accepted(out/"accepted_curve.csv",header,resume);
    double area_sum0=0;
    for(size_t i=0;i<ns;++i)area_sum0+=sp.Vol_[i]*RI*RI/(T*X[i].squaredNorm());
    const double awscale=4*Pi*RI*RI/area_sum0, v0=4*Pi*RI*RI*RI/3.;
    uint64_t segment_start=h.step;double last_ke=0;int stable_samples=0;
    const double load_duration=1000.,min_hold=500.; // ns; rate-independent quasi-static acceleration, acceptance gated.
    std::vector<Vecd> inner_pos(ns),inner_normal(ns);
    while(h.stage<=10 && !stop_requested){
        double dt=std::min(pd_dt.parallel_exec(),fluid_dt.parallel_exec());
        if(!(dt>0)||!std::isfinite(dt)||dt<1e-8)throw std::runtime_error("Invalid/vanishing stable time step");
        double u=std::clamp((h.time-h.start_time)/load_duration,0.,1.);
        double blend=u*u*(3-2*u), deriv=6*u*(1-u)/load_duration;
        double target=R-4.8*(h.stage-1)-4.8*blend;
        double plate_velocity=-4.8*deriv,shift=target-h.top;
        for(size_t i=0;i<tp.total_real_particles_;++i){tp.pos_[i][1]+=shift;tp.vel_[i]=Vecd(0,plate_velocity,0);}
        h.top=target;
        top.updateCellLinkedList();shell.updateCellLinkedList();fluid.updateCellLinkedList();
        fluid_rel.updateConfiguration();shell_fluid.updateConfiguration();shell_top.updateConfiguration();shell_bottom.updateConfiguration();
        fluid_init.parallel_exec();average.initialize_displacement_.parallel_exec();
        for(size_t i=0;i<ns;++i)sp.acc_prior_[i]=Vecd::Zero();
        fsi.parallel_exec(dt);
        top_density.parallel_exec();
        std::vector<Vecd> before(ns);for(size_t i=0;i<ns;++i)before[i]=sp.acc_prior_[i];
        top_contact.parallel_exec(dt);force_top=0;
        for(size_t i=0;i<ns;++i){force_top-=(sp.acc_prior_[i][1]-before[i][1])*sp.mass_[i];before[i]=sp.acc_prior_[i];}
        bottom_density.parallel_exec();bottom_contact.parallel_exec(dt);force_bottom=0;
        for(size_t i=0;i<ns;++i)force_bottom+=(sp.acc_prior_[i][1]-before[i][1])*sp.mass_[i];
        pd_first.parallel_exec(dt);pd_second.parallel_exec(dt);hourglass.parallel_exec(dt);damping.parallel_exec(dt);pd_third.parallel_exec(dt);pd_fourth.parallel_exec(dt);
        average.update_averages_.parallel_exec(dt);
        // On a thin curved shell the virtual fluid normal must follow its inner face.
        for(size_t i=0;i<ns;++i){sp.n_[i]=-sp.F_[i].inverse().transpose()*X[i].normalized();sp.n_[i].normalize();}
        fluid_first.parallel_exec(dt);fluid_second.parallel_exec(dt);
        h.time+=dt;++h.step;GlobalStaticVariables::physical_time_=h.time;
        if(u>=1 && h.hold_start<0)h.hold_start=h.time;
        if(h.step%100==0 || (h.step<=100&&h.step%10==0)){
            double ke=0,ue=0,pmean=0,p2=0,vf=0,vsurf=0,eqp=0,rmax=0,crossing=0;
            uint64_t outside_box=0;Vecd box_min=Vecd::Constant(1e30),box_max=Vecd::Constant(-1e30);
            const double G=ES/(2*(1+NU)),K=ES/(3*(1-2*NU));
            for(size_t i=0;i<ns;++i){
                double J=sp.F_[i].determinant();if(!(J>0)||!std::isfinite(J))throw std::runtime_error("PD inverted/nonfinite deformation");
                ke+=.5*sp.mass_[i]*sp.vel_[i].squaredNorm();
                double tr=sp.stress_[i].trace();Matd dev=sp.stress_[i]-Matd::Identity()*(tr/3.);
                ue+=sp.Vol_[i]*(dev.squaredNorm()/(4*G)+tr*tr/(18*K));
                eqp=std::max(eqp,iso[i]);
                Vecd n=X[i].normalized();Vecd xi=sp.pos_[i]+sp.F_[i]*(RI*n-X[i]);
                Vecd da=J*sp.F_[i].inverse().transpose()*n;
                inner_pos[i]=xi;inner_normal[i]=-da.normalized();box_min=box_min.cwiseMin(xi);box_max=box_max.cwiseMax(xi);
                vsurf+=xi.dot(da)*sp.Vol_[i]*RI*RI/(T*X[i].squaredNorm())*awscale/3.;
            }
            for(size_t i=0;i<nf;++i){
                if(!(fp.rho_[i]>0)||!std::isfinite(fp.rho_[i]))throw std::runtime_error("Invalid liquid density");
                double v=fp.mass_[i]/fp.rho_[i];vf+=v;pmean+=v*fp.p_[i];p2+=v*fp.p_[i]*fp.p_[i];
                ke+=.5*fp.mass_[i]*fp.vel_[i].squaredNorm();ue+=.5*v*fp.p_[i]*fp.p_[i]/KF;
                rmax=std::max(rmax,fp.pos_[i].norm());
                auto&nb=fluid_rel.contact_configuration_[0][i];double best=1e30,margin=1e30;
                for(size_t n=0;n<nb.current_size_;++n){size_t j=nb.j_[n];Vecd delta=fp.pos_[i]-inner_pos[j];double q=delta.squaredNorm();if(q<best){best=q;margin=delta.dot(inner_normal[j]);}}
                if(best<1e29)crossing=std::max(crossing,-margin);
                if((fp.pos_[i].array()<box_min.array()-.2).any()||(fp.pos_[i].array()>box_max.array()+.2).any())++outside_box;
            }
            pmean/=vf;double cv=std::sqrt(std::max(0.,p2/vf-pmean*pmean))/std::max(std::abs(pmean),1e-9);
            double imb=std::abs(force_top-force_bottom)/std::max({std::abs(force_top),std::abs(force_bottom),1e-5});
            double kr=ke/std::max(ue,1e-12),wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-wall0).count();
            double strain=(R-h.top)/(2*R),pv=-KF*(vsurf-v0)/v0;
            double pverr=std::abs(pmean-pv)/std::max({std::abs(pmean),std::abs(pv),1e-5});
            trace.line(h.step,h.stage,strain,R-h.top,h.time*1e-9,force_top,force_bottom,imb,ke,ue,kr,pmean*1000,cv,vf,vsurf,pv*1000,eqp,rmax,dt,wall,crossing,outside_box,pverr);
            std::ofstream state(out/"status.tmp");state<<std::setprecision(15)<<"{\"state\":\"RUNNING_UNACCEPTED\",\"step\":"<<h.step<<",\"stage\":"<<h.stage<<",\"strain\":"<<strain<<",\"force_mN\":"<<force_top<<",\"kinetic_ratio\":"<<kr<<",\"force_imbalance\":"<<imb<<",\"wall_seconds\":"<<wall<<"}\n";state.close();fs::rename(out/"status.tmp",out/"status.json");
            std::cout<<"step="<<h.step<<" strain="<<strain<<" F_mN="<<force_top<<" imbalance="<<imb<<" KE/U="<<kr<<" elapsed="<<wall<<std::endl;
            bool held=h.hold_start>=0 && h.time-h.hold_start>=min_hold;
            double drift=std::abs(force_top-h.previous_force)/std::max(std::abs(force_top),1e-5);
            bool settled=held&&imb<.01&&kr<.01&&drift<.005&&cv<.05;
            stable_samples=settled?stable_samples+1:0;h.previous_force=force_top;last_ke=kr;
            if(stable_samples>=5){
                if(pverr>.05 || crossing>.2 || outside_box){checkpoint();throw std::runtime_error("Model gate failed: pressure-volume relation or resolved inner-surface crossing; no curve point accepted");}
                accepted.line(h.step,h.stage,strain,R-h.top,h.time*1e-9,force_top,force_bottom,imb,ke,ue,kr,pmean*1000,cv,vf,vsurf,pv*1000,eqp,rmax,dt,wall,crossing,outside_box,pverr);
                ++h.stage;h.start_time=h.time;h.hold_start=-1;stable_samples=0;checkpoint();
            }
            if(held&&h.time-h.hold_start>10000.){checkpoint();throw std::runtime_error("Stage did not equilibrate; no curve point accepted");}
        }
        if(h.step%1000==0)checkpoint();
        if(max_steps && h.step-segment_start>=max_steps)break;
    }
    checkpoint();std::ofstream f(out/"exit.json");f<<"{\"completed\":"<<(h.stage>10?"true":"false")<<",\"accepted_points\":"<<h.stage-1<<",\"step\":"<<h.step<<",\"kinetic_ratio\":"<<last_ke<<",\"checkpoint_preserved\":true}\n";
    return 0;
 }catch(const std::exception&e){std::cerr<<"FAILED: "<<e.what()<<std::endl;return 2;}
}

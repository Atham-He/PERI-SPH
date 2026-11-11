/** NOSB-PD volume shell + total-Lagrangian solid-SPH DPC core.
 * SI units. A working numerical case, not experimentally calibrated.
 * Three/five Gaussian thickness layers retain real shell thickness and mass.
 * They are volume quadrature, NOT a validated reduced-order PD shell theory.
 */
#include "capsule_core_fixture.h"
#include "all_peridynamics.h"
#include "pd_quadratic_correspondence.h"
#include "pd_solid_reference_coupling.h"
#include "pd_principal_stress_onset.h"
#include <sstream>
#include <map>

struct HybridOptions : PDSolidCapsuleOptions
{
    bool preflight_only=false;
    bool quadratic_kinematics=false;
    size_t shell_layers=3;
    unsigned interface_order=1;
    Real interface_core_support=2.5,interface_shell_support=2.5;
    std::string shell_material="elastic",fracture_model="none",fracture_scope="none";
    Real shell_yield=0,shell_hardening=0,fracture_stress=0;
    Real pd_spacing=0,pd_cfl=.1,pd_hourglass=.5;
    Real ramp_duration=2.8132387884274494e-6,max_stage_duration=1.4066193942137247e-5;
    Real damping_rate=.002/2.8132387884274494e-10;
    void parseHybrid(int argc,char **argv)
    {
        std::vector<std::string> base{argv[0]};
        std::function<void(const std::string &)> consume=[&](const std::string &arg)
        {
            if(arg=="--preflight-only"){preflight_only=true;return;}
            auto eq=arg.find('=');if(eq==std::string::npos) {base.push_back(arg);return;}
            std::string key=arg.substr(2,eq-2),value=arg.substr(eq+1);
            if(key=="config")
            {
                std::ifstream f(value);if(!f)throw std::invalid_argument("Cannot read hybrid config");
                std::string line;while(std::getline(f,line))
                {line=line.substr(0,line.find('#'));auto a=line.find_first_not_of(" \t\r");
                 if(a==std::string::npos)continue;line=line.substr(a,line.find_last_not_of(" \t\r")-a+1);consume("--"+line);}return;
            }
            if(key=="pd-kinematics"){if(value!="linear" && value!="quadratic")throw std::invalid_argument("pd-kinematics must be linear or quadratic");quadratic_kinematics=value=="quadratic";return;}
            if(key=="interface-order"){if(value!="1" && value!="2")throw std::invalid_argument("interface-order must be 1 or 2");interface_order=unsigned(value[0]-'0');return;}
            if(key=="core-backend"){if(value!="solid_sph")throw std::invalid_argument("Core requires solid_sph,not a WCSPH fluid");return;}
            if(key=="shell-backend"){if(value!="nosb_pd")throw std::invalid_argument("Shell requires nosb_pd");return;}
            if(key=="shell-material"){shell_material=value;return;}
            if(key=="fracture-model"){fracture_model=value;return;}
            if(key=="fracture-scope"){fracture_scope=value;return;}
            if(key=="shell-layers"){size_t used=0;shell_layers=std::stoul(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid layer count");return;}
            Real *target=nullptr;
            if(key=="shell-yield-stress")target=&shell_yield;
            if(key=="shell-hardening-modulus")target=&shell_hardening;
            if(key=="fracture-stress")target=&fracture_stress;
            if(key=="interface-core-support-ratio")target=&interface_core_support;
            if(key=="interface-shell-support-ratio")target=&interface_shell_support;
            if(key=="pd-spacing")target=&pd_spacing;
            if(key=="pd-cfl")target=&pd_cfl;
            if(key=="pd-hourglass")target=&pd_hourglass;
            if(key=="ramp-duration")target=&ramp_duration;
            if(key=="max-stage-duration")target=&max_stage_duration;
            if(key=="damping-rate")target=&damping_rate;
            if(target){size_t used=0;*target=std::stod(value,&used);if(used!=value.size() || !std::isfinite(*target))throw std::invalid_argument("Invalid hybrid number");return;}
            base.push_back(arg);
        };
        for(int i=1;i<argc;++i)consume(argv[i]);
        std::vector<char*> ptr;for(auto &a:base)ptr.push_back(a.data());
        PDSolidCapsuleOptions::parse(int(ptr.size()),ptr.data());
        if(pd_spacing==0)pd_spacing=spacing;
        if(shell_material!="elastic"&&shell_material!="j2")throw std::invalid_argument("shell-material must be elastic or j2");
        if(shell_yield<0||shell_hardening<0||(shell_material=="j2"&&shell_yield<=0)||
           (shell_material=="elastic"&&(shell_yield!=0||shell_hardening!=0)))
            throw std::invalid_argument("J2 requires an explicit positive shell-yield-stress; elastic forbids plastic parameters");
        if(fracture_model=="none")
        {if(fracture_scope!="none"||fracture_stress!=0)throw std::invalid_argument("No-fracture mode forbids active fracture parameters");}
        else if(fracture_model!="principal_stress"||fracture_scope!="first_crack"||fracture_stress<=0)
            throw std::invalid_argument("Only principal_stress with explicit threshold and fracture-scope=first_crack is supported; post-fracture reconstruction/contact are unqualified");
        if(!coated || !quasistatic || material!="dpc" || plate_contact!="outer_surface" ||
            patch_test || rest_test || smoke_test || (shell_layers!=3 && shell_layers!=5) ||
            pd_spacing<=0 || pd_cfl<=0 || pd_cfl>.2 || pd_hourglass<0 ||
            ramp_duration<=0 || max_stage_duration<=ramp_duration || damping_rate<0 ||
            interface_core_support<=1 || interface_shell_support<=1)
            throw std::invalid_argument("Invalid PD-solid configuration or unsupported dynamic/fracture path");
    }
};

struct ThicknessQuadrature
{
    std::vector<Real> xi,weight;
    explicit ThicknessQuadrature(size_t layers)
    {
        if(layers==3){xi={-std::sqrt(Real(3)/5),0,std::sqrt(Real(3)/5)};weight={Real(5)/9,Real(8)/9,Real(5)/9};}
        else if(layers==5){xi={-.906179845938664,-.538469310105683,0,.538469310105683,.906179845938664};
            weight={.236926885056189,.478628670499366,.568888888888889,.478628670499366,.236926885056189};}
        else throw std::invalid_argument("Thickness quadrature requires3or5 layers");
    }
};

namespace SPH
{
class LayeredPDCapsuleGenerator;
template<> class ParticleGenerator<BaseParticles,LayeredPDCapsuleGenerator> : public ParticleGenerator<BaseParticles>
{
    Real inner_,thickness_,spacing_;size_t layers_;
  public:
    ParticleGenerator(SPHBody &body,BaseParticles &particles,Real inner,Real thickness,Real spacing,size_t layers)
        : ParticleGenerator<BaseParticles>(body,particles),inner_(inner),thickness_(thickness),spacing_(spacing),layers_(layers){}
    void prepareGeometricData() override
    {
        auto directions=sphereDirections(inner_+thickness_,spacing_);ThicknessQuadrature rule(layers_);
        Real mid=inner_+.5*thickness_,solid_angle=4*Pi/directions.size();
        for(size_t l=0;l<rule.xi.size();++l)
        {
            Real r=mid+.5*thickness_*rule.xi[l];
            Real volume=solid_angle*r*r*.5*thickness_*rule.weight[l];
            for(const auto &n:directions)addPositionAndVolumetricMeasure(r*n,volume);
        }
    }
};
}

/** Retain the native PD kick/work recording; apply the SAME relaxation rate
 * before the final drift as the SPH core. This is numerical relaxation only. */
class RelaxedPDFourthStep : public solid_dynamics::NosbPDFourthStep
{
    Real rate_;Real *loss_;
  public:
    RelaxedPDFourthStep(SPHBody &body,Real rate):NosbPDFourthStep(body),rate_(rate),
        loss_(body.getBaseParticles().registerStateVariableData<Real>("PDGlobalRelaxationLoss")){}
    void update(size_t i,Real dt=0)
    {
        Vecd midpoint=pos_[i];NosbPDFourthStep::update(i,dt);
        Real before=.5*mass_[i]*vel_[i].squaredNorm(),factor=std::exp(-rate_*dt);
        vel_[i]*=factor;loss_[i]+=before*(1-factor*factor);pos_[i]=midpoint+.5*dt*vel_[i];
    }
};

int runHybrid(const HybridOptions &o)
{
    using namespace peridynamics;
    fs::path out=fs::absolute(o.output);
    if(fs::exists(out) && !fs::is_empty(out))throw std::runtime_error("Output must be new or empty;no automatic history overwrite");
    fs::create_directories(out/"runtime");
    struct Cwd {fs::path old=fs::current_path();explicit Cwd(const fs::path &p){fs::current_path(p);}~Cwd(){fs::current_path(old);}} cwd(out/"runtime");
    Real R=o.diameter/2,rc=R-o.thickness;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-2*R),Vecd::Constant(2*R)),o.spacing,o.threads);
    IO::getEnvironment().resetOutputFolder((out/"vtp").string(),true);
    SolidBody core(system,makeShared<GeometricShapeBall>(Vecd::Zero(),rc,"Core"));
    core.defineAdaptationRatios(1.15,1);
    core.defineMatterMaterial<DruckerPragerCapSolid>(o.core_rho,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
    auto &dpc_material=DynamicCast<DruckerPragerCapSolid>(&core,core.getMatterMaterial());
    core.generateParticles<BaseParticles,CapsuleCoreGenerator>(rc,o.spacing);
    InnerRelation core_inner(core);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> core_correction(core_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> core_first(core_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> core_second(core_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> core_dt(core,o.cfl);

    PDBody shell(system,makeShared<DefaultShape>("PDShell"));
    shell.defineAdaptation<PDAdaptation>(1.5075,o.spacing/o.pd_spacing);
    if(o.shell_material=="elastic")shell.defineMatterMaterial<HughesWingetSolid>(o.shell_rho,o.shell_E,o.shell_nu);
    else shell.defineMatterMaterial<J2PlasticityforPD>(o.shell_rho,o.shell_E,o.shell_nu,o.shell_yield,
        PDIsotropicHardening::Linear(o.shell_hardening),0);
    shell.generateParticles<BaseParticles,LayeredPDCapsuleGenerator>(rc,o.thickness,o.pd_spacing,o.shell_layers);
    InnerRelation shell_inner(shell);
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> pd_initialize(shell,makeShared<Gravity>(Vecd::Zero()));
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> pd_shape(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> pd_first(shell);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> pd_linear_second(shell_inner);
    solid_dynamics::PDQuadraticReference quadratic_reference(shell_inner);
    InteractionWithUpdate<solid_dynamics::NosbPDQuadraticSecondStep> pd_quadratic_second(shell_inner,quadratic_reference);
    UniquePtr<BaseDynamics<void>> plastic_second;
    if(o.shell_material=="j2")
    {
        if(o.quadratic_kinematics)plastic_second=makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDQuadraticSecondStepPlastic>>(shell_inner,quadratic_reference);
        else plastic_second=makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic>>(shell_inner);
    }
    auto pd_second=[&](Real dt){if(plastic_second)plastic_second->exec(dt);else if(o.quadratic_kinematics)pd_quadratic_second.exec(dt);else pd_linear_second.exec(dt);};
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> pd_linear_third(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDQuadraticThirdStep> pd_quadratic_third(shell_inner,quadratic_reference);
    auto pd_third=[&](Real dt){if(o.quadratic_kinematics)pd_quadratic_third.exec(dt);else pd_linear_third.exec(dt);};
    InteractionDynamics<solid_dynamics::ProjectedPairStabilization> pd_stabilize(shell_inner,o.pd_hourglass,.1);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> pd_pair_damping(shell_inner);
    SimpleDynamics<RelaxedPDFourthStep> pd_fourth(shell,o.damping_rate);
    ReduceDynamics<solid_dynamics::PDTimeStep> pd_dt(shell,o.pd_cfl);
    State c(core),s(shell);
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();core_correction.exec();pd_shape.exec();
    if(o.quadratic_kinematics)quadratic_reference.build();
    pd_second(0);
    SolidForcePort core_port(core,SolidForcePort::Storage::ForcePrior);
    SolidForcePort pd_port(shell,SolidForcePort::Storage::PDAccelerationPrior);
    auto interface_directions=sphereDirections(R,o.pd_spacing);
    auto contact_directions=sphereDirections(R,o.contact_quadrature_scale*o.spacing);
    StdVec<Vecd> inner_points,outer_points;
    for(const auto &n:interface_directions)inner_points.push_back(rc*n);
    for(const auto &n:contact_directions)outer_points.push_back(R*n);
    ReferenceSolidSurface core_surface(core_port,inner_points,o.spacing,o.interface_core_support,o.interface_order);
    ReferenceSolidSurface pd_inner(pd_port,inner_points,o.pd_spacing,o.interface_shell_support,o.interface_order);
    ReferenceSolidSurface pd_outer(pd_port,outer_points,o.pd_spacing);
    Real Eref=o.constraint_E>0?o.constraint_E:o.core_E;
    StdVec<Real> inner_area(inner_points.size(),4*Pi*rc*rc/inner_points.size());
    StdVec<Real> outer_area(outer_points.size(),4*Pi*R*R/outer_points.size());
    BondedSolidInterface tie(pd_inner,core_surface,inner_area,o.tie_factor*Eref/o.spacing,R);
    SolidSurfacePlaneContact contact(pd_outer,outer_area,o.contact_factor*Eref/o.spacing,R);
    Real constraint_dt=1/std::sqrt(1/(tie.stable_dt*tie.stable_dt)+1/(contact.stable_dt*contact.stable_dt));
    ThicknessQuadrature quadrature(o.shell_layers);
    Real radial_spacing=o.thickness;
    for(size_t i=1;i<quadrature.xi.size();++i)radial_spacing=std::min(radial_spacing,.5*o.thickness*(quadrature.xi[i]-quadrature.xi[i-1]));
    auto &pd_material=DynamicCast<HughesWingetSolid>(&shell,shell.getMatterMaterial());
    Real radial_dt=o.pd_cfl*std::min(radial_spacing,o.pd_spacing)/pd_material.ReferenceSoundSpeed();
    auto &cp=core.getBaseParticles();auto &sp=shell.getBaseParticles();
    auto *pd_acc=sp.getVariableDataByName<Vecd>("PDAcceleration");
    auto *pd_internal=sp.getVariableDataByName<Vecd>("PDInternalForce");
    auto *pd_K=sp.getVariableDataByName<Matd>("ShapeTensorK");
    auto *pd_stress=sp.getVariableDataByName<Matd>("PDStress");
    auto *failure=sp.getVariableDataByName<int>("NumericalFailureReason");
    auto *damage=sp.getVariableDataByName<Real>("Damage");
    UniquePtr<solid_dynamics::PrincipalStressBondOnset> onset;
    if(o.fracture_model!="none")onset=makeUnique<solid_dynamics::PrincipalStressBondOnset>(shell_inner,o.fracture_stress);
    auto *dpc_r=cp.getVariableDataByName<Real>("DPCRelativeDensity");
    auto *dpc_kappa=cp.getVariableDataByName<Real>("DPCPlasticCompaction");
    auto *dpc_residual=cp.getVariableDataByName<Real>("DPCYieldResidual");
    auto *dpc_work=cp.getVariableDataByName<Real>("DPCPlasticStressWork");
    auto *dpc_surface=cp.getVariableDataByName<Real>("DPCActiveSurface");
    std::vector<bool> ever_cap(c.n,false);
    Real shape_condition=0,initial_F_error=0;
    for(size_t i=0;i<s.n;++i)
    {
        Eigen::SelfAdjointEigenSolver<Matd> eig(pd_K[i]);
        if(eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff()<=0)throw std::runtime_error("Invalid physical-thickness PD shape tensor");
        shape_condition=std::max(shape_condition,eig.eigenvalues().maxCoeff()/eig.eigenvalues().minCoeff());
        initial_F_error=std::max(initial_F_error,(s.F[i]-Matd::Identity()).norm());
    }
    Real expected_core_mass=o.core_rho*4*Pi*rc*rc*rc/3;
    Real expected_shell_mass=o.shell_rho*4*Pi*(R*R*R-rc*rc*rc)/3;
    if(std::abs(c.totalMass()/expected_core_mass-1)>1e-10 || std::abs(s.totalMass()/expected_shell_mass-1)>1e-10 || initial_F_error>1e-9)
        throw std::runtime_error("Physical mass or initial PD affine kinematics failed");
    std::ofstream config(out/"parameters.txt");config<<std::setprecision(17)
        <<"status=CONFIGURED_PD_SOLID_CASE_NOT_CALIBRATED\ncore_backend=total_lagrangian_solid_SPH\nshell_backend=NOSB_PD_volume\nfluid_backend=none\n"
        <<"shell_constitutive="<<(o.shell_material=="elastic"?"HughesWinget_elastic":"J2_HW_plastic")
        <<"\nfracture_model="<<o.fracture_model<<"\nfracture_scope="<<o.fracture_scope<<"\nfracture_stress="<<o.fracture_stress
        <<"\nshell_yield_stress="<<o.shell_yield<<"\nshell_hardening_modulus="<<o.shell_hardening<<"\ncore_material=DPC\n"
        <<"shell_E="<<o.shell_E<<"\nshell_nu="<<o.shell_nu<<"\nshell_reference_density="<<o.shell_rho
        <<"\ninterface_law=perfect_bond_vector_penalty\n"
        <<"core_particles="<<c.n<<"\nshell_particles="<<s.n<<"\nshell_layers="<<o.shell_layers
        <<"\ncore_diameter="<<2*rc<<"\nouter_diameter="<<2*R<<"\nshell_thickness="<<o.thickness
        <<"\ncore_spacing="<<o.spacing<<"\npd_tangential_spacing="<<o.pd_spacing<<"\nminimum_PD_radial_spacing="<<radial_spacing
        <<"\nPD_horizon="<<shell.getSPHAdaptation().getKernel()->CutOffRadius()
        <<"\ncore_mass="<<c.totalMass()<<"\nshell_mass="<<s.totalMass()<<"\nexpected_shell_mass="<<expected_shell_mass
        <<"\nPD_shape_condition_max="<<shape_condition<<"\ninitial_PD_F_error="<<initial_F_error
        <<"\ninterface_points="<<inner_points.size()<<"\ncontact_points="<<outer_points.size()
        <<"\ncore_surface_affine_error="<<core_surface.affine_error<<"\nPD_inner_affine_error="<<pd_inner.affine_error<<"\nPD_outer_affine_error="<<pd_outer.affine_error
        <<"\ninterface_order="<<o.interface_order<<"\ninterface_core_support_ratio="<<o.interface_core_support<<"\ninterface_shell_support_ratio="<<o.interface_shell_support
        <<"\ncore_surface_polynomial_error="<<core_surface.polynomial_error<<"\nPD_inner_polynomial_error="<<pd_inner.polynomial_error
        <<"\ninterface_K="<<o.tie_factor*Eref/o.spacing<<"\ncontact_K="<<o.contact_factor*Eref/o.spacing
        <<"\ninterface_dt_bound="<<tie.stable_dt<<"\ncontact_dt_bound="<<contact.stable_dt
        <<"\ncombined_constraint_dt="<<constraint_dt<<"\nPD_radial_acoustic_dt="<<radial_dt
        <<"\nmaximum_dt_cap="<<o.max_dt<<"\nnumerical_ramp_duration="<<o.ramp_duration<<"\nnumerical_damping_rate="<<o.damping_rate
        <<"\nload_steps="<<o.load_steps<<"\nstrain="<<o.strain<<"\nthreads="<<o.threads
        <<"\npd_kinematics="<<(o.quadratic_kinematics?"quadratic_reproducing":"linear_correspondence")
        <<"\nquadratic_moment_condition_max="<<quadratic_reference.max_moment_condition<<"\nquadratic_reproduction_error="<<quadratic_reference.max_reproduction_error
        <<"\nquadratic_initial_elastic_dt_bound="<<quadratic_reference.elastic_dt_bound
        <<"\ndpc_shear_hardening="<<dpc_material.ShearHardeningModulus()
        <<"\nsource_sha256="<<CAPSULE_SOURCE_SHA256<<"\ndpc_material_sha256="<<DPC_MATERIAL_SHA256<<'\n';config.flush();
    if(o.preflight_only)
    {
        // Nonzero manufactured loads: check the implemented force ports and
        // virtual work against a finite difference of the actual tie/contact energy.
        Matd affine=Matd::Identity();affine(0,0)=.998;affine(1,1)=1.001;affine(2,2)=.999;
        for(size_t i=0;i<c.n;++i){c.x[i]=affine*c.x0[i]+Vecd(5e-8,0,0);c.v[i]=Vecd(.2*std::sin(37*c.x0[i][0]/R),-.1,.3*c.x0[i][2]/R);}
        for(size_t i=0;i<s.n;++i){s.x[i]=affine*s.x0[i];s.v[i]=Vecd(-.1,.2*std::cos(19*s.x0[i][1]/R),-.2*s.x0[i][2]/R);}
        pd_second(0);Real affine_PD_error=0;
        for(size_t i=0;i<s.n;++i)affine_PD_error=std::max(affine_PD_error,(s.F[i]-affine).norm());
        auto energy_and_forces=[&](){core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+2e-6,R-2e-6);return tie.energy+contact.energy;};
        energy_and_forces();Real force_power=0;
        for(size_t i=0;i<c.n;++i)force_power+=core_port.external[i].dot(c.v[i]);
        for(size_t i=0;i<s.n;++i)force_power+=s.mass[i]*pd_port.external[i].dot(s.v[i]);
        const Real epsilon=1e-12;
        for(size_t i=0;i<c.n;++i)c.x[i]+=epsilon*c.v[i];
        for(size_t i=0;i<s.n;++i)s.x[i]+=epsilon*s.v[i];
        Real plus=energy_and_forces();
        for(size_t i=0;i<c.n;++i)c.x[i]-=2*epsilon*c.v[i];
        for(size_t i=0;i<s.n;++i)s.x[i]-=2*epsilon*s.v[i];
        Real minus=energy_and_forces(),energy_rate=(plus-minus)/(2*epsilon);
        Real finite_difference_error=std::abs(energy_rate+force_power)/std::max({std::abs(force_power),std::abs(energy_rate),Real(1e-20)});
        if(affine_PD_error>1e-9 || finite_difference_error>1e-5 ||
            std::max({tie.worst_force,tie.worst_torque,tie.worst_power,contact.worst_force,contact.worst_torque,contact.worst_power})>1e-9)
            throw std::runtime_error("Manufactured affine PD/loaded interface virtual-work check failed");
        std::ofstream f(out/"preflight.json");f<<std::setprecision(17)
            <<"{\"status\":\"INITIAL_GEOMETRY_AFFINE_AND_LOADED_FORCE_MAPS_PASSED\",\"core_particles\":"<<c.n<<",\"PD_particles\":"<<s.n
            <<",\"PD_shape_condition_max\":"<<shape_condition<<",\"initial_F_error\":"<<initial_F_error
            <<",\"manufactured_PD_affine_error\":"<<affine_PD_error<<",\"loaded_energy_gradient_relative_error\":"<<finite_difference_error
            <<",\"interface_force_residual\":"<<tie.worst_force<<",\"interface_torque_residual\":"<<tie.worst_torque<<",\"interface_power_residual\":"<<tie.worst_power
            <<",\"contact_force_residual\":"<<contact.worst_force<<",\"contact_torque_residual\":"<<contact.worst_torque<<",\"contact_power_residual\":"<<contact.worst_power
            <<",\"full_bending_validation\":false,\"physical_history_points\":0}\n";
        std::cout<<"PREFLIGHT geometry and affine reference maps passed;no compression history evolved."<<std::endl;return 0;
    }
    for(auto *body:{static_cast<RealBody*>(&core),static_cast<RealBody*>(&shell)})
    {body->getBaseParticles().addVariableToWrite<Vecd>("Velocity");body->getBaseParticles().addVariableToWrite<Matd>("DeformationGradient");}
    sp.addVariableToWrite<Matd>("PDStress");
    sp.addVariableToWrite<Real>("Damage");sp.addVariableToWrite<int>("BondNumber");sp.addVariableToWrite<int>("ParticleLive");
    cp.addVariableToWrite<Matd>("DPCCauchyStress");
    BodyStatesRecordingToVtp vtp(system);
    std::ofstream eq(out/"equilibria.csv"),history(out/"material_history.csv"),progress(out/"progress.csv");
    eq<<std::setprecision(17)<<"load_step,relaxation_steps,compression_m,strain,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,max_tie_gap_m,min_core_J,min_PD_J,interface_force_residual,interface_torque_residual,interface_power_residual\n";
    history<<std::setprecision(17)<<"load_step,strain,core_volume_ratio,min_relative_density,max_relative_density,min_hardening_strain,max_hardening_strain,ever_cap_fraction,plastic_stress_work_J,max_yield_residual\n";
    progress<<std::setprecision(17)<<"step,stage,stage_time_s,dt_s,compression_m,upper_N,lower_N,kinetic_J,interface_J,contact_J,min_core_J,min_PD_J\n";
    Real &time=*system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t stage=1,stage_step=0,step=0,consecutive=0;Real stage_time=0,max_gap=0,min_core_J=1,min_pd_J=1,work=0,core_relax_loss=0;
    auto delta_at=[&](Real t){Real a=std::clamp(t/o.ramp_duration,Real(0),Real(1));
        return o.strain*o.diameter*(Real(stage-1)+.5*(1-std::cos(Pi*a)))/o.load_steps;};
    auto apply=[&](Real delta)
    {core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+.5*delta,R-.5*delta);};
    auto write_fields=[&](size_t frame)
    {
        vtp.writeToFile(frame);
        fs::create_directories(out/"interface_surface");std::ofstream f(out/"interface_surface"/("interface_"+std::to_string(frame)+".csv"));
        f<<std::setprecision(17)<<"x0_m,y0_m,z0_m,PD_x_m,PD_y_m,PD_z_m,core_x_m,core_y_m,core_z_m,fx_N,fy_N,fz_N,area_m2\n";
        for(size_t i=0;i<inner_points.size();++i)f<<inner_points[i][0]<<','<<inner_points[i][1]<<','<<inner_points[i][2]<<','<<tie.x_a[i][0]<<','<<tie.x_a[i][1]<<','<<tie.x_a[i][2]<<','
            <<tie.x_b[i][0]<<','<<tie.x_b[i][1]<<','<<tie.x_b[i][2]<<','<<tie.force_a[i][0]<<','<<tie.force_a[i][1]<<','<<tie.force_a[i][2]<<','<<inner_area[i]<<'\n';
    };
    auto inspect=[&](State &a,Real &kinetic,Real &minimum,bool is_core)
    {
        for(size_t i=0;i<a.n;++i)
        {
            Real J=a.F[i].determinant();if(!a.x[i].allFinite() || !a.v[i].allFinite() || !std::isfinite(J) || J<=.2)throw std::runtime_error("Invalid solid deformation state");
            minimum=std::min(minimum,J);kinetic+=.5*a.mass[i]*a.v[i].squaredNorm();
            if(is_core)ever_cap[i]=ever_cap[i] || dpc_surface[i]==2;
            else if(failure[i]!=0 || damage[i]!=0 || !pd_stress[i].allFinite())throw std::runtime_error("PD failure in an intact prefracture case");
        }
    };
    while(true)
    {
        Real delta=delta_at(stage_time);apply(delta);
        if(onset && step%200==0 && onset->scan())
        {
            // A genuine reciprocal bond-state change is written, but NO force,
            // constitutive or MLS update is evaluated on the broken topology.
            // This explicit boundary prevents a first-crack diagnostic from
            // masquerading as a qualified post-fracture capsule simulation.
            onset->commitFirstBreak();onset->writeBonds((out/"first_crack_bonds.csv").string());write_fields(stage);
            std::ofstream event(out/"first_crack.json");event<<std::setprecision(17)
                <<"{\"status\":\"FIRST_CRACK_DETECTED_STOPPED_BEFORE_POSTFRACTURE_UPDATE\",\"step\":"<<step
                <<",\"stage\":"<<stage<<",\"nominal_strain\":"<<delta/o.diameter<<",\"compression_m\":"<<delta
                <<",\"tensile_threshold_Pa\":"<<onset->threshold()<<",\"maximum_particle_tension_Pa\":"<<onset->maximum_particle_tension
                <<",\"maximum_triggered_pair_tension_Pa\":"<<onset->maximum_triggered_pair_tension
                <<",\"broken_undirected_bonds\":"<<onset->candidates.size()<<",\"maximum_damage\":"<<onset->maximum_damage
                <<",\"reaction_before_break_N\":"<<.5*(contact.upper+contact.lower)
                <<",\"postfracture_evolution\":false,\"postfracture_equilibrium\":false,\"fracture_energy_calibrated\":false}\n";
            std::cout<<"FIRST_CRACK strain="<<delta/o.diameter<<" broken_bonds="<<onset->candidates.size()<<"; diagnostic stopped before post-fracture evolution"<<std::endl;
            return 0;
        }
        Real kinetic=0,jc=1,js=1;inspect(c,kinetic,jc,true);inspect(s,kinetic,js,false);
        min_core_J=std::min(min_core_J,jc);min_pd_J=std::min(min_pd_J,js);max_gap=std::max(max_gap,tie.max_gap);
        Real dpc_dt=std::numeric_limits<Real>::infinity();
        for(size_t i=0;i<c.n;++i)dpc_dt=std::min(dpc_dt,o.cfl*1.15*o.spacing/(dpc_material.acousticBound(c.F[i].determinant())+c.v[i].norm()));
        Real quadratic_dt=o.quadratic_kinematics?quadratic_reference.elastic_dt_bound:std::numeric_limits<Real>::infinity();
        Real dt=std::min({o.max_dt,core_dt.exec(),pd_dt.exec(),constraint_dt,radial_dt,dpc_dt,quadratic_dt})*o.dt_scale;
        if(!std::isfinite(dt) || dt<1e-15 || step>=o.max_steps)throw std::runtime_error("Time step/count guard");
        if(step==0)write_fields(0);
        if(step%1000==0)
        {progress<<step<<','<<stage<<','<<stage_time<<','<<dt<<','<<delta<<','<<contact.upper<<','<<contact.lower<<','<<kinetic<<','<<tie.energy<<','<<contact.energy<<','<<jc<<','<<js<<'\n';progress.flush();}
        if(stage_time>=o.ramp_duration && stage_step%200==0)
        {
            Real squared=0;
            for(size_t i=0;i<c.n;++i)squared+=(c.force[i]+c.prior[i]).squaredNorm();
            // Internal PD force plus the CURRENT configuration's mapped loads.
            // Include native stabilization/damping held from the last midpoint.
            auto *hg=sp.getVariableDataByName<Vecd>("HourglassForce");auto *damp=sp.getVariableDataByName<Vecd>("PDNumericalDampingForce");
            for(size_t i=0;i<s.n;++i)squared+=(pd_internal[i]+hg[i]+damp[i]+s.mass[i]*pd_port.external[i]).squaredNorm();
            Real mean=.5*(contact.upper+contact.lower),residual=std::sqrt(squared)/std::max(mean,Real(1e-20));
            Real imbalance=std::abs(contact.upper-contact.lower)/std::max(mean,Real(1e-20));
            Real kr=kinetic/std::max(mean*delta,Real(1e-30));
            bool good=mean>1e-10 && residual<=o.force_tolerance && imbalance<=o.balance_tolerance && kr<=o.kinetic_tolerance;
            consecutive=good?consecutive+1:0;
            if(stage_step%10000==0)std::cout<<"load="<<stage<<" relax="<<stage_step<<" residual="<<residual<<" imbalance="<<imbalance<<" kinetic="<<kr<<std::endl;
            if(consecutive>=3)
            {
                eq<<stage<<','<<stage_step<<','<<delta<<','<<delta/o.diameter<<','<<contact.upper<<','<<contact.lower<<','<<residual<<','<<imbalance<<','<<kr<<','<<tie.max_gap<<','<<jc<<','<<js<<','<<tie.force_residual<<','<<tie.torque_residual<<','<<tie.power_residual<<'\n';eq.flush();
                Real V=0,V0=0,W=0,yr=0;size_t cap=0;
                for(size_t i=0;i<c.n;++i){V+=c.measure[i]*c.F[i].determinant();V0+=c.measure[i];W+=c.measure[i]*dpc_work[i];yr=std::max(yr,dpc_residual[i]);cap+=ever_cap[i];}
                history<<stage<<','<<delta/o.diameter<<','<<V/V0<<','<<*std::min_element(dpc_r,dpc_r+c.n)<<','<<*std::max_element(dpc_r,dpc_r+c.n)<<','
                    <<*std::min_element(dpc_kappa,dpc_kappa+c.n)<<','<<*std::max_element(dpc_kappa,dpc_kappa+c.n)<<','<<Real(cap)/c.n<<','<<W<<','<<yr<<'\n';history.flush();
                write_fields(stage);
                std::cout<<"ACCEPTED load="<<stage<<'/'<<o.load_steps<<" strain="<<delta/o.diameter<<" force="<<mean<<" residual="<<residual<<" imbalance="<<imbalance<<std::endl;
                if(stage==o.load_steps)break;
                ++stage;stage_time=0;stage_step=0;consecutive=0;continue;
            }
        }
        if(stage_time>o.max_stage_duration)throw std::runtime_error("Hybrid load step did not reach equilibrium;no complete result written");
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){core_first.initialization(i,dt);});
        pd_first.exec(dt);
        apply(delta_at(stage_time+.5*dt));
        // PD first-half copied old prior acceleration before drifting. Replace
        // it here with the common-midpoint interface/contact loads exactly once.
        for(size_t i=0;i<s.n;++i)pd_acc[i]=pd_port.external[i];
        pd_second(dt);pd_third(dt);pd_stabilize.exec(dt);pd_pair_damping.exec(dt);
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){core_first.interaction(i,dt);});
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){core_first.update(i,dt);});
        Real factor=std::exp(-o.damping_rate*dt);
        for(size_t i=0;i<c.n;++i){Real before=.5*c.mass[i]*c.v[i].squaredNorm();c.v[i]*=factor;core_relax_loss+=before*(1-factor*factor);}
        work+=.5*(contact.upper+contact.lower)*(delta_at(stage_time+dt)-delta);
        core_second.exec(dt);pd_fourth.exec(dt);
        time+=dt;stage_time+=dt;++step;++stage_step;
    }
    std::ofstream summary(out/"summary.json");summary<<std::setprecision(17)
        <<"{\n\"status\":\"COMPLETED_CONFIGURED_PD_SOLID_QUASISTATIC_CASE\",\n\"equilibrated_load_steps\":"<<stage
        <<",\n\"steps\":"<<step<<",\n\"core_particles\":"<<c.n<<",\n\"PD_shell_particles\":"<<s.n
        <<",\n\"core_min_J\":"<<min_core_J<<",\n\"shell_min_J\":"<<min_pd_J<<",\n\"largest_interface_gap_m\":"<<max_gap
        <<",\n\"max_interface_force_mapping_residual\":"<<tie.worst_force<<",\n\"max_interface_torque_mapping_residual\":"<<tie.worst_torque<<",\n\"max_interface_power_mapping_residual\":"<<tie.worst_power
        <<",\n\"max_contact_force_mapping_residual\":"<<contact.worst_force<<",\n\"max_contact_torque_mapping_residual\":"<<contact.worst_torque<<",\n\"max_contact_power_mapping_residual\":"<<contact.worst_power
        <<",\n\"plate_work_quadrature_J\":"<<work<<",\n\"core_global_relaxation_loss_J\":"<<core_relax_loss
        <<",\n\"material_calibrated\":false,\n\"thin_layer_spatial_convergence\":false,\n\"full_stored_energy_audit\":false\n}\n";
    return 0;
}

int main(int argc,char **argv)
{
    try{HybridOptions options;options.parseHybrid(argc,argv);return runHybrid(options);}
    catch(const std::exception &e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
}

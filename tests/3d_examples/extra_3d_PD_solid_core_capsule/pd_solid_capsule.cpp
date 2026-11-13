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
#include "pd_fragment_contact.h"
#include "pd_energy_crack_band.h"
#include "pd_layered_correspondence.h"
#include "pd_layered_surface.h"
#include "solid_residual_integration.h"
#include <sstream>
#include <map>
#include <tbb/task_arena.h>

struct HybridOptions : PDSolidCapsuleOptions
{
    std::string solid_stabilization="none";
    bool preflight_only=false;
    size_t postfracture_checkpoint_steps=0;
    bool quadratic_kinematics=false,layered_kinematics=false;
    unsigned plate_surface_order=1;
    size_t shell_layers=3;
    unsigned interface_order=1,interface_core_order=0; // zero inherits the common order
    Real interface_core_support=2.5,interface_shell_support=2.5;
    std::string shell_material="elastic",fracture_model="none",fracture_scope="none";
    std::string fragment_policy="error",core_contact_mesh;
    Real shell_yield=0,shell_hardening=0,fracture_stress=0,fracture_energy=0;
    Real pd_spacing=0,pd_cfl=.1,pd_hourglass=.5;
    Real ramp_duration=2.8132387884274494e-6,max_stage_duration=1.4066193942137247e-5;
    Real damping_rate=.002/2.8132387884274494e-10;
    Real plate_friction=0,plate_tangent_factor=1;
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
            if(key=="pd-kinematics"){if(value!="linear"&&value!="quadratic"&&value!="layered_quadratic")throw std::invalid_argument("pd-kinematics must be linear, quadratic or layered_quadratic");quadratic_kinematics=value!="linear";layered_kinematics=value=="layered_quadratic";return;}
            if(key=="plate-surface-order"){if(value!="1"&&value!="2")throw std::invalid_argument("plate-surface-order must be 1 or 2");plate_surface_order=unsigned(value[0]-'0');return;}
            if(key=="interface-order"){if(value!="1" && value!="2")throw std::invalid_argument("interface-order must be 1 or 2");interface_order=unsigned(value[0]-'0');return;}
            if(key=="interface-core-order")
            {if(value=="inherit")interface_core_order=0;else if(value=="1"||value=="2")interface_core_order=unsigned(value[0]-'0');
             else throw std::invalid_argument("interface-core-order must be inherit, 1 or 2");return;}
            if(key=="solid-stabilization"){if(value!="none"&&value!="quadratic_residual")throw std::invalid_argument("solid-stabilization must be none or quadratic_residual");solid_stabilization=value;return;}
            if(key=="core-backend"){if(value!="solid_sph")throw std::invalid_argument("Core requires solid_sph,not a WCSPH fluid");return;}
            if(key=="shell-backend"){if(value!="nosb_pd")throw std::invalid_argument("Shell requires nosb_pd");return;}
            if(key=="shell-material"){shell_material=value;return;}
            if(key=="fracture-model"){fracture_model=value;return;}
            if(key=="fracture-scope"){fracture_scope=value;return;}
            if(key=="fragment-policy"){fragment_policy=value;return;}
            if(key=="core-contact-mesh"){core_contact_mesh=value;return;}
            if(key=="postfracture-checkpoint-steps")
            {size_t used=0;postfracture_checkpoint_steps=std::stoull(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid checkpoint count");return;}
            if(key=="shell-layers"){size_t used=0;shell_layers=std::stoul(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid layer count");return;}
            Real *target=nullptr;
            if(key=="shell-yield-stress")target=&shell_yield;
            if(key=="shell-hardening-modulus")target=&shell_hardening;
            if(key=="fracture-stress")target=&fracture_stress;
            if(key=="fracture-energy")target=&fracture_energy;
            if(key=="interface-core-support-ratio")target=&interface_core_support;
            if(key=="interface-shell-support-ratio")target=&interface_shell_support;
            if(key=="pd-spacing")target=&pd_spacing;
            if(key=="pd-cfl")target=&pd_cfl;
            if(key=="pd-hourglass")target=&pd_hourglass;
            if(key=="ramp-duration")target=&ramp_duration;
            if(key=="max-stage-duration")target=&max_stage_duration;
            if(key=="damping-rate")target=&damping_rate;
            if(key=="plate-friction")target=&plate_friction;
            if(key=="plate-tangent-factor")target=&plate_tangent_factor;
            if(target){size_t used=0;*target=std::stod(value,&used);if(used!=value.size() || !std::isfinite(*target))throw std::invalid_argument("Invalid hybrid number");return;}
            base.push_back(arg);
        };
        for(int i=1;i<argc;++i)consume(argv[i]);
        std::vector<char*> ptr;for(auto &a:base)ptr.push_back(a.data());
        if(fracture_model=="crack_band_energy"&&fracture_scope=="growth")maximum_supported_strain=.5;
        else if(fracture_model=="principal_stress"&&fracture_scope=="growth")maximum_supported_strain=.25;
        PDSolidCapsuleOptions::parse(int(ptr.size()),ptr.data());
        if(pd_spacing==0)pd_spacing=spacing;
        if(interface_core_order==0)interface_core_order=interface_order;
        if(shell_material!="elastic"&&shell_material!="j2"&&shell_material!="corotated_biot"&&shell_material!="finite_volume_biot")throw std::invalid_argument("Unknown shell material");
        if(shell_yield<0||shell_hardening<0||(shell_material=="j2"&&shell_yield<=0)||
           (shell_material!="j2"&&(shell_yield!=0||shell_hardening!=0)))
            throw std::invalid_argument("J2 requires an explicit positive shell-yield-stress; elastic forbids plastic parameters");
        if(fracture_model=="none")
        {if(fracture_scope!="none"||fracture_stress!=0||fracture_energy!=0)throw std::invalid_argument("No-fracture mode forbids active fracture parameters");}
        else if(fracture_model=="crack_band_energy")
        {if(fracture_scope!="growth"||(shell_material!="corotated_biot"&&shell_material!="finite_volume_biot")||fracture_stress<=0||fracture_energy<=0||!quadratic_kinematics||interface_order!=2)
            throw std::invalid_argument("crack_band_energy requires corotated_biot, positive strength/Gc, quadratic PD and interface, growth scope");}
        else if(fracture_model!="principal_stress"||(fracture_scope!="first_crack"&&fracture_scope!="growth")||fracture_stress<=0||fracture_energy!=0)
            throw std::invalid_argument("principal_stress requires explicit threshold and first_crack or growth scope");
        if((shell_material=="corotated_biot"||shell_material=="finite_volume_biot")&&fracture_model!="crack_band_energy")throw std::invalid_argument("Biot energy material currently requires its verified crack-band assembly");
        if(fragment_policy!="error")throw std::invalid_argument("Only fragment-policy=error is allowed: numerical erosion is not a physical fracture parameter");
        if(fracture_scope=="growth"&&fracture_model=="principal_stress"&&(!quadratic_kinematics||interface_order!=2||core_contact_mesh.empty()))
            throw std::invalid_argument("Fracture growth requires quadratic PD, quadratic interface and a core-contact triangulation");
        if(fracture_model=="principal_stress"&&fracture_scope=="growth"&&interface_core_order!=interface_order)
            throw std::invalid_argument("Mixed core/PD reconstruction orders are not qualified for sharp-fracture growth");
        if(postfracture_checkpoint_steps&&(fracture_scope!="growth"||preflight_only||postfracture_checkpoint_steps>max_steps))throw std::invalid_argument("Invalid fracture-growth diagnostic checkpoint");
        if(postfracture_checkpoint_steps&&fracture_model=="crack_band_energy")throw std::invalid_argument("Use energy history outputs for continuous damage diagnostics");
        if(layered_kinematics&&(interface_order!=2||plate_surface_order!=2||fracture_model!="crack_band_energy"))throw std::invalid_argument("The qualified layered volume path requires quadratic inner/outer traces and energy crack-band damage");
        if(!coated || !quasistatic || material!="dpc" || plate_contact!="outer_surface" ||
            patch_test || rest_test || smoke_test || (shell_layers!=3 && shell_layers!=5) ||
            pd_spacing<=0 || pd_cfl<=0 || pd_cfl>.2 || pd_hourglass<0 ||
            ramp_duration<=0 || max_stage_duration<=ramp_duration || damping_rate<0 ||
            interface_core_support<=1 || interface_shell_support<=1)
            throw std::invalid_argument("Invalid PD-solid configuration or unsupported dynamic/fracture path");
        if(plate_friction<0||plate_tangent_factor<=0)
            throw std::invalid_argument("Coulomb friction must be nonnegative and tangent factor positive");
        if(plate_friction>0&&fracture_model=="principal_stress"&&fracture_scope=="growth")
            throw std::invalid_argument("Platen friction is not qualified for sharp-fragment contact; use the energy-damage path");
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
    UniquePtr<solid_dynamics::SolidQuadraticReference> core_residual_reference;
    UniquePtr<solid_dynamics::SolidResidualIntegration1stHalf> core_residual_first;
    if(o.solid_stabilization=="quadratic_residual")
    {core_residual_reference=makeUnique<solid_dynamics::SolidQuadraticReference>(core_inner,o.core_E/(2*(1+o.core_nu)),o.core_E/(1-2*o.core_nu));core_residual_first=makeUnique<solid_dynamics::SolidResidualIntegration1stHalf>(core_inner,*core_residual_reference);
     auto&p=core.getBaseParticles();p.addVariableToWrite<Real>("SolidResidualEnergyDensity");p.addVariableToWrite<Vecd>("SolidResidualForce");}
    Dynamics1Level<solid_dynamics::Integration2ndHalf> core_second(core_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> core_dt(core,o.cfl);

    PDBody shell(system,makeShared<DefaultShape>("PDShell"));
    shell.defineAdaptation<PDAdaptation>(1.5075,o.spacing/o.pd_spacing);
    if(o.shell_material!="j2")shell.defineMatterMaterial<HughesWingetSolid>(o.shell_rho,o.shell_E,o.shell_nu);
    else shell.defineMatterMaterial<J2PlasticityforPD>(o.shell_rho,o.shell_E,o.shell_nu,o.shell_yield,
        PDIsotropicHardening::Linear(o.shell_hardening),0);
    shell.generateParticles<BaseParticles,LayeredPDCapsuleGenerator>(rc,o.thickness,o.pd_spacing,o.shell_layers);
    InnerRelation shell_inner(shell);
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> pd_initialize(shell,makeShared<Gravity>(Vecd::Zero()));
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> pd_shape(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> pd_first(shell);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> pd_linear_second(shell_inner);
    bool energy_growth=o.fracture_model=="crack_band_energy";
    if(energy_growth&&o.shell_material!="finite_volume_biot")throw std::invalid_argument("Use shell-material=finite_volume_biot for the actual-volume corrected crack-band model; legacy Biot remains in its archived executable");
    bool growth=o.fracture_scope=="growth"&&!energy_growth;
    UniquePtr<LayeredParticleGeometry> layered_geometry;
    if(o.layered_kinematics)
    {
        auto&p=shell.getBaseParticles();size_t columns_count=p.TotalRealParticles()/o.shell_layers;auto*x=p.getVariableDataByName<Vecd>("Position");
        StdVec<StdVec<size_t>> columns(columns_count);StdVec<Vecd>directors(columns_count);
        for(size_t c=0;c<columns_count;++c){for(size_t l=0;l<o.shell_layers;++l)columns[c].push_back(l*columns_count+c);directors[c]=x[c].normalized();}
        ThicknessQuadrature normal_rule(o.shell_layers);layered_geometry=makeUnique<LayeredParticleGeometry>(p,columns,directors,normal_rule.xi);
    }
    UniquePtr<solid_dynamics::PDQuadraticReference> reference_owner;
    solid_dynamics::PDFractureReference *fracture_reference=nullptr;
    if(growth){auto r=makeUnique<solid_dynamics::PDFractureReference>(shell_inner);fracture_reference=r.get();reference_owner=std::move(r);}
    else if(layered_geometry)reference_owner=makeUnique<solid_dynamics::PDLayeredQuadraticReference>(shell_inner,*layered_geometry);
    else reference_owner=makeUnique<solid_dynamics::PDQuadraticReference>(shell_inner);
    auto &quadratic_reference=*reference_owner;
    InteractionWithUpdate<solid_dynamics::NosbPDQuadraticSecondStep> pd_quadratic_second(shell_inner,quadratic_reference);
    UniquePtr<BaseDynamics<void>> plastic_second;
    if(o.shell_material=="j2")
    {
        if(o.quadratic_kinematics)plastic_second=makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDQuadraticSecondStepPlastic>>(shell_inner,quadratic_reference);
        else plastic_second=makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic>>(shell_inner);
    }
    UniquePtr<solid_dynamics::PDEnergyCrackBand> energy_model;
    auto pd_second=[&](Real dt){if(energy_model)particle_for(execution::ParallelPolicy(),shell.LoopRange(),[&](size_t i){energy_model->evaluate(i,dt>0);});
        else if(plastic_second)plastic_second->exec(dt);else if(o.quadratic_kinematics)pd_quadratic_second.exec(dt);else pd_linear_second.exec(dt);};
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> pd_linear_third(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDQuadraticThirdStep> pd_quadratic_third(shell_inner,quadratic_reference);
    auto pd_third=[&](Real dt){if(energy_model)particle_for(execution::ParallelPolicy(),shell.LoopRange(),[&](size_t i){energy_model->forces(i);});
        else if(o.quadratic_kinematics)pd_quadratic_third.exec(dt);else pd_linear_third.exec(dt);};
    InteractionDynamics<solid_dynamics::ProjectedPairStabilization> pd_stabilize(shell_inner,o.pd_hourglass,.1);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> pd_pair_damping(shell_inner);
    SimpleDynamics<RelaxedPDFourthStep> pd_fourth(shell,o.damping_rate);
    ReduceDynamics<solid_dynamics::PDTimeStep> pd_dt(shell,o.pd_cfl);
    if(energy_growth)energy_model=makeUnique<solid_dynamics::PDEnergyCrackBand>(shell_inner,quadratic_reference,o.shell_E,o.shell_nu,o.fracture_stress,o.fracture_energy,o.pd_hourglass,
        shell.getSPHAdaptation().getKernel()->CutOffRadius());
    auto pd_auxiliary=[&](Real dt){if(!energy_model){pd_stabilize.exec(dt);pd_pair_damping.exec(dt);}};
    State c(core),s(shell);
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();core_correction.exec();if(core_residual_reference)core_residual_reference->build();pd_shape.exec();
    if(growth)fracture_reference->initialize();else if(o.quadratic_kinematics)quadratic_reference.build();
    if(energy_model)energy_model->initialize();
    pd_second(0);
    SolidForcePort core_port(core,SolidForcePort::Storage::ForcePrior);
    SolidForcePort pd_port(shell,SolidForcePort::Storage::PDAccelerationPrior);
    auto interface_directions=sphereDirections(R,o.pd_spacing);
    auto contact_directions=sphereDirections(R,o.contact_quadrature_scale*o.spacing);
    StdVec<Vecd> inner_points,outer_points;
    for(const auto &n:interface_directions)inner_points.push_back(rc*n);
    for(const auto &n:contact_directions)outer_points.push_back(R*n);
    ReferenceSolidSurface core_surface(core_port,inner_points,o.spacing,o.interface_core_support,o.interface_core_order);
    UniquePtr<ReferenceSolidSurface> inner_owner,outer_owner;
    FractureAwareSurface *fracture_inner=nullptr,*fracture_outer=nullptr;
    if(growth)
    {auto a=makeUnique<FractureAwareSurface>(pd_port,inner_points,o.pd_spacing,o.interface_shell_support,o.interface_order,*fracture_reference);fracture_inner=a.get();inner_owner=std::move(a);
     auto b=makeUnique<FractureAwareSurface>(pd_port,outer_points,o.pd_spacing,2.5,1,*fracture_reference);fracture_outer=b.get();outer_owner=std::move(b);}
    else if(layered_geometry)
    {inner_owner=makeUnique<LayeredSolidSurface>(pd_port,inner_points,*layered_geometry,-1,o.pd_spacing,o.interface_shell_support);
     outer_owner=makeUnique<LayeredSolidSurface>(pd_port,outer_points,*layered_geometry,1,o.pd_spacing,2.5);}
    else {inner_owner=makeUnique<ReferenceSolidSurface>(pd_port,inner_points,o.pd_spacing,o.interface_shell_support,o.interface_order);outer_owner=makeUnique<ReferenceSolidSurface>(pd_port,outer_points,o.pd_spacing,2.5,o.plate_surface_order);}
    auto &pd_inner=*inner_owner;auto &pd_outer=*outer_owner;
    Real Eref=o.constraint_E>0?o.constraint_E:o.core_E;
    StdVec<Real> inner_area(inner_points.size(),4*Pi*rc*rc/inner_points.size());
    StdVec<Real> outer_area(outer_points.size(),4*Pi*R*R/outer_points.size());
    BondedSolidInterface tie(pd_inner,core_surface,inner_area,o.tie_factor*Eref/o.spacing,R);
    SolidSurfacePlaneContact contact(pd_outer,outer_area,o.contact_factor*Eref/o.spacing,R,o.plate_friction,o.plate_tangent_factor);
    Real constraint_dt=1/std::sqrt(1/(tie.stable_dt*tie.stable_dt)+1/(contact.stable_dt*contact.stable_dt));
    ThicknessQuadrature quadrature(o.shell_layers);
    Real radial_spacing=o.thickness;
    for(size_t i=1;i<quadrature.xi.size();++i)radial_spacing=std::min(radial_spacing,.5*o.thickness*(quadrature.xi[i]-quadrature.xi[i-1]));
    auto &pd_material=DynamicCast<HughesWingetSolid>(&shell,shell.getMatterMaterial());
    Real radial_dt=o.pd_cfl*std::min(radial_spacing,o.pd_spacing)/pd_material.ReferenceSoundSpeed();
    UniquePtr<ReferenceSolidSurface> exposed_surface;UniquePtr<ExposedCorePlaneContact> exposed_contact;UniquePtr<PDFragmentContact> fragment_contact;
    if(growth)
    {
        StdVec<Vecd> q;for(const auto &n:contact_directions)q.push_back(rc*n);
        exposed_surface=makeUnique<ReferenceSolidSurface>(core_port,q,o.spacing);
        StdVec<Real> area(q.size(),4*Pi*rc*rc/q.size());exposed_contact=makeUnique<ExposedCorePlaneContact>(*exposed_surface,pd_outer,area,o.contact_factor*Eref/o.spacing);
        StdVec<Matd> metric(s.n);
        for(size_t i=0;i<s.n;++i){Real r=s.x0[i].norm();Vecd n=s.x0[i]/r;
            Real ar=std::min(.5*radial_spacing,.75*std::min(r-rc,R-r)),at=.5*o.pd_spacing;
            metric[i]=at*at*Matd::Identity()+(ar*ar-at*at)*n*n.transpose();}
        fragment_contact=makeUnique<PDFragmentContact>(pd_port,*fracture_reference,core_surface,tie.x_b,metric,o.contact_factor*Eref/o.spacing,o.core_contact_mesh);
    }
    auto &cp=core.getBaseParticles();auto &sp=shell.getBaseParticles();
    std::ofstream residual_history;
    if(core_residual_reference){residual_history.open(out/"solid_stabilization.csv");residual_history<<std::setprecision(17)<<"load_step,strain,residual_energy_J,core_elastic_energy_J,core_plastic_work_J,residual_fraction_of_material_energy\n";}
    auto *pd_acc=sp.getVariableDataByName<Vecd>("PDAcceleration");
    auto *pd_internal=sp.getVariableDataByName<Vecd>("PDInternalForce");
    auto *pd_K=sp.getVariableDataByName<Matd>("ShapeTensorK");
    auto *pd_stress=sp.getVariableDataByName<Matd>("PDStress");
    auto *failure=sp.getVariableDataByName<int>("NumericalFailureReason");
    auto *damage=sp.getVariableDataByName<Real>("Damage");
    UniquePtr<solid_dynamics::PrincipalStressBondOnset> onset;
    if(o.fracture_model=="principal_stress")onset=makeUnique<solid_dynamics::PrincipalStressBondOnset>(shell_inner,o.fracture_stress,growth);
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
        <<"shell_constitutive="<<(energy_growth?"finite_volume_Biot_energy_damage":o.shell_material=="elastic"?"HughesWinget_elastic":"J2_HW_plastic")
        <<"\nfracture_model="<<o.fracture_model<<"\nfracture_scope="<<o.fracture_scope<<"\nfracture_stress="<<o.fracture_stress
        <<"\nfragment_policy="<<o.fragment_policy<<"\ncore_contact_mesh="<<o.core_contact_mesh<<"\nfracture_energy_regularized="<<energy_growth
        <<"\nfracture_energy_J_m2="<<o.fracture_energy<<"\nfracture_representation="<<(energy_growth?"diffuse_energy_damage":"sharp_bond_failure")
        <<"\npostfracture_checkpoint_steps="<<o.postfracture_checkpoint_steps
        <<"\nshell_yield_stress="<<o.shell_yield<<"\nshell_hardening_modulus="<<o.shell_hardening<<"\ncore_material=DPC\n"
        <<"shell_E="<<o.shell_E<<"\nshell_nu="<<o.shell_nu<<"\nshell_reference_density="<<o.shell_rho
        <<"\nsolid_stabilization="<<o.solid_stabilization<<"\nsolid_residual_reference_dt_bound="<<(core_residual_reference?core_residual_reference->native_dt_bound:0)
        <<"\ninterface_law=perfect_bond_vector_penalty\n"
        <<"core_particles="<<c.n<<"\nshell_particles="<<s.n<<"\nshell_layers="<<o.shell_layers
        <<"\ncore_diameter="<<2*rc<<"\nouter_diameter="<<2*R<<"\nshell_thickness="<<o.thickness
        <<"\ncore_spacing="<<o.spacing<<"\npd_tangential_spacing="<<o.pd_spacing<<"\nminimum_PD_radial_spacing="<<radial_spacing
        <<"\nPD_horizon="<<shell.getSPHAdaptation().getKernel()->CutOffRadius()
        <<"\ncore_mass="<<c.totalMass()<<"\nshell_mass="<<s.totalMass()<<"\nexpected_shell_mass="<<expected_shell_mass
        <<"\nPD_shape_condition_max="<<shape_condition<<"\ninitial_PD_F_error="<<initial_F_error
        <<"\ninterface_points="<<inner_points.size()<<"\ncontact_points="<<outer_points.size()
        <<"\ncore_surface_affine_error="<<core_surface.affine_error<<"\nPD_inner_affine_error="<<pd_inner.affine_error<<"\nPD_outer_affine_error="<<pd_outer.affine_error
        <<"\nplate_surface_order="<<o.plate_surface_order<<"\nshell_surface_reconstruction="<<(layered_geometry?"normal_Lagrange_then_quadratic_surface_MLS":"volume_MLS")<<"\ninterface_order="<<o.interface_order<<"\ninterface_core_support_ratio="<<o.interface_core_support<<"\ninterface_shell_support_ratio="<<o.interface_shell_support
        <<"\ninterface_core_order="<<o.interface_core_order<<"\ninterface_shell_order="<<(layered_geometry?2:o.interface_order)
        <<"\ncore_surface_polynomial_error_applicable="<<(o.interface_core_order==2)
        <<"\ncore_surface_polynomial_error="<<core_surface.polynomial_error<<"\nPD_inner_polynomial_error="<<pd_inner.polynomial_error
        <<"\ninterface_K="<<o.tie_factor*Eref/o.spacing<<"\ncontact_K="<<o.contact_factor*Eref/o.spacing
        <<"\nplate_contact_law="<<(o.plate_friction>0?"normal_penalty_coulomb":"normal_penalty_frictionless")
        <<"\nplate_friction_mu="<<o.plate_friction<<"\nplate_tangent_factor="<<o.plate_tangent_factor
        <<"\nfriction_history_time_layer=midpoint_commit\ninterface_dt_bound="<<tie.stable_dt<<"\ncontact_dt_bound="<<contact.stable_dt
        <<"\ncombined_constraint_dt="<<constraint_dt<<"\nPD_radial_acoustic_dt="<<radial_dt
        <<"\nmaximum_dt_cap="<<o.max_dt<<"\nnumerical_ramp_duration="<<o.ramp_duration<<"\nnumerical_damping_rate="<<o.damping_rate
        <<"\nload_steps="<<o.load_steps<<"\nstrain="<<o.strain<<"\nthreads="<<o.threads
        <<"\npd_kinematics="<<(o.layered_kinematics?"layered_quadratic_reproducing":o.quadratic_kinematics?"quadratic_reproducing":"linear_correspondence")
        <<"\nquadratic_moment_condition_max="<<quadratic_reference.max_moment_condition<<"\nquadratic_reproduction_error="<<quadratic_reference.max_reproduction_error
        <<"\nquadratic_initial_elastic_dt_bound="<<quadratic_reference.elastic_dt_bound
        <<"\ndpc_shear_hardening="<<dpc_material.ShearHardeningModulus()
        <<"\nsource_sha256="<<CAPSULE_SOURCE_SHA256<<"\ndpc_material_sha256="<<DPC_MATERIAL_SHA256<<'\n';
    if(energy_model)config<<"crack_band_width_m="<<energy_model->fracture.width<<"\ncrack_band_kappa0="<<energy_model->fracture.k0<<"\ncrack_band_kappaf="<<energy_model->fracture.kf
        <<"\nenergy_stabilizer_dt="<<energy_model->stabilizer_dt<<"\nenergy_stabilizer=quadratic_LS_residual\nlegacy_pair_damping=disabled_in_energy_branch\nstrength_definition=uniaxial_nominal_exact_lateral_equilibrium\nvolumetric_energy=K/4*(J^2-1-2*log(J))\nvolumetric_split=actual_J\nfinite_strain_tangent_bound=enabled\n";
    config.flush();
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
        Real friction_force_sum=0,friction_cone_ratio=0,friction_stored=0;
        if(o.plate_friction>0)
        {
            // The preceding conservative FD audit uses uncommitted, zero
            // tangential histories. Sliding is checked separately, not falsely
            // treated as the gradient of a conservative scalar potential.
            core_port.clear();pd_port.clear();contact.apply(-R+2e-6,R-2e-6,true);
            for(size_t i=0;i<c.n;++i)c.x[i]+=Vecd(2e-8,-1e-8,0);
            for(size_t i=0;i<s.n;++i)s.x[i]+=Vecd(2e-8,-1e-8,0);
            auto commits=contact.history_commits;Real prior_D=contact.slip_dissipation;
            core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+2e-6,R-2e-6);
            friction_force_sum=contact.tangential_force_sum;friction_cone_ratio=contact.max_cone_ratio;
            friction_stored=contact.tangent_energy;
            core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+2e-6,R-2e-6);
            if(friction_force_sum<=0||friction_stored<=0||friction_cone_ratio>1+1e-12||
               contact.history_commits!=commits||contact.slip_dissipation!=prior_D||
               std::max({contact.worst_force,contact.worst_torque,contact.worst_power})>1e-9)
                throw std::runtime_error("Nonzero Coulomb force-map or const-trial check failed");
            core_port.clear();pd_port.clear();contact.apply(-R+2e-6,R-2e-6,true);
            if(contact.slip_dissipation<0||contact.release_energy<0)
                throw std::runtime_error("Negative Coulomb dissipation");
        }
        std::ofstream f(out/"preflight.json");f<<std::setprecision(17)
            <<"{\"status\":\"INITIAL_GEOMETRY_AFFINE_AND_LOADED_FORCE_MAPS_PASSED\",\"core_particles\":"<<c.n<<",\"PD_particles\":"<<s.n
            <<",\"PD_shape_condition_max\":"<<shape_condition<<",\"initial_F_error\":"<<initial_F_error
            <<",\"manufactured_PD_affine_error\":"<<affine_PD_error<<",\"loaded_energy_gradient_relative_error\":"<<finite_difference_error
            <<",\"interface_force_residual\":"<<tie.worst_force<<",\"interface_torque_residual\":"<<tie.worst_torque<<",\"interface_power_residual\":"<<tie.worst_power
            <<",\"contact_force_residual\":"<<contact.worst_force<<",\"contact_torque_residual\":"<<contact.worst_torque<<",\"contact_power_residual\":"<<contact.worst_power
            <<",\"friction_mu\":"<<o.plate_friction<<",\"nonzero_friction_force_sum_N\":"<<friction_force_sum
            <<",\"friction_cone_ratio\":"<<friction_cone_ratio<<",\"friction_stored_J\":"<<friction_stored
            <<",\"loaded_energy_FD_includes_friction_history\":false,\"friction_trial_is_const\":true"
            <<",\"full_bending_validation\":false,\"physical_history_points\":0}\n";
        std::cout<<"PREFLIGHT geometry and affine reference maps passed;no compression history evolved."<<std::endl;return 0;
    }
    for(auto *body:{static_cast<RealBody*>(&core),static_cast<RealBody*>(&shell)})
    {body->getBaseParticles().addVariableToWrite<Vecd>("Velocity");body->getBaseParticles().addVariableToWrite<Matd>("DeformationGradient");}
    sp.addVariableToWrite<Matd>("PDStress");
    sp.addVariableToWrite<Real>("Damage");sp.addVariableToWrite<int>("BondNumber");sp.addVariableToWrite<int>("ParticleLive");
    if(growth){sp.addVariableToWrite<int>("PDReproductionOrder");sp.addVariableToWrite<int>("ResolutionEroded");}
    if(energy_model){for(const auto&name:{"FractureDissipationDensity","ShellStoredEnergyDensity","ShellStabilizationEnergyDensity","FractureHistoryKappa","ShellTangentModulusBound"})sp.addVariableToWrite<Real>(name);
        sp.addVariableToWrite<int>("MaterialTensileFailure");sp.addVariableToWrite<Vecd>("CrackNormalReference");}
    cp.addVariableToWrite<Matd>("DPCCauchyStress");
    // Keep material histories available for auditing coupled load paths.
    for(const auto &name:{"DPCRelativeDensity","DPCPlasticCompaction","DPCAccumulatedPlasticShear","DPCEffectiveCohesion","DPCPlasticStressWork","DPCElasticEnergyDensity","DPCYieldResidual","DPCActiveSurface"})cp.addVariableToWrite<Real>(name);
    cp.addVariableToWrite<Matd>("DPCInversePlasticMetric");
    BodyStatesRecordingToVtp vtp(system);
    std::ofstream eq(out/"equilibria.csv"),history(out/"material_history.csv"),progress(out/"progress.csv");
    std::ofstream fracture_events(out/"fracture_events.csv"),fracture_history(out/"fracture_history.csv");
    std::ofstream energy_history,energy_events;
    const std::string energy_columns="step,stage,strain,initiated_particles,tensile_failed_particles,max_damage,stored_J,stabilization_J,fracture_dissipation_J,equivalent_crack_area_m2\n";
    if(energy_model){energy_history.open(out/"energy_history.csv");energy_events.open(out/"energy_events.csv");energy_history<<energy_columns;energy_events<<energy_columns;}
    fracture_events<<"event,step,stage,strain,new_physical_broken_bonds,total_physical_broken_bonds,resolution_removed_bonds,ballistic_particles,affine_gradient_rows,inactive_inner_points,inactive_outer_points,topology_version,gradient_dt_bound,tie_dt_bound,plate_dt_bound,fragment_dt_bound\n";
    fracture_history<<"stage,strain,total_physical_broken_bonds,maximum_damage,ballistic_particles,resolution_removed_bonds,affine_gradient_rows,inactive_inner_points,inactive_outer_points,self_contacts,core_contacts,debris_plate_contacts,fragment_force_residual,fragment_torque_residual,fragment_power_residual\n";
    eq<<std::setprecision(17)<<"load_step,relaxation_steps,compression_m,strain,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,max_tie_gap_m,min_core_J,min_PD_J,interface_force_residual,interface_torque_residual,interface_power_residual\n";
    history<<std::setprecision(17)<<"load_step,strain,core_volume_ratio,min_relative_density,max_relative_density,min_hardening_strain,max_hardening_strain,ever_cap_fraction,plastic_stress_work_J,max_yield_residual\n";
    progress<<std::setprecision(17)<<"step,stage,stage_time_s,dt_s,compression_m,upper_N,lower_N,kinetic_J,interface_J,contact_J,min_core_J,min_PD_J\n";
    Real &time=*system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t stage=1,stage_step=0,step=0,consecutive=0;Real stage_time=0,max_gap=0,min_core_J=1,min_pd_J=1,work=0,core_relax_loss=0;
    auto delta_at=[&](Real t){Real a=std::clamp(t/o.ramp_duration,Real(0),Real(1));
        return o.strain*o.diameter*(Real(stage-1)+.5*(1-std::cos(Pi*a)))/o.load_steps;};
    size_t fracture_event_count=0,first_fracture_step=0;Real maximum_fragment_force_error=0,maximum_fragment_torque_error=0,maximum_fragment_power_error=0;
    auto update_constraint_bound=[&]()
    {Real rate=1/(tie.stable_dt*tie.stable_dt)+1/(contact.stable_dt*contact.stable_dt);
     if(growth&&fracture_outer->inactive_points)rate+=1/(exposed_contact->stable_dt*exposed_contact->stable_dt);constraint_dt=1/std::sqrt(rate);};
    auto plate_speed=[&](Real t){Real a=std::clamp(t/o.ramp_duration,Real(0),Real(1));
        return t<o.ramp_duration?o.strain*o.diameter*.5*Pi*std::sin(Pi*a)/(o.load_steps*o.ramp_duration):0;};
    auto apply=[&](Real delta,Real speed=0,bool commit_friction=false)
    {
        core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+.5*delta,R-.5*delta,commit_friction);
        if(growth&&fracture_reference->topology_version)
        {exposed_contact->apply(-R+.5*delta,R-.5*delta);fragment_contact->apply(-R+.5*delta,R-.5*delta,.5*speed,-.5*speed);
         contact.upper+=exposed_contact->upper+fragment_contact->upper;contact.lower+=exposed_contact->lower+fragment_contact->lower;
         contact.energy+=exposed_contact->energy+fragment_contact->energy;
         maximum_fragment_force_error=std::max(maximum_fragment_force_error,fragment_contact->force_residual);
         maximum_fragment_torque_error=std::max(maximum_fragment_torque_error,fragment_contact->torque_residual);
         maximum_fragment_power_error=std::max(maximum_fragment_power_error,fragment_contact->power_residual);}
    };
    std::ofstream friction_output(out/"friction_history.csv");
    friction_output<<std::setprecision(17)<<"kind,step,stage,strain,normal_penalty_J,tangent_penalty_J,slip_dissipation_J,contact_release_J,active_points,sliding_points,max_cone_ratio,history_commits\n";
    auto write_friction=[&](const char *kind,Real delta)
    {friction_output<<kind<<','<<step<<','<<stage<<','<<delta/o.diameter<<','<<contact.energy-contact.tangent_energy<<','
        <<contact.tangent_energy<<','<<contact.slip_dissipation<<','<<contact.release_energy<<','<<contact.active_friction_points<<','
        <<contact.sliding_points<<','<<contact.max_cone_ratio<<','<<contact.history_commits<<'\n';friction_output.flush();};
    auto write_fields=[&](size_t frame)
    {
        vtp.writeToFile(frame);
        if(o.plate_friction>0)
        {
            fs::create_directories(out/"platen_surface");std::ofstream q(out/"platen_surface"/("platen_"+std::to_string(frame)+".csv"));
            q<<std::setprecision(17)<<"x_m,y_m,z_m,fx_N,fy_N,fz_N,area_m2,contact_side\n";
            for(size_t i=0;i<contact.position.size();++i)q<<contact.position[i][0]<<','<<contact.position[i][1]<<','<<contact.position[i][2]<<','
                <<contact.force[i][0]<<','<<contact.force[i][1]<<','<<contact.force[i][2]<<','<<outer_area[i]<<','<<contact.friction_history[i].side<<'\n';
        }
        fs::create_directories(out/"interface_surface");std::ofstream f(out/"interface_surface"/("interface_"+std::to_string(frame)+".csv"));
        f<<std::setprecision(17)<<"x0_m,y0_m,z0_m,PD_x_m,PD_y_m,PD_z_m,core_x_m,core_y_m,core_z_m,fx_N,fy_N,fz_N,area_m2\n";
        for(size_t i=0;i<inner_points.size();++i)f<<inner_points[i][0]<<','<<inner_points[i][1]<<','<<inner_points[i][2]<<','<<tie.x_a[i][0]<<','<<tie.x_a[i][1]<<','<<tie.x_a[i][2]<<','
            <<tie.x_b[i][0]<<','<<tie.x_b[i][1]<<','<<tie.x_b[i][2]<<','<<tie.force_a[i][0]<<','<<tie.force_a[i][1]<<','<<tie.force_a[i][2]<<','<<inner_area[i]<<'\n';
    };
    size_t last_energy_initiated=0,last_energy_failed=0,last_energy_vtp_failed=0,energy_event_count=0;
    auto write_energy=[&](std::ostream&file,Real delta)
    {auto a=energy_model->statistics();file<<std::setprecision(17)<<step<<','<<stage<<','<<delta/o.diameter<<','<<a.initiated<<','<<a.failed<<','<<a.max_damage<<','
        <<a.stored<<','<<a.stabilization<<','<<a.dissipated<<','<<a.dissipated/o.fracture_energy<<'\n';file.flush();return a;};
    std::ofstream stability_trace(out/"stability_trace.csv");
    stability_trace<<std::setprecision(17)<<"step,stage,strain,particle,J,damage,lambda_max,lambda_mid,lambda_min,trace_Biot,HG_J_m3,x0,y0,z0,xx,xy,xz,yx,yy,yz,zx,zy,zz\n";
    auto dump_energy_failure=[&](const std::string &reason)
    {
        std::ofstream f(out/"energy_update_failure.txt");f<<std::setprecision(17)<<"step="<<step<<" stage="<<stage<<" strain="<<delta_at(stage_time)/o.diameter<<" reason="<<reason<<'\n';
        for(size_t i=0;i<s.n;++i)if(!std::isfinite(s.F[i].determinant())||s.F[i].determinant()<=.2)
            f<<"shell_particle="<<i<<" J="<<s.F[i].determinant()<<" damage="<<damage[i]<<" F="<<s.F[i]<<'\n';
        write_fields(999998);
    };
    auto inspect=[&](State &a,Real &kinetic,Real &minimum,bool is_core)
    {
        for(size_t i=0;i<a.n;++i)
        {
            Real J=a.F[i].determinant();if(!a.x[i].allFinite() || !a.v[i].allFinite() || !std::isfinite(J) || J<=.2)
            {
                std::ofstream f(out/"invalid_state.json");f<<std::setprecision(17)
                    <<"{\"body\":\""<<(is_core?"core":"shell")<<"\",\"particle\":"<<i<<",\"step\":"<<step
                    <<",\"stage\":"<<stage<<",\"strain\":"<<delta_at(stage_time)/o.diameter<<",\"J\":"<<J;
                if(!is_core&&energy_model)f<<",\"damage\":"<<damage[i]<<",\"kappa\":"<<energy_model->history[i].kappa;
                f<<",\"F_row_major\":[";for(int k=0;k<9;++k){if(k)f<<',';f<<a.F[i](k/3,k%3);}f<<"]}\n";f.flush();
                write_fields(999999);throw std::runtime_error("Invalid solid deformation state; exact particle saved in invalid_state.json");
            }
            minimum=std::min(minimum,J);kinetic+=.5*a.mass[i]*a.v[i].squaredNorm();
            if(is_core)ever_cap[i]=ever_cap[i] || dpc_surface[i]==2;
            else if(failure[i]!=0 || (!growth&&!energy_growth&&damage[i]!=0) || !std::isfinite(damage[i])||damage[i]<0||damage[i]>1 || !pd_stress[i].allFinite())throw std::runtime_error("Invalid PD damage/material state");
        }
    };
    while(true)
    {
        Real delta=delta_at(stage_time);apply(delta,plate_speed(stage_time));
        if(onset && step%200==0 && onset->scan())
        {
            if(growth)
            {
                onset->commitFirstBreak();size_t newly=onset->candidates.size();fracture_reference->updateAfterFailure();
                onset->maximum_damage=std::max(onset->maximum_damage,*std::max_element(damage,damage+s.n));
                fracture_inner->rebuildChanged();fracture_outer->rebuildChanged();
                tie.stable_dt=ReferenceSolidSurface::stableBound(pd_inner,&core_surface,tie.stiffness);
                contact.refreshTimeStepBound();update_constraint_bound();
                pd_second(0);for(size_t i=0;i<s.n;++i)pd_acc[i].setZero();pd_third(0);pd_auxiliary(0);
                apply(delta,plate_speed(stage_time));for(size_t i=0;i<s.n;++i)pd_acc[i]+=pd_port.external[i];consecutive=0;
                if(fracture_event_count==0)first_fracture_step=step;++fracture_event_count;
                fracture_events<<std::setprecision(17)<<fracture_event_count<<','<<step<<','<<stage<<','<<delta/o.diameter<<','<<newly<<','<<onset->total_physical_broken_bonds<<','
                    <<fracture_reference->resolution_removed_bonds<<','<<fracture_reference->ballistic_particles<<','<<fracture_reference->affine_rows<<','
                    <<fracture_inner->inactive_points<<','<<fracture_outer->inactive_points<<','<<fracture_reference->topology_version<<','<<quadratic_reference.elastic_dt_bound<<','
                    <<tie.stable_dt<<','<<contact.stable_dt<<','<<fragment_contact->stable_dt<<'\n';fracture_events.flush();
                if(fracture_event_count==1 || fracture_event_count%10==0)write_fields(100000+fracture_event_count);
                std::cout<<"FRACTURE event="<<fracture_event_count<<" strain="<<delta/o.diameter<<" new_bonds="<<newly<<" ballistic="<<fracture_reference->ballistic_particles<<std::endl;
            }
            else
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
        }
        Real kinetic=0,jc=1,js=1;inspect(c,kinetic,jc,true);inspect(s,kinetic,js,false);
        if(energy_model&&step%200==0)
        {auto a=energy_model->statistics();if(a.initiated!=last_energy_initiated||a.failed!=last_energy_failed)
         {bool frame=(last_energy_initiated==0&&a.initiated>0)||(last_energy_failed==0&&a.failed>0)||a.failed>=last_energy_vtp_failed+std::max(size_t(1),s.n/100);
          write_energy(energy_events,delta);++energy_event_count;if(frame){write_fields(100000+energy_event_count);last_energy_vtp_failed=a.failed;}
          if(last_energy_initiated==0||last_energy_failed==0&&a.failed>0)std::cout<<"ENERGY_DAMAGE strain="<<delta/o.diameter<<" initiated="<<a.initiated<<" failed="<<a.failed<<" dissipation="<<a.dissipated<<std::endl;
          last_energy_initiated=a.initiated;last_energy_failed=a.failed;}}
        min_core_J=std::min(min_core_J,jc);min_pd_J=std::min(min_pd_J,js);max_gap=std::max(max_gap,tie.max_gap);
        if(energy_model && stage>=6 && step%200==0)
        {
            size_t j=0;for(size_t i=1;i<s.n;++i)if(s.F[i].determinant()<s.F[j].determinant())j=i;
            Eigen::JacobiSVD<Matd> svd(s.F[j]);auto l=svd.singularValues();
            stability_trace<<step<<','<<stage<<','<<delta/o.diameter<<','<<j<<','<<s.F[j].determinant()<<','<<damage[j]<<','<<l[0]<<','<<l[1]<<','<<l[2]<<','<<l.sum()-3<<','<<energy_model->hg_density[j];
            for(int k=0;k<3;++k)stability_trace<<','<<s.x0[j][k];for(int k=0;k<9;++k)stability_trace<<','<<s.F[j](k/3,k%3);stability_trace<<'\n';stability_trace.flush();
        }

        if(o.postfracture_checkpoint_steps&&fracture_event_count&&step-first_fracture_step>=o.postfracture_checkpoint_steps)
        {
            write_fields(200000);std::ofstream f(out/"continuation_checkpoint.json");f<<std::setprecision(17)
                <<"{\"status\":\"DIAGNOSTIC_POSTFRACTURE_CONTINUATION_NOT_EQUILIBRATED\",\"step\":"<<step<<",\"first_fracture_step\":"<<first_fracture_step
                <<",\"postfracture_steps\":"<<step-first_fracture_step<<",\"events\":"<<fracture_event_count<<",\"physical_broken_bonds\":"<<onset->total_physical_broken_bonds
                <<",\"resolution_removed_bonds\":"<<fracture_reference->resolution_removed_bonds<<",\"ballistic_particles\":"<<fracture_reference->ballistic_particles
                <<",\"mass_preserved\":"<<(std::abs(s.totalMass()/expected_shell_mass-1)<1e-10&&std::abs(c.totalMass()/expected_core_mass-1)<1e-10?"true":"false")
                <<",\"all_inspected_states_finite\":true,\"min_core_J\":"<<min_core_J<<",\"min_PD_J\":"<<min_pd_J
                <<",\"max_fragment_force_residual\":"<<maximum_fragment_force_error<<",\"max_fragment_torque_residual\":"<<maximum_fragment_torque_error
                <<",\"max_fragment_power_residual\":"<<maximum_fragment_power_error<<",\"quasistatic_equilibrium_claimed\":false,\"experimental_curve_claimed\":false}\n";
            std::cout<<"CONTINUATION_CHECKPOINT passed; diagnostic only, no equilibrated curve claimed"<<std::endl;return 0;
        }
        Real dpc_dt=std::numeric_limits<Real>::infinity();
        for(size_t i=0;i<c.n;++i)dpc_dt=std::min(dpc_dt,o.cfl*1.15*o.spacing/(dpc_material.acousticBound(c.F[i].determinant())+c.v[i].norm()));
        Real quadratic_dt=o.quadratic_kinematics?(energy_model?energy_model->finiteStrainTimeStep(quadratic_reference.elastic_dt_bound):quadratic_reference.elastic_dt_bound):std::numeric_limits<Real>::infinity();
        Real fragment_dt=growth?fragment_contact->stable_dt:std::numeric_limits<Real>::infinity();
        Real energy_dt=energy_model?energy_model->stabilizer_dt:std::numeric_limits<Real>::infinity();
        Real dt=std::min({o.max_dt,core_dt.exec(),pd_dt.exec(),constraint_dt,radial_dt,dpc_dt,quadratic_dt,fragment_dt,energy_dt,core_residual_reference?core_residual_reference->native_dt_bound:std::numeric_limits<Real>::infinity()})*o.dt_scale;
        if(!std::isfinite(dt) || dt<1e-15 || step>=o.max_steps)throw std::runtime_error("Time step/count guard");
        if(step==0)write_fields(0);
        if(step%1000==0)
        {progress<<step<<','<<stage<<','<<stage_time<<','<<dt<<','<<delta<<','<<contact.upper<<','<<contact.lower<<','<<kinetic<<','<<tie.energy<<','<<contact.energy<<','<<jc<<','<<js<<'\n';progress.flush();write_friction("progress",delta);}
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
                if(growth){fracture_history<<std::setprecision(17)<<stage<<','<<delta/o.diameter<<','<<onset->total_physical_broken_bonds<<','<<onset->maximum_damage<<','
                    <<fracture_reference->ballistic_particles<<','<<fracture_reference->resolution_removed_bonds<<','<<fracture_reference->affine_rows<<','<<fracture_inner->inactive_points<<','<<fracture_outer->inactive_points<<','
                    <<fragment_contact->self_contacts<<','<<fragment_contact->core_contacts<<','<<fragment_contact->plane_contacts<<','<<fragment_contact->force_residual<<','<<fragment_contact->torque_residual<<','<<fragment_contact->power_residual<<'\n';fracture_history.flush();}
                if(energy_model)write_energy(energy_history,delta);
                if(core_residual_reference)
                {Real U=0,W=0;auto*p=cp.getVariableDataByName<Real>("DPCElasticEnergyDensity");for(size_t i=0;i<c.n;++i){U+=c.measure[i]*p[i];W+=c.measure[i]*dpc_work[i];}Real S=core_residual_reference->totalEnergy();residual_history<<stage<<','<<delta/o.diameter<<','<<S<<','<<U<<','<<W<<','<<S/std::max(U+W,Real(1e-30))<<'\n';residual_history.flush();}
                write_friction("equilibrium",delta);
                write_fields(stage);
                std::cout<<"ACCEPTED load="<<stage<<'/'<<o.load_steps<<" strain="<<delta/o.diameter<<" force="<<mean<<" residual="<<residual<<" imbalance="<<imbalance<<std::endl;
                if(stage==o.load_steps)break;
                ++stage;stage_time=0;stage_step=0;consecutive=0;continue;
            }
        }
        if(stage_time>o.max_stage_duration)throw std::runtime_error("Hybrid load step did not reach equilibrium;no complete result written");
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){if(core_residual_first)core_residual_first->initialization(i,dt);else core_first.initialization(i,dt);});
        pd_first.exec(dt);
        apply(delta_at(stage_time+.5*dt),plate_speed(stage_time+.5*dt),true);
        // PD first-half copied old prior acceleration before drifting. Replace
        // it here with the common-midpoint interface/contact loads exactly once.
        for(size_t i=0;i<s.n;++i)pd_acc[i]=pd_port.external[i];
        try{pd_second(dt);pd_third(dt);pd_auxiliary(dt);}
        catch(const std::exception&e){if(energy_model)dump_energy_failure(e.what());throw;}
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){if(core_residual_first)core_residual_first->interaction(i,dt);else core_first.interaction(i,dt);});
        particle_for(execution::ParallelPolicy(),core.LoopRange(),[&](size_t i){if(core_residual_first)core_residual_first->update(i,dt);else core_first.update(i,dt);});
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
        <<",\n\"plate_friction_mu\":"<<o.plate_friction<<",\n\"plate_tangent_factor\":"<<o.plate_tangent_factor
        <<",\n\"friction_slip_dissipation_J\":"<<contact.slip_dissipation<<",\n\"friction_contact_release_J\":"<<contact.release_energy
        <<",\n\"friction_tangent_penalty_J\":"<<contact.tangent_energy<<",\n\"friction_history_commits\":"<<contact.history_commits
        <<",\n\"material_calibrated\":false,\n\"thin_layer_spatial_convergence\":false,\n\"full_stored_energy_audit\":false\n}\n";
    if(growth){std::ofstream f(out/"fracture_summary.json");f<<std::setprecision(17)
        <<"{\"status\":\"COMPLETED_STRENGTH_CONTROLLED_FRACTURE_GROWTH\",\"events\":"<<fracture_event_count
        <<",\"physical_broken_bonds\":"<<onset->total_physical_broken_bonds<<",\"resolution_removed_bonds\":"<<fracture_reference->resolution_removed_bonds
        <<",\"ballistic_particles\":"<<fracture_reference->ballistic_particles<<",\"inactive_inner_points\":"<<fracture_inner->inactive_points
        <<",\"inactive_outer_points\":"<<fracture_outer->inactive_points<<",\"shell_mass_preserved\":"<<(std::abs(s.totalMass()/expected_shell_mass-1)<1e-10?"true":"false")
        <<",\"max_fragment_force_residual\":"<<maximum_fragment_force_error<<",\"max_fragment_torque_residual\":"<<maximum_fragment_torque_error
        <<",\"max_fragment_power_residual\":"<<maximum_fragment_power_error<<",\"fracture_energy_regularized\":false,\"experimentally_calibrated_fracture\":false}\n";}
    if(energy_model){auto a=energy_model->statistics();std::ofstream f(out/"energy_summary.json");f<<std::setprecision(17)
        <<"{\"status\":\"COMPLETED_NOSB_PD_ENERGY_DAMAGE_COMPRESSION\",\"ft_Pa\":"<<o.fracture_stress<<",\"Gc_J_m2\":"<<o.fracture_energy
        <<",\"band_width_m\":"<<energy_model->fracture.width<<",\"initiated_particles\":"<<a.initiated<<",\"tensile_failed_particles\":"<<a.failed
        <<",\"maximum_damage\":"<<a.max_damage<<",\"stored_J\":"<<a.stored<<",\"stabilization_J\":"<<a.stabilization<<",\"fracture_dissipated_J\":"<<a.dissipated
        <<",\"equivalent_fracture_area_m2\":"<<a.dissipated/o.fracture_energy<<",\"numerical_erosion_bonds\":0,\"interface_releases\":0,\"sharp_crack_geometry_claimed\":false,\"core_fracture_implemented\":false,\"material_calibrated\":false}\n";}
    return 0;
}

int main(int argc,char **argv)
{
    try{HybridOptions options;options.parseHybrid(argc,argv);
        // Match the arena concurrency to the configured worker limit. The
        // implicit hardware-sized arena can retain affinity proxy tasks when
        // oneTBB global_control limits the workers (upstream issue #1403).
        tbb::task_arena arena(options.threads);
        return arena.execute([&](){return runHybrid(options);});}
    catch(const std::exception &e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
}

#pragma once
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
#include <sstream>
#include <map>
#include <tbb/task_arena.h>

struct HybridOptions : PDSolidCapsuleOptions
{
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


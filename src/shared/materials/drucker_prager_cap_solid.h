/* Adapter for the classic SPHinXsys total-Lagrangian solid stress path.
 * Not the isochoric DecomposedPlasticIntegration1stHalf path.
 * Constitutive updates must be called once per accepted integration midpoint.
 */
#ifndef SPH_DRUCKER_PRAGER_CAP_SOLID_H
#define SPH_DRUCKER_PRAGER_CAP_SOLID_H
#include "elastic_solid.h"
#include "base_particles.hpp"
#include "drucker_prager_cap_kernel.h"

namespace SPH {
class DruckerPragerCapSolid : public SaintVenantKirchhoffSolid
{
    dpc::Model model_;
    Real acoustic_modulus_bound_ = 0;
    Matd *inverse_cp_ = nullptr, *cauchy_ = nullptr;
    Real *kappa_ = nullptr, *relative_density_ = nullptr, *plastic_work_ = nullptr;
    Real *elastic_energy_ = nullptr, *yield_residual_ = nullptr, *active_surface_ = nullptr;
    Real *alpha_ = nullptr, *effective_d_ = nullptr;
  public:
    DruckerPragerCapSolid(Real rho0,const dpc::Model &model)
        : SaintVenantKirchhoffSolid(rho0,model.parameters(model.r0()).E,model.parameters(model.r0()).nu),model_(model)
    {
        material_type_name_="DruckerPragerCapSolid";
        // The table is immutable for this material instance. Cache exactly
        // the same global modulus maximum used by acousticBound instead of
        // scanning every table node for every particle at every time step.
        for(const auto &row:model_.nodes())
        {
            const auto p=model_.parameters(row.r);
            acoustic_modulus_bound_=std::max(acoustic_modulus_bound_,p.K+4*p.G/3);
        }
    }
    void initializeLocalParameters(BaseParticles *particles) override
    {
        SaintVenantKirchhoffSolid::initializeLocalParameters(particles);
        inverse_cp_=particles->registerStateVariableData<Matd>("DPCInversePlasticMetric",[](size_t)->Matd{return Matd::Identity();});
        cauchy_=particles->registerStateVariableData<Matd>("DPCCauchyStress");
        kappa_=particles->registerStateVariableData<Real>("DPCPlasticCompaction");
        relative_density_=particles->registerStateVariableData<Real>("DPCRelativeDensity",model_.r0());
        plastic_work_=particles->registerStateVariableData<Real>("DPCPlasticStressWork");
        elastic_energy_=particles->registerStateVariableData<Real>("DPCElasticEnergyDensity");
        yield_residual_=particles->registerStateVariableData<Real>("DPCYieldResidual");
        active_surface_=particles->registerStateVariableData<Real>("DPCActiveSurface");
        alpha_=particles->registerStateVariableData<Real>("DPCAccumulatedPlasticShear");
        effective_d_=particles->registerStateVariableData<Real>("DPCEffectiveCohesion",model_.parameters(model_.r0()).d);
        particles->addEvolvingVariable<Matd>("DPCInversePlasticMetric");
        particles->addEvolvingVariable<Real>("DPCPlasticCompaction");
        particles->addEvolvingVariable<Real>("DPCPlasticStressWork");
        particles->addEvolvingVariable<Real>("DPCAccumulatedPlasticShear");
    }
    Matd StressPK2(Matd &F,size_t i) override
    {
        if(Dimensions!=3)throw std::runtime_error("DPC classic adapter currently qualified for 3-D only");
        dpc::State old;old.inverse_Cp=inverse_cp_[i];old.kappa=kappa_[i];old.alpha=alpha_[i];
        const auto out=model_.update(F,old);
        inverse_cp_[i]=out.state.inverse_Cp;kappa_[i]=out.state.kappa;cauchy_[i]=out.sigma;
        relative_density_[i]=out.relative_density;plastic_work_[i]+=out.plastic_stress_work;
        elastic_energy_[i]=out.elastic_energy;yield_residual_[i]=out.yield_residual;active_surface_[i]=out.surface;
        alpha_[i]=out.state.alpha;effective_d_[i]=out.effective_d;
        const Matd invF=F.inverse();return F.determinant()*invF*out.sigma*invF.transpose();
    }
    Matd StressCauchy(Matd &,size_t) override
    { throw std::logic_error("DPC requires the total-deformation PK1/PK2 integration path; strain-only and shell adapters are unsupported"); }
    Real VolumetricKirchhoff(Real) override
    { throw std::logic_error("DPC must not use an elastic/isochoric decomposed stress integration path"); }
    Matd DeviatoricKirchhoff(const Matd &) override
    { throw std::logic_error("DPC must not use a separated elastic deviatoric stress update"); }
    Real acousticBound(Real J) const
    {
        return std::sqrt(acoustic_modulus_bound_*J/rho0_);
    }
    Real ShearHardeningModulus() const {return model_.shearHardening();}
};
} // namespace SPH
#endif

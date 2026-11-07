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
    Matd *inverse_cp_ = nullptr, *cauchy_ = nullptr;
    Real *kappa_ = nullptr, *relative_density_ = nullptr, *plastic_work_ = nullptr;
    Real *elastic_energy_ = nullptr, *yield_residual_ = nullptr, *active_surface_ = nullptr;
  public:
    DruckerPragerCapSolid(Real rho0,const dpc::Model &model)
        : SaintVenantKirchhoffSolid(rho0,model.parameters(model.r0()).E,model.parameters(model.r0()).nu),model_(model)
    { material_type_name_="DruckerPragerCapSolid"; }
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
        particles->addEvolvingVariable<Matd>("DPCInversePlasticMetric");
        particles->addEvolvingVariable<Real>("DPCPlasticCompaction");
        particles->addEvolvingVariable<Real>("DPCPlasticStressWork");
    }
    Matd StressPK2(Matd &F,size_t i) override
    {
        if(Dimensions!=3)throw std::runtime_error("DPC classic adapter currently qualified for 3-D only");
        dpc::State old;old.inverse_Cp=inverse_cp_[i];old.kappa=kappa_[i];
        const auto out=model_.update(F,old);
        inverse_cp_[i]=out.state.inverse_Cp;kappa_[i]=out.state.kappa;cauchy_[i]=out.sigma;
        relative_density_[i]=out.relative_density;plastic_work_[i]+=out.plastic_stress_work;
        elastic_energy_[i]=out.elastic_energy;yield_residual_[i]=out.yield_residual;active_surface_[i]=out.surface;
        const Matd invF=F.inverse();return F.determinant()*invF*out.sigma*invF.transpose();
    }
    Real acousticBound(Real J) const
    {
        Real maxmod=0;for(const auto &row:model_.nodes())
        {auto p=model_.parameters(row.r);maxmod=std::max(maxmod,p.K+4*p.G/3);}
        return std::sqrt(maxmod*J/rho0_);
    }
};
} // namespace SPH
#endif

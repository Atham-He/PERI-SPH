#include "pd_gel_material.h"
#include "pd_gel_kernel.h"
#include "base_particles.hpp"
#include <stdexcept>

namespace SPH::peridynamics
{
void requireGelStatus(GelKernelStatus status)
{
    switch (status)
    {
    case GelKernelStatus::OK:
        return;
    case GelKernelStatus::InvalidDeformation:
        throw std::runtime_error("gel deformation determinant below 0.05 or non-finite");
    case GelKernelStatus::InvalidHistory:
        throw std::runtime_error("invalid gel material history");
    case GelKernelStatus::InvalidDissipation:
        throw std::runtime_error("invalid gel dissipation history");
    case GelKernelStatus::InvalidIncrement:
        throw std::invalid_argument("gel material time increment must be finite and nonnegative");
    case GelKernelStatus::InvalidRelaxation:
        throw std::invalid_argument("gel relaxation increment outside representable range");
    case GelKernelStatus::InvalidDissipationIntegral:
        throw std::runtime_error("gel dissipation integral invalid; no clipping applied");
    case GelKernelStatus::NonfiniteResponse:
        throw std::runtime_error("non-finite gel response");
    }
    throw std::logic_error("unknown gel material status");
}

GelParameters GelParameters::fromYoungPoisson(Real young, Real poisson, Real shear_ratio, Real tau)
{
    if (!std::isfinite(young) || young <= 0 || !std::isfinite(poisson) || poisson <= -1 || poisson >= 0.5 ||
        !std::isfinite(shear_ratio) || shear_ratio < 0)
        throw std::invalid_argument("invalid equilibrium gel moduli or relaxing shear ratio");
    const Real shear = young / (2 * (1 + poisson));
    GelParameters parameters{young / (3 * (1 - 2 * poisson)), shear, shear * shear_ratio, tau};
    parameters.validate();
    return parameters;
}

void GelParameters::validate() const
{
    if (Dimensions != 3)
        throw std::invalid_argument("the effective gel kernel requires three dimensions");
    if (!std::isfinite(bulk) || bulk <= 0 || !std::isfinite(equilibrium_shear) || equilibrium_shear <= 0 ||
        !std::isfinite(branch_shear) || branch_shear < 0 || !std::isfinite(relaxation_time) || relaxation_time <= 0)
        throw std::invalid_argument("gel requires positive finite bulk/equilibrium shear/tau and nonnegative branch shear");
}

Real GelParameters::instantaneousYoung() const
{
    const Real shear = equilibrium_shear + branch_shear;
    return 9 * bulk * shear / (3 * bulk + shear);
}

Real GelParameters::instantaneousPoisson() const
{
    const Real shear = equilibrium_shear + branch_shear;
    return (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear));
}

Matd gelDeviatoricGreen(const Matd &F)
{
    if (!validGelDeformation(F, F.determinant()))
        requireGelStatus(GelKernelStatus::InvalidDeformation);
    return deviatoricGelGreen(F);
}

GelResponse gelResponse(const Matd &F, const GelParameters &p, const GelHistory &history)
{
    p.validate();
    GelResponse result;
    requireGelStatus(gelResponseKernel(F, p, history, result));
    return result;
}

GelHistory advanceGel(const Matd &F, Real dt, const GelParameters &p, const GelHistory &old)
{
    p.validate();
    GelHistory next;
    requireGelStatus(advanceGelKernel(F, dt, p, old, next));
    return next;
}
} // namespace SPH::peridynamics

namespace SPH
{
GelViscoelasticSolid::GelViscoelasticSolid(Real density, Real young, Real poisson, Real ratio, Real tau)
    : HughesWingetSolid(density,
        peridynamics::GelParameters::fromYoungPoisson(young, poisson, ratio, tau).instantaneousYoung(),
        peridynamics::GelParameters::fromYoungPoisson(young, poisson, ratio, tau).instantaneousPoisson()),
      parameters_(peridynamics::GelParameters::fromYoungPoisson(young, poisson, ratio, tau))
{
    material_type_name_ = "GelViscoelasticSolid";
    // Use the initial longitudinal modulus, including the relaxing branch.
    c0_ = std::sqrt((parameters_.bulk + 4 * (parameters_.equilibrium_shear + parameters_.branch_shear) / 3) / density);
    setContactStiffness(density, c0_);
}

void GelViscoelasticSolid::initializeLocalParameters(BaseParticles *particles)
{
    HughesWingetSolid::initializeLocalParameters(particles);
    for (const char *name : {"GelGreenStrain", "GelViscousStrain", "GelMidpointGreenStrain", "GelMidpointViscousStrain", "GelAcceptedStress"})
    {
        particles->registerStateVariableData<Matd>(name);
        particles->addEvolvingVariable<Matd>(name);
    }
    for (const char *name : {"GelDissipationDensity", "GelMidpointDissipationDensity", "GelEnergyDensity", "GelHistoryTime"})
    {
        particles->registerStateVariableData<Real>(name);
        particles->addEvolvingVariable<Real>(name);
    }
}
} // namespace SPH

#include "pd_gel_material.h"
#include "base_particles.hpp"
#include <stdexcept>

namespace SPH::peridynamics
{
namespace
{
Real requireDeformation(const Matd &F)
{
    const Real J = F.determinant();
    if (!F.allFinite() || !std::isfinite(J) || J < 0.05)
        throw std::runtime_error("gel deformation determinant below 0.05 or non-finite");
    return J;
}

void requireHistory(const GelHistory &history)
{
    for (const Matd *tensor : {&history.green, &history.viscous})
        if (!tensor->allFinite() || (tensor->transpose() - *tensor).norm() > 1e-12 * (1 + tensor->norm()) ||
            std::abs(tensor->trace()) > 1e-12 * (1 + tensor->norm()))
            throw std::runtime_error("invalid gel material history");
    if (!std::isfinite(history.dissipation) || history.dissipation < 0)
        throw std::runtime_error("invalid gel dissipation history");
}
} // namespace

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
    requireDeformation(F);
    const Matd E = 0.5 * (F.transpose() * F - Matd::Identity());
    return E - (E.trace() / Dimensions) * Matd::Identity();
}

GelResponse gelResponse(const Matd &F, const GelParameters &p, const GelHistory &history)
{
    p.validate();
    requireHistory(history);
    const Real J = requireDeformation(F);
    const Real isochoric = std::pow(J, -2.0 / 3.0);
    const Real trace_C = F.squaredNorm();
    const Matd inverse_transpose = F.inverse().transpose();
    const Matd q = gelDeviatoricGreen(F) - history.viscous;
    GelResponse result;
    result.free_energy = 0.5 * p.equilibrium_shear * (isochoric * trace_C - 3) +
        0.25 * p.bulk * (J * J - 1 - 2 * std::log(J)) + p.branch_shear * q.squaredNorm();
    result.pk1 = p.equilibrium_shear * isochoric * (F - trace_C / 3 * inverse_transpose) +
        0.5 * p.bulk * (J * J - 1) * inverse_transpose + 2 * p.branch_shear * F * q;
    result.cauchy = result.pk1 * F.transpose() / J;
    if (!result.pk1.allFinite() || !result.cauchy.allFinite() || !std::isfinite(result.free_energy))
        throw std::runtime_error("non-finite gel response");
    return result;
}

GelHistory advanceGel(const Matd &F, Real dt, const GelParameters &p, const GelHistory &old)
{
    p.validate();
    requireHistory(old);
    const Matd green = gelDeviatoricGreen(F);
    if (!std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("gel material time increment must be finite and nonnegative");
    if (dt == 0)
        return old;
    const Real z = dt / p.relaxation_time;
    if (!std::isfinite(z) || z <= 0)
        throw std::invalid_argument("gel relaxation increment outside representable range");
    const Real decay = std::exp(-z), one_minus = -std::expm1(-z);
    const Matd q0 = old.green - old.viscous, change = green - old.green;
    const Matd q1 = decay * q0 + (one_minus / z) * change;
    // Integral of (1-exp(-s))^2 / z^2. The series avoids cancellation.
    const Real w22 = z < 1e-3 ? z * (1.0 / 3 + z * (-1.0 / 4 + z * (7.0 / 60 +
        z * (-1.0 / 24 + z * (31.0 / 2520 - z / 320))))) :
        (1 - (one_minus + 0.5 * one_minus * one_minus) / z) / z;
    const Real increment = 2 * p.branch_shear * (
        (-0.5 * std::expm1(-2 * z)) * q0.squaredNorm() +
        (one_minus * one_minus / z) * q0.cwiseProduct(change).sum() + w22 * change.squaredNorm());
    if (!std::isfinite(increment) || increment < 0)
        throw std::runtime_error("gel dissipation integral invalid; no clipping applied");
    GelHistory next{green, green - q1, old.dissipation + increment};
    requireHistory(next);
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

#ifndef PD_ENERGY_AUDIT_H
#define PD_ENERGY_AUDIT_H

#include "base_particles.hpp"
#include "base_local_dynamics.h"
#include <cmath>
#include <stdexcept>

namespace SPH
{
namespace peridynamics
{
/** Specific barotropic storage for p=c^2(rho-rho0), with zero at rho0.
 * Integrates de=p/rho^2 drho. This is mechanical compression energy, not heat.
 */
inline Real linearEOSSpecificEnergy(Real rho, Real rho0, Real c)
{
    if (!(rho > 0) || !(rho0 > 0) || !(c > 0))
        throw std::invalid_argument("EOS energy requires positive density and sound speed");
    Real x = rho / rho0 - 1;
    // Avoid cancellation near rest, where the energy is O(x^2).
    if (std::abs(x) < 1e-4)
        return c * c * x * x * (0.5 + x * (-2.0 / 3.0 + x * (0.75 - 0.8 * x)));
    return c * c * (std::log1p(x) - x / (1 + x));
}

/** Records work of the actual fluid velocity kick. Run immediately after
 * pressure.exec(dt), before any additional damper changes Velocity.
 * The midpoint velocity reconstructs the exact kinetic-energy increment.
 */
class FluidKineticWork : public LocalDynamics
{
  public:
    explicit FluidKineticWork(SPHBody &body) : LocalDynamics(body),
        mass_(particles_->getVariableDataByName<Real>("Mass")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        force_(particles_->getVariableDataByName<Vecd>("Force")),
        prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
        viscous_(particles_->getVariableDataByName<Vecd>("ViscousForce")),
        pressure_work_(particles_->registerStateVariableData<Real>("FluidPressureWork")),
        prior_work_(particles_->registerStateVariableData<Real>("FluidPriorWork")),
        viscous_work_(particles_->registerStateVariableData<Real>("FluidViscousWork"))
    {
        for (const char *name : {"FluidPressureWork", "FluidPriorWork", "FluidViscousWork"})
            particles_->addEvolvingVariable<Real>(name);
    }
    void update(size_t i, Real dt)
    {
        Vecd midpoint = vel_[i] - 0.5 * dt * (force_[i] + prior_[i]) / mass_[i];
        pressure_work_[i] += dt * force_[i].dot(midpoint);
        prior_work_[i] += dt * prior_[i].dot(midpoint);
        viscous_work_[i] += dt * viscous_[i].dot(midpoint);
    }
  protected:
    Real *mass_;
    Vecd *vel_, *force_, *prior_, *viscous_;
    Real *pressure_work_, *prior_work_, *viscous_work_;
};
} // namespace peridynamics
} // namespace SPH
#endif

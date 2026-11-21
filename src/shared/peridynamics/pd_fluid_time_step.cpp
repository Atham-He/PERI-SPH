#include "pd_fluid_time_step.h"
#include "weakly_compressible_fluid.h"

#include <cmath>
#include <typeinfo>

namespace SPH::fluid_dynamics
{
namespace
{
peridynamics::ParticleVector3 packVector(const Vecd &v)
{
    peridynamics::ParticleVector3 result{};
    for (int axis = 0; axis < Dimensions; ++axis)
        result.component[axis] = v[axis];
    return result;
}
}

FluidAcousticTimeStepExecution::FluidAcousticTimeStepExecution(
    SPHBody &body, peridynamics::FluidCflBatch *backend, Real cfl)
    : Base(body, cfl), backend_(backend), constant_sound_speed_(typeid(fluid_) == typeid(WeaklyCompressibleFluid)) {}

Real FluidAcousticTimeStepExecution::exec(Real dt)
{
    if (!backend_)
        return Base::exec(dt);

    setupDynamics(dt);
    sites_.resize(particles_->TotalRealParticles());
    const Real uniform_sound = constant_sound_speed_ ? fluid_.ReferenceSoundSpeed() : 0.0;
    for (std::size_t i = 0; i < sites_.size(); ++i)
    {
        const Real sound = constant_sound_speed_ ? uniform_sound : fluid_.getSoundSpeed(p_[i], rho_[i]);
        sites_[i] = {packVector(vel_[i]), sound};
    }

    const Real maximum = backend_->acousticMaximum(sites_, Reference());
    return std::isfinite(maximum) ? outputResult(maximum) : Base::exec(dt);
}

template <class AdvectionCriterion>
FluidAdvectionTimeStepExecution<AdvectionCriterion>::FluidAdvectionTimeStepExecution(
    SPHBody &body, peridynamics::FluidCflBatch *backend, Real reference_speed, Real cfl)
    : Base(body, reference_speed, cfl), backend_(backend) {}

template <class AdvectionCriterion>
Real FluidAdvectionTimeStepExecution<AdvectionCriterion>::exec(Real dt)
{
    if (!backend_)
        return Base::exec(dt);

    this->setupDynamics(dt);
    sites_.resize(this->particles_->TotalRealParticles());
    for (std::size_t i = 0; i < sites_.size(); ++i)
        sites_[i] = {packVector(this->vel_[i]), packVector(this->force_[i]),
                     packVector(this->force_prior_[i]), this->mass_[i]};

    const Real maximum = backend_->advectionMaximum(sites_, this->Reference(), this->h_min_);
    return std::isfinite(maximum) ? this->outputResult(maximum) : Base::exec(dt);
}

template class FluidAdvectionTimeStepExecution<AdvectionTimeStep>;
template class FluidAdvectionTimeStepExecution<AdvectionViscousTimeStep>;
}

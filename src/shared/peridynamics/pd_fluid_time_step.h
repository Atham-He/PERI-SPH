#pragma once

#include "dynamics_algorithms.h"
#include "fluid_time_step.h"
#include "pd_fluid_cfl_batch.h"

namespace SPH::fluid_dynamics
{
class FluidAcousticTimeStepExecution : public ReduceDynamics<AcousticTimeStep>
{
    using Base = ReduceDynamics<AcousticTimeStep>;

  public:
    FluidAcousticTimeStepExecution(SPHBody &body, peridynamics::FluidCflBatch *backend, Real cfl = 0.6);
    Real exec(Real dt = 0.0) override;

  private:
    peridynamics::FluidCflBatch *backend_;
    bool constant_sound_speed_;
    std::vector<peridynamics::FluidAcousticCflSite> sites_;
};

template <class AdvectionCriterion>
class FluidAdvectionTimeStepExecution : public ReduceDynamics<AdvectionCriterion>
{
    using Base = ReduceDynamics<AdvectionCriterion>;

  public:
    FluidAdvectionTimeStepExecution(SPHBody &body, peridynamics::FluidCflBatch *backend,
                                    Real reference_speed, Real cfl = 0.25);
    Real exec(Real dt = 0.0) override;

  private:
    peridynamics::FluidCflBatch *backend_;
    std::vector<peridynamics::FluidAdvectionCflSite> sites_;
};

extern template class FluidAdvectionTimeStepExecution<AdvectionTimeStep>;
extern template class FluidAdvectionTimeStepExecution<AdvectionViscousTimeStep>;
}

#pragma once

#include "dynamics_algorithms.h"
#include "pd_numerics.h"
#include "pd_solid_cfl_batch.h"

namespace SPH::solid_dynamics
{
class PDTimeStepExecution : public ReduceDynamics<PDTimeStep>
{
    using Base = ReduceDynamics<PDTimeStep>;

  public:
    PDTimeStepExecution(SPHBody &body, peridynamics::PDCflBatch *backend, Real cfl = 0.1)
        : Base(body, cfl), backend_(backend) {}

    Real exec(Real dt = 0.0) override
    {
        if (!backend_)
            return Base::exec(dt);

        setupDynamics(dt);
        sites_.resize(particles_->TotalRealParticles());
        for (std::size_t i = 0; i < sites_.size(); ++i)
        {
            auto &site = sites_[i];
            site = {};
            for (int axis = 0; axis < Dimensions; ++axis)
            {
                site.velocity.component[axis] = velocity_[i][axis];
                site.acceleration.component[axis] = acceleration_[i][axis];
            }
        }

        const auto inverse_step = backend_->maximumInverseStep(sites_, Reference(), spacing_, sound_speed_);
        return std::isfinite(inverse_step) ? outputResult(inverse_step) : Base::exec(dt);
    }

  private:
    peridynamics::PDCflBatch *backend_;
    std::vector<peridynamics::PDCflSite> sites_;
};
}

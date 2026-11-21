#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
enum class SolidKickInput
{
    Acceleration,
    Force
};

struct SolidDriftInput
{
    ParticleVector3 displacement, velocity, drive;
};

struct SolidDriftResult
{
    ParticleVector3 displacement, half_velocity;
    double damped_kinetic;
};

struct SolidFinishInput
{
    ParticleVector3 half_velocity, force, additional_force;
};

struct SolidFinishResult
{
    ParticleVector3 velocity, acceleration;
    double undamped_kinetic;
};

class SolidTimeStepBatch
{
  public:
    virtual ~SolidTimeStepBatch() = default;
    // No accepted state is changed by a trial drift. Callers can reject it
    // and retry from the same old fields with a different time step.
    virtual void predict(const std::vector<SolidDriftInput> &input, double dt, double damping_factor,
                         SolidKickInput drive, std::vector<SolidDriftResult> &output) = 0;
    // The additional force is included in the stored acceleration. Its
    // impulse during the drift is owned by the calling operator.
    virtual void finish(const std::vector<SolidFinishInput> &input, double dt, double damping_factor,
                        SolidKickInput drive, std::vector<SolidFinishResult> &output) = 0;
};

std::unique_ptr<SolidTimeStepBatch> makeCudaSolidTimeStep(const std::vector<double> &mass);
}

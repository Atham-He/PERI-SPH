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

// NOSB-PD uses a drift-kick-drift update with work evaluated at the
// midpoint velocity; the solid SPH integrator above uses kick-drift-kick.
struct PositionVerletStart
{
    ParticleVector3 position, velocity;
};

struct PositionVerletFinish
{
    ParticleVector3 position, velocity, acceleration;
    ParticleVector3 internal, hourglass, damping, viscous, friction, acoustic;
    double mass, boundary_power, friction_boundary_power;
};

struct PositionVerletResult
{
    ParticleVector3 position, velocity;
    double internal_work, external_work, hourglass_work, damping_work;
    double viscous_work, acoustic_work, boundary_work, friction_loss;
};

class PositionVerletBatch
{
  public:
    virtual ~PositionVerletBatch() = default;
    virtual void start(const std::vector<PositionVerletStart> &state, double dt,
                       std::vector<ParticleVector3> &position) = 0;
    virtual void finish(const std::vector<PositionVerletFinish> &state, double dt,
                        std::vector<PositionVerletResult> &result) = 0;
};

std::unique_ptr<PositionVerletBatch> makeCudaPositionVerlet();
std::unique_ptr<SolidTimeStepBatch> makeCudaSolidTimeStep(const std::vector<double> &mass);
}

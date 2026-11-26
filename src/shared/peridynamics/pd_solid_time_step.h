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

struct DeviceCoupledMaterialForces
{
    DeviceParticleVectorView first_adjoint, second_adjoint, second_stabilizing_force;

    explicit operator bool() const noexcept
    {
        return first_adjoint && second_adjoint && second_stabilizing_force;
    }
};

struct SolidMotion
{
    ParticleVector3 displacement, velocity;
};

struct SolidCoupledFinishResult
{
    ParticleVector3 velocity;
    double undamped_kinetic;
};

struct SolidImplicitCorrection
{
    ParticleVector3 displacement;
    double error_squared, increment_squared;
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

    // Coupled operators publish forces on the same device and stream. The
    // accepted motion is retained across trials and is changed only on finish.
    virtual const double *deviceMass() const noexcept { return nullptr; }
    virtual void setCoupledDrive(DeviceParticleVectorView, SolidKickInput,
                                 DeviceParticleVectorView, SolidKickInput,
                                 const std::vector<ParticleVector3> &, DeviceCoupledMaterialForces,
                                 const std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This solid integrator does not consume device forces");
    }
    virtual void readCoupledDrive(std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This solid integrator does not store device forces");
    }
    virtual void configureCoupledReference(const std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This solid integrator does not store device reference geometry");
    }
    virtual DeviceSolidKinematics coupledKinematics(bool)
    {
        throw std::logic_error("This solid integrator does not expose device kinematics");
    }
    virtual void initializeCoupledMotion(const std::vector<SolidMotion> &)
    {
        throw std::logic_error("This solid integrator does not retain device motion");
    }
    virtual void predictCoupled(double, double, bool, std::vector<SolidDriftResult> &)
    {
        throw std::logic_error("This solid integrator does not predict coupled motion");
    }
    virtual void correctCoupledTrial(double, std::vector<SolidImplicitCorrection> &)
    {
        throw std::logic_error("This solid integrator does not correct coupled motion");
    }
    virtual void advanceCoupledTrial()
    {
        throw std::logic_error("This solid integrator does not support implicit coupled trials");
    }
    virtual void readCoupledAcceleration(std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This solid integrator does not store coupled acceleration");
    }
    virtual void finishCoupled(double, double, std::vector<SolidCoupledFinishResult> &)
    {
        throw std::logic_error("This solid integrator does not finish coupled motion");
    }
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

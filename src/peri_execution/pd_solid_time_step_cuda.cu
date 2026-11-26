#include "pd_solid_time_step.h"
#include "pd_cuda_buffer.cuh"

#include <cmath>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

__device__ double normSquared(const ParticleVector3 &vector)
{
    const auto &x = vector.component;
    return x[0]*x[0] + (x[1]*x[1] + x[2]*x[2]);
}

__device__ SolidDriftResult predictMotion(const ParticleVector3 &displacement,
                                          const ParticleVector3 &velocity, const ParticleVector3 &old_drive,
                                          double mass, double dt, double factor, SolidKickInput drive)
{
    SolidDriftResult next;
    ParticleVector3 damped_velocity;
    for (int axis = 0; axis < 3; ++axis)
    {
        damped_velocity.component[axis] = factor * velocity.component[axis];
        double kick = (.5 * dt) * old_drive.component[axis];
        if (drive == SolidKickInput::Force)
            kick /= mass;
        next.half_velocity.component[axis] = damped_velocity.component[axis] + kick;
        next.displacement.component[axis] = displacement.component[axis] + dt * next.half_velocity.component[axis];
    }
    next.damped_kinetic = .5 * mass * normSquared(damped_velocity);
    return next;
}

__global__ void predictSolidStep(const double *mass, const SolidDriftInput *input, std::size_t count,
                                 double dt, double factor, SolidKickInput drive, SolidDriftResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
    {
        const auto &old = input[i];
        output[i] = predictMotion(old.displacement, old.velocity, old.drive, mass[i], dt, factor, drive);
    }
}

__device__ SolidFinishResult finishMotion(const ParticleVector3 &half_velocity, const ParticleVector3 &force,
                                           const ParticleVector3 &additional_force, double mass,
                                           double dt, double factor, SolidKickInput drive)
{
    SolidFinishResult next;
    ParticleVector3 undamped_velocity;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double value = force.component[axis];
        const double additional = additional_force.component[axis];
        next.acceleration.component[axis] = (value + additional) / mass;
        // Both forms have the same discrete impulse. Keep the arithmetic
        // order used by the corresponding CPU stepping operator.
        const double kick = drive == SolidKickInput::Force ? (.5 * dt) * value / mass
            : (.5 * dt) * (value / mass);
        undamped_velocity.component[axis] = half_velocity.component[axis] + kick;
        next.velocity.component[axis] = factor * undamped_velocity.component[axis];
    }
    next.undamped_kinetic = .5 * mass * normSquared(undamped_velocity);
    return next;
}

__global__ void finishSolidStep(const double *mass, const SolidFinishInput *input, std::size_t count,
                                double dt, double factor, SolidKickInput drive, SolidFinishResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
    {
        const auto &trial = input[i];
        output[i] = finishMotion(trial.half_velocity, trial.force, trial.additional_force,
                                 mass[i], dt, factor, drive);
    }
}

__global__ void assembleCoupledDrive(const double *mass, const ParticleVector3 *internal,
                                      DeviceCoupledMaterialForces material,
                                      DeviceParticleVectorView first, SolidKickInput first_units,
                                      DeviceParticleVectorView second, SolidKickInput second_units,
                                      ParticleVector3 *drive)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= first.count + second.count)
        return;
    const bool on_first = i < first.count;
    const auto &external = on_first ? first.data[i] : second.data[i - first.count];
    const auto units = on_first ? first_units : second_units;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double value = external.component[axis];
        double material_force;
        if (internal)
            material_force = internal[i].component[axis];
        else if (on_first)
            material_force = material.first_adjoint.data[i].component[axis];
        else
        {
            const auto j = i - first.count;
            material_force = material.second_adjoint.data[j].component[axis] -
                             material.second_stabilizing_force.data[j].component[axis];
        }
        drive[i].component[axis] = (units == SolidKickInput::Force ? value : mass[i] * value) - material_force;
    }
}

__global__ void predictCoupledStep(const double *mass, const SolidMotion *motion, const ParticleVector3 *drive,
                                   std::size_t count, double dt, double factor, SolidDriftResult *output,
                                   ParticleVector3 *trial)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
    {
        const auto &old = motion[i];
        const auto next = predictMotion(old.displacement, old.velocity, drive[i], mass[i], dt, factor, SolidKickInput::Force);
        output[i] = next;
        trial[i] = next.displacement;
    }
}

__global__ void correctCoupledStep(const double *mass, const SolidMotion *old, const SolidDriftResult *drift,
                                   const ParticleVector3 *old_force, const ParticleVector3 *force,
                                   const ParticleVector3 *trial, std::size_t count, double dt,
                                   SolidImplicitCorrection *result, ParticleVector3 *candidate)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    ParticleVector3 next, error, increment;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double reference = old[i].displacement.component[axis];
        const double change = (force[i].component[axis] - old_force[i].component[axis]) / mass[i];
        next.component[axis] = reference + dt * (drift[i].half_velocity.component[axis] + (.25 * dt) * change);
        error.component[axis] = next.component[axis] - trial[i].component[axis];
        increment.component[axis] = trial[i].component[axis] - reference;
    }
    candidate[i] = next;
    result[i] = {next, normSquared(error), normSquared(increment)};
}

template <bool Trial>
__global__ void coupledGeometry(const ParticleVector3 *reference, const SolidMotion *motion,
                                 const ParticleVector3 *trial, const SolidDriftResult *drift,
                                 std::size_t count, ParticleVector3 *displacement,
                                 ParticleVector3 *position, ParticleVector3 *velocity)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    const ParticleVector3 d = Trial ? trial[i] : motion[i].displacement;
    displacement[i] = d;
    velocity[i] = Trial ? drift[i].half_velocity : motion[i].velocity;
    for (int axis = 0; axis < 3; ++axis)
        position[i].component[axis] = reference[i].component[axis] + d.component[axis];
}

__global__ void finishCoupledStep(const double *mass, SolidMotion *motion, const SolidDriftResult *drift,
                                  const ParticleVector3 *force, const ParticleVector3 *trial,
                                  std::size_t count, double dt, double factor,
                                  SolidCoupledFinishResult *output, ParticleVector3 *acceleration)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
    {
        const auto next = finishMotion(drift[i].half_velocity, force[i], ParticleVector3{}, mass[i],
                                       dt, factor, SolidKickInput::Force);
        output[i] = {next.velocity, next.undamped_kinetic};
        acceleration[i] = next.acceleration;
        motion[i] = {trial[i], next.velocity};
    }
}

__device__ double dot(const ParticleVector3 &a, const ParticleVector3 &b)
{
    return a.component[0] * b.component[0] +
        (a.component[1] * b.component[1] + a.component[2] * b.component[2]);
}

__global__ void startPositionVerlet(const PositionVerletStart *state, std::size_t count, double dt,
                                    ParticleVector3 *position)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    for (int axis = 0; axis < 3; ++axis)
        position[i].component[axis] = state[i].position.component[axis] +
            state[i].velocity.component[axis] * dt * 0.5;
}

__global__ void finishPositionVerlet(const PositionVerletFinish *state, std::size_t count, double dt,
                                     PositionVerletResult *result)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &old = state[i];
    auto &next = result[i];
    ParticleVector3 displacement, external;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double acceleration = old.acceleration.component[axis];
        const double velocity = old.velocity.component[axis];
        displacement.component[axis] = dt * (velocity + 0.5 * dt * acceleration);
        external.component[axis] = old.mass * acceleration - old.internal.component[axis] -
            old.hourglass.component[axis] - old.damping.component[axis];
        next.velocity.component[axis] = velocity + acceleration * dt;
        next.position.component[axis] = old.position.component[axis] + next.velocity.component[axis] * dt * 0.5;
    }
    next.internal_work = dot(old.internal, displacement);
    next.external_work = dot(external, displacement);
    next.hourglass_work = dot(old.hourglass, displacement);
    next.damping_work = dot(old.damping, displacement);
    next.viscous_work = dot(old.viscous, displacement);
    next.acoustic_work = dot(old.acoustic, displacement);
    next.boundary_work = dt * old.boundary_power;
    next.friction_loss = dt * old.friction_boundary_power - dot(old.friction, displacement);
}

class CudaPositionVerlet final : public PositionVerletBatch
{
  public:
    CudaPositionVerlet()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select position Verlet device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read position Verlet grid limit");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void start(const std::vector<PositionVerletStart> &state, double dt,
               std::vector<ParticleVector3> &position) override
    {
        const auto blocks = grid(state.size(), dt);
        position.resize(state.size());
        if (!blocks)
            return;
        upload(start_, state, "position Verlet start state");
        position_.reserve(state.size());
        startPositionVerlet<<<blocks, threads_>>>(start_.get(), state.size(), dt, position_.get());
        check(cudaGetLastError(), "start position Verlet");
        download(position, position_, "position Verlet starting positions");
    }

    void finish(const std::vector<PositionVerletFinish> &state, double dt,
                std::vector<PositionVerletResult> &result) override
    {
        const auto blocks = grid(state.size(), dt);
        result.resize(state.size());
        if (!blocks)
            return;
        upload(finish_, state, "position Verlet finish state");
        result_.reserve(state.size());
        finishPositionVerlet<<<blocks, threads_>>>(finish_.get(), state.size(), dt, result_.get());
        check(cudaGetLastError(), "finish position Verlet");
        download(result, result_, "position Verlet finished positions and work");
    }

  private:
    unsigned grid(std::size_t count, double dt) const
    {
        if (!std::isfinite(dt) || dt < 0)
            throw std::invalid_argument("Invalid position Verlet time step");
        const auto blocks = count ? (count - 1) / threads_ + 1 : 0;
        if (blocks > max_blocks_)
            throw std::length_error("Position Verlet batch exceeds the device grid");
        return static_cast<unsigned>(blocks);
    }

    static constexpr unsigned threads_ = 128;
    std::size_t max_blocks_ = 0;
    DeviceBuffer<PositionVerletStart> start_;
    DeviceBuffer<PositionVerletFinish> finish_;
    DeviceBuffer<ParticleVector3> position_;
    DeviceBuffer<PositionVerletResult> result_;
};

class CudaSolidTimeStep final : public SolidTimeStepBatch
{
  public:
    explicit CudaSolidTimeStep(const std::vector<double> &mass) : count_(mass.size())
    {
        for (const double value : mass)
            if (!std::isfinite(value) || value <= 0)
                throw std::invalid_argument("Solid explicit time stepping requires positive finite particle masses");
        int device = 0;
        int maximum = 0;
        check(cudaGetDevice(&device), "select solid time stepping device");
        check(cudaDeviceGetAttribute(&maximum, cudaDevAttrMaxGridDimX, device), "read solid stepping grid limit");
        const auto blocks = count_ == 0 ? 0 : (count_ - 1) / threads + 1;
        if (blocks > static_cast<std::size_t>(maximum))
            throw std::length_error("Solid time stepping batch exceeds the device grid");
        blocks_ = static_cast<unsigned>(blocks);
        upload(mass_, mass, "solid time stepping masses");
    }

    void predict(const std::vector<SolidDriftInput> &input, double dt, double damping_factor,
                 SolidKickInput drive, std::vector<SolidDriftResult> &output) override
    {
        validate(input.size(), dt, damping_factor);
        output.resize(count_);
        if (count_ == 0)
            return;
        upload(drift_input_, input, "solid time stepping drift inputs");
        drift_output_.reserve(count_);
        predictSolidStep<<<blocks_, threads>>>(mass_.get(), drift_input_.get(), count_, dt, damping_factor,
                                               drive, drift_output_.get());
        check(cudaGetLastError(), "predict solid time step");
        download(output, drift_output_, "solid time stepping drift results");
    }

    void finish(const std::vector<SolidFinishInput> &input, double dt, double damping_factor,
                SolidKickInput drive, std::vector<SolidFinishResult> &output) override
    {
        validate(input.size(), dt, damping_factor);
        output.resize(count_);
        if (count_ == 0)
            return;
        upload(finish_input_, input, "solid time stepping final kick inputs");
        finish_output_.reserve(count_);
        finishSolidStep<<<blocks_, threads>>>(mass_.get(), finish_input_.get(), count_, dt, damping_factor,
                                              drive, finish_output_.get());
        check(cudaGetLastError(), "finish solid time step");
        download(output, finish_output_, "solid time stepping final kick results");
    }

    const double *deviceMass() const noexcept override { return mass_.get(); }

    void setCoupledDrive(DeviceParticleVectorView first, SolidKickInput first_units,
                         DeviceParticleVectorView second, SolidKickInput second_units,
                         const std::vector<ParticleVector3> &host_second, DeviceCoupledMaterialForces material,
                         const std::vector<ParticleVector3> &internal) override
    {
        if (!first || first.count > count_ ||
            (second ? second.count != count_ - first.count : host_second.size() != count_ - first.count) ||
            (material ? material.first_adjoint.count != first.count ||
                        material.second_adjoint.count != count_ - first.count ||
                        material.second_stabilizing_force.count != count_ - first.count
                      : internal.size() != count_))
            throw std::invalid_argument("Coupled force fields do not match the integrator");
        if (!second)
        {
            upload(host_external_, host_second, "host contact fallback forces");
            second = {host_external_.get(), host_second.size()};
        }
        if (!material)
            upload(coupled_internal_, internal, "coupled material forces");
        coupled_drive_.reserve(count_);
        assembleCoupledDrive<<<blocks_, threads>>>(mass_.get(), material ? nullptr : coupled_internal_.get(),
                                                     material, first, first_units, second, second_units, coupled_drive_.get());
        check(cudaGetLastError(), "assemble coupled forces");
    }

    void readCoupledDrive(std::vector<ParticleVector3> &force) override
    {
        force.resize(count_);
        download(force, coupled_drive_, "read coupled driving force");
    }

    void configureCoupledReference(const std::vector<ParticleVector3> &reference) override
    {
        if (reference.size() != count_)
            throw std::invalid_argument("Coupled reference geometry does not match the integrator");
        upload(coupled_reference_, reference, "coupled reference geometry");
    }

    DeviceSolidKinematics coupledKinematics(bool trial) override
    {
        if (!coupled_reference_.get() || !coupled_motion_.get() || (trial && !currentTrial()))
            throw std::logic_error("Coupled geometry requested without its motion state");
        geometry_displacement_.reserve(count_);
        geometry_position_.reserve(count_);
        geometry_velocity_.reserve(count_);
        if (trial)
            coupledGeometry<true><<<blocks_, threads>>>(
                coupled_reference_.get(), coupled_motion_.get(), currentTrial(), coupled_drift_.get(),
                count_, geometry_displacement_.get(), geometry_position_.get(), geometry_velocity_.get());
        else
            coupledGeometry<false><<<blocks_, threads>>>(
                coupled_reference_.get(), coupled_motion_.get(), nullptr, nullptr,
                count_, geometry_displacement_.get(), geometry_position_.get(), geometry_velocity_.get());
        check(cudaGetLastError(), "build coupled geometry on the device");
        return {{geometry_displacement_.get(), count_}, {geometry_position_.get(), count_},
                {geometry_velocity_.get(), count_}};
    }

    void initializeCoupledMotion(const std::vector<SolidMotion> &motion) override
    {
        if (motion.size() != count_)
            throw std::invalid_argument("Coupled motion does not match the integrator");
        upload(coupled_motion_, motion, "initial coupled motion");
    }

    void predictCoupled(double dt, double factor, bool implicit, std::vector<SolidDriftResult> &output) override
    {
        validate(count_, dt, factor);
        if (!coupled_drive_.get() || !coupled_motion_.get())
            throw std::logic_error("Coupled prediction needs accepted motion and a driving force");
        trial_first_ = true;
        saved_implicit_force_ = implicit;
        coupled_drift_.reserve(count_);
        trial_a_.reserve(count_);
        if (implicit)
        {
            initial_drive_.reserve(count_);
            check(cudaMemcpyAsync(initial_drive_.get(), coupled_drive_.get(), count_ * sizeof(ParticleVector3),
                                  cudaMemcpyDeviceToDevice), "preserve implicit initial force");
        }
        predictCoupledStep<<<blocks_, threads>>>(mass_.get(), coupled_motion_.get(), coupled_drive_.get(),
                                                   count_, dt, factor, coupled_drift_.get(), trial_a_.get());
        check(cudaGetLastError(), "predict coupled motion");
        output.resize(count_);
        download(output, coupled_drift_, "read trial geometry for material and contact search");
    }

    void correctCoupledTrial(double dt, std::vector<SolidImplicitCorrection> &output) override
    {
        validate(count_, dt, 1);
        if (!saved_implicit_force_)
            throw std::logic_error("Coupled correction needs the force at the start of the trial");
        trial_b_.reserve(count_);
        correction_.reserve(count_);
        auto *candidate = trial_first_ ? trial_b_.get() : trial_a_.get();
        correctCoupledStep<<<blocks_, threads>>>(
            mass_.get(), coupled_motion_.get(), coupled_drift_.get(), initial_drive_.get(), coupled_drive_.get(),
            currentTrial(), count_, dt, correction_.get(), candidate);
        check(cudaGetLastError(), "correct coupled implicit trial");
        output.resize(count_);
        download(output, correction_, "read corrected trial geometry");
    }

    void advanceCoupledTrial() override { trial_first_ = !trial_first_; }

    void readCoupledAcceleration(std::vector<ParticleVector3> &acceleration) override
    {
        acceleration.resize(count_);
        download(acceleration, coupled_acceleration_, "read coupled acceleration for checkpoint");
    }

    void finishCoupled(double dt, double factor, std::vector<SolidCoupledFinishResult> &output) override
    {
        validate(count_, dt, factor);
        coupled_finish_.reserve(count_);
        coupled_acceleration_.reserve(count_);
        finishCoupledStep<<<blocks_, threads>>>(mass_.get(), coupled_motion_.get(), coupled_drift_.get(),
                                                  coupled_drive_.get(), currentTrial(), count_, dt, factor,
                                                  coupled_finish_.get(), coupled_acceleration_.get());
        check(cudaGetLastError(), "finish coupled motion");
        output.resize(count_);
        download(output, coupled_finish_, "read accepted velocities for host contact geometry");
    }

  private:
    const ParticleVector3 *currentTrial() const { return trial_first_ ? trial_a_.get() : trial_b_.get(); }

    void validate(std::size_t count, double dt, double factor) const
    {
        if (count != count_ || !std::isfinite(dt) || dt < 0 || !std::isfinite(factor) || factor < 0)
            throw std::invalid_argument("Solid time stepping inputs do not match the configured state");
    }

    static constexpr unsigned threads = 128;
    std::size_t count_;
    unsigned blocks_ = 0;
    DeviceBuffer<double> mass_;
    DeviceBuffer<SolidDriftInput> drift_input_;
    DeviceBuffer<SolidDriftResult> drift_output_;
    DeviceBuffer<SolidFinishInput> finish_input_;
    DeviceBuffer<SolidFinishResult> finish_output_;
    DeviceBuffer<SolidMotion> coupled_motion_;
    DeviceBuffer<SolidCoupledFinishResult> coupled_finish_;
    DeviceBuffer<ParticleVector3> coupled_acceleration_;
    DeviceBuffer<ParticleVector3> coupled_reference_, geometry_displacement_, geometry_position_, geometry_velocity_;
    DeviceBuffer<SolidDriftResult> coupled_drift_;
    DeviceBuffer<ParticleVector3> coupled_drive_, initial_drive_, coupled_internal_, host_external_, trial_a_, trial_b_;
    DeviceBuffer<SolidImplicitCorrection> correction_;
    bool trial_first_ = true, saved_implicit_force_ = false;
};
}

std::unique_ptr<PositionVerletBatch> makeCudaPositionVerlet()
{
    return std::make_unique<CudaPositionVerlet>();
}

std::unique_ptr<SolidTimeStepBatch> makeCudaSolidTimeStep(const std::vector<double> &mass)
{
    return std::make_unique<CudaSolidTimeStep>(mass);
}
}

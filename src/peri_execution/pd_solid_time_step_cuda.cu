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

__global__ void predictSolidStep(const double *mass, const SolidDriftInput *input, std::size_t count,
                                 double dt, double factor, SolidKickInput drive, SolidDriftResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &old = input[i];
    auto &next = output[i];
    ParticleVector3 damped_velocity;
    for (int axis = 0; axis < 3; ++axis)
    {
        damped_velocity.component[axis] = factor * old.velocity.component[axis];
        double kick = (.5 * dt) * old.drive.component[axis];
        if (drive == SolidKickInput::Force)
            kick /= mass[i];
        next.half_velocity.component[axis] = damped_velocity.component[axis] + kick;
        next.displacement.component[axis] = old.displacement.component[axis] + dt * next.half_velocity.component[axis];
    }
    next.damped_kinetic = .5 * mass[i] * normSquared(damped_velocity);
}

__global__ void finishSolidStep(const double *mass, const SolidFinishInput *input, std::size_t count,
                                double dt, double factor, SolidKickInput drive, SolidFinishResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &trial = input[i];
    auto &next = output[i];
    ParticleVector3 undamped_velocity;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double force = trial.force.component[axis];
        const double additional = trial.additional_force.component[axis];
        next.acceleration.component[axis] = (force + additional) / mass[i];
        // Both forms have the same discrete impulse. Keep the arithmetic
        // order used by the corresponding CPU stepping operator.
        const double kick = drive == SolidKickInput::Force ? (.5 * dt) * force / mass[i]
            : (.5 * dt) * (force / mass[i]);
        undamped_velocity.component[axis] = trial.half_velocity.component[axis] + kick;
        next.velocity.component[axis] = factor * undamped_velocity.component[axis];
    }
    next.undamped_kinetic = .5 * mass[i] * normSquared(undamped_velocity);
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

  private:
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

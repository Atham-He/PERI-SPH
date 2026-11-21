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

std::unique_ptr<SolidTimeStepBatch> makeCudaSolidTimeStep(const std::vector<double> &mass)
{
    return std::make_unique<CudaSolidTimeStep>(mass);
}
}

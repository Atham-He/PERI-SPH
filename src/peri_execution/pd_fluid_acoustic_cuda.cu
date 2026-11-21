#include "pd_fluid_acoustic_batch.h"
#include "pd_cuda_buffer.cuh"

#include <algorithm>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

struct DeviceStencil
{
    DeviceBuffer<std::size_t> offsets;
    DeviceBuffer<FluidAcousticNeighbor> neighbors;

    void uploadRows(const FluidAcousticStencil &stencil)
    {
        upload(offsets, stencil.offsets, "fluid acoustic offsets");
        upload(neighbors, stencil.neighbors, "fluid acoustic neighbors");
    }
};

__device__ double dot(const ParticleVector3 &a, const ParticleVector3 &b)
{
    return a.component[0] * b.component[0] +
        (a.component[1] * b.component[1] + a.component[2] * b.component[2]);
}

__device__ ParticleVector3 difference(const ParticleVector3 &a, const ParticleVector3 &b)
{
    ParticleVector3 result;
    for (int axis = 0; axis < 3; ++axis)
        result.component[axis] = a.component[axis] - b.component[axis];
    return result;
}

__device__ void addForce(ParticleVector3 &force, double scale, const ParticleVector3 &direction)
{
    for (int axis = 0; axis < 3; ++axis)
        force.component[axis] += scale * direction.component[axis];
}

__device__ double dissipativePressure(double jump, FluidAcousticCoefficients coefficients)
{
    // BaseAcousticRiemannSolver with its default TruncatedLinear limiter.
    const double scaled = 3.0 * (coefficients.inverse_sound_speed * (jump > 0.0 ? jump : 0.0));
    return coefficients.impedance * jump * (scaled < 1.0 ? scaled : 1.0);
}

template <bool pressure>
__global__ void evaluateAcoustic(const FluidAcousticSite *sites, const FluidAcousticWall *walls,
                                 const std::size_t *inner_offsets, const FluidAcousticNeighbor *inner,
                                 const std::size_t *wall_offsets, const FluidAcousticNeighbor *wall,
                                 std::size_t count, FluidAcousticCoefficients coefficients,
                                 FluidAcousticResult *result)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &source = sites[i];
    FluidAcousticResult row{};
    for (auto n = inner_offsets[i]; n < inner_offsets[i + 1]; ++n)
    {
        const auto &pair = inner[n];
        const auto &target = sites[pair.target];
        const double gradient = pair.kernel_derivative * target.volume;
        if constexpr (pressure)
        {
            addForce(row.inner_force, -(source.pressure + target.pressure) * gradient, pair.direction);
            row.inner_density_rate += ((source.pressure - target.pressure) * coefficients.inverse_impedance) * gradient;
        }
        else
        {
            const double jump = dot(difference(source.velocity, target.velocity), pair.direction);
            row.inner_density_rate += jump * gradient;
            addForce(row.inner_force, dissipativePressure(jump, coefficients) * gradient, pair.direction);
        }
    }

    for (auto n = wall_offsets[i]; n < wall_offsets[i + 1]; ++n)
    {
        const auto &pair = wall[n];
        const auto &target = walls[pair.target];
        const double gradient = pair.kernel_derivative * target.volume;
        if constexpr (pressure)
        {
            ParticleVector3 acceleration, direction_to_wall;
            for (int axis = 0; axis < 3; ++axis)
            {
                acceleration.component[axis] = source.force_prior.component[axis] / source.mass -
                    target.average_acceleration.component[axis];
                direction_to_wall.component[axis] = -pair.direction.component[axis];
            }
            const double face_acceleration = dot(acceleration, direction_to_wall);
            const double wall_pressure = source.pressure + target.acceleration_factor * source.density * pair.distance *
                (face_acceleration < 0.0 ? 0.0 : face_acceleration);
            addForce(row.wall_force, -(source.pressure + wall_pressure) * gradient, pair.direction);
            row.wall_density_rate += ((source.pressure - wall_pressure) * coefficients.inverse_impedance) * gradient;
        }
        else
        {
            const double facing = dot(pair.direction, target.normal);
            const double sign = (0.0 < facing) - (facing < 0.0);
            ParticleVector3 ghost_velocity, normal;
            for (int axis = 0; axis < 3; ++axis)
            {
                ghost_velocity.component[axis] = 2.0 * target.average_velocity.component[axis] - source.velocity.component[axis];
                normal.component[axis] = sign * target.normal.component[axis];
            }
            row.wall_density_rate += dot(difference(source.velocity, ghost_velocity), pair.direction) * gradient;
            const double jump = 2.0 * dot(difference(source.velocity, target.average_velocity), normal);
            addForce(row.wall_force, dissipativePressure(jump, coefficients) * gradient, normal);
        }
    }

    for (int axis = 0; axis < 3; ++axis)
    {
        row.inner_force.component[axis] *= source.volume;
        row.wall_force.component[axis] *= source.volume;
    }
    row.inner_density_rate *= source.density;
    row.wall_density_rate *= source.density;
    result[i] = row;
}

bool validStencil(const FluidAcousticStencil &stencil, std::size_t count, std::size_t target_count)
{
    return stencil.offsets.size() == count + 1 && stencil.offsets.front() == 0 &&
        stencil.offsets.back() == stencil.neighbors.size() &&
        std::is_sorted(stencil.offsets.begin(), stencil.offsets.end()) &&
        std::all_of(stencil.neighbors.begin(), stencil.neighbors.end(),
                    [=](const auto &neighbor) { return neighbor.target < target_count; });
}

class CudaFluidAcoustic final : public FluidAcousticBatch
{
  public:
    CudaFluidAcoustic()
    {
        int device = 0, maximum = 0;
        check(cudaGetDevice(&device), "select fluid acoustic device");
        check(cudaDeviceGetAttribute(&maximum, cudaDevAttrMaxGridDimX, device), "read fluid acoustic grid limit");
        max_blocks_ = static_cast<std::size_t>(maximum);
    }

    void evaluate(FluidAcousticHalf half, std::size_t count, const std::vector<FluidAcousticSite> &sites,
                  const std::vector<FluidAcousticWall> &walls, const FluidAcousticStencil &inner,
                  const FluidAcousticStencil &wall, FluidAcousticCoefficients coefficients,
                  std::vector<FluidAcousticResult> &result) override
    {
        if (count > sites.size() || !validStencil(inner, count, sites.size()) || !validStencil(wall, count, walls.size()))
            throw std::invalid_argument("Invalid fluid acoustic stencil");
        result.resize(count);
        if (count == 0)
            return;

        constexpr unsigned threads = 128;
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("Fluid acoustic batch exceeds the device grid");
        upload(sites_, sites, "fluid acoustic sites");
        upload(walls_, walls, "fluid acoustic walls");
        inner_.uploadRows(inner);
        wall_.uploadRows(wall);
        result_.reserve(count);

        if (half == FluidAcousticHalf::Pressure)
            evaluateAcoustic<true><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), walls_.get(),
                inner_.offsets.get(), inner_.neighbors.get(), wall_.offsets.get(), wall_.neighbors.get(),
                count, coefficients, result_.get());
        else
            evaluateAcoustic<false><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), walls_.get(),
                inner_.offsets.get(), inner_.neighbors.get(), wall_.offsets.get(), wall_.neighbors.get(),
                count, coefficients, result_.get());

        check(cudaGetLastError(), "evaluate fluid acoustic interactions");
        download(result, result_, "fluid acoustic results");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceStencil inner_, wall_;
    DeviceBuffer<FluidAcousticSite> sites_;
    DeviceBuffer<FluidAcousticWall> walls_;
    DeviceBuffer<FluidAcousticResult> result_;
};
}

std::unique_ptr<FluidAcousticBatch> makeCudaFluidAcoustic()
{
    return std::make_unique<CudaFluidAcoustic>();
}
}

#include "pd_fluid_cfl_batch.h"
#include "pd_solid_cfl_batch.h"
#include "pd_cuda_buffer.cuh"

#include <algorithm>
#include <cub/block/block_reduce.cuh>
#include <math_constants.h>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::upload;

constexpr unsigned block_size = 128;
using BlockReduce = cub::BlockReduce<double, block_size>;

struct Maximum
{
    __device__ double operator()(double lhs, double rhs) const { return lhs > rhs ? lhs : rhs; }
};

__device__ double squaredNorm(const ParticleVector3 &v)
{
    return v.component[0] * v.component[0] +
        (v.component[1] * v.component[1] + v.component[2] * v.component[2]);
}

__device__ double candidate(const FluidAcousticCflSite &site, double)
{
    const double speed = site.sound_speed + sqrt(squaredNorm(site.velocity));
    return isfinite(speed) ? speed : CUDART_INF;
}

__device__ double candidate(const FluidAdvectionCflSite &site, double h)
{
    ParticleVector3 total_force;
    for (int axis = 0; axis < 3; ++axis)
        total_force.component[axis] = site.force.component[axis] + site.force_prior.component[axis];

    const double speed_squared = squaredNorm(site.velocity);
    const double acceleration_scale = 4.0 * h * sqrt(squaredNorm(total_force)) / site.mass;
    return isfinite(speed_squared) && isfinite(acceleration_scale) ?
        Maximum{}(speed_squared, acceleration_scale) : CUDART_INF;
}

struct PDWaveParameters
{
    double spacing, sound_speed;
};

__device__ double candidate(const PDCflSite &site, PDWaveParameters parameters)
{
    const auto wave = (parameters.sound_speed + sqrt(squaredNorm(site.velocity))) / parameters.spacing;
    const auto acceleration = sqrt(sqrt(squaredNorm(site.acceleration)) / parameters.spacing);
    return isfinite(wave) && isfinite(acceleration) ? Maximum{}(wave, acceleration) : CUDART_INF;
}

template <class Site, class Parameters>
__global__ void reduceSiteSpeeds(const Site *sites, std::size_t count, double initial, Parameters parameters,
                                 double *block_values)
{
    std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto stride = std::size_t(gridDim.x) * blockDim.x;
    double maximum = initial;
    for (; i < count; i += stride)
        maximum = Maximum{}(maximum, candidate(sites[i], parameters));

    __shared__ BlockReduce::TempStorage temp;
    const auto block_maximum = BlockReduce(temp).Reduce(maximum, Maximum{});
    if (threadIdx.x == 0)
        block_values[blockIdx.x] = block_maximum;
}

__global__ void reduceBlockSpeeds(const double *blocks, unsigned count, double initial, double *result)
{
    double maximum = initial;
    for (unsigned i = threadIdx.x; i < count; i += blockDim.x)
        maximum = Maximum{}(maximum, blocks[i]);

    __shared__ BlockReduce::TempStorage temp;
    const auto aggregate = BlockReduce(temp).Reduce(maximum, Maximum{});
    if (threadIdx.x == 0)
        *result = aggregate;
}

class CudaCflReduction
{
  public:
    template <class Site, class Parameters>
    double evaluate(const std::vector<Site> &sites, DeviceBuffer<Site> &device,
                    double initial, Parameters parameters)
    {
        if (sites.empty())
            return initial;

        // Grid-stride rows keep the final reduction to one block, even for
        // large bodies. This cap limits only temporary storage.
        constexpr std::size_t max_blocks = 512;
        const auto blocks = static_cast<unsigned>(std::min(max_blocks, 1 + (sites.size() - 1) / block_size));
        upload(device, sites, "CFL inputs");
        partial_.reserve(blocks);

        reduceSiteSpeeds<<<blocks, block_size>>>(device.get(), sites.size(), initial, parameters, partial_.get());
        check(cudaGetLastError(), "reduce CFL per-particle bounds");
        const double *result = partial_.get();
        if (blocks > 1)
        {
            aggregate_.reserve(1);
            reduceBlockSpeeds<<<1, block_size>>>(partial_.get(), blocks, initial, aggregate_.get());
            check(cudaGetLastError(), "reduce CFL block bounds");
            result = aggregate_.get();
        }
        double maximum = initial;
        check(cudaMemcpy(&maximum, result, sizeof(maximum), cudaMemcpyDeviceToHost), "read CFL maximum");
        return maximum;
    }

  private:
    DeviceBuffer<double> partial_, aggregate_;
};

class CudaFluidCfl final : public FluidCflBatch
{
  public:
    double acousticMaximum(const std::vector<FluidAcousticCflSite> &sites, double initial) override
    {
        return reduction_.evaluate(sites, acoustic_sites_, initial, 0.0);
    }

    double advectionMaximum(const std::vector<FluidAdvectionCflSite> &sites, double initial,
                            double smoothing_length) override
    {
        return reduction_.evaluate(sites, advection_sites_, initial, smoothing_length);
    }

  private:
    CudaCflReduction reduction_;
    DeviceBuffer<FluidAcousticCflSite> acoustic_sites_;
    DeviceBuffer<FluidAdvectionCflSite> advection_sites_;
};

class CudaPDCfl final : public PDCflBatch
{
  public:
    double maximumInverseStep(const std::vector<PDCflSite> &sites, double initial,
                              double spacing, double sound_speed) override
    {
        return reduction_.evaluate(sites, sites_, initial, PDWaveParameters{spacing, sound_speed});
    }

  private:
    CudaCflReduction reduction_;
    DeviceBuffer<PDCflSite> sites_;
};
}

std::unique_ptr<FluidCflBatch> makeCudaFluidCfl()
{
    return std::make_unique<CudaFluidCfl>();
}
std::unique_ptr<PDCflBatch> makeCudaPDCfl()
{
    return std::make_unique<CudaPDCfl>();
}
}

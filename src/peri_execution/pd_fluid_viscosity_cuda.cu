#include "pd_fluid_viscosity_batch.h"
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

__device__ ParticleVector3 sumViscosity(const FluidViscousSite &source, const FluidViscousSite *sites,
                                        const FluidViscousPair *pairs, std::size_t begin, std::size_t end)
{
    ParticleVector3 sum{};
    for (auto n = begin; n < end; ++n)
    {
        const auto &pair = pairs[n];
        const auto &target = sites[pair.target];
        for (int axis = 0; axis < 3; ++axis)
        {
            const double difference = source.velocity.component[axis] - target.velocity.component[axis];
            const double derivative = (pair.velocity_factor * difference) / pair.distance_with_offset;
            sum.component[axis] += pair.force_factor * derivative * pair.kernel_derivative * target.volume;
        }
    }
    return sum;
}

__global__ void evaluateViscosity(const FluidViscousSite *sites, const std::size_t *offsets,
                                  const std::size_t *inner_ends, const FluidViscousPair *pairs,
                                  std::size_t count, ParticleVector3 *force)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &source = sites[i];
    const auto inner = sumViscosity(source, sites, pairs, offsets[i], inner_ends[i]);
    const auto contact = sumViscosity(source, sites, pairs, inner_ends[i], offsets[i + 1]);
    for (int axis = 0; axis < 3; ++axis)
        force[i].component[axis] = inner.component[axis] * source.volume + contact.component[axis] * source.volume;
}

class CudaFluidViscosity final : public FluidViscosityBatch
{
  public:
    CudaFluidViscosity()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select fluid viscosity device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read fluid viscosity grid limit");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void evaluate(const FluidViscousRows &rows, std::vector<ParticleVector3> &force) override
    {
        const auto count = rows.inner_ends.size();
        if (rows.offsets.size() != count + 1 || rows.offsets.front() != 0 ||
            rows.offsets.back() != rows.pairs.size() || count > rows.sites.size() ||
            std::any_of(rows.pairs.begin(), rows.pairs.end(),
                        [&](const auto &pair) { return pair.target >= rows.sites.size(); }))
            throw std::invalid_argument("Invalid fluid viscosity stencil");
        for (std::size_t i = 0; i < count; ++i)
            if (rows.offsets[i] > rows.inner_ends[i] || rows.inner_ends[i] > rows.offsets[i + 1])
                throw std::invalid_argument("Invalid inner/contact split in fluid viscosity stencil");

        force.resize(count);
        if (count == 0)
            return;

        constexpr unsigned threads = 128;
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("Fluid viscosity batch exceeds the device grid");
        upload(sites_, rows.sites, "fluid viscosity sites");
        upload(offsets_, rows.offsets, "fluid viscosity offsets");
        upload(inner_ends_, rows.inner_ends, "fluid viscosity inner offsets");
        upload(pairs_, rows.pairs, "fluid viscosity pairs");
        force_.reserve(count);
        evaluateViscosity<<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), offsets_.get(),
            inner_ends_.get(), pairs_.get(), count, force_.get());
        check(cudaGetLastError(), "evaluate fluid viscosity");
        download(force, force_, "fluid viscosity forces");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<FluidViscousSite> sites_;
    DeviceBuffer<std::size_t> offsets_, inner_ends_;
    DeviceBuffer<FluidViscousPair> pairs_;
    DeviceBuffer<ParticleVector3> force_;
};
}

std::unique_ptr<FluidViscosityBatch> makeCudaFluidViscosity()
{
    return std::make_unique<CudaFluidViscosity>();
}
}

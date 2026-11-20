#include "pd_spatial_derivatives.h"
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

template <bool second_order>
__global__ void evaluateSpatialDerivatives(const std::size_t *offsets, const std::size_t *neighbors,
                                           const std::size_t *origins, const ParticleVector3 *gradients,
                                           const ShapeHessian3 *hessians, const ParticleVector3 *particles,
                                           bool add_identity, std::size_t rows, SpatialDerivatives3 *result)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto row = index / 3;
    const auto axis = index % 3;
    if (row >= rows)
        return;

    double gradient[3]{}, hessian[3][3]{};
    if (add_identity)
        gradient[axis] = 1.0;
    const double origin = origins ? particles[origins[row]].component[axis] : 0.0;
    for (auto i = offsets[row]; i < offsets[row + 1]; ++i)
    {
        const double value = particles[neighbors[i]].component[axis] - origin;
        for (int j = 0; j < 3; ++j)
        {
            gradient[j] += value * gradients[i].component[j];
            if constexpr (second_order)
                for (int k = 0; k < 3; ++k)
                    hessian[j][k] += value * hessians[i].derivative[j][k];
        }
    }
    for (int j = 0; j < 3; ++j)
    {
        result[row].gradient[axis][j] = gradient[j];
        for (int k = 0; k < 3; ++k)
            result[row].hessian[axis][j][k] = hessian[j][k];
    }
}

class CudaSpatialDerivatives final : public SpatialDerivativeBatch
{
  public:
    CudaSpatialDerivatives()
    {
        int device = 0;
        check(cudaGetDevice(&device), "select spatial derivative device");
        cudaDeviceProp properties{};
        check(cudaGetDeviceProperties(&properties, device), "get spatial derivative device properties");
        max_blocks_ = static_cast<std::size_t>(properties.maxGridSize[0]);
    }

    void configure(const SpatialDerivativeMap &map) override
    {
        const auto entries = map.neighbors.size();
        if (map.offsets.empty() || map.offsets.front() != 0 || map.offsets.back() != entries ||
            !std::is_sorted(map.offsets.begin(), map.offsets.end()) || map.gradients.size() != entries ||
            (!map.hessians.empty() && map.hessians.size() != entries))
            throw std::invalid_argument("Invalid spatial derivative stencil");
        const auto rows = map.offsets.size() - 1;
        if (!map.origins.empty() && map.origins.size() != rows)
            throw std::invalid_argument("Invalid spatial derivative origins");
        for (auto node : map.neighbors)
            if (node >= map.particle_count)
                throw std::invalid_argument("Spatial derivative neighbor is outside the particle field");
        for (auto node : map.origins)
            if (node >= map.particle_count)
                throw std::invalid_argument("Spatial derivative origin is outside the particle field");

        configured_ = false;
        upload(offsets_, map.offsets, "spatial derivative offsets");
        upload(neighbors_, map.neighbors, "spatial derivative neighbors");
        upload(origins_, map.origins, "spatial derivative origins");
        upload(gradients_, map.gradients, "reference gradients");
        upload(hessians_, map.hessians, "reference Hessians");
        particle_count_ = map.particle_count;
        rows_ = rows;
        use_origins_ = !map.origins.empty();
        use_hessians_ = !map.hessians.empty();
        add_identity_ = map.add_identity;
        configured_ = true;
    }

    void evaluate(const std::vector<ParticleVector3> &particles,
                  std::vector<SpatialDerivatives3> &result) override
    {
        if (!configured_ || particles.size() != particle_count_)
            throw std::invalid_argument("Particle field does not match the spatial derivative stencil");
        result.resize(rows_);
        if (rows_ == 0)
            return;

        constexpr unsigned int threads = 128;
        const std::size_t blocks = (rows_ * 3 + threads - 1) / threads;
        if (blocks > max_blocks_)
            throw std::length_error("Spatial derivative batch exceeds the device grid");
        upload(particles_, particles, "spatial derivative particle field");
        result_.reserve(rows_);
        const auto *origins = use_origins_ ? origins_.get() : nullptr;
        if (use_hessians_)
            evaluateSpatialDerivatives<true><<<static_cast<unsigned int>(blocks), threads>>>(
                offsets_.get(), neighbors_.get(), origins, gradients_.get(), hessians_.get(),
                particles_.get(), add_identity_, rows_, result_.get());
        else
            evaluateSpatialDerivatives<false><<<static_cast<unsigned int>(blocks), threads>>>(
                offsets_.get(), neighbors_.get(), origins, gradients_.get(), nullptr,
                particles_.get(), add_identity_, rows_, result_.get());
        check(cudaGetLastError(), "spatial derivative evaluation");
        download(result, result_, "spatial derivative results");
    }

  private:
    std::size_t particle_count_ = 0, rows_ = 0, max_blocks_ = 0;
    bool use_origins_ = false, use_hessians_ = false, add_identity_ = false, configured_ = false;
    DeviceBuffer<std::size_t> offsets_, neighbors_, origins_;
    DeviceBuffer<ParticleVector3> gradients_, particles_;
    DeviceBuffer<ShapeHessian3> hessians_;
    DeviceBuffer<SpatialDerivatives3> result_;
};
}

std::unique_ptr<SpatialDerivativeBatch> makeCudaSpatialDerivatives()
{
    return std::make_unique<CudaSpatialDerivatives>();
}
}

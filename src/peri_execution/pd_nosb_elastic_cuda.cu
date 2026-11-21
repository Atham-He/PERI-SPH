#include "pd_nosb_elastic_batch.h"
#include "pd_hughes_winget.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

template <int dimensions>
using Tensor = Eigen::Matrix<double, dimensions, dimensions>;

template <int dimensions>
__device__ Tensor<dimensions> load(const ShapeHessian3 &entries)
{
    Tensor<dimensions> tensor;
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            tensor(i, j) = entries.derivative[i][j];
    return tensor;
}

template <int dimensions>
__device__ ShapeHessian3 pack(const Tensor<dimensions> &tensor)
{
    ShapeHessian3 entries{};
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            entries.derivative[i][j] = tensor(i, j);
    return entries;
}

template <int dimensions>
__device__ bool allFinite(const Tensor<dimensions> &tensor)
{
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            if (!isfinite(tensor(i, j)))
                return false;
    return true;
}

template <int dimensions>
__global__ void updateNosbPDElasticity(const NosbPDElasticSite *sites, std::size_t count,
                                       double lambda, double shear, NosbPDElasticResult *results)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &site = sites[i];
    NosbPDElasticResult result{};
    if (!site.active)
    {
        result.status = NosbPDElasticStatus::Inactive;
        results[i] = result;
        return;
    }

    const auto shape_inverse = load<dimensions>(site.shape_inverse);
    const Tensor<dimensions> full = load<dimensions>(site.full) * shape_inverse;
    const Tensor<dimensions> increment = load<dimensions>(site.increment) * shape_inverse;
    const Tensor<dimensions> midpoint = load<dimensions>(site.midpoint) * shape_inverse;
    result.deformation = pack(full);
    result.deformation_delta = pack(increment);
    result.deformation_half = pack(midpoint);

    if (!allFinite(full) || !allFinite(midpoint))
    {
        result.status = NosbPDElasticStatus::NonFiniteGradient;
        results[i] = result;
        return;
    }
    const double determinant = full.determinant(), midpoint_determinant = midpoint.determinant();
    if (!isfinite(determinant) || !isfinite(midpoint_determinant) ||
        determinant < nosb_minimum_deformation_determinant || midpoint_determinant < nosb_minimum_deformation_determinant)
    {
        result.status = NosbPDElasticStatus::InvalidDeterminant;
        results[i] = result;
        return;
    }

    const Tensor<dimensions> inverse = full.inverse(), inverse_half = midpoint.inverse();
    const Tensor<dimensions> gradient = increment * inverse_half;
    const Tensor<dimensions> stress = hughesWingetStress(gradient, load<dimensions>(site.stress), lambda, shear);
    const Tensor<dimensions> piola = determinant * inverse * stress;
    const Tensor<dimensions> piola_transpose = piola.transpose();
    const Tensor<dimensions> force_state = piola_transpose * shape_inverse;
    result.inverse = pack(inverse);
    result.inverse_half = pack(inverse_half);
    result.stress = pack(stress);
    result.piola = pack(piola_transpose);
    result.force_state = pack(force_state);
    result.status = NosbPDElasticStatus::Accepted;
    results[i] = result;
}

class CudaNosbPDElastic final : public NosbPDElasticBatch
{
  public:
    CudaNosbPDElastic()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select NOSB-PD elastic device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read NOSB-PD elastic grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void update(const std::vector<NosbPDElasticSite> &sites, double lambda, double shear, int dimensions,
                std::vector<NosbPDElasticResult> &results) override
    {
        if (dimensions != 2 && dimensions != 3)
            throw std::invalid_argument("NOSB-PD elasticity requires two or three dimensions");
        results.resize(sites.size());
        if (sites.empty())
            return;

        constexpr unsigned threads = 128;
        const auto count = sites.size();
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("NOSB-PD elastic update exceeds the device grid");
        upload(sites_, sites, "NOSB-PD elastic sites");
        results_.reserve(count);
        if (dimensions == 2)
            updateNosbPDElasticity<2><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, lambda, shear, results_.get());
        else
            updateNosbPDElasticity<3><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, lambda, shear, results_.get());
        check(cudaGetLastError(), "update NOSB-PD elasticity");
        download(results, results_, "NOSB-PD elastic results");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<NosbPDElasticSite> sites_;
    DeviceBuffer<NosbPDElasticResult> results_;
};
}

std::unique_ptr<NosbPDElasticBatch> makeCudaNosbPDElastic()
{
    return std::make_unique<CudaNosbPDElastic>();
}
}

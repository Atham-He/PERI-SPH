#include "drucker_prager_cap_kernel.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::dpc
{
namespace
{
using peridynamics::cuda_detail::check;
using peridynamics::cuda_detail::DeviceBuffer;
using peridynamics::cuda_detail::download;
using peridynamics::cuda_detail::upload;

__global__ void integrateMaterialPoints(PointInputModel model, const PointInput *input,
                                        std::size_t count, PointOutput *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &point = input[i];
    State history;
    history.kappa = point.history.kappa;
    history.alpha = point.history.alpha;
    Tensor deformation;
    for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
        {
            deformation(row, col) = point.deformation[row][col];
            history.inverse_Cp(row, col) = point.history.inverse_Cp[row][col];
        }

    Result material;
    const auto failure = PointModel(model).integrate(deformation, history, material);
    auto &result = output[i];
    result.failure = failure;
    if (failure != PointFailure::None && failure != PointFailure::BadReturn)
        return;

    for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
        {
            result.sigma[row][col] = material.sigma(row, col);
            result.history.inverse_Cp[row][col] = material.state.inverse_Cp(row, col);
        }
    result.history.kappa = material.state.kappa;
    result.history.alpha = material.state.alpha;
    result.p = material.p;
    result.q = material.q;
    result.relative_density = material.relative_density;
    result.plastic_shear = material.plastic_shear;
    result.plastic_stress_work = material.plastic_stress_work;
    result.elastic_energy = material.elastic_energy;
    result.yield_residual = material.yield_residual;
    result.effective_d = material.effective_d;
    result.surface = material.surface;
}

class CudaPointBatch final : public PointBatch
{
  public:
    explicit CudaPointBatch(const Model &model) : model_(model.pointInput())
    {
        int device = 0;
        int max_blocks = 0;
        check(cudaGetDevice(&device), "select DPC material device");
        check(cudaDeviceGetAttribute(&max_blocks, cudaDevAttrMaxGridDimX, device), "read DPC grid limit");
        max_blocks_ = static_cast<std::size_t>(max_blocks);
        upload(nodes_, model.nodes(), "DPC material table");
        model_.nodes = nodes_.get();
    }

    void evaluate(const std::vector<PointInput> &input, std::vector<PointOutput> &output) override
    {
        output.resize(input.size());
        if (input.empty())
            return;
        constexpr unsigned threads = 128;
        const auto blocks = (input.size() - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("DPC material batch exceeds the device grid");
        upload(inputs_, input, "DPC material inputs");
        outputs_.reserve(input.size());
        integrateMaterialPoints<<<static_cast<unsigned>(blocks), threads>>>(model_,
            inputs_.get(), input.size(), outputs_.get());
        check(cudaGetLastError(), "integrate DPC material points");
        download(output, outputs_, "DPC material outputs");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<Node> nodes_;
    DeviceBuffer<PointInput> inputs_;
    DeviceBuffer<PointOutput> outputs_;
    PointInputModel model_;
};
}

std::unique_ptr<PointBatch> makeCudaPointBatch(const Model &model)
{
    return std::make_unique<CudaPointBatch>(model);
}
}

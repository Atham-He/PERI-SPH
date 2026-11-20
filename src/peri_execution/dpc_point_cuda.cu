#include "dpc_snsni_adapter.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::dpc
{
namespace
{
using peridynamics::cuda_detail::check;
using peridynamics::cuda_detail::DeviceBuffer;
using peridynamics::cuda_detail::download;
using peridynamics::cuda_detail::upload;

__device__ void loadPoint(const PointInput &point, Tensor &deformation, State &history)
{
    history.kappa = point.history.kappa;
    history.alpha = point.history.alpha;
    for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
        {
            deformation(row, col) = point.deformation[row][col];
            history.inverse_Cp(row, col) = point.history.inverse_Cp[row][col];
        }

}

__device__ void saveResult(const Result &material, PointFailure failure, PointOutput &result)
{
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

__global__ void integrateMaterialPoints(PointInputModel model, const PointInput *input,
                                        std::size_t count, PointOutput *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    State history;
    Tensor deformation;
    loadPoint(input[i], deformation, history);
    Result material;
    const auto failure = PointModel(model).integrate(deformation, history, material);
    saveResult(material, failure, output[i]);
}

__global__ void integrateSnsniPoints(PointInputModel model, const SnsniPointInput *input,
                                     std::size_t count, double rho0, double cfl_length,
                                     SnsniPointOutput *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &point=input[i];
    auto &result=output[i];
    State history;
    Tensor deformation;
    loadPoint(point.material, deformation, history);
    if (!(deformation.determinant()>0))
    {
        result.material.failure=PointFailure::Deformation;
        return;
    }

    TensorGradient gradient_F, gradient_Cp;
    Eigen::Vector3d gradient_kappa;
    for (int direction=0;direction<3;++direction)
    {
        gradient_kappa[direction]=point.kappa_gradient[direction];
        for (int row=0;row<3;++row)
            for (int col=0;col<3;++col)
            {
                gradient_F[direction](row,col)=point.deformation_gradient[direction][row][col];
                gradient_Cp[direction](row,col)=point.inverse_Cp_gradient[direction][row][col];
            }
    }

    Spectrum elastic_spectrum;
    computeElasticSpectrum(deformation, history, elastic_spectrum);
    const PointModel point_model(model);
    TrialField trial;
    auto failure=evaluateHenckyTrial(point_model,deformation,history,gradient_F,gradient_Cp,gradient_kappa,
                                     trial,&elastic_spectrum);
    if (failure!=PointFailure::None)
    {
        result.material.failure=failure;
        return;
    }

    Result material;
    failure=point_model.integrate(deformation,history,material,&elastic_spectrum);
    saveResult(material,failure,result.material);
    if (failure!=PointFailure::None)
        return;

    Parameters parameters;
    failure=point_model.parameters(material.relative_density, material.state.alpha, parameters);
    if (failure!=PointFailure::None)
    {
        result.material.failure=failure;
        return;
    }
    for (int row=0;row<3;++row)
        for (int col=0;col<3;++col)
        {
            result.trial_sigma[row][col]=trial.sigma(row,col);
            for (int direction=0;direction<3;++direction)
                result.trial_sigma_gradient[direction][row][col]=trial.gradient[direction](row,col);
        }
    result.acoustic_dt=acousticStep(deformation,deformation.inverse(),parameters,rho0,cfl_length);
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
        const unsigned blocks = blocksFor(input.size());
        upload(inputs_, input, "DPC material inputs");
        outputs_.reserve(input.size());
        integrateMaterialPoints<<<blocks, threads>>>(model_,
            inputs_.get(), input.size(), outputs_.get());
        check(cudaGetLastError(), "integrate DPC material points");
        download(output, outputs_, "DPC material outputs");
    }

    void evaluateSnsni(const std::vector<SnsniPointInput> &input, double rho0, double cfl_length,
                        std::vector<SnsniPointOutput> &output) override
    {
        output.resize(input.size());
        if (input.empty())
            return;
        const unsigned blocks=blocksFor(input.size());
        upload(snsni_inputs_,input,"S-NSNI material inputs");
        snsni_outputs_.reserve(input.size());
        integrateSnsniPoints<<<blocks,threads>>>(model_,snsni_inputs_.get(),input.size(),rho0,cfl_length,snsni_outputs_.get());
        check(cudaGetLastError(),"integrate S-NSNI material points");
        download(output,snsni_outputs_,"S-NSNI material outputs");
    }

  private:
    static constexpr unsigned threads=128;

    unsigned blocksFor(std::size_t count) const
    {
        const auto blocks=(count-1)/threads+1;
        if (blocks>max_blocks_)
            throw std::length_error("DPC material batch exceeds the device grid");
        return static_cast<unsigned>(blocks);
    }

    std::size_t max_blocks_ = 0;
    DeviceBuffer<Node> nodes_;
    DeviceBuffer<PointInput> inputs_;
    DeviceBuffer<PointOutput> outputs_;
    DeviceBuffer<SnsniPointInput> snsni_inputs_;
    DeviceBuffer<SnsniPointOutput> snsni_outputs_;
    PointInputModel model_;
};
}

std::unique_ptr<PointBatch> makeCudaPointBatch(const Model &model)
{
    return std::make_unique<CudaPointBatch>(model);
}
}

#include "pd_finite_volume_biot_cuda.cuh"
#include "pd_cuda_buffer.cuh"

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

__global__ void evaluateBiot(const ShapeHessian3 *input,std::size_t count,double bulk,double shear,
                            FiniteVolumeBiotPointResult *result)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i<count)
        result[i]=cuda_detail::evaluateFiniteVolumeBiotPoint(input[i].derivative,bulk,shear);
}

class CudaFiniteVolumeBiot final : public FiniteVolumeBiotBatch
{
  public:
    CudaFiniteVolumeBiot()
    {
        int device=0;
        check(cudaGetDevice(&device),"select finite-volume Biot device");
        check(cudaDeviceGetAttribute(&max_blocks_,cudaDevAttrMaxGridDimX,device),"read finite-volume Biot grid limit");
    }
    void evaluate(const std::vector<ShapeHessian3> &deformation,double bulk,double shear,
                  std::vector<FiniteVolumeBiotPointResult> &result) override
    {
        const std::size_t count=deformation.size();
        result.resize(count);
        if (count==0) return;
        const auto blocks=(count-1)/threads_+1;
        if (blocks>static_cast<std::size_t>(max_blocks_))
            throw std::length_error("Finite-volume Biot batch exceeds the CUDA grid");
        upload(input_,deformation,"finite-volume Biot deformation");
        result_.reserve(count);
        evaluateBiot<<<static_cast<unsigned>(blocks),threads_>>>(input_.get(),count,bulk,shear,result_.get());
        check(cudaGetLastError(),"evaluate finite-volume Biot material");
        download(result,result_,"finite-volume Biot stresses and energy");
    }

  private:
    static constexpr unsigned threads_=128;
    int max_blocks_=0;
    DeviceBuffer<ShapeHessian3> input_;
    DeviceBuffer<FiniteVolumeBiotPointResult> result_;
};
}

std::unique_ptr<FiniteVolumeBiotBatch> makeCudaFiniteVolumeBiot()
{
    return std::make_unique<CudaFiniteVolumeBiot>();
}
}

#include "pd_hertz_sphere_contact.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

__global__ void evaluateHertzPairs(const HertzSphereSite<double> *particles,const HertzSpherePair *candidates,
                                  std::size_t count,double effective_modulus,HertzSphereContact<double> *result)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i<count)
    {
        const auto pair=candidates[i];
        result[i]=evaluateHertzSpherePair(particles[pair.first],particles[pair.second],effective_modulus);
    }
}

class CudaHertzSpherePairs final : public HertzSpherePairBatch
{
  public:
    CudaHertzSpherePairs()
    {
        int device=0;
        check(cudaGetDevice(&device),"select Hertz sphere contact device");
        check(cudaDeviceGetAttribute(&max_blocks_,cudaDevAttrMaxGridDimX,device),"read Hertz sphere contact grid limit");
    }

    void evaluate(const std::vector<HertzSphereSite<double>> &particles,const std::vector<HertzSpherePair> &candidates,
                  double effective_modulus,std::vector<HertzSphereContact<double>> &result) override
    {
        result.resize(candidates.size());
        if (candidates.empty()) return;
        const auto blocks=(candidates.size()-1)/threads_+1;
        if (blocks>static_cast<std::size_t>(max_blocks_))
            throw std::length_error("Hertz sphere contact batch exceeds the CUDA grid");
        for (const auto &pair:candidates)
            if (pair.first>=particles.size()||pair.second>=particles.size())
                throw std::out_of_range("Hertz sphere pair refers to an unknown grain");
        upload(particles_,particles,"Hertz sphere positions and radii");
        upload(candidates_,candidates,"Hertz sphere candidates");
        result_.reserve(candidates.size());
        evaluateHertzPairs<<<static_cast<unsigned>(blocks),threads_>>>(particles_.get(),candidates_.get(),candidates.size(),
                                                                     effective_modulus,result_.get());
        check(cudaGetLastError(),"evaluate Hertz sphere pairs");
        download(result,result_,"Hertz sphere pair results");
    }

  private:
    static constexpr unsigned threads_=128;
    int max_blocks_=0;
    DeviceBuffer<HertzSphereSite<double>> particles_;
    DeviceBuffer<HertzSpherePair> candidates_;
    DeviceBuffer<HertzSphereContact<double>> result_;
};
}

std::unique_ptr<HertzSpherePairBatch> makeCudaHertzSpherePairs()
{
    return std::make_unique<CudaHertzSpherePairs>();
}
}

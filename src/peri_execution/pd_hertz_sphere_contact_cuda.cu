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

__global__ void evaluateFrictionalPairs(const FrictionalSphereSite<double> *particles,
                                        const FrictionalSphereCandidate<double> *candidates,
                                        std::size_t count,FrictionalSphereParameters<double> parameters,
                                        FrictionalSphereResult<double> *result)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<count)
    {
        const auto &entry=candidates[i];
        result[i]=evaluateFrictionalSpherePair(particles[entry.indices.first],particles[entry.indices.second],
                                               entry.old,parameters);
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
        const auto blocks=grid(candidates.size());
        for (const auto &pair:candidates)
            validatePair(pair,particles.size());
        upload(particles_,particles,"Hertz sphere positions and radii");
        upload(candidates_,candidates,"Hertz sphere candidates");
        result_.reserve(candidates.size());
        evaluateHertzPairs<<<blocks,threads_>>>(particles_.get(),candidates_.get(),candidates.size(),
                                                                     effective_modulus,result_.get());
        check(cudaGetLastError(),"evaluate Hertz sphere pairs");
        download(result,result_,"Hertz sphere pair results");
    }

    bool supportsFriction() const noexcept override {return true;}

    void evaluateFrictional(const std::vector<FrictionalSphereSite<double>> &particles,
                            const std::vector<FrictionalSphereCandidate<double>> &candidates,
                            FrictionalSphereParameters<double> parameters,
                            std::vector<FrictionalSphereResult<double>> &result) override
    {
        result.resize(candidates.size());
        if(candidates.empty()) return;
        const auto blocks=grid(candidates.size());
        for(const auto &pair:candidates)
            validatePair(pair.indices,particles.size());
        upload(friction_sites_,particles,"grain positions, velocities and spin");
        upload(friction_candidates_,candidates,"grain contact indices and committed history");
        friction_result_.reserve(candidates.size());
        evaluateFrictionalPairs<<<blocks,threads_>>>(friction_sites_.get(),friction_candidates_.get(),
                                                      candidates.size(),parameters,friction_result_.get());
        check(cudaGetLastError(),"evaluate frictional sphere pairs");
        download(result,friction_result_,"frictional sphere pair results");
    }

  private:
    unsigned grid(std::size_t count) const
    {
        const auto blocks=(count-1)/threads_+1;
        if(blocks>static_cast<std::size_t>(max_blocks_))
            throw std::length_error("Hertz sphere contact batch exceeds the CUDA grid");
        return static_cast<unsigned>(blocks);
    }

    static void validatePair(HertzSpherePair pair,std::size_t count)
    {
        if(pair.first>=count||pair.second>=count)
            throw std::out_of_range("Hertz sphere pair refers to an unknown grain");
    }

    static constexpr unsigned threads_=128;
    int max_blocks_=0;
    DeviceBuffer<HertzSphereSite<double>> particles_;
    DeviceBuffer<HertzSpherePair> candidates_;
    DeviceBuffer<HertzSphereContact<double>> result_;
    DeviceBuffer<FrictionalSphereSite<double>> friction_sites_;
    DeviceBuffer<FrictionalSphereCandidate<double>> friction_candidates_;
    DeviceBuffer<FrictionalSphereResult<double>> friction_result_;
};
}

std::unique_ptr<HertzSpherePairBatch> makeCudaHertzSpherePairs()
{
    return std::make_unique<CudaHertzSpherePairs>();
}
}

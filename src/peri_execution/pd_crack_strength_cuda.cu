#include "pd_crack_strength.h"
#include "pd_cuda_buffer.cuh"

#include <stdexcept>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

__global__ void evaluateCrackStrength(const CrackStrengthCandidate *candidates, const CrackStrengthPair *pairs,
                                     const CrackStrengthState *nodes, const unsigned char *active,
                                     std::size_t count, double strength, CrackStrengthResult *result)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    CrackStrengthResult maximum{0, true};
    if (!active[i])
    {
        const auto &candidate = candidates[i];
        for (auto k = candidate.begin; k < candidate.end; ++k)
        {
            const auto value = crackStrength(nodes[pairs[k].first], nodes[pairs[k].second],
                                             candidate.normal, strength);
            if (!value.valid)
            {
                maximum.valid = false;
                break;
            }
            if (value.ratio > maximum.ratio)
                maximum.ratio = value.ratio;
        }
    }
    result[i] = maximum;
}

class CudaCrackStrength final : public CrackStrengthBatch
{
  public:
    CudaCrackStrength()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select crack strength device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read crack strength grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void configure(std::size_t node_count, const std::vector<CrackStrengthCandidate> &candidates,
                   const std::vector<CrackStrengthPair> &pairs) override
    {
        std::size_t previous_end = 0;
        for (const auto &candidate : candidates)
        {
            if (candidate.begin != previous_end || candidate.end <= candidate.begin || candidate.end > pairs.size())
                throw std::invalid_argument("Crack strength candidate sample offsets are inconsistent");
            previous_end = candidate.end;
        }
        if (previous_end != pairs.size())
            throw std::invalid_argument("Crack strength candidates do not cover all samples");
        for (const auto &pair : pairs)
            if (pair.first >= node_count || pair.second >= node_count)
                throw std::out_of_range("Crack strength sample index");

        configured_ = false;
        upload(candidates_, candidates, "crack strength candidate bank");
        upload(pairs_, pairs, "crack strength samples");
        count_ = candidates.size();
        node_count_ = node_count;
        configured_ = true;
    }

    void evaluate(const std::vector<CrackStrengthState> &nodes, const std::vector<unsigned char> &active,
                  double strength, std::vector<CrackStrengthResult> &result) override
    {
        if (!configured_ || nodes.size() != node_count_ || active.size() != count_)
            throw std::invalid_argument("Crack strength state does not match its candidate bank");
        result.resize(count_);
        if (count_ == 0)
            return;

        constexpr unsigned threads = 128;
        const auto blocks = (count_ - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("Crack strength bank exceeds the device grid");
        upload(nodes_, nodes, "crack strength nodal state");
        upload(active_, active, "crack strength active mask");
        result_.reserve(count_);
        evaluateCrackStrength<<<static_cast<unsigned>(blocks), threads>>>(
            candidates_.get(), pairs_.get(), nodes_.get(), active_.get(), count_, strength, result_.get());
        check(cudaGetLastError(), "evaluate crack strength");
        download(result, result_, "crack strength ratios");
    }

  private:
    std::size_t node_count_ = 0, count_ = 0, max_blocks_ = 0;
    bool configured_ = false;
    DeviceBuffer<CrackStrengthCandidate> candidates_;
    DeviceBuffer<CrackStrengthPair> pairs_;
    DeviceBuffer<CrackStrengthState> nodes_;
    DeviceBuffer<unsigned char> active_;
    DeviceBuffer<CrackStrengthResult> result_;
};
}

std::unique_ptr<CrackStrengthBatch> makeCudaCrackStrength()
{
    return std::make_unique<CudaCrackStrength>();
}
}

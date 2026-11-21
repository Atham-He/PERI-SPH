#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace SPH::peridynamics
{
struct CrackStrengthPair
{
    std::size_t first, second;
};

struct CrackStrengthCandidate
{
    double normal[3];
    std::size_t begin, end;
};

struct CrackStrengthState
{
    double deformation[3][3], piola[3][3];
};

struct CrackStrengthResult
{
    double ratio;
    bool valid;
};

class CrackStrengthBatch
{
  public:
    virtual ~CrackStrengthBatch() = default;
    virtual void configure(std::size_t node_count, const std::vector<CrackStrengthCandidate> &candidates,
                           const std::vector<CrackStrengthPair> &pairs) = 0;
    virtual void evaluate(const std::vector<CrackStrengthState> &nodes,
                          const std::vector<unsigned char> &active, double strength,
                          std::vector<CrackStrengthResult> &result) = 0;
};

std::unique_ptr<CrackStrengthBatch> makeCudaCrackStrength();
}

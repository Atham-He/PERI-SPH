#pragma once

#include <memory>
#include <vector>

namespace SPH::peridynamics
{
struct CrackProjectionInput
{
    // Four master corners, followed by the slave point.
    double coordinates[5][3];
};

struct CrackProjectionResult
{
    double s, t;
    int status;
};

class CrackProjectionBatch
{
  public:
    virtual ~CrackProjectionBatch() = default;
    virtual void project(const std::vector<CrackProjectionInput> &input,
                         std::vector<CrackProjectionResult> &output) = 0;
};

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection();
}

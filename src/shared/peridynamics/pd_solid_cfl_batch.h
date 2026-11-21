#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
struct PDCflSite
{
    ParticleVector3 velocity, acceleration;
};

class PDCflBatch
{
  public:
    virtual ~PDCflBatch() = default;
    // Non-finite results are handed back to the CPU criterion so its
    // particle-specific error reporting is preserved.
    virtual double maximumInverseStep(const std::vector<PDCflSite> &sites, double initial,
                                      double spacing, double sound_speed) = 0;
};

std::unique_ptr<PDCflBatch> makeCudaPDCfl();
}

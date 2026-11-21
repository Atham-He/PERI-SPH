#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
struct FluidAcousticCflSite
{
    ParticleVector3 velocity;
    double sound_speed;
};

struct FluidAdvectionCflSite
{
    ParticleVector3 velocity, force, force_prior;
    double mass;
};

class FluidCflBatch
{
  public:
    virtual ~FluidCflBatch() = default;

    // Non-finite input is reported through a non-finite result so the
    // integration can retain its ordinary CPU diagnostic path.
    virtual double acousticMaximum(const std::vector<FluidAcousticCflSite> &sites, double initial) = 0;
    virtual double advectionMaximum(const std::vector<FluidAdvectionCflSite> &sites, double initial,
                                    double smoothing_length) = 0;
};

std::unique_ptr<FluidCflBatch> makeCudaFluidCfl();
}

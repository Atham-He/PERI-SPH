#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
struct FluidViscousSite
{
    ParticleVector3 velocity;
    double volume;
};

struct FluidViscousPair
{
    std::size_t target;
    double distance_with_offset, velocity_factor, force_factor, kernel_derivative;
};

struct FluidViscousRows
{
    std::vector<std::size_t> offsets, inner_ends;
    std::vector<FluidViscousPair> pairs;
};

class FluidViscosityStencil
{
  public:
    virtual ~FluidViscosityStencil() = default;
    virtual void configure(const FluidViscousRows &rows, std::size_t site_count) = 0;
};

class FluidViscosityBatch
{
  public:
    virtual ~FluidViscosityBatch() = default;
    virtual std::unique_ptr<FluidViscosityStencil> makeStencil() = 0;
    // Source sites occupy the first entries. The remaining sites supply
    // referenced ghosts and neighbors on the opposite body.
    virtual void evaluate(const FluidViscosityStencil &stencil, const std::vector<FluidViscousSite> &sites,
                          std::vector<ParticleVector3> &force) = 0;
};

std::unique_ptr<FluidViscosityBatch> makeCudaFluidViscosity();
}

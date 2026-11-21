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
    // Source sites occupy the first offsets.size() - 1 entries. Other entries
    // hold referenced ghost particles and sites on the opposite body.
    std::vector<FluidViscousSite> sites;
    std::vector<std::size_t> offsets, inner_ends;
    std::vector<FluidViscousPair> pairs;
};

class FluidViscosityBatch
{
  public:
    virtual ~FluidViscosityBatch() = default;
    virtual void evaluate(const FluidViscousRows &rows, std::vector<ParticleVector3> &force) = 0;
};

std::unique_ptr<FluidViscosityBatch> makeCudaFluidViscosity();
}

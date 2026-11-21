#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
struct FluidTractionSolid
{
    ParticleVector3 normal, average_velocity, average_acceleration;
    double volume;
};

struct FluidTractionFluid
{
    ParticleVector3 velocity, force_prior;
    double mass, density, pressure, divergence, volume;
    double acoustic_impedance, inverse_sound_speed;
};

struct FluidTractionNeighbor
{
    std::size_t fluid;
    ParticleVector3 direction;
    double distance, kernel_value, kernel_derivative;
};

struct FluidTractionResult
{
    ParticleVector3 pressure_force, acoustic_force;
    double pressure, kernel_weight;
};

class FluidTractionBatch
{
  public:
    virtual ~FluidTractionBatch() = default;
    // Each row is one solid particle. Neighbor order is shared with the CPU
    // relation and the referenced fluid state is shared by all incident rows.
    virtual void configure(std::size_t solid_sites, std::size_t fluid_sites,
                           const std::vector<std::size_t> &offsets,
                           const std::vector<FluidTractionNeighbor> &neighbors) = 0;
    virtual void evaluate(const std::vector<FluidTractionSolid> &solid,
                          const std::vector<FluidTractionFluid> &fluid,
                          double acceleration_factor,
                          std::vector<FluidTractionResult> &result) = 0;
};

std::unique_ptr<FluidTractionBatch> makeCudaFluidTraction();
}

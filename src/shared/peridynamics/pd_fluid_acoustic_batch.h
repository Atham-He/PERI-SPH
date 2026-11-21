#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
enum class FluidAcousticStage { Pressure, Density, Divergence, Damping };

struct FluidAcousticCoefficients
{
    double inverse_impedance, impedance, inverse_sound_speed;
};

struct FluidAcousticSite
{
    ParticleVector3 velocity, force_prior;
    double volume, density, pressure, mass, divergence;
};

struct FluidAcousticWall
{
    ParticleVector3 average_velocity, average_acceleration, normal;
    double volume, acceleration_factor;
};

struct FluidAcousticNeighbor
{
    std::size_t target;
    ParticleVector3 direction;
    double distance, kernel_derivative;
};

struct FluidAcousticStencil
{
    std::vector<std::size_t> offsets;
    std::vector<FluidAcousticNeighbor> neighbors;
};

// Pressure and density use scaled increments. Divergence returns the two
// unscaled rate sums, and damping returns the two unscaled vector sums.
struct FluidAcousticResult
{
    ParticleVector3 inner_force, wall_force;
    double inner_density_rate, wall_density_rate;
};

class FluidAcousticBatch
{
  public:
    virtual ~FluidAcousticBatch() = default;
    // The first count fluid sites are source rows; the rest are referenced
    // ghosts. The neighbor geometry stays fixed until the next configure call.
    virtual void configure(std::size_t count, std::size_t fluid_sites, std::size_t wall_sites,
                           const FluidAcousticStencil &inner, const FluidAcousticStencil &wall) = 0;
    virtual void evaluate(FluidAcousticStage stage, const std::vector<FluidAcousticSite> &sites,
                          const std::vector<FluidAcousticWall> &walls, FluidAcousticCoefficients coefficients,
                          std::vector<FluidAcousticResult> &result) = 0;
};

std::unique_ptr<FluidAcousticBatch> makeCudaFluidAcoustic();
}

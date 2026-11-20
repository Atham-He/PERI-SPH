#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace SPH::peridynamics
{
struct ParticleVector3
{
    double component[3];
};

struct ShapeHessian3
{
    double derivative[3][3];
};

// Reference gradients and optional Hessians are stored in neighbor order.
// An origin per row requests differences, as in a live PD support.
// Weighted least-squares fit of quadratic differences on each stencil row.
struct WeightedDifferenceEntry
{
    double basis[9], coefficient[9];
    double energy_weight, force_weight;
};

struct SpatialDerivativeMap
{
    std::size_t particle_count;
    std::vector<std::size_t> offsets, neighbors, origins;
    std::vector<ParticleVector3> gradients;
    std::vector<ShapeHessian3> hessians;
    std::vector<WeightedDifferenceEntry> weighted_differences;
    bool add_identity = false;
};

struct SpatialDerivatives3
{
    double gradient[3][3];
    double hessian[3][3][3];
};

struct WeightedDifferenceResult
{
    std::vector<double> energy_density;
    std::vector<ParticleVector3> forces;
};

class SpatialDerivativeBatch
{
  public:
    virtual ~SpatialDerivativeBatch() = default;
    virtual void configure(const SpatialDerivativeMap &map) = 0;
    virtual void evaluate(const std::vector<ParticleVector3> &particles,
                          std::vector<SpatialDerivatives3> &result,
                          WeightedDifferenceResult *weighted_difference = nullptr) = 0;
    // Transpose of the linear part of evaluate. The constant identity has
    // zero derivative. A difference stencil contributes with both signs.
    virtual void assembleAdjoint(const std::vector<SpatialDerivatives3> &coefficients,
                                 std::vector<ParticleVector3> &result) = 0;
};

std::unique_ptr<SpatialDerivativeBatch> makeCudaSpatialDerivatives();
}

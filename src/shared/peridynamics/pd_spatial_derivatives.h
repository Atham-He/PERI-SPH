#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
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

// Weighted least-squares fit of quadratic differences on each stencil row.
struct WeightedDifferenceEntry
{
    double basis[9], coefficient[9];
    double energy_weight, force_weight;
};

struct SmoothedQuadratureSite
{
    ShapeHessian3 moment;
    double volume;
};

struct SmoothedStressField
{
    ShapeHessian3 stress, trial_gradient[3], correction_gradient[3];
};

// Reference gradients and optional Hessians are stored in neighbor order.
// An origin per row requests differences, as in a live PD support.
struct SpatialDerivativeMap
{
    std::size_t particle_count;
    std::vector<std::size_t> offsets, neighbors, origins;
    std::vector<ParticleVector3> gradients;
    std::vector<ShapeHessian3> hessians;
    std::vector<WeightedDifferenceEntry> weighted_differences;
    bool add_identity = false;
    bool build_adjoint = true;
    std::vector<SmoothedQuadratureSite> smoothed_sites;
    std::vector<double> biot_point_volumes;
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

struct FiniteVolumeBiotPointResult;

class SpatialDerivativeBatch
{
  public:
    virtual ~SpatialDerivativeBatch() = default;
    virtual void configure(const SpatialDerivativeMap &map) = 0;
    virtual void evaluate(const std::vector<ParticleVector3> &particles,
                          std::vector<SpatialDerivatives3> &result,
                          WeightedDifferenceResult *weighted_difference = nullptr) = 0;
    // Each scalar channel uses the same first derivative stencil. Inputs are
    // indexed by particle then channel, and outputs by row then channel.
    // Origins and differences are retained; no identity term is added.
    virtual void evaluateChannelGradients(const std::vector<double> &particles, std::size_t channels,
                                          std::vector<ParticleVector3> &result) = 0;
    // Transpose of the linear part of evaluate. The constant identity has
    // zero derivative. A difference stencil contributes with both signs.
    virtual void assembleAdjoint(const std::vector<SpatialDerivatives3> &coefficients,
                                 std::vector<ParticleVector3> &result) = 0;
    // Use the deformation derivatives from the preceding evaluate call.
    // Coefficients are also returned for checkpoints and virtual-work audits.
    virtual void assembleSmoothedAdjoint(const std::vector<SmoothedStressField> &stress,
                                         std::vector<SpatialDerivatives3> &coefficients,
                                         std::vector<ParticleVector3> &result) = 0;

    virtual bool supportsResidentBiot() const noexcept { return false; }
    // Evaluates Biot from the latest deformation gradients. A false result
    // means the caller must complete the rejected material points on the host.
    virtual bool evaluateResidentBiot(double, double, std::vector<FiniteVolumeBiotPointResult> &)
    {
        throw std::logic_error("This derivative backend does not evaluate resident Biot material");
    }
    // Reuses accepted Biot stresses with the configured point volumes.
    virtual void assembleResidentBiotAdjoint(std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This derivative backend does not assemble resident Biot forces");
    }
};

std::unique_ptr<SpatialDerivativeBatch> makeCudaSpatialDerivatives();
}

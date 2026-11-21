#pragma once

#include "pd_spatial_derivatives.h"

#include <cstdint>

namespace SPH::peridynamics
{
struct NosbPDBondGeometry
{
    std::size_t particle;
    ParticleVector3 direction;
    double reference_distance, kernel_weight, kernel_derivative;
};

struct NosbPDKinematicSite
{
    ParticleVector3 position, velocity;
    double volume;
    bool active;
};

struct NosbPDForceSite
{
    ShapeHessian3 force_state;
    double volume;
    bool active;
};

struct NosbPDDeformationTensors
{
    ShapeHessian3 full, increment, midpoint;
};

struct NosbPDCorrectionSite
{
    ParticleVector3 position, velocity;
    ShapeHessian3 deformation;
    double volume, damage;
    bool active;
};

struct NosbPDStabilizationParameters
{
    double stiffness, damage_cutoff, collapsed_length_ratio;
};

struct NosbPDDampingParameters
{
    double density, sound_speed, smoothing_length, factor, inverse_kernel_at_origin;
    int dimensions;
};

struct NosbPDCorrectionForce
{
    ParticleVector3 density;
    bool collapsed;
};

class NosbPDBondBatch
{
  public:
    virtual ~NosbPDBondBatch() = default;
    // Each offset row is a real particle, followed by optional referenced
    // ghosts in the separate site field. Bonds retain the CPU neighbor order.
    virtual void configure(std::size_t sites, const std::vector<std::size_t> &offsets,
                           const std::vector<NosbPDBondGeometry> &bonds) = 0;
    virtual void evaluateDeformation(const std::vector<NosbPDKinematicSite> &sites,
                          const std::vector<std::uint8_t> &bond_live, double dt,
                                     std::vector<NosbPDDeformationTensors> &tensors) = 0;
    virtual void assembleInternalForce(const std::vector<NosbPDForceSite> &sites,
                                       const std::vector<std::uint8_t> &bond_live,
                                       std::vector<ParticleVector3> &force_density) = 0;
    virtual void assembleStabilization(const std::vector<NosbPDCorrectionSite> &sites,
                                       const std::vector<std::uint8_t> &bond_live,
                                       NosbPDStabilizationParameters parameters,
                                       std::vector<NosbPDCorrectionForce> &force) = 0;
    virtual void assembleDamping(const std::vector<NosbPDCorrectionSite> &sites,
                                 const std::vector<std::uint8_t> &bond_live,
                                 NosbPDDampingParameters parameters,
                                 std::vector<NosbPDCorrectionForce> &force) = 0;
};

std::unique_ptr<NosbPDBondBatch> makeCudaNosbPDBonds();
}

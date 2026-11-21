#pragma once

#include "pd_spatial_derivatives.h"

#include <cstdint>

namespace SPH::peridynamics
{
struct NosbPDBondGeometry
{
    std::size_t particle;
    ParticleVector3 reference;
    double kernel_weight;
};

struct NosbPDKinematicSite
{
    ParticleVector3 position, velocity;
    double volume;
    bool active;
};

struct NosbPDDeformationTensors
{
    ShapeHessian3 full, increment, midpoint;
};

class NosbPDDeformationBatch
{
  public:
    virtual ~NosbPDDeformationBatch() = default;
    // Each offset row is a real particle, followed by optional referenced
    // ghosts in the separate site field. Bonds retain the CPU neighbor order.
    virtual void configure(std::size_t sites, const std::vector<std::size_t> &offsets,
                           const std::vector<NosbPDBondGeometry> &bonds) = 0;
    virtual void evaluate(const std::vector<NosbPDKinematicSite> &sites,
                          const std::vector<std::uint8_t> &bond_live, double dt,
                          std::vector<NosbPDDeformationTensors> &tensors) = 0;
};

std::unique_ptr<NosbPDDeformationBatch> makeCudaNosbPDDeformation();
}

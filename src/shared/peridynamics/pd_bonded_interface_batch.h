#pragma once

#include "pd_particle_trace_map.h"
#include <memory>
#include <vector>

namespace SPH::peridynamics
{
using BondedTraceEntry = ParticleTraceEntry;
using BondedTraceMap = ParticleTraceMap;

struct BondedInterfaceNode
{
    double position[3], velocity[3], external[3];
    double force_divisor;
};

struct BondedInterfacePoint
{
    double position_a[3], position_b[3], force_a[3], sum[3], torque[3];
    double energy, gap, power, force_scale, power_scale;
};

class BondedInterfaceBatch
{
  public:
    virtual ~BondedInterfaceBatch() = default;
    virtual void configure(const BondedTraceMap &a, const BondedTraceMap &b,
                           const std::vector<double> &stiffness) = 0;
    virtual void apply(std::vector<BondedInterfaceNode> &a,
                       std::vector<BondedInterfaceNode> &b,
                       std::vector<BondedInterfacePoint> &points) = 0;
};

std::unique_ptr<BondedInterfaceBatch> makeCudaBondedInterface();
}

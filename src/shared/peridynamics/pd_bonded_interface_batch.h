#pragma once

#include "pd_coulomb_friction.h"
#include "pd_particle_trace_map.h"
#include <memory>
#include <vector>

namespace SPH::peridynamics
{
using BondedTraceEntry = ParticleTraceEntry;
using BondedTraceMap = ParticleTraceMap;

struct TraceForceNode
{
    double position[3], velocity[3], external[3];
    double force_divisor;
};

using BondedInterfaceNode = TraceForceNode;

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

struct SurfacePlaneContactPoint
{
    double position[3], force[3], sum[3], torque[3];
    double upper, lower, energy, penetration, power, force_scale, power_scale;
};

struct SurfaceFrictionState
{
    double anchor[2], stored;
    int side;
};

struct SurfaceFrictionTrial
{
    double anchor[2], elastic_slip[2];
    CoulombSlipUpdate<double> update;
};

struct SurfaceFrictionTransfer
{
    std::vector<SurfaceFrictionState> committed;
    std::vector<SurfaceFrictionTrial> trials;
};

class SurfacePlaneContactBatch
{
  public:
    virtual ~SurfacePlaneContactBatch() = default;
    virtual void configure(const ParticleTraceMap &trace, const std::vector<double> &stiffness,
                           double friction = 0, double tangent_factor = 1) = 0;
    virtual void apply(double bottom, double top, std::vector<TraceForceNode> &nodes,
                       std::vector<SurfacePlaneContactPoint> &points,
                       SurfaceFrictionTransfer *friction = nullptr) = 0;
};

std::unique_ptr<BondedInterfaceBatch> makeCudaBondedInterface();
std::unique_ptr<SurfacePlaneContactBatch> makeCudaSurfacePlaneContact();
}

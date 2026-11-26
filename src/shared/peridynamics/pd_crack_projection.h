#pragma once

#include "pd_particle_trace_map.h"
#include "pd_spatial_derivatives.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace SPH::peridynamics
{
struct CrackProjectionInput
{
    // Four master corners, followed by the slave point.
    double coordinates[5][3];
};

// Current geometry is shared by every candidate referring to that face or point.
struct CrackProjectionMaster
{
    double corners[4][3];
};

struct CrackProjectionPoint
{
    double position[3];
};

struct CrackProjectionPair
{
    std::size_t master, point;
};

struct CrackProjectionResult
{
    enum Status { Accepted, Outside, InvalidGeometry };
    double s, t;
    int status;
};

// Quantities at a selected slave point and its four master traces.
struct CrackLocalContactInput
{
    double tangent_u[3], tangent_v[3], mixed_tangent[3], normal[3];
    double shape[4], derivative_u[4], derivative_v[4];
    double mass_weight[5];
    double gap, stiffness;
};

struct CrackLocalContactResult
{
    double row_bound[5][3], force_on_slave[3];
    bool valid;
};

struct CrackTraceAddress
{
    std::size_t rows[5];
};

struct CrackTraceParticleState
{
    double mass, external[3];
};

struct CrackTraceParticleResult
{
    double external[3], rate[3];
};

class CrackProjectionBatch
{
  public:
    virtual ~CrackProjectionBatch() = default;
    virtual void project(const std::vector<CrackProjectionInput> &input,
                         std::vector<CrackProjectionResult> &output) = 0;
    virtual void projectIndexed(const std::vector<CrackProjectionMaster> &masters,
                                const std::vector<CrackProjectionPoint> &points,
                                const std::vector<CrackProjectionPair> &pairs,
                                std::vector<CrackProjectionResult> &output)
    {
        // Keep projection-only implementations compatible with indexed callers.
        std::vector<CrackProjectionInput> input(pairs.size());
        for (std::size_t i = 0; i < pairs.size(); ++i)
        {
            const auto &master = masters.at(pairs[i].master);
            const auto &point = points.at(pairs[i].point);
            for (std::size_t corner = 0; corner < 4; ++corner)
                std::copy_n(master.corners[corner], 3, input[i].coordinates[corner]);
            std::copy_n(point.position, 3, input[i].coordinates[4]);
        }
        project(input, output);
    }
    virtual bool assemblesLocalContact() const noexcept { return false; }
    virtual void assembleLocalContact(const std::vector<CrackLocalContactInput> &,
                                      std::vector<CrackLocalContactResult> &)
    {
        throw std::logic_error("This crack backend only implements projection");
    }
    virtual bool assemblesTracePullback() const noexcept { return false; }
    virtual std::uint64_t contactTraceGeneration() const noexcept { return 0; }
    virtual std::uint64_t configureContactTrace(const ParticleTraceMap &)
    {
        throw std::logic_error("This crack backend does not implement trace assembly");
    }
    // Assembles the most recent local response batch in the configured trace
    // generation. The host can retain its reference path for invalid tangents.
    virtual void assembleContactTrace(std::uint64_t,
                                      const std::vector<CrackTraceAddress> &,
                                      const std::vector<CrackTraceParticleState> &,
                                      bool,
                                      std::vector<CrackTraceParticleResult> &)
    {
        throw std::logic_error("This crack backend does not implement trace assembly");
    }

    virtual DeviceParticleVectorView assembleContactTraceOnDevice(
        std::uint64_t, const std::vector<CrackTraceAddress> &, DeviceParticleVectorView,
        const double *, bool, double &)
    {
        throw std::logic_error("This crack backend does not assemble device forces");
    }
    virtual void copyDeviceForces(DeviceParticleVectorView, std::vector<ParticleVector3> &)
    {
        throw std::logic_error("This crack backend does not own device forces");
    }
};

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection(bool assemble_local_contact = false);
}

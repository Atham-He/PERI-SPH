#include "pd_bonded_interface_batch.h"
#include "pd_cuda_buffer.cuh"

#include <cmath>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

struct TraceView
{
    const std::size_t *offsets;
    const BondedTraceEntry *entries;
    const TraceForceNode *nodes;
};

__device__ void gather(TraceView trace, std::size_t point, double *position, double *velocity)
{
    for (auto k = trace.offsets[point]; k < trace.offsets[point + 1]; ++k)
    {
        const auto entry = trace.entries[k];
        const auto &node = trace.nodes[entry.node];
        for (int d = 0; d < 3; ++d)
        {
            position[d] += entry.weight * node.position[d];
            velocity[d] += entry.weight * node.velocity[d];
        }
    }
}

template <class Point>
__device__ void audit(TraceView trace, std::size_t point, double sign, const double *force,
                      double force_norm, Point &result)
{
    for (auto k = trace.offsets[point]; k < trace.offsets[point + 1]; ++k)
    {
        const auto entry = trace.entries[k];
        const auto &node = trace.nodes[entry.node];
        const double w = sign * entry.weight;
        const double fx = w * force[0];
        const double fy = w * force[1];
        const double fz = w * force[2];
        const double *x = node.position;
        const double *v = node.velocity;
        result.sum[0] += fx;
        result.sum[1] += fy;
        result.sum[2] += fz;
        result.torque[0] += x[1] * fz - x[2] * fy;
        result.torque[1] += x[2] * fx - x[0] * fz;
        result.torque[2] += x[0] * fy - x[1] * fx;
        const double power = (fx * v[0] + fy * v[1]) + fz * v[2];
        result.power += power;
        result.force_scale += fabs(w) * force_norm;
        result.power_scale += fabs(power);
    }
}

__global__ void evaluateInterface(TraceView a, TraceView b, const double *stiffness, std::size_t count,
                                  BondedInterfacePoint *points)
{
    const std::size_t point = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (point >= count)
        return;
    BondedInterfacePoint result{};
    double v_a[3]{}, v_b[3]{};
    gather(b, point, result.position_b, v_b);
    const bool active = a.offsets[point] != a.offsets[point + 1] &&
                        b.offsets[point] != b.offsets[point + 1];
    if (active)
    {
        gather(a, point, result.position_a, v_a);
        double gap[3];
        const double k = stiffness[point];
        for (int d = 0; d < 3; ++d)
        {
            gap[d] = result.position_a[d] - result.position_b[d];
            result.force_a[d] = -k * gap[d];
        }
        const double gap_squared = (gap[0] * gap[0] + gap[1] * gap[1]) + gap[2] * gap[2];
        result.energy = .5 * k * gap_squared;
        result.gap = sqrt(gap_squared);
        // Each nodal contribution is a scalar multiple of this point force.
        const double force_norm = k * result.gap;
        audit(a, point, 1, result.force_a, fabs(force_norm), result);
        audit(b, point, -1, result.force_a, fabs(force_norm), result);
        const double energy_rate = -((result.force_a[0] * (v_a[0] - v_b[0]) +
                                      result.force_a[1] * (v_a[1] - v_b[1])) +
                                      result.force_a[2] * (v_a[2] - v_b[2]));
        result.power += energy_rate;
        result.power_scale += fabs(energy_rate);
    }
    else
    {
        for (int d = 0; d < 3; ++d)
            result.position_a[d] = result.position_b[d];
    }
    points[point] = result;
}

__global__ void evaluatePlaneContact(TraceView trace, const double *stiffness, std::size_t count,
                                     double bottom, double top, SurfacePlaneContactPoint *points)
{
    const std::size_t point = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (point >= count)
        return;

    SurfacePlaneContactPoint result{};
    if (trace.offsets[point] != trace.offsets[point + 1])
    {
        double velocity[3]{};
        gather(trace, point, result.position, velocity);
        const double lower = bottom - result.position[2];
        const double upper = result.position[2] - top;
        const double pb = lower < 0 ? 0 : lower;
        const double pt = upper < 0 ? 0 : upper;
        if (pb != 0 || pt != 0)
        {
            const double k = stiffness[point];
            const double normal_force = k * (pb - pt);
            result.force[2] = normal_force;
            result.lower = k * pb;
            result.upper = k * pt;
            result.energy = .5 * k * (pb * pb + pt * pt);
            result.penetration = fmax(pb, pt);
            result.sum[2] = -normal_force;
            result.torque[0] = -result.position[1] * normal_force;
            result.torque[1] = result.position[0] * normal_force;
            result.power = -normal_force * velocity[2];
            result.power_scale = fabs(result.power);
            audit(trace, point, 1, result.force, fabs(normal_force), result);
        }
    }
    points[point] = result;
}

__device__ const double *forceComponents(const BondedInterfacePoint &point) { return point.force_a; }
__device__ const double *forceComponents(const SurfacePlaneContactPoint &point) { return point.force; }

constexpr unsigned int force_lanes = 8;
constexpr unsigned int warp_lanes = 32;

template <bool divide_by_mass, class Point>
__device__ void sumParticleTraceContributions(const ParticleTraceContribution *entries, const Point *points,
                                     std::size_t begin, std::size_t end, double divisor, double *force)
{
    for (auto k = begin; k < end; k += force_lanes)
    {
        const auto entry = entries[k];
        const auto *force_components = forceComponents(points[entry.row]);
        for (int d = 0; d < 3; ++d)
        {
            const double contribution = entry.weight * force_components[d];
            if constexpr (divide_by_mass)
                force[d] += contribution / divisor;
            else
                force[d] += contribution;
        }
    }
}

template <class Point>
__global__ void accumulateNodalForce(const std::size_t *__restrict__ offsets,
                                     const ParticleTraceContribution *__restrict__ entries,
                                     const Point *__restrict__ points, std::size_t count,
                                     TraceForceNode *__restrict__ nodes)
{
    const std::size_t thread = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const std::size_t node = thread / force_lanes;
    unsigned int lane = threadIdx.x % force_lanes;
    if (node >= count)
        return;

    const auto begin = offsets[node], end = offsets[node + 1];
    if (begin == end)
        return;

    const double divisor = nodes[node].force_divisor;
    double force[3]{};
    if (divisor == 1.0)
        sumParticleTraceContributions<false>(entries, points, begin + lane, end, divisor, force);
    else
        sumParticleTraceContributions<true>(entries, points, begin + lane, end, divisor, force);

    // One group owns each output node. The fixed tree keeps the force sum
    // reproducible between GPU runs, including after trace reconstruction.
    const unsigned int warp_lane = threadIdx.x % warp_lanes;
    const unsigned int mask = ((1u << force_lanes) - 1u) << (warp_lane - lane);
    for (unsigned int offset = force_lanes / 2; offset > 0; offset /= 2)
        for (int d = 0; d < 3; ++d)
            force[d] += __shfl_down_sync(mask, force[d], offset, force_lanes);

    if (lane == 0)
        for (int d = 0; d < 3; ++d)
            nodes[node].external[d] += force[d];
}

struct DeviceTrace
{
    DeviceBuffer<std::size_t> offsets, transpose_offsets;
    DeviceBuffer<BondedTraceEntry> entries;
    DeviceBuffer<ParticleTraceContribution> contributions;
    DeviceBuffer<TraceForceNode> nodes;
    std::size_t count = 0;

    template <class IncludeRow>
    void configure(const ParticleTraceMap &map, IncludeRow include_row, double sign)
    {
        count = map.node_count;
        upload(offsets, map.offsets);
        upload(entries, map.entries);

        const auto transpose = transposeParticleTrace(map, include_row, sign);
        upload(transpose_offsets, transpose.offsets);
        upload(contributions, transpose.entries);
        nodes.reserve(count);
    }

    TraceView view() const { return {offsets.get(), entries.get(), nodes.get()}; }
};

constexpr unsigned int threads = 128;

unsigned int launchBlocks(std::size_t count, int maximum)
{
    const auto blocks = (count + threads - 1) / threads;
    if (blocks > static_cast<std::size_t>(maximum))
        throw std::length_error("CUDA surface trace exceeds the device grid");
    return static_cast<unsigned int>(blocks);
}

template <class Point>
void accumulate(DeviceTrace &trace, const DeviceBuffer<Point> &points, int maximum)
{
    if (trace.count == 0)
        return;
    static_assert(threads % warp_lanes == 0);
    accumulateNodalForce<<<launchBlocks(trace.count * force_lanes, maximum), threads>>>(trace.transpose_offsets.get(),
        trace.contributions.get(), points.get(), trace.count, trace.nodes.get());
    check(cudaGetLastError(), "surface trace force accumulation");
}

class CudaBondedInterface final : public BondedInterfaceBatch
{
  public:
    CudaBondedInterface()
    {
        int device = 0;
        check(cudaGetDevice(&device), "bonded interface device selection");
        check(cudaGetDeviceProperties(&properties_, device), "bonded interface device properties");
    }

    void configure(const BondedTraceMap &a, const BondedTraceMap &b,
                   const std::vector<double> &stiffness) override
    {
        validateParticleTrace(a, stiffness.size());
        validateParticleTrace(b, stiffness.size());
        a_.configure(a, [&](std::size_t row) { return b.offsets[row] != b.offsets[row + 1]; }, 1);
        b_.configure(b, [&](std::size_t row) { return a.offsets[row] != a.offsets[row + 1]; }, -1);
        upload(stiffness_, stiffness);
        point_count_ = stiffness.size();
        points_.reserve(point_count_);
    }

    void apply(std::vector<BondedInterfaceNode> &a, std::vector<BondedInterfaceNode> &b,
               std::vector<BondedInterfacePoint> &points) override
    {
        if (a.size() != a_.count || b.size() != b_.count)
            throw std::invalid_argument("CUDA bonded interface node counts changed without map rebuild");
        points.resize(point_count_);
        if (point_count_ == 0)
            return;
        upload(a_.nodes, a);
        upload(b_.nodes, b);
        evaluateInterface<<<launchBlocks(point_count_, properties_.maxGridSize[0]), threads>>>(a_.view(), b_.view(), stiffness_.get(), point_count_, points_.get());
        check(cudaGetLastError(), "bonded interface evaluation");
        accumulate(a_, points_, properties_.maxGridSize[0]);
        accumulate(b_, points_, properties_.maxGridSize[0]);
        download(points, points_);
        download(a, a_.nodes);
        download(b, b_.nodes);
    }

  private:
    DeviceTrace a_, b_;
    DeviceBuffer<double> stiffness_;
    DeviceBuffer<BondedInterfacePoint> points_;
    std::size_t point_count_ = 0;
    cudaDeviceProp properties_{};

};

class CudaSurfacePlaneContact final : public SurfacePlaneContactBatch
{
  public:
    CudaSurfacePlaneContact()
    {
        int device = 0;
        check(cudaGetDevice(&device), "plane contact device selection");
        check(cudaDeviceGetAttribute(&max_grid_, cudaDevAttrMaxGridDimX, device), "plane contact device grid");
    }

    void configure(const ParticleTraceMap &trace, const std::vector<double> &stiffness) override
    {
        validateParticleTrace(trace, stiffness.size());
        trace_.configure(trace, [](std::size_t) { return true; }, 1);
        upload(stiffness_, stiffness);
        point_count_ = stiffness.size();
        points_.reserve(point_count_);
    }

    void apply(double bottom, double top, std::vector<TraceForceNode> &nodes,
               std::vector<SurfacePlaneContactPoint> &points) override
    {
        if (nodes.size() != trace_.count)
            throw std::invalid_argument("CUDA plane contact node count changed without map rebuild");
        points.resize(point_count_);
        if (point_count_ == 0)
            return;
        upload(trace_.nodes, nodes);
        evaluatePlaneContact<<<launchBlocks(point_count_, max_grid_), threads>>>(
            trace_.view(), stiffness_.get(), point_count_, bottom, top, points_.get());
        check(cudaGetLastError(), "plane contact evaluation");
        accumulate(trace_, points_, max_grid_);
        download(points, points_);
        download(nodes, trace_.nodes);
    }

  private:
    DeviceTrace trace_;
    DeviceBuffer<double> stiffness_;
    DeviceBuffer<SurfacePlaneContactPoint> points_;
    std::size_t point_count_ = 0;
    int max_grid_ = 0;
};
}

std::unique_ptr<BondedInterfaceBatch> makeCudaBondedInterface()
{
    return std::make_unique<CudaBondedInterface>();
}

std::unique_ptr<SurfacePlaneContactBatch> makeCudaSurfacePlaneContact()
{
    return std::make_unique<CudaSurfacePlaneContact>();
}
}

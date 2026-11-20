#include "pd_bonded_interface_batch.h"
#include "pd_cuda_buffer.cuh"

#include <cmath>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;

struct NodeContribution
{
    std::size_t point;
    double weight;
};

struct TraceView
{
    const std::size_t *offsets;
    const BondedTraceEntry *entries;
    const BondedInterfaceNode *nodes;
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

__device__ void audit(TraceView trace, std::size_t point, double sign, BondedInterfacePoint &result)
{
    for (auto k = trace.offsets[point]; k < trace.offsets[point + 1]; ++k)
    {
        const auto entry = trace.entries[k];
        const auto &node = trace.nodes[entry.node];
        const double w = sign * entry.weight;
        const double fx = w * result.force_a[0];
        const double fy = w * result.force_a[1];
        const double fz = w * result.force_a[2];
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
        result.force_scale += sqrt((fx * fx + fy * fy) + fz * fz);
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
        audit(a, point, 1, result);
        audit(b, point, -1, result);
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

__global__ void accumulateNodalForce(const std::size_t *offsets, const NodeContribution *entries,
                                     const BondedInterfacePoint *points, std::size_t count,
                                     BondedInterfaceNode *nodes)
{
    const std::size_t node = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (node >= count || offsets[node] == offsets[node + 1])
        return;
    auto &state = nodes[node];
    // Each node sums its entries in the original trace order.
    for (auto k = offsets[node]; k < offsets[node + 1]; ++k)
    {
        const auto entry = entries[k];
        const auto &point = points[entry.point];
        for (int d = 0; d < 3; ++d)
            state.external[d] += (entry.weight * point.force_a[d]) / state.force_divisor;
    }
}

template <class T>
void upload(DeviceBuffer<T> &device, const std::vector<T> &host)
{
    device.reserve(host.size());
    if (!host.empty())
        check(cudaMemcpy(device.get(), host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice),
              "bonded interface upload");
}

template <class T>
void download(std::vector<T> &host, const DeviceBuffer<T> &device)
{
    if (!host.empty())
        check(cudaMemcpy(host.data(), device.get(), host.size() * sizeof(T), cudaMemcpyDeviceToHost),
              "bonded interface download");
}

struct DeviceTrace
{
    DeviceBuffer<std::size_t> offsets, transpose_offsets;
    DeviceBuffer<BondedTraceEntry> entries;
    DeviceBuffer<NodeContribution> contributions;
    DeviceBuffer<BondedInterfaceNode> nodes;
    std::size_t count = 0;

    void configure(const BondedTraceMap &map, const BondedTraceMap &other, double sign)
    {
        count = map.node_count;
        upload(offsets, map.offsets);
        upload(entries, map.entries);

        std::vector<std::size_t> transpose(count + 1, 0);
        for (std::size_t point = 0; point + 1 < map.offsets.size(); ++point)
            if (other.offsets[point] != other.offsets[point + 1])
                for (auto k = map.offsets[point]; k < map.offsets[point + 1]; ++k)
                    ++transpose[map.entries[k].node + 1];
        for (std::size_t i = 1; i < transpose.size(); ++i)
            transpose[i] += transpose[i - 1];

        std::vector<NodeContribution> sorted(transpose.back());
        auto cursor = transpose;
        for (std::size_t point = 0; point + 1 < map.offsets.size(); ++point)
            if (other.offsets[point] != other.offsets[point + 1])
                for (auto k = map.offsets[point]; k < map.offsets[point + 1]; ++k)
                {
                    const auto entry = map.entries[k];
                    sorted[cursor[entry.node]++] = {point, sign * entry.weight};
                }

        upload(transpose_offsets, transpose);
        upload(contributions, sorted);
        nodes.reserve(count);
    }

    TraceView view() const { return {offsets.get(), entries.get(), nodes.get()}; }
};

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
        validate(a, stiffness.size());
        validate(b, stiffness.size());
        a_.configure(a, b, 1);
        b_.configure(b, a, -1);
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
        evaluateInterface<<<grid(point_count_), threads>>>(a_.view(), b_.view(), stiffness_.get(), point_count_, points_.get());
        check(cudaGetLastError(), "bonded interface evaluation");
        accumulate(a_);
        accumulate(b_);
        download(points, points_);
        download(a, a_.nodes);
        download(b, b_.nodes);
    }

  private:
    static constexpr unsigned int threads = 128;
    DeviceTrace a_, b_;
    DeviceBuffer<double> stiffness_;
    DeviceBuffer<BondedInterfacePoint> points_;
    std::size_t point_count_ = 0;
    cudaDeviceProp properties_{};

    static void validate(const BondedTraceMap &map, std::size_t points)
    {
        if (map.offsets.size() != points + 1 || map.offsets.front() != 0 ||
            map.offsets.back() != map.entries.size())
            throw std::invalid_argument("Invalid CUDA bonded trace offsets");
        for (std::size_t i = 1; i < map.offsets.size(); ++i)
            if (map.offsets[i] < map.offsets[i - 1])
                throw std::invalid_argument("Unordered CUDA bonded trace offsets");
        for (const auto &entry : map.entries)
            if (entry.node >= map.node_count)
                throw std::invalid_argument("CUDA bonded trace contains a node outside its port");
    }

    unsigned int grid(std::size_t count) const
    {
        const auto blocks = (count + threads - 1) / threads;
        if (blocks > static_cast<std::size_t>(properties_.maxGridSize[0]))
            throw std::length_error("CUDA bonded interface exceeds the device grid");
        return static_cast<unsigned int>(blocks);
    }

    void accumulate(DeviceTrace &trace)
    {
        if (trace.count == 0)
            return;
        accumulateNodalForce<<<grid(trace.count), threads>>>(trace.transpose_offsets.get(),
            trace.contributions.get(), points_.get(), trace.count, trace.nodes.get());
        check(cudaGetLastError(), "bonded interface force accumulation");
    }
};
}

std::unique_ptr<BondedInterfaceBatch> makeCudaBondedInterface()
{
    return std::make_unique<CudaBondedInterface>();
}
}

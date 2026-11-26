#include "pd_bilinear_contact_projection.h"
#include "pd_cuda_buffer.cuh"

#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

using crack_projection_detail::Vector3;

__global__ void projectBatch(const CrackProjectionInput *packed,
                             const CrackProjectionMaster *masters, const CrackProjectionPoint *points,
                             const CrackProjectionPair *pairs, CrackProjectionResult *out, std::size_t count)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto *corners = pairs ? masters[pairs[i].master].corners : packed[i].coordinates;
    const auto &point = pairs ? points[pairs[i].point].position : packed[i].coordinates[4];
    out[i] = projectCrackContact(corners, point);
}

__device__ Vector3 loadVector(const double (&x)[3])
{
    return {x[0], x[1], x[2]};
}

__device__ CrackLocalContactResult localContact(const CrackLocalContactInput &input)
{
    CrackLocalContactResult response{};
    const Vector3 u = loadVector(input.tangent_u);
    const Vector3 v = loadVector(input.tangent_v);
    const Vector3 normal = loadVector(input.normal);
    const double penetration = fmin(0.0, input.gap);
    const Vector3 residual = normal * penetration;
    const double a = u.dot(u);
    const double b = u.dot(v) - residual.dot(loadVector(input.mixed_tangent));
    const double c = v.dot(v);
    const double largest = .5 * (a + c) + hypot(.5 * (a - c), b);
    const double determinant = a * c - b * b;

    // det / largest is the smaller eigenvalue of the symmetric tangent.
    // The host eigensolver handles borderline geometry and reports failures.
    if (!isfinite(determinant) || !isfinite(largest) || largest <= 0.0 ||
        determinant <= 1.e-13 * largest * largest)
        return response;

    const double kinv00 = c / determinant;
    const double kinv01 = -b / determinant;
    const double kinv11 = a / determinant;
    double coefficient[5] = {1.0, -input.shape[0], -input.shape[1], -input.shape[2], -input.shape[3]};
    double coupling_u[15], coupling_v[15];
    for (int i = 0; i < 5; ++i)
    {
        const Vector3 pu = i == 0 ? u * (-coefficient[i]) :
            u * (-coefficient[i]) - residual * input.derivative_u[i - 1];
        const Vector3 pv = i == 0 ? v * (-coefficient[i]) :
            v * (-coefficient[i]) - residual * input.derivative_v[i - 1];
        coupling_u[3 * i] = pu.x;
        coupling_u[3 * i + 1] = pu.y;
        coupling_u[3 * i + 2] = pu.z;
        coupling_v[3 * i] = pv.x;
        coupling_v[3 * i + 1] = pv.y;
        coupling_v[3 * i + 2] = pv.z;
    }

    for (int row = 0; row < 15; ++row)
    {
        const double inverse_u = -coupling_u[row] * kinv00 - coupling_v[row] * kinv01;
        const double inverse_v = -coupling_u[row] * kinv01 - coupling_v[row] * kinv11;
        double bound = 0.0;
        for (int column = 0; column < 15; ++column)
        {
            // Condense out the motion of the two closest-point coordinates.
            double h = inverse_u * coupling_u[column] + inverse_v * coupling_v[column];
            if (row % 3 == column % 3)
                h += coefficient[row / 3] * coefficient[column / 3];
            bound += fabs(input.stiffness * h) * input.mass_weight[column / 3];
        }
        if (!isfinite(bound))
            return response;
        response.row_bound[row / 3][row % 3] = bound;
    }
    if (input.gap < 0.0)
    {
        const Vector3 force = normal * (-input.stiffness * input.gap);
        response.force_on_slave[0] = force.x;
        response.force_on_slave[1] = force.y;
        response.force_on_slave[2] = force.z;
    }
    response.valid = true;
    return response;
}

__global__ void localContactBatch(const CrackLocalContactInput *input, CrackLocalContactResult *output,
                                 std::size_t count)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        output[i] = localContact(input[i]);
}

struct ContactTraceUse
{
    std::size_t contact;
    unsigned int leg;
};

struct TraceLoad
{
    double force[3], rate[3];
};

__global__ void collectTraceLoads(const CrackLocalContactInput *input,
                                  const CrackLocalContactResult *response,
                                  const std::size_t *offsets, const ContactTraceUse *uses,
                                  std::size_t trace_count, TraceLoad *output)
{
    const auto trace = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (trace >= trace_count)
        return;

    TraceLoad load{};
    for (auto i = offsets[trace]; i < offsets[trace + 1]; ++i)
    {
        const auto use = uses[i];
        const auto &contact = response[use.contact];
        const double coefficient = use.leg == 0 ? 1.0 : -input[use.contact].shape[use.leg - 1];
        for (int d = 0; d < 3; ++d)
        {
            load.force[d] += coefficient * contact.force_on_slave[d];
            load.rate[d] += contact.row_bound[use.leg][d];
        }
    }
    output[trace] = load;
}

__device__ CrackTraceParticleResult crackTraceForce(
    const std::size_t *node_offsets, const ParticleTraceContribution *entries, const TraceLoad *traces,
    std::size_t i, const double *external, double mass, bool acceleration_port)
{
    const double divisor = acceleration_port ? mass : 1.0;
    const double sqrt_mass = sqrt(mass);
    CrackTraceParticleResult state{{external[0], external[1], external[2]}, {}};
    for (auto k = node_offsets[i]; k < node_offsets[i + 1]; ++k)
    {
        const auto entry = entries[k];
        const auto &trace = traces[entry.row];
        for (int d = 0; d < 3; ++d)
        {
            state.external[d] += (entry.weight * trace.force[d]) / divisor;
            state.rate[d] += (fabs(entry.weight) * trace.rate[d]) / sqrt_mass;
        }
    }
    return state;
}

__global__ void pullBackCrackTraces(const std::size_t *node_offsets,
                                    const ParticleTraceContribution *entries,
                                    const TraceLoad *traces,
                                    const CrackTraceParticleState *particles, bool acceleration_port,
                                    std::size_t count, CrackTraceParticleResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        output[i] = crackTraceForce(node_offsets, entries, traces, i, particles[i].external,
                                   particles[i].mass, acceleration_port);
}

__global__ void pullBackDeviceCrackTraces(const std::size_t *node_offsets,
                                          const ParticleTraceContribution *entries,
                                          const TraceLoad *traces, const double *mass,
                                          DeviceParticleVectorView base, bool acceleration_port,
                                          ParticleVector3 *output, double *block_rates)
{
    extern __shared__ double local_rates[];
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    double rate = 0;
    if (i < base.count)
    {
        const auto state = crackTraceForce(node_offsets, entries, traces, i,
                                           base.data[i].component, mass[i], acceleration_port);
        output[i] = {{state.external[0], state.external[1], state.external[2]}};
        rate = fmax(state.rate[0], fmax(state.rate[1], state.rate[2]));
    }
    local_rates[threadIdx.x] = rate;
    __syncthreads();
    for (unsigned stride = blockDim.x / 2; stride != 0; stride /= 2)
    {
        if (threadIdx.x < stride)
            local_rates[threadIdx.x] = fmax(local_rates[threadIdx.x], local_rates[threadIdx.x + stride]);
        __syncthreads();
    }
    if (threadIdx.x == 0)
        block_rates[blockIdx.x] = local_rates[0];
}

class CudaCrackProjection final : public CrackProjectionBatch
{
  public:
    explicit CudaCrackProjection(bool assemble_local_contact) : assemble_local_contact_(assemble_local_contact)
    {
        int device = 0;
        check(cudaGetDevice(&device), "device selection");
        check(cudaGetDeviceProperties(&properties_, device), "device properties");
    }

    void project(const std::vector<CrackProjectionInput> &input,
                 std::vector<CrackProjectionResult> &output) override
    {
        const std::size_t count = input.size();
        output.resize(count);
        if (count == 0)
            return;

        const auto blocks = grid(count);
        input_.reserve(count);
        output_.reserve(count);
        check(cudaMemcpy(input_.get(), input.data(), count * sizeof(input[0]),
                         cudaMemcpyHostToDevice), "input transfer");
        projectBatch<<<blocks, threads>>>(input_.get(), nullptr, nullptr, nullptr, output_.get(), count);
        check(cudaGetLastError(), "kernel launch");
        check(cudaMemcpy(output.data(), output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "output transfer");
    }

    void projectIndexed(const std::vector<CrackProjectionMaster> &masters,
                        const std::vector<CrackProjectionPoint> &points,
                        const std::vector<CrackProjectionPair> &pairs,
                        std::vector<CrackProjectionResult> &output) override
    {
        output.resize(pairs.size());
        if (pairs.empty())
            return;

        const auto blocks = grid(pairs.size());
        stageProjectionPairs(masters.size(), points.size(), pairs);
        upload(masters_, masters, "crack projection master faces");
        upload(points_, points, "crack projection integration points");
        output_.reserve(pairs.size());
        projectBatch<<<blocks, threads>>>(nullptr, masters_.get(), points_.get(), pairs_.get(),
                                         output_.get(), pairs.size());
        check(cudaGetLastError(), "project indexed crack faces");
        download(output, output_, "indexed crack projection results");
    }

    bool assemblesLocalContact() const noexcept override { return assemble_local_contact_; }

    void assembleLocalContact(const std::vector<CrackLocalContactInput> &input,
                              std::vector<CrackLocalContactResult> &output) override
    {
        if (!assemble_local_contact_)
            throw std::logic_error("Local contact assembly is disabled for this projection backend");
        const auto count = input.size();
        local_contact_count_ = count;
        local_contact_valid_ = count == 0;
        output.resize(count);
        if (count == 0)
            return;
        const auto blocks = grid(count);
        contact_input_.reserve(count);
        contact_output_.reserve(count);
        check(cudaMemcpy(contact_input_.get(), input.data(), count * sizeof(input[0]),
                         cudaMemcpyHostToDevice), "crack tangent input transfer");
        localContactBatch<<<blocks, threads>>>(contact_input_.get(), contact_output_.get(), count);
        check(cudaGetLastError(), "crack tangent kernel launch");
        check(cudaMemcpy(output.data(), contact_output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "crack tangent output transfer");
        local_contact_valid_ = std::all_of(output.begin(), output.end(),
                                           [](const auto &value) { return value.valid; });
    }

    bool assemblesTracePullback() const noexcept override { return assemble_local_contact_; }
    std::uint64_t contactTraceGeneration() const noexcept override { return trace_ready_ ? trace_generation_ : 0; }

    std::uint64_t configureContactTrace(const ParticleTraceMap &map) override
    {
        if (!assemble_local_contact_)
            throw std::logic_error("Trace assembly is disabled for this projection backend");
        if (map.offsets.empty())
            throw std::invalid_argument("Crack trace offsets are missing");
        const auto rows = map.offsets.size() - 1;
        validateParticleTrace(map, rows);
        const auto transpose = transposeParticleTrace(map, [](std::size_t) { return true; });
        trace_ready_ = false;
        cached_contact_address_.clear();
        upload(node_trace_offsets_, transpose.offsets, "crack node trace offsets");
        upload(node_trace_entries_, transpose.entries, "crack node trace entries");
        trace_count_ = rows;
        particle_count_ = map.node_count;
        if (++trace_generation_ == 0)
            ++trace_generation_;
        trace_ready_ = true;
        return trace_generation_;
    }

    void assembleContactTrace(std::uint64_t generation, const std::vector<CrackTraceAddress> &address,
                              const std::vector<CrackTraceParticleState> &particles, bool acceleration_port,
                              std::vector<CrackTraceParticleResult> &output) override
    {
        if (generation == 0 || generation != contactTraceGeneration() || !local_contact_valid_ ||
            particles.size() != particle_count_ ||
            address.size() != local_contact_count_)
            throw std::invalid_argument("Crack contact batch does not match the active trace generation");

        for (const auto &particle : particles)
            if (!std::isfinite(particle.mass) || particle.mass <= 0.0)
                throw std::invalid_argument("Crack trace assembly requires positive finite particle masses");
        output.resize(particles.size());
        if (address.empty())
        {
            for (std::size_t i = 0; i < particles.size(); ++i)
            {
                const double *external = particles[i].external;
                output[i] = {{external[0], external[1], external[2]}, {}};
            }
            return;
        }
        if (particles.empty())
            return;
        collectDeviceTraceLoads(address);
        upload(trace_particle_input_, particles, "crack trace particle state");
        trace_particle_output_.reserve(particle_count_);
        pullBackCrackTraces<<<grid(particle_count_), threads>>>(node_trace_offsets_.get(), node_trace_entries_.get(),
            trace_loads_.get(), trace_particle_input_.get(), acceleration_port, particle_count_, trace_particle_output_.get());
        check(cudaGetLastError(), "pull back crack trace loads");
        download(output, trace_particle_output_, "crack trace particle forces and rates");
    }

    DeviceParticleVectorView assembleContactTraceOnDevice(
        std::uint64_t generation, const std::vector<CrackTraceAddress> &address,
        DeviceParticleVectorView base, const double *mass, bool acceleration_port, double &max_rate) override
    {
        if (generation == 0 || generation != contactTraceGeneration() || !local_contact_valid_ ||
            !base || !mass || base.count != particle_count_ || address.size() != local_contact_count_)
            throw std::invalid_argument("Device crack forces do not match the active trace generation");
        max_rate = 0;
        if (address.empty())
            return base;

        collectDeviceTraceLoads(address);
        trace_device_forces_.reserve(base.count);
        const auto blocks = grid(base.count);
        trace_block_rates_.reserve(blocks);
        pullBackDeviceCrackTraces<<<blocks, threads, threads * sizeof(double)>>>(
            node_trace_offsets_.get(), node_trace_entries_.get(), trace_loads_.get(), mass,
            base, acceleration_port, trace_device_forces_.get(), trace_block_rates_.get());
        check(cudaGetLastError(), "pull back crack forces into device state");
        host_block_rates_.resize(blocks);
        download(host_block_rates_, trace_block_rates_, "crack contact stability reduction");
        max_rate = *std::max_element(host_block_rates_.begin(), host_block_rates_.end());
        return {trace_device_forces_.get(), base.count};
    }

    void copyDeviceForces(DeviceParticleVectorView force, std::vector<ParticleVector3> &output) override
    {
        output.resize(force.count);
        if (force.count)
            check(cudaMemcpy(output.data(), force.data, force.count * sizeof(ParticleVector3), cudaMemcpyDeviceToHost),
                  "read device forces for host crack fallback");
    }

  private:
    void collectDeviceTraceLoads(const std::vector<CrackTraceAddress> &address)
    {
        stageContactAddress(address);
        trace_loads_.reserve(trace_count_);
        if (trace_count_ > 0)
        {
            collectTraceLoads<<<grid(trace_count_), threads>>>(contact_input_.get(), contact_output_.get(),
                trace_use_device_offsets_.get(), trace_use_device_entries_.get(), trace_count_, trace_loads_.get());
            check(cudaGetLastError(), "collect crack trace loads");
        }
    }

    void stageProjectionPairs(std::size_t master_count, std::size_t point_count,
                              const std::vector<CrackProjectionPair> &pairs)
    {
        const auto same_pair = [](const CrackProjectionPair &a, const CrackProjectionPair &b)
        {
            return a.master == b.master && a.point == b.point;
        };
        if (master_count == cached_master_count_ && point_count == cached_point_count_ &&
            std::equal(pairs.begin(), pairs.end(), cached_pairs_.begin(), cached_pairs_.end(), same_pair))
            return;

        for (const auto &pair : pairs)
            if (pair.master >= master_count || pair.point >= point_count)
                throw std::out_of_range("Crack projection pair refers to missing geometry");
        cached_pairs_.clear();
        upload(pairs_, pairs, "crack projection candidate pairs");
        cached_pairs_ = pairs;
        cached_master_count_ = master_count;
        cached_point_count_ = point_count;
    }

    void stageContactAddress(const std::vector<CrackTraceAddress> &address)
    {
        const auto same_address = [](const CrackTraceAddress &a, const CrackTraceAddress &b)
        {
            return std::equal(a.rows, a.rows + 5, b.rows);
        };
        if (std::equal(address.begin(), address.end(), cached_contact_address_.begin(),
                       cached_contact_address_.end(), same_address))
            return;

        trace_use_offsets_.assign(trace_count_ + 1, 0);
        for (const auto &contact : address)
            for (auto trace : contact.rows)
            {
                if (trace >= trace_count_)
                    throw std::invalid_argument("Contact references a trace outside the active map");
                ++trace_use_offsets_[trace + 1];
            }
        for (std::size_t i = 1; i < trace_use_offsets_.size(); ++i)
            trace_use_offsets_[i] += trace_use_offsets_[i - 1];
        trace_uses_.resize(trace_use_offsets_.back());
        trace_use_cursor_ = trace_use_offsets_;
        for (std::size_t i = 0; i < address.size(); ++i)
            for (unsigned int leg = 0; leg < 5; ++leg)
                trace_uses_[trace_use_cursor_[address[i].rows[leg]]++] = {i, leg};
        cached_contact_address_.clear();
        upload(trace_use_device_offsets_, trace_use_offsets_, "crack trace use offsets");
        upload(trace_use_device_entries_, trace_uses_, "crack trace uses");
        cached_contact_address_ = address;
    }

    unsigned int grid(std::size_t count) const
    {
        const auto blocks = (count + threads - 1) / threads;
        if (blocks > static_cast<std::size_t>(properties_.maxGridSize[0]))
            throw std::length_error("CUDA crack batch exceeds the device grid");
        return static_cast<unsigned int>(blocks);
    }

    static constexpr unsigned int threads = 128;
    bool assemble_local_contact_;
    cudaDeviceProp properties_{};
    DeviceBuffer<CrackProjectionInput> input_;
    DeviceBuffer<CrackProjectionMaster> masters_;
    DeviceBuffer<CrackProjectionPoint> points_;
    DeviceBuffer<CrackProjectionPair> pairs_;
    std::vector<CrackProjectionPair> cached_pairs_;
    std::size_t cached_master_count_ = 0, cached_point_count_ = 0;
    DeviceBuffer<CrackProjectionResult> output_;
    DeviceBuffer<CrackLocalContactInput> contact_input_;
    DeviceBuffer<CrackLocalContactResult> contact_output_;
    std::size_t local_contact_count_ = 0, trace_count_ = 0, particle_count_ = 0;
    std::uint64_t trace_generation_ = 0;
    bool trace_ready_ = false, local_contact_valid_ = false;
    std::vector<std::size_t> trace_use_offsets_, trace_use_cursor_;
    std::vector<ContactTraceUse> trace_uses_;
    std::vector<CrackTraceAddress> cached_contact_address_;
    DeviceBuffer<std::size_t> node_trace_offsets_, trace_use_device_offsets_;
    DeviceBuffer<ParticleTraceContribution> node_trace_entries_;
    DeviceBuffer<ContactTraceUse> trace_use_device_entries_;
    DeviceBuffer<TraceLoad> trace_loads_;
    DeviceBuffer<CrackTraceParticleState> trace_particle_input_;
    DeviceBuffer<CrackTraceParticleResult> trace_particle_output_;
    DeviceBuffer<ParticleVector3> trace_device_forces_;
    DeviceBuffer<double> trace_block_rates_;
    std::vector<double> host_block_rates_;
};
}

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection(bool assemble_local_contact)
{
    return std::make_unique<CudaCrackProjection>(assemble_local_contact);
}
}

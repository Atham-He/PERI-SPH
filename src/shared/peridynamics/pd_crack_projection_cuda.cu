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

__global__ void projectBatch(const CrackProjectionInput *in, CrackProjectionResult *out,std::size_t n) {
    const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = projectCrackContact(in[i]);
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

__global__ void pullBackCrackTraces(const std::size_t *node_offsets,
                                    const ParticleTraceContribution *entries,
                                    const TraceLoad *traces,
                                    const CrackTraceParticleState *particles, bool acceleration_port,
                                    std::size_t count, CrackTraceParticleResult *output)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &particle = particles[i];
    const double divisor = acceleration_port ? particle.mass : 1.0;
    const double sqrt_mass = sqrt(particle.mass);
    CrackTraceParticleResult state{{particle.external[0], particle.external[1], particle.external[2]}, {}};
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
    output[i] = state;
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
        projectBatch<<<blocks, threads>>>(input_.get(), output_.get(), count);
        check(cudaGetLastError(), "kernel launch");
        check(cudaMemcpy(output.data(), output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "output transfer");
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

        for (const auto &particle : particles)
            if (!std::isfinite(particle.mass) || particle.mass <= 0.0)
                throw std::invalid_argument("Crack trace assembly requires positive finite particle masses");
        output.resize(particles.size());
        if (particles.empty())
            return;

        trace_loads_.reserve(trace_count_);
        if (trace_count_ > 0)
        {
            upload(trace_use_device_offsets_, trace_use_offsets_, "crack trace use offsets");
            upload(trace_use_device_entries_, trace_uses_, "crack trace uses");
            collectTraceLoads<<<grid(trace_count_), threads>>>(contact_input_.get(), contact_output_.get(),
                trace_use_device_offsets_.get(), trace_use_device_entries_.get(), trace_count_, trace_loads_.get());
            check(cudaGetLastError(), "collect crack trace loads");
        }

        upload(trace_particle_input_, particles, "crack trace particle state");
        trace_particle_output_.reserve(particle_count_);
        pullBackCrackTraces<<<grid(particle_count_), threads>>>(node_trace_offsets_.get(), node_trace_entries_.get(),
            trace_loads_.get(), trace_particle_input_.get(), acceleration_port, particle_count_, trace_particle_output_.get());
        check(cudaGetLastError(), "pull back crack trace loads");
        download(output, trace_particle_output_, "crack trace particle forces and rates");
    }

  private:
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
    DeviceBuffer<CrackProjectionResult> output_;
    DeviceBuffer<CrackLocalContactInput> contact_input_;
    DeviceBuffer<CrackLocalContactResult> contact_output_;
    std::size_t local_contact_count_ = 0, trace_count_ = 0, particle_count_ = 0;
    std::uint64_t trace_generation_ = 0;
    bool trace_ready_ = false, local_contact_valid_ = false;
    std::vector<std::size_t> trace_use_offsets_, trace_use_cursor_;
    std::vector<ContactTraceUse> trace_uses_;
    DeviceBuffer<std::size_t> node_trace_offsets_, trace_use_device_offsets_;
    DeviceBuffer<ParticleTraceContribution> node_trace_entries_;
    DeviceBuffer<ContactTraceUse> trace_use_device_entries_;
    DeviceBuffer<TraceLoad> trace_loads_;
    DeviceBuffer<CrackTraceParticleState> trace_particle_input_;
    DeviceBuffer<CrackTraceParticleResult> trace_particle_output_;
};
}

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection(bool assemble_local_contact)
{
    return std::make_unique<CudaCrackProjection>(assemble_local_contact);
}
}

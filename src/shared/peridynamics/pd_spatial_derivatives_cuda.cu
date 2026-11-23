#include "pd_spatial_derivatives.h"
#include "pd_finite_volume_biot_cuda.cuh"
#include "pd_smoothed_force_projection.h"
#include "pd_cuda_buffer.cuh"

#include <algorithm>
#include <limits>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

template <bool second_order>
__global__ void evaluateSpatialDerivatives(const std::size_t *offsets, const std::size_t *neighbors,
                                           const std::size_t *origins, const ParticleVector3 *gradients,
                                           const ShapeHessian3 *hessians, const ParticleVector3 *particles,
                                           bool add_identity, std::size_t rows, SpatialDerivatives3 *result)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto row = index / 3;
    const auto axis = index % 3;
    if (row >= rows)
        return;

    double gradient[3]{}, hessian[3][3]{};
    if (add_identity)
        gradient[axis] = 1.0;
    const double origin = origins ? particles[origins[row]].component[axis] : 0.0;
    for (auto i = offsets[row]; i < offsets[row + 1]; ++i)
    {
        const double value = particles[neighbors[i]].component[axis] - origin;
        for (int j = 0; j < 3; ++j)
        {
            gradient[j] += value * gradients[i].component[j];
            if constexpr (second_order)
                for (int k = 0; k < 3; ++k)
                    hessian[j][k] += value * hessians[i].derivative[j][k];
        }
    }
    for (int j = 0; j < 3; ++j)
    {
        result[row].gradient[axis][j] = gradient[j];
        for (int k = 0; k < 3; ++k)
            result[row].hessian[axis][j][k] = hessian[j][k];
    }
}

__global__ void evaluateChannelGradientsKernel(const std::size_t *offsets, const std::size_t *neighbors,
                                               const std::size_t *origins, const ParticleVector3 *weights,
                                               const double *particles, std::size_t channels, std::size_t entries,
                                               ParticleVector3 *result)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= entries)
        return;
    const auto row = index / channels;
    const auto channel = index % channels;
    const double origin = origins ? particles[origins[row] * channels + channel] : 0.0;
    ParticleVector3 gradient{};
    for (auto i = offsets[row]; i < offsets[row + 1]; ++i)
    {
        const double value = particles[neighbors[i] * channels + channel] - origin;
        for (int axis = 0; axis < 3; ++axis)
            gradient.component[axis] += value * weights[i].component[axis];
    }
    result[index] = gradient;
}

struct AdjointContribution
{
    std::size_t row, entry;
    int sign;
};

__global__ void evaluateWeightedDifferences(const std::size_t *offsets, const std::size_t *neighbors,
                                            const std::size_t *origins, const WeightedDifferenceEntry *weights,
                                            const ParticleVector3 *particles, std::size_t rows,
                                            ParticleVector3 *edge_forces, ParticleVector3 *row_energies)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto row = index / 3;
    const auto axis = index % 3;
    if (row >= rows)
        return;

    const double origin = particles[origins[row]].component[axis];
    double fit[9]{};
    for (auto entry = offsets[row]; entry < offsets[row + 1]; ++entry)
    {
        const double difference = particles[neighbors[entry]].component[axis] - origin;
        for (int term = 0; term < 9; ++term)
            fit[term] += difference * weights[entry].coefficient[term];
    }

    double energy = 0.0;
    for (auto entry = offsets[row]; entry < offsets[row + 1]; ++entry)
    {
        double prediction = 0.0;
        for (int term = 0; term < 9; ++term)
            prediction += fit[term] * weights[entry].basis[term];

        const double residual = (particles[neighbors[entry]].component[axis] - origin) - prediction;
        energy += weights[entry].energy_weight * residual * residual;
        edge_forces[entry].component[axis] = weights[entry].force_weight * residual;
    }
    row_energies[row].component[axis] = energy;
}

__global__ void gatherWeightedDifferenceForces(const std::size_t *offsets, const AdjointContribution *entries,
                                               const ParticleVector3 *edge_forces, std::size_t particle_count,
                                               ParticleVector3 *forces)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto particle = index / 3;
    const auto axis = index % 3;
    if (particle >= particle_count)
        return;

    double force = 0.0;
    for (auto i = offsets[particle]; i < offsets[particle + 1]; ++i)
    {
        const auto entry = entries[i];
        force -= entry.sign * edge_forces[entry.entry].component[axis];
    }
    forces[particle].component[axis] = force;
}

template <bool second_order>
__global__ void assembleSpatialAdjoint(const std::size_t *offsets, const AdjointContribution *entries,
                                       const ParticleVector3 *gradients, const ShapeHessian3 *hessians,
                                       const SpatialDerivatives3 *coefficients, std::size_t particle_count,
                                       ParticleVector3 *result)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto particle = index / 3;
    const auto axis = index % 3;
    if (particle >= particle_count)
        return;

    double force = 0.0;
    for (auto i = offsets[particle]; i < offsets[particle + 1]; ++i)
    {
        const auto entry = entries[i];
        const auto &coefficient = coefficients[entry.row];
        const auto &gradient = gradients[entry.entry];
        const double first = (coefficient.gradient[axis][0] * gradient.component[0] +
                              coefficient.gradient[axis][1] * gradient.component[1]) +
                              coefficient.gradient[axis][2] * gradient.component[2];
        force += entry.sign * first;
        if constexpr (second_order)
        {
            double second = 0.0;
            const auto &hessian = hessians[entry.entry];
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    second += coefficient.hessian[axis][j][k] * hessian.derivative[j][k];
            force += entry.sign * second;
        }
    }
    result[particle].component[axis] = force;
}

__device__ SmoothedForceTensor loadTensor(const double (&value)[3][3])
{
    SmoothedForceTensor result;
    for (int i=0;i<3;++i)
        for (int j=0;j<3;++j)
            result(i,j)=value[i][j];
    return result;
}

__global__ void projectSmoothedCoefficients(const SmoothedQuadratureSite *quadrature,
                                            const SmoothedStressField *stress_field,
                                            const SpatialDerivatives3 *kinematics,
                                            std::size_t rows, SpatialDerivatives3 *result)
{
    const auto l = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (l >= rows)
        return;

    const auto &state=kinematics[l];
    const auto &stress=stress_field[l];
    const auto &site=quadrature[l];
    const SmoothedForceTensor F=loadTensor(state.gradient), inverse_F=F.inverse();
    SmoothedForceGradient trial_gradient, correction_gradient, displacement_hessian, nodal_hessian;
    for (int a=0;a<3;++a)
    {
        trial_gradient[a]=loadTensor(stress.trial_gradient[a].derivative);
        correction_gradient[a]=loadTensor(stress.correction_gradient[a].derivative);
        displacement_hessian[a]=loadTensor(state.hessian[a]);
    }
    SmoothedForceTensor nodal_gradient;
    projectSmoothedForce(F,inverse_F,site.volume,loadTensor(site.moment.derivative),
        loadTensor(stress.stress.derivative),trial_gradient,correction_gradient,displacement_hessian,
        nodal_gradient,nodal_hessian);
    for (int a=0;a<3;++a)
        for (int i=0;i<3;++i)
        {
            result[l].gradient[a][i]=nodal_gradient(a,i);
            for (int j=0;j<3;++j)
                result[l].hessian[a][i][j]=nodal_hessian[a](i,j);
        }
}

__global__ void evaluateResidentBiotKernel(const SpatialDerivatives3 *deformation,std::size_t count,
                                              double bulk,double shear,FiniteVolumeBiotPointResult *result)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i<count)
        result[i]=cuda_detail::evaluateFiniteVolumeBiotPoint(deformation[i].gradient,bulk,shear);
}

__global__ void weightBiotAdjoint(const FiniteVolumeBiotPointResult *material,const double *volume,
                                  std::size_t count,SpatialDerivatives3 *coefficients)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i>=count) return;
    for (int a=0;a<3;++a)
        for (int b=0;b<3;++b)
            coefficients[i].gradient[a][b]=volume[i]*material[i].piola.derivative[a][b];
}

class CudaSpatialDerivatives final : public SpatialDerivativeBatch
{
  public:
    CudaSpatialDerivatives()
    {
        int device = 0;
        check(cudaGetDevice(&device), "select spatial derivative device");
        cudaDeviceProp properties{};
        check(cudaGetDeviceProperties(&properties, device), "get spatial derivative device properties");
        max_blocks_ = static_cast<std::size_t>(properties.maxGridSize[0]);
    }

    void configure(const SpatialDerivativeMap &map) override
    {
        const auto entries = map.neighbors.size();
        if (map.offsets.empty() || map.offsets.front() != 0 || map.offsets.back() != entries ||
            !std::is_sorted(map.offsets.begin(), map.offsets.end()) || map.gradients.size() != entries ||
            (!map.hessians.empty() && map.hessians.size() != entries))
            throw std::invalid_argument("Invalid spatial derivative stencil");
        const auto rows = map.offsets.size() - 1;
        if (!map.smoothed_sites.empty() &&
            (map.smoothed_sites.size()!=rows || map.hessians.empty() || !map.build_adjoint || !map.add_identity))
            throw std::invalid_argument("Smoothed force projection requires matching integration sites, deformation Hessians and an adjoint stencil");
        if (!map.biot_point_volumes.empty() && (map.biot_point_volumes.size()!=rows || !map.build_adjoint))
            throw std::invalid_argument("Biot force assembly requires row volumes and an adjoint stencil");
        if (!map.weighted_differences.empty() &&
            (map.weighted_differences.size() != entries || map.origins.empty() || !map.build_adjoint))
            throw std::invalid_argument("Weighted differences require entry fits, row origins and the adjoint stencil");
        if (!map.origins.empty() && map.origins.size() != rows)
            throw std::invalid_argument("Invalid spatial derivative origins");
        for (auto node : map.neighbors)
            if (node >= map.particle_count)
                throw std::invalid_argument("Spatial derivative neighbor is outside the particle field");
        for (auto node : map.origins)
            if (node >= map.particle_count)
                throw std::invalid_argument("Spatial derivative origin is outside the particle field");

        configured_ = false;
        derivative_result_ready_ = false;
        resident_biot_ready_ = false;
        upload(offsets_, map.offsets, "spatial derivative offsets");
        upload(neighbors_, map.neighbors, "spatial derivative neighbors");
        upload(origins_, map.origins, "spatial derivative origins");
        upload(gradients_, map.gradients, "reference gradients");
        upload(hessians_, map.hessians, "reference Hessians");
        upload(smoothed_sites_, map.smoothed_sites, "smoothed integration sites");
        upload(biot_point_volumes_,map.biot_point_volumes,"Biot integration volumes");
        has_biot_volumes_ = !map.biot_point_volumes.empty();
        has_smoothed_sites_ = !map.smoothed_sites.empty();
        upload(weighted_differences_, map.weighted_differences, "reference weighted differences");
        if (map.build_adjoint)
            configureAdjoint(map);
        has_adjoint_ = map.build_adjoint;
        particle_count_ = map.particle_count;
        rows_ = rows;
        use_origins_ = !map.origins.empty();
        use_hessians_ = !map.hessians.empty();
        add_identity_ = map.add_identity;
        weighted_differences_count_ = map.weighted_differences.size();
        configured_ = true;
    }

    void evaluate(const std::vector<ParticleVector3> &particles,
                  std::vector<SpatialDerivatives3> &result,
                  WeightedDifferenceResult *weighted_difference) override
    {
        if (!configured_ || particles.size() != particle_count_)
            throw std::invalid_argument("Particle field does not match the spatial derivative stencil");
        if (weighted_difference && weighted_differences_count_ == 0)
            throw std::invalid_argument("Stencil has no weighted difference fit");
        derivative_result_ready_ = false;
        resident_biot_ready_ = false;
        result.resize(rows_);
        if (rows_ == 0)
            return;

        const auto blocks = grid(rows_);
        upload(particles_, particles, "spatial derivative particle field");
        result_.reserve(rows_);
        const auto *origins = use_origins_ ? origins_.get() : nullptr;
        if (use_hessians_)
            evaluateSpatialDerivatives<true><<<blocks, threads>>>(
                offsets_.get(), neighbors_.get(), origins, gradients_.get(), hessians_.get(),
                particles_.get(), add_identity_, rows_, result_.get());
        else
            evaluateSpatialDerivatives<false><<<blocks, threads>>>(
                offsets_.get(), neighbors_.get(), origins, gradients_.get(), nullptr,
                particles_.get(), add_identity_, rows_, result_.get());
        check(cudaGetLastError(), "spatial derivative evaluation");
        if (weighted_difference)
            evaluateWeightedDifference(*weighted_difference);
        download(result, result_, "spatial derivative results");
        derivative_result_ready_ = true;
    }

    void evaluateChannelGradients(const std::vector<double> &particles, std::size_t channels,
                                  std::vector<ParticleVector3> &result) override
    {
        constexpr auto maximum = std::numeric_limits<std::size_t>::max();
        if (!configured_ || channels == 0 || particle_count_ > maximum / channels ||
            particles.size() != particle_count_ * channels)
            throw std::invalid_argument("Channel field does not match the spatial derivative stencil");
        if (rows_ > maximum / channels)
            throw std::length_error("Spatial derivative channel result exceeds its index range");
        const auto entries = rows_ * channels;
        const auto blocks = grid(entries, 1);
        result.resize(entries);
        if (entries == 0)
            return;
        upload(channel_field_, particles, "spatial derivative channels");
        channel_gradients_.reserve(entries);
        const auto *origins = use_origins_ ? origins_.get() : nullptr;
        evaluateChannelGradientsKernel<<<blocks, threads>>>(offsets_.get(), neighbors_.get(), origins,
            gradients_.get(), channel_field_.get(), channels, entries, channel_gradients_.get());
        check(cudaGetLastError(), "spatial derivative channel evaluation");
        download(result, channel_gradients_, "spatial derivative channel results");
    }

    void assembleAdjoint(const std::vector<SpatialDerivatives3> &coefficients,
                         std::vector<ParticleVector3> &result) override
    {
        if (!configured_ || !has_adjoint_ || coefficients.size() != rows_)
            throw std::invalid_argument("Coefficient field does not match an adjoint stencil");
        upload(adjoint_coefficients_, coefficients, "spatial adjoint coefficients");
        gatherAdjoint(result,use_hessians_);
    }

    void assembleSmoothedAdjoint(const std::vector<SmoothedStressField> &stress,
                                 std::vector<SpatialDerivatives3> &coefficients,
                                 std::vector<ParticleVector3> &result) override
    {
        if (!configured_ || !has_smoothed_sites_ || !derivative_result_ready_ || stress.size()!=rows_)
            throw std::invalid_argument("Smoothed force projection requires a current deformation and matching stress field");

        upload(smoothed_stress_,stress,"smoothed stress fields");
        adjoint_coefficients_.reserve(rows_);
        projectSmoothedCoefficients<<<grid(rows_,1),threads>>>(smoothed_sites_.get(),smoothed_stress_.get(),
            result_.get(),rows_,adjoint_coefficients_.get());
        check(cudaGetLastError(),"project smoothed force coefficients");

        gatherAdjoint(result,true);
        coefficients.resize(rows_);
        download(coefficients,adjoint_coefficients_,"smoothed force coefficients");
    }

    bool supportsResidentBiot() const noexcept override { return true; }

    bool evaluateResidentBiot(double bulk,double shear,std::vector<FiniteVolumeBiotPointResult> &result) override
    {
        if (!derivative_result_ready_)
            throw std::logic_error("Resident Biot material requires a current deformation gradient");
        resident_biot_ready_=false;
        result.resize(rows_);
        biot_material_.reserve(rows_);
        evaluateResidentBiotKernel<<<grid(rows_,1),threads>>>(result_.get(),rows_,bulk,shear,biot_material_.get());
        check(cudaGetLastError(),"evaluate Biot material from resident gradients");
        download(result,biot_material_,"resident Biot material results");
        resident_biot_ready_=std::all_of(result.begin(),result.end(),
                                         [](const auto &point){return point.accepted;});
        return resident_biot_ready_;
    }

    void assembleResidentBiotAdjoint(std::vector<ParticleVector3> &result) override
    {
        if (!has_biot_volumes_ || !resident_biot_ready_)
            throw std::logic_error("Resident Biot forces require accepted stresses and integration volumes");
        adjoint_coefficients_.reserve(rows_);
        weightBiotAdjoint<<<grid(rows_,1),threads>>>(biot_material_.get(),biot_point_volumes_.get(),rows_,
                                                    adjoint_coefficients_.get());
        check(cudaGetLastError(),"weight resident Biot Piola stresses");
        gatherAdjoint(result,false);
    }

  private:
    void gatherAdjoint(std::vector<ParticleVector3> &result,bool second_order)
    {
        result.resize(particle_count_);
        if (particle_count_ == 0)
            return;
        const auto blocks=grid(particle_count_);
        adjoint_result_.reserve(particle_count_);
        if (second_order)
            assembleSpatialAdjoint<true><<<blocks, threads>>>(adjoint_offsets_.get(), adjoint_entries_.get(),
                gradients_.get(), hessians_.get(), adjoint_coefficients_.get(), particle_count_, adjoint_result_.get());
        else
            assembleSpatialAdjoint<false><<<blocks, threads>>>(adjoint_offsets_.get(), adjoint_entries_.get(),
                gradients_.get(), nullptr, adjoint_coefficients_.get(), particle_count_, adjoint_result_.get());
        check(cudaGetLastError(), "assemble spatial adjoint");
        download(result, adjoint_result_, "spatial adjoint particle result");
    }

    static constexpr unsigned int threads = 128;

    unsigned int grid(std::size_t sites, std::size_t components = 3) const
    {
        if (sites == 0)
            return 0;
        if (sites > std::numeric_limits<std::size_t>::max() / components)
            throw std::length_error("Spatial derivative batch exceeds its index range");
        const std::size_t blocks = (sites * components - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("Spatial derivative batch exceeds the device grid");
        return static_cast<unsigned int>(blocks);
    }

    void evaluateWeightedDifference(WeightedDifferenceResult &result)
    {
        edge_forces_.reserve(weighted_differences_count_);
        row_energies_.reserve(rows_);
        weighted_forces_.reserve(particle_count_);
        evaluateWeightedDifferences<<<grid(rows_), threads>>>(
            offsets_.get(), neighbors_.get(), origins_.get(), weighted_differences_.get(), particles_.get(),
            rows_, edge_forces_.get(), row_energies_.get());
        check(cudaGetLastError(), "evaluate weighted difference fit");
        gatherWeightedDifferenceForces<<<grid(particle_count_), threads>>>(
            adjoint_offsets_.get(), adjoint_entries_.get(), edge_forces_.get(), particle_count_, weighted_forces_.get());
        check(cudaGetLastError(), "gather weighted difference forces");

        host_energies_.resize(rows_);
        download(host_energies_, row_energies_, "weighted difference energies");
        result.energy_density.resize(rows_);
        for (std::size_t row = 0; row < rows_; ++row)
        {
            const auto &energy = host_energies_[row].component;
            result.energy_density[row] = energy[0] + energy[1] + energy[2];
        }
        result.forces.resize(particle_count_);
        download(result.forces, weighted_forces_, "weighted difference forces");
    }

    void configureAdjoint(const SpatialDerivativeMap &map)
    {
        // The traversal is stable in reference stencil order, including
        // negative contributions at difference origins.
        auto for_each_contribution = [&](auto apply)
        {
            for (std::size_t row = 0; row + 1 < map.offsets.size(); ++row)
                for (auto entry = map.offsets[row]; entry < map.offsets[row + 1]; ++entry)
                {
                    apply(map.neighbors[entry], AdjointContribution{row, entry, 1});
                    if (!map.origins.empty())
                        apply(map.origins[row], AdjointContribution{row, entry, -1});
                }
        };
        std::vector<std::size_t> offsets(map.particle_count + 1, 0);
        for_each_contribution([&](std::size_t node, const AdjointContribution &) { ++offsets[node + 1]; });
        for (std::size_t i = 1; i < offsets.size(); ++i)
            offsets[i] += offsets[i - 1];

        std::vector<AdjointContribution> entries(offsets.back());
        auto cursor = offsets;
        for_each_contribution([&](std::size_t node, const AdjointContribution &entry)
        {
            entries[cursor[node]++] = entry;
        });
        upload(adjoint_offsets_, offsets, "spatial adjoint offsets");
        upload(adjoint_entries_, entries, "spatial adjoint entries");
    }

    std::size_t particle_count_ = 0, rows_ = 0, max_blocks_ = 0;
    bool use_origins_ = false, use_hessians_ = false, add_identity_ = false, configured_ = false;
    bool has_adjoint_ = false, has_smoothed_sites_ = false, derivative_result_ready_ = false;
    bool has_biot_volumes_ = false, resident_biot_ready_ = false;
    DeviceBuffer<double> biot_point_volumes_;
    DeviceBuffer<FiniteVolumeBiotPointResult> biot_material_;
    DeviceBuffer<std::size_t> offsets_, neighbors_, origins_;
    DeviceBuffer<ParticleVector3> gradients_, particles_, channel_gradients_;
    DeviceBuffer<double> channel_field_;
    DeviceBuffer<ShapeHessian3> hessians_;
    DeviceBuffer<SmoothedQuadratureSite> smoothed_sites_;
    DeviceBuffer<SmoothedStressField> smoothed_stress_;
    DeviceBuffer<SpatialDerivatives3> result_, adjoint_coefficients_;
    DeviceBuffer<std::size_t> adjoint_offsets_;
    DeviceBuffer<AdjointContribution> adjoint_entries_;
    DeviceBuffer<ParticleVector3> adjoint_result_;
    std::size_t weighted_differences_count_ = 0;
    DeviceBuffer<WeightedDifferenceEntry> weighted_differences_;
    DeviceBuffer<ParticleVector3> edge_forces_, row_energies_, weighted_forces_;
    std::vector<ParticleVector3> host_energies_;
};
}

std::unique_ptr<SpatialDerivativeBatch> makeCudaSpatialDerivatives()
{
    return std::make_unique<CudaSpatialDerivatives>();
}
}

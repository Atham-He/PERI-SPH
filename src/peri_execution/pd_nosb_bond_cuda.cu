#include "pd_nosb_bond_batch.h"
#include "pd_cuda_buffer.cuh"

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

__global__ void evaluateBondTensors(const std::size_t *offsets, const NosbPDBondGeometry *bonds,
                                    const NosbPDKinematicSite *sites, const std::uint8_t *live,
                                    std::size_t rows, double dt, NosbPDDeformationTensors *tensors)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto i = index / 3;
    const auto axis = index % 3;
    if (i >= rows)
        return;

    double full[3]{}, increment[3]{}, midpoint[3]{};
    const auto end = sites[i].active ? offsets[i + 1] : offsets[i];
    for (auto entry = offsets[i]; entry < end; ++entry)
    {
        const auto &bond = bonds[entry];
        const auto &target = sites[bond.particle];
        double eta = (-bond.reference_distance * bond.direction.component[axis]);
        double du = 0.0;
        double eta_half = eta;
        if (live[entry])
        {
            eta = target.position.component[axis] - sites[i].position.component[axis];
            du = (target.velocity.component[axis] - sites[i].velocity.component[axis]) * dt;
            eta_half = eta - du * 0.5;
        }

        const double weight = target.volume * bond.kernel_weight;
        for (int k = 0; k < 3; ++k)
        {
            const double reference = -bond.reference_distance * bond.direction.component[k];
            full[k] += weight * (eta * reference);
            increment[k] += weight * (du * reference);
            midpoint[k] += weight * (eta_half * reference);
        }
    }

    for (int k = 0; k < 3; ++k)
    {
        tensors[i].full.derivative[axis][k] = full[k];
        tensors[i].increment.derivative[axis][k] = increment[k];
        tensors[i].midpoint.derivative[axis][k] = midpoint[k];
    }
}

template <bool Reference>
__global__ void assemblePositionMoments(const std::size_t *offsets, const NosbPDBondGeometry *bonds,
                                        const NosbPDPositionSite *sites, const double *reference_volume,
                                        std::size_t rows, ShapeHessian3 *moments)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto i = index / 3;
    const auto axis = index % 3;
    if (i >= rows)
        return;

    double row[3]{};
    for (auto entry = offsets[i]; entry < offsets[i + 1]; ++entry)
    {
        const auto &bond = bonds[entry];
        const auto j = bond.particle;
        const double weight = bond.kernel_weight * (Reference ? reference_volume[j] : sites[j].volume);
        const double bond_axis = -bond.reference_distance * bond.direction.component[axis];
        const double displacement = Reference ? 0 : sites[j].position.component[axis] - sites[i].position.component[axis];
        for (int k = 0; k < 3; ++k)
        {
            const double reference = -bond.reference_distance * bond.direction.component[k];
            if constexpr (Reference)
                row[k] += weight * (bond_axis * reference);
            else
                row[k] += (weight * displacement) * reference;
        }
    }
    for (int k = 0; k < 3; ++k)
        moments[i].derivative[axis][k] = row[k];
}

__device__ double projectBond(const ShapeHessian3 &state, std::size_t axis, double weight,
                              const ParticleVector3 &reference, double sign)
{
    const auto &row = state.derivative[axis];
    return (weight * row[0]) * (sign * reference.component[0]) +
        ((weight * row[1]) * (sign * reference.component[1]) +
         (weight * row[2]) * (sign * reference.component[2]));
}

__global__ void assembleInternalBondForce(const std::size_t *offsets, const NosbPDBondGeometry *bonds,
                                          const NosbPDForceSite *sites, const std::uint8_t *live,
                                          std::size_t rows, ParticleVector3 *force_density)
{
    const auto index = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto i = index / 3;
    const auto axis = index % 3;
    if (i >= rows)
        return;

    const auto &source = sites[i];
    double force = 0.0;
    const auto end = source.active ? offsets[i + 1] : offsets[i];
    for (auto entry = offsets[i]; entry < end; ++entry)
    {
        if (!live[entry])
            continue;

        const auto &bond = bonds[entry];
        const auto &target = sites[bond.particle];
        ParticleVector3 reference;
        for (int k = 0; k < 3; ++k)
            reference.component[k] = -bond.reference_distance * bond.direction.component[k];
        const double from = projectBond(source.force_state, axis, bond.kernel_weight, reference, 1.0);
        const double to = projectBond(target.force_state, axis, bond.kernel_weight, reference, -1.0);
        force += (from - to) * target.volume;
    }
    force_density[i].component[axis] = force;
}

__device__ ParticleVector3 difference(const ParticleVector3 &lhs, const ParticleVector3 &rhs)
{
    ParticleVector3 result;
    for (int k = 0; k < 3; ++k)
        result.component[k] = lhs.component[k] - rhs.component[k];
    return result;
}

__device__ double dot(const ParticleVector3 &a, const ParticleVector3 &b)
{
    return a.component[0] * b.component[0] +
        (a.component[1] * b.component[1] + a.component[2] * b.component[2]);
}

__device__ ShapeHessian3 averageDeformation(const NosbPDCorrectionSite &source, const NosbPDCorrectionSite &target)
{
    ShapeHessian3 average;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            average.derivative[a][b] = 0.5 * (source.deformation.derivative[a][b] + target.deformation.derivative[a][b]);
    return average;
}

__device__ bool correctionApplies(NosbPDStabilizationParameters parameters, const NosbPDCorrectionSite &source)
{
    return source.active && !(source.damage >= parameters.damage_cutoff);
}

__device__ bool correctionApplies(NosbPDDampingParameters, const NosbPDCorrectionSite &source)
{
    return source.active;
}

__device__ bool accumulateCorrection(NosbPDStabilizationParameters parameters,
                                      const NosbPDCorrectionSite &source, const NosbPDCorrectionSite &target,
                                      const NosbPDBondGeometry &bond, ParticleVector3 &force)
{
    if (!target.active || target.damage >= parameters.damage_cutoff)
        return true;

    const auto current = difference(target.position, source.position);
    const double length = sqrt(dot(current, current));
    if (!isfinite(length) || length <= parameters.collapsed_length_ratio * bond.reference_distance)
        return false;

    const auto deformation = averageDeformation(source, target);
    ParticleVector3 reference, error;
    for (int k = 0; k < 3; ++k)
        reference.component[k] = -bond.reference_distance * bond.direction.component[k];
    for (int k = 0; k < 3; ++k)
        error.component[k] = projectBond(deformation, k, 1.0, reference, 1.0) - current.component[k];
    const double projection = dot(error, current);
    const double scale = parameters.stiffness * projection / (bond.reference_distance * length) * target.volume;
    for (int k = 0; k < 3; ++k)
        force.component[k] -= scale * current.component[k];
    return true;
}

__device__ bool accumulateCorrection(NosbPDDampingParameters parameters,
                                      const NosbPDCorrectionSite &source, const NosbPDCorrectionSite &target,
                                      const NosbPDBondGeometry &bond, ParticleVector3 &force)
{
    const double inverse_distance = double(parameters.dimensions) / bond.reference_distance;
    const auto position_jump = difference(source.position, target.position);
    const auto velocity_jump = difference(source.velocity, target.velocity);
    const double strain_rate = inverse_distance * inverse_distance * dot(position_jump, velocity_jump);
    const double stress_scale = 0.5 * parameters.density * parameters.sound_speed * strain_rate * parameters.smoothing_length;
    const double weight = bond.kernel_weight * parameters.inverse_kernel_at_origin;
    const double gradient = bond.kernel_derivative * target.volume;
    auto stress = averageDeformation(source, target);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            stress.derivative[a][b] = gradient * ((parameters.factor * weight) * (stress.derivative[a][b] * stress_scale));
    for (int k = 0; k < 3; ++k)
        force.component[k] += projectBond(stress, k, 1.0, bond.direction, 1.0);
    return true;
}

template <class Parameters>
__global__ void assembleCorrectionForce(const std::size_t *offsets, const NosbPDBondGeometry *bonds,
                                        const NosbPDCorrectionSite *sites, const std::uint8_t *live,
                                        std::size_t rows, Parameters parameters, NosbPDCorrectionForce *force)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= rows)
        return;

    const auto &source = sites[i];
    NosbPDCorrectionForce row{};
    const auto end = correctionApplies(parameters, source) ? offsets[i + 1] : offsets[i];
    for (auto entry = offsets[i]; entry < end; ++entry)
    {
        if (!live[entry])
            continue;
        const auto &bond = bonds[entry];
        if (!accumulateCorrection(parameters, source, sites[bond.particle], bond, row.density))
        {
            row.collapsed = true;
            break;
        }
    }
    force[i] = row;
}

// Extreme eigenvalues of a finite symmetric tensor, scaled before forming
// the characteristic invariants. The host resolves decisions near a limit.
__device__ void principalExtremes(double (&matrix)[3][3], int dim, double scale,
                                 double &minimum, double &maximum)
{
    if (scale == 0)
    {
        minimum = maximum = 0;
        return;
    }
    for (int a = 0; a < dim; ++a)
        for (int b = 0; b < dim; ++b)
            matrix[a][b] /= scale;

    if (dim == 2)
    {
        const double center = (matrix[0][0] + matrix[1][1]) * .5;
        const double radius = hypot((matrix[0][0] - matrix[1][1]) * .5, matrix[0][1]);
        minimum = scale * (center - radius);
        maximum = scale * (center + radius);
        return;
    }

    const double center = (matrix[0][0] + matrix[1][1] + matrix[2][2]) / 3;
    double squared = 0;
    for (int a = 0; a < 3; ++a)
    {
        matrix[a][a] -= center;
        for (int b = 0; b < 3; ++b)
            squared += matrix[a][b] * matrix[a][b];
    }
    const double radius = sqrt(squared / 6);
    if (radius == 0)
    {
        minimum = maximum = center * scale;
        return;
    }
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            matrix[a][b] /= radius;
    const auto &a = matrix;
    const double determinant = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                               a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                               a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    const double angle = acos(fmax(-1.0, fmin(1.0, 0.5 * determinant))) / 3;
    maximum = scale * (center + 2 * radius * cos(angle));
    minimum = scale * (center + 2 * radius * cos(angle + 2.0943951023931954923));
}

__device__ NosbPDBondDecision screenFailure(const NosbPDBondFailureSite &source, const NosbPDBondFailureSite &target,
                                            double reference_distance, NosbPDBondFailureParameters parameters)
{
    using Decision = NosbPDBondDecision;
    using Law = NosbPDBondFailureLaw;
    constexpr double epsilon = 2.2204460492503131e-16;
    const bool plastic = parameters.law == Law::PrincipalPlasticStrain;
    const auto delta = difference(target.position, source.position);
    const double stretch = (sqrt(dot(delta, delta)) - reference_distance) / reference_distance;
    if (!isfinite(stretch))
        return Decision::EvaluateOnHost;
    const double stretch_limit = plastic ? 1.0 : 0.5;
    if (fabs(stretch - stretch_limit) <= 32 * epsilon * fmax(1.0, fabs(stretch)))
        return Decision::EvaluateOnHost;
    if (plastic && stretch > stretch_limit)
        return Decision::Broken;

    const auto dim = parameters.dimensions;
    double tensor[3][3]{};
    double scale = 0, trace = 0;
    for (int a = 0; a < dim; ++a)
        for (int b = 0; b < dim; ++b)
        {
            const double ab = .5 * (source.tensor.derivative[a][b] + target.tensor.derivative[a][b]);
            const double ba = .5 * (source.tensor.derivative[b][a] + target.tensor.derivative[b][a]);
            tensor[a][b] = .5 * (ab + ba);
            if (!isfinite(tensor[a][b]))
                return Decision::EvaluateOnHost;
            scale = fmax(scale, fabs(tensor[a][b]));
            if (a == b)
                trace += ab;
        }

    double minimum, maximum;
    principalExtremes(tensor, dim, scale, minimum, maximum);
    const bool shear = parameters.law == Law::StressShear;
    const bool pressure = parameters.law == Law::StressPressure;
    const double secondary = shear ? (maximum - minimum) * .5 : -trace / dim;
    // Near repeated eigenvalues the inverse cosine can lose order sqrt(eps).
    // Delegate these threshold decisions to the original host eigensolver.
    const double margin = 64 * sqrt(epsilon) * scale;
    if (!isfinite(maximum) || !isfinite(minimum) ||
        fabs(maximum - parameters.limit) <= margin ||
        ((shear || pressure) && (!isfinite(secondary) || fabs(secondary - parameters.secondary_limit) <= margin)))
        return Decision::EvaluateOnHost;
    const bool fails = maximum > parameters.limit || stretch > stretch_limit ||
                       ((shear || pressure) && secondary > parameters.secondary_limit);
    return fails ? Decision::Broken : Decision::Live;
}

__global__ void screenBondFailures(const std::size_t *offsets, const NosbPDBondGeometry *bonds,
                                   const NosbPDBondFailureSite *sites, const std::uint8_t *live,
                                   NosbPDBondFailureParameters parameters, NosbPDBondDecision *decisions)
{
    const std::size_t i = blockIdx.x;
    for (auto entry = offsets[i] + threadIdx.x; entry < offsets[i + 1]; entry += blockDim.x)
    {
        auto decision = live[entry] ? NosbPDBondDecision::Live : NosbPDBondDecision::Broken;
        if (sites[i].active && live[entry])
        {
            const auto &bond = bonds[entry];
            decision = screenFailure(sites[i], sites[bond.particle], bond.reference_distance, parameters);
        }
        decisions[entry] = decision;
    }
}

class CudaNosbPDBonds final : public NosbPDBondBatch
{
  public:
    CudaNosbPDBonds()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select NOSB-PD deformation device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read NOSB-PD deformation grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void configure(std::size_t sites, const std::vector<std::size_t> &offsets,
                   const std::vector<NosbPDBondGeometry> &bonds) override
    {
        configured_ = false;
        if (offsets.empty() || offsets.front() != 0 || offsets.back() != bonds.size() ||
            offsets.size() - 1 > sites || !std::is_sorted(offsets.begin(), offsets.end()) ||
            std::any_of(bonds.begin(), bonds.end(), [sites](const auto &bond) { return bond.particle >= sites; }))
            throw std::invalid_argument("Invalid NOSB-PD deformation stencil");

        rows_ = offsets.size() - 1;
        if (rows_ > max_blocks_ * threads_ / 3)
            throw std::length_error("NOSB-PD deformation stencil exceeds the device grid");
        upload(offsets_, offsets, "NOSB-PD deformation offsets");
        upload(bonds_, bonds, "NOSB-PD reference bonds");
        site_count_ = sites;
        bond_count_ = bonds.size();
        configured_ = true;
    }

    void evaluateDeformation(const std::vector<NosbPDKinematicSite> &sites, const std::vector<std::uint8_t> &bond_live,
                  double dt, std::vector<NosbPDDeformationTensors> &tensors) override
    {
        validate(sites.size(), bond_live.size());
        tensors.resize(rows_);
        if (rows_ == 0)
            return;

        const auto blocks = static_cast<unsigned>((3 * rows_ - 1) / threads_ + 1);
        upload(sites_, sites, "NOSB-PD kinematic sites");
        upload(live_, bond_live, "NOSB-PD live bond flags");
        tensors_.reserve(rows_);
        evaluateBondTensors<<<blocks, threads_>>>(offsets_.get(), bonds_.get(), sites_.get(), live_.get(),
            rows_, dt, tensors_.get());
        check(cudaGetLastError(), "evaluate NOSB-PD bond tensors");
        download(tensors, tensors_, "NOSB-PD deformation tensors");
    }

    bool assemblesReferenceShape() const noexcept override { return true; }

    void assembleReferenceShape(const std::vector<double> &volume, std::vector<ShapeHessian3> &moments) override
    {
        validate(volume.size(), bond_count_);
        moments.resize(rows_);
        if (rows_ == 0)
            return;
        upload(reference_volume_, volume, "NOSB-PD reference particle volumes");
        accumulatePositionMoments<true>(nullptr, reference_volume_.get(), moments);
    }

    bool assemblesIntactPositionMoments() const noexcept override { return true; }

    void assembleIntactPositionMoments(const std::vector<NosbPDPositionSite> &sites,
                                      std::vector<ShapeHessian3> &moments) override
    {
        validate(sites.size(), bond_count_);
        moments.resize(rows_);
        if (rows_ == 0)
            return;
        upload(position_sites_, sites, "NOSB-PD positions and volumes");
        accumulatePositionMoments<false>(position_sites_.get(), nullptr, moments);
    }

    bool screensBondFailure() const noexcept override { return true; }

    void screenBondFailure(const std::vector<NosbPDBondFailureSite> &sites,
                           const std::vector<std::uint8_t> &bond_live, NosbPDBondFailureParameters parameters,
                           std::vector<NosbPDBondDecision> &decisions) override
    {
        validate(sites.size(), bond_live.size());
        if (parameters.dimensions != 2 && parameters.dimensions != 3)
            throw std::invalid_argument("NOSB-PD failure screening requires two or three dimensions");
        decisions.resize(bond_count_);
        if (bond_count_ == 0)
            return;
        if (rows_ > max_blocks_)
            throw std::length_error("NOSB-PD failure stencil exceeds the device grid");
        upload(failure_sites_, sites, "NOSB-PD failure sites");
        upload(live_, bond_live, "NOSB-PD live flags before failure screening");
        failure_decisions_.reserve(bond_count_);
        screenBondFailures<<<static_cast<unsigned>(rows_), threads_>>>(
            offsets_.get(), bonds_.get(), failure_sites_.get(), live_.get(), parameters, failure_decisions_.get());
        check(cudaGetLastError(), "screen NOSB-PD bond failures");
        download(decisions, failure_decisions_, "NOSB-PD bond decisions");
    }

    void assembleInternalForce(const std::vector<NosbPDForceSite> &sites,
                               const std::vector<std::uint8_t> &bond_live,
                               std::vector<ParticleVector3> &force_density) override
    {
        validate(sites.size(), bond_live.size());
        force_density.resize(rows_);
        if (rows_ == 0)
            return;

        const auto blocks = static_cast<unsigned>((3 * rows_ - 1) / threads_ + 1);
        upload(force_sites_, sites, "NOSB-PD force sites");
        upload(live_, bond_live, "NOSB-PD live bond flags before internal force");
        force_density_.reserve(rows_);
        assembleInternalBondForce<<<blocks, threads_>>>(offsets_.get(), bonds_.get(), force_sites_.get(),
            live_.get(), rows_, force_density_.get());
        check(cudaGetLastError(), "assemble NOSB-PD internal bond force");
        download(force_density, force_density_, "NOSB-PD internal force density");
    }

    void assembleStabilization(const std::vector<NosbPDCorrectionSite> &sites,
                               const std::vector<std::uint8_t> &bond_live,
                               NosbPDStabilizationParameters parameters,
                               std::vector<NosbPDCorrectionForce> &force) override
    {
        assembleCorrection(sites, bond_live, parameters, force);
    }

    void assembleDamping(const std::vector<NosbPDCorrectionSite> &sites,
                         const std::vector<std::uint8_t> &bond_live,
                         NosbPDDampingParameters parameters,
                         std::vector<NosbPDCorrectionForce> &force) override
    {
        assembleCorrection(sites, bond_live, parameters, force);
    }

  private:
    template <bool Reference>
    void accumulatePositionMoments(const NosbPDPositionSite *positions, const double *reference_volume,
                                   std::vector<ShapeHessian3> &moments)
    {
        const auto blocks = static_cast<unsigned>((3 * rows_ - 1) / threads_ + 1);
        position_moments_.reserve(rows_);
        assemblePositionMoments<Reference><<<blocks, threads_>>>(
            offsets_.get(), bonds_.get(), positions, reference_volume, rows_, position_moments_.get());
        check(cudaGetLastError(), "assemble NOSB-PD position moments");
        download(moments, position_moments_, "NOSB-PD position moments");
    }

    template <class Parameters>
    void assembleCorrection(const std::vector<NosbPDCorrectionSite> &sites,
                            const std::vector<std::uint8_t> &bond_live, Parameters parameters,
                            std::vector<NosbPDCorrectionForce> &force)
    {
        validate(sites.size(), bond_live.size());
        force.resize(rows_);
        if (rows_ == 0)
            return;

        const auto blocks = static_cast<unsigned>((rows_ - 1) / threads_ + 1);
        upload(correction_sites_, sites, "NOSB-PD correction sites");
        upload(live_, bond_live, "NOSB-PD live bond flags before force correction");
        correction_force_.reserve(rows_);
        assembleCorrectionForce<<<blocks, threads_>>>(offsets_.get(), bonds_.get(), correction_sites_.get(),
            live_.get(), rows_, parameters, correction_force_.get());
        check(cudaGetLastError(), "assemble NOSB-PD force correction");
        download(force, correction_force_, "NOSB-PD correction forces");
    }

    void validate(std::size_t sites, std::size_t bonds) const
    {
        if (!configured_ || sites != site_count_ || bonds != bond_count_)
            throw std::invalid_argument("NOSB-PD inputs do not match the reference stencil");
    }

    static constexpr unsigned threads_ = 128;
    std::size_t max_blocks_ = 0, rows_ = 0, site_count_ = 0, bond_count_ = 0;
    bool configured_ = false;
    DeviceBuffer<std::size_t> offsets_;
    DeviceBuffer<NosbPDBondGeometry> bonds_;
    DeviceBuffer<NosbPDKinematicSite> sites_;
    DeviceBuffer<NosbPDPositionSite> position_sites_;
    DeviceBuffer<double> reference_volume_;
    DeviceBuffer<ShapeHessian3> position_moments_;
    DeviceBuffer<NosbPDForceSite> force_sites_;
    DeviceBuffer<NosbPDBondFailureSite> failure_sites_;
    DeviceBuffer<NosbPDBondDecision> failure_decisions_;
    DeviceBuffer<ParticleVector3> force_density_;
    DeviceBuffer<NosbPDCorrectionSite> correction_sites_;
    DeviceBuffer<NosbPDCorrectionForce> correction_force_;
    DeviceBuffer<std::uint8_t> live_;
    DeviceBuffer<NosbPDDeformationTensors> tensors_;
};
}

std::unique_ptr<NosbPDBondBatch> makeCudaNosbPDBonds()
{
    return std::make_unique<CudaNosbPDBonds>();
}
}

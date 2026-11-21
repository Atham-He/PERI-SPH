#include "pd_nosb_elastic_batch.h"
#include "pd_nosb_plastic_batch.h"
#include "pd_drucker_prager_return.h"
#include "pd_hughes_winget.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;

template <int dimensions>
using Tensor = Eigen::Matrix<double, dimensions, dimensions>;

template <int dimensions>
__device__ Tensor<dimensions> load(const ShapeHessian3 &entries)
{
    Tensor<dimensions> tensor;
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            tensor(i, j) = entries.derivative[i][j];
    return tensor;
}

template <int dimensions>
__device__ ShapeHessian3 pack(const Tensor<dimensions> &tensor)
{
    ShapeHessian3 entries{};
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            entries.derivative[i][j] = tensor(i, j);
    return entries;
}

template <int dimensions>
__device__ bool allFinite(const Tensor<dimensions> &tensor)
{
    for (int i = 0; i < dimensions; ++i)
        for (int j = 0; j < dimensions; ++j)
            if (!isfinite(tensor(i, j)))
                return false;
    return true;
}

template <int dimensions>
__global__ void updateNosbPDElasticity(const NosbPDElasticSite *sites, std::size_t count,
                                       double lambda, double shear, NosbPDElasticResult *results)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &site = sites[i];
    NosbPDElasticResult result{};
    if (!site.active)
    {
        result.status = NosbPDElasticStatus::Inactive;
        results[i] = result;
        return;
    }

    const auto shape_inverse = load<dimensions>(site.shape_inverse);
    const Tensor<dimensions> full = load<dimensions>(site.full) * shape_inverse;
    const Tensor<dimensions> increment = load<dimensions>(site.increment) * shape_inverse;
    const Tensor<dimensions> midpoint = load<dimensions>(site.midpoint) * shape_inverse;
    result.deformation = pack(full);
    result.deformation_delta = pack(increment);
    result.deformation_half = pack(midpoint);

    if (!allFinite(full) || !allFinite(midpoint))
    {
        result.status = NosbPDElasticStatus::NonFiniteGradient;
        results[i] = result;
        return;
    }
    const double determinant = full.determinant(), midpoint_determinant = midpoint.determinant();
    if (!isfinite(determinant) || !isfinite(midpoint_determinant) ||
        determinant < nosb_minimum_deformation_determinant || midpoint_determinant < nosb_minimum_deformation_determinant)
    {
        result.status = NosbPDElasticStatus::InvalidDeterminant;
        results[i] = result;
        return;
    }

    const Tensor<dimensions> inverse = full.inverse(), inverse_half = midpoint.inverse();
    const Tensor<dimensions> gradient = increment * inverse_half;
    const Tensor<dimensions> stress = hughesWingetStress(gradient, load<dimensions>(site.stress), lambda, shear);
    const Tensor<dimensions> piola = determinant * inverse * stress;
    const Tensor<dimensions> piola_transpose = piola.transpose();
    const Tensor<dimensions> force_state = piola_transpose * shape_inverse;
    result.inverse = pack(inverse);
    result.inverse_half = pack(inverse_half);
    result.stress = pack(stress);
    result.piola = pack(piola_transpose);
    result.force_state = pack(force_state);
    result.status = NosbPDElasticStatus::Accepted;
    results[i] = result;
}

__device__ double j2Hardening(NosbPDPlasticParameters parameters, double alpha)
{
    const double argument = -parameters.hardening_rate * alpha;
    const double saturation = parameters.legacy_hardening ? 1.0 - exp(argument) : -expm1(argument);
    return parameters.hardening_linear * alpha + parameters.hardening_saturation * saturation;
}

struct J2Residual
{
    NosbPDPlasticParameters parameters;
    double norm, resistance, alpha, root_two_thirds;

    __device__ double operator()(double gamma) const
    {
        return norm - resistance * gamma - root_two_thirds *
            (parameters.yield + j2Hardening(parameters, alpha + root_two_thirds * gamma));
    }
};

template <int dimensions>
__global__ void updateNosbPDPlastic(const NosbPDPlasticSite *sites, std::size_t count, NosbPDPlasticParameters parameters,
                               double dt, NosbPDPlasticResult *results)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &site = sites[i];
    NosbPDPlasticResult result{};
    if (!site.active)
    {
        result.status = NosbPDPlasticStatus::Inactive;
        results[i] = result;
        return;
    }
    if (!parameters.drucker_prager && (!isfinite(site.alpha) || site.alpha < 0))
    {
        result.status = NosbPDPlasticStatus::InvalidHistory;
        results[i] = result;
        return;
    }

    using Matrix = Tensor<dimensions>;
    const Matrix gradient = load<dimensions>(site.gradient);
    const Matrix symmetric = (gradient + gradient.transpose()) * 0.5;
    const Matrix skew = (gradient - gradient.transpose()) * 0.5;
    const Matrix system = Matrix::Identity() - skew * 0.5;
    const Matrix rotation = Matrix::Identity() + system.inverse() * skew;
    const Matrix stress_increment =
        parameters.lambda * symmetric.trace() * Matrix::Identity() + 2.0 * parameters.shear * symmetric;
    Matrix stress = rotation * load<dimensions>(site.stress) * rotation.transpose() + stress_increment;
    const Matrix rotated_plastic = rotation * load<dimensions>(site.plastic_strain) * rotation.transpose();
    Matrix plastic = rotated_plastic;
    Matrix kinematic = load<dimensions>(site.kinematic);
    double alpha = site.alpha;
    if (parameters.drucker_prager)
    {
        const DruckerPragerHardening<double> hardening{
            parameters.hardening_linear, parameters.hardening_saturation, parameters.hardening_rate};
        const auto status = returnDruckerPrager(stress, plastic, alpha, parameters.shear, parameters.bulk,
                                               parameters.pressure_sensitivity, parameters.yield,
                                               hardening, parameters.yield_tolerance);
        if (status != PlasticMultiplierStatus::Accepted)
        {
            result.status = NosbPDPlasticStatus::ReturnFailed;
            result.return_status = status;
            results[i] = result;
            return;
        }
    }
    else
    {
        kinematic = rotation * kinematic * rotation.transpose();
        const Matrix backstress = parameters.kinematic * kinematic;
        const Matrix relative = stress - backstress;
        const Matrix deviator = relative - relative.trace() * (1.0 / dimensions) * Matrix::Identity();
        const double norm = deviator.norm();
        const double root_two_thirds = sqrt(2.0 / 3.0);
        const double yield_value = norm - root_two_thirds * (j2Hardening(parameters, site.alpha) + parameters.yield);

        if (yield_value > parameters.yield_tolerance)
        {
            if (parameters.viscosity > 0 && (!isfinite(dt) || dt <= 0))
            {
                result.status = NosbPDPlasticStatus::InvalidViscousStep;
                results[i] = result;
                return;
            }
            const double viscosity_resistance = parameters.viscosity > 0 ? parameters.viscosity / (3.0 * dt) : 0.0;
            const double flow_resistance = parameters.shear + parameters.hardening_linear / 3.0 + parameters.kinematic / 2.0;
            const double resistance = 2.0 * (flow_resistance + viscosity_resistance) - 2.0 * parameters.hardening_linear / 3.0;
            const J2Residual residual{parameters, norm, resistance, alpha, root_two_thirds};
            double gamma = 0;
            const double scale = norm > parameters.yield ? norm : parameters.yield;
            const auto status = bracketPlasticMultiplier(residual, 0.0, yield_value / resistance, scale, gamma);
            if (status != PlasticMultiplierStatus::Accepted)
            {
                result.status = NosbPDPlasticStatus::ReturnFailed;
                result.return_status = status;
                results[i] = result;
                return;
            }
            const Matrix direction = deviator / norm;
            stress -= 2.0 * parameters.shear * gamma * direction;
            plastic += gamma * direction;
            alpha += root_two_thirds * gamma;
            kinematic += gamma * direction;
        }
    }

    result.status = NosbPDPlasticStatus::Accepted;
    result.stress = pack(stress);
    result.plastic_strain = pack(plastic);
    result.rotated_plastic_strain = pack(rotated_plastic);
    result.kinematic = pack(kinematic);
    result.alpha = alpha;
    results[i] = result;
}

class CudaNosbPDPlastic final : public NosbPDPlasticBatch
{
  public:
    CudaNosbPDPlastic()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select NOSB-PD J2 device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read NOSB-PD J2 grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void update(const std::vector<NosbPDPlasticSite> &sites, NosbPDPlasticParameters parameters,
                double dt, int dimensions, std::vector<NosbPDPlasticResult> &results) override
    {
        if (dimensions != 2 && dimensions != 3)
            throw std::invalid_argument("NOSB-PD plastic requires two or three dimensions");
        const auto count = sites.size();
        results.resize(count);
        if (count == 0)
            return;
        constexpr unsigned threads = 128;
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("NOSB-PD plastic update exceeds the device grid");

        upload(sites_, sites, "NOSB-PD plastic sites");
        results_.reserve(count);
        if (dimensions == 2)
            updateNosbPDPlastic<2><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, parameters, dt, results_.get());
        else
            updateNosbPDPlastic<3><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, parameters, dt, results_.get());
        check(cudaGetLastError(), "update NOSB-PD J2");
        download(results, results_, "NOSB-PD plastic results");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<NosbPDPlasticSite> sites_;
    DeviceBuffer<NosbPDPlasticResult> results_;
};

class CudaNosbPDElastic final : public NosbPDElasticBatch
{
  public:
    CudaNosbPDElastic()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select NOSB-PD elastic device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read NOSB-PD elastic grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void update(const std::vector<NosbPDElasticSite> &sites, double lambda, double shear, int dimensions,
                std::vector<NosbPDElasticResult> &results) override
    {
        if (dimensions != 2 && dimensions != 3)
            throw std::invalid_argument("NOSB-PD elasticity requires two or three dimensions");
        results.resize(sites.size());
        if (sites.empty())
            return;

        constexpr unsigned threads = 128;
        const auto count = sites.size();
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("NOSB-PD elastic update exceeds the device grid");
        upload(sites_, sites, "NOSB-PD elastic sites");
        results_.reserve(count);
        if (dimensions == 2)
            updateNosbPDElasticity<2><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, lambda, shear, results_.get());
        else
            updateNosbPDElasticity<3><<<static_cast<unsigned>(blocks), threads>>>(sites_.get(), count, lambda, shear, results_.get());
        check(cudaGetLastError(), "update NOSB-PD elasticity");
        download(results, results_, "NOSB-PD elastic results");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<NosbPDElasticSite> sites_;
    DeviceBuffer<NosbPDElasticResult> results_;
};
}

std::unique_ptr<NosbPDElasticBatch> makeCudaNosbPDElastic()
{
    return std::make_unique<CudaNosbPDElastic>();
}
std::unique_ptr<NosbPDPlasticBatch> makeCudaNosbPDPlastic()
{
    return std::make_unique<CudaNosbPDPlastic>();
}

}

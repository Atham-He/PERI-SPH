#include "pd_nosb_bond_batch.h"
#include "pd_cuda_buffer.cuh"

#include <algorithm>

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
        double eta = bond.reference.component[axis];
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
            const double reference = bond.reference.component[k];
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
        const double from = projectBond(source.force_state, axis, bond.kernel_weight, bond.reference, 1.0);
        const double to = projectBond(target.force_state, axis, bond.kernel_weight, bond.reference, -1.0);
        force += (from - to) * target.volume;
    }
    force_density[i].component[axis] = force;
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

  private:
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
    DeviceBuffer<NosbPDForceSite> force_sites_;
    DeviceBuffer<ParticleVector3> force_density_;
    DeviceBuffer<std::uint8_t> live_;
    DeviceBuffer<NosbPDDeformationTensors> tensors_;
};
}

std::unique_ptr<NosbPDBondBatch> makeCudaNosbPDBonds()
{
    return std::make_unique<CudaNosbPDBonds>();
}
}

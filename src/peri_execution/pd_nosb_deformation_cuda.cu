#include "pd_nosb_deformation_batch.h"
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

class CudaNosbPDDeformation final : public NosbPDDeformationBatch
{
  public:
    CudaNosbPDDeformation()
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

    void evaluate(const std::vector<NosbPDKinematicSite> &sites, const std::vector<std::uint8_t> &bond_live,
                  double dt, std::vector<NosbPDDeformationTensors> &tensors) override
    {
        if (!configured_ || sites.size() != site_count_ || bond_live.size() != bond_count_)
            throw std::invalid_argument("NOSB-PD deformation inputs do not match the reference stencil");
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

  private:
    static constexpr unsigned threads_ = 128;
    std::size_t max_blocks_ = 0, rows_ = 0, site_count_ = 0, bond_count_ = 0;
    bool configured_ = false;
    DeviceBuffer<std::size_t> offsets_;
    DeviceBuffer<NosbPDBondGeometry> bonds_;
    DeviceBuffer<NosbPDKinematicSite> sites_;
    DeviceBuffer<std::uint8_t> live_;
    DeviceBuffer<NosbPDDeformationTensors> tensors_;
};
}

std::unique_ptr<NosbPDDeformationBatch> makeCudaNosbPDDeformation()
{
    return std::make_unique<CudaNosbPDDeformation>();
}
}

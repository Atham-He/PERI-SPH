#include "pd_fluid_traction_batch.h"
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

__device__ double dot(const ParticleVector3 &a, const ParticleVector3 &b)
{
    return a.component[0]*b.component[0]+(a.component[1]*b.component[1]+a.component[2]*b.component[2]);
}

__global__ void evaluateTraction(const FluidTractionSolid *solids, const FluidTractionFluid *fluids,
                                 const std::size_t *offsets, const FluidTractionNeighbor *neighbors,
                                 std::size_t count, double acceleration_factor, FluidTractionResult *output)
{
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i>=count)
        return;

    const auto &solid=solids[i];
    FluidTractionResult row{};
    for (auto n=offsets[i];n<offsets[i+1];++n)
    {
        const auto &neighbor=neighbors[n];
        const auto &fluid=fluids[neighbor.fluid];
        const double gradient=neighbor.kernel_derivative*fluid.volume;
        const double facing=dot(neighbor.direction,solid.normal);
        const double sign=facing<0 ? 1. : (facing>0 ? -1. : 0.);
        ParticleVector3 acceleration, normal, relative_velocity;
        for (int axis=0;axis<3;++axis)
        {
            acceleration.component[axis]=fluid.force_prior.component[axis]/fluid.mass-solid.average_acceleration.component[axis];
            normal.component[axis]=sign*solid.normal.component[axis];
            relative_velocity.component[axis]=fluid.velocity.component[axis]-solid.average_velocity.component[axis];
        }
        const double face_acceleration=dot(acceleration,neighbor.direction);
        const double wall_pressure=fluid.pressure+acceleration_factor*fluid.density*neighbor.distance*
            (0.>face_acceleration ? 0. : face_acceleration);
        const double velocity_jump=2.*dot(relative_velocity,normal);
        double pressure_jump=0.;
        if (fluid.acoustic_impedance!=0.)
        {
            const double scaled_jump=3.*(fluid.inverse_sound_speed*(velocity_jump>0. ? velocity_jump : 0.));
            const double limiter=scaled_jump<1. ? scaled_jump : 1.;
            pressure_jump=fluid.acoustic_impedance*velocity_jump*limiter;
        }
        for (int axis=0;axis<3;++axis)
        {
            row.pressure_force.component[axis]-=(pressure_jump*normal.component[axis]+
                (wall_pressure+fluid.pressure)*neighbor.direction.component[axis])*solid.volume*gradient;
            row.acoustic_force.component[axis]+=(fluid.divergence*2.)*solid.volume*gradient*neighbor.direction.component[axis];
        }
        row.pressure+=fluid.pressure*neighbor.kernel_value;
        row.kernel_weight+=neighbor.kernel_value;
    }
    output[i]=row;
}

class CudaFluidTraction final : public FluidTractionBatch
{
  public:
    CudaFluidTraction()
    {
        int device=0, maximum=0;
        check(cudaGetDevice(&device),"select fluid traction device");
        check(cudaDeviceGetAttribute(&maximum,cudaDevAttrMaxGridDimX,device),"read fluid traction grid limit");
        max_blocks_=static_cast<std::size_t>(maximum);
    }

    void evaluate(const std::vector<FluidTractionSolid> &solid,
                  const std::vector<FluidTractionFluid> &fluid,
                  const std::vector<std::size_t> &offsets,
                  const std::vector<FluidTractionNeighbor> &neighbors,
                  double acceleration_factor,
                  std::vector<FluidTractionResult> &result) override
    {
        if (offsets.size()!=solid.size()+1 || offsets.front()!=0 || offsets.back()!=neighbors.size() ||
            !std::is_sorted(offsets.begin(),offsets.end()) ||
            std::any_of(neighbors.begin(),neighbors.end(),[&](const auto &n){return n.fluid>=fluid.size();}))
            throw std::invalid_argument("Invalid PD fluid traction stencil");
        result.resize(solid.size());
        if (solid.empty())
            return;
        constexpr unsigned threads=128;
        const auto blocks=(solid.size()-1)/threads+1;
        if (blocks>max_blocks_)
            throw std::length_error("PD fluid traction batch exceeds the device grid");
        upload(solid_,solid,"fluid traction solid particles");
        upload(fluid_,fluid,"fluid traction fluid particles");
        upload(offsets_,offsets,"fluid traction offsets");
        upload(neighbors_,neighbors,"fluid traction neighbors");
        result_.reserve(solid.size());
        evaluateTraction<<<static_cast<unsigned>(blocks),threads>>>(solid_.get(),fluid_.get(),offsets_.get(),
            neighbors_.get(),solid.size(),acceleration_factor,result_.get());
        check(cudaGetLastError(),"evaluate PD fluid traction");
        download(result,result_,"PD fluid traction results");
    }

  private:
    std::size_t max_blocks_=0;
    DeviceBuffer<FluidTractionSolid> solid_;
    DeviceBuffer<FluidTractionFluid> fluid_;
    DeviceBuffer<std::size_t> offsets_;
    DeviceBuffer<FluidTractionNeighbor> neighbors_;
    DeviceBuffer<FluidTractionResult> result_;
};
}

std::unique_ptr<FluidTractionBatch> makeCudaFluidTraction()
{
    return std::make_unique<CudaFluidTraction>();
}
}

#pragma once

#include "dpc_point_batch.h"
#include "pd_bonded_interface_batch.h"
#include "pd_crack_projection.h"
#include "pd_crack_strength_batch.h"
#include "pd_finite_volume_biot_batch.h"
#include "pd_fluid_acoustic_batch.h"
#include "pd_fluid_cfl_batch.h"
#include "pd_fluid_traction_batch.h"
#include "pd_fluid_viscosity_batch.h"
#include "pd_gel_batch.h"
#include "pd_hertz_sphere_contact.h"
#include "pd_nosb_bond_batch.h"
#include "pd_nosb_elastic_batch.h"
#include "pd_nosb_plastic_batch.h"
#include "pd_solid_time_step.h"
#include "pd_solid_cfl_batch.h"
#include "pd_spatial_derivatives.h"

#include <string_view>

namespace SPH::peridynamics
{
enum class PeriExecutionMode
{
    Default,
    CPU,
    GPU
};

PeriExecutionMode parsePeriExecution(std::string_view input);
PeriExecutionMode takePeriExecutionOption(int &argc, char **argv);

class PeriExecution
{
  public:
    // Partial phase selection remains available to existing diagnostic runs.
    struct LegacyGpuSelection
    {
        bool crack_projection;
        bool bonded_interface;
    };

    explicit PeriExecution(PeriExecutionMode mode,
                           LegacyGpuSelection legacy = {false, false});

    static bool hasCudaBackend() noexcept;
    const char *backendName() const noexcept;
    bool bondedInterfaceOnDevice() const noexcept { return bonded_interface_; }

    std::unique_ptr<CrackProjectionBatch> makeCrackProjection() const;
    std::unique_ptr<CrackStrengthBatch> makeCrackStrength() const;
    std::unique_ptr<BondedInterfaceBatch> makeBondedInterface() const;
    std::unique_ptr<SurfacePlaneContactBatch> makeSurfacePlaneContact() const;
    std::unique_ptr<HertzSpherePairBatch> makeHertzSpherePairs() const;
    std::unique_ptr<SpatialDerivativeBatch> makeSpatialDerivatives() const;
    std::unique_ptr<FiniteVolumeBiotBatch> makeFiniteVolumeBiot() const;
    std::unique_ptr<dpc::PointBatch> makeDpcPoints(const dpc::Model &model) const;
    std::unique_ptr<SolidTimeStepBatch> makeSolidTimeStep(const std::vector<double> &mass) const;
    std::unique_ptr<PositionVerletBatch> makePositionVerlet() const;
    std::unique_ptr<NosbPDBondBatch> makeNosbPDBonds() const;
    std::unique_ptr<NosbPDElasticBatch> makeNosbPDElastic() const;
    std::unique_ptr<NosbPDPlasticBatch> makeNosbPDPlastic() const;
    std::unique_ptr<GelPointBatch> makeGelPoints() const;
    std::unique_ptr<FluidTractionBatch> makeFluidTraction() const;
    std::unique_ptr<FluidAcousticBatch> makeFluidAcoustic() const;
    std::unique_ptr<FluidViscosityBatch> makeFluidViscosity() const;
    std::unique_ptr<FluidCflBatch> makeFluidCfl() const;
    std::unique_ptr<PDCflBatch> makePDCfl() const;

  private:
    PeriExecutionMode mode_;
    bool crack_projection_;
    bool bonded_interface_;
};
}

#pragma once

#include "dpc_point_batch.h"
#include "pd_bonded_interface_batch.h"
#include "pd_crack_projection.h"
#include "pd_fluid_acoustic_batch.h"
#include "pd_fluid_traction_batch.h"
#include "pd_solid_time_step.h"
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

    const char *backendName() const noexcept;
    bool bondedInterfaceOnDevice() const noexcept { return bonded_interface_; }

    std::unique_ptr<CrackProjectionBatch> makeCrackProjection() const;
    std::unique_ptr<BondedInterfaceBatch> makeBondedInterface() const;
    std::unique_ptr<SpatialDerivativeBatch> makeSpatialDerivatives() const;
    std::unique_ptr<dpc::PointBatch> makeDpcPoints(const dpc::Model &model) const;
    std::unique_ptr<SolidTimeStepBatch> makeSolidTimeStep(const std::vector<double> &mass) const;
    std::unique_ptr<FluidTractionBatch> makeFluidTraction() const;
    std::unique_ptr<FluidAcousticBatch> makeFluidAcoustic() const;

  private:
    PeriExecutionMode mode_;
    bool crack_projection_;
    bool bonded_interface_;
};
}

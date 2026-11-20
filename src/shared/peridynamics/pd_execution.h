#pragma once

#include "pd_bonded_interface_batch.h"
#include "pd_crack_projection.h"

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

  private:
    PeriExecutionMode mode_;
    bool crack_projection_;
    bool bonded_interface_;
};
}

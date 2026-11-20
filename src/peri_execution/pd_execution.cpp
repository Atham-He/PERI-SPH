#include "pd_execution.h"

#include <stdexcept>

namespace SPH::peridynamics
{
PeriExecutionMode parsePeriExecution(std::string_view input)
{
    if (input == "cpu")
        return PeriExecutionMode::CPU;
    if (input == "gpu")
        return PeriExecutionMode::GPU;
    throw std::invalid_argument("PERI-SPH execution must be cpu or gpu");
}

PeriExecution::PeriExecution(PeriExecutionMode mode, LegacyGpuSelection legacy)
    : mode_(mode), crack_projection_(mode == PeriExecutionMode::GPU || legacy.crack_projection),
      bonded_interface_(mode == PeriExecutionMode::GPU || legacy.bonded_interface)
{
    if (mode == PeriExecutionMode::CPU && (crack_projection_ || bonded_interface_))
        throw std::invalid_argument("execution=cpu conflicts with GPU diagnostic flags");
#if !defined(PERI_SPH_HAS_CUDA)
    if (crack_projection_ || bonded_interface_)
        throw std::invalid_argument("GPU execution requires a build with PERI_SPH_ENABLE_CUDA=ON");
#endif
}

const char *PeriExecution::backendName() const noexcept
{
    return crack_projection_ || bonded_interface_ ? "cuda" : "cpu";
}

std::unique_ptr<CrackProjectionBatch> PeriExecution::makeCrackProjection() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (crack_projection_)
        return makeCudaCrackProjection(mode_ == PeriExecutionMode::GPU);
#endif
    return nullptr;
}

std::unique_ptr<BondedInterfaceBatch> PeriExecution::makeBondedInterface() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (bonded_interface_)
        return makeCudaBondedInterface();
#endif
    return nullptr;
}

std::unique_ptr<SpatialDerivativeBatch> PeriExecution::makeSpatialDerivatives() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaSpatialDerivatives();
#endif
    return nullptr;
}

}

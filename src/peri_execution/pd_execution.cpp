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

PeriExecutionMode takePeriExecutionOption(int &argc, char **argv)
{
    PeriExecutionMode mode=PeriExecutionMode::Default;
    int remaining=1;
    for (int i=1;i<argc;++i)
    {
        const std::string_view argument(argv[i]);
        if (argument.substr(0,12)=="--execution=")
            mode=parsePeriExecution(argument.substr(12));
        else
            argv[remaining++]=argv[i];
    }
    argc=remaining;
    argv[argc]=nullptr;
    return mode;
}

PeriExecution::PeriExecution(PeriExecutionMode mode, LegacyGpuSelection legacy)
    : mode_(mode), crack_projection_(mode == PeriExecutionMode::GPU || legacy.crack_projection),
      bonded_interface_(mode == PeriExecutionMode::GPU || legacy.bonded_interface)
{
    if (mode == PeriExecutionMode::CPU && (crack_projection_ || bonded_interface_))
        throw std::invalid_argument("execution=cpu conflicts with GPU diagnostic flags");
    if ((crack_projection_ || bonded_interface_) && !hasCudaBackend())
        throw std::invalid_argument("GPU execution requires a build with PERI_SPH_ENABLE_CUDA=ON");
}

bool PeriExecution::hasCudaBackend() noexcept
{
#if defined(PERI_SPH_HAS_CUDA)
    return true;
#else
    return false;
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

std::unique_ptr<dpc::PointBatch> PeriExecution::makeDpcPoints(const dpc::Model &model) const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return dpc::makeCudaPointBatch(model);
#else
    (void)model;
#endif
    return nullptr;
}

std::unique_ptr<SolidTimeStepBatch> PeriExecution::makeSolidTimeStep(const std::vector<double> &mass) const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaSolidTimeStep(mass);
#else
    (void)mass;
#endif
    return nullptr;
}

std::unique_ptr<NosbPDBondBatch> PeriExecution::makeNosbPDBonds() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaNosbPDBonds();
#endif
    return nullptr;
}

std::unique_ptr<NosbPDElasticBatch> PeriExecution::makeNosbPDElastic() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaNosbPDElastic();
#endif
    return nullptr;
}

std::unique_ptr<NosbPDJ2Batch> PeriExecution::makeNosbPDJ2() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaNosbPDJ2();
#endif
    return nullptr;
}

std::unique_ptr<PositionVerletBatch> PeriExecution::makePositionVerlet() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaPositionVerlet();
#endif
    return nullptr;
}

std::unique_ptr<FluidTractionBatch> PeriExecution::makeFluidTraction() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaFluidTraction();
#endif
    return nullptr;
}

std::unique_ptr<FluidAcousticBatch> PeriExecution::makeFluidAcoustic() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaFluidAcoustic();
#endif
    return nullptr;
}

std::unique_ptr<FluidViscosityBatch> PeriExecution::makeFluidViscosity() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaFluidViscosity();
#endif
    return nullptr;
}

std::unique_ptr<FluidCflBatch> PeriExecution::makeFluidCfl() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaFluidCfl();
#endif
    return nullptr;
}
std::unique_ptr<PDCflBatch> PeriExecution::makePDCfl() const
{
#if defined(PERI_SPH_HAS_CUDA)
    if (mode_ == PeriExecutionMode::GPU)
        return makeCudaPDCfl();
#endif
    return nullptr;
}

}

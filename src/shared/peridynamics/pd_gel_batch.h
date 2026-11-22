#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
enum class GelKernelStatus
{
    OK, InvalidDeformation, InvalidHistory, InvalidDissipation,
    InvalidIncrement, InvalidRelaxation, InvalidDissipationIntegral, NonfiniteResponse
};

struct GelPointParameters
{
    double bulk, equilibrium_shear, branch_shear, relaxation_time;
};

struct GelPointState
{
    ShapeHessian3 deformation, green, viscous;
    double dissipation;
};

struct GelPointResult
{
    ShapeHessian3 green, viscous, pk1, cauchy;
    double dissipation, free_energy;
    GelKernelStatus status;
};

class GelPointBatch
{
  public:
    virtual ~GelPointBatch() = default;
    virtual void evaluate(const std::vector<GelPointState> &points, GelPointParameters parameters,
                          double dt, std::vector<GelPointResult> &result) = 0;
};

std::unique_ptr<GelPointBatch> makeCudaGelPoints();
}

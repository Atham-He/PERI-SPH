#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
struct FiniteVolumeBiotPointResult
{
    ShapeHessian3 piola;
    double jacobian, elastic_energy, tangent_bound;
    bool accepted;
};

class FiniteVolumeBiotBatch
{
  public:
    virtual ~FiniteVolumeBiotBatch() = default;
    virtual void evaluate(const std::vector<ShapeHessian3> &deformation, double bulk, double shear,
                          std::vector<FiniteVolumeBiotPointResult> &result) = 0;
};

std::unique_ptr<FiniteVolumeBiotBatch> makeCudaFiniteVolumeBiot();
}

#pragma once

#include "pd_spatial_derivatives.h"

namespace SPH::peridynamics
{
inline constexpr double nosb_minimum_deformation_determinant = 0.05;

enum class NosbPDElasticStatus
{
    Accepted,
    Inactive,
    NonFiniteGradient,
    InvalidDeterminant
};

struct NosbPDElasticSite
{
    ShapeHessian3 shape_inverse, full, increment, midpoint, stress;
    bool active;
};

struct NosbPDElasticResult
{
    ShapeHessian3 deformation, deformation_delta, deformation_half;
    ShapeHessian3 inverse, inverse_half, stress, piola, force_state;
    NosbPDElasticStatus status;
};

class NosbPDElasticBatch
{
  public:
    virtual ~NosbPDElasticBatch() = default;
    virtual void update(const std::vector<NosbPDElasticSite> &sites,
                        double lambda, double shear, int dimensions,
                        std::vector<NosbPDElasticResult> &results) = 0;
};

std::unique_ptr<NosbPDElasticBatch> makeCudaNosbPDElastic();
}

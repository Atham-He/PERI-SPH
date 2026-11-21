#pragma once

#include "pd_nosb_elastic_batch.h"
#include "pd_plastic_multiplier_status.h"

namespace SPH::peridynamics
{
struct NosbPDJ2Parameters
{
    double lambda, shear, yield, hardening_linear, hardening_saturation, hardening_rate;
    double kinematic, viscosity, yield_tolerance;
    bool legacy_hardening;
};

struct NosbPDJ2Site
{
    ShapeHessian3 gradient, stress, plastic_strain, kinematic;
    double alpha;
    bool active;
};

enum class NosbPDJ2Status { Accepted, Inactive, InvalidHistory, InvalidViscousStep, ReturnFailed };

struct NosbPDJ2Result
{
    ShapeHessian3 stress, plastic_strain, rotated_plastic_strain, kinematic;
    double alpha;
    NosbPDJ2Status status;
    PlasticMultiplierStatus return_status;
};

class NosbPDJ2Batch
{
  public:
    virtual ~NosbPDJ2Batch() = default;
    virtual void update(const std::vector<NosbPDJ2Site> &sites, NosbPDJ2Parameters parameters,
                        double dt, int dimensions, std::vector<NosbPDJ2Result> &results) = 0;
};

std::unique_ptr<NosbPDJ2Batch> makeCudaNosbPDJ2();
}

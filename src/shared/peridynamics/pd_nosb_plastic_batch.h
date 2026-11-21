#pragma once

#include "pd_nosb_elastic_batch.h"
#include "pd_plastic_multiplier_status.h"

namespace SPH::peridynamics
{
struct NosbPDPlasticParameters
{
    double lambda, shear, yield, hardening_linear, hardening_saturation, hardening_rate;
    double kinematic, viscosity, yield_tolerance;
    bool legacy_hardening;
    double pressure_sensitivity = 0, bulk = 0;
    bool drucker_prager = false;
};

struct NosbPDPlasticSite
{
    ShapeHessian3 gradient, stress, plastic_strain, kinematic;
    double alpha;
    bool active;
};

enum class NosbPDPlasticStatus { Accepted, Inactive, InvalidHistory, InvalidViscousStep, ReturnFailed };

struct NosbPDPlasticResult
{
    ShapeHessian3 stress, plastic_strain, rotated_plastic_strain, kinematic;
    double alpha;
    NosbPDPlasticStatus status;
    PlasticMultiplierStatus return_status;
};

class NosbPDPlasticBatch
{
  public:
    virtual ~NosbPDPlasticBatch() = default;
    virtual void update(const std::vector<NosbPDPlasticSite> &sites, NosbPDPlasticParameters parameters,
                        double dt, int dimensions, std::vector<NosbPDPlasticResult> &results) = 0;
};

std::unique_ptr<NosbPDPlasticBatch> makeCudaNosbPDPlastic();
}

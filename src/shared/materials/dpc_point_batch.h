#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace SPH::dpc
{
struct Node
{
    double r, E, nu, d, beta, R, pb;
};

struct PointInputModel
{
    const Node *nodes;
    std::size_t count;
    double initial_density, shear_dilation, shear_hardening;
};

enum class PointFailure
{
    None,
    DensityDomain,
    ShearHistory,
    HardeningDomain,
    Deformation,
    ElasticMetric,
    HydrostaticBound,
    ReturnDomain,
    ZeroDilationApex,
    ApexFlow,
    BadReturn
};

struct PointHistory
{
    double inverse_Cp[3][3];
    double kappa, alpha;
};

struct PointInput
{
    double deformation[3][3];
    PointHistory history;
};

struct PointOutput
{
    double sigma[3][3];
    PointHistory history;
    double p, q, relative_density, plastic_shear;
    double plastic_stress_work, elastic_energy, yield_residual, effective_d;
    int surface;
    PointFailure failure;
};

class PointBatch
{
  public:
    virtual ~PointBatch() = default;
    // Outputs trial material histories. The caller owns acceptance and commit.
    // Numerical fields are populated for successful points and BadReturn diagnostics.
    virtual void evaluate(const std::vector<PointInput> &input, std::vector<PointOutput> &output) = 0;
};

class Model;
std::unique_ptr<PointBatch> makeCudaPointBatch(const Model &model);
}

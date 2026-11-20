#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

namespace SPH::peridynamics
{
struct CrackProjectionInput
{
    // Four master corners, followed by the slave point.
    double coordinates[5][3];
};

struct CrackProjectionResult
{
    double s, t;
    int status;
};

// Quantities at a selected slave point and its four master traces.
struct CrackLocalContactInput
{
    double tangent_u[3], tangent_v[3], mixed_tangent[3], normal[3];
    double shape[4], derivative_u[4], derivative_v[4];
    double mass_weight[5];
    double gap, stiffness;
};

struct CrackLocalContactResult
{
    double row_bound[5][3], force_on_slave[3];
    bool valid;
};

class CrackProjectionBatch
{
  public:
    virtual ~CrackProjectionBatch() = default;
    virtual void project(const std::vector<CrackProjectionInput> &input,
                         std::vector<CrackProjectionResult> &output) = 0;
    virtual bool assemblesLocalContact() const noexcept { return false; }
    virtual void assembleLocalContact(const std::vector<CrackLocalContactInput> &,
                                      std::vector<CrackLocalContactResult> &)
    {
        throw std::logic_error("This crack backend only implements projection");
    }
};

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection(bool assemble_local_contact = false);
}

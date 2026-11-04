/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_GEL_MATERIAL_H
#define PD_GEL_MATERIAL_H

#include "pd_materials.h"
#include <cmath>

namespace SPH::peridynamics
{
/** Reference-volume free energy: isochoric neo-Hookean equilibrium + convex
 * volumetric energy + one Green-strain Maxwell branch. This is an effective
 * recoverable gel model; plasticity, damage and solvent transport are absent. */
struct GelParameters
{
    Real bulk, equilibrium_shear, branch_shear, relaxation_time;
    static GelParameters fromYoungPoisson(Real young, Real poisson, Real shear_ratio, Real tau);
    void validate() const;
    Real instantaneousYoung() const;
    Real instantaneousPoisson() const;
};

struct GelHistory
{
    Matd green = Matd::Zero(); // Accepted deviatoric Green strain.
    Matd viscous = Matd::Zero(); // Symmetric deviatoric reference internal strain.
    Real dissipation = 0; // Accumulated J / reference m^3.
};

struct GelResponse
{
    Matd pk1, cauchy;
    Real free_energy;
};

Matd gelDeviatoricGreen(const Matd &deformation);
GelResponse gelResponse(const Matd &deformation, const GelParameters &parameters, const GelHistory &history);
/** Exact Maxwell evolution and dissipation for a linear segment in deviatoric
 * Green strain. dt=0 leaves every history value unchanged. */
GelHistory advanceGel(const Matd &deformation, Real dt, const GelParameters &parameters, const GelHistory &old);
} // namespace SPH::peridynamics

namespace SPH
{
class GelViscoelasticSolid : public HughesWingetSolid
{
  public:
    GelViscoelasticSolid(Real density, Real equilibrium_young, Real equilibrium_poisson,
                        Real relaxing_shear_ratio, Real relaxation_time);
    void initializeLocalParameters(BaseParticles *particles) override;
    Real PairNumericalDamping(Real, Real) override { return 0; }
    const peridynamics::GelParameters &gelParameters() const { return parameters_; }

  private:
    peridynamics::GelParameters parameters_;
};
} // namespace SPH
#endif

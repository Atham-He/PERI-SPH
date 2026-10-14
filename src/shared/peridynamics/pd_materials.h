/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_materials.h
 * @brief 	Materials for peridynamic bodies.
 * @details	Ported from Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU).
 * 			The Hughes-Winget incremental constitutive update supports
 * 			finite strain and finite rotation within the NOSB-PD framework.
 */

#ifndef PD_MATERIALS_H
#define PD_MATERIALS_H

#include "elastic_solid.h"

namespace SPH
{
/**
 * @class HughesWingetSolid
 * @brief Incremental linear-elastic constitutive model for NOSB-PD
 * with finite strain and finite rotation (Hughes & Winget, 1980).
 * @details This class registers all particle variables required by the
 * NOSB-PD algorithm. The actual incremental update is performed by the
 * peridynamics dynamics classes in nosbpd_dynamics.h.
 */
class HughesWingetSolid : public LinearElasticSolid
{
  public:
    explicit HughesWingetSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
    virtual ~HughesWingetSolid() {};

    virtual void registerLocalParametersToReload(BaseParticles *base_particles) override;
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    /** Hughes-Winget incremental update of the Cauchy stress
     * with the Eulerian velocity gradient G (finite rotation objective). */
    virtual Matd StressHW(Matd &velocity_gradient, Matd &stress_old, size_t particle_index_i);
    /** Kelvin-Voigt type damping stress for bond-pair interaction. */
    virtual Real PairNumericalDamping(Real strain_rate, Real smoothing_length);
    /** bulk modulus used by the Littlewood hourglass control. */
    virtual Real getBulkModulusforHGC(const int &ndim);
};

/**
 * @class PlasticSolidforPD
 * @brief Abstract elastoplastic solid for PD with combined
 * isotropic-kinematic hardening (rate-independent).
 */
class PlasticSolidforPD : public HughesWingetSolid
{
  public:
    explicit PlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                               Real yield_stress, Real isotropic_hardening_modulus = 0.0,
                               Real kinematic_hardening_modulus = 0.0);
    virtual ~PlasticSolidforPD() {};

    virtual void registerLocalParametersToReload(BaseParticles *base_particles) override;
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    Real YieldStress() { return yield_stress_; };
    Real IsotropicHardeningModulus() { return isotropic_hardening_modulus_; };
    Real KinematicHardeningModulus() { return kinematic_hardening_modulus_; };

    /** yield function value */
    virtual Real YieldFunc(const Matd &stress, const Real &iso_q = 0.0, const Matd &kin_q = Matd::Zero()) = 0;
    /** convex cutting-plane return mapping with Hughes-Winget objectivity */
    virtual Matd PlasticConstitutiveRelation(const Matd &vel_grad, const Matd &stress_old,
                                             Matd &epsilon_p, size_t index_i, Real dt = 0.0) = 0;

  protected:
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    Real yield_stress_;
    Real isotropic_hardening_modulus_;
    Real kinematic_hardening_modulus_;
    /** internal state variables registered as particle variables */
    Real *iso_hardening_ISV1_;
    Matd *kin_hardening_ISV2_;
};

/**
 * @class J2PlasticityforPD
 * @brief Associated J2 plasticity for PD with isotropic-kinematic hardening.
 */
class J2PlasticityforPD : public PlasticSolidforPD
{
  public:
    explicit J2PlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
                               Real isotropic_hardening_modulus = 0.0, Real kinematic_hardening_modulus = 0.0);
    virtual ~J2PlasticityforPD() {};

    virtual Real YieldFunc(const Matd &stress, const Real &iso_q = 0.0, const Matd &kin_q = Matd::Zero()) override;
    inline Real IsoHardeningFunc(const Real &ISV1)
    {
        Real exponent = exp(-index_kesi_ * ISV1);
        return isotropic_hardening_modulus_ * ISV1 + (sigmaY_infinite_ - sigmaY_0_) * (1.0 - exponent);
    };
    inline Matd KinHardeningFunc(const Matd &ISV2) { return kinematic_hardening_modulus_ * ISV2; };
    virtual Matd PlasticConstitutiveRelation(const Matd &vel_grad, const Matd &stress_old,
                                             Matd &epsilon_p, size_t index_i, Real dt = 0.0) override;

  protected:
    /** denominator of the J2 plastic corrector (flow resistance).
     * Rate-independent J2: G + H_iso/3 + H_kin/2. Rate-dependent
     * extensions (e.g. Perzyna) add a viscosity term. */
    virtual Real J2FlowResistance(Real dt) const;
    Real sigmaY_infinite_, sigmaY_0_, index_kesi_;
};

/**
 * @class ViscoPlasticSolidforPD
 * @brief Perzyna-type visco-plastic solid for PD (rate-dependent J2).
 * @details A scalar bracketed solve enforces the nonlinear consistency relation
 * with saturation hardening. Viscosity adds 2*eta*gamma/(3*dt) to the
 * J2 overstress residual. The kinematic variable is defined by backstress
 * = H_kin * internal_strain, giving resistance G + H_iso/3 + H_kin/2
 * in the linear-hardening limit.
 * The model reduces to rate-independent J2 as viscosity -> 0 and to
 * purely elastic response as viscosity -> inf; the characteristic
 * overstress relaxation time is eta/(3G). Suitable for alginate-type
 * hydrogel shells with rate-dependent yielding.
 */
class ViscoPlasticSolidforPD : public J2PlasticityforPD
{
  public:
    explicit ViscoPlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                    Real yield_stress, Real viscosity,
                                    Real isotropic_hardening_modulus = 0.0,
                                    Real kinematic_hardening_modulus = 0.0);
    virtual ~ViscoPlasticSolidforPD() {};

    Real Viscosity() { return viscosity_; };

  protected:
    virtual Real J2FlowResistance(Real dt) const override;
    Real viscosity_; ///< Perzyna viscosity parameter [Pa s]
};

/**
 * @class DruckerPragerPlasticityforPD
 * @brief Associated Drucker-Prager plasticity for PD, isotropic hardening only.
 * Suitable for pressure-sensitive quasi-brittle materials (ice, rock, sand).
 */
class DruckerPragerPlasticityforPD : public PlasticSolidforPD
{
  public:
    explicit DruckerPragerPlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                          Real yield_stress, Real alpha_phi,
                                          Real isotropic_hardening_modulus = 0.0,
                                          Real kinematic_hardening_modulus = 0.0);
    virtual ~DruckerPragerPlasticityforPD() {};

    virtual Real YieldFunc(const Matd &stress, const Real &iso_q = 0.0, const Matd &kin_q = Matd::Zero()) override;
    inline Real IsoHardeningFunc(const Real &ISV1)
    {
        Real exponent = exp(-index_kesi_ * ISV1);
        return isotropic_hardening_modulus_ * ISV1 + (sigmaY_infinite_ - sigmaY_0_) * (1.0 - exponent);
    };
    inline Matd KinHardeningFunc(const Matd &ISV2) { return kinematic_hardening_modulus_ * ISV2; };
    virtual Matd PlasticConstitutiveRelation(const Matd &vel_grad, const Matd &stress_old,
                                             Matd &epsilon_p, size_t index_i, Real dt = 0.0) override;

  protected:
    const Real sqrt2_ = sqrt(2.0);
    Real alpha_phi_;
    Real k_phi_;
    Real sigmaY_infinite_, sigmaY_0_, index_kesi_;
};
} // namespace SPH
#endif // PD_MATERIALS_H

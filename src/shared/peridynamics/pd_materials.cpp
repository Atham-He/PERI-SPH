/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_materials.cpp
 * @brief 	Materials for peridynamic bodies.
 */

#include "pd_materials.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
HughesWingetSolid::HughesWingetSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
    : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio)
{
    material_type_name_ = "HughesWingetSolid";
}
//=================================================================================================//
void HughesWingetSolid::registerLocalParametersToReload(BaseParticles *base_particles)
{
    /** NOTE: keep this EMPTY on purpose. The framework calls this method
     * BEFORE initializeLocalParameters even in non-reload mode, and the
     * variable registry keeps the first registration (with default zero
     * initialization), which would override the physical initial values
     * (e.g. ParticleLive = 1, identity matrices) assigned in
     * initializeLocalParameters. Restart registration of the PD variables
     * is left for future work. */
}
//=================================================================================================//
void HughesWingetSolid::registerLocalParametersFromReload(BaseParticles *base_particles)
{
    /** see the note in registerLocalParametersToReload. */
}
//=================================================================================================//
void HughesWingetSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    LinearElasticSolid::initializeLocalParameters(base_particles);

    /** bond state: 1 = particle effective, 0 = free particle (all bonds broken) */
    base_particles->registerStateVariableData<int>("ParticleLive", [&](size_t i) -> int
                                                   { return 1; });
    base_particles->registerStateVariableData<int>("BondNumber", [&](size_t i) -> int
                                                   { return 0; });
    base_particles->registerStateVariableData<Real>("Damage", [&](size_t i) -> Real
                                                    { return 0.0; });
    /** shape tensor and its inverse in the reference configuration */
    base_particles->registerStateVariableData<Matd>("ShapeTensorK");
    base_particles->registerStateVariableData<Matd>("ShapeTensorKInv", IdentityMatrix<Matd>::value);
    /** variables for the Hughes-Winget incremental algorithm */
    base_particles->registerStateVariableData<Matd>("DeformationTensorN", IdentityMatrix<Matd>::value);
    base_particles->registerStateVariableData<Matd>("VelocityTensorNdeltaU");
    base_particles->registerStateVariableData<Matd>("DeformationTensorNHalf", IdentityMatrix<Matd>::value);
    base_particles->registerStateVariableData<Matd>("DeformationGradientHalf", IdentityMatrix<Matd>::value);
    base_particles->registerStateVariableData<Matd>("VelocityGradientFdelta");
    base_particles->registerStateVariableData<Matd>("DeformationGradientInv", IdentityMatrix<Matd>::value);
    base_particles->registerStateVariableData<Matd>("DeformationGradientHalfInv", IdentityMatrix<Matd>::value);
    /** Cauchy stress, first Piola-Kirchhoff stress and peridynamic force state */
    base_particles->registerStateVariableData<Matd>("PDStress");
    base_particles->registerStateVariableData<Matd>("PDPK1");
    base_particles->registerStateVariableData<Matd>("PDForceStateT0");
    /** generic solid variables required by standard reduce classes, e.g. AcousticTimeStep */
    base_particles->registerStateVariableData<Vecd>("Force");
    base_particles->registerStateVariableData<Vecd>("ForcePrior");
    /** wall normal for FSI coupling, initialized by NormalDirectionFromBodyShape */
    base_particles->registerStateVariableData<Vecd>("NormalDirection");

    /** variables evolving with the motion, e.g. for particle sorting */
    base_particles->addEvolvingVariable<int>("ParticleLive");
    base_particles->addEvolvingVariable<int>("BondNumber");
    base_particles->addEvolvingVariable<Real>("Damage");
    base_particles->addEvolvingVariable<Matd>("PDStress");
}
//=================================================================================================//
Matd HughesWingetSolid::StressHW(Matd &G, Matd &stress_old, size_t particle_index_i)
{
    /** symmetric part of G: rate of deformation tensor */
    Matd Gsymm = (G + G.transpose()) * 0.5;
    /** antisymmetric (skew) part of G: rate of rotation tensor */
    Matd Gskew = (G - G.transpose()) * 0.5;

    Matd delta_sigma = lambda0_ * Gsymm.trace() * Matd::Identity() + 2.0 * G0_ * Gsymm;

    Matd R = Matd::Identity() - Gskew * 0.5;
    Matd Q = Matd::Identity() + R.inverse() * Gskew;

    Matd trial_sigma = Q * stress_old * Q.transpose();

    return delta_sigma + trial_sigma;
}
//=================================================================================================//
Real HughesWingetSolid::PairNumericalDamping(Real strain_rate, Real smoothing_length)
{
    return 0.5 * rho0_ * c0_ * strain_rate * smoothing_length;
}
//=================================================================================================//
Real HughesWingetSolid::getBulkModulusforHGC(const int &ndim)
{
    return ndim == 2 ? E0_ / 2.0 / (1.0 - nu_) : K0_;
}
//=================================================================================================//
PlasticSolidforPD::
    PlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                      Real yield_stress, Real isotropic_hardening_modulus,
                      Real kinematic_hardening_modulus)
    : HughesWingetSolid(rho0, youngs_modulus, poisson_ratio),
      yield_stress_(yield_stress),
      isotropic_hardening_modulus_(isotropic_hardening_modulus),
      kinematic_hardening_modulus_(kinematic_hardening_modulus),
      iso_hardening_ISV1_(nullptr), kin_hardening_ISV2_(nullptr)
{
    material_type_name_ = "PlasticSolidforPD";
}
//=================================================================================================//
void PlasticSolidforPD::registerLocalParametersToReload(BaseParticles *base_particles)
{
    HughesWingetSolid::registerLocalParametersToReload(base_particles);
}
//=================================================================================================//
void PlasticSolidforPD::registerLocalParametersFromReload(BaseParticles *base_particles)
{
    HughesWingetSolid::registerLocalParametersFromReload(base_particles);
}
//=================================================================================================//
void PlasticSolidforPD::initializeLocalParameters(BaseParticles *base_particles)
{
    HughesWingetSolid::initializeLocalParameters(base_particles);
    base_particles->registerStateVariableData<Matd>("PlasticStrain");
    iso_hardening_ISV1_ = base_particles->registerStateVariableData<Real>("IsotropicHardeningParam");
    kin_hardening_ISV2_ = base_particles->registerStateVariableData<Matd>("KinematicHardeningParam");

    base_particles->addEvolvingVariable<Matd>("PlasticStrain");
    base_particles->addEvolvingVariable<Real>("IsotropicHardeningParam");
    base_particles->addEvolvingVariable<Matd>("KinematicHardeningParam");
}
//=================================================================================================//
J2PlasticityforPD::
    J2PlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
                      Real isotropic_hardening_modulus, Real kinematic_hardening_modulus)
    : PlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                        isotropic_hardening_modulus, kinematic_hardening_modulus),
      sigmaY_infinite_(1.589 * yield_stress), sigmaY_0_(yield_stress), index_kesi_(16.93)
{
    material_type_name_ = "J2PlasticityforPD";
}
//=================================================================================================//
Real J2PlasticityforPD::YieldFunc(const Matd &stress, const Real &iso_q, const Matd &kin_q)
{
    return 0.0;
}
//=================================================================================================//
Real J2PlasticityforPD::J2FlowResistance(Real dt) const
{
    return G0_ + (isotropic_hardening_modulus_ + kinematic_hardening_modulus_) / 3.0;
}
//=================================================================================================//
Matd J2PlasticityforPD::PlasticConstitutiveRelation(const Matd &G, const Matd &stress_old_const,
                                                    Matd &epsilon_p, size_t index_i, Real dt)
{
    Matd stress_old = stress_old_const;
    /** Hughes-Winget incremental objectivity */
    Matd Gsymm = (G + G.transpose()) * 0.5;
    Matd Gskew = (G - G.transpose()) * 0.5;
    Matd R = Matd::Identity() - Gskew * 0.5;
    Matd Q = Matd::Identity() + R.inverse() * Gskew;

    /** elastic predictor */
    Matd delta_sigma = lambda0_ * Gsymm.trace() * Matd::Identity() + 2.0 * G0_ * Gsymm;

    Matd trial_sigma = Q * stress_old * Q.transpose() + delta_sigma;
    Matd trial_epsilon_p = Q * epsilon_p * Q.transpose();
    /** predict the internal state variables */
    Real trial_isoHardening_ISV1 = iso_hardening_ISV1_[index_i];
    Matd trial_kinHardening_ISV2 = Q * kin_hardening_ISV2_[index_i] * Q.transpose();
    /** hardening variables */
    Real trial_isoHardening_q = IsoHardeningFunc(trial_isoHardening_ISV1);
    Matd trial_kinHardening_q = KinHardeningFunc(trial_kinHardening_ISV2);

    /** intermediate variables */
    Matd trial_relative_stress = trial_sigma - trial_kinHardening_q;
    Real norm1 = trial_relative_stress.trace() * OneOverDimensions;
    Matd dev_eta = trial_relative_stress - norm1 * Matd::Identity();
    Real dev_eta_norm = dev_eta.norm();

    /** yield function */
    Real trial_function = dev_eta_norm - sqrt_2_over_3_ * (trial_isoHardening_q + yield_stress_);

    if (trial_function > TinyReal)
    {
        /** plastic corrector */
        Real relax_increment = 0.5 * trial_function / J2FlowResistance(dt);
        Matd flow_direction = dev_eta / dev_eta_norm;
        trial_sigma -= 2.0 * G0_ * relax_increment * flow_direction;
        trial_epsilon_p += relax_increment * flow_direction;
        trial_isoHardening_ISV1 += sqrt_2_over_3_ * relax_increment;
        trial_kinHardening_ISV2 += relax_increment * flow_direction;
    }

    epsilon_p = trial_epsilon_p;
    iso_hardening_ISV1_[index_i] = trial_isoHardening_ISV1;
    kin_hardening_ISV2_[index_i] = trial_kinHardening_ISV2;

    return trial_sigma;
}
//=================================================================================================//
ViscoPlasticSolidforPD::
    ViscoPlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                           Real yield_stress, Real viscosity,
                           Real isotropic_hardening_modulus, Real kinematic_hardening_modulus)
    : J2PlasticityforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                        isotropic_hardening_modulus, kinematic_hardening_modulus),
      viscosity_(viscosity)
{
    material_type_name_ = "ViscoPlasticSolidforPD";
}
//=================================================================================================//
Real ViscoPlasticSolidforPD::J2FlowResistance(Real dt) const
{
    /** Perzyna overstress, backward Euler: d_gamma = 0.5 f / (G + H/3 + eta/(3 dt)).
     * dt <= 0 (rate information unavailable) falls back to rate-independent J2. */
    Real viscous_resistance = (viscosity_ > 0.0 && dt > 0.0) ? viscosity_ / (3.0 * dt) : 0.0;
    return J2PlasticityforPD::J2FlowResistance(dt) + viscous_resistance;
}
//=================================================================================================//
DruckerPragerPlasticityforPD::
    DruckerPragerPlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                 Real yield_stress, Real alpha_phi,
                                 Real isotropic_hardening_modulus, Real kinematic_hardening_modulus)
    : PlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                        isotropic_hardening_modulus, kinematic_hardening_modulus),
      k_phi_(yield_stress), alpha_phi_(alpha_phi),
      sigmaY_infinite_(1.589 * yield_stress), sigmaY_0_(yield_stress), index_kesi_(16.93)
{
    material_type_name_ = "DruckerPragerPlasticityforPD";
}
//=================================================================================================//
Real DruckerPragerPlasticityforPD::YieldFunc(const Matd &stress, const Real &iso_q, const Matd &kin_q)
{
    return 0.0;
}
//=================================================================================================//
Matd DruckerPragerPlasticityforPD::PlasticConstitutiveRelation(const Matd &G, const Matd &stress_old_const,
                                                               Matd &epsilon_p, size_t index_i, Real dt)
{
    Matd stress_old = stress_old_const;
    /** Hughes-Winget incremental objectivity */
    Matd Gsymm = (G + G.transpose()) * 0.5;
    Matd Gskew = (G - G.transpose()) * 0.5;
    Matd R = Matd::Identity() - Gskew * 0.5;
    Matd Q = Matd::Identity() + R.inverse() * Gskew;

    /** elastic predictor */
    Matd delta_sigma = lambda0_ * Gsymm.trace() * Matd::Identity() + 2.0 * G0_ * Gsymm;

    Matd trial_sigma = Q * stress_old * Q.transpose() + delta_sigma;
    Matd trial_epsilon_p = Q * epsilon_p * Q.transpose();
    /** predict the internal state variables, hardening in the "k" level */
    Real trial_isoHardening_ISV1 = iso_hardening_ISV1_[index_i];
    Real trial_isoHardening_q = IsoHardeningFunc(trial_isoHardening_ISV1);
    /** intermediate variables */
    Real norm1 = trial_sigma.trace();
    Matd dev_sigma = trial_sigma - norm1 * OneOverDimensions * Matd::Identity();
    Real dev_sigma_norm = dev_sigma.norm();

    /** yield function */
    Real trial_function = alpha_phi_ * norm1 + dev_sigma_norm / sqrt2_ -
                          (trial_isoHardening_q + k_phi_);

    if (trial_function > TinyReal)
    {
        /** plastic corrector */
        Real relax_increment = trial_function /
                               (9.0 * K0_ * alpha_phi_ * alpha_phi_ + G0_ + isotropic_hardening_modulus_);
        Matd eta = dev_sigma / dev_sigma_norm;
        Matd flow_direction = alpha_phi_ * Matd::Identity() + eta / sqrt2_;
        trial_sigma -= relax_increment * (3.0 * K0_ * alpha_phi_ * Matd::Identity() + sqrt2_ * G0_ * eta);
        trial_epsilon_p += relax_increment * flow_direction;
        trial_isoHardening_ISV1 += relax_increment;
    }

    epsilon_p = trial_epsilon_p;
    iso_hardening_ISV1_[index_i] = trial_isoHardening_ISV1;

    return trial_sigma;
}
//=================================================================================================//
} // namespace SPH

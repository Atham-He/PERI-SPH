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
#include "pd_numerics.h"
#include <functional>

namespace SPH
{
namespace
{
Real solvePlasticMultiplier(const std::function<Real(Real)> &residual, Real lower, Real upper, Real scale)
{
    for (int i = 0; residual(upper) > 0.0 && i < 64; ++i) upper = 2.0 * upper + 1e-16;
    if (!std::isfinite(residual(upper)) || residual(upper) > 0.0)
        throw std::runtime_error("PD constitutive return mapping: failed to bracket plastic multiplier");
    for (int i = 0; i < 100; ++i)
    {
        Real mid = 0.5 * (lower + upper);
        Real value = residual(mid);
        if (!std::isfinite(value)) throw std::runtime_error("PD constitutive return mapping: non-finite residual");
        if (std::abs(value) <= 1e-11 * std::max(scale, Real(1))) return mid;
        if (value > 0.0) lower = mid; else upper = mid;
    }
    throw std::runtime_error("PD constitutive return mapping: local solve did not converge");
}
} // namespace
//=================================================================================================//
HughesWingetSolid::HughesWingetSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
    : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio)
{
    if (!std::isfinite(rho0) || rho0 <= 0 || !std::isfinite(youngs_modulus) || youngs_modulus <= 0 ||
        !std::isfinite(poisson_ratio) || poisson_ratio <= -1 || poisson_ratio >= 0.5)
        throw std::invalid_argument("HughesWingetSolid: positive density/modulus and -1 < Poisson ratio < 0.5 required");
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

    base_particles->registerStateVariableData<int>("NumericalFailureReason", 0);
    base_particles->addVariableToWrite<int>("NumericalFailureReason");
    base_particles->addEvolvingVariable<int>("NumericalFailureReason");
    for (const char *name : {"PDInternalForce", "PDNumericalDampingForce", "HourglassForce",
                             "ContactFrictionForce", "ViscousForceFromFluid", "AcousticDampingForceFromFluid"})
        base_particles->registerStateVariableData<Vecd>(name);
    for (const char *name : {"PDInternalWork", "PDExternalWork", "PDHourglassWork", "PDDampingWork",
                             "ShellViscousWork", "ContactBoundaryWork", "ContactFrictionLoss",
                             "ContactBoundaryPower", "ContactFrictionBoundaryPower", "PlasticWork", "ShellAcousticWallWork"})
        base_particles->registerStateVariableData<Real>(name);
    /** bond state: 1 = particle effective, 0 = free particle (all bonds broken) */
    base_particles->registerStateVariableData<int>("ParticleLive", [&](size_t i) -> int
                                                   { return 1; });
    base_particles->registerStateVariableData<int>("BondNumber", [&](size_t i) -> int
                                                   { return 0; });
    base_particles->registerStateVariableData<Real>("Damage", [&](size_t i) -> Real
                                                    { return 0.0; });
    base_particles->addVariableToWrite<int>("ParticleLive");
    base_particles->addVariableToWrite<int>("BondNumber");
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
    if (!std::isfinite(yield_stress) || yield_stress <= 0.0 ||
        !std::isfinite(isotropic_hardening_modulus) || isotropic_hardening_modulus < 0.0 ||
        !std::isfinite(kinematic_hardening_modulus) || kinematic_hardening_modulus < 0.0)
        throw std::invalid_argument("PlasticSolidforPD: yield stress must be positive and hardening non-negative");
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
    : J2PlasticityforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                       PDIsotropicHardening::Legacy(yield_stress, isotropic_hardening_modulus),
                       kinematic_hardening_modulus)
{}
//=================================================================================================//
J2PlasticityforPD::
    J2PlasticityforPD(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
                      const PDIsotropicHardening &hardening, Real kinematic_hardening_modulus)
    : PlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                        hardening.linearModulus(), kinematic_hardening_modulus),
      hardening_(hardening)
{
    material_type_name_ = "J2PlasticityforPD";
}
//=================================================================================================//
Real J2PlasticityforPD::YieldFunc(const Matd &stress, const Real &iso_q, const Matd &kin_q)
{
    Matd relative = stress - kin_q;
    Matd deviator = relative - relative.trace() * OneOverDimensions * Matd::Identity();
    return deviator.norm() - sqrt_2_over_3_ * (yield_stress_ + iso_q);
}
//=================================================================================================//
Real J2PlasticityforPD::J2FlowResistance(Real dt) const
{
    return G0_ + isotropic_hardening_modulus_ / 3.0 + kinematic_hardening_modulus_ / 2.0;
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
        // q(alpha) contains the complete isotropic law, including H*alpha.
        // Remove that linear contribution here so it is counted exactly once.
        Real linear_resistance = 2.0 * J2FlowResistance(dt) - 2.0 * isotropic_hardening_modulus_ / 3.0;
        auto residual = [&](Real gamma)
        {
            return dev_eta_norm - linear_resistance * gamma - sqrt_2_over_3_ *
                   (yield_stress_ + IsoHardeningFunc(trial_isoHardening_ISV1 + sqrt_2_over_3_ * gamma));
        };
        Real relax_increment = solvePlasticMultiplier(residual, 0.0, trial_function / linear_resistance,
                                                       std::max(dev_eta_norm, yield_stress_));
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
    : ViscoPlasticSolidforPD(rho0, youngs_modulus, poisson_ratio, yield_stress, viscosity,
                            PDIsotropicHardening::Legacy(yield_stress, isotropic_hardening_modulus),
                            kinematic_hardening_modulus)
{}
//=================================================================================================//
ViscoPlasticSolidforPD::
    ViscoPlasticSolidforPD(Real rho0, Real youngs_modulus, Real poisson_ratio,
                           Real yield_stress, Real viscosity, const PDIsotropicHardening &hardening,
                           Real kinematic_hardening_modulus)
    : J2PlasticityforPD(rho0, youngs_modulus, poisson_ratio, yield_stress,
                        hardening, kinematic_hardening_modulus),
      viscosity_(viscosity)
{
    if (!std::isfinite(viscosity) || viscosity < 0.0)
        throw std::invalid_argument("ViscoPlasticSolidforPD: viscosity must be finite and non-negative");
    material_type_name_ = "ViscoPlasticSolidforPD";
}
//=================================================================================================//
Real ViscoPlasticSolidforPD::J2FlowResistance(Real dt) const
{
    /** The viscous resistance enters the nonlinear consistency residual.
     * A rate-dependent plastic update always requires a positive time increment. */
    if (viscosity_ > 0.0 && (!std::isfinite(dt) || dt <= 0.0))
        throw std::invalid_argument("ViscoPlasticSolidforPD requires a positive time increment");
    Real viscous_resistance = viscosity_ > 0.0 ? viscosity_ / (3.0 * dt) : 0.0;
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
    if (!std::isfinite(alpha_phi) || alpha_phi < 0.0 || kinematic_hardening_modulus != 0.0)
        throw std::invalid_argument("DruckerPragerPlasticityforPD: non-negative alpha and isotropic hardening only");
    material_type_name_ = "DruckerPragerPlasticityforPD";
}
//=================================================================================================//
Real DruckerPragerPlasticityforPD::YieldFunc(const Matd &stress, const Real &iso_q, const Matd &kin_q)
{
    Matd deviator = stress - stress.trace() * OneOverDimensions * Matd::Identity();
    return alpha_phi_ * stress.trace() + deviator.norm() / sqrt2_ - (k_phi_ + iso_q);
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
        Real q_trial = dev_sigma_norm / sqrt2_;
        Real gamma_apex = q_trial / G0_;
        auto smooth_residual = [&](Real gamma)
        {
            return alpha_phi_ * norm1 + q_trial - (9.0 * K0_ * alpha_phi_ * alpha_phi_ + G0_) * gamma -
                   k_phi_ - IsoHardeningFunc(trial_isoHardening_ISV1 + gamma);
        };
        Real relax_increment;
        if (smooth_residual(gamma_apex) <= 0.0)
        {
            relax_increment = solvePlasticMultiplier(smooth_residual, 0.0, gamma_apex,
                                                      std::max(std::abs(alpha_phi_ * norm1) + q_trial, k_phi_));
            Matd eta = dev_sigma / dev_sigma_norm;
            trial_sigma -= relax_increment * (3.0 * K0_ * alpha_phi_ * Matd::Identity() + sqrt2_ * G0_ * eta);
            trial_epsilon_p += relax_increment * (alpha_phi_ * Matd::Identity() + eta / sqrt2_);
        }
        else
        {
            // Cone apex: deviatoric stress is exhausted. Continue on the hydrostatic branch.
            auto apex_residual = [&](Real gamma)
            {
                return alpha_phi_ * norm1 - 9.0 * K0_ * alpha_phi_ * alpha_phi_ * gamma -
                       k_phi_ - IsoHardeningFunc(trial_isoHardening_ISV1 + gamma);
            };
            Real upper = gamma_apex + trial_function / (9.0 * K0_ * alpha_phi_ * alpha_phi_);
            relax_increment = solvePlasticMultiplier(apex_residual, gamma_apex, upper,
                                                      std::max(std::abs(alpha_phi_ * norm1), k_phi_));
            trial_sigma = (norm1 * OneOverDimensions - 3.0 * K0_ * alpha_phi_ * relax_increment) * Matd::Identity();
            trial_epsilon_p += dev_sigma / (2.0 * G0_) + alpha_phi_ * relax_increment * Matd::Identity();
        }
        trial_isoHardening_ISV1 += relax_increment;
    }

    epsilon_p = trial_epsilon_p;
    iso_hardening_ISV1_[index_i] = trial_isoHardening_ISV1;

    return trial_sigma;
}
//=================================================================================================//
} // namespace SPH

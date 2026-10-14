/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	nosbpd_dynamics.cpp
 * @brief 	Non-ordinary state-based peridynamics (NOSB-PD) dynamics.
 */

#include "nosbpd_dynamics.h"

#include "adaptation.h"
#include "base_particles.hpp"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
PDTimeStepInitialization::
    PDTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : LocalDynamics(sph_body), gravity_ptr_(gravity_ptr),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      acc_prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")) {}
//=================================================================================================//
void PDTimeStepInitialization::update(size_t index_i, Real dt)
{
    /** gravity applies to all particles, including fully debonded (free) ones:
     * shell fragments remain ordinary debris with mass. */
    acc_prior_[index_i] = gravity_ptr_->InducedAcceleration(pos_[index_i]);
}
//=================================================================================================//
NosbPDShapeMatrix::
    NosbPDShapeMatrix(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      shape_K_(particles_->getVariableDataByName<Matd>("ShapeTensorK")),
      shape_K_1_(particles_->getVariableDataByName<Matd>("ShapeTensorKInv")) {}
//=================================================================================================//
void NosbPDShapeMatrix::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        /** reference bond vector r_j - r_i */
        Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration += inner_neighborhood.W_ij_[n] * Vol_[index_j] * (r_ij * r_ij.transpose());
    }
    shape_K_[index_i] += local_configuration;
}
//=================================================================================================//
void NosbPDShapeMatrix::update(size_t index_i, Real dt)
{
    shape_K_1_[index_i] = shape_K_[index_i].inverse();
}
//=================================================================================================//
NosbPDFirstStep::NosbPDFirstStep(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, sph_body.getMatterMaterial())),
      rho_(particles_->registerStateVariableData<Real>("Density", hw_solid_.ReferenceDensity())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_old_(particles_->registerStateVariableData<Vecd>("PDAccelerationOld")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      acc_prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value))
{
    rho0_ = hw_solid_.ReferenceDensity();
}
//=================================================================================================//
void NosbPDFirstStep::update(size_t index_i, Real dt)
{
    acc_old_[index_i] = acc_[index_i];
    acc_[index_i] = acc_prior_[index_i];
    rho_[index_i] = rho0_ / F_[index_i].determinant();
    /** NaN fuse: a particle with non-finite velocity is disabled and frozen
     * before its position is poisoned and the contamination spreads through
     * the bond network (and into the fluid via AverageVelocity). */
    if (!vel_[index_i].allFinite())
    {
        std::cout << "Particle_index = " << index_i << " disabled (non-finite velocity)\n";
        particleLive_[index_i] = 0;
        vel_[index_i].setZero();
        acc_[index_i].setZero();
        return;
    }
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
NosbPDSecondStep::
    NosbPDSecondStep(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      shape_K_1_(particles_->getVariableDataByName<Matd>("ShapeTensorKInv")),
      N_(particles_->getVariableDataByName<Matd>("DeformationTensorN")),
      N_deltaU_(particles_->getVariableDataByName<Matd>("VelocityTensorNdeltaU")),
      N_half_(particles_->getVariableDataByName<Matd>("DeformationTensorNHalf")),
      F_half_(particles_->getVariableDataByName<Matd>("DeformationGradientHalf")),
      F_delta_(particles_->getVariableDataByName<Matd>("VelocityGradientFdelta")),
      F_1_(particles_->getVariableDataByName<Matd>("DeformationGradientInv")),
      F_1_half_(particles_->getVariableDataByName<Matd>("DeformationGradientHalfInv")),
      PK1_(particles_->getVariableDataByName<Matd>("PDPK1")),
      T0_(particles_->getVariableDataByName<Matd>("PDForceStateT0")),
      stress_(particles_->getVariableDataByName<Matd>("PDStress")) {}
//=================================================================================================//
void NosbPDSecondStep::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        Matd N = Matd::Zero();
        Matd N_deltaU = Matd::Zero();
        Matd N_half = Matd::Zero();

        Vecd eta = Vecd::Zero();
        Vecd eta_U = Vecd::Zero();
        Vecd eta_half = Vecd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            /** reference bond vector r_j - r_i */
            Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                eta = pos_[index_j] - pos_[index_i];
                eta_U = (vel_[index_j] - vel_[index_i]) * dt;
                eta_half = eta - eta_U * 0.5;
            }
            else
            {
                eta = r_ij;
                eta_U = Vecd::Zero();
                eta_half = r_ij;
            }
            Matd dN = eta * r_ij.transpose();
            Matd dN_deltaU = eta_U * r_ij.transpose();
            Matd dN_half = eta_half * r_ij.transpose();
            N += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN;
            N_deltaU += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_deltaU;
            N_half += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_half;
        }
        N_[index_i] = N;
        N_deltaU_[index_i] = N_deltaU;
        N_half_[index_i] = N_half;
    }
}
//=================================================================================================//
void NosbPDSecondStep::update(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        F_[index_i] = N_[index_i] * shape_K_1_[index_i];
        F_delta_[index_i] = N_deltaU_[index_i] * shape_K_1_[index_i];
        F_half_[index_i] = N_half_[index_i] * shape_K_1_[index_i];

        /** NaN isolation: a non-finite F marks a blown-up particle; disable it
         * and zero its force state so the contamination cannot propagate to
         * the neighbors through the bond forces. */
        if (!F_[index_i].allFinite() || !F_half_[index_i].allFinite())
        {
            std::cout << "Particle_index = " << index_i << " disabled (non-finite F)\n";
            particleLive_[index_i] = 0;
            stress_[index_i] = Matd::Zero();
            T0_[index_i] = Matd::Zero();
            return;
        }
        Real detF = F_[index_i].determinant();
        Real detF_half = F_half_[index_i].determinant();
        /** guard BEFORE inverting: a nearly singular F (volume collapsed to a
         * few percent, beyond the resolution limit) produces an exploding
         * inverse and NaN contamination through the bond network. */
        if (detF < 0.05 || detF_half < 0.05)
        {
            /** disable the particle instead of aborting the whole run */
            std::cout << "Particle_index = " << index_i << " has been disabled since det F < 0.05"
                      << "\n dt: " << dt << "\n";
            particleLive_[index_i] = 0;
            stress_[index_i] = Matd::Zero();
            T0_[index_i] = Matd::Zero();
            return;
        }
        F_1_[index_i] = F_[index_i].inverse();
        F_1_half_[index_i] = F_half_[index_i].inverse();
        /** G: Eulerian velocity gradient tensor */
        Matd G = F_delta_[index_i] * F_1_half_[index_i];
        /** update the Cauchy stress by the Hughes-Winget algorithm */
        stress_[index_i] = hw_solid_.StressHW(G, stress_[index_i], index_i);
        /** PK1: the first Piola-Kirchhoff stress tensor */
        Matd PK1 = detF * F_1_[index_i] * stress_[index_i];
        PK1_[index_i] = PK1.transpose();
        /** T0_: peridynamic force state */
        T0_[index_i] = PK1_[index_i] * shape_K_1_[index_i];
    }
    else
    {
        stress_[index_i] = Matd::Zero();
    }
}
//=================================================================================================//
NosbPDSecondStepPlastic::
    NosbPDSecondStepPlastic(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      plastic_solid_(DynamicCast<PlasticSolidforPD>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      shape_K_1_(particles_->getVariableDataByName<Matd>("ShapeTensorKInv")),
      N_(particles_->getVariableDataByName<Matd>("DeformationTensorN")),
      N_deltaU_(particles_->getVariableDataByName<Matd>("VelocityTensorNdeltaU")),
      N_half_(particles_->getVariableDataByName<Matd>("DeformationTensorNHalf")),
      F_half_(particles_->getVariableDataByName<Matd>("DeformationGradientHalf")),
      F_delta_(particles_->getVariableDataByName<Matd>("VelocityGradientFdelta")),
      F_1_(particles_->getVariableDataByName<Matd>("DeformationGradientInv")),
      F_1_half_(particles_->getVariableDataByName<Matd>("DeformationGradientHalfInv")),
      PK1_(particles_->getVariableDataByName<Matd>("PDPK1")),
      T0_(particles_->getVariableDataByName<Matd>("PDForceStateT0")),
      stress_(particles_->getVariableDataByName<Matd>("PDStress")),
      plastic_strain_(particles_->registerStateVariableData<Matd>("PlasticStrain")) {}
//=================================================================================================//
void NosbPDSecondStepPlastic::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        Matd N = Matd::Zero();
        Matd N_deltaU = Matd::Zero();
        Matd N_half = Matd::Zero();

        Vecd eta = Vecd::Zero();
        Vecd eta_U = Vecd::Zero();
        Vecd eta_half = Vecd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            /** reference bond vector r_j - r_i */
            Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                eta = pos_[index_j] - pos_[index_i];
                eta_U = (vel_[index_j] - vel_[index_i]) * dt;
                eta_half = eta - eta_U * 0.5;
            }
            else
            {
                eta = r_ij;
                eta_U = Vecd::Zero();
                eta_half = r_ij;
            }
            Matd dN = eta * r_ij.transpose();
            Matd dN_deltaU = eta_U * r_ij.transpose();
            Matd dN_half = eta_half * r_ij.transpose();
            N += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN;
            N_deltaU += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_deltaU;
            N_half += Vol_[index_j] * inner_neighborhood.W_ij_[n] * dN_half;
        }
        N_[index_i] = N;
        N_deltaU_[index_i] = N_deltaU;
        N_half_[index_i] = N_half;
    }
}
//=================================================================================================//
void NosbPDSecondStepPlastic::update(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        F_[index_i] = N_[index_i] * shape_K_1_[index_i];
        F_delta_[index_i] = N_deltaU_[index_i] * shape_K_1_[index_i];
        F_half_[index_i] = N_half_[index_i] * shape_K_1_[index_i];

        /** NaN isolation: see NosbPDSecondStep::update. */
        if (!F_[index_i].allFinite() || !F_half_[index_i].allFinite())
        {
            std::cout << "Particle_index = " << index_i << " disabled (non-finite F)\n";
            particleLive_[index_i] = 0;
            stress_[index_i] = Matd::Zero();
            T0_[index_i] = Matd::Zero();
            return;
        }
        Real detF = F_[index_i].determinant();
        Real detF_half = F_half_[index_i].determinant();
        /** guard BEFORE inverting: see NosbPDSecondStep::update. */
        if (detF < 0.05 || detF_half < 0.05)
        {
            /** disable the particle instead of aborting the whole run */
            std::cout << "Particle_index = " << index_i << " has been disabled since det F < 0.05"
                      << "\n dt: " << dt << "\n";
            particleLive_[index_i] = 0;
            stress_[index_i] = Matd::Zero();
            T0_[index_i] = Matd::Zero();
            return;
        }
        F_1_[index_i] = F_[index_i].inverse();
        F_1_half_[index_i] = F_half_[index_i].inverse();
        /** G: Eulerian velocity gradient tensor */
        Matd G = F_delta_[index_i] * F_1_half_[index_i];
        /** return-mapping (convex cutting-plane) with Hughes-Winget objectivity;
         * dt passed so that rate-dependent (visco-plastic) materials can
         * evaluate the overstress term. */
        stress_[index_i] = plastic_solid_.
            PlasticConstitutiveRelation(G, stress_[index_i], plastic_strain_[index_i], index_i, dt);
        /** PK1: the first Piola-Kirchhoff stress tensor */
        Matd PK1 = detF * F_1_[index_i] * stress_[index_i];
        PK1_[index_i] = PK1.transpose();
        /** T0_: peridynamic force state */
        T0_[index_i] = PK1_[index_i] * shape_K_1_[index_i];
    }
    else
    {
        stress_[index_i] = Matd::Zero();
    }
}
//=================================================================================================//
NosbPDCheckBondLive::
    NosbPDCheckBondLive(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      critical_value_(MaxReal),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")) {}
//=================================================================================================//
BondBreakByPrinStress::
    BondBreakByPrinStress(BaseInnerRelation &inner_relation, Real critical_value)
    : NosbPDCheckBondLive(inner_relation),
      bond_(particles_->getVariableDataByName<int>("BondNumber")),
      damage_(particles_->getVariableDataByName<Real>("Damage")),
      stress_(particles_->getVariableDataByName<Matd>("PDStress"))
{
    critical_value_ = critical_value;
}
//=================================================================================================//
bool BondBreakByPrinStress::checkBondLive(Matd &stress_eq, Real &stretch_rate)
{
    Real val_eq = stress_eq.eigenvalues().real().maxCoeff();
    return !(val_eq > critical_value_ || stretch_rate > 0.5);
}
//=================================================================================================//
void BondBreakByPrinStress::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        size_t bondcount = 0;

        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                Vecd eta = pos_[index_j] - pos_[index_i];
                Real stretchRate = (eta.norm() - inner_neighborhood.r_ij_[n]) / inner_neighborhood.r_ij_[n];
                Matd stress_avg = (stress_[index_i] + stress_[index_j]) * 0.5;

                inner_neighborhood.bondLive_[n] = checkBondLive(stress_avg, stretchRate);
                if (inner_neighborhood.bondLive_[n]) bondcount++;
            }
        }
        bond_[index_i] = int(bondcount);
        damage_[index_i] = 1.0 - Real(bondcount) / Real(inner_neighborhood.current_size_);
        if (bondcount < 1)
        {
            particleLive_[index_i] = 0;
            std::cout << "Particle_index_i = " << index_i << " becomes a FREE Particle !" << "\n";
        }
    }
}
//=================================================================================================//
BondBreakBySigma1andSigma3::
    BondBreakBySigma1andSigma3(BaseInnerRelation &inner_relation, Real max_tension, Real max_shear)
    : BondBreakByPrinStress(inner_relation, max_tension),
      max_tension_(max_tension), max_shear_(max_shear) {}
//=================================================================================================//
bool BondBreakBySigma1andSigma3::checkBondLive(Matd &stress_eq, Real &stretch_rate)
{
    Vecd prin_stress = stress_eq.eigenvalues().real();
    Real sigma1 = prin_stress.maxCoeff();
    Real sigma3 = prin_stress.minCoeff();
    Real max_shear = (sigma1 - sigma3) / 2.0;

    return !(sigma1 > max_tension_ || max_shear > max_shear_ || stretch_rate > 0.5);
}
//=================================================================================================//
BondBreakBySigma1andNorm1::
    BondBreakBySigma1andNorm1(BaseInnerRelation &inner_relation, Real max_tension, Real max_pressure)
    : BondBreakByPrinStress(inner_relation, max_tension),
      max_tension_(max_tension), max_pressure_(max_pressure) {}
//=================================================================================================//
bool BondBreakBySigma1andNorm1::checkBondLive(Matd &stress_eq, Real &stretch_rate)
{
    Real sigma1 = stress_eq.eigenvalues().real().maxCoeff();
    Real pressure = -OneOverDimensions * stress_eq.trace();

    return !(sigma1 > max_tension_ || pressure > max_pressure_ || stretch_rate > 0.5);
}
//=================================================================================================//
BondBreakByPlasticStrain::
    BondBreakByPlasticStrain(BaseInnerRelation &inner_relation, Real critical_value)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      critical_value_(critical_value),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      bond_(particles_->getVariableDataByName<int>("BondNumber")),
      damage_(particles_->getVariableDataByName<Real>("Damage")),
      plastic_strain_(particles_->registerStateVariableData<Matd>("PlasticStrain")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")) {}
//=================================================================================================//
bool BondBreakByPlasticStrain::checkBondLive(Matd &plastic_strain_eq, Real &stretch_rate)
{
    Real val_eq = plastic_strain_eq.eigenvalues().real().maxCoeff();
    return !(val_eq > critical_value_ || stretch_rate > 1.0);
}
//=================================================================================================//
void BondBreakByPlasticStrain::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        size_t bondcount = 0;

        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                Vecd eta = pos_[index_j] - pos_[index_i];
                Real stretchRate = (eta.norm() - inner_neighborhood.r_ij_[n]) / inner_neighborhood.r_ij_[n];
                Matd plastic_strain_avg = (plastic_strain_[index_i] + plastic_strain_[index_j]) * 0.5;

                inner_neighborhood.bondLive_[n] = checkBondLive(plastic_strain_avg, stretchRate);
                if (inner_neighborhood.bondLive_[n]) bondcount++;
            }
        }
        bond_[index_i] = int(bondcount);
        damage_[index_i] = 1.0 - Real(bondcount) / Real(inner_neighborhood.current_size_);
        if (bondcount < 1)
        {
            particleLive_[index_i] = 0;
            std::cout << "Particle_index_i = " << index_i << " becomes a FREE Particle !" << "\n";
        }
    }
}
//=================================================================================================//
NosbPDThirdStep::
    NosbPDThirdStep(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      T0_(particles_->getVariableDataByName<Matd>("PDForceStateT0"))
{
    rho0_ = hw_solid_.ReferenceDensity();
}
//=================================================================================================//
void NosbPDThirdStep::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        Vecd acceleration = Vecd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                /** force state T0 projected onto the bond; T_ij: force from j to i */
                Vecd T_ij = inner_neighborhood.W_ij_[n] * T0_[index_i] * r_ij;
                Vecd T_ji = inner_neighborhood.W_ij_[n] * T0_[index_j] * -r_ij;
                acceleration += (T_ij - T_ji) * Vol_[index_j];
            }
        }
        acc_[index_i] += acceleration / rho0_;
    }
}
//=================================================================================================//
NosbPDFourthStep::NosbPDFourthStep(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")) {}
//=================================================================================================//
void NosbPDFourthStep::update(size_t index_i, Real dt)
{
    vel_[index_i] += acc_[index_i] * dt;
    /** NaN fuse (second fence): velocity is checked immediately after the
     * integration so a poisoned value can never reach the position or the
     * AverageVelocity read by the fluid side. */
    if (!vel_[index_i].allFinite())
    {
        std::cout << "Particle_index = " << index_i << " disabled (non-finite velocity)\n";
        particleLive_[index_i] = 0;
        vel_[index_i].setZero();
        acc_[index_i].setZero();
        return;
    }
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
LittleWoodHourGlassControl::
    LittleWoodHourGlassControl(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      damage_(particles_->getVariableDataByName<Real>("Damage")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value))
{
    rho0_ = hw_solid_.ReferenceDensity();
    Chg_ = 0.5;
    Kbulk_ = hw_solid_.getBulkModulusforHGC(Dimensions);
    horizon_ = sph_body_->getSPHAdaptation().getKernel()->CutOffRadius();
}
//=================================================================================================//
void LittleWoodHourGlassControl::interaction(size_t index_i, Real dt)
{
    /** skip damaged particles: the hourglass force measures the deviation
     * between the F-predicted and the actual neighbor positions. In damaged
     * (crushed) zones F is no longer representative, so the deviation is
     * physical (fracture) rather than a zero-energy mode — and the
     * horizon^-3 stiffness would blow the velocity field up. */
    if (particleLive_[index_i] == 1 && damage_[index_i] < 0.1)
    {
        Real bond_const = 12.0 * Kbulk_ / Pi / pow(horizon_, 3);
        Vecd dhg_force = Vecd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            /** reference bond vector r_j - r_i */
            Vecd r_ij = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                Vecd r_ij_predict = F_[index_i] * r_ij;
                Vecd pos_j_predict = pos_[index_i] + r_ij_predict;
                Vecd hour_dir = pos_j_predict - pos_[index_j];
                Vecd eta = pos_[index_j] - pos_[index_i];
                Real hproj = hour_dir.transpose() * eta;

                dhg_force += Chg_ * bond_const * hproj / r_ij.norm() / eta.norm() * Vol_[index_j] * eta;
            }
        }
        acc_[index_i] -= dhg_force / rho0_;
    }
}
//=================================================================================================//
PairNumericalDampingforPD::
    PairNumericalDampingforPD(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value))
{
    rho0_ = hw_solid_.ReferenceDensity();
    c0_ = hw_solid_.ReferenceSoundSpeed();
    numerical_dissipation_factor_ = 0.25;
    smoothing_length_ = sph_body_->getSPHAdaptation().ReferenceSmoothingLength();
    inv_W0_ = 1.0 / sph_body_->getSPHAdaptation().getKernel()->W0(Vecd());
}
//=================================================================================================//
void PairNumericalDampingforPD::interaction(size_t index_i, Real dt)
{
    if (particleLive_[index_i] == 1)
    {
        Vecd damping_force = Vecd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd e_ij = inner_neighborhood.e_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            if (inner_neighborhood.bondLive_[n])
            {
                Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
                Real dim_r_ij_1 = Real(Dimensions) / r_ij;
                Vecd pos_jump = pos_[index_i] - pos_[index_j];
                Vecd vel_jump = vel_[index_i] - vel_[index_j];
                Real strain_rate = dim_r_ij_1 * dim_r_ij_1 * pos_jump.dot(vel_jump);
                Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
                Matd numerical_stress_ij =
                    0.5 * (F_[index_i] + F_[index_j]) * hw_solid_.PairNumericalDamping(strain_rate, smoothing_length_);
                damping_force += dW_ijV_j * (numerical_dissipation_factor_ * weight * numerical_stress_ij) * e_ij;
            }
        }
        acc_[index_i] += damping_force / rho0_;
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH

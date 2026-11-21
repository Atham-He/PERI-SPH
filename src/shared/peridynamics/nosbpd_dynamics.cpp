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
#include "pd_numerics.h"

#include "adaptation.h"
#include "base_particles.hpp"

namespace SPH
{
namespace solid_dynamics
{
namespace
{
Vecd principalValues(const Matd &tensor)
{
    if (!tensor.allFinite()) throw std::runtime_error("PD fracture criterion received a non-finite tensor");
    Eigen::SelfAdjointEigenSolver<Matd> solver;
    solver.computeDirect(0.5 * (tensor + tensor.transpose()), Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success) throw std::runtime_error("PD principal-value calculation failed");
    return solver.eigenvalues();
}
}

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
    Eigen::SelfAdjointEigenSolver<Matd> eig(shape_K_[index_i]);
    if (!shape_K_[index_i].allFinite() || eig.info() != Eigen::Success ||
        eig.eigenvalues().minCoeff() <= 1e-12 * eig.eigenvalues().maxCoeff())
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "singular ShapeTensorK");
    shape_K_1_[index_i] = shape_K_[index_i].inverse();
}
//=================================================================================================//
namespace
{
peridynamics::ParticleVector3 packPDVector(const Vecd &value)
{
    peridynamics::ParticleVector3 packed{};
    for (int axis = 0; axis < Dimensions; ++axis)
        packed.component[axis] = value[axis];
    return packed;
}

Vecd unpackPDVector(const peridynamics::ParticleVector3 &value)
{
    Vecd result;
    for (int axis = 0; axis < Dimensions; ++axis)
        result[axis] = value.component[axis];
    return result;
}

peridynamics::ShapeHessian3 packNosbPDTensor(const Matd &value)
{
    peridynamics::ShapeHessian3 tensor{};
    for (int i = 0; i < Dimensions; ++i)
        for (int j = 0; j < Dimensions; ++j)
            tensor.derivative[i][j] = value(i, j);
    return tensor;
}

Matd unpackNosbPDTensor(const peridynamics::ShapeHessian3 &tensor)
{
    Matd result;
    for (int i = 0; i < Dimensions; ++i)
        for (int j = 0; j < Dimensions; ++j)
            result(i, j) = tensor.derivative[i][j];
    return result;
}

void checkPositionVerletVelocity(SPHBody &body, BaseParticles &particles, std::size_t i, const Vecd &velocity)
{
    if (!velocity.allFinite())
    {
        particles.getVariableDataByName<int>("NumericalFailureReason")[i] = 1;
        peridynamics::numericalFailure(body.Name(), i, "non-finite Velocity");
    }
}
}

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
void NosbPDFirstStep::prepare(std::size_t i)
{
    acc_old_[i] = acc_[i];
    acc_[i] = acc_prior_[i];
    rho_[i] = rho0_ / F_[i].determinant();
    checkPositionVerletVelocity(*sph_body_, *particles_, i, vel_[i]);
}

void NosbPDFirstStep::update(size_t index_i, Real dt)
{
    prepare(index_i);
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}

void NosbPDFirstStep::evaluateDevice(peridynamics::PositionVerletBatch &backend, Real dt)
{
    const auto count = particles_->TotalRealParticles();
    device_input_.resize(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        prepare(i);
        device_input_[i] = {packPDVector(pos_[i]), packPDVector(vel_[i])};
    }
    backend.start(device_input_, dt, device_output_);
    for (std::size_t i = 0; i < count; ++i)
        pos_[i] = unpackPDVector(device_output_[i]);
}
//=================================================================================================//
void NosbPDBondWorkspace::configure(std::size_t count, std::size_t bound)
{
    using namespace peridynamics;
    const auto missing = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> ghost_indices(bound - count, missing);
    std::vector<NosbPDBondGeometry> bonds;
    offsets_.clear();
    ghost_sources_.clear();
    std::vector<std::size_t> offsets(1, 0);
    offsets.reserve(count + 1);
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto &row = relation_.inner_configuration_[i];
        for (std::size_t n = 0; n < row.current_size_; ++n)
        {
            const auto j = row.j_[n];
            auto target = j;
            if (j >= count)
            {
                auto &slot = ghost_indices.at(j - count);
                if (slot == missing)
                {
                    slot = count + ghost_sources_.size();
                    ghost_sources_.push_back(j);
                }
                target = slot;
            }
            bonds.push_back({target, packPDVector(row.e_ij_[n]), row.r_ij_[n], row.W_ij_[n], row.dW_ij_[n]});
        }
        offsets.push_back(bonds.size());
    }
    bond_live_.resize(bonds.size());
    backend_->configure(count + ghost_sources_.size(), offsets, bonds);
    offsets_ = std::move(offsets);
    bound_ = bound;
}

void NosbPDBondWorkspace::refreshLiveBonds(std::size_t count, std::size_t bound)
{
    // Pack the frozen geometry on first use, after checkpoint restore.
    if (offsets_.size() != count + 1 || bound_ != bound)
        configure(count, bound);
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto &row = relation_.inner_configuration_[i];
        if (row.current_size_ != offsets_[i + 1] - offsets_[i])
            throw std::logic_error("NOSB-PD reference row changed after bond workspace initialization");
        for (std::size_t n = 0; n < row.current_size_; ++n)
            bond_live_[offsets_[i] + n] = row.bondLive_[n];
    }
}

void NosbPDBondWorkspace::evaluateDeformation(Real dt, Matd *full, Matd *increment, Matd *midpoint)
{
    auto &particles = relation_.getSPHBody().getBaseParticles();
    const auto count = particles.TotalRealParticles();
    refreshLiveBonds(count, particles.ParticlesBound());
    const auto *position = particles.getVariableDataByName<Vecd>("Position");
    const auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    const auto *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const auto *active = particles.getVariableDataByName<int>("ParticleLive");

    sites_.resize(count + ghost_sources_.size());
    for (std::size_t i = 0; i < count; ++i)
        sites_[i] = {packPDVector(position[i]), packPDVector(velocity[i]), volume[i], active[i] == 1};
    for (std::size_t n = 0; n < ghost_sources_.size(); ++n)
    {
        const auto j = ghost_sources_[n];
        sites_[count + n] = {packPDVector(position[j]), packPDVector(velocity[j]), volume[j], false};
    }

    backend_->evaluateDeformation(sites_, bond_live_, dt, tensors_);
    for (std::size_t i = 0; i < count; ++i)
    {
        if (active[i] != 1)
            continue;
        full[i] = unpackNosbPDTensor(tensors_[i].full);
        increment[i] = unpackNosbPDTensor(tensors_[i].increment);
        midpoint[i] = unpackNosbPDTensor(tensors_[i].midpoint);
    }
}

void NosbPDBondWorkspace::evaluateInternalForce(Real reference_density, Vecd *acceleration, Vecd *force)
{
    auto &particles = relation_.getSPHBody().getBaseParticles();
    const auto count = particles.TotalRealParticles();
    refreshLiveBonds(count, particles.ParticlesBound());
    const auto *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const auto *active = particles.getVariableDataByName<int>("ParticleLive");
    const auto *force_state = particles.getVariableDataByName<Matd>("PDForceStateT0");

    force_sites_.resize(count + ghost_sources_.size());
    for (std::size_t i = 0; i < count; ++i)
        force_sites_[i] = {packNosbPDTensor(force_state[i]), volume[i], active[i] == 1};
    for (std::size_t n = 0; n < ghost_sources_.size(); ++n)
    {
        const auto j = ghost_sources_[n];
        force_sites_[count + n] = {packNosbPDTensor(force_state[j]), volume[j], false};
    }

    backend_->assembleInternalForce(force_sites_, bond_live_, force_density_);
    for (std::size_t i = 0; i < count; ++i)
    {
        force[i].setZero();
        if (active[i] != 1)
            continue;
        const auto density = unpackPDVector(force_density_[i]);
        acceleration[i] += density / reference_density;
        force[i] = density * volume[i];
    }
}

void NosbPDBondWorkspace::gatherCorrectionSites(bool damping)
{
    auto &particles = relation_.getSPHBody().getBaseParticles();
    const auto count = particles.TotalRealParticles();
    refreshLiveBonds(count, particles.ParticlesBound());
    const auto *position = particles.getVariableDataByName<Vecd>("Position");
    const auto *velocity = damping ? particles.getVariableDataByName<Vecd>("Velocity") : nullptr;
    const auto *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    const auto *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const auto *damage = damping ? nullptr : particles.getVariableDataByName<Real>("Damage");
    const auto *active = particles.getVariableDataByName<int>("ParticleLive");

    const auto pack = [&](std::size_t i)
    {
        return peridynamics::NosbPDCorrectionSite{
            packPDVector(position[i]), velocity ? packPDVector(velocity[i]) : peridynamics::ParticleVector3{},
            packNosbPDTensor(deformation[i]), volume[i], damage ? damage[i] : Real(0), active[i] == 1};
    };
    correction_sites_.resize(count + ghost_sources_.size());
    for (std::size_t i = 0; i < count; ++i)
        correction_sites_[i] = pack(i);
    for (std::size_t n = 0; n < ghost_sources_.size(); ++n)
        correction_sites_[count + n] = pack(ghost_sources_[n]);
}

void NosbPDBondWorkspace::publishCorrection(bool stabilization, Real damage_cutoff, Real reference_density,
                                            Vecd *acceleration, Vecd *force)
{
    for (std::size_t i = 0; i < correction_force_.size(); ++i)
    {
        force[i].setZero();
        const auto &site = correction_sites_[i];
        if (!site.active || (stabilization && site.damage >= damage_cutoff))
            continue;
        const auto &result = correction_force_[i];
        if (result.collapsed)
            peridynamics::numericalFailure(relation_.getSPHBody().Name(), i, "collapsed stabilization bond");
        const auto density = unpackPDVector(result.density);
        if (stabilization && !density.allFinite())
            peridynamics::numericalFailure(relation_.getSPHBody().Name(), i, "HourglassForce");
        acceleration[i] += density / reference_density;
        force[i] = density * site.volume;
    }
}

void NosbPDBondWorkspace::evaluateStabilization(peridynamics::NosbPDStabilizationParameters parameters,
                                               Real reference_density, Vecd *acceleration, Vecd *force)
{
    gatherCorrectionSites(false);
    backend_->assembleStabilization(correction_sites_, bond_live_, parameters, correction_force_);
    publishCorrection(true, parameters.damage_cutoff, reference_density, acceleration, force);
}

void NosbPDBondWorkspace::evaluateDamping(peridynamics::NosbPDDampingParameters parameters,
                                         Real reference_density, Vecd *acceleration, Vecd *force)
{
    gatherCorrectionSites(true);
    backend_->assembleDamping(correction_sites_, bond_live_, parameters, correction_force_);
    publishCorrection(false, 0, reference_density, acceleration, force);
}

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
            particles_->getVariableDataByName<int>("NumericalFailureReason")[index_i] = 2;
            peridynamics::numericalFailure(sph_body_->Name(), index_i, "non-finite DeformationGradient");
        }
        Real detF = F_[index_i].determinant();
        Real detF_half = F_half_[index_i].determinant();
        /** guard BEFORE inverting: a nearly singular F (volume collapsed to a
         * few percent, beyond the resolution limit) produces an exploding
         * inverse and NaN contamination through the bond network. */
        if (!std::isfinite(detF) || !std::isfinite(detF_half) || detF < 0.05 || detF_half < 0.05)
        {
            particles_->getVariableDataByName<int>("NumericalFailureReason")[index_i] = 3;
            peridynamics::numericalFailure(sph_body_->Name(), index_i, "DeformationGradient determinant below 0.05");
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
      plastic_strain_(particles_->registerStateVariableData<Matd>("PlasticStrain")),
      plastic_work_(particles_->getVariableDataByName<Real>("PlasticWork")) {}
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
            particles_->getVariableDataByName<int>("NumericalFailureReason")[index_i] = 2;
            peridynamics::numericalFailure(sph_body_->Name(), index_i, "non-finite DeformationGradient");
        }
        Real detF = F_[index_i].determinant();
        Real detF_half = F_half_[index_i].determinant();
        /** guard BEFORE inverting: see NosbPDSecondStep::update. */
        if (!std::isfinite(detF) || !std::isfinite(detF_half) || detF < 0.05 || detF_half < 0.05)
        {
            particles_->getVariableDataByName<int>("NumericalFailureReason")[index_i] = 3;
            peridynamics::numericalFailure(sph_body_->Name(), index_i, "DeformationGradient determinant below 0.05");
        }
        F_1_[index_i] = F_[index_i].inverse();
        F_1_half_[index_i] = F_half_[index_i].inverse();
        /** G: Eulerian velocity gradient tensor */
        Matd G = F_delta_[index_i] * F_1_half_[index_i];
        /** return-mapping (convex cutting-plane) with Hughes-Winget objectivity;
         * dt passed so that rate-dependent (visco-plastic) materials can
         * evaluate the overstress term. */
        Matd skew = 0.5 * (G - G.transpose());
        Matd Q = Matd::Identity() + (Matd::Identity() - 0.5 * skew).inverse() * skew;
        Matd rotated_plastic = Q * plastic_strain_[index_i] * Q.transpose();
        stress_[index_i] = plastic_solid_.
            PlasticConstitutiveRelation(G, stress_[index_i], plastic_strain_[index_i], index_i, dt);
        // End-point stress plastic work includes hardening and rate dependence;
        // it is not labelled as heat or as a calibrated fracture energy.
        plastic_work_[index_i] += Vol_[index_i] * detF *
            stress_[index_i].cwiseProduct(plastic_strain_[index_i] - rotated_plastic).sum();
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
    Real val_eq = principalValues(stress_eq).maxCoeff();
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
    Vecd prin_stress = principalValues(stress_eq);
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
    Real sigma1 = principalValues(stress_eq).maxCoeff();
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
    if (stretch_rate > 1.0) return false;
    if (plastic_strain_eq.squaredNorm() <= critical_value_ * critical_value_) return true;
    return principalValues(plastic_strain_eq).maxCoeff() <= critical_value_;
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
      T0_(particles_->getVariableDataByName<Matd>("PDForceStateT0")),
      internal_force_(particles_->getVariableDataByName<Vecd>("PDInternalForce"))
{
    rho0_ = hw_solid_.ReferenceDensity();
}
//=================================================================================================//
void NosbPDThirdStep::interaction(size_t index_i, Real dt)
{
    internal_force_[index_i].setZero();
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
        internal_force_[index_i] = acceleration * Vol_[index_i];
    }
}
//=================================================================================================//
NosbPDFourthStep::NosbPDFourthStep(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      prior_(particles_->getVariableDataByName<Vecd>("PDAccelerationPrior")),
      internal_(particles_->getVariableDataByName<Vecd>("PDInternalForce")),
      hg_(particles_->getVariableDataByName<Vecd>("HourglassForce")),
      damping_(particles_->getVariableDataByName<Vecd>("PDNumericalDampingForce")),
      viscous_(particles_->getVariableDataByName<Vecd>("ViscousForceFromFluid")),
      friction_(particles_->getVariableDataByName<Vecd>("ContactFrictionForce")),
      acoustic_(particles_->getVariableDataByName<Vecd>("AcousticDampingForceFromFluid")),
      internal_work_(particles_->getVariableDataByName<Real>("PDInternalWork")),
      external_work_(particles_->getVariableDataByName<Real>("PDExternalWork")),
      hg_work_(particles_->getVariableDataByName<Real>("PDHourglassWork")),
      damping_work_(particles_->getVariableDataByName<Real>("PDDampingWork")),
      viscous_work_(particles_->getVariableDataByName<Real>("ShellViscousWork")),
      acoustic_work_(particles_->getVariableDataByName<Real>("ShellAcousticWallWork")),
      boundary_power_(particles_->getVariableDataByName<Real>("ContactBoundaryPower")),
      friction_boundary_power_(particles_->getVariableDataByName<Real>("ContactFrictionBoundaryPower")),
      boundary_work_(particles_->getVariableDataByName<Real>("ContactBoundaryWork")),
      friction_loss_(particles_->getVariableDataByName<Real>("ContactFrictionLoss")) {}
//=================================================================================================//
void NosbPDFourthStep::update(size_t index_i, Real dt)
{
    Vecd displacement = dt * (vel_[index_i] + 0.5 * dt * acc_[index_i]);
    internal_work_[index_i] += internal_[index_i].dot(displacement);
    external_work_[index_i] += (mass_[index_i] * acc_[index_i] - internal_[index_i] - hg_[index_i] - damping_[index_i]).dot(displacement);
    hg_work_[index_i] += hg_[index_i].dot(displacement);
    damping_work_[index_i] += damping_[index_i].dot(displacement);
    viscous_work_[index_i] += viscous_[index_i].dot(displacement);
    acoustic_work_[index_i] += acoustic_[index_i].dot(displacement);
    boundary_work_[index_i] += dt * boundary_power_[index_i];
    friction_loss_[index_i] += dt * friction_boundary_power_[index_i] - friction_[index_i].dot(displacement);
    vel_[index_i] += acc_[index_i] * dt;
    checkPositionVerletVelocity(*sph_body_, *particles_, index_i, vel_[index_i]);
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}

void NosbPDFourthStep::evaluateDevice(peridynamics::PositionVerletBatch &backend, Real dt)
{
    const auto count = particles_->TotalRealParticles();
    device_input_.resize(count);
    for (std::size_t i = 0; i < count; ++i)
        device_input_[i] = {packPDVector(pos_[i]), packPDVector(vel_[i]), packPDVector(acc_[i]),
                            packPDVector(internal_[i]), packPDVector(hg_[i]), packPDVector(damping_[i]),
                            packPDVector(viscous_[i]), packPDVector(friction_[i]), packPDVector(acoustic_[i]),
                            mass_[i], boundary_power_[i], friction_boundary_power_[i]};

    backend.finish(device_input_, dt, device_output_);
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto &result = device_output_[i];
        internal_work_[i] += result.internal_work;
        external_work_[i] += result.external_work;
        hg_work_[i] += result.hourglass_work;
        damping_work_[i] += result.damping_work;
        viscous_work_[i] += result.viscous_work;
        acoustic_work_[i] += result.acoustic_work;
        boundary_work_[i] += result.boundary_work;
        friction_loss_[i] += result.friction_loss;
        vel_[i] = unpackPDVector(result.velocity);
        checkPositionVerletVelocity(*sph_body_, *particles_, i, vel_[i]);
        pos_[i] = unpackPDVector(result.position);
    }
}
//=================================================================================================//
ProjectedPairStabilization::
    ProjectedPairStabilization(BaseInnerRelation &inner_relation, Real coefficient, Real damage_cutoff)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, inner_relation.getSPHBody().getMatterMaterial())),
      particleLive_(particles_->getVariableDataByName<int>("ParticleLive")),
      damage_(particles_->getVariableDataByName<Real>("Damage")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      acc_(particles_->registerStateVariableData<Vecd>("PDAcceleration")),
      hourglass_force_(particles_->registerStateVariableData<Vecd>("HourglassForce")),
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value))
{
    if (!std::isfinite(coefficient) || coefficient < 0 || !std::isfinite(damage_cutoff) || damage_cutoff <= 0 || damage_cutoff > 1)
        throw std::invalid_argument("PD stabilization requires a non-negative coefficient and 0 < damage cutoff <= 1");
    particles_->addVariableToWrite<Vecd>("HourglassForce");
    rho0_ = hw_solid_.ReferenceDensity();
    Chg_ = coefficient;
    damage_cutoff_ = damage_cutoff;
    Kbulk_ = hw_solid_.getBulkModulusforHGC(Dimensions);
    horizon_ = sph_body_->getSPHAdaptation().getKernel()->CutOffRadius();
}
//=================================================================================================//
void ProjectedPairStabilization::interaction(size_t index_i, Real dt)
{
    hourglass_force_[index_i].setZero();
    if (particleLive_[index_i] != 1 || damage_[index_i] >= damage_cutoff_ || Chg_ == 0) return;
    // Scale-consistent projected pair penalty. In d dimensions, K/delta^(d+2)
    // times a length error and neighbor measure produces force density (Pa/m).
    Real bond_const = 12.0 * Kbulk_ / (Pi * std::pow(horizon_, Dimensions + 2));
    Vecd force_density = Vecd::Zero();
    const Neighborhood &neighbors = inner_configuration_[index_i];
    for (size_t n = 0; n < neighbors.current_size_; ++n)
    {
        size_t j = neighbors.j_[n];
        // Eligibility and endpoint averaging are symmetric; only this row's
        // acceleration is written, so parallel execution does not scatter data.
        if (!neighbors.bondLive_[n] || particleLive_[j] != 1 || damage_[j] >= damage_cutoff_) continue;
        Vecd reference = -neighbors.r_ij_[n] * neighbors.e_ij_[n];
        Vecd current = pos_[j] - pos_[index_i];
        Real length = current.norm();
        if (!std::isfinite(length) || length <= SqrtEps * neighbors.r_ij_[n])
            peridynamics::numericalFailure(sph_body_->Name(), index_i, "collapsed stabilization bond");
        Vecd error = (0.5 * (F_[index_i] + F_[j])) * reference - current;
        Real projection = error.dot(current);
        force_density -= Chg_ * bond_const * projection / (neighbors.r_ij_[n] * length) * Vol_[j] * current;
    }
    if (!force_density.allFinite())
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "HourglassForce");
    acc_[index_i] += force_density / rho0_;
    hourglass_force_[index_i] = force_density * Vol_[index_i];
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
      F_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      damping_force_(particles_->getVariableDataByName<Vecd>("PDNumericalDampingForce"))
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
    damping_force_[index_i].setZero();
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
        damping_force_[index_i] = damping_force * Vol_[index_i];
    }
}
//=================================================================================================//
void NosbPDStabilizationWithExecution::runMainStep(Real dt)
{
    if (!workspace_.onDevice() || Chg_ == 0)
    {
        Base::runMainStep(dt);
        return;
    }
    const Real bond_constant = 12.0 * Kbulk_ / (Pi * std::pow(horizon_, Dimensions + 2));
    workspace_.evaluateStabilization({Chg_ * bond_constant, damage_cutoff_, SqrtEps}, rho0_, acc_, hourglass_force_);
}

NosbPDDampingWithExecution::NosbPDDampingWithExecution(NosbPDBondWorkspace &workspace)
    : Base(workspace.relation()), workspace_(workspace)
{
    // These exact types use the linear pair law of HughesWingetSolid.
    // User-supplied virtual pair damping laws continue through the CPU operator.
    const auto &material = typeid(hw_solid_);
    supported_ = material == typeid(HughesWingetSolid) || material == typeid(J2PlasticityforPD) ||
        material == typeid(ViscoPlasticSolidforPD) || material == typeid(DruckerPragerPlasticityforPD);
}

void NosbPDDampingWithExecution::runMainStep(Real dt)
{
    if (onDevice())
        workspace_.evaluateDamping({rho0_, c0_, smoothing_length_, numerical_dissipation_factor_, inv_W0_, Dimensions},
                                   rho0_, acc_, damping_force_);
    else
        Base::runMainStep(dt);
}
} // namespace solid_dynamics
} // namespace SPH

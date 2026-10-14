/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_contact.cpp
 * @brief 	Solid-wall contact for peridynamic bodies.
 */

#include "pd_contact.h"

#include "adaptation.h"
#include "base_particles.hpp"

#include <cmath>

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
ContactRepulsionFactorforPD::
    ContactRepulsionFactorforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      repulsion_factor_(particles_->registerStateVariableData<Real>("RepulsionFactor"))
{
    Kernel *source_kernel = sph_body_->getSPHAdaptation().getKernel();
    Real dp = sph_body_->getSPHAdaptation().ReferenceSpacing();
    for (auto *body : contact_bodies_)
    {
        Kernel *target_kernel = body->getSPHAdaptation().getKernel();
        Kernel *kernel = target_kernel->SmoothingLength() > source_kernel->SmoothingLength() ? target_kernel : source_kernel;
        offset_W_ij_.push_back(kernel->W(dp, ZeroVecd));
    }
    particles_->addVariableToWrite<Real>("RepulsionFactor");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->registerStateVariableData<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void ContactRepulsionFactorforPD::interaction(size_t index_i, Real dt)
{
    Real sigma = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], Real(0));
            sigma += corrected_W_ij * Vol_k[contact_neighborhood.j_[n]];
        }
    }
    repulsion_factor_[index_i] = sigma;
}
//=================================================================================================//
ContactRepulsionFromWallforPD::
    ContactRepulsionFromWallforPD(BaseContactRelation &contact_relation, Real friction_mu, Real friction_veps)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, contact_relation.getSPHBody().getMatterMaterial())),
      contact_stiffness_(0.0),
      friction_mu_(friction_mu), friction_veps_(friction_veps),
      vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      repulsion_factor_(particles_->getVariableDataByName<Real>("RepulsionFactor")),
      acc_prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")),
      repulsion_force_(particles_->registerStateVariableData<Vecd>("ContactRepulsionFromWall"))
{
    if (!std::isfinite(friction_mu_) || friction_mu_ < 0.0 || !std::isfinite(friction_veps_) || friction_veps_ <= 0.0)
        throw std::invalid_argument("PD contact: finite non-negative friction and positive regularization speed required");
    for (size_t k = 0; k < contact_particles_.size(); ++k)
    {
        contact_velocity_.push_back(DynamicCast<Solid>(this, contact_bodies_[k]->getMatterMaterial()).AverageVelocity(contact_particles_[k]));
        contact_normal_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("NormalDirection"));
    }
    contact_stiffness_ = hw_solid_.ContactStiffness();
    particles_->addVariableToWrite<Vecd>("ContactRepulsionFromWall");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->registerStateVariableData<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void ContactRepulsionFromWallforPD::interaction(size_t index_i, Real dt)
{
    Real pressure = repulsion_factor_[index_i] * contact_stiffness_;
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        const Neighborhood &neighbors = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != neighbors.current_size_; ++n)
        {
            size_t j = neighbors.j_[n];
            Vecd normal_force = -2.0 * pressure * Vol_[index_i] * neighbors.e_ij_[n] * neighbors.dW_ij_[n] * contact_Vol_[k][j];
            force += pdContactPairForce(normal_force, vel_[index_i] - contact_velocity_[k][j], contact_normal_[k][j], friction_mu_, friction_veps_);
        }
    }
    repulsion_force_[index_i] = force;
    acc_prior_[index_i] += force / (hw_solid_.ReferenceDensity() * Vol_[index_i]);
}
//=================================================================================================//
ContactRepulsionToWallforPD::
    ContactRepulsionToWallforPD(BaseContactRelation &contact_relation, Real friction_mu, Real friction_veps)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      Vol_(particles_->registerStateVariableData<Real>("VolumetricMeasure")),
      repulsion_force_(particles_->registerStateVariableData<Vecd>("ContactRepulsionToWall")),
      velocity_(DynamicCast<Solid>(this, contact_relation.getSPHBody().getMatterMaterial()).AverageVelocity(particles_)),
      normal_(particles_->registerStateVariableData<Vecd>("NormalDirection")),
      friction_mu_(friction_mu), friction_veps_(friction_veps)
{
    particles_->addVariableToWrite<Vecd>("ContactRepulsionToWall");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_velocity_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
        contact_solids_.push_back(&DynamicCast<Solid>(this, contact_bodies_[k]->getMatterMaterial()));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_repulsion_factor_.push_back(contact_particles_[k]->getVariableDataByName<Real>("RepulsionFactor"));
    }
}
//=================================================================================================//
void ContactRepulsionToWallforPD::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Solid *solid_k = contact_solids_[k];
        Real *repulsion_factor_k = contact_repulsion_factor_[k];

        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            /** contact pressure from the PD-side factor, following RepulsionForce<Wall, Contact<>> */
            Real p_star = repulsion_factor_k[index_j] * solid_k->ContactStiffness();
            Vecd shell_normal_force = 2.0 * p_star * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * Vol_[index_i];
            force -= pdContactPairForce(shell_normal_force, contact_velocity_[k][index_j] - velocity_[index_i], normal_[index_i], friction_mu_, friction_veps_);
        }
    }
    repulsion_force_[index_i] = force;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH

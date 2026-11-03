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
      repulsion_force_(particles_->registerStateVariableData<Vecd>("ContactRepulsionFromWall")),
      friction_force_(particles_->getVariableDataByName<Vecd>("ContactFrictionForce")),
      boundary_power_(particles_->getVariableDataByName<Real>("ContactBoundaryPower")),
      friction_boundary_power_(particles_->getVariableDataByName<Real>("ContactFrictionBoundaryPower"))
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
    friction_force_[index_i].setZero();
    boundary_power_[index_i] = friction_boundary_power_[index_i] = 0;
    Real pressure = repulsion_factor_[index_i] * contact_stiffness_;
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        const Neighborhood &neighbors = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != neighbors.current_size_; ++n)
        {
            size_t j = neighbors.j_[n];
            Vecd normal_force = -2.0 * pressure * Vol_[index_i] * neighbors.e_ij_[n] * neighbors.dW_ij_[n] * contact_Vol_[k][j];
            Vecd pair = pdContactPairForce(normal_force, vel_[index_i] - contact_velocity_[k][j], contact_normal_[k][j], friction_mu_, friction_veps_);
            force += pair;
            friction_force_[index_i] += pair - normal_force;
            boundary_power_[index_i] += pair.dot(contact_velocity_[k][j]);
            friction_boundary_power_[index_i] += (pair - normal_force).dot(contact_velocity_[k][j]);
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

RigidBoxContactforPD::RigidBoxContactforPD(SPHBody &body, Real friction_mu, Real friction_veps)
    : LocalDynamics(body), spacing_(body.getSPHAdaptation().ReferenceSpacing()),
      stiffness_(DynamicCast<HughesWingetSolid>(this, body.getMatterMaterial()).ContactStiffness()),
      rho0_(DynamicCast<HughesWingetSolid>(this, body.getMatterMaterial()).ReferenceDensity()),
      mu_(friction_mu), veps_(friction_veps),
      position_(particles_->getVariableDataByName<Vecd>("Position")),
      velocity_(particles_->registerStateVariableData<Vecd>("Velocity")),
      prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")),
      force_(particles_->registerStateVariableData<Vecd>("ContactRepulsionFromWall")),
      friction_(particles_->getVariableDataByName<Vecd>("ContactFrictionForce")),
      volume_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      energy_(particles_->registerStateVariableData<Real>("AnalyticContactEnergy")),
      boundary_power_(particles_->getVariableDataByName<Real>("ContactBoundaryPower")),
      friction_boundary_power_(particles_->getVariableDataByName<Real>("ContactFrictionBoundaryPower"))
{
    acceleration_ = particles_->registerStateVariableData<Vecd>("PDAcceleration");
    if (!std::isfinite(mu_) || mu_ < 0 || !std::isfinite(veps_) || veps_ <= 0)
        throw std::invalid_argument("RigidBoxContactforPD: invalid friction parameters");
}

size_t RigidBoxContactforPD::addBox(const std::string &name, const Vecd &center, const Vecd &half_size, const Vecd &velocity)
{
    if (!center.allFinite() || !half_size.allFinite() || !velocity.allFinite() || half_size.minCoeff() <= 0)
        throw std::invalid_argument("RigidBoxContactforPD: invalid box geometry or motion");
    auto *reaction = particles_->registerStateVariableData<Vecd>(name + "Reaction");
    boxes_.push_back({center, half_size, velocity, reaction});
    return boxes_.size()-1;
}

void RigidBoxContactforPD::setKinematics(size_t box, const Vecd &center, const Vecd &velocity)
{
    if (!center.allFinite() || !velocity.allFinite()) throw std::invalid_argument("non-finite rigid-box motion");
    boxes_.at(box).center = center; boxes_.at(box).velocity = velocity;
}

RigidBoxContactResponse RigidBoxContactforPD::evaluate(size_t particle, size_t box_index) const
{
    const auto &box = boxes_.at(box_index);
    const auto geometry = rigidBoxDistance(position_[particle], box.center, box.half_size);
    RigidBoxContactResponse response;
    response.penetration = std::max(0.5 * spacing_ - geometry.first, Real(0));
    if (response.penetration == 0)
        return response;

    const Real stiffness = stiffness_ * volume_[particle] / (spacing_ * spacing_);
    const Vecd normal_force = stiffness * response.penetration * geometry.second;
    response.force = pdContactPairForce(normal_force, velocity_[particle] - box.velocity,
                                       geometry.second, mu_, veps_);
    response.friction = response.force - normal_force;
    response.potential_energy = 0.5 * stiffness * response.penetration * response.penetration;
    response.boundary_power = response.force.dot(box.velocity);
    response.friction_boundary_power = response.friction.dot(box.velocity);
    return response;
}

StdVec<RigidBoxContactObservation> RigidBoxContactforPD::observe() const
{
    StdVec<RigidBoxContactObservation> observations(boxes_.size());
    for (size_t i = 0; i < particles_->TotalRealParticles(); ++i)
        for (size_t k = 0; k < boxes_.size(); ++k)
        {
            const auto response = evaluate(i, k);
            auto &total = observations[k];
            total.reaction -= response.force;
            total.potential_energy += response.potential_energy;
            total.boundary_power += response.boundary_power;
            total.max_penetration = std::max(total.max_penetration, response.penetration);
            total.contacting_particles += response.penetration > 0;
        }
    return observations;
}

void RigidBoxContactforPD::update(size_t i, Real)
{
    force_[i].setZero(); friction_[i].setZero();
    energy_[i] = boundary_power_[i] = friction_boundary_power_[i] = 0;
    for (size_t k = 0; k < boxes_.size(); ++k)
    {
        const auto response = evaluate(i, k);
        force_[i] += response.force;
        friction_[i] += response.friction;
        boxes_[k].reaction[i] = -response.force;
        energy_[i] += response.potential_energy;
        boundary_power_[i] += response.boundary_power;
        friction_boundary_power_[i] += response.friction_boundary_power;
    }
    if (apply_acceleration_)
        (current_acceleration_ ? acceleration_[i] : prior_[i]) += force_[i] / (rho0_ * volume_[i]);
}
} // namespace solid_dynamics
} // namespace SPH

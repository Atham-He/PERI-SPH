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
    /** The contact neighborhoods are built with the largest smoothing length
     * among source and contact bodies (NeighborBuilder::chooseKernel),
     * so the W offset must be evaluated with the same kernel. */
    Kernel *contact_kernel = contact_relation.getSPHBody().getSPHAdaptation().getKernel();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Kernel *kernel_k = contact_bodies_[k]->getSPHAdaptation().getKernel();
        if (kernel_k->SmoothingLength() > contact_kernel->SmoothingLength())
            contact_kernel = kernel_k;
    }
    Real dp = contact_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing();
    offset_W_ij_ = contact_kernel->W(dp, ZeroVecd);
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
            Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_, Real(0));
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
    /** note: fully debonded (free) particles are NOT skipped — shell
     * fragments must still be pushed out of the wall as ordinary debris. */

    /** contact pressure scaled by the kernel-summation factor,
     * following RepulsionForce<Contact<Wall>> */
    Real p_i = repulsion_factor_[index_i] * contact_stiffness_;
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            force -= 2.0 * p_i * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    repulsion_force_[index_i] = force * Vol_[index_i];
    /** smooth-Coulomb tangential friction: f_t = -mu |f_n| tanh(|v_t|/v_eps) v_t_hat,
     * with the tangential relative velocity measured against the wall (zero tangential motion).
     * Calibration need: without friction a soft intact shell squirts out sideways
     * between the plates at ~40% stroke (watermelon-seed effect). */
    if (friction_mu_ > 0.0)
    {
        Vecd f_n_vec = repulsion_force_[index_i];
        Real f_n = f_n_vec.norm();
        if (f_n > 0.0)
        {
            Vecd n_hat = f_n_vec / f_n;
            Vecd v_i = vel_[index_i];
            Vecd v_t = v_i - n_hat * v_i.dot(n_hat);
            Real vt_mag = v_t.norm();
            if (vt_mag > 0.0)
            {
                Vecd f_t = -friction_mu_ * f_n * std::tanh(vt_mag / friction_veps_) * (v_t / vt_mag);
                repulsion_force_[index_i] += f_t;
                acc_prior_[index_i] += f_t / (hw_solid_.ReferenceDensity() * Vol_[index_i]);
            }
        }
    }
    /** additive contribution onto the baseline set by PDTimeStepInitialization */
    acc_prior_[index_i] += repulsion_force_[index_i] / (hw_solid_.ReferenceDensity() * Vol_[index_i]);
}
//=================================================================================================//
ContactRepulsionToWallforPD::
    ContactRepulsionToWallforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      Vol_(particles_->registerStateVariableData<Real>("VolumetricMeasure")),
      repulsion_force_(particles_->registerStateVariableData<Vecd>("ContactRepulsionToWall"))
{
    particles_->addVariableToWrite<Vecd>("ContactRepulsionToWall");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
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
            force -= 2.0 * p_star * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    repulsion_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH

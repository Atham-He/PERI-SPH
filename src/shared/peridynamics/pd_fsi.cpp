/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_fsi.cpp
 * @brief 	Fluid-structure interaction for peridynamic bodies.
 */

#include "pd_fsi.h"

#include "base_particles.hpp"
#include "cell_linked_list.hpp"

namespace SPH
{
//=================================================================================================//
ContactRelationSelectKernel::
    ContactRelationSelectKernel(SPHBody &sph_body, RealBodyVector contact_bodies,
                                bool use_source_kernel)
    : ContactRelationCrossResolution(sph_body, contact_bodies)
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_contact_neighbors_.push_back(
            neighbor_builder_contacts_keeper_.createPtr<NeighborBuilderContactSelectKernel>(
                sph_body_, *contact_bodies_[k], use_source_kernel));
    }
}
//=================================================================================================//
void ContactRelationSelectKernel::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
VelocityDivergenceInner::
    VelocityDivergenceInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      u_div_(particles_->registerStateVariableData<Real>("VelocityDivergence")) {}
//=================================================================================================//
void VelocityDivergenceInner::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        density_change_rate += (vel_[index_i] - vel_[index_j]).dot(inner_neighborhood.e_ij_[n]) *
                               inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }
    u_div_[index_i] = -density_change_rate;
}
//=================================================================================================//
VelocityDivergenceWithWallforPD::
    VelocityDivergenceWithWallforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      u_div_(particles_->registerStateVariableData<Real>("VelocityDivergence")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_ave_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("AverageVelocity"));
    }
}
//=================================================================================================//
void VelocityDivergenceWithWallforPD::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd vel_in_wall = 2.0 * contact_vel_ave_[k][index_j] - vel_[index_i];
            density_change_rate += (vel_[index_i] - vel_in_wall).dot(contact_neighborhood.e_ij_[n]) *
                                   contact_neighborhood.dW_ij_[n] * contact_Vol_[k][index_j];
        }
    }
    u_div_[index_i] -= density_change_rate;
}
//=================================================================================================//
AcousticDamperInnerforPD::
    AcousticDamperInnerforPD(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      u_div_(particles_->getVariableDataByName<Real>("VelocityDivergence")),
      force_(particles_->registerStateVariableData<Vecd>("Force")) {}
//=================================================================================================//
void AcousticDamperInnerforPD::interaction(size_t index_i, Real dt)
{
    Vecd acoustic_damper = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        acoustic_damper += (u_div_[index_i] + u_div_[index_j]) *
                           inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    }
    force_[index_i] += coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
}
//=================================================================================================//
AcousticDamperWithWallforPD::
    AcousticDamperWithWallforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      u_div_(particles_->getVariableDataByName<Real>("VelocityDivergence")),
      force_(particles_->registerStateVariableData<Vecd>("Force"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void AcousticDamperWithWallforPD::interaction(size_t index_i, Real dt)
{
    Vecd acoustic_damper = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            acoustic_damper += (2.0 * u_div_[index_i]) *
                               contact_neighborhood.dW_ij_[n] * contact_Vol_[k][index_j] *
                               contact_neighborhood.e_ij_[n];
        }
    }
    force_[index_i] += coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics

namespace solid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
PressureForceFromFluidforPD<RiemannSolverType>::
    PressureForceFromFluidforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, contact_relation.getSPHBody().getMatterMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      vel_ave_(particles_->registerStateVariableData<Vecd>("AverageVelocity")),
      acc_prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")),
      acc_ave_(particles_->registerStateVariableData<Vecd>("AverageAcceleration")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      force_from_fluid_(particles_->registerStateVariableData<Vecd>("PressureForceFromFluid")),
      pressure_from_fluid_(particles_->registerStateVariableData<Real>("PressureFromFluid"))
{
    particles_->addVariableToWrite<Vecd>("PressureForceFromFluid");
    particles_->addVariableToWrite<Real>("PressureFromFluid");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_fluids_.push_back(&DynamicCast<Fluid>(this, contact_bodies_[k]->getMatterMaterial()));
        contact_rho_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Density"));
        contact_p_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Pressure"));
        contact_u_div_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VelocityDivergence"));
        contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Velocity"));
        contact_force_prior_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("ForcePrior"));
        riemann_solvers_.push_back(RiemannSolverType(*contact_fluids_[k], *contact_fluids_[k]));
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void PressureForceFromFluidforPD<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Vecd acoustic_damper = Vecd::Zero();
    Real p_from_f = 0.0;
    Real sum_w = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_Vol_[k][index_j];

            Real face_wall_external_acceleration =
                (contact_force_prior_[k][index_j] / contact_mass_[k][index_j] - acc_ave_[index_i]).dot(e_ij);
            Real p_in_wall = contact_p_[k][index_j] +
                             contact_rho_[k][index_j] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            p_from_f += contact_p_[k][index_j] * contact_neighborhood.W_ij_[n];
            sum_w += contact_neighborhood.W_ij_[n];

            acoustic_damper += (contact_u_div_[k][index_j] * 2.0) * Vol_[index_i] * dW_ijV_j * e_ij;
            /** Interface dissipation consistent with the SPHinXsys standard FSI/wall
             * convention (fluid_structure_interaction.hpp): face_to_fluid_n points
             * from the fluid neighbor toward the shell, so the sign of u_jump and the
             * force term match the fluid-side Riemann solver exactly. */
            Vecd face_to_fluid_n = -SGN(e_ij.dot(n_[index_i])) * n_[index_i];
            Real u_jump = 2.0 * (contact_vel_[k][index_j] - vel_ave_[index_i]).dot(face_to_fluid_n);

            force -= (riemann_solvers_[k].DissipativePJump(u_jump) * face_to_fluid_n +
                      (p_in_wall + contact_p_[k][index_j]) * e_ij) *
                     Vol_[index_i] * dW_ijV_j;
        }
    }
    force_from_fluid_[index_i] = force + coeff_acoustic_damper_ * acoustic_damper;
    /** isolate a non-finite fluid force (e.g. from a blown-up fluid
     * neighborhood) so it cannot poison the PD acceleration. */
    if (force_from_fluid_[index_i].allFinite())
    {
        /** additive contribution: the baseline of PDAccelerationPrior (e.g. gravity)
         * must be set by PDTimeStepInitialization executed before this dynamics. */
        acc_prior_[index_i] += force_from_fluid_[index_i] / (hw_solid_.ReferenceDensity() * Vol_[index_i]);
    }
    else
    {
        force_from_fluid_[index_i].setZero();
    }
    sum_w = (sum_w == 0.0 ? 1.0 : sum_w);
    pressure_from_fluid_[index_i] = p_from_f / sum_w;
}
//=================================================================================================//
template class PressureForceFromFluidforPD<NoRiemannSolver>;
template class PressureForceFromFluidforPD<AcousticRiemannSolver>;
} // namespace solid_dynamics
} // namespace SPH

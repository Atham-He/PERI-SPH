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
#include "pd_numerics.h"

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
PressureWallTermsforPD::PressureWallTermsforPD(BaseContactRelation &relation, Real factor)
    : Integration1stHalf<Contact<Wall>, AcousticRiemannSolver, NoKernelCorrection>(relation)
{
    if (!std::isfinite(factor) || factor < 0 || factor > 1)
        throw std::invalid_argument("wall acceleration pressure factor must be in [0,1]");
    for (auto *body : contact_bodies_)
        acceleration_factors_.push_back(dynamic_cast<HughesWingetSolid *>(&body->getMatterMaterial()) ? factor : 1);
}
void PressureWallTermsforPD::interaction(size_t i, Real)
{
    Vecd force=Vecd::Zero(); Real rho_dissipation=0;
    for(size_t k=0;k<contact_configuration_.size();++k)
    {
        const Neighborhood &neighbors=(*contact_configuration_[k])[i];
        for(size_t n=0;n<neighbors.current_size_;++n)
        {
            size_t j=neighbors.j_[n]; Vecd direction=neighbors.e_ij_[n];
            Real gradient=neighbors.dW_ij_[n]*wall_Vol_[k][j];
            Real acceleration=(force_prior_[i]/mass_[i]-wall_acc_ave_[k][j]).dot(-direction);
            Real wall_pressure=p_[i]+acceleration_factors_[k]*rho_[i]*neighbors.r_ij_[n]*std::max(acceleration,Real(0));
            force-=(p_[i]+wall_pressure)*gradient*direction;
            rho_dissipation+=riemann_solver_.DissipativeUJump(p_[i]-wall_pressure)*gradient;
        }
    }
    force_[i]+=force*Vol_[i]; drho_dt_[i]+=rho_dissipation*rho_[i];
}

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
        contact_vel_ave_.push_back(DynamicCast<Solid>(this, contact_bodies_[k]->getMatterMaterial()).AverageVelocity(contact_particles_[k]));
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
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      work_(particles_->registerStateVariableData<Real>("AcousticDampingWorkInner")),
      force_(particles_->registerStateVariableData<Vecd>("AcousticDampingForceInner")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    particles_->addEvolvingVariable<Real>("AcousticDampingWorkInner");
}
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
    force_[index_i] = mass_[index_i] * coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
    if (!force_[index_i].allFinite() || !std::isfinite(dt) || dt < 0.0)
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "AcousticDampingForce");
    work_[index_i] += dt * force_[index_i].dot(vel_[index_i] + 0.5 * dt * force_[index_i] / mass_[index_i]);
    vel_[index_i] += dt * force_[index_i] / mass_[index_i];
}
//=================================================================================================//
AcousticDamperWithWallforPD::
    AcousticDamperWithWallforPD(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      u_div_(particles_->getVariableDataByName<Real>("VelocityDivergence")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      work_(particles_->registerStateVariableData<Real>("AcousticDampingWorkWall")),
      force_(particles_->registerStateVariableData<Vecd>("AcousticDampingForceWall")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    particles_->addEvolvingVariable<Real>("AcousticDampingWorkWall");
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
    force_[index_i] = mass_[index_i] * coeff_acoustic_damper_ * acoustic_damper / rho_[index_i];
    if (!force_[index_i].allFinite() || !std::isfinite(dt) || dt < 0.0)
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "AcousticDampingForce");
    work_[index_i] += dt * force_[index_i].dot(vel_[index_i] + 0.5 * dt * force_[index_i] / mass_[index_i]);
    vel_[index_i] += dt * force_[index_i] / mass_[index_i];
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
      acoustic_force_from_fluid_(particles_->getVariableDataByName<Vecd>("AcousticDampingForceFromFluid")),
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
        Real impedance=0., inverse_sound=0.;
        if constexpr (std::is_same_v<RiemannSolverType,AcousticRiemannSolver>)
        {
            // This relation supplies the same fluid material on each side
            // of the Riemann interface.
            const Real rho0=contact_fluids_[k]->ReferenceDensity();
            const Real rho_c0=rho0*contact_fluids_[k]->ReferenceSoundSpeed();
            const Real impedance_sum=rho_c0+rho_c0;
            impedance=2.*rho_c0*rho_c0*(1./impedance_sum);
            inverse_sound=.5*(rho0+rho0)*(impedance_sum/(rho_c0*rho_c0+rho_c0*rho_c0));
        }
        acoustic_coefficients_.emplace_back(impedance,inverse_sound);
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
                             acceleration_factor_ * contact_rho_[k][index_j] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
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
    publishTraction(index_i,force,acoustic_damper,p_from_f,sum_w);
}

namespace
{
peridynamics::ParticleVector3 packTractionVector(const Vecd &input)
{
    peridynamics::ParticleVector3 x{};
    for (int k=0;k<Dimensions;++k)
        x.component[k]=input[k];
    return x;
}

Vecd unpackTractionVector(const peridynamics::ParticleVector3 &input)
{
    Vecd x;
    for (int k=0;k<Dimensions;++k)
        x[k]=input.component[k];
    return x;
}
}

template<class RiemannSolverType>
void PressureForceFromFluidforPD<RiemannSolverType>::evaluateDevice(peridynamics::FluidTractionBatch &backend)
{
    using namespace peridynamics;
    const size_t count=particles_->TotalRealParticles();
    solid_rows_.resize(count);
    row_offsets_.assign(count+1,0);
    traction_neighbors_.clear();
    fluid_sites_.clear();
    fluid_indices_.resize(contact_particles_.size());
    constexpr size_t missing=std::numeric_limits<size_t>::max();
    for (size_t k=0;k<contact_particles_.size();++k)
        fluid_indices_[k].assign(contact_particles_[k]->ParticlesBound(),missing);

    for (size_t i=0;i<count;++i)
    {
        solid_rows_[i]={packTractionVector(n_[i]),packTractionVector(vel_ave_[i]),
                        packTractionVector(acc_ave_[i]),Vol_[i]};
        for (size_t k=0;k<contact_configuration_.size();++k)
        {
            const auto &neighbors=(*contact_configuration_[k])[i];
            for (size_t n=0;n<neighbors.current_size_;++n)
            {
                const size_t j=neighbors.j_[n];
                auto &fluid_id=fluid_indices_[k].at(j);
                if (fluid_id==missing)
                {
                    fluid_id=fluid_sites_.size();
                    const auto [impedance,inverse_sound]=acoustic_coefficients_[k];
                    fluid_sites_.push_back({packTractionVector(contact_vel_[k][j]),
                        packTractionVector(contact_force_prior_[k][j]),contact_mass_[k][j],contact_rho_[k][j],
                        contact_p_[k][j],contact_u_div_[k][j],contact_Vol_[k][j],impedance,inverse_sound});
                }
                traction_neighbors_.push_back({fluid_id,packTractionVector(neighbors.e_ij_[n]),
                    neighbors.r_ij_[n],neighbors.W_ij_[n],neighbors.dW_ij_[n]});
            }
        }
        row_offsets_[i+1]=traction_neighbors_.size();
    }
    backend.evaluate(solid_rows_,fluid_sites_,row_offsets_,traction_neighbors_,acceleration_factor_,traction_result_);
    if (traction_result_.size()!=count)
        throw std::logic_error("PD fluid traction backend returned an incomplete force field");
    for (size_t i=0;i<count;++i)
    {
        const auto &result=traction_result_[i];
        publishTraction(i,unpackTractionVector(result.pressure_force),unpackTractionVector(result.acoustic_force),
                        result.pressure,result.kernel_weight);
    }
}

template<class RiemannSolverType>
void PressureForceFromFluidforPD<RiemannSolverType>::publishTraction(
    size_t index_i,const Vecd &force,const Vecd &acoustic_damper,Real p_from_f,Real sum_w)
{
    acoustic_force_from_fluid_[index_i] = coeff_acoustic_damper_ * acoustic_damper;
    force_from_fluid_[index_i] = force + acoustic_force_from_fluid_[index_i];
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
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "PressureForceFromFluid");
    }
    sum_w = (sum_w == 0.0 ? 1.0 : sum_w);
    pressure_from_fluid_[index_i] = p_from_f / sum_w;
}
//=================================================================================================//
template class PressureForceFromFluidforPD<NoRiemannSolver>;
template class PressureForceFromFluidforPD<AcousticRiemannSolver>;
} // namespace solid_dynamics
} // namespace SPH

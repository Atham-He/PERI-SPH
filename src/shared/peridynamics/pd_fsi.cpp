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
namespace
{
peridynamics::ParticleVector3 packFluidVector(const Vecd &input)
{
    peridynamics::ParticleVector3 x{};
    for (int k = 0; k < Dimensions; ++k)
        x.component[k] = input[k];
    return x;
}

Vecd unpackFluidVector(const peridynamics::ParticleVector3 &input)
{
    Vecd x;
    for (int k = 0; k < Dimensions; ++k)
        x[k] = input.component[k];
    return x;
}
Real pdWallAccelerationFactor(SPHBody &body, Real factor)
{
    return dynamic_cast<HughesWingetSolid *>(&body.getMatterMaterial()) ? factor : 1;
}

}

namespace peridynamics
{
FluidNeighborConfiguration::Scope::Scope(FluidNeighborConfiguration &configuration)
    : configuration_(configuration)
{
    if (configuration_.depth_++ == 0)
        ++configuration_.generation_;
}

FluidNeighborConfiguration::Scope::~Scope()
{
    --configuration_.depth_;
}

void SPHViscosityWorkspace::addTarget(Form form, const ParticleConfiguration &configuration,
                                     const Vecd *velocity, const Real *volume, std::size_t bound,
                                     Real viscosity, Real smoothing_length)
{
    targets_.push_back({form, &configuration, velocity, volume, bound, viscosity, smoothing_length});
}

void SPHViscosityWorkspace::configure(FluidViscosityBatch &backend, std::size_t count)
{
    cache_.invalidate();
    constexpr auto missing = std::numeric_limits<std::size_t>::max();
    target_sources_.clear();
    rows_.offsets.assign(count + 1, 0);
    rows_.inner_ends.resize(count);
    rows_.pairs.clear();
    site_ids_.resize(targets_.size());

    for (std::size_t k = 0; k < targets_.size(); ++k)
    {
        const auto &target = targets_[k];
        const auto shared = target.form == Form::Inner ? count : 0;
        site_ids_[k].assign(target.bound - shared, missing);
    }

    for (std::size_t i = 0; i < count; ++i)
    {
        rows_.inner_ends[i] = rows_.pairs.size();
        for (std::size_t k = 0; k < targets_.size(); ++k)
        {
            const auto &target = targets_[k];
            const auto &neighbors = (*target.configuration)[i];
            const bool inner = target.form == Form::Inner;
            const auto shared = inner ? count : 0;
            const Real velocity_factor = inner ? 1.0 : 2.0;
            const Real offset = 0.01 * target.smoothing_length;
            for (std::size_t n = 0; n < neighbors.current_size_; ++n)
            {
                const auto j = neighbors.j_[n];
                auto site = j;
                if (j >= shared)
                {
                    auto &id = site_ids_[k].at(j - shared);
                    if (id == missing)
                    {
                        id = count + target_sources_.size();
                        target_sources_.emplace_back(k, j);
                    }
                    site = id;
                }

                const auto &e = neighbors.e_ij_[n];
                Real force_factor = 2.0 * target.viscosity;
                if (inner)
                    force_factor = e.dot(2.0 * e) * target.viscosity;
                else if (target.form == Form::Wall)
                    force_factor = 2.0 * e.dot(e) * target.viscosity;

                rows_.pairs.push_back({site, neighbors.r_ij_[n] + offset, velocity_factor,
                                       force_factor, neighbors.dW_ij_[n]});
            }
            if (inner)
                rows_.inner_ends[i] = rows_.pairs.size();
        }
        rows_.offsets[i + 1] = rows_.pairs.size();
    }

    if (!stencil_ || configured_backend_ != &backend)
        stencil_ = backend.makeStencil();
    stencil_->configure(rows_, count + target_sources_.size());
    configured_backend_ = &backend;
    cache_.configured(count);
}

void SPHViscosityWorkspace::evaluate(FluidViscosityBatch &backend, std::size_t count,
                                    const Vecd *velocity, const Real *volume, Vecd *force)
{
    if (cache_.needsConfiguration(count) || configured_backend_ != &backend)
        configure(backend, count);

    sites_.resize(count + target_sources_.size());
    for (std::size_t i = 0; i < count; ++i)
        sites_[i] = {packFluidVector(velocity[i]), volume[i]};
    for (std::size_t n = 0; n < target_sources_.size(); ++n)
    {
        const auto [body, particle] = target_sources_[n];
        const auto &target = targets_[body];
        sites_[count + n] = {packFluidVector(target.velocity[particle]), target.volume[particle]};
    }

    backend.evaluate(*stencil_, sites_, forces_);
    for (std::size_t i = 0; i < count; ++i)
        force[i] = unpackFluidVector(forces_[i]);
}
}

namespace fluid_dynamics
{
void ViscousWallTermsforPD::addDeviceTargets(peridynamics::SPHViscosityWorkspace &workspace)
{
    for (std::size_t k = 0; k < contact_configuration_.size(); ++k)
        workspace.addTarget(peridynamics::SPHViscosityWorkspace::Form::Wall, *contact_configuration_[k],
                            wall_vel_ave_[k], wall_Vol_[k], contact_particles_[k]->ParticlesBound(),
                            mu_(0, 0), smoothing_length_);
}

void ViscousForceWithWallforPD::evaluateDevice(peridynamics::FluidViscosityBatch &backend,
                                              peridynamics::SPHViscosityWorkspace &workspace)
{
    workspace.clearTargets();
    workspace.addTarget(peridynamics::SPHViscosityWorkspace::Form::Inner, inner_configuration_,
                        vel_, Vol_, particles_->ParticlesBound(), mu_(0, 0), smoothing_length_);
    walls_.addDeviceTargets(workspace);
    workspace.evaluate(backend, particles_->TotalRealParticles(), vel_, Vol_, viscous_force_);
}

PressureWallTermsforPD::PressureWallTermsforPD(BaseContactRelation &relation, Real factor)
    : Integration1stHalf<Contact<Wall>, AcousticRiemannSolver, NoKernelCorrection>(relation)
{
    if (!std::isfinite(factor) || factor < 0 || factor > 1)
        throw std::invalid_argument("wall acceleration pressure factor must be in [0,1]");
    for (auto *body : contact_bodies_)
        acceleration_factors_.push_back(pdWallAccelerationFactor(*body, factor));
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

PDFluidAcousticWorkspace::PDFluidAcousticWorkspace(BaseInnerRelation &inner, BaseContactRelation &wall,
                                                   std::unique_ptr<peridynamics::FluidAcousticBatch> backend, Real acceleration_factor,
                                                   const peridynamics::FluidNeighborConfiguration *configuration)
    : inner_(inner), wall_(wall), acceleration_factor_(acceleration_factor), backend_(std::move(backend)),
      neighbor_cache_(configuration)
{
    if (&inner.getSPHBody() != &wall.getSPHBody())
        throw std::invalid_argument("Fluid acoustic relations must share their source body");
    if (!std::isfinite(acceleration_factor) || acceleration_factor < 0 || acceleration_factor > 1)
        throw std::invalid_argument("wall acceleration pressure factor must be in [0,1]");
    if (!backend_)
        return;

    auto &fluid = dynamic_cast<Fluid &>(inner.getSPHBody().getMatterMaterial());
    const Real density = fluid.ReferenceDensity();
    const Real z = density * fluid.ReferenceSoundSpeed();
    const Real z_squared = math::pow(z, Real(2));
    const Real inverse_impedance = (z + z) / (z_squared + z_squared);
    coefficients_ = {inverse_impedance, 2.0 * z * z * (1.0 / (z + z)),
                     0.5 * (density + density) * inverse_impedance};
    for (auto *body : wall.contact_bodies_)
        wall_acceleration_factors_.push_back(pdWallAccelerationFactor(*body, acceleration_factor));
}

PDFluidPressureHalf::PDFluidPressureHalf(PDFluidAcousticWorkspace &workspace)
    : Base(workspace.inner_, workspace.wall_, workspace.acceleration_factor_), workspace_(workspace) {}

void PDFluidPressureHalf::runMainStep(Real dt)
{
    if (workspace_.onDevice())
        workspace_.evaluateHalf(peridynamics::FluidAcousticStage::Pressure);
    else
        Base::runMainStep(dt);
}

PDFluidDensityHalf::PDFluidDensityHalf(PDFluidAcousticWorkspace &workspace)
    : Base(workspace.inner_, workspace.wall_), workspace_(workspace) {}

void PDFluidDensityHalf::runMainStep(Real dt)
{
    if (workspace_.onDevice())
        workspace_.evaluateHalf(peridynamics::FluidAcousticStage::Density);
    else
        Base::runMainStep(dt);
}

void PDFluidAcousticWorkspace::configure(std::size_t count, std::size_t bound)
{
    using peridynamics::FluidAcousticStencil;
    constexpr auto missing = std::numeric_limits<std::size_t>::max();
    neighbor_cache_.invalidate();
    ghost_indices_.assign(bound - count, missing);
    ghost_sources_.clear();
    wall_sources_.clear();
    const auto resetStencil = [count](FluidAcousticStencil &stencil)
    {
        stencil.offsets.assign(count + 1, 0);
        stencil.neighbors.clear();
    };
    resetStencil(inner_stencil_);
    resetStencil(wall_stencil_);

    wall_indices_.resize(wall_.contact_particles_.size());
    for (std::size_t k = 0; k < wall_indices_.size(); ++k)
        wall_indices_[k].assign(wall_.contact_particles_[k]->ParticlesBound(), missing);

    for (std::size_t i = 0; i < count; ++i)
    {
        const auto &neighbors = inner_.inner_configuration_[i];
        for (std::size_t n = 0; n < neighbors.current_size_; ++n)
        {
            const auto j = neighbors.j_[n];
            auto target = j;
            if (j >= count)
            {
                auto &slot = ghost_indices_.at(j - count);
                if (slot == missing)
                {
                    slot = count + ghost_sources_.size();
                    ghost_sources_.push_back(j);
                }
                target = slot;
            }
            inner_stencil_.neighbors.push_back({target, packFluidVector(neighbors.e_ij_[n]),
                                               neighbors.r_ij_[n], neighbors.dW_ij_[n]});
        }
        inner_stencil_.offsets[i + 1] = inner_stencil_.neighbors.size();

        for (std::size_t k = 0; k < wall_indices_.size(); ++k)
        {
            const auto &contact = wall_.contact_configuration_[k][i];
            for (std::size_t n = 0; n < contact.current_size_; ++n)
            {
                const auto j = contact.j_[n];
                auto &target = wall_indices_[k].at(j);
                if (target == missing)
                {
                    target = wall_sources_.size();
                    wall_sources_.emplace_back(k, j);
                }
                wall_stencil_.neighbors.push_back({target, packFluidVector(contact.e_ij_[n]),
                                                   contact.r_ij_[n], contact.dW_ij_[n]});
            }
        }
        wall_stencil_.offsets[i + 1] = wall_stencil_.neighbors.size();
    }
    backend_->configure(count, count + ghost_sources_.size(), wall_sources_.size(), inner_stencil_, wall_stencil_);
    neighbor_cache_.configured(count);
}

void PDFluidAcousticWorkspace::compute(peridynamics::FluidAcousticStage stage)
{
    using namespace peridynamics;
    const bool pressure = stage == FluidAcousticStage::Pressure;
    const bool density_stage = stage == FluidAcousticStage::Density;
    const bool damping = stage == FluidAcousticStage::Damping;
    const bool gather_velocity = density_stage || stage == FluidAcousticStage::Divergence;
    auto &particles = inner_.getSPHBody().getBaseParticles();
    const auto count = particles.TotalRealParticles();
    if (neighbor_cache_.needsConfiguration(count))
        configure(count, particles.ParticlesBound());

    const auto *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const auto *density = pressure || density_stage ? particles.getVariableDataByName<Real>("Density") : nullptr;
    const auto *p = pressure ? particles.getVariableDataByName<Real>("Pressure") : nullptr;
    const auto *mass = pressure ? particles.getVariableDataByName<Real>("Mass") : nullptr;
    const auto *velocity = gather_velocity ? particles.getVariableDataByName<Vecd>("Velocity") : nullptr;
    const auto *prior = pressure ? particles.getVariableDataByName<Vecd>("ForcePrior") : nullptr;
    const auto *divergence = damping ? particles.getVariableDataByName<Real>("VelocityDivergence") : nullptr;

    fluid_sites_.resize(count + ghost_sources_.size());
    for (std::size_t i = 0; i < count; ++i)
        fluid_sites_[i] = {velocity ? packFluidVector(velocity[i]) : ParticleVector3{},
                           prior ? packFluidVector(prior[i]) : ParticleVector3{}, volume[i],
                           density ? density[i] : Real(0), p ? p[i] : Real(0),
                           mass ? mass[i] : Real(0), divergence ? divergence[i] : Real(0)};
    for (std::size_t n = 0; n < ghost_sources_.size(); ++n)
    {
        const auto j = ghost_sources_[n];
        fluid_sites_[count + n] = {velocity ? packFluidVector(velocity[j]) : ParticleVector3{},
                                   {}, volume[j], 0, p ? p[j] : Real(0), 0,
                                   divergence ? divergence[j] : Real(0)};
    }

    struct WallFields
    {
        const Real *volume;
        const Vecd *velocity, *acceleration, *normal;
    };
    StdVec<WallFields> wall_fields;
    wall_fields.reserve(wall_.contact_particles_.size());
    for (std::size_t k = 0; k < wall_.contact_particles_.size(); ++k)
    {
        auto *wall_particles = wall_.contact_particles_[k];
        auto &solid = dynamic_cast<Solid &>(wall_.contact_bodies_[k]->getMatterMaterial());
        wall_fields.push_back({wall_particles->getVariableDataByName<Real>("VolumetricMeasure"),
                               gather_velocity ? solid.AverageVelocity(wall_particles) : nullptr,
                               pressure ? solid.AverageAcceleration(wall_particles) : nullptr,
                               density_stage ? wall_particles->getVariableDataByName<Vecd>("NormalDirection") : nullptr});
    }
    wall_sites_.resize(wall_sources_.size());
    for (std::size_t n = 0; n < wall_sources_.size(); ++n)
    {
        const auto [k, j] = wall_sources_[n];
        const auto &fields = wall_fields[k];
        wall_sites_[n] = {fields.velocity ? packFluidVector(fields.velocity[j]) : ParticleVector3{},
                          fields.acceleration ? packFluidVector(fields.acceleration[j]) : ParticleVector3{},
                          fields.normal ? packFluidVector(fields.normal[j]) : ParticleVector3{},
                          fields.volume[j], wall_acceleration_factors_[k]};
    }
    backend_->evaluate(stage, fluid_sites_, wall_sites_, coefficients_, results_);
}

void PDFluidAcousticWorkspace::evaluateHalf(peridynamics::FluidAcousticStage stage)
{
    compute(stage);
    auto &particles = inner_.getSPHBody().getBaseParticles();
    auto *force = particles.getVariableDataByName<Vecd>("Force");
    auto *density_rate = particles.getVariableDataByName<Real>("DensityChangeRate");
    const bool pressure = stage == peridynamics::FluidAcousticStage::Pressure;
    for (std::size_t i = 0; i < results_.size(); ++i)
    {
        const auto &row = results_[i];
        if (pressure)
        {
            force[i] += unpackFluidVector(row.inner_force);
            density_rate[i] = row.inner_density_rate;
        }
        else
        {
            force[i] = unpackFluidVector(row.inner_force);
            density_rate[i] += row.inner_density_rate;
        }
        force[i] += unpackFluidVector(row.wall_force);
        density_rate[i] += row.wall_density_rate;
    }
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
    applyDamping(index_i, acoustic_damper, dt);
}

void AcousticDamperInnerforPD::applyDamping(std::size_t index_i, const Vecd &sum, Real dt)
{
    force_[index_i] = mass_[index_i] * coeff_acoustic_damper_ * sum / rho_[index_i];
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
    applyDamping(index_i, acoustic_damper, dt);
}

void AcousticDamperWithWallforPD::applyDamping(std::size_t index_i, const Vecd &sum, Real dt)
{
    force_[index_i] = mass_[index_i] * coeff_acoustic_damper_ * sum / rho_[index_i];
    if (!force_[index_i].allFinite() || !std::isfinite(dt) || dt < 0.0)
        peridynamics::numericalFailure(sph_body_->Name(), index_i, "AcousticDampingForce");
    work_[index_i] += dt * force_[index_i].dot(vel_[index_i] + 0.5 * dt * force_[index_i] / mass_[index_i]);
    vel_[index_i] += dt * force_[index_i] / mass_[index_i];
}
//=================================================================================================//
PDFluidVelocityDivergence::PDFluidVelocityDivergence(PDFluidAcousticWorkspace &workspace)
    : workspace_(workspace), inner_(workspace.inner_), wall_(workspace.wall_) {}

void PDFluidVelocityDivergence::exec(Real dt)
{
    if (!workspace_.onDevice())
    {
        inner_.exec(dt);
        wall_.exec(dt);
        return;
    }

    workspace_.compute(peridynamics::FluidAcousticStage::Divergence);
    auto &body = workspace_.inner_.getSPHBody();
    auto *divergence = body.getBaseParticles().getVariableDataByName<Real>("VelocityDivergence");
    body.setNewlyUpdated();
    for (std::size_t i = 0; i < workspace_.results_.size(); ++i)
    {
        const auto &row = workspace_.results_[i];
        divergence[i] = -row.inner_density_rate - row.wall_density_rate;
    }
}

PDFluidAcousticDamping::PDFluidAcousticDamping(PDFluidAcousticWorkspace &workspace)
    : workspace_(workspace), inner_(workspace.inner_), wall_(workspace.wall_) {}

void PDFluidAcousticDamping::setCoeffAcousticDamper(Real rho0, Real c0, Real h)
{
    inner_.setCoeffAcousticDamper(rho0, c0, h);
    wall_.setCoeffAcousticDamper(rho0, c0, h);
}

void PDFluidAcousticDamping::exec(Real dt)
{
    if (!workspace_.onDevice())
    {
        inner_.exec(dt);
        wall_.exec(dt);
        return;
    }

    workspace_.compute(peridynamics::FluidAcousticStage::Damping);
    workspace_.inner_.getSPHBody().setNewlyUpdated();
    for (std::size_t i = 0; i < workspace_.results_.size(); ++i)
    {
        const auto &row = workspace_.results_[i];
        inner_.applyDamping(i, unpackFluidVector(row.inner_force), dt);
        wall_.applyDamping(i, unpackFluidVector(row.wall_force), dt);
    }
}
} // namespace fluid_dynamics

namespace solid_dynamics
{
void ViscousForceFromFluidWithExecution::evaluateDevice(peridynamics::FluidViscosityBatch &backend,
                                                        peridynamics::SPHViscosityWorkspace &workspace)
{
    workspace.clearTargets();
    for (std::size_t k = 0; k < contact_configuration_.size(); ++k)
        workspace.addTarget(peridynamics::SPHViscosityWorkspace::Form::ReverseWall, *contact_configuration_[k],
                            contact_vel_[k], contact_Vol_[k], contact_particles_[k]->ParticlesBound(),
                            mu_[k], smoothing_length_[k]);

    workspace.evaluate(backend, particles_->TotalRealParticles(), vel_ave_, Vol_, force_from_fluid_);
}

//=================================================================================================//
template <class RiemannSolverType>
PressureForceFromFluidforPD<RiemannSolverType>::
    PressureForceFromFluidforPD(BaseContactRelation &contact_relation,
                                const peridynamics::FluidNeighborConfiguration *configuration)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      hw_solid_(DynamicCast<HughesWingetSolid>(this, contact_relation.getSPHBody().getMatterMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      vel_ave_(particles_->registerStateVariableData<Vecd>("AverageVelocity")),
      acc_prior_(particles_->registerStateVariableData<Vecd>("PDAccelerationPrior")),
      acc_ave_(particles_->registerStateVariableData<Vecd>("AverageAcceleration")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      force_from_fluid_(particles_->registerStateVariableData<Vecd>("PressureForceFromFluid")),
      acoustic_force_from_fluid_(particles_->getVariableDataByName<Vecd>("AcousticDampingForceFromFluid")),
      pressure_from_fluid_(particles_->registerStateVariableData<Real>("PressureFromFluid")),
      neighbor_cache_(configuration)
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

template<class RiemannSolverType>
void PressureForceFromFluidforPD<RiemannSolverType>::configureDevice(
    peridynamics::FluidTractionBatch &backend, std::size_t count)
{
    neighbor_cache_.invalidate();
    row_offsets_.assign(count + 1, 0);
    traction_neighbors_.clear();
    fluid_sources_.clear();
    fluid_indices_.resize(contact_particles_.size());
    constexpr auto missing = std::numeric_limits<std::size_t>::max();
    for (std::size_t k = 0; k < contact_particles_.size(); ++k)
        fluid_indices_[k].assign(contact_particles_[k]->ParticlesBound(), missing);

    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const auto &neighbors = (*contact_configuration_[k])[i];
            for (std::size_t n = 0; n < neighbors.current_size_; ++n)
            {
                const auto j = neighbors.j_[n];
                auto &fluid_id = fluid_indices_[k].at(j);
                if (fluid_id == missing)
                {
                    fluid_id = fluid_sources_.size();
                    fluid_sources_.emplace_back(k, j);
                }
                traction_neighbors_.push_back({fluid_id, packFluidVector(neighbors.e_ij_[n]),
                    neighbors.r_ij_[n], neighbors.W_ij_[n], neighbors.dW_ij_[n]});
            }
        }
        row_offsets_[i + 1] = traction_neighbors_.size();
    }

    backend.configure(count, fluid_sources_.size(), row_offsets_, traction_neighbors_);
    configured_backend_ = &backend;
    neighbor_cache_.configured(count);
}

template<class RiemannSolverType>
void PressureForceFromFluidforPD<RiemannSolverType>::evaluateDevice(peridynamics::FluidTractionBatch &backend)
{
    const auto count = particles_->TotalRealParticles();
    if (neighbor_cache_.needsConfiguration(count) || configured_backend_ != &backend)
        configureDevice(backend, count);

    solid_rows_.resize(count);
    for (std::size_t i = 0; i < count; ++i)
        solid_rows_[i] = {packFluidVector(n_[i]), packFluidVector(vel_ave_[i]),
                          packFluidVector(acc_ave_[i]), Vol_[i]};

    fluid_sites_.resize(fluid_sources_.size());
    for (std::size_t n = 0; n < fluid_sources_.size(); ++n)
    {
        const auto [k, j] = fluid_sources_[n];
        const auto [impedance, inverse_sound] = acoustic_coefficients_[k];
        fluid_sites_[n] = {packFluidVector(contact_vel_[k][j]), packFluidVector(contact_force_prior_[k][j]),
                           contact_mass_[k][j], contact_rho_[k][j], contact_p_[k][j],
                           contact_u_div_[k][j], contact_Vol_[k][j], impedance, inverse_sound};
    }
    backend.evaluate(solid_rows_, fluid_sites_, acceleration_factor_, traction_result_);
    if (traction_result_.size()!=count)
        throw std::logic_error("PD fluid traction backend returned an incomplete force field");
    for (size_t i=0;i<count;++i)
    {
        const auto &result=traction_result_[i];
        publishTraction(i,unpackFluidVector(result.pressure_force),unpackFluidVector(result.acoustic_force),
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

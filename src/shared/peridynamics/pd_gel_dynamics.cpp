#include "pd_gel_dynamics.h"
#include "base_particles.hpp"
#include "pd_numerics.h"

namespace SPH::solid_dynamics
{
GelConstitutiveUpdate::GelConstitutiveUpdate(BaseInnerRelation &relation, Stage stage)
    : LocalDynamics(relation.getSPHBody()), DataDelegateInner(relation), stage_(stage),
      material_(DynamicCast<GelViscoelasticSolid>(this, relation.getSPHBody().getMatterMaterial())),
      position_(particles_->getVariableDataByName<Vecd>("Position")),
      volume_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      density_(particles_->registerStateVariableData<Real>("Density")),
      dissipation_(particles_->getVariableDataByName<Real>("GelDissipationDensity")),
      midpoint_dissipation_(particles_->getVariableDataByName<Real>("GelMidpointDissipationDensity")),
      energy_(particles_->getVariableDataByName<Real>("GelEnergyDensity")),
      history_time_(particles_->getVariableDataByName<Real>("GelHistoryTime")),
      inverse_shape_(particles_->getVariableDataByName<Matd>("ShapeTensorKInv")),
      deformation_(particles_->registerStateVariableData<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      stress_(particles_->getVariableDataByName<Matd>("PDStress")),
      pk1_(particles_->getVariableDataByName<Matd>("PDPK1")),
      force_state_(particles_->getVariableDataByName<Matd>("PDForceStateT0")),
      accepted_stress_(particles_->getVariableDataByName<Matd>("GelAcceptedStress")),
      green_(particles_->getVariableDataByName<Matd>("GelGreenStrain")),
      viscous_(particles_->getVariableDataByName<Matd>("GelViscousStrain")),
      midpoint_green_(particles_->getVariableDataByName<Matd>("GelMidpointGreenStrain")),
      midpoint_viscous_(particles_->getVariableDataByName<Matd>("GelMidpointViscousStrain")) {}

void GelConstitutiveUpdate::interaction(size_t i, Real dt)
{
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("staged gel integration requires a positive time increment");
    Matd moment = Matd::Zero();
    const auto &neighbors = inner_configuration_[i];
    for (size_t n = 0; n < neighbors.current_size_; ++n)
    {
        if (!neighbors.bondLive_[n])
            peridynamics::numericalFailure(sph_body_->Name(), i, "gel control does not support broken reference bonds");
        const size_t j = neighbors.j_[n];
        const Vecd reference = -neighbors.r_ij_[n] * neighbors.e_ij_[n];
        moment += neighbors.W_ij_[n] * volume_[j] * (position_[j] - position_[i]) * reference.transpose();
    }
    const Matd F = moment * inverse_shape_[i];
    const bool midpoint = stage_ == Stage::MidpointForce;
    const peridynamics::GelHistory old = midpoint ? peridynamics::GelHistory{green_[i], viscous_[i], dissipation_[i]} :
        peridynamics::GelHistory{midpoint_green_[i], midpoint_viscous_[i], midpoint_dissipation_[i]};
    const auto next = peridynamics::advanceGel(F, 0.5 * dt, material_.gelParameters(), old);
    const auto response = peridynamics::gelResponse(F, material_.gelParameters(), next);
    deformation_[i] = F;
    if (midpoint)
    {
        midpoint_green_[i] = next.green;
        midpoint_viscous_[i] = next.viscous;
        midpoint_dissipation_[i] = next.dissipation;
        stress_[i] = response.cauchy;
        pk1_[i] = response.pk1;
        force_state_[i] = response.pk1 * inverse_shape_[i];
    }
    else
    {
        green_[i] = next.green;
        viscous_[i] = next.viscous;
        dissipation_[i] = next.dissipation;
        energy_[i] = response.free_energy;
        accepted_stress_[i] = response.cauchy;
        density_[i] = material_.ReferenceDensity() / F.determinant();
        history_time_[i] += dt;
    }
}
} // namespace SPH::solid_dynamics

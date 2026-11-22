#include "pd_gel_dynamics.h"
#include "base_particles.hpp"
#include "pd_numerics.h"
#include <algorithm>

namespace SPH::solid_dynamics
{
void canonicalizeInitialGelReference(BaseInnerRelation &relation, Real time)
{
    auto &particles = relation.getSPHBody().getBaseParticles();
    if (time != 0)
        throw std::invalid_argument("gel reference ordering is only initialized at fresh time zero");
    const size_t count = particles.TotalRealParticles();
    const Real *history_time = particles.getVariableDataByName<Real>("GelHistoryTime");
    const Matd *green = particles.getVariableDataByName<Matd>("GelGreenStrain");
    const Matd *viscous = particles.getVariableDataByName<Matd>("GelViscousStrain");
    const Matd *shape = particles.getVariableDataByName<Matd>("ShapeTensorK");
    for (size_t i = 0; i < count; ++i)
        if (history_time[i] != 0 || green[i].squaredNorm() != 0 || viscous[i].squaredNorm() != 0 || shape[i].squaredNorm() != 0)
            throw std::invalid_argument("cannot reorder an evolved or already constructed gel reference");
    struct Packet { size_t j; Real W, dW, r; Vecd e; bool live; };
    for (size_t i = 0; i < count; ++i)
    {
        auto &row = relation.inner_configuration_[i];
        StdVec<Packet> packets;
        packets.reserve(row.current_size_);
        for (size_t n = 0; n < row.current_size_; ++n)
            packets.push_back({row.j_[n], row.W_ij_[n], row.dW_ij_[n], row.r_ij_[n], row.e_ij_[n], row.bondLive_[n]});
        std::sort(packets.begin(), packets.end(), [](const auto &a, const auto &b) { return a.j < b.j; });
        for (size_t n = 0; n < packets.size(); ++n)
        {
            if (packets[n].j >= count || !packets[n].live || (n && packets[n-1].j == packets[n].j))
                throw std::invalid_argument("invalid initial gel reference packet");
            const auto &p = packets[n];
            row.j_[n] = p.j; row.W_ij_[n] = p.W; row.dW_ij_[n] = p.dW;
            row.r_ij_[n] = p.r; row.e_ij_[n] = p.e; row.bondLive_[n] = p.live;
        }
    }
}

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

namespace
{
void validateStagedTime(Real dt)
{
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("staged gel integration requires a positive time increment");
}

peridynamics::ShapeHessian3 packGelTensor(const Matd &tensor)
{
    peridynamics::ShapeHessian3 packed{};
    for (int a = 0; a < Dimensions; ++a)
        for (int b = 0; b < Dimensions; ++b)
            packed.derivative[a][b] = tensor(a,b);
    return packed;
}

Matd unpackGelTensor(const peridynamics::ShapeHessian3 &tensor)
{
    Matd value;
    for (int a = 0; a < Dimensions; ++a)
        for (int b = 0; b < Dimensions; ++b)
            value(a,b) = tensor.derivative[a][b];
    return value;
}
}

Matd GelConstitutiveUpdate::deformationAt(size_t i) const
{
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
    return moment * inverse_shape_[i];
}

peridynamics::GelHistory GelConstitutiveUpdate::historyAt(size_t i) const
{
    if (stage_ == Stage::MidpointForce)
        return {green_[i], viscous_[i], dissipation_[i]};
    return {midpoint_green_[i], midpoint_viscous_[i], midpoint_dissipation_[i]};
}

void GelConstitutiveUpdate::interaction(size_t i, Real dt)
{
    validateStagedTime(dt);
    const Matd F = deformationAt(i);
    const auto next = peridynamics::advanceGel(F, 0.5 * dt, material_.gelParameters(), historyAt(i));
    const auto response = peridynamics::gelResponse(F, material_.gelParameters(), next);
    publish(i, dt, F, next, response);
}

void GelConstitutiveUpdate::evaluateDevice(peridynamics::GelPointBatch &backend, Real dt)
{
    using namespace peridynamics;
    validateStagedTime(dt);
    const size_t count = particles_->TotalRealParticles();
    device_points_.resize(count);
    tbb::parallel_for(size_t(0), count, [&](size_t i)
    {
        const auto old = historyAt(i);
        device_points_[i] = {packGelTensor(deformationAt(i)), packGelTensor(old.green),
                             packGelTensor(old.viscous), old.dissipation};
    });
    const auto &p = material_.gelParameters();
    backend.evaluate(device_points_, {p.bulk, p.equilibrium_shear, p.branch_shear, p.relaxation_time},
                     0.5 * dt, device_results_);
    if (device_results_.size() != count)
        throw std::logic_error("Gel material backend returned an incomplete batch");
    for (const auto &result : device_results_)
        requireGelStatus(result.status);

    tbb::parallel_for(size_t(0), count, [&](size_t i)
    {
        const auto &result = device_results_[i];
        const GelHistory next{unpackGelTensor(result.green), unpackGelTensor(result.viscous), result.dissipation};
        const GelResponse response{unpackGelTensor(result.pk1), unpackGelTensor(result.cauchy), result.free_energy};
        publish(i, dt, unpackGelTensor(device_points_[i].deformation), next, response);
    });
}

void GelConstitutiveUpdate::publish(size_t i, Real dt, const Matd &F, const peridynamics::GelHistory &next,
                                    const peridynamics::GelResponse &response)
{
    deformation_[i] = F;
    if (stage_ == Stage::MidpointForce)
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

/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_ELASTIC_NORMAL_STEP_H
#define PD_ELASTIC_NORMAL_STEP_H

#include "pd_normal_projection.h"

namespace SPH::peridynamics
{
/** Reduced diagonal-spring research step. No PD body or constitutive history.
 * A=M+dt^2 K/2 is an effective metric, never a change to physical mass M.
 * Forces use midpoint displacement; velocity/position are accepted endpoints.
 * A nonconverged result must not be committed by a caller.
 */
struct ElasticNormalStepResult
{
    NormalMetricProjectionResult projection;
    StdVec<Real> elastic_force;
    Real kinetic_old = 0, kinetic_new = 0, elastic_old = 0, elastic_new = 0;
    Real algorithmic_loss = 0, wall_work = 0;
    Real energy_residual = 0, maximum_momentum_residual = 0;
};

inline ElasticNormalStepResult advanceElasticNormal(const StdVec<Real> &mass,
    const StdVec<Real> &stiffness, const StdVec<Real> &position, const StdVec<Real> &velocity,
    Real dt, const StdVec<NormalContactConstraint> &constraints, const NormalProjectionSettings &settings = {})
{
    const size_t n = mass.size();
    if (!n || n > 128 || stiffness.size()!=n || position.size()!=n || velocity.size()!=n ||
        !std::isfinite(dt) || dt<=0)
        throw std::invalid_argument("elastic normal step: invalid dimensions or dt");
    StdVec<Real> metric(n), free_velocity(n);
    for (size_t i=0; i<n; ++i)
    {
        if (!std::isfinite(mass[i]) || mass[i]<=0 || !std::isfinite(stiffness[i]) || stiffness[i]<0 ||
            !std::isfinite(position[i]) || !std::isfinite(velocity[i]))
            throw std::invalid_argument("elastic normal step: invalid physical state");
        metric[i] = mass[i] + 0.5*dt*dt*stiffness[i];
        free_velocity[i] = stiffness[i] == 0 ? velocity[i] :
            (mass[i]*velocity[i]-dt*stiffness[i]*position[i])/metric[i];
    }
    ElasticNormalStepResult r;
    r.projection = projectNormalMetric(metric, position, free_velocity, dt, constraints, settings);
    r.elastic_force.resize(n);
    StdVec<Real> residual(n);
    for (size_t i=0; i<n; ++i)
    {
        const Real q = r.projection.position[i], v = r.projection.velocity[i], dv = v-velocity[i];
        r.elastic_force[i] = -0.5*stiffness[i]*(position[i]+q);
        residual[i] = mass[i]*dv-dt*r.elastic_force[i];
        r.kinetic_old += 0.5*mass[i]*velocity[i]*velocity[i];
        r.kinetic_new += 0.5*mass[i]*v*v;
        r.elastic_old += 0.5*stiffness[i]*position[i]*position[i];
        r.elastic_new += 0.5*stiffness[i]*q*q;
        r.algorithmic_loss += 0.5*mass[i]*dv*dv;
    }
    for (size_t s=0; s<constraints.size(); ++s)
    {
        const auto &row = constraints[s];
        Real old_gap = -row.old_offset;
        for (size_t j=0; j<row.targets.size(); ++j)
        {
            const size_t i=row.targets[j];
            old_gap += row.normal*row.weights[j]*position[i];
            residual[i] -= row.normal*row.weights[j]*r.projection.normal_impulse[s];
        }
        const Real impulse=r.projection.normal_impulse[s];
        r.wall_work += impulse*(row.new_offset-row.old_offset)/dt;
        r.algorithmic_loss += impulse*old_gap/dt;
    }
    for (Real value : residual) r.maximum_momentum_residual=std::max(r.maximum_momentum_residual,std::abs(value));
    r.energy_residual = r.wall_work-(r.kinetic_new-r.kinetic_old+r.elastic_new-r.elastic_old)-r.algorithmic_loss;
    for (Real value : {r.kinetic_old,r.kinetic_new,r.elastic_old,r.elastic_new,r.algorithmic_loss,
                       r.wall_work,r.energy_residual,r.maximum_momentum_residual})
        if (!std::isfinite(value)) throw std::runtime_error("elastic normal step: non-finite physical ledger");
    return r;
}
} // namespace SPH::peridynamics
#endif

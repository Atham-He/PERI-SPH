/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_NORMAL_PROJECTION_H
#define PD_NORMAL_PROJECTION_H

#include "base_data_type.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Candidate normal-only contact mathematics. It is NOT connected to the PD
 * integrator. y_new=y_old+dt*v_new, C*v_new>=b. Plate normal is +/-1; signed
 * reference interpolation weights are retained, with sum(weights)=1. */
struct NormalContactConstraint
{
    StdVec<size_t> targets;
    StdVec<Real> weights;
    Real normal = 1, old_offset = 0, new_offset = 0;
};

struct NormalProjectionSettings
{
    Real velocity_scale = 1;
    Real relative_tolerance = 1e-12;
    size_t max_sweeps = 10000;
};

// Quadratic forms below use the supplied mathematical metric A. They are
// physical kinetic energies ONLY when the caller explicitly supplies mass M.
struct NormalMetricProjectionResult
{
    bool converged = false;
    size_t sweeps = 0;
    StdVec<Real> velocity, position, normal_impulse;
    Real primal_velocity_error = 0, natural_velocity_error = 0;
    Real stationarity_impulse_error = 0, complementarity_energy_error = 0;
    Real free_quadratic = 0, accepted_quadratic = 0, wall_work = 0;
    Real metric_step_loss = 0, energy_residual = 0, stationarity_sum_residual = 0;
};

/** Solve min 1/2*(v-v_free)' A (v-v_free), C v >= b by projected dual
 * coordinate descent. The exact primal minimizer is unique, but redundant
 * contacts need not have a unique distribution of dual impulses.
 * All inputs are immutable. A nonconverged result MUST NOT be committed.
 * This bounded normal-only prototype contains no friction/restitution fit,
 * warm-start history, stiffness penalty, position clamp or production fallback.
 */
inline NormalMetricProjectionResult projectNormalMetric(const StdVec<Real> &metric, const StdVec<Real> &position,
    const StdVec<Real> &free_velocity, Real dt, const StdVec<NormalContactConstraint> &constraints,
    const NormalProjectionSettings &settings = {})
{
    const size_t count = metric.size(), rows = constraints.size();
    if (count == 0 || count > 128 || rows > 1024 || position.size() != count || free_velocity.size() != count ||
        !std::isfinite(dt) || dt <= 0 || !std::isfinite(settings.velocity_scale) || settings.velocity_scale <= 0 ||
        !std::isfinite(settings.relative_tolerance) || settings.relative_tolerance <= 0 || settings.relative_tolerance > 1e-3 ||
        settings.max_sweeps == 0 || settings.max_sweeps > 100000)
        throw std::invalid_argument("normal contact prototype: invalid dimensions, step or settings");
    Real total_metric = 0;
    for (size_t i = 0; i < count; ++i)
    {
        if (!std::isfinite(metric[i]) || metric[i] <= 0 || !std::isfinite(position[i]) || !std::isfinite(free_velocity[i]))
            throw std::invalid_argument("normal contact prototype: invalid body state");
        total_metric += metric[i];
    }
    const Real scale_v = settings.velocity_scale, scale_p = total_metric * scale_v, scale_e = scale_p * scale_v;
    if (!std::isfinite(total_metric) || !std::isfinite(scale_p) || !std::isfinite(scale_e) || scale_e <= 0)
        throw std::invalid_argument("normal contact prototype: invalid metric scales");
    StdVec<Real> m(count), velocity(count), lambda(rows, 0), diagonal(rows, 0), bound(rows), old_gap(rows), wall_speed(rows);
    for (size_t i = 0; i < count; ++i)
    {
        m[i] = metric[i] / total_metric; velocity[i] = free_velocity[i] / scale_v;
        if (m[i] <= 0 || !std::isfinite(velocity[i]))
            throw std::invalid_argument("normal contact prototype: ill-scaled metric or velocity");
    }
    auto rowValue = [&](size_t s, const StdVec<Real> &values)
    {
        Real value = 0;
        const auto &row = constraints[s];
        for (size_t j = 0; j < row.targets.size(); ++j) value += row.normal * row.weights[j] * values[row.targets[j]];
        return value;
    };
    for (size_t s = 0; s < rows; ++s)
    {
        const auto &row = constraints[s];
        if (row.targets.empty() || row.targets.size() != row.weights.size() || !std::isfinite(row.normal) ||
            std::abs(row.normal) != 1 || !std::isfinite(row.old_offset) || !std::isfinite(row.new_offset))
            throw std::invalid_argument("normal contact prototype: invalid constraint");
        Real partition = 0;
        for (size_t j = 0; j < row.targets.size(); ++j)
        {
            if (row.targets[j] >= count || !std::isfinite(row.weights[j]) ||
                std::find(row.targets.begin(), row.targets.begin() + j, row.targets[j]) != row.targets.begin() + j)
                throw std::invalid_argument("normal contact prototype: invalid/repeated target");
            partition += row.weights[j];
            diagonal[s] += row.weights[j] * row.weights[j] / m[row.targets[j]];
        }
        if (std::abs(partition - 1) > 1e-10 || !std::isfinite(diagonal[s]) || diagonal[s] <= 0)
            throw std::invalid_argument("normal contact prototype: partition or metric rank failure");
        const Real old_normal_position = rowValue(s, position);
        old_gap[s] = (old_normal_position - row.old_offset) / (dt * scale_v);
        bound[s] = (row.new_offset - old_normal_position) / (dt * scale_v);
        wall_speed[s] = (row.new_offset - row.old_offset) / (dt * scale_v);
        if (!std::isfinite(old_gap[s]) || !std::isfinite(bound[s]) || !std::isfinite(wall_speed[s]) ||
            old_gap[s] < -settings.relative_tolerance)
            throw std::invalid_argument("normal contact prototype: non-finite or initially penetrating geometry");
    }
    NormalMetricProjectionResult result;
    result.converged = rows == 0;
    for (size_t sweep = 0; rows && sweep < settings.max_sweeps; ++sweep)
    {
        for (size_t s = 0; s < rows; ++s)
        {
            const Real trial = lambda[s] + (bound[s] - rowValue(s, velocity)) / diagonal[s];
            if (!std::isfinite(trial)) throw std::runtime_error("normal contact prototype: non-finite dual iterate");
            const Real next = std::max(Real(0), trial), change = next - lambda[s];
            lambda[s] = next;
            const auto &row = constraints[s];
            for (size_t j = 0; j < row.targets.size(); ++j)
            {
                const size_t i = row.targets[j];
                velocity[i] += change * row.normal * row.weights[j] / m[i];
                if (!std::isfinite(velocity[i])) throw std::runtime_error("normal contact prototype: non-finite primal iterate");
            }
        }
        Real primal = 0, natural = 0, complementarity = 0, energy_scale = 1;
        for (size_t s = 0; s < rows; ++s)
        {
            const Real gap = rowValue(s, velocity) - bound[s];
            primal = std::max(primal, -gap);
            natural = std::max(natural, std::abs(std::min(lambda[s] * diagonal[s], gap)));
            complementarity += std::abs(lambda[s] * gap);
            energy_scale += std::abs(lambda[s] * wall_speed[s]);
        }
        StdVec<Real> stationarity(count);
        for (size_t i = 0; i < count; ++i) stationarity[i] = m[i] * (velocity[i] - free_velocity[i] / scale_v);
        for (size_t s = 0; s < rows; ++s)
            for (size_t j = 0; j < constraints[s].targets.size(); ++j)
                stationarity[constraints[s].targets[j]] -= lambda[s] * constraints[s].normal * constraints[s].weights[j];
        Real stationarity_error = 0;
        for (Real value : stationarity) stationarity_error = std::max(stationarity_error, std::abs(value));
        result.sweeps = sweep + 1;
        result.primal_velocity_error = primal * scale_v;
        result.natural_velocity_error = natural * scale_v;
        result.stationarity_impulse_error = stationarity_error * scale_p;
        result.complementarity_energy_error = complementarity * scale_e;
        result.converged = primal <= settings.relative_tolerance && natural <= settings.relative_tolerance &&
            stationarity_error <= settings.relative_tolerance && complementarity <= settings.relative_tolerance * energy_scale;
        if (result.converged) break;
    }
    result.velocity.resize(count); result.position.resize(count); result.normal_impulse.resize(rows);
    for (size_t i = 0; i < count; ++i)
    {
        result.velocity[i] = velocity[i] * scale_v;
        result.position[i] = position[i] + dt * result.velocity[i];
        result.free_quadratic += 0.5 * metric[i] * free_velocity[i] * free_velocity[i];
        result.accepted_quadratic += 0.5 * metric[i] * result.velocity[i] * result.velocity[i];
        const Real change = result.velocity[i] - free_velocity[i];
        result.metric_step_loss += 0.5 * metric[i] * change * change;
        result.stationarity_sum_residual += metric[i] * change;
    }
    for (size_t s = 0; s < rows; ++s)
    {
        const Real impulse = lambda[s] * scale_p;
        result.normal_impulse[s] = impulse;
        result.wall_work += impulse * wall_speed[s] * scale_v;
        result.metric_step_loss += impulse * old_gap[s] * scale_v;
        result.stationarity_sum_residual -= constraints[s].normal * impulse;
    }
    result.energy_residual = result.wall_work - (result.accepted_quadratic - result.free_quadratic) - result.metric_step_loss;
    for (Real value : result.position) if (!std::isfinite(value)) throw std::runtime_error("normal contact prototype: non-finite accepted position");
    for (Real value : {result.free_quadratic, result.accepted_quadratic, result.wall_work, result.metric_step_loss,
                       result.energy_residual, result.stationarity_sum_residual})
        if (!std::isfinite(value)) throw std::runtime_error("normal contact prototype: non-finite energy or impulse");
    return result;
}

/** Physical-mass specialization. Keeps the original mass/kinetic-energy API;
 * implicit elastic callers must use projectNormalMetric and a separate ledger. */
struct NormalProjectionResult
{
    bool converged = false;
    size_t sweeps = 0;
    StdVec<Real> velocity, position, normal_impulse;
    Real primal_velocity_error = 0, natural_velocity_error = 0;
    Real stationarity_impulse_error = 0, complementarity_energy_error = 0;
    Real free_kinetic = 0, accepted_kinetic = 0, wall_work = 0;
    Real contact_step_dissipation = 0, energy_residual = 0, momentum_residual = 0;
};

inline NormalProjectionResult projectNormalContact(const StdVec<Real> &mass, const StdVec<Real> &position,
    const StdVec<Real> &free_velocity, Real dt, const StdVec<NormalContactConstraint> &constraints,
    const NormalProjectionSettings &settings = {})
{
    const auto r = projectNormalMetric(mass, position, free_velocity, dt, constraints, settings);
    NormalProjectionResult out;
    out.converged=r.converged; out.sweeps=r.sweeps;
    out.velocity=r.velocity; out.position=r.position; out.normal_impulse=r.normal_impulse;
    out.primal_velocity_error=r.primal_velocity_error; out.natural_velocity_error=r.natural_velocity_error;
    out.stationarity_impulse_error=r.stationarity_impulse_error;
    out.complementarity_energy_error=r.complementarity_energy_error;
    out.free_kinetic=r.free_quadratic; out.accepted_kinetic=r.accepted_quadratic;
    out.wall_work=r.wall_work; out.contact_step_dissipation=r.metric_step_loss;
    out.energy_residual=r.energy_residual; out.momentum_residual=r.stationarity_sum_residual;
    return out;
}

} // namespace SPH::peridynamics
#endif

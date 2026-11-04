/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_DENSE_METRIC_CONTACT_H
#define PD_DENSE_METRIC_CONTACT_H

#include "pd_normal_projection.h"

namespace SPH::peridynamics
{
using DenseContactMatrix = StdVec<StdVec<Real>>;

/** A row C may represent signed reference interpolation and any flattened
 * normal. The identical row is used for gather and transpose impulse scatter.
 * Offsets have the same units as C*q. No diagonal replacement is permitted. */
struct DenseContactConstraint
{
    StdVec<Real> coefficients;
    Real old_offset = 0, new_offset = 0;
};

/** Bounded dense SPD factorization for research controls, not a scalable PD
 * solver. Normalization affects only arithmetic, never the physical operator.
 * Non-SPD or nonsymmetric inputs fail; there is no diagonal shift/fallback. */
class DensePositiveMetric
{
  public:
    struct Solution { StdVec<Real> value; Real relative_residual; };

    explicit DensePositiveMetric(const DenseContactMatrix &matrix) : matrix_(matrix)
    {
        const size_t n = matrix.size();
        if (!n || n > 128) throw std::invalid_argument("dense metric: invalid dimensions");
        long double trace = 0, norm = 0;
        for (size_t i = 0; i < n; ++i)
        {
            if (matrix[i].size() != n) throw std::invalid_argument("dense metric: non-square matrix");
            long double row = 0;
            for (Real value : matrix[i])
            {
                if (!std::isfinite(value)) throw std::invalid_argument("dense metric: non-finite matrix");
                row += std::abs(static_cast<long double>(value));
            }
            norm = std::max(norm, row); trace += matrix[i][i];
        }
        if (!(trace > 0) || !(norm > 0) || !std::isfinite(trace) || !std::isfinite(norm))
            throw std::invalid_argument("dense metric: nonpositive trace/norm");
        scale_ = trace; norm_ = norm;
        lower_.assign(n, StdVec<long double>(n, 0));
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j <= i; ++j)
            {
                if (std::abs(static_cast<long double>(matrix[i][j]) - matrix[j][i]) > 1e-14L * norm)
                    throw std::invalid_argument("dense metric: nonsymmetric matrix");
                long double value = matrix[i][j] / scale_;
                for (size_t k = 0; k < j; ++k) value -= lower_[i][k] * lower_[j][k];
                if (i == j)
                {
                    if (!(value > 1e-14L * norm_ / scale_))
                        throw std::invalid_argument("dense metric: non-SPD or numerically rank-deficient pivot");
                    lower_[i][j] = std::sqrt(value);
                }
                else lower_[i][j] = value / lower_[j][j];
            }
    }

    Solution solve(const StdVec<Real> &rhs) const
    {
        const size_t n = matrix_.size();
        if (rhs.size() != n) throw std::invalid_argument("dense metric: rhs size mismatch");
        StdVec<long double> value(n);
        for (size_t i = 0; i < n; ++i)
        {
            if (!std::isfinite(rhs[i])) throw std::invalid_argument("dense metric: non-finite rhs");
            value[i] = rhs[i] / scale_;
            for (size_t j = 0; j < i; ++j) value[i] -= lower_[i][j] * value[j];
            value[i] /= lower_[i][i];
        }
        for (size_t i = n; i-- > 0;)
        {
            for (size_t j = i + 1; j < n; ++j) value[i] -= lower_[j][i] * value[j];
            value[i] /= lower_[i][i];
        }
        Solution result; result.value.resize(n);
        long double max_x = 0, max_b = 0, error = 0;
        for (size_t i = 0; i < n; ++i)
        {
            result.value[i] = value[i];
            if (!std::isfinite(result.value[i])) throw std::runtime_error("dense metric: non-finite solution");
            max_x = std::max(max_x, std::abs(static_cast<long double>(result.value[i])));
            max_b = std::max(max_b, std::abs(static_cast<long double>(rhs[i])));
        }
        for (size_t i = 0; i < n; ++i)
        {
            long double residual = -rhs[i];
            for (size_t j = 0; j < n; ++j) residual += static_cast<long double>(matrix_[i][j]) * result.value[j];
            error = std::max(error, std::abs(residual));
        }
        const long double denominator = norm_ * max_x + max_b;
        result.relative_residual = denominator > 0 ? error / denominator : 0;
        if (result.relative_residual > 1e-12) throw std::runtime_error("dense metric: linear residual gate failed");
        return result;
    }

    long double scale() const { return scale_; }

  private:
    DenseContactMatrix matrix_;
    StdVec<StdVec<long double>> lower_;
    long double scale_, norm_;
};

inline Real denseContactDot(const StdVec<Real> &a, const StdVec<Real> &b)
{
    long double value = 0;
    for (size_t i = 0; i < a.size(); ++i) value += static_cast<long double>(a[i]) * b[i];
    return value;
}

struct DenseMetricContactResult
{
    bool converged = false;
    size_t sweeps = 0;
    StdVec<Real> velocity, normal_impulse;
    Real linear_relative_residual = 0, primal_velocity_error = 0, natural_velocity_error = 0;
    Real stationarity_impulse_error = 0, complementarity_energy_error = 0;
};

inline DenseMetricContactResult projectDenseMetricContact(const DenseContactMatrix &metric,
    const StdVec<Real> &position, const StdVec<Real> &free_velocity, Real dt,
    const StdVec<DenseContactConstraint> &constraints, const NormalProjectionSettings &settings = {})
{
    const DensePositiveMetric factor(metric);
    const size_t n = metric.size(), rows = constraints.size();
    const Real speed = settings.velocity_scale;
    if (position.size() != n || free_velocity.size() != n || rows > 1024 || !std::isfinite(dt) || dt <= 0 ||
        !std::isfinite(speed) || speed <= 0 || !std::isfinite(settings.relative_tolerance) ||
        settings.relative_tolerance <= 0 || settings.relative_tolerance > 1e-3 || !settings.max_sweeps || settings.max_sweeps > 100000)
        throw std::invalid_argument("dense contact: invalid state/settings");
    const Real impulse_scale = factor.scale() * speed, energy_scale = impulse_scale * speed;
    if (!std::isfinite(impulse_scale) || !std::isfinite(energy_scale) || energy_scale <= 0)
        throw std::invalid_argument("dense contact: invalid scales");
    StdVec<Real> velocity(n), lambda(rows, 0), diagonal(rows), bound(rows), wall_speed(rows);
    StdVec<StdVec<Real>> compliance(rows);
    DenseMetricContactResult result;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(position[i]) || !std::isfinite(free_velocity[i])) throw std::invalid_argument("dense contact: non-finite state");
        velocity[i] = free_velocity[i] / speed;
    }
    for (size_t s = 0; s < rows; ++s)
    {
        const auto &row = constraints[s];
        if (row.coefficients.size() != n || !std::isfinite(row.old_offset) || !std::isfinite(row.new_offset))
            throw std::invalid_argument("dense contact: invalid row");
        const auto solved = factor.solve(row.coefficients);
        result.linear_relative_residual = std::max(result.linear_relative_residual, solved.relative_residual);
        compliance[s] = solved.value;
        for (auto &value : compliance[s]) value *= factor.scale();
        diagonal[s] = denseContactDot(row.coefficients, compliance[s]);
        const Real previous = denseContactDot(row.coefficients, position);
        const Real old_gap = (previous - row.old_offset) / (dt * speed);
        if (!std::isfinite(diagonal[s]) || diagonal[s] <= 0 || !std::isfinite(old_gap) || old_gap < -settings.relative_tolerance)
            throw std::invalid_argument("dense contact: rank failure or initial penetration");
        bound[s] = (row.new_offset - previous) / (dt * speed);
        wall_speed[s] = (row.new_offset - row.old_offset) / (dt * speed);
        if (!std::isfinite(bound[s]) || !std::isfinite(wall_speed[s])) throw std::invalid_argument("dense contact: non-finite geometry");
    }
    result.converged = rows == 0;
    for (size_t sweep = 0; rows && sweep < settings.max_sweeps; ++sweep)
    {
        for (size_t s = 0; s < rows; ++s)
        {
            const Real trial = lambda[s] + (bound[s] - denseContactDot(constraints[s].coefficients, velocity)) / diagonal[s];
            if (!std::isfinite(trial)) throw std::runtime_error("dense contact: non-finite dual iterate");
            const Real next = std::max(Real(0), trial), change = next - lambda[s];
            lambda[s] = next;
            for (size_t i = 0; i < n; ++i) velocity[i] += compliance[s][i] * change;
        }
        Real primal = 0, natural = 0, complementarity = 0, work_scale = 1;
        for (size_t s = 0; s < rows; ++s)
        {
            const Real gap = denseContactDot(constraints[s].coefficients, velocity) - bound[s];
            primal = std::max(primal, -gap);
            natural = std::max(natural, std::abs(std::min(lambda[s] * diagonal[s], gap)));
            complementarity += std::abs(lambda[s] * gap);
            work_scale += std::abs(lambda[s] * wall_speed[s]);
        }
        Real stationarity = 0;
        for (size_t i = 0; i < n; ++i)
        {
            long double value = 0;
            for (size_t j = 0; j < n; ++j)
                value += metric[i][j] / factor.scale() * (velocity[j] - free_velocity[j] / speed);
            for (size_t s = 0; s < rows; ++s) value -= constraints[s].coefficients[i] * lambda[s];
            stationarity = std::max(stationarity, std::abs(Real(value)));
        }
        result.sweeps = sweep + 1;
        result.primal_velocity_error = primal * speed; result.natural_velocity_error = natural * speed;
        result.stationarity_impulse_error = stationarity * impulse_scale;
        result.complementarity_energy_error = complementarity * energy_scale;
        result.converged = primal <= settings.relative_tolerance && natural <= settings.relative_tolerance &&
            stationarity <= settings.relative_tolerance && complementarity <= settings.relative_tolerance * work_scale;
        if (result.converged) break;
    }
    result.velocity.resize(n); result.normal_impulse.resize(rows);
    for (size_t i = 0; i < n; ++i)
    {
        result.velocity[i] = velocity[i] * speed;
        if (!std::isfinite(result.velocity[i])) throw std::runtime_error("dense contact: non-finite candidate velocity");
    }
    for (size_t s = 0; s < rows; ++s) result.normal_impulse[s] = lambda[s] * impulse_scale;
    return result;
}

/** Validate convex quadratic stiffness without changing its entries. A null
 * pivot in a PSD matrix must have a null remaining row. No eigenvalue clipping. */
inline void validateDenseStiffness(const DenseContactMatrix &matrix)
{
    const size_t n = matrix.size();
    if (!n || n > 128) throw std::invalid_argument("dense stiffness: invalid dimensions");
    StdVec<StdVec<long double>> schur(n, StdVec<long double>(n));
    long double norm = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (matrix[i].size() != n) throw std::invalid_argument("dense stiffness: non-square");
        for (size_t j = 0; j < n; ++j)
        {
            if (!std::isfinite(matrix[i][j])) throw std::invalid_argument("dense stiffness: non-finite");
            schur[i][j] = matrix[i][j]; norm = std::max(norm, std::abs(schur[i][j]));
        }
    }
    const long double tolerance = 1e-14L * norm;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (std::abs(schur[i][j] - schur[j][i]) > tolerance) throw std::invalid_argument("dense stiffness: nonsymmetric");
    for (size_t i = 0; i < n; ++i)
    {
        const long double pivot = schur[i][i];
        if (pivot < -tolerance) throw std::invalid_argument("dense stiffness: non-PSD");
        if (std::abs(pivot) <= tolerance)
        {
            for (size_t j = i + 1; j < n; ++j)
                if (std::abs(schur[j][i]) > tolerance) throw std::invalid_argument("dense stiffness: non-PSD null row");
            continue;
        }
        for (size_t j = i + 1; j < n; ++j)
            for (size_t k = i + 1; k < n; ++k) schur[j][k] -= schur[j][i] * schur[k][i] / pivot;
    }
}

struct DenseElasticContactResult
{
    DenseMetricContactResult projection;
    StdVec<Real> position, internal_force, momentum_residual;
    Real kinetic_old = 0, kinetic_new = 0, elastic_old = 0, elastic_new = 0;
    Real algorithmic_loss = 0, normal_contact_loss = 0, wall_work = 0, energy_residual = 0;
};

inline DenseElasticContactResult advanceDenseElasticContact(const StdVec<Real> &mass,
    const DenseContactMatrix &stiffness, const StdVec<Real> &position, const StdVec<Real> &velocity,
    Real dt, const StdVec<DenseContactConstraint> &constraints, const NormalProjectionSettings &settings = {})
{
    const size_t n = mass.size();
    if (!n || n > 128 || stiffness.size() != n || position.size() != n || velocity.size() != n || !std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("dense elastic: invalid state dimensions/dt");
    validateDenseStiffness(stiffness);
    DenseContactMatrix metric = stiffness; StdVec<Real> rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(mass[i]) || mass[i] <= 0 || !std::isfinite(position[i]) || !std::isfinite(velocity[i]))
            throw std::invalid_argument("dense elastic: invalid physical state");
        for (size_t j = 0; j < n; ++j) metric[i][j] *= dt * dt / 2;
        metric[i][i] += mass[i];
        rhs[i] = mass[i] * velocity[i] - dt * denseContactDot(stiffness[i], position);
    }
    const auto free = DensePositiveMetric(metric).solve(rhs);
    DenseElasticContactResult r;
    r.projection = projectDenseMetricContact(metric, position, free.value, dt, constraints, settings);
    r.projection.linear_relative_residual = std::max(r.projection.linear_relative_residual, free.relative_residual);
    r.position.resize(n); r.internal_force.resize(n); r.momentum_residual.resize(n);
    StdVec<Real> average(n);
    for (size_t i = 0; i < n; ++i)
    {
        const Real v = r.projection.velocity[i], dv = v - velocity[i];
        r.position[i] = position[i] + dt * v; average[i] = (position[i] + r.position[i]) / 2;
        r.kinetic_old += mass[i] * velocity[i] * velocity[i] / 2;
        r.kinetic_new += mass[i] * v * v / 2; r.algorithmic_loss += mass[i] * dv * dv / 2;
    }
    for (size_t i = 0; i < n; ++i)
    {
        r.internal_force[i] = -denseContactDot(stiffness[i], average);
        r.momentum_residual[i] = mass[i] * (r.projection.velocity[i] - velocity[i]) - dt * r.internal_force[i];
        r.elastic_old += position[i] * denseContactDot(stiffness[i], position) / 2;
        r.elastic_new += r.position[i] * denseContactDot(stiffness[i], r.position) / 2;
    }
    for (size_t s = 0; s < constraints.size(); ++s)
    {
        const Real impulse = r.projection.normal_impulse[s]; const auto &row = constraints[s];
        r.wall_work += impulse * (row.new_offset - row.old_offset) / dt;
        r.normal_contact_loss += impulse * (denseContactDot(row.coefficients, position) - row.old_offset) / dt;
        for (size_t i = 0; i < n; ++i) r.momentum_residual[i] -= row.coefficients[i] * impulse;
    }
    r.energy_residual = r.wall_work - (r.kinetic_new - r.kinetic_old + r.elastic_new - r.elastic_old) - r.algorithmic_loss - r.normal_contact_loss;
    for (Real value : {r.kinetic_old, r.kinetic_new, r.elastic_old, r.elastic_new, r.algorithmic_loss,
                       r.normal_contact_loss, r.wall_work, r.energy_residual})
        if (!std::isfinite(value)) throw std::runtime_error("dense elastic: non-finite physical ledger");
    return r;
}
} // namespace SPH::peridynamics
#endif

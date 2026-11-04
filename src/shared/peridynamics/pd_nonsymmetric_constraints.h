/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_NONSYMMETRIC_CONSTRAINTS_H
#define PD_NONSYMMETRIC_CONSTRAINTS_H

#include "pd_dense_metric_contact.h"
#include <functional>
#include <string>

namespace SPH::peridynamics
{
/** Scaled partial-pivot Gaussian solve. The full supplied matrix is retained;
 * row permutations are algebra, not symmetrization or a physical mass shift. */
struct GeneralLinearSolution
{
    StdVec<Real> value;
    Real relative_residual = 0, minimum_relative_pivot = 1;
};

inline GeneralLinearSolution solveGeneralContactSystem(const DenseContactMatrix &matrix, const StdVec<Real> &rhs)
{
    const size_t n = rhs.size();
    if (!n || n > 132 || matrix.size() != n) throw std::invalid_argument("general contact linear solve: dimensions");
    StdVec<StdVec<long double>> a(n, StdVec<long double>(n + 1));
    long double norm = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (matrix[i].size() != n || !std::isfinite(rhs[i])) throw std::invalid_argument("general contact linear solve: invalid row");
        long double row = 0;
        for (size_t j = 0; j < n; ++j)
        {
            if (!std::isfinite(matrix[i][j])) throw std::invalid_argument("general contact linear solve: nonfinite");
            a[i][j] = matrix[i][j]; row += std::abs(a[i][j]);
        }
        a[i][n] = rhs[i]; norm = std::max(norm, row);
    }
    if (!(norm > 0)) throw std::runtime_error("general contact linear solve: singular matrix");
    GeneralLinearSolution result;
    for (size_t k = 0; k < n; ++k)
    {
        size_t pivot = k;
        for (size_t i = k + 1; i < n; ++i) if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
        const long double relative = std::abs(a[pivot][k]) / norm;
        if (relative <= 1e-14L) throw std::runtime_error("general contact linear solve: pivot/conditioning guard");
        result.minimum_relative_pivot = std::min(result.minimum_relative_pivot, Real(relative));
        std::swap(a[k], a[pivot]);
        for (size_t i = k + 1; i < n; ++i)
        {
            const long double multiplier = a[i][k] / a[k][k]; a[i][k] = 0;
            for (size_t j = k + 1; j <= n; ++j) a[i][j] -= multiplier * a[k][j];
        }
    }
    StdVec<long double> value(n, 0); result.value.resize(n);
    for (size_t i = n; i-- > 0;)
    {
        value[i] = a[i][n];
        for (size_t j = i + 1; j < n; ++j) value[i] -= a[i][j] * value[j];
        value[i] /= a[i][i]; result.value[i] = value[i];
        if (!std::isfinite(result.value[i])) throw std::runtime_error("general contact linear solve: nonfinite solution");
    }
    long double error = 0, max_x = 0, max_b = 0;
    for (size_t i = 0; i < n; ++i)
    {
        max_x = std::max(max_x, std::abs(static_cast<long double>(result.value[i])));
        max_b = std::max(max_b, std::abs(static_cast<long double>(rhs[i])));
        long double residual = -rhs[i];
        for (size_t j = 0; j < n; ++j) residual += static_cast<long double>(matrix[i][j]) * result.value[j];
        error = std::max(error, std::abs(residual));
    }
    const long double denominator = norm * max_x + max_b;
    result.relative_residual = denominator > 0 ? error / denominator : 0;
    if (result.relative_residual > 1e-12) throw std::runtime_error("general contact linear solve: residual gate");
    return result;
}

struct NonlinearContactEquation
{
    StdVec<Real> momentum;
    DenseContactMatrix jacobian;
};

struct PrimalDualSettings
{
    Real velocity_scale = 1, relative_tolerance = 1e-12;
    size_t max_newton_steps = 32, max_line_search = 24;
};

struct PrimalDualTrace
{
    size_t active_mask, iteration, domain_rejections;
    Real merit, accepted_step, linear_residual;
};

struct PrimalDualResult
{
    bool converged = false;
    std::string status = "no_KKT_feasible_active_set";
    StdVec<Real> velocity, normal_impulse;
    StdVec<PrimalDualTrace> trace;
    size_t active_mask = 0, dependent_sets_skipped = 0;
    Real normalized_KKT = 0, maximum_linear_residual = 0;
};

inline bool independentContactRows(const StdVec<DenseContactConstraint> &rows, const StdVec<size_t> &active)
{
    StdVec<StdVec<long double>> orthogonal;
    for (size_t index : active)
    {
        const auto &row = rows[index].coefficients;
        StdVec<long double> residual(row.begin(), row.end()); long double before = 0;
        for (auto value : residual) before += value * value;
        for (const auto &basis : orthogonal)
        {
            long double projection = 0;
            for (size_t i = 0; i < row.size(); ++i) projection += residual[i] * basis[i];
            for (size_t i = 0; i < row.size(); ++i) residual[i] -= projection * basis[i];
        }
        long double after = 0; for (auto value : residual) after += value * value;
        if (after <= 1e-24L * before) return false;
        for (auto &value : residual) value /= std::sqrt(after);
        orthogonal.push_back(residual);
    }
    return true;
}

/** Bounded active-set equality-root search. Only algebraically independent
 * sets are solved; all original rows remain in the final KKT check. A singular
 * global solve, failed line search or Newton limit aborts the whole call.
 * Trial vectors are cold iterates, never committed physical resets/history.
 */
inline PrimalDualResult solveNonsymmetricContact(const StdVec<Real> &old_position, Real dt, Real physical_mass_scale,
    const StdVec<DenseContactConstraint> &rows,
    const std::function<NonlinearContactEquation(const StdVec<Real> &)> &evaluate,
    const std::function<bool(const StdVec<Real> &)> &in_domain,
    const PrimalDualSettings &settings = {})
{
    const size_t n = old_position.size(), count = rows.size();
    const Real speed = settings.velocity_scale, impulse_scale = physical_mass_scale * speed;
    if (!n || n > 128 || count > 4 || !std::isfinite(dt) || dt <= 0 || !std::isfinite(physical_mass_scale) || physical_mass_scale <= 0 ||
        !std::isfinite(speed) || speed <= 0 || !std::isfinite(impulse_scale) || !std::isfinite(settings.relative_tolerance) ||
        settings.relative_tolerance <= 0 || settings.relative_tolerance > 1e-3 || !settings.max_newton_steps || settings.max_newton_steps > 64 ||
        !settings.max_line_search || settings.max_line_search > 32)
        throw std::invalid_argument("nonsymmetric contact: bounded dimensions/settings");
    StdVec<Real> bounds(count);
    for (Real value : old_position) if (!std::isfinite(value)) throw std::invalid_argument("nonsymmetric contact: nonfinite old position");
    for (size_t s = 0; s < count; ++s)
    {
        const auto &row = rows[s];
        if (row.coefficients.size() != n || !std::isfinite(row.old_offset) || !std::isfinite(row.new_offset))
            throw std::invalid_argument("nonsymmetric contact: invalid row");
        Real norm = 0;
        for (Real value : row.coefficients) { if (!std::isfinite(value)) throw std::invalid_argument("nonsymmetric contact: nonfinite row"); norm += value * value; }
        const Real previous = denseContactDot(row.coefficients, old_position);
        if (!(norm > 0) || !std::isfinite(norm) || !std::isfinite(previous) || (previous - row.old_offset) / (dt * speed) < -settings.relative_tolerance)
            throw std::invalid_argument("nonsymmetric contact: zero row/initial penetration");
        bounds[s] = (row.new_offset - previous) / dt;
        if (!std::isfinite(bounds[s])) throw std::invalid_argument("nonsymmetric contact: nonfinite bound");
    }
    PrimalDualResult result;
    for (size_t cardinality = 0; cardinality <= count; ++cardinality)
    for (size_t mask = 0; mask < (size_t(1) << count); ++mask)
    {
        StdVec<size_t> active;
        for (size_t s = 0; s < count; ++s) if (mask & (size_t(1) << s)) active.push_back(s);
        if (active.size() != cardinality) continue;
        if (!independentContactRows(rows, active)) { ++result.dependent_sets_skipped; continue; }
        StdVec<Real> velocity(n, 0), lambda(active.size(), 0);
        auto residual = [&](const NonlinearContactEquation &equation, const StdVec<Real> &v, const StdVec<Real> &l)
        {
            if (equation.momentum.size() != n || equation.jacobian.size() != n) throw std::runtime_error("nonsymmetric contact: callback dimensions");
            StdVec<Real> values(n + active.size());
            for (size_t i = 0; i < n; ++i)
            {
                if (equation.jacobian[i].size() != n || !std::isfinite(equation.momentum[i])) throw std::runtime_error("nonsymmetric contact: callback row");
                for (Real value : equation.jacobian[i]) if (!std::isfinite(value)) throw std::runtime_error("nonsymmetric contact: nonfinite Jacobian");
                values[i] = equation.momentum[i];
                for (size_t j = 0; j < active.size(); ++j) values[i] -= rows[active[j]].coefficients[i] * l[j];
                values[i] /= impulse_scale;
            }
            for (size_t j = 0; j < active.size(); ++j)
                values[n + j] = (denseContactDot(rows[active[j]].coefficients, v) - bounds[active[j]]) / speed;
            return values;
        };
        auto merit = [](const StdVec<Real> &values)
        {
            Real maximum = 0;
            for (Real value : values) { if (!std::isfinite(value)) throw std::runtime_error("nonsymmetric contact: nonfinite residual"); maximum = std::max(maximum, std::abs(value)); }
            return maximum;
        };
        bool root = false; Real root_merit = 0;
        for (size_t iteration = 0; iteration <= settings.max_newton_steps; ++iteration)
        {
            if (!in_domain(velocity)) throw std::runtime_error("nonsymmetric contact: invalid initial iterate domain");
            const auto equation = evaluate(velocity);
            const auto values = residual(equation, velocity, lambda); const Real current = merit(values);
            if (current <= settings.relative_tolerance)
            {
                result.trace.push_back({mask, iteration, 0, current, 0, 0}); root = true; root_merit = current; break;
            }
            if (iteration == settings.max_newton_steps)
            {
                // Retain the terminal residual even though no physical state
                // is committed. Frozen v11 controls predate this diagnostic fix.
                result.trace.push_back({mask, iteration, 0, current, 0, 0});
                result.status = "Newton_iteration_limit"; return result;
            }
            const size_t size = n + active.size();
            DenseContactMatrix system(size, StdVec<Real>(size, 0)); StdVec<Real> rhs(size);
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < n; ++j) system[i][j] = equation.jacobian[i][j] * speed / impulse_scale;
                for (size_t j = 0; j < active.size(); ++j) system[i][n + j] = -rows[active[j]].coefficients[i];
            }
            for (size_t j = 0; j < active.size(); ++j)
                for (size_t i = 0; i < n; ++i) system[n + j][i] = rows[active[j]].coefficients[i];
            for (size_t i = 0; i < size; ++i) rhs[i] = -values[i];
            const auto direction = solveGeneralContactSystem(system, rhs);
            result.maximum_linear_residual = std::max(result.maximum_linear_residual, direction.relative_residual);
            bool accepted = false; Real alpha = 1; size_t domain_rejections = 0;
            for (size_t search = 0; search < settings.max_line_search; ++search, alpha *= .5)
            {
                auto v = velocity, l = lambda;
                for (size_t i = 0; i < n; ++i) v[i] += alpha * speed * direction.value[i];
                for (size_t j = 0; j < active.size(); ++j) l[j] += alpha * impulse_scale * direction.value[n + j];
                if (!in_domain(v)) { ++domain_rejections; continue; }
                const Real next = merit(residual(evaluate(v), v, l));
                if (next <= (1 - 1e-4 * alpha) * current)
                {
                    velocity = std::move(v); lambda = std::move(l); accepted = true; break;
                }
            }
            result.trace.push_back({mask, iteration, domain_rejections, current, accepted ? alpha : 0, direction.relative_residual});
            if (!accepted) { result.status = "line_search_limit"; return result; }
        }
        if (!root) throw std::runtime_error("nonsymmetric contact: missing equality root state");
        StdVec<Real> impulses(count, 0);
        bool dual = true;
        for (size_t j = 0; j < active.size(); ++j) { impulses[active[j]] = lambda[j]; if (lambda[j] < 0) dual = false; }
        Real primal = 0, complementarity = 0;
        for (size_t s = 0; s < count; ++s)
        {
            const Real gap = denseContactDot(rows[s].coefficients, velocity) - bounds[s];
            primal = std::max(primal, -gap / speed); complementarity += std::abs(impulses[s] * gap) / (impulse_scale * speed);
        }
        if (!dual || primal > settings.relative_tolerance || complementarity > settings.relative_tolerance) continue;
        result.converged = true; result.status = "KKT_root_accepted"; result.active_mask = mask;
        result.velocity = velocity; result.normal_impulse = impulses;
        result.normalized_KKT = std::max(root_merit, std::max(primal, complementarity)); return result;
    }
    return result;
}
} // namespace SPH::peridynamics
#endif

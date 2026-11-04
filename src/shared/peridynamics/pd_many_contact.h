/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_MANY_CONTACT_H
#define PD_MANY_CONTACT_H

#include "pd_nonsymmetric_constraints.h"
#include <Eigen/SVD>

namespace SPH::peridynamics
{
/** A single algebraic method for full-rank and redundant Newton systems.
 * Singular values below the declared cutoff are omitted, but the ORIGINAL
 * unmodified matrix must still satisfy the residual gate. This is not a
 * regularized physical operator, mass change or runtime fallback. */
struct ConsistentLinearSolution
{
    StdVec<Real> value;
    size_t rank = 0;
    Real relative_residual = 0;
};

inline ConsistentLinearSolution solveConsistentContactSystem(const DenseContactMatrix &matrix, const StdVec<Real> &rhs)
{
    const size_t n = rhs.size();
    if (!n || n > 384 || matrix.size() != n) throw std::invalid_argument("consistent contact solve: dimensions");
    Eigen::MatrixXd a(n, n); Eigen::VectorXd b(n);
    Real norm = 0, norm_b = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (matrix[i].size() != n || !std::isfinite(rhs[i])) throw std::invalid_argument("consistent contact solve: row");
        b[i] = rhs[i]; norm_b = std::max(norm_b, std::abs(rhs[i])); Real row_norm = 0;
        for (size_t j = 0; j < n; ++j)
        {
            if (!std::isfinite(matrix[i][j])) throw std::invalid_argument("consistent contact solve: nonfinite");
            a(i,j) = matrix[i][j]; row_norm += std::abs(matrix[i][j]);
        }
        norm = std::max(norm, row_norm);
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1e-13);
    const Eigen::VectorXd x = svd.solve(b);
    ConsistentLinearSolution result; result.rank = svd.rank(); result.value.resize(n);
    long double norm_x = 0, error = 0;
    for (size_t j = 0; j < n; ++j)
    {
        if (!std::isfinite(x[j])) throw std::runtime_error("consistent contact solve: nonfinite solution");
        result.value[j] = x[j]; norm_x = std::max(norm_x, std::abs(static_cast<long double>(x[j])));
    }
    for (size_t i = 0; i < n; ++i)
    {
        long double residual = -rhs[i];
        for (size_t j = 0; j < n; ++j) residual += static_cast<long double>(matrix[i][j]) * result.value[j];
        error = std::max(error, std::abs(residual));
    }
    const long double scale = norm * norm_x + norm_b;
    result.relative_residual = scale > 0 ? error / scale : 0;
    if (result.relative_residual > 1e-12) throw std::runtime_error("consistent contact solve: inconsistent or inaccurate original linear system");
    return result;
}

struct ComplementarityValue
{
    Real value, impulse_derivative, gap_derivative;
};

inline ComplementarityValue contactFischerBurmeister(Real impulse, Real gap)
{
    if (!std::isfinite(impulse) || !std::isfinite(gap)) throw std::invalid_argument("contact complementarity: nonfinite");
    const Real radius = std::hypot(impulse, gap), sum = impulse + gap;
    // The rational form retains small nonzero products in the feasible quadrant.
    const Real phi = sum > 0 ? -2 * impulse * (gap / (radius + sum)) : radius - sum;
    if (radius == 0)
    {
        const Real derivative = std::sqrt(Real(.5)) - 1;
        return {0, derivative, derivative};
    }
    return {phi, impulse / radius - 1, gap / radius - 1};
}

struct ManyContactSettings
{
    Real velocity_scale = 1, relative_tolerance = 1e-12;
    size_t max_newton_steps = 48, max_line_search = 32;
};

struct ManyContactTrace
{
    size_t iteration = 0, linear_rank = 0, domain_rejections = 0;
    Real merit = 0, squared_norm = 0, KKT = 0, accepted_step = 0, linear_residual = 0;
    StdVec<Real> velocity, normal_impulse, normalized_direction;
    StdVec<Real> normalized_residual, normalized_gaps;
    DenseContactMatrix generalized_jacobian;
};

struct ManyContactResult
{
    bool converged = false;
    std::string status;
    // Empty unless the original momentum and every original KKT row pass.
    StdVec<Real> velocity, normal_impulse;
    StdVec<ManyContactTrace> trace;
    Real normalized_KKT = 0;
};

/** Bounded semismooth Newton control. Every original contact row is retained;
 * there is no active-set subset enumeration, matrix shift or line-search
 * gradient fallback. Failure returns no accepted physical state. Callbacks
 * must read immutable accepted material/history; only the caller may commit.
 */
inline ManyContactResult solveManyNonsymmetricContact(const StdVec<Real> &old_position, Real dt, Real mass_scale,
    const StdVec<DenseContactConstraint> &rows,
    const std::function<NonlinearContactEquation(const StdVec<Real> &)> &evaluate,
    const std::function<bool(const StdVec<Real> &)> &in_domain,
    const ManyContactSettings &settings = {})
{
    const size_t n = old_position.size(), count = rows.size(), size = n + count;
    const Real speed = settings.velocity_scale, impulse_scale = mass_scale * speed;
    if (!n || n > 128 || count > 256 || !std::isfinite(dt) || dt <= 0 || !std::isfinite(mass_scale) || mass_scale <= 0 ||
        !std::isfinite(speed) || speed <= 0 || !std::isfinite(impulse_scale) || !std::isfinite(settings.relative_tolerance) ||
        settings.relative_tolerance <= 0 || settings.relative_tolerance > 1e-3 || !settings.max_newton_steps || settings.max_newton_steps > 48 ||
        !settings.max_line_search || settings.max_line_search > 32)
        throw std::invalid_argument("many contact: bounded dimensions/settings");
    for (Real x : old_position) if (!std::isfinite(x)) throw std::invalid_argument("many contact: nonfinite position");
    StdVec<Real> bounds(count);
    for (size_t s = 0; s < count; ++s)
    {
        const auto &row = rows[s];
        if (row.coefficients.size() != n || !std::isfinite(row.old_offset) || !std::isfinite(row.new_offset)) throw std::invalid_argument("many contact: row");
        Real norm = 0;
        for (Real c : row.coefficients) { if (!std::isfinite(c)) throw std::invalid_argument("many contact: nonfinite row"); norm += c*c; }
        const Real old_gather = denseContactDot(row.coefficients, old_position);
        if (!(norm > 0) || !std::isfinite(norm) || !std::isfinite(old_gather) ||
            (old_gather - row.old_offset) / (dt * speed) < -settings.relative_tolerance)
            throw std::invalid_argument("many contact: zero row/initial penetration");
        bounds[s] = (row.new_offset - old_gather) / dt;
        if (!std::isfinite(bounds[s])) throw std::invalid_argument("many contact: nonfinite bound");
    }
    struct Values { StdVec<Real> residual, gaps; DenseContactMatrix jacobian; Real merit = 0, squared_norm = 0, KKT = 0; };
    auto values = [&](const StdVec<Real> &v, const StdVec<Real> &lambda, bool tangent)
    {
        const auto equation = evaluate(v); Values out; out.residual.resize(size); out.gaps.resize(count);
        if (equation.momentum.size() != n || equation.jacobian.size() != n) throw std::runtime_error("many contact: callback dimensions");
        if (tangent) out.jacobian.assign(size, StdVec<Real>(size, 0));
        for (size_t i = 0; i < n; ++i)
        {
            if (equation.jacobian[i].size() != n || !std::isfinite(equation.momentum[i])) throw std::runtime_error("many contact: callback row");
            Real r = equation.momentum[i];
            for (size_t s = 0; s < count; ++s) r -= rows[s].coefficients[i] * lambda[s];
            out.residual[i] = r / impulse_scale;
            for (size_t j = 0; j < n; ++j)
            {
                if (!std::isfinite(equation.jacobian[i][j])) throw std::runtime_error("many contact: nonfinite tangent");
                if (tangent) out.jacobian[i][j] = equation.jacobian[i][j] * speed / impulse_scale;
            }
            if (tangent) for (size_t s = 0; s < count; ++s) out.jacobian[i][n+s] = -rows[s].coefficients[i];
        }
        Real product_sum = 0;
        for (size_t s = 0; s < count; ++s)
        {
            const Real a = lambda[s] / impulse_scale, b = (denseContactDot(rows[s].coefficients, v) - bounds[s]) / speed;
            out.gaps[s] = b;
            const auto fb = contactFischerBurmeister(a, b); out.residual[n+s] = fb.value;
            out.KKT = std::max(out.KKT, std::max(-a, -b)); product_sum += std::abs(a*b);
            if (tangent)
            {
                for (size_t i = 0; i < n; ++i) out.jacobian[n+s][i] = fb.gap_derivative * rows[s].coefficients[i];
                out.jacobian[n+s][n+s] = fb.impulse_derivative;
            }
        }
        out.KKT = std::max(out.KKT, product_sum);
        for (Real r : out.residual)
        {
            if (!std::isfinite(r)) throw std::runtime_error("many contact: nonfinite residual");
            out.merit = std::max(out.merit, std::abs(r)); out.squared_norm += r*r;
        }
        if (!std::isfinite(out.squared_norm) || !std::isfinite(out.KKT)) throw std::runtime_error("many contact: nonfinite merit/KKT");
        return out;
    };
    ManyContactResult result; StdVec<Real> velocity(n, 0), lambda(count, 0);
    if (!in_domain(velocity)) throw std::invalid_argument("many contact: invalid initial iterate domain");
    for (size_t iteration = 0; iteration <= settings.max_newton_steps; ++iteration)
    {
        const auto current = values(velocity, lambda, true);
        ManyContactTrace trace; trace.iteration = iteration; trace.merit = current.merit; trace.squared_norm = current.squared_norm;
        trace.KKT = current.KKT; trace.velocity = velocity; trace.normal_impulse = lambda;
        trace.normalized_residual = current.residual; trace.normalized_gaps = current.gaps;
        trace.generalized_jacobian = current.jacobian;
        if (std::max(current.merit, current.KKT) <= settings.relative_tolerance)
        {
            result.trace.push_back(std::move(trace)); result.converged = true; result.status = "all_rows_KKT_root_accepted";
            result.velocity = velocity; result.normal_impulse = lambda; result.normalized_KKT = std::max(current.merit, current.KKT); return result;
        }
        if (iteration == settings.max_newton_steps)
        { result.trace.push_back(std::move(trace)); result.status = "Newton_iteration_limit"; return result; }
        StdVec<Real> rhs(size); for (size_t i = 0; i < size; ++i) rhs[i] = -current.residual[i];
        ConsistentLinearSolution direction;
        try { direction = solveConsistentContactSystem(current.jacobian, rhs); }
        catch (const std::runtime_error &)
        { result.trace.push_back(std::move(trace)); result.status = "original_linear_system_rejected"; return result; }
        trace.linear_rank = direction.rank; trace.linear_residual = direction.relative_residual; trace.normalized_direction = direction.value;
        Real alpha = 1; bool accepted = false;
        for (size_t search = 0; search < settings.max_line_search; ++search, alpha *= .5)
        {
            auto v = velocity, l = lambda;
            for (size_t i = 0; i < n; ++i) v[i] += alpha * speed * direction.value[i];
            for (size_t s = 0; s < count; ++s) l[s] += alpha * impulse_scale * direction.value[n+s];
            if (!in_domain(v)) { ++trace.domain_rejections; continue; }
            const auto trial = values(v, l, false);
            if (trial.squared_norm <= (1 - 1e-4 * alpha) * current.squared_norm)
            { trace.accepted_step = alpha; velocity = std::move(v); lambda = std::move(l); accepted = true; break; }
        }
        result.trace.push_back(std::move(trace));
        if (!accepted) { result.status = "line_search_limit"; return result; }
    }
    throw std::logic_error("many contact: unreachable");
}
} // namespace SPH::peridynamics
#endif

/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_GEL_AFFINE_PATH_H
#define PD_GEL_AFFINE_PATH_H

#include "pd_gel_material.h"
#include <array>

namespace SPH::peridynamics
{
inline Matd gelPathDeviator(const Matd &value)
{
    return value - value.trace() / 3 * Matd::Identity();
}

/** Coefficients of X(s)=exp(-hs)X0+a1(s)G1+a2(s)G2 for a quadratic
 * Green path E(s)=E0+sG1+s²G2. This is not a linear-Green approximation. */
inline std::array<long double, 3> gelPathCoefficients(long double h, long double s)
{
    const long double z = h * s;
    if (z < 0.25L)
    {
        long double a0 = 1, a1 = 1, a2 = 1, term = 1;
        for (int k = 1; k < 32; ++k)
        {
            term *= -z / k; a0 += term;
            a1 += term / (k + 1); a2 += 2 * term / ((k + 1) * (k + 2));
        }
        return {a0, s * a1, s * s * a2};
    }
    const long double decay = std::exp(-z), a1 = -std::expm1(-z) / h;
    return {decay, a1, 2 * (s - a1) / h};
}

/** Exact analytic Gram moments int a_i(s)*a_j(s) ds, with positive integral
 * dissipation. Series avoid cancellation for small h; no clipping is applied. */
inline std::array<std::array<long double, 3>, 3> gelPathGram(long double h)
{
    std::array<std::array<long double, 3>, 3> gram{};
    if (h < 0.25L)
    {
        std::array<long double, 34> factorial{}; factorial[0] = 1;
        for (size_t i = 1; i < factorial.size(); ++i) factorial[i] = factorial[i - 1] * i;
        std::array<std::array<long double, 32>, 3> coefficient{};
        long double power = 1;
        for (size_t k = 0; k < 32; ++k)
        {
            for (size_t i = 0; i < 3; ++i) coefficient[i][k] = (i == 2 ? 2 : 1) * power / factorial[k + i];
            power *= -h;
        }
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                for (size_t k = 0; k < 32; ++k)
                    for (size_t l = 0; l < 32; ++l)
                        gram[i][j] += coefficient[i][k] * coefficient[j][l] / (k + l + i + j + 1);
    }
    else
    {
        const long double decay = std::exp(-h), e0 = -std::expm1(-h) / h;
        const long double e00 = -std::expm1(-2 * h) / (2 * h);
        const long double e1 = (1 - (1 + h) * decay) / (h * h), square = 1 - 2 * e0 + e00;
        gram[0][0] = e00; gram[0][1] = gram[1][0] = (e0 - e00) / h;
        gram[0][2] = gram[2][0] = 2 * e1 / h - 2 * (e0 - e00) / (h * h);
        gram[1][1] = square / (h * h);
        gram[1][2] = gram[2][1] = (1 - 2 * e1) / (h * h) - 2 * square / (h * h * h);
        gram[2][2] = 4 / (h * h) * (1.L / 3 - (1 - 2 * e1) / h + square / (h * h));
    }
    return gram;
}

inline const StdVec<std::pair<Real, Real>> &gelPathGauss32()
{
    static const auto rule = []
    {
        StdVec<std::pair<Real, Real>> values;
        for (int i = 0; i < 16; ++i)
        {
            long double x = std::cos(std::acos(-1.L) * (i + 0.75L) / 32.5L), derivative = 0;
            bool converged = false;
            for (int iteration = 0; iteration < 32; ++iteration)
            {
                long double previous = 1, current = x;
                for (int degree = 2; degree <= 32; ++degree)
                {
                    const long double next = ((2 * degree - 1) * x * current - (degree - 1) * previous) / degree;
                    previous = current; current = next;
                }
                derivative = 32 * (x * current - previous) / (x * x - 1);
                const long double next = x - current / derivative;
                if (std::abs(next - x) <= 2e-19L) { x = next; converged = true; break; }
                x = next;
            }
            if (!converged) throw std::runtime_error("gel affine path: fixed quadrature construction failed");
            const long double weight = 1 / ((1 - x * x) * derivative * derivative);
            values.emplace_back((1 - x) / 2, weight); values.emplace_back((1 + x) / 2, weight);
        }
        return values;
    }();
    return rule;
}

inline Real gelAffineMinimumDeterminant(const Matd &initial, const Matd &increment)
{
    std::array<long double, 4> coefficients{};
    for (unsigned mask = 0; mask < 8; ++mask)
    {
        Matd matrix = initial; size_t order = 0;
        for (size_t column = 0; column < 3; ++column)
            if (mask & (1u << column)) { matrix.col(column) = increment.col(column); ++order; }
        coefficients[order] += matrix.determinant();
    }
    Real minimum = std::min(initial.determinant(), (initial + increment).determinant());
    auto consider = [&](long double s)
    {
        if (s > 0 && s < 1) minimum = std::min(minimum, (initial + Real(s) * increment).determinant());
    };
    const long double a = 3 * coefficients[3], b = 2 * coefficients[2], c = coefficients[1];
    if (a == 0) { if (b != 0) consider(-c / b); }
    else
    {
        const long double discriminant = b * b - 4 * a * c;
        if (discriminant >= 0)
        {
            consider((-b + std::sqrt(discriminant)) / (2 * a));
            consider((-b - std::sqrt(discriminant)) / (2 * a));
        }
    }
    return minimum;
}

struct GelAffinePathResult
{
    GelHistory history;
    Matd average_pk1 = Matd::Zero();
    std::array<Matd, 9> endpoint_tangent; // Column k maps row-major delta F1[k] to delta average P.
    Real stored_old = 0, stored_new = 0, physical_dissipation = 0, work = 0, energy_residual = 0;
    Real minimum_determinant = 0, tangent_asymmetry = 0;
    bool tangent_available = true;
};

/** Path-consistent material research candidate. Immutable accepted history;
 * no correction to force/work, no parameter fit and no production integration.
 * The consistent tangent may be nonsymmetric: do not feed it to an SPD solver
 * by symmetrizing or replacing its diagonal without a separate derivation. */
inline GelAffinePathResult advanceGelAffinePath(const Matd &initial, const Matd &final, Real dt,
    const GelParameters &p, const GelHistory &old)
{
    p.validate();
    if (!initial.allFinite() || !final.allFinite() || !std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("gel affine path: invalid interval/state");
    const auto initial_response = gelResponse(initial, p, old);
    const Matd e0 = gelDeviatoricGreen(initial), delta = final - initial;
    if ((old.green - e0).norm() > 1e-12 * (1 + e0.norm()))
        throw std::invalid_argument("gel affine path: inconsistent initial history");
    GelAffinePathResult r; for (auto &column : r.endpoint_tangent) column.setZero();
    r.stored_old = initial_response.free_energy;
    r.minimum_determinant = gelAffineMinimumDeterminant(initial, delta);
    if (!std::isfinite(r.minimum_determinant) || r.minimum_determinant < 0.05)
        throw std::runtime_error("gel affine path: determinant below 0.05 within interval");
    if (dt == 0)
    {
        if (delta.norm() != 0) throw std::invalid_argument("gel affine path: nonzero path at zero dt");
        r.history = old; r.stored_new = r.stored_old; r.average_pk1 = initial_response.pk1;
        r.tangent_available = false; return r;
    }
    const long double h = static_cast<long double>(dt) / p.relaxation_time;
    const Matd g1 = gelPathDeviator(0.5 * (initial.transpose() * delta + delta.transpose() * initial));
    const Matd g2 = gelPathDeviator(0.5 * delta.transpose() * delta), x0 = e0 - old.viscous;
    const std::array<Matd, 3> history_basis{x0, g1, g2};
    const auto end = gelPathCoefficients(h, 1);
    const auto gram = gelPathGram(h);
    const Matd x1 = Real(end[0]) * x0 + Real(end[1]) * g1 + Real(end[2]) * g2;
    long double integral = 0;
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t row = 0; row < 3; ++row)
                for (size_t column = 0; column < 3; ++column)
                    integral += gram[i][j] * static_cast<long double>(history_basis[i](row, column)) * history_basis[j](row, column);
    r.physical_dissipation = 2 * p.branch_shear * h * integral;
    if (!std::isfinite(r.physical_dissipation) || r.physical_dissipation < 0)
        throw std::runtime_error("gel affine path: invalid physical dissipation, no clipping");
    r.history.green = gelDeviatoricGreen(final); r.history.viscous = r.history.green - x1;
    r.history.dissipation = old.dissipation + r.physical_dissipation;
    r.stored_new = gelResponse(final, p, r.history).free_energy;
    std::array<long double, 9> stress_sum{};
    std::array<std::array<long double, 9>, 9> tangent_sum{};
    for (const auto &[s, weight] : gelPathGauss32())
    {
        const Matd f = initial + s * delta, inverse_transpose = f.inverse().transpose();
        const Real determinant = f.determinant(), trace = f.squaredNorm(), isochoric = std::pow(determinant, -2.L / 3);
        const auto coefficient = gelPathCoefficients(h, s);
        const Matd x = Real(coefficient[0]) * x0 + Real(coefficient[1]) * g1 + Real(coefficient[2]) * g2;
        const Matd deviator = f - trace / 3 * inverse_transpose;
        const Matd pk1 = p.equilibrium_shear * isochoric * deviator +
            0.5 * p.bulk * (determinant * determinant - 1) * inverse_transpose + 2 * p.branch_shear * f * x;
        for (size_t a = 0; a < 9; ++a) stress_sum[a] += static_cast<long double>(weight) * pk1(a / 3, a % 3);
        for (size_t b = 0; b < 9; ++b)
        {
            Matd direction = Matd::Zero(); direction(b / 3, b % 3) = 1;
            const Matd df = s * direction;
            const Real jacobian_rate = inverse_transpose.cwiseProduct(df).sum();
            const Matd dinverse = -inverse_transpose * df.transpose() * inverse_transpose;
            const Real dtrace = 2 * f.cwiseProduct(df).sum();
            const Matd dg1 = gelPathDeviator(0.5 * (initial.transpose() * direction + direction.transpose() * initial));
            const Matd dg2 = gelPathDeviator(0.5 * (delta.transpose() * direction + direction.transpose() * delta));
            const Matd dx = Real(coefficient[1]) * dg1 + Real(coefficient[2]) * dg2;
            const Matd dp = p.equilibrium_shear * isochoric * (df - dtrace / 3 * inverse_transpose - trace / 3 * dinverse -
                2 * jacobian_rate / 3 * deviator) + 0.5 * p.bulk * (2 * determinant * determinant * jacobian_rate * inverse_transpose +
                (determinant * determinant - 1) * dinverse) + 2 * p.branch_shear * (df * x + f * dx);
            for (size_t a = 0; a < 9; ++a) tangent_sum[b][a] += static_cast<long double>(weight) * dp(a / 3, a % 3);
        }
    }
    long double work = 0, asymmetry = 0, norm = 0;
    for (size_t a = 0; a < 9; ++a)
    {
        r.average_pk1(a / 3, a % 3) = stress_sum[a];
        work += stress_sum[a] * delta(a / 3, a % 3);
        for (size_t b = 0; b < 9; ++b)
        {
            r.endpoint_tangent[b](a / 3, a % 3) = tangent_sum[b][a];
            norm += tangent_sum[b][a] * tangent_sum[b][a];
            const long double difference = tangent_sum[b][a] - tangent_sum[a][b];
            asymmetry += difference * difference;
        }
    }
    r.work = work; r.energy_residual = r.work - (r.stored_new - r.stored_old + r.physical_dissipation);
    r.tangent_asymmetry = norm > 0 ? std::sqrt(asymmetry / norm) : 0;
    if (!r.average_pk1.allFinite() || !std::isfinite(r.energy_residual) || !std::isfinite(r.tangent_asymmetry))
        throw std::runtime_error("gel affine path: non-finite work/stress/tangent");
    return r;
}
} // namespace SPH::peridynamics
#endif

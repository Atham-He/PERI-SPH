/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_COULOMB_DISK_H
#define PD_COULOMB_DISK_H

#include "base_data_type.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SPH::peridynamics
{
using ContactTangentVector = Eigen::Matrix<Real, 2, 1>;
using ContactTangentMatrix = Eigen::Matrix<Real, 2, 2>;

enum class CoulombDiskBranch { StickInterior, SlipExterior, ZeroRadiusOrigin };

/** Dimensionless disk natural map t-Proj_{mu*max(l,0)}(t-w).
 * l is normal impulse/I, t is the two-component tangent impulse/I, and w is
 * relative tangent velocity/speed. This is NOT an associated full-cone law.
 * The normal Signorini complementarity equation remains a separate equation.
 */
struct CoulombDiskValue
{
    ContactTangentVector residual = ContactTangentVector::Zero();
    ContactTangentVector normal_derivative = ContactTangentVector::Zero();
    ContactTangentMatrix impulse_derivative = ContactTangentMatrix::Zero();
    ContactTangentMatrix velocity_derivative = ContactTangentMatrix::Zero();
    ContactTangentVector residual_roundoff_bound = ContactTangentVector::Zero();
    Real maximum_derivative_roundoff_bound = 0;
    Real friction_loss = 0, friction_loss_roundoff_bound = 0;
    Real cone_excess_upper = 0, maximum_dissipation_error_upper = 0;
    bool representable = true, within_budget = true;
    CoulombDiskBranch branch = CoulombDiskBranch::ZeroRadiusOrigin;
};

namespace coulomb_disk_detail
{
using Wide = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<
    256, boost::multiprecision::backends::digit_base_2>>;
static_assert(std::numeric_limits<Wide>::digits == 256);
static_assert(std::numeric_limits<Real>::is_iec559 && std::numeric_limits<Real>::digits == 53);

inline Wide absolute(const Wide &x) { return boost::multiprecision::abs(x); }

/** Each expression has fewer than 64 elementary operations. Input binary64
 * values are exact in Wide. Squared norms contain no cancellation; close
 * binary64 subtraction is exact in Wide, while distant subtraction has
 * relative error below its Wide rounding unit. The 2^-240 envelope leaves
 * 16 guard bits over 256-bit arithmetic for projection/norm/derivative/work
 * expressions and branch-boundary continuity. Binary64 conversion error is
 * added explicitly. This conditional budget must still be independently
 * checked over the frozen input family before physical qualification.
 */
inline Wide arithmeticEnvelope(const Wide &magnitude)
{
    return boost::multiprecision::ldexp(Wide(1), -240) * (1 + magnitude);
}

inline Real outwardPositive(const Wide &x, CoulombDiskValue &result)
{
    if (x < 0) throw std::logic_error("Coulomb disk: negative error enclosure");
    if (x >= Wide(std::numeric_limits<Real>::max()))
    {
        result.representable = result.within_budget = false;
        return std::numeric_limits<Real>::max();
    }
    const Real rounded = x.convert_to<Real>();
    return x == 0 ? Real(0) : std::nextafter(rounded, std::numeric_limits<Real>::infinity());
}

inline Real convert(const Wide &x, const Wide &envelope, Real &bound, CoulombDiskValue &result)
{
    if (absolute(x) > Wide(std::numeric_limits<Real>::max()))
    {
        result.representable = result.within_budget = false;
        bound = std::numeric_limits<Real>::max();
        return 0; // No caller may consume this candidate when within_budget is false.
    }
    const Real rounded = x.convert_to<Real>();
    bound = outwardPositive(absolute(Wide(rounded) - x) + envelope, result);
    result.within_budget &= bound <= 1e-12;
    return rounded;
}
} // namespace coulomb_disk_detail

inline CoulombDiskValue coulombDiskNaturalMap(Real normal_impulse,
    const ContactTangentVector &impulse, const ContactTangentVector &relative_velocity, Real mu)
{
    using namespace coulomb_disk_detail;
    if (!std::isfinite(normal_impulse) || !impulse.allFinite() || !relative_velocity.allFinite() ||
        !std::isfinite(mu) || mu < 0)
        throw std::invalid_argument("Coulomb disk: nonfinite input or negative coefficient");

    CoulombDiskValue result;
    const Wide l(normal_impulse), coefficient(mu);
    const Wide radius = coefficient * (l > 0 ? l : Wide(0));
    Wide t[2], w[2], y[2], dt[2][2]{}, dw[2][2]{}, dl[2]{}, value[2];
    for (size_t i = 0; i < 2; ++i)
    {
        t[i] = Wide(impulse[i]); w[i] = Wide(relative_velocity[i]); y[i] = t[i] - w[i];
    }
    const Wide size = boost::multiprecision::sqrt(y[0]*y[0] + y[1]*y[1]);
    const Wide scalar_magnitude = absolute(t[0]) + absolute(t[1]) + absolute(w[0]) + absolute(w[1]) + radius;
    const Wide scalar_envelope = arithmeticEnvelope(scalar_magnitude);
    if (size < radius)
    {
        result.branch = CoulombDiskBranch::StickInterior;
        for (size_t i = 0; i < 2; ++i) { value[i] = w[i]; dw[i][i] = 1; }
    }
    else if (size == 0)
    {
        // A valid limiting derivative from l<0, radius=0. No division by zero.
        result.branch = CoulombDiskBranch::ZeroRadiusOrigin;
        for (size_t i = 0; i < 2; ++i) { value[i] = t[i]; dt[i][i] = 1; }
    }
    else
    {
        result.branch = CoulombDiskBranch::SlipExterior;
        const Wide h[2]{y[0]/size, y[1]/size}, fraction = radius/size;
        for (size_t i = 0; i < 2; ++i)
        {
            value[i] = t[i] - radius*h[i];
            dl[i] = l >= 0 ? -coefficient*h[i] : Wide(0);
            for (size_t j = 0; j < 2; ++j)
            {
                dw[i][j] = fraction * (Wide(i == j) - h[i]*h[j]);
                dt[i][j] = Wide(i == j) - dw[i][j];
            }
        }
    }
    const Wide derivative_envelope = arithmeticEnvelope(1 + coefficient);
    auto derivative = [&](const Wide &x)
    {
        Real bound = 0; const Real rounded = convert(x, derivative_envelope, bound, result);
        result.maximum_derivative_roundoff_bound = std::max(result.maximum_derivative_roundoff_bound, bound);
        return rounded;
    };
    for (size_t i = 0; i < 2; ++i)
    {
        result.residual[i] = convert(value[i], scalar_envelope, result.residual_roundoff_bound[i], result);
        result.normal_derivative[i] = derivative(dl[i]);
        for (size_t j = 0; j < 2; ++j)
        {
            result.impulse_derivative(i,j) = derivative(dt[i][j]);
            result.velocity_derivative(i,j) = derivative(dw[i][j]);
        }
    }
    const Wide tnorm = boost::multiprecision::sqrt(t[0]*t[0] + t[1]*t[1]);
    const Wide wnorm = boost::multiprecision::sqrt(w[0]*w[0] + w[1]*w[1]);
    const Wide product = t[0]*w[0] + t[1]*w[1];
    const Wide work_envelope = arithmeticEnvelope(absolute(t[0]*w[0]) + absolute(t[1]*w[1]) + radius*wnorm);
    result.friction_loss = convert(-product, work_envelope, result.friction_loss_roundoff_bound, result);
    result.cone_excess_upper = outwardPositive((tnorm > radius ? tnorm-radius : Wide(0)) + scalar_envelope, result);
    result.maximum_dissipation_error_upper = outwardPositive(absolute(product + radius*wnorm) + work_envelope, result);
    return result;
}
} // namespace SPH::peridynamics
#endif

/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_EXACT_CONTACT_BALANCE_H
#define PD_EXACT_CONTACT_BALANCE_H

#include "base_data_type.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace SPH::peridynamics
{
namespace exact_contact_detail
{
using Integer = boost::multiprecision::cpp_int;
static_assert(sizeof(Real)==8 && std::numeric_limits<Real>::is_iec559 &&
    std::numeric_limits<Real>::radix==2 && std::numeric_limits<Real>::digits==53,
    "Exact contact balance requires IEEE binary64 Real");

struct BinaryNumber { Integer mantissa; int exponent; };

inline BinaryNumber decode(Real x)
{
    if (!std::isfinite(x)) throw std::invalid_argument("exact contact balance: nonfinite input");
    uint64_t bits; std::memcpy(&bits,&x,sizeof(bits));
    const uint64_t field=(bits>>52)&2047, fraction=bits&((uint64_t(1)<<52)-1);
    Integer mantissa=field ? ((uint64_t(1)<<52)|fraction) : fraction;
    if (bits>>63) mantissa=-mantissa;
    return {mantissa,field ? int(field)-1023-52 : -1074};
}

/** Round an exact positive integer ratio times2^exponent to binary64,
 * ties-to-even. All rounding decisions use integer quotient/remainder.
 * ldexp only scales an exactly representable <=53bit integer to another
 * exactly representable binary64 number; no current rounding mode is used. */
struct RoundedRatio { Real value=0, absolute_error_bound=0; };
inline RoundedRatio roundRatio(Integer numerator, const Integer &denominator, int exponent)
{
    if (denominator<=0) throw std::invalid_argument("exact contact balance: denominator");
    if (numerator==0) return {};
    const bool negative=numerator<0; if(negative)numerator=-numerator;
    int highest=int(boost::multiprecision::msb(numerator))-int(boost::multiprecision::msb(denominator))+exponent;
    const int comparison_shift=highest-exponent;
    const bool at_least=comparison_shift>=0 ? numerator>=(denominator<<comparison_shift) :
        (numerator<<(-comparison_shift))>=denominator;
    if(!at_least)--highest;
    if(highest>1023)throw std::overflow_error("exact contact balance: normalized overflow");
    const int unit=std::max(highest-52,-1074), shift=exponent-unit;
    Integer dividend=numerator, divisor=denominator;
    if(shift>=0)dividend<<=shift;else divisor<<=(-shift);
    Integer quotient=dividend/divisor, remainder=dividend%divisor;
    const Integer twice=remainder<<1;
    if(twice>divisor || (twice==divisor && bool(quotient&1)))++quotient;
    if(quotient>(Integer(1)<<53))throw std::logic_error("exact contact balance: rounding range");
    const Real rounded=std::ldexp(Real(quotient.convert_to<uint64_t>()),unit);
    if(!std::isfinite(rounded))throw std::overflow_error("exact contact balance: rounded overflow");
    // A half-unit enclosure is exact except at the subnormal bottom, where
    // one minimum subnormal is the smallest outward binary64 enclosure.
    const Real bound=remainder==0 ? 0 : std::ldexp(Real(1),std::max(unit-1,-1074));
    return {negative ? -rounded : rounded,bound};
}
} // namespace exact_contact_detail

struct ExactContactBalance
{
    Real normalized_value=0, normalized_roundoff_bound=0;
    bool within_budget=true;
};

/** Exact supplied R - sum C_s lambda_s, then exact ratio normalization.
 * Each finite binary64 product has exponent>=-2148, so a fixed2^-2148
 * integer accumulator represents every product and their sum exactly.
 * This changes arithmetic only, never R, C, lambda or the physical mass.
 * The enclosure is ONLY for this accumulation/normalization of supplied
 * floating inputs, not a certificate of material or tangent evaluation. */
inline ExactContactBalance exactContactBalance(Real momentum, const StdVec<Real> &coefficients,
    const StdVec<Real> &impulses, Real impulse_scale)
{
    using namespace exact_contact_detail;
    if(coefficients.size()!=impulses.size() || coefficients.size()>256 || !std::isfinite(impulse_scale) || impulse_scale<=0)
        throw std::invalid_argument("exact contact balance: dimensions/scale");
    const auto r=decode(momentum), scale=decode(impulse_scale);
    Integer sum=r.mantissa<<(r.exponent+2148);
    for(size_t s=0;s<coefficients.size();++s)
    {
        const auto c=decode(coefficients[s]), lambda=decode(impulses[s]);
        if(c.mantissa!=0 && lambda.mantissa!=0)
            sum-=(c.mantissa*lambda.mantissa)<<(c.exponent+lambda.exponent+2148);
    }
    const auto rounded=roundRatio(sum,scale.mantissa,-2148-scale.exponent);
    return {rounded.value,rounded.absolute_error_bound,rounded.absolute_error_bound<=Real(1e-12)};
}
} // namespace SPH::peridynamics
#endif

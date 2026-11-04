/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_EARLY_TANGENT_CYCLE_H
#define PD_EARLY_TANGENT_CYCLE_H
#include <array>
#include <cmath>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Prospective diagnostic displacement, independent of the normal motion.
 * 0 -> +1 -> -1 -> 0 over 2/4/2 nominal time steps, then held at zero.
 * Each join has zero prescribed velocity. It does not assert that any contact
 * stays active: the physical trajectory must prove same-row reversal itself.
 * Kind2 uses this diagnostic; the original kind0/1 motion is unchanged.
 */
struct EarlyTangentCycle
{
    double nominal_dt=.01;

    std::array<double,3> ends() const
    { return {2*nominal_dt,6*nominal_dt,8*nominal_dt}; }

    void validate() const
    {
        if(!std::isfinite(nominal_dt) || nominal_dt<=0 || !std::isfinite(8*nominal_dt))
            throw std::invalid_argument("early tangent cycle: time scale");
    }

    double at(double time) const
    {
        validate();
        if(!std::isfinite(time) || time<0)
            throw std::invalid_argument("early tangent cycle: time cursor");
        const auto t=ends();
        if(time>=t[2])return 0;
        const size_t phase=time<t[0] ? 0 : (time<t[1] ? 1 : 2);
        const double begin=phase==0 ? 0 : t[phase-1];
        const double u=(time-begin)/(t[phase]-begin),s=u*u*(3-2*u);
        return phase==0 ? s : (phase==1 ? 1-2*s : -1+s);
    }
};
}
#endif

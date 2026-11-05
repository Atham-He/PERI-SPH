/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_PHYSICAL_PROTOCOL_H
#define PD_BOUNDED_PHYSICAL_PROTOCOL_H
#include "pd_bounded_friction_protocol.h"
#include "pd_physical_time_schedule.h"

namespace SPH::peridynamics
{
struct BoundedSessionInterval : BoundedFrictionInterval
{
    uint64_t old_tick=0,new_tick=0,time_denominator=0;
};

/** Explicit v2 physical protocol. The integer grid is reconstructed in at
 * most six event segments, independent of the accepted step count. Nothing
 * here constructs reference/neighbor/material state or advances a solver.
 */
struct BoundedPhysicalTimeProtocol
{
    PhysicalTimeDefinition physical;
    uint64_t nominal_ticks=5000;
    static constexpr uint64_t maximum_steps=98304;
    std::vector<uint64_t> events() const
    {
        std::vector<uint64_t> e(physical.normal_ends.begin(),physical.normal_ends.end());
        e.insert(e.end(),physical.tangent_ends.begin(),physical.tangent_ends.end());
        std::sort(e.begin(),e.end());e.erase(std::unique(e.begin(),e.end()),e.end());return e;
    }
    uint64_t steps() const
    {
        (void)PhysicalTimeSchedule{physical};
        if(nominal_ticks==0 || nominal_ticks>PhysicalTimeSchedule::maximum_ticks)
            throw std::invalid_argument("physical PD protocol: nominal tick range");
        uint64_t begin=0,total=0;
        for(uint64_t end:events())
        {
            if(end<=begin)throw std::invalid_argument("physical PD protocol: event order");
            const uint64_t count=1+(end-begin-1)/nominal_ticks;
            if(count>maximum_steps-total)throw std::invalid_argument("physical PD protocol: accepted step budget");
            total+=count;begin=end;
        }
        return total;
    }
    void validate() const {(void)steps();}
    uint64_t tickAtStep(uint64_t step) const
    {
        validate();
        if(step>steps())throw std::invalid_argument("physical PD protocol: accepted step cursor");
        uint64_t begin=0;
        for(uint64_t end:events())
        {
            const uint64_t count=1+(end-begin-1)/nominal_ticks;
            if(step<=count)return begin+(step==count?end-begin:step*nominal_ticks);
            step-=count;begin=end;
        }
        throw std::invalid_argument("physical PD protocol: missing terminal cursor");
    }
    Real timeAtStep(uint64_t step) const {return PhysicalTimeSchedule(physical).time(tickAtStep(step));}
    Real nominalDt() const {return Real(nominal_ticks)/Real(physical.denominator);}
    Real endTime() const {return PhysicalTimeSchedule(physical).time(physical.normal_ends.back());}
    ContactTangentVector tangentStroke() const
    {return ContactTangentVector(physical.tangent_stroke[0],physical.tangent_stroke[1]);}
    std::array<Real,3> phaseEnds() const
    {
        const PhysicalTimeSchedule clock(physical);
        return {clock.time(physical.normal_ends[0]),clock.time(physical.normal_ends[1]),clock.time(physical.normal_ends[2])};
    }
    std::array<uint64_t,3> phaseStepCounts() const
    {
        validate();std::array<uint64_t,3> result{0,0,0};uint64_t begin=0;
        for(uint64_t end:events())
        {
            const size_t phase=begin<physical.normal_ends[0]?0:(begin<physical.normal_ends[1]?1:2);
            result[phase]+=1+(end-begin-1)/nominal_ticks;begin=end;
        }
        return result;
    }
    BoundedSessionInterval interval(uint64_t step,Real accepted_time) const
    {
        const PhysicalTimeSchedule clock(physical);
        const auto v=clock.interval(tickAtStep(step),nominal_ticks);
        const auto same=[](Real a,Real b){return std::memcmp(&a,&b,sizeof(Real))==0;};
        if(!same(accepted_time,v.old_time) || !same(accepted_time+v.actual_dt,v.new_time))
            throw std::invalid_argument("physical PD protocol: exact accepted binary clock");
        BoundedSessionInterval r;r.step=step;r.phase=v.normal_phase;r.old_tick=v.old_tick;r.new_tick=v.new_tick;
        r.old_time=v.old_time;r.new_time=v.new_time;r.actual_dt=v.actual_dt;r.time_denominator=physical.denominator;
        r.old_offset[1]=v.old_wall.normal_offset;r.new_offset[1]=v.new_wall.normal_offset;
        r.normal_cubic[1].coefficients=v.normal_cubic;
        for(size_t a=0;a<2;++a)
        {
            r.old_tangent_displacement[1][a]=v.old_wall.tangent_displacement[a];
            r.new_tangent_displacement[1][a]=v.new_wall.tangent_displacement[a];
            r.wall_tangent_velocity[1][a]=v.tangent_interval_velocity[a];
        }
        return r;
    }
    uint64_t phaseAtStep(uint64_t step) const
    {
        const auto t=tickAtStep(step);
        return t==physical.normal_ends.back()?3:(t<physical.normal_ends[0]?0:(t<physical.normal_ends[1]?1:2));
    }
};
} // namespace SPH::peridynamics
#endif

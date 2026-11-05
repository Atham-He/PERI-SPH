/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_PHYSICAL_TIME_SCHEDULE_H
#define PD_PHYSICAL_TIME_SCHEDULE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace SPH::peridynamics
{
/** A bounded rational physical clock. Time is tick/denominator, not step*dt.
 * The denominator is explicit input, not an implicit quantization or rounding
 * of a requested double. Unsupported tick ranges are rejected. A future PD
 * driver must bind this entire schedule AND numerical step in its CP contract.
 * This independent component does not confer physical restart qualification.
 */
struct PhysicalTimeDefinition
{
    uint64_t denominator=1000000;
    std::array<uint64_t,3> normal_ends{160000,240000,400000};
    std::array<uint64_t,3> tangent_ends{20000,60000,80000};
    double normal_stroke=.01;
    std::array<double,2> tangent_stroke{.0006,-.0004};
};

struct PhysicalPlateSample
{
    double normal_offset=-1,normal_velocity=0;
    std::array<double,2> tangent_displacement{0,0},tangent_velocity{0,0};
};

struct PhysicalTimeInterval
{
    uint64_t old_tick=0,new_tick=0,nominal_ticks=0;
    unsigned normal_phase=0;
    double old_time=0,new_time=0,actual_dt=0;
    PhysicalPlateSample old_wall,new_wall;
    // Polynomial in unit interval s, for whole-interval gap observation.
    std::array<double,4> normal_cubic{-1,0,0,0};
    std::array<double,2> tangent_interval_velocity{0,0};
};

class PhysicalTimeSchedule
{
  public:
    static constexpr uint64_t maximum_ticks=uint64_t(1)<<48;
    explicit PhysicalTimeSchedule(const PhysicalTimeDefinition &definition={}) : definition_(definition)
    {
        check(definition_.denominator>0 && definition_.denominator<=maximum_ticks,"denominator");
        check(increasing(definition_.normal_ends) && increasing(definition_.tangent_ends),"event order");
        check(definition_.normal_ends.back()<=maximum_ticks &&
              definition_.tangent_ends.back()<=definition_.normal_ends.back(),"event range");
        check(std::isfinite(definition_.normal_stroke) && definition_.normal_stroke>0 &&
              definition_.normal_stroke<.1,"normal stroke");
        double norm2=0;
        for(double x:definition_.tangent_stroke){check(std::isfinite(x),"tangent stroke");norm2+=x*x;}
        check(norm2<.01,"tangent stroke range");
    }

    const PhysicalTimeDefinition &definition() const noexcept { return definition_; }
    uint64_t endTick() const noexcept { return definition_.normal_ends.back(); }
    double time(uint64_t tick) const
    {
        check(tick<=endTick(),"time cursor range");
        return double(tick)/double(definition_.denominator);
    }

    PhysicalPlateSample sample(uint64_t tick) const
    {
        check(tick<=endTick(),"sample cursor range");
        const auto normal=shape(tick,definition_.normal_ends,{0,1,1,0});
        const auto tangent=shape(tick,definition_.tangent_ends,{0,1,-1,0});
        PhysicalPlateSample result;
        result.normal_offset=-1+definition_.normal_stroke*normal[0];
        result.normal_velocity=definition_.normal_stroke*normal[1];
        for(size_t k=0;k<2;++k)
        {
            result.tangent_displacement[k]=definition_.tangent_stroke[k]*tangent[0];
            result.tangent_velocity[k]=definition_.tangent_stroke[k]*tangent[1];
        }
        return result;
    }

    PhysicalTimeInterval interval(uint64_t accepted_tick,uint64_t nominal_ticks) const
    {
        check(accepted_tick<endTick(),"completed or invalid cursor");
        check(nominal_ticks>0 && nominal_ticks<=maximum_ticks,"nominal tick range");
        uint64_t boundary=endTick();
        for(const auto &events:{definition_.normal_ends,definition_.tangent_ends})
            for(uint64_t event:events)if(event>accepted_tick)boundary=std::min(boundary,event);
        PhysicalTimeInterval result;
        result.old_tick=accepted_tick;
        result.new_tick=accepted_tick+std::min(nominal_ticks,boundary-accepted_tick);
        result.nominal_ticks=nominal_ticks;
        result.old_time=time(result.old_tick);result.new_time=time(result.new_tick);
        // Tell the physical integrator the ACTUAL binary interval. Both exact
        // rational tick endpoints are retained for identity and independent QA.
        result.actual_dt=result.new_time-result.old_time;
        check(std::isfinite(result.actual_dt) && result.actual_dt>0,"unrepresentable interval");
        result.old_wall=sample(result.old_tick);result.new_wall=sample(result.new_tick);
        const auto &events=definition_.normal_ends;
        result.normal_phase=accepted_tick<events[0]?0:(accepted_tick<events[1]?1:2);
        const unsigned phase=result.normal_phase;
        const uint64_t begin=phase==0?0:events[phase-1];
        const double duration=double(events[phase]-begin);
        const double u=double(accepted_tick-begin)/duration;
        const double h=double(result.new_tick-result.old_tick)/duration;
        result.normal_cubic[0]=result.old_wall.normal_offset;
        if(phase!=1)
        {
            const double amplitude=phase==0?definition_.normal_stroke:-definition_.normal_stroke;
            result.normal_cubic[1]=amplitude*(6*u-6*u*u)*h;
            result.normal_cubic[2]=amplitude*(3-6*u)*h*h;
            result.normal_cubic[3]=-2*amplitude*h*h*h;
        }
        for(size_t k=0;k<2;++k)
            result.tangent_interval_velocity[k]=
                (result.new_wall.tangent_displacement[k]-result.old_wall.tangent_displacement[k])/result.actual_dt;
        return result;
    }

    /** Canonical schedule-only packet, not a complete PD checkpoint. Exact
     * expected bytes bind all events/amplitudes AND the numerical step. Decode
     * stages privately and rejects before assignment; a rehashed alteration is
     * still a contract mismatch. Hash authentication belongs to the outer CP.
     */
    std::vector<uint8_t> packet(uint64_t nominal_ticks,uint64_t accepted_tick) const
    {
        check(nominal_ticks>0 && nominal_ticks<=maximum_ticks && accepted_tick<=endTick(),"packet cursor");
        std::vector<uint8_t> result{'P','D','T','I','M','E',1,0};
        append(result,definition_.denominator);
        for(auto x:definition_.normal_ends)append(result,x);
        for(auto x:definition_.tangent_ends)append(result,x);
        appendDouble(result,definition_.normal_stroke);
        for(auto x:definition_.tangent_stroke)appendDouble(result,x);
        append(result,nominal_ticks);append(result,accepted_tick);
        return result;
    }

    static uint64_t restoreCursor(const std::vector<uint8_t> &bytes,
                                  const PhysicalTimeSchedule &expected,uint64_t nominal_ticks)
    {
        check(bytes.size()==104,"packet length");
        const std::vector<uint8_t> magic{'P','D','T','I','M','E',1,0};
        check(std::equal(magic.begin(),magic.end(),bytes.begin()),"packet magic/version");
        size_t at=8;PhysicalTimeDefinition d;d.denominator=take(bytes,at);
        for(auto &x:d.normal_ends)x=take(bytes,at);
        for(auto &x:d.tangent_ends)x=take(bytes,at);
        d.normal_stroke=takeDouble(bytes,at);
        for(auto &x:d.tangent_stroke)x=takeDouble(bytes,at);
        const auto stored_step=take(bytes,at),cursor=take(bytes,at);
        const PhysicalTimeSchedule staged(d);
        check(stored_step==nominal_ticks,"numerical step contract mismatch");
        check(staged.packet(stored_step,cursor)==expected.packet(nominal_ticks,cursor),"physical schedule contract mismatch");
        return cursor;
    }

  private:
    PhysicalTimeDefinition definition_;
    static void check(bool ok,const char *why)
    { if(!ok)throw std::invalid_argument(std::string("physical time schedule: ")+why); }
    static bool increasing(const std::array<uint64_t,3> &v)
    { return v[0]>0 && v[0]<v[1] && v[1]<v[2]; }
    std::array<double,2> shape(uint64_t tick,const std::array<uint64_t,3> &ends,
                               const std::array<double,4> &values) const
    {
        if(tick>=ends.back())return {values.back(),0};
        const unsigned phase=tick<ends[0]?0:(tick<ends[1]?1:2);
        const uint64_t begin=phase==0?0:ends[phase-1];
        const double duration=double(ends[phase]-begin);
        const double u=double(tick-begin)/duration,amplitude=values[phase+1]-values[phase];
        return {values[phase]+amplitude*u*u*(3-2*u),
                amplitude*6*u*(1-u)*double(definition_.denominator)/duration};
    }
    static void append(std::vector<uint8_t> &bytes,uint64_t x)
    { for(unsigned k=0;k<8;++k)bytes.push_back(uint8_t(x>>(8*k))); }
    static uint64_t take(const std::vector<uint8_t> &bytes,size_t &at)
    {
        uint64_t x=0;for(unsigned k=0;k<8;++k)x|=uint64_t(bytes.at(at++))<<(8*k);return x;
    }
    static void appendDouble(std::vector<uint8_t> &bytes,double x)
    { uint64_t bits=0;static_assert(sizeof(x)==sizeof(bits));std::memcpy(&bits,&x,8);append(bytes,bits); }
    static double takeDouble(const std::vector<uint8_t> &bytes,size_t &at)
    { const uint64_t bits=take(bytes,at);double x;std::memcpy(&x,&bits,8);return x; }
};
} // namespace SPH::peridynamics
#endif

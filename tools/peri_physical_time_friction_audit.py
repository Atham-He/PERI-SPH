"""Explicit v2 physical-time front end to the unchanged full PD math audit.

Only clock/header admission is replaced. Original every-trial material tangent,
non-symmetric R/J, true mass, contact and physical/algorithm ledger assessment
is inherited verbatim. No native solver is called by this reader.
"""
import struct
from decimal import Decimal as D
from fractions import Fraction as F
from peri_bounded_friction_audit import FrictionInspector,inspect_stream

def schedule_grid(case):
    p=case['physical_time'];den=p['denominator'];n=p['normal_event_ticks'];t=p['tangent_event_ticks'];dt=p['nominal_ticks']
    assert isinstance(den,int) and 0<den<=2**48 and isinstance(dt,int) and 0<dt<=2**48
    for v in (n,t):assert len(v)==3 and all(isinstance(x,int) for x in v) and 0<v[0]<v[1]<v[2]<=2**48
    assert t[-1]<=n[-1]
    events=sorted(set(n+t));grid=[0]
    while grid[-1]<n[-1]:
        assert len(grid)<=98304
        boundary=next(e for e in events if e>grid[-1]);grid.append(min(grid[-1]+dt,boundary))
    return grid

def decimal(x):
    if isinstance(x,F):return D(x.numerator)/D(x.denominator)
    return D.from_float(x) if isinstance(x,float) else D(x)

def shape(t,events,values):
    if t>=events[-1]:return F(values[-1])
    k=next(k for k,e in enumerate(events) if t<e);begin=0 if k==0 else events[k-1]
    u=F(t-begin,events[k]-begin)
    return F(values[k])+F(values[k+1]-values[k])*(3*u*u-2*u*u*u)

def protocol(case,step,time):
    p=case['physical_time'];grid=schedule_grid(case)
    assert 0<=step<len(grid)-1
    old,new=grid[step:step+2];den=p['denominator'];normal=p['normal_event_ticks'];tangent=p['tangent_event_ticks']
    old_time,new_time=float(F(old,den)),float(F(new,den));dt=new_time-old_time
    assert struct.pack('<d',time)==struct.pack('<d',old_time)
    assert struct.pack('<d',old_time+dt)==struct.pack('<d',new_time) and dt>0
    phase=next(k for k,e in enumerate(normal) if old<e);begin=0 if phase==0 else normal[phase-1]
    u=F(old-begin,normal[phase]-begin);h=F(new-old,normal[phase]-begin)
    a,b=shape(old,normal,(0,1,1,0)),shape(new,normal,(0,1,1,0));stroke=F(case['stroke'])
    cubic=[-1+stroke*a,F(0),F(0),F(0)]
    if phase!=1:
        amplitude=stroke*(1 if phase==0 else -1)
        cubic[1:]=[amplitude*(6*u-6*u*u)*h,amplitude*(3-6*u)*h*h,-2*amplitude*h*h*h]
    ta,tb=shape(old,tangent,(0,1,-1,0)),shape(new,tangent,(0,1,-1,0))
    strokes=list(map(F,case['tangent_stroke']))
    return dict(phase=phase,step=step,old_time=old_time,new_time=new_time,actual_dt=dt,
                old_tick=old,new_tick=new,time_denominator=den,
                old_offset=[D(-1),decimal(-1+stroke*a)],new_offset=[D(-1),decimal(-1+stroke*b)],
                normal_cubic=[[D(-1),D(0),D(0),D(0)],list(map(decimal,cubic))],
                old_tangent_displacement=[[D(0),D(0)],[decimal(x*ta) for x in strokes]],
                new_tangent_displacement=[[D(0),D(0)],[decimal(x*tb) for x in strokes]])

class PhysicalTimeFrictionInspector(FrictionInspector):
    def expected_end_time(self):
        p=self.case['physical_time'];return decimal(F(p['normal_event_ticks'][-1],p['denominator']))

    def validate_protocol_header(self,case,header):
        p=case['physical_time'];grid=schedule_grid(case);den=p['denominator']
        assert header['format']==2 and header['protocol_kind']==case['protocol_kind']==2
        assert header['physical_clock_profile']=='explicit_rational_events_v2'
        for h,c in [('time_denominator','denominator'),('nominal_ticks','nominal_ticks'),
                    ('normal_event_ticks','normal_event_ticks'),('tangent_event_ticks','tangent_event_ticks')]:
            assert header[h]==p[c],('physical clock header',h)
        assert case['steps']==header['steps_requested']==len(grid)-1
        assert case['dt']==header['nominal_dt']==float(F(p['nominal_ticks'],den))
        ends=p['normal_event_ticks'];counts=[0,0,0]
        for t in grid[:-1]:counts[next(k for k,e in enumerate(ends) if t<e)]+=1
        assert header['phase_step_counts']==counts
        assert header['phase_end_times']==[float(F(t,den)) for t in ends]
        assert header['end_time']==float(F(ends[-1],den))
        assert header['tangent_cycle_cursor']=='fixed_rational_C1_early_cycle_then_zero'
        assert header['tangent_cycle_factors']==[0,1,-1,0]
        assert header['tangent_cycle_end_times']==[float(F(t,den)) for t in p['tangent_event_ticks']]
        cp=case['checkpoint_steps'];assert cp==header['checkpoint_steps'] and 0<len(cp)<=64
        assert cp==sorted(set(cp)) and cp[0]>0 and cp[-1]==len(grid)-1

    def expected_protocol(self,step,time):return protocol(self.case,step,time)

    def interval(self,interval,old):
        expected=self.expected_protocol(old['step'],old['time'])
        for field in ('old_tick','new_tick','time_denominator'):
            assert interval[field]==expected[field],('physical clock accepted cursor',field)
        return super().interval(interval,old)

def inspect_physical_stream(case,path,full=True,capture=None):
    assert 'physical_time' in case and case['protocol_kind']==2
    return inspect_stream(case,path,0,full,capture,inspector_type=PhysicalTimeFrictionInspector)

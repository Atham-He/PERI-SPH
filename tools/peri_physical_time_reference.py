"""Exact-rational, zero-physics assessment of a NEW frozen clock component.

This is not a PD material/trajectory/complete-checkpoint qualification.
No executable is run by this module; its controller supplies recorded packets.
"""
import json
import math
import struct
from fractions import Fraction as F

DEN=1000000
NORMAL=(160000,240000,400000)
TANGENT=(20000,60000,80000)
STROKES=(.01,.0006,-.0004)

def smooth(t,events,values):
    t=F(t)
    if t>=events[-1]:
        return F(values[-1]),F(0)
    k=next(k for k,e in enumerate(events) if t<e)
    begin=0 if k==0 else events[k-1]
    duration=F(events[k]-begin)
    u=(t-begin)/duration
    a=F(values[k+1]-values[k])
    return F(values[k])+a*(3*u*u-2*u*u*u),a*(6*u-6*u*u)*DEN/duration

def wall(t):
    n,nv=smooth(t,NORMAL,(0,1,1,0))
    u,uv=smooth(t,TANGENT,(0,1,-1,0))
    sn,s0,s1=map(F.from_float,STROKES)
    return [-1+sn*n,sn*nv,s0*u,s1*u,s0*uv,s1*uv]

def expected_interval(t,nominal):
    boundary=min(e for e in (*NORMAL,*TANGENT) if e>t)
    end=t+min(nominal,boundary-t)
    phase=next(k for k,e in enumerate(NORMAL) if t<e)
    begin=0 if phase==0 else NORMAL[phase-1]
    duration=F(NORMAL[phase]-begin)
    u=F(t-begin)/duration
    h=F(end-t)/duration
    old,new=wall(t),wall(end)
    a=F.from_float(.01)*(1 if phase==0 else -1)
    cubic=[old[0],F(0),F(0),F(0)] if phase==1 else [old[0],a*(6*u-6*u*u)*h,a*(3-6*u)*h*h,-2*a*h*h*h]
    return dict(old_tick=t,new_tick=end,phase=phase,old_time=F(t,DEN),new_time=F(end,DEN),
                actual_dt=F(end-t,DEN),old_wall=old,new_wall=new,normal_cubic=cubic,
                tangent_interval_velocity=[(new[k]-old[k])/F(end-t,DEN) for k in (2,3)])

def normalized(x,y):
    assert isinstance(x,(int,float)) and math.isfinite(x),('nonfinite/wrong type',x)
    return float(abs(F(x)-y)/max(F(1),abs(y)))

def assess_grid(raw,nominal,start=0):
    records=[json.loads(x) for x in raw.splitlines()]
    assert records
    t=start
    maximum=0.
    for r in records:
        expect=expected_interval(t,nominal)
        assert set(r)==set(expect),('packet fields',set(r)^set(expect))
        for k in ('old_tick','new_tick','phase'):
            assert r[k]==expect[k],(k,r[k],expect[k])
        for k in ('old_time','new_time','actual_dt'):
            maximum=max(maximum,normalized(r[k],expect[k]))
        assert r['actual_dt']==r['new_time']-r['old_time']>0
        for k in ('old_wall','new_wall','normal_cubic','tangent_interval_velocity'):
            assert len(r[k])==len(expect[k])
            maximum=max(maximum,*(normalized(a,b) for a,b in zip(r[k],expect[k])))
        # Check the polynomial on eleven rational interior coordinates as well
        # as endpoints, against an independently evaluated full physical path.
        for j in range(11):
            x=F(j,10)
            polynomial=sum(F(r['normal_cubic'][i])*x**i for i in range(4))
            true=wall(F(t)+(expect['new_tick']-t)*x)[0]
            maximum=max(maximum,float(abs(polynomial-true)/max(F(1),abs(true))))
        t=expect['new_tick']
    assert t==NORMAL[-1],('missing physical tail',t)
    assert maximum<=1e-12,('clock normalized error',maximum)
    return dict(intervals=len(records),last_tick=t,last_time=records[-1]['new_time'],
                maximum_normalized_error=maximum,all_whole_normal_cubics_checked=True,
                nominal_ticks=nominal,start_tick=start)

def packet(nominal,tick=20000):
    return b'PDTIME\x01\x00'+struct.pack('<7Q3d2Q',DEN,*NORMAL,*TANGENT,*STROKES,nominal,tick)

def negative_packets():
    base=packet(10000)
    def u64(offset,value):
        b=bytearray(base);struct.pack_into('<Q',b,offset,value);return bytes(b)
    def f64(offset,value):
        b=bytearray(base);struct.pack_into('<d',b,offset,value);return bytes(b)
    return [
        ('zero_denominator',u64(8,0),'denominator'),
        ('event_order',u64(16,0),'event order'),
        ('normal_time_changed',u64(16,160001),'physical schedule contract mismatch'),
        ('tangent_time_changed',u64(40,20001),'physical schedule contract mismatch'),
        ('normal_stroke_changed',f64(64,.0101),'physical schedule contract mismatch'),
        ('tangent_stroke_changed',f64(72,.0007),'physical schedule contract mismatch'),
        ('nominal_step_changed',u64(88,5000),'numerical step contract mismatch'),
        ('cursor_outside',u64(96,400001),'packet cursor'),
        ('cursor_overflow',u64(96,2**64-1),'packet cursor'),
        ('normal_nan',f64(64,float('nan')),'normal stroke'),
        ('tangent_inf',f64(72,float('inf')),'tangent stroke'),
        ('tangent_event_after_end',u64(56,400001),'event range'),
        ('denominator_outside',u64(8,2**48+1),'denominator'),
        ('wrong_version',base[:6]+b'\x02'+base[7:],'packet magic/version'),
        ('short_packet',base[:-1],'packet length'),
        ('extra_tail',base+b'\x00','packet length'),
    ]

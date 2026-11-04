"""Independent high-precision S2 reader. No native solver/controller imports.

Decimal closed-form scalar active/box solutions and exponential convolution,
plus continuous quintic-wall and piecewise-Coulomb references. Every physical
ledger is recomputed using actual M. CSV rows never become the reference state.
"""
import csv
from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import struct


def rows(path):
    with Path(path).open() as stream:
        return [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]


def scalar_interval(q,x,qn,dt,k0,kv,tau):
    h=dt/tau; decay=(-h).exp(); alpha=(1-decay)/h
    dx=qn-q; steady=tau*dx/dt; transient=x-steady
    xn=steady+transient*decay
    average=k0*(q+qn)/2+kv*(steady+transient*alpha)
    old=(k0*q*q+kv*x*x)/2; new=(k0*qn*qn+kv*xn*xn)/2
    # Independent dimensional exponential antiderivative, unlike native series.
    integral=steady*steady*dt+2*steady*transient*tau*(1-decay)+transient*transient*tau*(1-decay*decay)/2
    diss=kv/tau*integral
    return xn,average,old,new,diss


def wall(kind,t,total,displacement,dt):
    if kind==0:
        return displacement*t if t<=1 else displacement if t<=2 else displacement*(3-t) if t<=3 else D(0)
    if kind==1: return displacement
    if kind==4: return D(0)
    if kind==5: return -max(D(0),t-8*dt)
    s=t/total
    return displacement*(10*s**3-15*s**4+6*s**5)


def tangent_wall(kind,t,boost):
    first,second=(D('.2'),D('.5')) if kind==4 else (D('.6'),D('1.2'))
    return boost+(D('.05') if t<first else D('-.025') if t<second else D('.1'))


def continuum(t,total,displacement,mass,k0,kv,tau):
    """Polynomial particular solution of x'+x/tau=w', with x(0)=0."""
    coefficients={3:10*displacement/total**3,4:-15*displacement/total**4,5:6*displacement/total**5}
    def derivative(at,order):
        return sum((coefficient*D(math.factorial(n)//math.factorial(n-order))*(at**(n-order) if n>order else D(1))
                    for n,coefficient in coefficients.items() if n>=order),D(0))
    def particular(at):
        return sum(((-1)**(order-1)*tau**order*derivative(at,order) for order in range(1,6)),D(0))
    x=particular(t)-particular(D(0))*(-t/tau).exp()
    q=derivative(t,0); velocity=derivative(t,1)
    force=k0*q+kv*x+mass*derivative(t,2)
    # Integrated force from 0, including the exact convolution integral.
    area=sum((c*t**(n+1)/D(n+1) for n,c in coefficients.items()),D(0))
    impulse=k0*area+kv*tau*(q-x)+mass*velocity
    return q,velocity,x,force,impulse


def assess(case,trace):
    with localcontext() as context:
        context.prec=65
        values=list(map(D,map(str,case['config'])))
        dt,total,mn,mt,k0,kv,tau,displacement,mu,boost,load,q,x,vn,vt=values
        kind=case['kind']; vt+=boost; qt=D(0)
        data=rows(trace); assert len(data)==int(total/dt)
        alpha=(1-(-dt/tau).exp())/(dt/tau)
        metric=mn+dt*(k0*dt/2+kv*tau*(1-alpha))
        worst=0.0; losses=[]; works=[]; material_losses=[]; friction_losses=[]; modes=set()
        force_errors=[[],[]]; force_refs=[[],[]]; velocity_errors=[[],[]]; velocity_refs=[[],[]]
        vt_exact=vt; previous_integral=D(0)
        for index,row in enumerate(data,1):
            t0=(index-1)*dt; t1=index*dt
            w0=wall(kind,t0,total,displacement,dt); w1=wall(kind,t1,total,displacement,dt)
            uw=tangent_wall(kind,(t0+t1)/2,boost)
            if kind<=1:
                qn=w1; vnn=(w1-w0)/dt; vtn=vt; qtn=qt; jn=jt=D(0)
                kinetic0=kinetic1=algorithm=normal_loss=friction_loss=bodywork=D(0)
                slip=D(0); mode=0
            else:
                free=(mn*vn+dt*(load-k0*q-kv*alpha*x))/metric
                bound=(w1-q)/dt
                jn=max(D(0),metric*(bound-free)); vnn=free+jn/metric; qn=q+dt*vnn
                trial=mt*(uw-vt); limit=mu*jn
                jt=min(limit,max(-limit,trial)); vtn=vt+jt/mt; qtn=qt+dt*vtn
                slip=vtn-uw
                mode=0 if jn==0 else 1 if abs(trial)<=limit else 2 if slip>0 else 3
                kinetic0=(mn*vn*vn+mt*vt*vt)/2
                kinetic1=(mn*vnn*vnn+mt*vtn*vtn)/2
                algorithm=(mn*(vnn-vn)**2+mt*(vtn-vt)**2)/2
                normal_loss=jn*(q-w0)/dt; friction_loss=-jt*slip; bodywork=load*(qn-q)
                assert jn>=0 and qn-w1>=D('-1e-54') and normal_loss>=D('-1e-54') and friction_loss>=D('-1e-54')
                assert abs(jn*(vnn-bound))<D('1e-54')
            xn,average,stored0,stored1,physical=scalar_interval(q,x,qn,dt,k0,kv,tau)
            assert physical>=0
            work=average*(qn-q) if kind<=1 else jn*(w1-w0)/dt+jt*uw
            assert abs(work+bodywork-(kinetic1-kinetic0+stored1-stored0)-physical-friction_loss-normal_loss-algorithm)<D('1e-54')
            expected=dict(step=D(index),t0=t0,t1=t1,dt=dt,q0=q,vn0=vn,x0=x,qt0=qt,vt0=vt,
                q1=qn,vn1=vnn,x1=xn,qt1=qtn,vt1=vtn,wall0=w0,wall1=w1,wall_ut=uw,
                jn=jn,jt=jt,fn=jn/dt,ft=jt/dt,material_force=average,U0=stored0,U1=stored1,
                K0=kinetic0,K1=kinetic1,Wwall=work,Wbody=bodywork,Dmaterial=physical,
                Dfriction=friction_loss,Dnormal=normal_loss,Dalgorithm=algorithm,Eres=D(0),Pnres=D(0),Ptres=D(0),
                primal=D(0),natural=D(0),stationarity=D(0),complementarity=D(0),cone=D(0),friction_natural=D(0),slip=slip)
            for field,value in expected.items():
                observed=row[field]; assert math.isfinite(observed),(case['name'],index,field,'nonfinite')
                if field=='step': scale=1
                elif field in ('t0','t1','dt'): scale=float(total)
                elif field in ('q0','q1','x0','x1','qt0','qt1','wall0','wall1'): scale=.2
                elif field in ('vn0','vn1','vt0','vt1','wall_ut','primal','natural','slip'): scale=.1
                elif field in ('jn','jt','Pnres','Ptres','stationarity','cone','friction_natural'): scale=.001
                elif field in ('fn','ft','material_force'): scale=.2
                else: scale=.0001 if kind in (4,5) else .04
                error=abs(observed-float(value))/scale
                worst=max(worst,error)
                assert error<=1e-10,(case['name'],index,field,error)
            if int(row['mode'])!=mode:
                assert kind>1 and abs(row['slip'])<1e-11 and abs(abs(row['jt'])-float(mu)*row['jn'])<1e-13, (case['name'],index,'mode')
            modes.add(mode); losses.append(row['Dalgorithm']); material_losses.append(row['Dmaterial']); friction_losses.append(row['Dfriction'])
            works.append(max(0,row['Wwall']+row['Wbody']))
            if kind in (2,3):
                assert k0/2>6*mn/total**2 and kv>=0, 'unproven maintained contact reference'
                _,vnref,_,_,integral=continuum(t1,total,displacement,mn,k0,kv,tau)
                _,_,_,fnref,_=continuum((t0+t1)/2,total,displacement,mn,k0,kv,tau)
                # U is constant on each interval; the integral is monotone because R>=0.
                delta=mu*(integral-previous_integral)/mt
                assert delta>=0
                vt_next=uw-max(D(0),uw-vt_exact-delta) if vt_exact<=uw else uw+max(D(0),vt_exact-uw-delta)
                ftref=mt*(vt_next-vt_exact)/dt
                for channel,field,reference in [(0,'fn',fnref),(1,'ft',ftref)]:
                    force_errors[channel].append(row[field]-float(reference)); force_refs[channel].append(float(reference))
                for channel,field,reference in [(0,'vn1',vnref),(1,'vt1',vt_next)]:
                    velocity_errors[channel].append(row[field]-float(reference)); velocity_refs[channel].append(float(reference))
                vt_exact=vt_next; previous_integral=integral
            elif kind==4:
                acceleration=mu*(-load)/mt
                vt_next=uw-max(D(0),uw-vt_exact-acceleration*dt) if vt_exact<=uw else uw+max(D(0),vt_exact-uw-acceleration*dt)
                ftref=mt*(vt_next-vt_exact)/dt
                force_errors[0].append(row['fn']-float(-load)); force_refs[0].append(float(-load))
                force_errors[1].append(row['ft']-float(ftref)); force_refs[1].append(float(ftref))
                velocity_errors[1].append(row['vt1']-float(vt_next)); velocity_refs[1].append(float(vt_next))
                vt_exact=vt_next
            elif kind==5:
                attached=index<=8
                fnref=-load if attached else D(0)
                vnref=D(0) if attached else load/mn*(t1-8*dt)
                delta=mu*fnref*dt/mt
                vt_next=uw-max(D(0),uw-vt_exact-delta) if vt_exact<=uw else uw+max(D(0),vt_exact-uw-delta)
                ftref=mt*(vt_next-vt_exact)/dt
                for channel,field,reference in [(0,'fn',fnref),(1,'ft',ftref)]:
                    force_errors[channel].append(row[field]-float(reference)); force_refs[channel].append(float(reference))
                for channel,field,reference in [(0,'vn1',vnref),(1,'vt1',vt_next)]:
                    velocity_errors[channel].append(row[field]-float(reference)); velocity_refs[channel].append(float(reference))
                vt_exact=vt_next
            q,x,vn,qt,vt=qn,xn,vnn,qtn,vtn
        def relative(errors,reference):
            denominator=math.fsum(v*v for v in reference)
            return math.sqrt(math.fsum(v*v for v in errors)/denominator) if denominator>0 else None
        positive_work=math.fsum(works)
        return dict(name=case['name'],kind=kind,dt=float(dt),steps=len(data),maximum_normalized_error=worst,
            force_relative_L2=[relative(a,b) for a,b in zip(force_errors,force_refs)],
            velocity_relative_L2=[relative(a,b) for a,b in zip(velocity_errors,velocity_refs)],
            algorithmic_loss=math.fsum(losses),positive_input_work=positive_work,
            algorithmic_fraction=math.fsum(losses)/positive_work if positive_work else None,
            material_dissipation=math.fsum(material_losses),friction_dissipation=math.fsum(friction_losses),modes=sorted(modes))


def check_checkpoint(path,case,trace):
    layout='<Qi15d5dQQ14d4di'
    raw=Path(path).read_bytes()
    assert len(raw)==struct.calcsize(layout)
    items=struct.unpack(layout,raw)
    assert items[0]==0x5356463200000001 and items[1]==case['kind']
    assert list(items[2:17])==list(map(float,case['config']))
    data=rows(trace); final=data[-1]
    assert items[17:22]==tuple(final[name] for name in ['q1','vn1','x1','qt1','vt1'])
    assert items[22]==len(data) and items[23]==Path(trace).stat().st_size
    for j,fields in enumerate([['Wwall','Wbody'],['Dmaterial'],['Dfriction'],['Dnormal'],['Dalgorithm'],['jn'],['jt']]):
        exact=math.fsum(math.fsum(row[field] for field in fields) for row in data)
        assert abs(items[24+2*j]-exact)<=1e-13*max(.001,abs(exact))
        accumulator=correction=0.0
        for row in data:
            value=sum(row[field] for field in fields)
            y=value-correction; next_value=accumulator+y
            correction=(next_value-accumulator)-y; accumulator=next_value
        assert items[24+2*j]==accumulator and items[25+2*j]==correction
    assert items[38:42]==(final['jn'],final['jt'],final['wall1'],final['wall_ut'])
    assert items[42]==int(final['mode'])
    return dict(size=len(raw),accepted_steps=len(data),full_state_and_ledger_checked=True)

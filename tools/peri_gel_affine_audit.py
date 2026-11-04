"""Independent 65-digit material-path reference.

Uses a polynomial-plus-exponential particular solution and dimensional closed
antiderivative for loss; native uses coefficient/Gram series. Averaged PK1 uses
independent48-point high-precision quadrature. Endpoint tangents are checked by
65-digit central differences of that independent integral, not native formulas.
"""
from decimal import Decimal as D, localcontext
from functools import lru_cache
import math


def identity_matrix(): return [D(int(i//3==i%3)) for i in range(9)]
def add(a,b): return [x+y for x,y in zip(a,b)]
def sub(a,b): return [x-y for x,y in zip(a,b)]
def scale(a,s): return [x*s for x in a]
def dot(a,b): return sum((x*y for x,y in zip(a,b)),D(0))
def transpose(a): return [a[j*3+i] for i in range(3) for j in range(3)]
def multiply(a,b): return [sum((a[i*3+k]*b[k*3+j] for k in range(3)),D(0)) for i in range(3) for j in range(3)]
def deviator(a):
    trace=(a[0]+a[4]+a[8])/3
    return [x-(trace if i//3==i%3 else 0) for i,x in enumerate(a)]
def green(f): return deviator(scale(sub(multiply(transpose(f),f),identity_matrix()),D('.5')))
def determinant(f):
    a,b,c,d,e,z,g,h,i=f
    return a*(e*i-z*h)-b*(d*i-z*g)+c*(d*h-e*g)
def inverse_transpose(f):
    a,b,c,d,e,z,g,h,i=f; j=determinant(f)
    return [x/j for x in [e*i-z*h,z*g-d*i,d*h-e*g,c*h-b*i,a*i-c*g,b*g-a*h,b*z-c*e,c*d-a*z,a*e-b*d]]


@lru_cache(None)
def rule48():
    values=[]
    for i in range(24):
        x=D(str(math.cos(math.pi*(i+.75)/48.5)))
        for _ in range(80):
            before=D(1); current=x
            for n in range(2,49): before,current=current,((2*n-1)*x*current-(n-1)*before)/n
            derivative=48*(x*current-before)/(x*x-1)
            next_value=x-current/derivative
            if abs(next_value-x)<D('1e-60'): x=next_value; break
            x=next_value
        else: raise AssertionError('reference quadrature root failed')
        weight=1/((1-x*x)*derivative*derivative)
        values.extend([((1-x)/2,weight),((1+x)/2,weight)])
    for power in range(13):
        assert abs(sum((w*(s**power if power else D(1)) for s,w in values),D(0))-D(1)/D(power+1))<D('1e-55')
    return tuple(values)


def minimum_det(f0,f1):
    df=sub(f1,f0); coefficients=[D(0)]*4
    for mask in range(8):
        matrix=[df[i] if mask & (1<<(i%3)) else f0[i] for i in range(9)]
        coefficients[mask.bit_count()]+=determinant(matrix)
    candidates=[D(0),D(1)]; a=3*coefficients[3]; b=2*coefficients[2]; c=coefficients[1]
    if a==0:
        if b!=0: candidates.append(-c/b)
    else:
        discriminant=b*b-4*a*c
        if discriminant>=0: candidates.extend([(-b+discriminant.sqrt())/(2*a),(-b-discriminant.sqrt())/(2*a)])
    return min(determinant(add(f0,scale(df,s))) for s in candidates if 0<=s<=1)


def stored(f,x,p):
    bulk,eq,branch,tau=p; j=determinant(f)
    return eq/2*((-D(2)/3*j.ln()).exp()*dot(f,f)-3)+bulk/4*(j*j-1-2*j.ln())+branch*dot(x,x)


def pk1(f,x,p):
    bulk,eq,branch,tau=p; j=determinant(f); inv=inverse_transpose(f)
    return add(add(scale(sub(f,scale(inv,dot(f,f)/3)),eq*(-D(2)/3*j.ln()).exp()),
                   scale(inv,bulk/2*(j*j-1))),scale(multiply(f,x),2*branch))


def basis(f0,f1,x0,h):
    df=sub(f1,f0)
    g1=deviator(scale(add(multiply(transpose(f0),df),multiply(transpose(df),f0)),D('.5')))
    g2=deviator(scale(multiply(transpose(df),df),D('.5')))
    a=sub(scale(g1,1/h),scale(g2,2/(h*h))); b=scale(g2,2/h); c=sub(x0,a)
    return a,b,c


def average(f0,f1,x0,dt,p):
    if dt==0: return pk1(f0,x0,p)
    h=dt/p[3]; a,b,c=basis(f0,f1,x0,h); df=sub(f1,f0); result=[D(0)]*9
    for s,w in rule48():
        x=add(add(a,scale(b,s)),scale(c,(-h*s).exp()))
        result=add(result,scale(pk1(add(f0,scale(df,s)),x,p),w))
    return result


def reference_step(f0,f1,z0,total_loss,dt,p,tangent=True):
    if p[0]<=0 or p[1]<=0 or p[2]<0 or p[3]<=0 or dt<0: raise ValueError('invalid parameters/step')
    min_j=minimum_det(f0,f1)
    if min_j<D('.05'): raise ValueError('invalid determinant inside path')
    e0=green(f0); x0=sub(e0,z0)
    if dt==0:
        if f0!=f1: raise ValueError('nonzero path at zero dt')
        xn=x0; diss=D(0)
    else:
        h=dt/p[3]; a,b,c=basis(f0,f1,x0,h); decay=(-h).exp()
        xn=add(add(a,b),scale(c,decay))
        e0h=(1-decay)/h; e1h=(1-(1+h)*decay)/(h*h); e00h=(1-decay*decay)/(2*h)
        integral=dot(a,a)+dot(a,b)+dot(b,b)/3+2*dot(a,c)*e0h+2*dot(b,c)*e1h+dot(c,c)*e00h
        diss=2*p[2]*h*integral
        assert diss>=D('-1e-45')
    e1=green(f1); z1=sub(e1,xn); avg=average(f0,f1,x0,dt,p)
    u0=stored(f0,x0,p); u1=stored(f1,xn,p); work=dot(avg,sub(f1,f0))
    assert abs(work-(u1-u0+diss))<D('1e-28')
    jacobian=[[D(0)]*9 for _ in range(9)]
    if tangent and dt:
        epsilon=D('1e-18')
        for column in range(9):
            plus=list(f1); minus=list(f1); plus[column]+=epsilon; minus[column]-=epsilon
            jacobian[column]=scale(sub(average(f0,plus,x0,dt,p),average(f0,minus,x0,dt,p)),1/(2*epsilon))
    norm=sum((v*v for col in jacobian for v in col),D(0))
    asym=sum(((jacobian[j][i]-jacobian[i][j])**2 for i in range(9) for j in range(9)),D(0))
    return dict(green=e1,viscous=z1,average_pk1=avg,tangent=jacobian,stored_old=u0,stored_new=u1,
        physical_dissipation=diss,total_dissipation=total_loss+diss,work=work,energy_residual=D(0),minimum_determinant=min_j,
        tangent_asymmetry=(asym/norm).sqrt() if norm else D(0))


def check(case,records,code):
    with localcontext() as context:
        context.prec=65
        assert code==case.get('expected_code',0),(case['name'],'return',code)
        p=list(map(lambda v:D(str(v)),case['parameters'])); f0=list(map(lambda v:D(str(v)),case['initial']))
        old_green=list(map(lambda v:D(str(v)),case['green'])); z0=list(map(lambda v:D(str(v)),case['viscous']))
        loss=D(str(case['old_dissipation'])); entries=[r for r in records if r['type']=='step']
        if code:
            invalid=False
            try:
                if any(abs(a-b)>D('1e-12') for a,b in zip(old_green,green(f0))): raise ValueError('history mismatch')
                for step in case['steps']:
                    ref=reference_step(f0,list(map(lambda v:D(str(v)),step['F'])),z0,loss,D(str(step['dt'])),p,False)
            except ValueError: invalid=True
            assert invalid and records[-1]['type']=='error' and not entries
            return dict(name=case['name'],status='invalid_material_path_rejected')
        assert records[-1]['completed'] and len(entries)==len(case['steps'])
        worst=0.0; max_asym=0.0; valid_tangents=0
        natural_scale=float(max(p[:3])); assert natural_scale>0
        for index,(step,row) in enumerate(zip(case['steps'],entries),1):
            f1=list(map(lambda v:D(str(v)),step['F'])); dt=D(str(step['dt']))
            ref=reference_step(f0,f1,z0,loss,dt,p)
            assert row['step']==index and abs(row['dt']-float(dt))<1e-14
            assert row['tangent_available']==bool(dt)
            for field in ['green','viscous','average_pk1','tangent','stored_old','stored_new','physical_dissipation',
                          'total_dissipation','work','energy_residual','minimum_determinant','tangent_asymmetry']:
                want=ref[field]; got=row[field]
                scale_value=1 if field in ['green','viscous','minimum_determinant','tangent_asymmetry'] else natural_scale
                if field=='tangent':
                    wanted=[v for column in want for v in column]; observed=[v for column in got for v in column]
                elif isinstance(want,list): wanted=want; observed=got
                else: wanted=[want]; observed=[got]
                assert len(wanted)==len(observed)
                for a,b in zip(observed,wanted):
                    error=abs(a-float(b))/scale_value
                    assert math.isfinite(error) and error<=1e-10,(case['name'],index,field,error)
                    worst=max(worst,error)
            assert row['physical_dissipation']>=-1e-10*natural_scale
            max_asym=max(max_asym,row['tangent_asymmetry']); valid_tangents+=bool(dt)
            f0=f1; z0=ref['viscous']; loss=ref['total_dissipation']
        return dict(name=case['name'],status='independent_affine_history_work_and_endpoint_tangent_passed',steps=len(entries),
                    maximum_normalized_error=worst,maximum_consistent_tangent_asymmetry=max_asym,
                    independently_checked_9x9_tangents=valid_tangents,
                    SPD_tangent_permitted=max_asym<=1e-10)

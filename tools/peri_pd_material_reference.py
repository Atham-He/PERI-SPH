"""Independent95digit path reference with cached exponentials, not native Gram."""
from decimal import Decimal as D, getcontext
from functools import lru_cache
from peri_gel_affine_audit import add,sub,scale,dot,green,basis,rule48,stored,pk1,minimum_det

@lru_cache(None)
def time_rule(h):
    assert getcontext().prec>=95
    return tuple((s,w,(-h*s).exp()) for s,w in rule48())

def average(f0,f1,x0,dt,p):
    h=dt/p[3];a,b,c=basis(f0,f1,x0,h);df=sub(f1,f0);result=[D(0)]*9
    for s,w,decay in time_rule(h):
        x=add(add(a,scale(b,s)),scale(c,decay))
        result=add(result,scale(pk1(add(f0,scale(df,s)),x,p),w))
    return result

def material_step(f0,f1,z0,loss0,dt,p,tangent):
    assert getcontext().prec>=95 and dt>0
    min_j=minimum_det(f0,f1);assert min_j>=D('.05')
    x0=sub(green(f0),z0);h=dt/p[3];a,b,c=basis(f0,f1,x0,h);decay=(-h).exp()
    x1=add(add(a,b),scale(c,decay));z1=sub(green(f1),x1)
    e0=(1-decay)/h;e1=(1-(1+h)*decay)/(h*h);e00=(1-decay*decay)/(2*h)
    integral=dot(a,a)+dot(a,b)+dot(b,b)/3+2*dot(a,c)*e0+2*dot(b,c)*e1+dot(c,c)*e00
    loss=2*p[2]*h*integral;assert loss>=D('-1e-45')
    avg=average(f0,f1,x0,dt,p);u0=stored(f0,x0,p);u1=stored(f1,x1,p);work=dot(avg,sub(f1,f0))
    assert abs(work-(u1-u0+loss))<D('1e-28')
    derivative=None
    if tangent:
        epsilon=D('1e-18');derivative=[]
        for column in range(9):
            plus=f1.copy();minus=f1.copy();plus[column]+=epsilon;minus[column]-=epsilon
            derivative.append(scale(sub(average(f0,plus,x0,dt,p),average(f0,minus,x0,dt,p)),1/(2*epsilon)))
    return dict(green=green(f1),viscous=z1,total_dissipation=loss0+loss,physical_dissipation=loss,
        stored_old=u0,stored_new=u1,work=work,average_pk1=avg,endpoint_tangent=derivative,minimum_determinant=min_j)

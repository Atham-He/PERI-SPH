"""Independent scalar exact-arithmetic and continuum assessor for S1.

Does not import the native controller or use an iterative dual QP. The scalar
active/inactive solution is computed with Fraction; actual physical M is used
in all momentum and energy checks. Force is an interval impulse/dt, compared
with the analytic force at interval midpoint. Velocity is an endpoint value.
"""
from fractions import Fraction as Q
import math


def assess(raw, model, dt_value, tolerance=1e-10):
    m, k, q, v, d, total = map(Q, model)
    dt = Q(str(dt_value))
    assert m > 0 and k >= 0 and total > 0 and dt > 0
    # For the frozen quintic, k*w+m*w'' = D*s*[k*s²*(10-15s+6s²)
    # +60*m/T²*(1-s)*(1-2s)]. On s<=1/2 both terms are >=0.
    # On s>=1/2 w>=D/2. Differentiation gives extrema at (3±sqrt(3))/6,
    # so |w''|<=10*sqrt(3)*D/(3T²)<6D/T². This conservative bound
    # proves nonnegative reaction for this entire prescribed interval.
    assert k/2 > 6*m/total**2, 'maintained-contact reference lacks a positivity proof'
    def wall(t):
        s=t/total
        return d*(10*s**3-15*s**4+6*s**5)
    def speed(t):
        s=t/total
        return d/total*(30*s**2-60*s**3+30*s**4)
    def reaction(t):
        s=t/total
        return k*wall(t)+m*d/total**2*(60*s-180*s**2+120*s**3)
    assert raw['completed'] is True and len(raw['steps'])==int(total/dt)
    worst=0.0; force_error=[]; force_ref=[]; speed_error=[]; speed_ref=[]
    works=[]; losses=[]; max_penetration=0.0
    for index,row in enumerate(raw['steps'],1):
        t0=(index-1)*dt; t1=index*dt
        a0=wall(t0); a1=wall(t1)
        metric=m+dt**2*k/2
        free=(m*v-dt*k*q)/metric
        bound=(a1-q)/dt
        impulse=max(Q(0), metric*(bound-free))
        vn=free+impulse/metric; qn=q+dt*vn
        kinetic0=m*v*v/2; kinetic1=m*vn*vn/2
        elastic0=k*q*q/2; elastic1=k*qn*qn/2
        work=impulse*(a1-a0)/dt
        loss=impulse*(q-a0)/dt+m*(vn-v)**2/2
        assert impulse>=0 and qn>=a1 and loss>=0
        assert impulse*(vn-bound)==0
        assert m*(vn-v)+dt*k*(q+qn)/2-impulse==0
        assert work-(kinetic1-kinetic0+elastic1-elastic0)-loss==0
        exact={'step':index,'time_old':t0,'time_new':t1,'dt':dt,
               'q_old':q,'v_old':v,'wall_old':a0,'wall_new':a1,
               'q_new':qn,'v_new':vn,'impulse':impulse,'reaction_interval':impulse/dt,
               'elastic_force':-k*(q+qn)/2,'kinetic_old':kinetic0,'kinetic_new':kinetic1,
               'elastic_old':elastic0,'elastic_new':elastic1,'wall_work':work,
               'algorithmic_loss':loss,'energy_residual':0,'maximum_momentum_residual':0,
               'primal_velocity_error':0,'natural_velocity_error':0,
               'stationarity_impulse_error':0,'complementarity_energy_error':0}
        for name,value in exact.items():
            observed=row[name]
            assert math.isfinite(observed), (index,name,'nonfinite')
            # Fixed physical scales, not the effective metric and not a zero net-force denominator.
            scale=1 if name=='step' else (
                float(total) if name.startswith('time') or name=='dt' else
                float(d) if name.startswith('q_') or name.startswith('wall_') and name!='wall_work' else
                float(d/total) if name.startswith('v_') or 'velocity_error' in name else
                float(m*d/total) if 'impulse' in name or 'momentum' in name else
                float(k*d) if name in ('reaction_interval','elastic_force') else
                float(k*d*d+m*(d/total)**2))
            error=abs(observed-float(value))/scale
            worst=max(worst,error)
            assert error<=tolerance,(index,name,error)
        assert row['sweeps']>0
        pen=max(0,row['wall_new']-row['q_new'])/float(d)
        max_penetration=max(max_penetration,pen)
        assert pen<=tolerance
        fr=float(reaction((t0+t1)/2)); vr=float(speed(t1))
        force_ref.append(fr); force_error.append(row['reaction_interval']-fr)
        speed_ref.append(vr); speed_error.append(row['v_new']-vr)
        works.append(max(0,row['wall_work'])); losses.append(row['algorithmic_loss'])
        q,v=qn,vn
    relative=lambda errors,reference: math.sqrt(math.fsum(x*x for x in errors)/math.fsum(x*x for x in reference))
    return dict(steps=len(raw['steps']),dt=float(dt),maximum_normalized_exact_error=worst,
                maximum_normalized_penetration=max_penetration,
                force_midpoint_relative_L2=relative(force_error,force_ref),
                velocity_endpoint_relative_L2=relative(speed_error,speed_ref),
                algorithmic_loss=math.fsum(losses),positive_wall_work=math.fsum(works),
                loss_fraction=math.fsum(losses)/math.fsum(works),
                terminal_q=raw['steps'][-1]['q_new'],terminal_v=raw['steps'][-1]['v_new'])

"""Independent95-digit witness for the one-step Newton-limit negative control.

This does not call the native solver. It reconstructs the single recorded full
Newton update from independent path-integral residual differences, then proves
the uncommitted next residual still exceeds the original stopping tolerance.
"""
from decimal import Decimal as D,localcontext
from peri_tetra_constraints_audit import gradient,mv,transpose_mv,dvalues
from peri_gel_affine_audit import average,green,sub


def reconstruct_limit_witness(case,raw):
    with localcontext() as context:
        context.prec=95
        assert case['max_newton']==1 and raw['status']=='Newton_iteration_limit'
        assert not raw['converged'] and raw['old_history_unchanged'] and 'velocity' not in raw
        assert len(raw['trace'])==1
        entry=raw['trace'][0]
        assert entry['mask']==0 and entry['iteration']==0 and entry['alpha']==1 and entry['domain_rejections']==0
        assert entry['linear_residual']<=1e-12
        b=gradient();q=dvalues(case['position']);v0=dvalues(case['velocity']);m=dvalues(case['node_mass'])
        dt=D(str(case['dt']));volume=D(str(case['volume']));p=dvalues(case['parameters']);z=dvalues(case['viscous'])
        f0=mv(b,q);x0=sub(green(f0),z);scale=sum(m)*D(str(case['velocity_scale']))
        def residual(v):
            f1=mv(b,[a+dt*u for a,u in zip(q,v)])
            force=transpose_mv(b,average(f0,f1,x0,dt,p))
            return [m[i//3]*(v[i]-v0[i])+dt*volume*force[i] for i in range(12)]
        cold=[D(0)]*12;r0=residual(cold);epsilon=D('1e-18');jac=[[D(0)]*12 for _ in range(12)]
        for j in range(12):
            plus=list(cold);minus=list(cold);plus[j]+=epsilon;minus[j]-=epsilon
            a=residual(plus);c=residual(minus)
            for i in range(12):jac[i][j]=(a[i]-c[i])/(2*epsilon)
        augmented=[row+[-value] for row,value in zip(jac,r0)]
        for k in range(12):
            pivot=max(range(k,12),key=lambda i:abs(augmented[i][k]))
            assert abs(augmented[pivot][k])>D('1e-30')
            augmented[k],augmented[pivot]=augmented[pivot],augmented[k]
            divisor=augmented[k][k];augmented[k]=[v/divisor for v in augmented[k]]
            for i in range(12):
                if i!=k:
                    multiplier=augmented[i][k];augmented[i]=[a-multiplier*c for a,c in zip(augmented[i],augmented[k])]
        step=[row[-1] for row in augmented]
        linear=max(abs(sum(jac[i][j]*step[j] for j in range(12))+r0[i]) for i in range(12))/scale
        before=max(map(abs,r0))/scale;after=max(map(abs,residual(step)))/scale
        tolerance=D(str(case['tolerance']))
        assert linear<D('1e-60')
        assert abs(float(before)-entry['merit'])<1e-10
        assert before>tolerance and after>tolerance and after<before*(1-D('1e-4'))
        return dict(status='independent_one_step_limit_witness_proved',precision_digits=95,
            normalized_initial_merit=str(before),normalized_uncommitted_next_merit=str(after),
            original_stopping_tolerance=str(tolerance),independent_linear_residual=str(linear),
            reconstructed_uncommitted_velocity=list(map(str,step)),native_solver_called=False)

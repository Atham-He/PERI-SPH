"""Independent95-digit residual/Jacobian/energy certificate for global controls.

The candidate root is checked as an unknown in independent equations, not fed
back as a prescribed reference state. Global derivatives use high-precision
differences of independent path integrals, not native local/assembled tangents.
"""
from decimal import Decimal as D,localcontext
import math
from peri_gel_affine_audit import average,reference_step,green,sub,dot,minimum_det
from peri_dense_contact_audit import exact as exact_quadratic


def gradient():
    gradients=[[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1]]
    return [[D(gradients[a][j] if i==k else 0) for a in range(4) for k in range(3)] for i in range(3) for j in range(3)]
def mv(matrix,vector): return [dot(row,vector) for row in matrix]
def transpose_mv(matrix,vector): return [sum((row[i]*v for row,v in zip(matrix,vector)),D(0)) for i in range(len(matrix[0]))]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def add(a,b): return [x+y for x,y in zip(a,b)]
def dvalues(values): return list(map(lambda v:D(str(v)),values))


def inspect(case,native,code):
    with localcontext() as context:
        context.prec=95
        assert code==case.get('expected_code',0),(case['name'],'return',code)
        if code==1:
            assert native.get('error') and case.get('invalid_reason')
            if case['invalid_reason']=='nonpositive_mass': assert any(m<=0 for m in case['node_mass'])
            elif case['invalid_reason']=='zero_constraint_row': assert any(not any(row['C']) for row in case['rows'])
            elif case['invalid_reason']=='invalid_trial_path':
                b=gradient(); q=dvalues(case['position']);v=dvalues(case['sample_velocity']);dt=D(str(case['dt']))
                assert minimum_det(mv(b,q),mv(b,[x+dt*u for x,u in zip(q,v)]))<D('.05')
            else: raise AssertionError('unknown independent invalidity proof')
            return dict(name=case['name'],status='declared_invalid_configuration_rejected',reason=case['invalid_reason'])
        if code==2:
            assert not native['converged'] and native['status']==case['expected_status']
            assert native['old_history_unchanged'] and 'velocity' not in native
            if case['expected_status']=='no_KKT_feasible_active_set':
                rows=case['rows']; a,b=rows[:2]
                assert all(D(str(x))==-D(str(y)) for x,y in zip(a['C'],b['C']))
                assert D(str(a['new']))+D(str(b['new']))>0
            else:
                from peri_tetra_limit_witness import reconstruct_limit_witness
                witness=reconstruct_limit_witness(case,native)
                return dict(name=case['name'],status='declared_terminal_rejection_without_commit',reason=native['status'],independent_limit_witness=witness)
            return dict(name=case['name'],status='declared_terminal_rejection_without_commit',reason=native['status'])
        if case['mode']==1:
            qcase=dict(name=case['name'],mode=1,matrix=case['stiffness'],mass=case['mass'],position=case['position'],
                velocity=case['velocity'],dt=case['dt'],rows=case['rows'])
            ref=exact_quadratic(qcase); assert ref is not None and native['converged']
            error=max(abs(v-float(w)) for v,w in zip(native['velocity'],ref['velocity']))
            assert error<=1e-10
            for l,w in zip(native['normal_impulse'],ref['normal_impulse']): assert abs(l-float(w))<=1e-10
            return dict(name=case['name'],status='independent_quadratic_limit_root_passed',maximum_normalized_error=error)
        b=gradient(); q0=dvalues(case['position']); oldv=dvalues(case['velocity']); mass=dvalues(case['node_mass'])
        dt=D(str(case['dt'])); volume=D(str(case['volume'])); p=dvalues(case['parameters'])
        z0=dvalues(case['viscous']); oldloss=D(str(case['old_dissipation'])); f0=mv(b,q0)
        assert all(abs(a-bb)<D('1e-12') for a,bb in zip(dvalues(case['green']),green(f0)))
        velocity=dvalues(native['velocity']); assert len(velocity)==12
        position=[q+dt*v for q,v in zip(q0,velocity)]; f1=mv(b,position)
        material=reference_step(f0,f1,z0,oldloss,dt,p,False)
        force=[-volume*v for v in transpose_mv(b,material['average_pk1'])]
        momentum=[mass[i//3]*(velocity[i]-oldv[i])-dt*force[i] for i in range(12)]
        x0=sub(green(f0),z0)
        def residual(v):
            end=mv(b,[q+dt*u for q,u in zip(q0,v)])
            stress=average(f0,end,x0,dt,p)
            internal=transpose_mv(b,stress)
            return [mass[i//3]*(v[i]-oldv[i])+dt*volume*internal[i] for i in range(12)]
        epsilon=D('1e-18'); jacobian=[[D(0)]*12 for _ in range(12)]
        for j in range(12):
            plus=list(velocity);minus=list(velocity);plus[j]+=epsilon;minus[j]-=epsilon
            a=residual(plus);c=residual(minus)
            for i in range(12): jacobian[i][j]=(a[i]-c[i])/(2*epsilon)
        norm=sum((v*v for row in jacobian for v in row),D(0))
        skew=sum(((jacobian[i][j]-jacobian[j][i])**2 for i in range(12) for j in range(12)),D(0))
        asym=(skew/norm).sqrt() if norm else D(0)
        total_mass=sum(mass); speed=D(str(case['velocity_scale'])); pscale=total_mass*speed
        escale=max(total_mass*speed*speed,volume*max(p)); worst=0.0
        def compare(field,wanted,scale_value):
            nonlocal worst
            got=native[field]
            if wanted and isinstance(wanted,list) and isinstance(wanted[0],list):
                got=[x for row in got for x in row]; wanted=[x for row in wanted for x in row]
            elif not isinstance(wanted,list): got=[got];wanted=[wanted]
            assert len(got)==len(wanted),(case['name'],field,'dimensions')
            for observed,expected in zip(got,wanted):
                error=abs(observed-float(expected))/float(scale_value)
                assert math.isfinite(error) and error<=1e-10,(case['name'],field,error)
                worst=max(worst,error)
        compare('B',b,D(1));compare('F1',f1,D(1));compare('internal_force',force,pscale/dt)
        compare('momentum_equation',momentum,pscale);compare('jacobian',jacobian,total_mass)
        compare('jacobian_asymmetry',asym,D(1))
        compare('average_pk1',material['average_pk1'],max(p));compare('green',material['green'],D(1));compare('viscous',material['viscous'],D(1))
        for field,source in [('material_stored_old','stored_old'),('material_stored_new','stored_new'),('material_loss','physical_dissipation'),
                             ('total_dissipation','total_dissipation'),('material_work','work')]: compare(field,material[source],max(p))
        compare('transpose_work_residual',-dt*dot(force,velocity)-volume*material['work'],escale)
        # Every row of B kills a rigid translation, and the reference tetra is reproduced.
        for row in b:
            for k in range(3): assert sum((row[3*a+k] for a in range(4)),D(0))==0
        if case['mode']==2:
            assert native['operator_probe_only'] and not native['converged']
            return dict(name=case['name'],status='independent_global_operator_and_tangent_passed',maximum_normalized_error=worst,
                        actual_jacobian_asymmetry=float(asym),accepted_physical_step=False)
        assert native['converged'] and native['old_history_unchanged']
        impulse=dvalues(native['normal_impulse']); assert len(impulse)==len(case['rows']) and all(l>=0 for l in impulse)
        scatter=[D(0)]*12; work=D(0); normal_loss=D(0); complementarity=D(0); min_gap=None
        for row,l in zip(case['rows'],impulse):
            c=dvalues(row['C']); old=D(str(row['old']));new=D(str(row['new'])); bound=(new-dot(c,q0))/dt
            gap=dot(c,velocity)-bound; min_gap=gap if min_gap is None else min(min_gap,gap)
            assert gap>=-D('1e-10')*speed
            complementarity+=abs(l*gap)
            scatter=add(scatter,[l*v for v in c]);work+=l*(new-old)/dt;normal_loss+=l*(dot(c,q0)-old)/dt
        assert complementarity<=D('1e-10')*pscale*speed
        equation=sub(momentum,scatter)
        assert max(map(abs,equation))<=D('1e-10')*pscale,(case['name'],'independent_momentum')
        ke0=sum((mass[i//3]*oldv[i]**2/2 for i in range(12)),D(0))
        ke1=sum((mass[i//3]*velocity[i]**2/2 for i in range(12)),D(0))
        algorithm=sum((mass[i//3]*(velocity[i]-oldv[i])**2/2 for i in range(12)),D(0))
        energy=work-(ke1-ke0+volume*(material['stored_new']-material['stored_old']+material['physical_dissipation']))-algorithm-normal_loss
        assert abs(energy)<=D('1e-10')*escale,(case['name'],'independent_energy')
        for field,value in [('kinetic_old',ke0),('kinetic_new',ke1),('algorithmic_loss',algorithm),('normal_contact_loss',normal_loss),
                            ('wall_work',work),('energy_residual',energy)]: compare(field,value,escale)
        compare('position',position,D(1));compare('contact_impulse',scatter,pscale);compare('momentum_residual',equation,pscale)
        totalp=[D(0)]*3; forcesum=[D(0)]*3; torque=[D(0)]*3; angular=[D(0)]*3
        for a in range(4):
            s=slice(3*a,3*a+3);dq=sub(velocity[s],oldv[s]);qnext=position[s]
            totalp=add(totalp,sub([mass[a]*u for u in dq],scatter[s]));forcesum=add(forcesum,force[s])
            torque=add(torque,cross(q0[s],force[s]))
            angular=add(angular,sub([mass[a]*v for v in sub(cross(qnext,velocity[s]),cross(q0[s],oldv[s]))],cross(q0[s],scatter[s])))
        balance=sub(angular,[dt*v for v in torque])
        for field,value,scale_value in [('total_momentum_residual',totalp,pscale),('internal_force_sum',forcesum,pscale/dt),
            ('internal_torque_old_cursor',torque,pscale/dt),('angular_change_minus_external',angular,pscale),('angular_balance_old_cursor',balance,pscale)]: compare(field,value,scale_value)
        assert max(map(abs,totalp))<=D('1e-10')*pscale
        assert max(map(abs,balance))<=D('1e-10')*pscale
        assert native['normalized_KKT']<=1e-10 and native['maximum_linear_residual']<=1e-12
        return dict(name=case['name'],status='independent_nonlinear_global_momentum_energy_Jacobian_KKT_passed',
            maximum_normalized_error=worst,actual_jacobian_asymmetry=float(asym),minimum_velocity_constraint_gap=float(min_gap) if min_gap is not None else None,
            angular_drift_minus_external=[float(x) for x in angular],internal_old_cursor_torque=[float(x) for x in torque],
            angular_conservation_qualified=False,one_step_algorithmic_loss=float(algorithm),positive_wall_work=float(max(D(0),work)),
            trajectory_accuracy_qualified=False)

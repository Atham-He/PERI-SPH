"""Exact-rational full-matrix oracle. Enumerates active sets, not native sweeps.

SPD/PSD is decided by exact principal determinants, not native Cholesky or
numerical Schur tolerances. Physical ledgers always use actual diagonal M.
"""
from fractions import Fraction as Q
from itertools import combinations
import math
from peri_normal_projection_audit import solve, dot


def determinant(matrix):
    a=[list(row) for row in matrix]; result=Q(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None: return Q(0)
        if pivot!=i: a[i],a[pivot]=a[pivot],a[i]; result=-result
        value=a[i][i]; result*=value
        for j in range(i+1,len(a)):
            multiplier=a[j][i]/value
            for k in range(i+1,len(a)): a[j][k]-=multiplier*a[i][k]
    return result


def exact(case):
    matrix=[[Q(str(v)) for v in row] for row in case['matrix']]
    mass,q,v=[[Q(str(value)) for value in case[key]] for key in ['mass','position','velocity']]
    n=len(q); dt=Q(str(case['dt']))
    if not n or n>8 or len(case['rows'])>4: raise ValueError('oracle budget')
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(n)): raise ValueError('nonsymmetric')
    if case['mode']==1:
        if any(m<=0 for m in mass): raise ValueError('nonpositive physical mass')
        for count in range(1,n+1):
            for subset in combinations(range(n),count):
                if determinant([[matrix[i][j] for j in subset] for i in subset])<0: raise ValueError('non-PSD stiffness')
        metric=[[dt*dt*matrix[i][j]/2+(mass[i] if i==j else 0) for j in range(n)] for i in range(n)]
        free=solve(metric,[mass[i]*v[i]-dt*dot(matrix[i],q) for i in range(n)])
    else: metric=matrix; free=v
    if any(determinant([row[:count] for row in metric[:count]])<=0 for count in range(1,n+1)):
        raise ValueError('non-SPD metric')
    C=[[Q(str(v)) for v in row['C']] for row in case['rows']]
    gaps=[dot(c,q)-Q(str(row['old'])) for c,row in zip(C,case['rows'])]
    if any(g<0 for g in gaps): raise ValueError('initial penetration')
    if any(not any(row) for row in C): raise ValueError('zero constraint row')
    bound=[(Q(str(row['new']))-dot(c,q))/dt for c,row in zip(C,case['rows'])]
    columns=[solve(metric,row) for row in C]
    for count in range(min(len(C),n)+1):
        for active in combinations(range(len(C)),count):
            gram=[[dot(C[i],columns[j]) for j in active] for i in active]
            multipliers=solve(gram,[bound[i]-dot(C[i],free) for i in active])
            if multipliers is None or any(v<0 for v in multipliers): continue
            impulse=[Q(0)]*len(C)
            for i,value in zip(active,multipliers): impulse[i]=value
            velocity=[free[i]+sum(columns[s][i]*impulse[s] for s in range(len(C))) for i in range(n)]
            if any(dot(c,velocity)<b for c,b in zip(C,bound)): continue
            result=dict(velocity=velocity,normal_impulse=impulse,C=C,bound=bound,metric=metric,free=free)
            if case['mode']==1:
                position=[old+dt*vn for old,vn in zip(q,velocity)]
                force=[-dot(row,[(a+b)/2 for a,b in zip(q,position)]) for row in matrix]
                work=sum((l*(Q(str(row['new']))-Q(str(row['old'])))/dt for l,row in zip(impulse,case['rows'])),Q(0))
                ke0=sum(m*u*u/2 for m,u in zip(mass,v)); ke1=sum(m*u*u/2 for m,u in zip(mass,velocity))
                pe0=dot(q,[dot(row,q) for row in matrix])/2
                pe1=dot(position,[dot(row,position) for row in matrix])/2
                loss=sum(m*(vn-vo)**2/2 for m,vn,vo in zip(mass,velocity,v))
                normal_loss=sum((l*g/dt for l,g in zip(impulse,gaps)),Q(0))
                assert work-(ke1-ke0+pe1-pe0)-loss-normal_loss==0
                result.update(position=position,internal_force=force,kinetic_old=ke0,kinetic_new=ke1,
                    elastic_old=pe0,elastic_new=pe1,algorithmic_loss=loss,normal_contact_loss=normal_loss,wall_work=work)
            return result
    return None


def check(case,native,code):
    expected=case.get('expected_code',0)
    assert code==expected,(case['name'],'return',code,expected)
    try: reference=exact(case)
    except ValueError as error:
        assert code==1 and native.get('error'),(case['name'],'invalid control')
        return dict(name=case['name'],status='invalid_operator_or_geometry_rejected',reason=str(error))
    if reference is None:
        assert code==2 and not native['converged']
        return dict(name=case['name'],status='infeasible_control_rejected',sweeps=native['sweeps'])
    assert code==0 and native['converged']
    V=case.get('velocity_scale',1); dt=case['dt']
    scale_p=sum(float(reference['metric'][i][i]) for i in range(len(case['mass'])))*V
    scale_e=scale_p*V; scale_q=max(dt*V,max(map(abs,case['position'])))
    errors=[]
    def compare(field,scale):
        wanted=reference[field]; got=native[field]
        if isinstance(wanted,list):
            assert len(wanted)==len(got)
            values=[abs(a-float(b))/scale for a,b in zip(got,wanted)]
        else: values=[abs(got-float(wanted))/scale]
        assert all(math.isfinite(x) and x<=1e-10 for x in values),(case['name'],field,values)
        errors.extend(values)
    compare('velocity',V)
    if case['mode']==1:
        compare('position',scale_q); compare('internal_force',scale_p/dt)
        for field in ['kinetic_old','kinetic_new','elastic_old','elastic_new','algorithmic_loss','normal_contact_loss','wall_work']:
            compare(field,scale_e)
        assert abs(native['energy_residual'])<=1e-10*scale_e
        assert max(map(abs,native['momentum_residual']))<=1e-10*scale_p
    C=[[float(value) for value in row] for row in reference['C']]
    impulse=native['normal_impulse']; velocity=native['velocity']
    assert len(impulse)==len(C) and all(math.isfinite(l) and l>=0 for l in impulse)
    stationarity=[]
    for i,row in enumerate(reference['metric']):
        value=math.fsum(float(a)*(b-float(f)) for a,b,f in zip(row,velocity,reference['free']))
        value-=math.fsum(c[i]*l for c,l in zip(C,impulse))
        stationarity.append(abs(value)/scale_p)
    gaps=[math.fsum(c*v for c,v in zip(row,velocity))-float(b) for row,b in zip(C,reference['bound'])]
    primal=max([0]+[-value/V for value in gaps])
    comp=math.fsum(abs(l*g) for l,g in zip(impulse,gaps))/scale_e
    assert max(stationarity+[primal,comp])<=1e-10,(case['name'],'KKT',stationarity,primal,comp)
    assert native['linear_relative_residual']<=1e-10
    errors.extend(stationarity+[primal,comp,native['linear_relative_residual']])
    return dict(name=case['name'],status='exact_full_matrix_primal_physical_energy_and_dual_KKT_passed',
                maximum_normalized_error=max(errors),DOF=len(velocity),rows=len(C),sweeps=native['sweeps'])

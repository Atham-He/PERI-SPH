"""Independent original-equation assessment of bounded simultaneous controls."""
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import math
from peri_coulomb_contact_audit import (
    frac,decimal,fdot,coefficients,operator_reference,original_reference,check_original,check_trace,close_tree,dec)


def analytic_particle(p):
    """A physical diagonal unit-mass impact; normalization mass is unrelated."""
    assert p['n']==3 and p['m']==1 and p['beta']==[0.]*3
    assert p['matrix']==[[float(i==j) for j in range(3)] for i in range(3)]
    c=p['contacts'][0];basis=[c['normal'],*c['tangent']];free=list(map(dec,p['load']))
    free_normal=sum((dec(a)*b for a,b in zip(basis[0],free)),D(0))
    boundary=(dec(c['new_offset'])-decimal(fdot(basis[0],p['old_position'])))/dec(p['dt'])
    normal=max(free_normal,boundary);l=normal-free_normal
    wall=list(map(dec,c['wall_tangent_velocity']))
    u=[sum((dec(a)*b for a,b in zip(t,free)),D(0))-w for t,w in zip(basis[1:],wall)]
    length=sum(x*x for x in u).sqrt();radius=dec(c['mu'])*l
    tau=[-x*min(D(1),radius/length) for x in u] if length else [D(0),D(0)]
    velocity=[free[i]+l*dec(basis[0][i])+sum((tau[k]*dec(basis[k+1][i]) for k in range(2)),D(0)) for i in range(3)]
    relative=[u[k]+tau[k] for k in range(2)]
    normal_loss=l*(normal-boundary);friction=-sum((t*w for t,w in zip(tau,relative)),D(0))
    wall_work=l*boundary+sum((t*w for t,w in zip(tau,wall)),D(0))
    kinetic_old=sum(x*x for x in free)/2;kinetic_new=sum(x*x for x in velocity)/2
    algorithmic=sum((a-b)**2 for a,b in zip(velocity,free))/2
    assert abs(kinetic_new-kinetic_old-wall_work+normal_loss+friction+algorithmic)<=D('1e-85')
    return velocity,[l,*tau],dict(kinetic_old=kinetic_old,kinetic_new=kinetic_new,
        wall_work=wall_work,normal_loss=normal_loss,friction_loss=friction,algorithmic_loss=algorithmic)


def check_header(p,h):
    assert h['type']=='header' and h['format']==1 and not h['actual_PD_body'] and h['manufactured_global_control']
    assert h['accepted_physical_steps']==0 and not h['S3_complete']
    for key in ['mode','fault','old_position','dt','mass_scale','velocity_scale','domain_limit','max_newton_steps','max_line_search',
                'contacts','trial_velocity','trial_impulses','linear_rhs']:
        assert h[key]==p[key],('input_header',key)
    for target,source in [('primal_dofs','n'),('contact_count','m'),('operator_matrix','matrix'),('operator_load','load'),('operator_cubic_beta','beta')]:
        assert h[target]==p[source],('input_header',target)
    assert h['impulse_scale']==p['mass_scale']*p['velocity_scale'] and h['relative_tolerance']==1e-12
    assert h['unknowns']==p['n']+3*p['m']
    c=coefficients(p['contacts'],p['n']);assert h['C_full_three_components']==c
    return c


def inspect(case,records,code):
    with localcontext() as ctx:
        ctx.prec=95
        p=case['input'];kind=case['kind'];n=p['n'];m=p['m'];size=n+3*m
        assert code==case['expected_code'],(case['name'],'return_code',code)
        assert records[0]['type']=='identity'
        if kind=='invalid':
            assert len(records)==2 and records[-1]['type']=='error'
            assert records[-1]['message']==case['expected_message'] and not records[-1]['accepted_physical_step']
            return dict(name=case['name'],status='declared_invalid_input_rejected',accepted_physical_steps=0)
        assert len(records)==4 and records[-1]['type']=='end'
        h,row,end=records[1:];c=check_header(p,h)
        assert end['accepted_physical_steps']==end['last_accepted_time']==0 and not end['full_S3_complete']
        assert end['completed']==(code==0) and not row['accepted_physical_step']
        I=h['impulse_scale'];V=h['velocity_scale'];maximum=0.
        material=lambda v,_:operator_reference(h,v)
        result=dict(name=case['name'],status='original_global_equations_independently_passed',accepted_physical_steps=0)
        if kind in ['sample','linear']:
            assert row['type']=='global_sample';R,J=material(p['trial_velocity'],None)
            reference,maximum=check_original(p['old_position'],p['contacts'],p['dt'],I,V,
                p['trial_velocity'],p['trial_impulses'],row['original'],R,J,c)
            if kind=='linear':
                assert row['linear_consistent'] and len(row['linear_solution'])==size
                assert row['linear_rank']==case['exact_rank'] and row['linear_backward_residual']<=1e-12
                a=row['original']['jacobian'];x=row['linear_solution'];rhs=p['linear_rhs']
                error=max(abs(fdot(r,x)-frac(b)) for r,b in zip(a,rhs))
                scale=max(sum(map(abs,map(frac,r)),F(0)) for r in a)*max(map(abs,map(frac,x)))+max(map(abs,map(frac,rhs)))
                assert error/scale<=frac(1e-12)
                maximum=max(maximum,close_tree(x[:n],case['known_primal']))
                # Repeated selector rows constrain every coordinate; full-row-rank
                # selector C has rank n, and [A,-C^T;D*C,0] has exact rank 2n.
                assert n==99 and m in [116,228] and case['exact_rank']==2*n
                assert all(sum(abs(z) for z in r)==1 for r in c)
                assert {tuple(r) for r in c}=={tuple(float(i==j) for i in range(n)) for j in range(n)}
                maximum=max(maximum,close_tree([decimal(fdot([r[i] for r in c],x[n:])) for i in range(n)],
                    [decimal(fdot([r[i] for r in c],case['known_dual'])) for i in range(n)]))
                result.update(original_matrix_dimension=size,exact_rank=case['exact_rank'],linear_backward=float(error/scale))
            if case.get('finite_difference'):
                epsilon=F(1,10**22);v=list(map(frac,p['trial_velocity']));imp=list(map(frac,p['trial_impulses']))
                error=D(0)
                for column in range(size):
                    vp,vm,pp,pm=v.copy(),v.copy(),imp.copy(),imp.copy()
                    if column<n:vp[column]+=epsilon*frac(V);vm[column]-=epsilon*frac(V)
                    else:pp[column-n]+=epsilon*frac(I);pm[column-n]-=epsilon*frac(I)
                    def value(v,p):
                        R,J=material(v,None)
                        return original_reference(h['old_position'],h['contacts'],h['dt'],I,V,v,p,R,J,c)['residual']
                    a,b=value(vp,pp),value(vm,pm)
                    error=max(error,max(abs((u-v)/D('2e-22')-reference['jacobian'][r][column]) for r,(u,v) in enumerate(zip(a,b))))
                assert error<=D('1e-10');result['joint_derivative_selfcheck']=float(error)
            if case['name']=='exact_original_gather_cancellation':
                assert reference['relative_velocity'][0][0]==1
            if case['name']=='exact_branch_inside_by_sub_95digit_gap':
                assert row['original']['disk'][0]['branch']=='stick_interior'
        else:
            assert row['type']=='global_solve';s=row['solve']
            root=kind=='root'
            maximum=check_trace(p['old_position'],p['contacts'],p['dt'],I,V,s,material,c,require_root=root)
            result['Newton_observations']=len(s['trace'])
            if root:
                if 'target_velocity' in case:
                    maximum=max(maximum,close_tree(s['velocity'],case['target_velocity']),close_tree(s['impulses'],case['target_impulses']))
                if case.get('physical_particle'):
                    v,imp,energy=analytic_particle(p)
                    maximum=max(maximum,close_tree(s['velocity'],v),close_tree(s['impulses'],imp))
                    result['independent_physical_unit_mass_energy']={key:float(value) for key,value in energy.items()}
                    fv=list(map(frac,s['velocity']));fp=list(map(frac,s['impulses']));free=list(map(frac,p['load']))
                    contact=p['contacts'][0];wall=list(map(frac,contact['wall_tangent_velocity']))
                    boundary=(frac(contact['new_offset'])-fdot(contact['normal'],p['old_position']))/frac(p['dt'])
                    kinetic_old=fdot(free,free)/2;kinetic_new=fdot(fv,fv)/2
                    relative=[fdot(t,fv)-w for t,w in zip(contact['tangent'],wall)]
                    W=fp[0]*boundary+fdot(fp[1:],wall)
                    Dnormal=fp[0]*(fdot(contact['normal'],fv)-boundary)
                    Dfriction=-fdot(fp[1:],relative);Dalgorithm=fdot([a-b for a,b in zip(fv,free)],[a-b for a,b in zip(fv,free)])/2
                    balance=abs(kinetic_new-kinetic_old-W+Dnormal+Dfriction+Dalgorithm)/max(F(1),abs(kinetic_old),abs(W))
                    assert balance<=frac(1e-10) and min(Dnormal,Dfriction,Dalgorithm)>=-frac(1e-10)
                    maximum=max(maximum,float(balance))
                    result['native_root_true_mass_energy']={key:float(x) for key,x in dict(
                        kinetic_old=kinetic_old,kinetic_new=kinetic_new,wall_work=W,normal_loss=Dnormal,
                        friction_loss=Dfriction,algorithmic_loss=Dalgorithm,normalized_balance=balance).items()}
            else:
                assert not s['converged'] and s['velocity']==s['impulses']==[]
                if kind=='precision':
                    w=s['precision_rejection'];assert w['present'] and not w['during_line_search']
                    assert s['status']=='original_contact_arithmetic_precision_budget' and not w['original']['precision_budget_ok']
                    r=w['original'];R,J=material(w['velocity'],None)
                    assert max(r['residual_roundoff_bound']+[r['maximum_jacobian_roundoff_bound']])>1e-12
                    assert abs(decimal(frac(R[0])/frac(I)))>D('1e140')
                    original=[(frac(a)-fdot([row[k] for row in c],w['impulses']))/frac(I) for k,a in enumerate(r['supplied_momentum'])]
                    assert any(abs(frac(float(x))-x)>frac(1e-12) for x in original)
                    assert all(a==float(b) for a,b in zip(r['residual'][:n],original))
                    result['status']='declared_unsupported_precision_rejected'
                elif kind=='evaluation_failure':
                    w=s['evaluation_failure'];assert w['present'] and not s['precision_rejection']['present']
                    assert s['status']=='original_material_or_contact_evaluation_failed'
                    assert len(w['velocity'])==n and len(w['impulses'])==3*m
                    assert w['during_line_search']==(p['fault']==2)
                    assert w['reason']==('Coulomb contact: material callback row' if p['fault']==3 else 'declared frozen material evaluation fault')
                    if w['during_line_search']:
                        assert w['alpha']>0 and any(w['velocity']) and len(s['trace'])==1
                        t=s['trace'][0];alpha=dec(w['alpha']);direction=t['normalized_direction']
                        maximum=max(maximum,close_tree(w['velocity'],[dec(a)+alpha*dec(V)*dec(b) for a,b in zip(t['velocity'],direction[:n])]),
                            close_tree(w['impulses'],[dec(a)+alpha*dec(I)*dec(b) for a,b in zip(t['impulses'],direction[n:])]))
                    result['status']='declared_callback_failure_preserved_uncommitted'
                else:
                    assert kind=='infeasible'
                    # Exact Farkas certificate: Cn0+Cn1=0 but b0+b1>0.
                    assert all(frac(a)+frac(b)==0 for a,b in zip(c[0],c[3]))
                    bounds=[(frac(h['contacts'][k]['new_offset'])-fdot(c[3*k],p['old_position']))/frac(p['dt']) for k in range(2)]
                    assert sum(bounds)>0
                    assert s['status'] in ['original_coupled_linear_system_rejected','line_search_limit',
                                           'Newton_iteration_limit','original_contact_arithmetic_precision_budget']
                    result.update(status='exact_infeasibility_certificate_no_accepted_root',native_rejection=s['status'])
        result.update(maximum_normalized_error=maximum,full_S3_complete=False,actual_PD_friction_qualified=False)
        return result

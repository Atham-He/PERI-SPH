"""Prospective algebra inputs; no native execution, no physical trajectory claim."""
import copy,math
from fractions import Fraction as F
from peri_coulomb_contact_audit import coefficients,fdot


def plane(weights=None,mu=.3):
    return dict(weights=[1.] if weights is None else weights,normal=[1.,0.,0.],tangent=[[0.,1.,0.],[0.,0.,1.]],
                old_offset=0.,new_offset=0.,mu=mu,wall_tangent_velocity=[0.,0.])


def problem(n=3,contacts=None):
    contacts=[plane()] if contacts is None else contacts
    return dict(mode=0,n=n,m=len(contacts),dt=1.,mass_scale=1.,velocity_scale=1.,domain_limit=1e6,
        max_newton_steps=48,max_line_search=32,fault=0,old_position=[0.]*n,
        matrix=[[float(i==j) for j in range(n)] for i in range(n)],load=[-1.,2.,1.]+[0.]*(n-3),
        beta=[0.]*n,contacts=contacts,trial_velocity=[0.]*n,trial_impulses=[0.]*(3*len(contacts)),
        linear_rhs=[0.]*(n+3*len(contacts)))


def tokens(p):
    v=[p[k] for k in ['mode','n','m','dt','mass_scale','velocity_scale','domain_limit','max_newton_steps','max_line_search','fault']]
    v+=p['old_position']+[x for row in p['matrix'] for x in row]+p['load']+p['beta']
    for c in p['contacts']:
        v+=c['weights']+c['normal']+[x for row in c['tangent'] for x in row]
        v+=[c[k] for k in ['old_offset','new_offset','mu']]+c['wall_tangent_velocity']
    v+=p['trial_velocity']+p['trial_impulses']+p['linear_rhs']
    return [repr(x) if isinstance(x,float) else str(x) for x in v]


def controls():
    cases=[]
    def add(name,p=None,code=0,kind='root',**kw):
        p=problem() if p is None else p
        cases.append(dict(name=name,input=p,tokens=tokens(p),expected_code=code,kind=kind,**kw))
    add('physical_slide_two_components',physical_particle=True)
    p=problem();p['load']=[-1.,-2.,-1.];add('physical_reverse_slide',p,physical_particle=True)
    p=problem(contacts=[plane(mu=3.)]);add('physical_strict_stick',p,physical_particle=True)
    p=problem(contacts=[plane(mu=0.)]);add('physical_zero_mu',p,physical_particle=True)
    p=problem();p['old_position'][0]=.1;p['load']=[1.,2.,1.];add('physical_separated',p,physical_particle=True)
    p=problem();p['contacts'][0]['new_offset']=.1;add('physical_moving_normal_wall',p,physical_particle=True)
    p=problem();p['contacts'][0]['wall_tangent_velocity']=[.4,-.2];p['load']=[-1.,.4,-.2]
    add('physical_common_wall_velocity',p,physical_particle=True)
    p=problem();p['mass_scale']=7.;p['velocity_scale']=2.5;add('normalization_preserves_physical_mass',p,physical_particle=True)
    p=problem();p['contacts'][0]['tangent']=[[0.,0.,1.],[0.,-1.,0.]]
    add('tangent_basis_covariance',p,physical_particle=True,comparison='physical_slide_two_components')
    p=problem();p['contacts'][0]['normal']=[0.,0.,1.];p['contacts'][0]['tangent']=[[1.,0.,0.],[0.,1.,0.]]
    p['load']=[2.,1.,-1.];add('coordinate_frame_covariance',p,physical_particle=True,comparison='physical_slide_two_components')
    p=problem(contacts=[plane(),plane()]);add('redundant_rows',p,comparison='physical_slide_two_components')
    matrix=[[2.,.2,-.1],[.05,1.4,.3],[-.12,.1,1.2]]
    target=[0.,.2,-.3];q=math.sqrt(.13);impulses=[.6,-.24*.2/q,.24*.3/q]
    for nonlinear in [False,True]:
        p=problem(contacts=[plane(mu=.4)]);p['matrix']=copy.deepcopy(matrix)
        if nonlinear:p['beta']=[.15,.25,.35]
        p['load']=[float(fdot(row,target)+F.from_float(b)*F.from_float(v)**3-F.from_float(t))
                   for row,b,v,t in zip(matrix,p['beta'],target,impulses)]
        add('nonsymmetric_'+('nonlinear' if nonlinear else 'linear')+'_root',p,target_velocity=target,target_impulses=impulses)
    p=problem(6,[plane([-.25,1.25],mu=.4)]);target=[0.,.3,.4,0.,.3,.4]
    impulses=[.7,-.168,-.224];c=coefficients(p['contacts'],6)
    p['beta']=[.1,.2,.3,.2,.3,.1]
    p['load']=[float(F.from_float(v)+F.from_float(b)*F.from_float(v)**3-fdot([row[i] for row in c],impulses))
               for i,(v,b) in enumerate(zip(target,p['beta']))]
    add('signed_two_node_nonlinear_root',p,target_velocity=target,target_impulses=impulses)
    p=problem();p['mode']=1;p['matrix']=matrix;p['beta']=[.1,.2,.3]
    p['trial_velocity']=[.1,.3,-.2];p['trial_impulses']=[.5,-.1,.04]
    add('full_joint_original_Jacobian',p,kind='sample',finite_difference=True)
    p=problem(contacts=[plane(mu=1.)]);p['mode']=1
    p['trial_velocity']=[0.,1e-100,0.];p['trial_impulses']=[1.,1.,0.]
    add('exact_branch_inside_by_sub_95digit_gap',p,kind='sample')
    p=problem(9,[plane([1e16,1.,-1e16],mu=1.)]);p['mode']=1;p['mass_scale']=1e16
    p['trial_velocity']=[0.,1.,0.]*3;p['trial_impulses']=[1e16,-1e16,0.]
    add('exact_original_gather_cancellation',p,kind='sample')
    for m in [116,228]:
        rows=[]
        for s in range(m):
            w=[float(i==s%33) for i in range(33)];rows.append(plane(w,mu=.5))
        p=problem(99,rows);p['mode']=2;p['matrix'][0][1]=.2;p['matrix'][1][0]=-.1
        p['trial_impulses']=[1.,0.,0.]*m
        # At strict stick, g=0 and lambda=1: original block is [A,-C^T;D*C,0].
        c=coefficients(rows,99);x=[(.25 if i%2 else -.125) for i in range(99)]
        dual=[(.125 if i%2 else -.25) for i in range(3*m)]
        p['linear_rhs']=[float(fdot(row,x)-fdot([r[i] for r in c],dual)) for i,row in enumerate(p['matrix'])]
        p['linear_rhs']+=[float(fdot(row,x)*(-1 if s%3==0 else 1)) for s,row in enumerate(c)]
        add('original_linear_'+str(99+3*m),p,kind='linear',known_primal=x,known_dual=dual,exact_rank=198)
    rows=[plane(),plane()];rows[1]['normal']=[-1.,0.,0.]
    rows[0]['new_offset']=rows[1]['new_offset']=.1
    add('infeasible_opposed_normals',problem(3,rows),2,'infeasible')
    p=problem();p['load']=[-1e150,1e150,0.];p['mass_scale']=3.
    add('precision_rejected_before_iteration',p,2,'precision')
    for name,change,message in [
        ('invalid_frame',lambda p:p['contacts'][0].update(tangent=[[1.,0.,0.],[0.,0.,1.]]),'Coulomb contact: invalid row/frame/coefficient'),
        ('invalid_partition',lambda p:p['contacts'][0].update(weights=[.9]),'Coulomb contact: signed partition mismatch'),
        ('invalid_mu',lambda p:p['contacts'][0].update(mu=-.1),'Coulomb contact: invalid row/frame/coefficient'),
        ('invalid_old_geometry',lambda p:p.update(old_position=[-.1,0.,0.]),'Coulomb contact: accepted initial penetration'),
        ('invalid_dt',lambda p:p.update(dt=0.),'Coulomb contact: bounded dimensions/settings'),
        ('invalid_speed',lambda p:p.update(velocity_scale=0.),'Coulomb contact: bounded dimensions/settings')]:
        p=problem();change(p);add(name,p,1,'invalid',expected_message=message)
    for fault,label in [(1,'callback_initial_failure'),(2,'callback_trial_failure'),(3,'callback_nonfinite_rejected')]:
        p=problem();p['fault']=fault;add(label,p,2,'evaluation_failure')
    return cases

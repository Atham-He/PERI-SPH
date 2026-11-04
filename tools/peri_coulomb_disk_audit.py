"""Independent Decimal disk reference and conversion-enclosure assessment."""
from decimal import Decimal as D,localcontext
import math

def dec(x):return D.from_float(x) if isinstance(x,float) else D(str(x))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),D(0))
def norm(x):return dot(x,x).sqrt()

def reference(values):
    l,tx,ty,wx,wy,mu=list(map(dec,values))
    assert mu>=0
    t=[tx,ty];w=[wx,wy];y=[a-b for a,b in zip(t,w)]
    radius=mu*max(l,D(0));size=norm(y);dl=[D(0)]*2
    jt=[[D(0)]*2 for _ in range(2)];jw=[[D(0)]*2 for _ in range(2)]
    if size<radius:
        r=w;branch='stick_interior'
        for i in range(2):jw[i][i]=1
    elif size==0:
        r=t;branch='zero_radius_origin'
        for i in range(2):jt[i][i]=1
    else:
        h=[x/size for x in y];r=[t[i]-radius*h[i] for i in range(2)];branch='slip_exterior'
        for i in range(2):
            dl[i]=-mu*h[i] if l>=0 else D(0)
            for j in range(2):
                jw[i][j]=radius/size*(D(i==j)-h[i]*h[j]);jt[i][j]=D(i==j)-jw[i][j]
    eps=D(2)**-240
    return dict(residual=r,normal_derivative=dl,impulse_derivative=jt,velocity_derivative=jw,branch=branch,
        friction_loss=-dot(t,w),cone_excess=max(D(0),norm(t)-radius),
        maximum_dissipation_error=abs(dot(t,w)+radius*norm(w)),
        scalar_envelope=eps*(1+sum(map(abs,[*t,*w]))+radius),
        derivative_envelope=eps*(2+mu),
        work_envelope=eps*(1+abs(tx*wx)+abs(ty*wy)+radius*norm(w)))

def flatten(x):
    if isinstance(x,(list,tuple)):
        return [v for item in x for v in flatten(item)]
    return [dec(x)]

def inspect(case,records,code):
    with localcontext() as ctx:
        ctx.prec=95
        assert code==case['expected_code'],(case['name'],'return_code',code)
        assert records[0]['type']=='identity'
        if code==1:
            assert len(records)==2 and records[-1]['type']=='error'
            assert records[-1]['message']==case['expected_message']
            assert not records[-1]['accepted_physical_step']
            return dict(name=case['name'],status='declared_invalid_input_rejected')
        assert len(records)==4 and [r['type'] for r in records[1:]]==['disk_input','disk_value','end']
        inp,r,end=records[1:];values=list(map(float,case['tokens']))
        assert [inp['normal_impulse'],*inp['tangent_impulse'],*inp['relative_velocity'],inp['mu']]==values
        assert r['kernel_only'] and not r['accepted_physical_step'] and end['accepted_steps']==end['last_accepted_time']==0
        assert not end['full_S3_complete'] and end['completed']
        assert r['within_budget']==end['arithmetic_budget_satisfied']==(code==0)
        assert all(len(r[field])==2 for field in ['residual','normal_derivative','residual_roundoff_bound'])
        assert all(len(r[field])==2 and all(len(row)==2 for row in r[field]) for field in ['impulse_derivative','velocity_derivative'])
        assert all(math.isfinite(float(x)) for field in ['residual','normal_derivative','impulse_derivative',
            'velocity_derivative','residual_roundoff_bound','maximum_derivative_roundoff_bound','friction_loss',
            'friction_loss_roundoff_bound','cone_excess_upper','maximum_dissipation_error_upper'] for x in flatten(r[field]))
        if code==2:
            assert not r['within_budget']
            assert not r['representable'] or max(r['residual_roundoff_bound']+[r['maximum_derivative_roundoff_bound'],r['friction_loss_roundoff_bound']])>1e-12
            policy=reference(values)
            assert max(policy['scalar_envelope'],policy['work_envelope'])>D('1e-12')
            if not r['representable']:
                assert abs(policy['friction_loss'])>dec(float.fromhex('0x1.fffffffffffffp+1023')) or policy['work_envelope']>dec(float.fromhex('0x1.fffffffffffffp+1023'))
            return dict(name=case['name'],status='declared_precision_range_rejected_no_physical_step',representable=r['representable'])
        assert r['representable'];expected=reference(values);assert r['branch']==expected['branch']
        maximum=D(0)
        for field in ['residual','normal_derivative','impulse_derivative','velocity_derivative','friction_loss']:
            for actual,wanted in zip(flatten(r[field]),flatten(expected[field])):
                error=abs(actual-wanted)/max(D(1),abs(wanted));maximum=max(maximum,error)
                assert error<=D('1e-10'),(case['name'],field,str(error))
        reference_noise=D('1e-85')*max(D(1),max(map(abs,flatten(values))))
        for i in range(2):
            bound=dec(r['residual_roundoff_bound'][i]);assert 0<=bound<=D('1e-12')
            assert abs(dec(r['residual'][i])-expected['residual'][i])<=bound+reference_noise
            nominal_error=abs(dec(r['residual'][i])-expected['residual'][i])+expected['scalar_envelope']
            assert abs(bound-nominal_error)<=reference_noise+4*dec(math.ulp(float(nominal_error)))
        bound=dec(r['maximum_derivative_roundoff_bound']);assert 0<=bound<=D('1e-12')
        derivative_errors=[abs(a-b) for field in ['normal_derivative','impulse_derivative','velocity_derivative']
            for a,b in zip(flatten(r[field]),flatten(expected[field]))]
        assert max(derivative_errors)<=bound+reference_noise
        nominal=max(derivative_errors)+expected['derivative_envelope']
        assert abs(bound-nominal)<=reference_noise+4*dec(math.ulp(float(nominal)))
        bound=dec(r['friction_loss_roundoff_bound']);assert 0<=bound<=D('1e-12')
        assert abs(dec(r['friction_loss'])-expected['friction_loss'])<=bound+reference_noise
        nominal=abs(dec(r['friction_loss'])-expected['friction_loss'])+expected['work_envelope']
        assert abs(bound-nominal)<=reference_noise+4*dec(math.ulp(float(nominal)))
        for emitted,field,envelope in [('cone_excess_upper','cone_excess','scalar_envelope'),
            ('maximum_dissipation_error_upper','maximum_dissipation_error','work_envelope')]:
            actual=dec(r[emitted]);wanted=expected[field]+expected[envelope]
            assert actual+reference_noise>=expected[field]
            assert abs(actual-wanted)<=reference_noise+4*dec(math.ulp(float(wanted)))
        fdmax=D(0)
        if case.get('finite_difference',False):
            epsilon=D('1e-20');v=list(map(dec,values))
            for column in range(5):
                plus=v.copy();minus=v.copy();plus[column]+=epsilon;minus[column]-=epsilon
                a=reference(plus)['residual'];b=reference(minus)['residual']
                for row in range(2):
                    numerical=(a[row]-b[row])/(2*epsilon)
                    derivative=expected['normal_derivative'][row] if column==0 else (
                        expected['impulse_derivative'][row][column-1] if column<3 else expected['velocity_derivative'][row][column-3])
                    fdmax=max(fdmax,abs(numerical-derivative))
            assert fdmax<=D('1e-10'),(case['name'],'independent_derivative_selfcheck',str(fdmax))
        if case.get('disk_root',False):
            assert values[0]>=0 and max(map(abs,expected['residual']))<=D('1e-12')
            assert r['cone_excess_upper']<=1e-12 and r['maximum_dissipation_error_upper']<=1e-12
            assert expected['friction_loss']>=-D('1e-12')
        return dict(name=case['name'],status='disk_map_and_derivative_enclosures_independently_passed',
            maximum_normalized_error=float(maximum),reference_derivative_selfcheck=float(fdmax),
            declared_disk_root=case.get('disk_root',False),physical_global_friction_qualified=False,S3_complete=False)

"""Independent decimal equations and exact positivity/rank certificates.

No native solver is imported or called. These controls are manufactured
strongly monotone operators, not material/energy or sphere qualification.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q

GATE = D('1e-10')

def decimal(value): return D.from_float(value) if isinstance(value,float) else D(str(value))
def dot(a,b): return sum((x*y for x,y in zip(a,b)),D(0))
def magnitude(v): return max(map(abs,v),default=D(0))
def normalized_error(a,b): return abs(a-b)/max(D(1),abs(a),abs(b))

def exact_rank(matrix):
    a=[[Q(str(x)) for x in row] for row in matrix]
    if not a: return 0
    row=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[row],a[pivot]=a[pivot],a[row];p=a[row][col]
        a[row]=[x/p for x in a[row]]
        for i in range(row+1,len(a)):
            factor=a[i][col]
            if factor:a[i]=[x-factor*y for x,y in zip(a[i],a[row])]
        row+=1
        if row==len(a): break
    return row

def strict_monotonicity_certificate(case):
    a=[[Q(str(x)) for x in row] for row in case['A']]
    n=len(a);margins=[]
    for i in range(n):
        margins.append(a[i][i]-sum((abs((a[i][j]+a[j][i])/2) for j in range(n) if j!=i),Q(0)))
    assert min(margins)>0 and all(Q(str(x))>=0 for x in case['beta'])
    # x^T Sym(A)x >= sum_i(diag_i-sum_j|off_ij|)*x_i^2.
    # The nonnegative componentwise cubic adds a monotone operator.
    return str(min(margins))

def equation(case,v):
    n=len(v);a=[[decimal(x) for x in row] for row in case['A']]
    beta=list(map(decimal,case['beta']));load=list(map(decimal,case['load']))
    r=[dot(a[i],v)+beta[i]*v[i]**3-load[i] for i in range(n)]
    j=[[a[i][k]+(3*beta[i]*v[i]**2 if i==k else 0) for k in range(n)] for i in range(n)]
    return r,j

def system(case,v,lam):
    n=len(v);count=len(lam);speed=decimal(case['speed']);impulse=decimal(case['mass_scale'])*speed
    r,j=equation(case,v);c=[[decimal(x) for x in row['C']] for row in case['rows']]
    q=list(map(decimal,case['position']));dt=decimal(case['dt'])
    values=[(r[i]-sum((c[s][i]*lam[s] for s in range(count)),D(0)))/impulse for i in range(n)]
    h=[[D(0)]*(n+count) for _ in range(n+count)]
    for i in range(n):
        for k in range(n):h[i][k]=j[i][k]*speed/impulse
        for s in range(count):h[i][n+s]=-c[s][i]
    gaps=[];kkt=D(0);products=D(0)
    for s,row in enumerate(case['rows']):
        bound=(decimal(row['new'])-dot(c[s],q))/dt
        a=lam[s]/impulse;b=(dot(c[s],v)-bound)/speed;gaps.append(b*speed)
        radius=(a*a+b*b).sqrt();values.append(radius-a-b)
        da=a/radius-1 if radius else D('.5').sqrt()-1
        db=b/radius-1 if radius else D('.5').sqrt()-1
        for i in range(n):h[n+s][i]=db*c[s][i]
        h[n+s][n+s]=da;kkt=max(kkt,-a,-b);products+=abs(a*b)
    return values,h,max(kkt,products),gaps

def inspect(case,raw,code):
    with localcontext() as ctx:
        ctx.prec=95
        assert code==case['expected_code'],('return_code',case['name'],code)
        assert raw['manufactured_operator_only'] is True
        maximum=D(0)
        def close(a,b,label):
            nonlocal maximum
            error=normalized_error(decimal(a),decimal(b));maximum=max(maximum,error)
            assert error<=GATE,(case['name'],label,str(error),str(a),str(b))
        if case['mode']==1:
            a=case['A'];rhs=case['rhs'];rank=exact_rank(a)
            if code:
                assert exact_rank([row+[b] for row,b in zip(a,rhs)])>rank
                assert 'inconsistent' in raw['error']
                return dict(name=case['name'],status='independently_certified_inconsistent_linear_rejection',rank=rank)
            x=list(map(decimal,raw['value']));assert len(x)==len(a)
            assert raw['linear_rank']==rank
            for row,b in zip(a,rhs):close(dot(list(map(decimal,row)),x),b,'original_linear_equation')
            for actual,wanted in zip(x,case['minimum_norm_reference']):close(actual,wanted,'minimum_norm_solution')
            assert raw['linear_residual']<=1e-12
            return dict(name=case['name'],status='independent_original_redundant_linear_system_passed',maximum_normalized_error=float(maximum),rank=rank)
        if code==1:
            label=case['invalid']
            if label=='zero_row':assert any(all(decimal(x)==0 for x in r['C']) for r in case['rows'])
            elif label=='initial_penetration':
                q=list(map(decimal,case['position']));assert any(dot(list(map(decimal,r['C'])),q)<decimal(r['old']) for r in case['rows'])
            elif label=='negative_cubic':assert any(decimal(x)<0 for x in case['beta'])
            elif label=='zero_mass_scale':assert decimal(case['mass_scale'])==0
            else:raise AssertionError(('unproved invalidity',label))
            assert raw['error']
            return dict(name=case['name'],status='independently_certified_invalid_input_rejection',reason=label)
        margin=strict_monotonicity_certificate(case)
        n=len(case['position']);count=len(case['rows']);trace=raw['trace'];assert trace
        speed=decimal(case['speed']);impulse=decimal(case['mass_scale'])*speed
        for index,t in enumerate(trace):
            assert t['iteration']==index
            v=list(map(decimal,t['velocity']));lam=list(map(decimal,t['normal_impulse']))
            assert len(v)==n and len(lam)==count
            f,h,kkt,_=system(case,v,lam)
            close(t['merit'],magnitude(f),'trace_residual');close(t['squared_norm'],dot(f,f),'trace_merit');close(t['KKT'],kkt,'trace_KKT')
            native_f=list(map(decimal,t['normalized_residual']))
            native_h=[list(map(decimal,row)) for row in t['generalized_jacobian']]
            native_gaps=list(map(decimal,t['normalized_gaps']))
            assert len(native_f)==len(native_h)==n+count and len(native_gaps)==count
            for actual,expected in zip(native_f,f):close(actual,expected,'recorded_original_residual')
            # Generalized derivatives at a nearly biactive corner are sensitive
            # to roundoff in the signed gap. Verify that measured gap against
            # exact geometry, then verify the selected scalar derivative using
            # that recorded floating-point argument. Never infer it from a
            # decimal-rendered coordinate that has already lost that cursor.
            _,_,_,physical_gaps=system(case,v,lam)
            for s,b in enumerate(native_gaps):
                close(b,physical_gaps[s]/speed,'recorded_signed_gap')
                a=lam[s]/impulse;radius=(a*a+b*b).sqrt()
                da=a/radius-1 if radius else D('.5').sqrt()-1
                db=b/radius-1 if radius else D('.5').sqrt()-1
                for i,c in enumerate(case['rows'][s]['C']):h[n+s][i]=db*decimal(c)
                h[n+s][n+s]=da
            for actual,expected in zip(native_h,h):
                assert len(actual)==n+count
                for x,y in zip(actual,expected):close(x,y,'full_generalized_Jacobian')
            direction=list(map(decimal,t['direction']))
            if direction:
                assert len(direction)==n+count and 0<t['linear_rank']<=n+count
                norm_h=max((sum(map(abs,row),D(0)) for row in native_h),default=D(0))
                residual=magnitude([dot(row,direction)+b for row,b in zip(native_h,native_f)])
                denominator=norm_h*magnitude(direction)+magnitude(native_f)
                independent_linear=residual/denominator if denominator else residual
                maximum=max(maximum,independent_linear);assert independent_linear<=D('1e-12')
                assert t['linear_residual']<=1e-12
            if index+1<len(trace):
                assert direction and 0<decimal(t['alpha'])<=1
                alpha=decimal(t['alpha']);new=trace[index+1]
                for i in range(n):close(new['velocity'][i],v[i]+alpha*speed*direction[i],'Newton_velocity_update')
                for s in range(count):close(new['normal_impulse'][s],lam[s]+alpha*impulse*direction[n+s],'Newton_impulse_update')
                next_f,_,_,_=system(case,list(map(decimal,new['velocity'])),list(map(decimal,new['normal_impulse'])))
                assert dot(next_f,next_f)<=(1-D('1e-4')*alpha)*dot(f,f)+D('1e-25')
        if code==2:
            assert not raw['converged'] and raw['velocity']==[] and raw['normal_impulse']==[]
            last=trace[-1];label=case['rejection']
            if label=='Newton_limit':
                assert raw['status']=='Newton_iteration_limit' and last['iteration']==case['max_newton']
                assert max(decimal(last['merit']),decimal(last['KKT']))>decimal(case['tolerance'])
            elif label=='domain_blocked':
                assert raw['status']=='line_search_limit' and last['domain_rejections']==case['max_line_search']
                for k in range(case['max_line_search']):
                    v=[decimal(x)+(D(2)**(-k))*speed*decimal(d) for x,d in zip(last['velocity'],last['direction'][:n])]
                    assert magnitude(v)>decimal(case['domain_radius'])
            elif label=='infeasible':
                a,b=case['rows'];assert all(decimal(x)==-decimal(y) for x,y in zip(a['C'],b['C']))
                assert decimal(a['new'])+decimal(b['new'])>0 and all(decimal(x)==0 for x in case['position'])
                assert raw['status'] in ['line_search_limit','Newton_iteration_limit','original_linear_system_rejected']
            else:raise AssertionError(('unproved rejection',label))
            return dict(name=case['name'],status='independent_terminal_rejection_no_accepted_state',reason=label,maximum_normalized_error=float(maximum))
        assert raw['converged'] and raw['status']=='all_rows_KKT_root_accepted'
        v=list(map(decimal,raw['velocity']));lam=list(map(decimal,raw['normal_impulse']))
        target=list(map(decimal,case['root_reference']));expected_lam=list(map(decimal,case['dual_witness']))
        assert len(v)==n and len(lam)==count and len(target)==n and len(expected_lam)==count
        f,_,kkt,gaps=system(case,v,lam);ref_f,_,ref_kkt,_=system(case,target,expected_lam)
        assert magnitude(ref_f)<=D('1e-50') and ref_kkt<=D('1e-50'),('invalid known root',case['name'])
        assert magnitude(f)<=GATE and kkt<=GATE and raw['normalized_KKT']<=1e-12
        for actual,wanted in zip(v,target):close(actual,wanted,'unique_velocity_root')
        r,j=equation(case,v)
        for i in range(n):
            close(raw['momentum_equation'][i],r[i],'full_operator')
            for k in range(n):
                # Independent central difference of the cubic equation, not
                # the native analytic tangent or its symmetrization.
                delta=D('1e-20');left=v.copy();right=v.copy();left[k]-=delta;right[k]+=delta
                if i==k:
                    numerical=(equation(case,right)[0][i]-equation(case,left)[0][i])/(2*delta)
                else:numerical=decimal(case['A'][i][k])
                close(raw['jacobian'][i][k],numerical,'global_nonsymmetric_Jacobian')
        c=[[decimal(x) for x in row['C']] for row in case['rows']]
        scatter=[sum((c[s][i]*lam[s] for s in range(count)),D(0)) for i in range(n)]
        expected=[sum((c[s][i]*expected_lam[s] for s in range(count)),D(0)) for i in range(n)]
        for actual,wanted in zip(scatter,expected):close(actual,wanted,'unique_transpose_scatter')
        active=[case['rows'][s]['C'] for s in range(count) if decimal(case['gap_reference'][s])==0]
        return dict(name=case['name'],status='independent_many_row_nonsymmetric_root_and_trace_passed',
            dofs=n,contact_rows=count,maximum_normalized_error=float(maximum),exact_strong_monotonicity_margin=margin,
            active_row_rank=exact_rank(active),active_rows=len(active),dual_uniqueness_claimed=False,
            minimum_gap=float(min(gaps,default=D(0))),minimum_impulse=float(min(lam,default=D(0))),
            iterations=trace[-1]['iteration'],scatter=list(map(float,scatter)),velocity=list(map(float,v)))

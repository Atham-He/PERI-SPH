"""Original exact-Fraction/95-digit simultaneous normal/tangent reference."""
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import math
from peri_pd_reference_audit import dec,dot,close_tree

def frac(x):
    if isinstance(x,F):return x
    if isinstance(x,float):return F.from_float(x)
    return F(x)
def decimal(x):
    if isinstance(x,F):return D(x.numerator)/D(x.denominator)
    return dec(x)
def fdot(a,b):return sum((frac(x)*frac(y) for x,y in zip(a,b)),F(0))
def absolute(x):return -x if x<0 else x
def fraction_disk(l,t,w,mu):
    """Branch sign is exact before converting exact differences to Decimal."""
    l,mu=frac(l),frac(mu);t=list(map(frac,t));w=list(map(frac,w))
    radius=mu*max(l,F(0));y=[a-b for a,b in zip(t,w)];squared=fdot(y,y)
    td,wd,yd=list(map(decimal,t)),list(map(decimal,w)),list(map(decimal,y));rd=decimal(radius)
    jt=[[D(0)]*2 for _ in range(2)];jw=[[D(0)]*2 for _ in range(2)];dl=[D(0)]*2
    if squared<radius*radius:
        value=wd;branch='stick_interior'
        for i in range(2):jw[i][i]=1
    elif squared==0:
        value=td;branch='zero_radius_origin'
        for i in range(2):jt[i][i]=1
    else:
        size=decimal(squared).sqrt();h=[a/size for a in yd];value=[td[i]-rd*h[i] for i in range(2)]
        branch='slip_exterior'
        for i in range(2):
            dl[i]=-decimal(mu)*h[i] if l>=0 else D(0)
            for j in range(2):
                jw[i][j]=rd/size*(D(i==j)-h[i]*h[j]);jt[i][j]=D(i==j)-jw[i][j]
    tn=decimal(fdot(t,t)).sqrt();wn=decimal(fdot(w,w)).sqrt();work=decimal(fdot(t,w))
    return dict(residual=value,normal_derivative=dl,impulse_derivative=jt,velocity_derivative=jw,branch=branch,
        friction_loss=-work,cone_excess=max(D(0),tn-rd),maximum_dissipation_error=abs(work+rd*wn))

def coefficients(contacts,n):
    rows=[]
    for c in contacts:
        assert len(c['weights'])*3==n and len(c['tangent'])==2
        assert abs(sum(map(frac,c['weights']),F(0))-1)<=frac(1e-12)
        basis=[c['normal'],*c['tangent']]
        for i in range(3):
            assert len(basis[i])==3
            for j in range(3):assert abs(fdot(basis[i],basis[j])-int(i==j))<=frac(1e-12)
        assert c['mu']>=0 and len(c['wall_tangent_velocity'])==2
        # The original discrete coefficient is the native binary64 product.
        for vector in basis:rows.append([float(w)*float(a) for w in c['weights'] for a in vector])
    return rows

def operator_reference(header,v):
    a=header['operator_matrix'];load=header['operator_load'];beta=header['operator_cubic_beta'];n=len(v)
    qv=list(map(frac,v));r=[];j=[]
    for i in range(n):
        r.append(decimal(fdot(a[i],qv)-frac(load[i])+frac(beta[i])*qv[i]**3))
        row=list(map(dec,a[i]));row[i]+=decimal(3*frac(beta[i])*qv[i]**2);j.append(row)
    return r,j

def original_reference(old,contacts,dt,I,speed,v,impulses,material_R,material_J,gather=None):
    n=len(v);m=len(contacts);size=n+3*m
    c=coefficients(contacts,n) if gather is None else gather
    qI,qV,qdt=frac(I),frac(speed),frac(dt);p=list(map(frac,impulses));fv=list(map(frac,v))
    residual=[D(0)]*size;jac=[[D(0)]*size for _ in range(size)];gaps=[];relative=[];disks=[]
    for i in range(n):
        residual[i]=(dec(material_R[i])-decimal(sum((frac(c[s][i])*p[s] for s in range(3*m)),F(0))))/decimal(qI)
        for k in range(n):jac[i][k]=dec(material_J[i][k])*decimal(qV/qI)
        for s in range(3*m):jac[i][n+s]=-dec(c[s][i])
    products=F(0);KKT=D(0);work_error=D(0)
    for s,contact in enumerate(contacts):
        row=n+3*s;l=p[3*s]/qI
        gap=(fdot(c[3*s],fv)+(fdot(c[3*s],old)-frac(contact['new_offset']))/qdt)/qV
        gaps.append(gap);products+=abs(l*gap);KKT=max(KKT,decimal(max(-l,-gap,F(0))))
        a,b=decimal(l),decimal(gap);radius=(a*a+b*b).sqrt()
        phi=-2*a*b/(radius+a+b) if a+b>0 else radius-a-b
        da=a/radius-1 if radius else D('.5').sqrt()-1
        db=b/radius-1 if radius else da
        residual[row]=phi;jac[row][row]=da
        for i in range(n):jac[row][i]=db*dec(c[3*s][i])
        t=[p[3*s+1+k]/qI for k in range(2)]
        w=[(fdot(c[3*s+1+k],fv)-frac(contact['wall_tangent_velocity'][k]))/qV for k in range(2)]
        relative.append(w);disk=fraction_disk(l,t,w,frac(contact['mu']));disks.append(disk)
        KKT=max(KKT,disk['cone_excess']);work_error+=disk['maximum_dissipation_error']
        for a in range(2):
            residual[row+1+a]=disk['residual'][a];jac[row+1+a][row]=disk['normal_derivative'][a]
            for k in range(2):jac[row+1+a][row+1+k]=disk['impulse_derivative'][a][k]
            for i in range(n):
                jac[row+1+a][i]=sum((disk['velocity_derivative'][a][k]*dec(c[3*s+1+k][i]) for k in range(2)),D(0))
    KKT=max(KKT,decimal(products),work_error)
    return dict(residual=residual,jacobian=jac,gaps=gaps,relative_velocity=relative,disks=disks,KKT=KKT)

def check_original(old,contacts,dt,I,speed,v,p,native,material_R,material_J,gather=None):
    """Check physical R/J and exact supplied-R accumulation separately."""
    n=len(v);m=len(contacts);size=n+3*m;c=coefficients(contacts,n) if gather is None else gather
    assert native['precision_budget_ok']
    assert len(native['residual'])==len(native['residual_roundoff_bound'])==len(native['jacobian'])==size
    assert all(len(row)==size for row in native['jacobian'])
    assert len(native['normalized_gaps'])==len(native['gap_roundoff_bound'])==len(native['disk'])==m
    assert len(native['normalized_relative_velocity'])==len(native['relative_velocity_roundoff_bound'])==m
    assert len(native['supplied_momentum'])==len(native['supplied_jacobian'])==n
    assert all(len(row)==n for row in native['supplied_jacobian'])
    maximum=close_tree(native['supplied_momentum'],material_R,dec(I))
    maximum=max(maximum,close_tree(native['supplied_jacobian'],material_J,dec(I)/dec(speed)))
    reference=original_reference(old,contacts,dt,I,speed,v,p,material_R,material_J,c)
    maximum=max(maximum,close_tree(native['residual'],reference['residual']),close_tree(native['jacobian'],reference['jacobian']))
    assert all(0<=b<=1e-12 for b in native['residual_roundoff_bound'])
    assert 0<=native['maximum_jacobian_roundoff_bound']<=1e-12
    supplied=original_reference(old,contacts,dt,I,speed,v,p,native['supplied_momentum'],native['supplied_jacobian'],c)
    for i in range(n):
        exact=(frac(native['supplied_momentum'][i])-sum((frac(c[s][i])*frac(p[s]) for s in range(3*m)),F(0)))/frac(I)
        assert native['residual'][i]==float(exact)
        assert abs(frac(native['residual'][i])-exact)<=frac(native['residual_roundoff_bound'][i])
    for row in range(n,size):
        noise=D('1e-80')*max(D(1),abs(supplied['residual'][row]))
        assert abs(dec(native['residual'][row])-supplied['residual'][row])<=dec(native['residual_roundoff_bound'][row])+noise
    for actual,wanted in zip(native['jacobian'],supplied['jacobian']):
        for a,b in zip(actual,wanted):
            noise=D('1e-80')*max(D(1),abs(b))
            assert abs(dec(a)-b)<=dec(native['maximum_jacobian_roundoff_bound'])+noise
    for s in range(m):
        assert 0<=native['gap_roundoff_bound'][s]<=1e-12
        assert abs(frac(native['normalized_gaps'][s])-reference['gaps'][s])<=frac(native['gap_roundoff_bound'][s])
        for k in range(2):
            bound=native['relative_velocity_roundoff_bound'][s][k];assert 0<=bound<=1e-12
            assert abs(frac(native['normalized_relative_velocity'][s][k])-reference['relative_velocity'][s][k])<=frac(bound)
        disk=native['disk'][s];expected=reference['disks'][s]
        assert disk['within_budget'] and disk['representable'] and disk['branch']==expected['branch']
        assert native['residual'][n+3*s+1:n+3*s+3]==disk['residual']
        assert native['residual_roundoff_bound'][n+3*s+1:n+3*s+3]==disk['residual_roundoff_bound']
        for field in ['residual','normal_derivative','impulse_derivative','velocity_derivative','friction_loss']:
            maximum=max(maximum,close_tree(disk[field],expected[field]))
        assert 0<=disk['maximum_derivative_roundoff_bound']<=1e-12
        for field in ['normal_derivative','impulse_derivative','velocity_derivative']:
            values=[disk[field]] if field=='normal_derivative' else disk[field]
            wants=[expected[field]] if field=='normal_derivative' else expected[field]
            for row,wanted in zip(values,wants):
                for a,b in zip(row,wanted):
                    assert abs(dec(a)-b)<=dec(disk['maximum_derivative_roundoff_bound'])+D('1e-80')*max(D(1),abs(b))
        assert 0<=disk['friction_loss_roundoff_bound']<=1e-12
        assert abs(dec(disk['friction_loss'])-expected['friction_loss'])<=dec(disk['friction_loss_roundoff_bound'])+D('1e-80')*max(D(1),abs(expected['friction_loss']))
        noise=D('1e-80')*max(D(1),abs(expected['friction_loss']))
        for emitted,field in [('cone_excess_upper','cone_excess'),('maximum_dissipation_error_upper','maximum_dissipation_error')]:
            assert dec(disk[emitted])+noise>=expected[field]
            maximum=max(maximum,close_tree(disk[emitted],expected[field],max(D(1),abs(expected[field]))))
    assert dec(native['KKT'])+D('1e-80')*max(D(1),abs(reference['KKT']))>=reference['KKT']
    maximum=max(maximum,close_tree(native['KKT'],reference['KKT'],max(D(1),abs(reference['KKT']))))
    expected_merit=max((abs(r)+b for r,b in zip(native['residual'],native['residual_roundoff_bound'])),default=0)
    assert native['merit']==expected_merit
    squared=sum((dec(r)**2 for r in native['residual']),D(0))
    maximum=max(maximum,close_tree(native['squared_norm'],squared,max(D(1),squared)))
    return reference,maximum

def check_trace(old,contacts,dt,I,speed,solve,material,gather=None,require_root=True):
    n=len(old);size=n+3*len(contacts);trace=solve['trace'];maximum=0
    if require_root:assert solve['converged'] and trace
    for index,t in enumerate(trace):
        assert t['iteration']==index
        R,J=material(t['velocity'],t['original'])
        reference,error=check_original(old,contacts,dt,I,speed,t['velocity'],t['impulses'],t['original'],R,J,gather)
        maximum=max(maximum,error)
        direction=t['normalized_direction'];native=t['original']
        if direction:
            assert len(direction)==size and 0<t['linear_rank']<=size and t['linear_residual']<=1e-12
            residual=max(abs(fdot(row,direction)+frac(rhs)) for row,rhs in zip(native['jacobian'],native['residual']))
            scale=max(sum(map(absolute,map(frac,row)),F(0)) for row in native['jacobian'])*max(map(absolute,map(frac,direction)))+max(map(absolute,map(frac,native['residual'])))
            assert residual/scale<=frac(1e-12) if scale else residual==0
        if index+1<len(trace):
            nxt=trace[index+1];alpha=t['alpha'];assert direction and 0<alpha<=1
            maximum=max(maximum,close_tree(nxt['velocity'],[dec(x)+dec(alpha)*dec(speed)*dec(d) for x,d in zip(t['velocity'],direction[:n])]))
            maximum=max(maximum,close_tree(nxt['impulses'],[dec(x)+dec(alpha)*dec(I)*dec(d) for x,d in zip(t['impulses'],direction[n:])],dec(I)))
            native_armijo=(1-1e-4*alpha)*native['squared_norm']
            assert nxt['original']['squared_norm']<=math.nextafter(native_armijo,math.inf)
    if require_root:
        assert not solve['precision_rejection']['present'] and not solve['evaluation_failure']['present']
        last=trace[-1]['original'];assert max(last['merit'],last['KKT'],solve['normalized_KKT'])<=1e-12
        maximum=max(maximum,close_tree(solve['velocity'],trace[-1]['velocity']),close_tree(solve['impulses'],trace[-1]['impulses'],dec(I)))
        assert max(map(abs,reference['residual']))<=D('1e-10') and reference['KKT']<=D('1e-10')
    return maximum

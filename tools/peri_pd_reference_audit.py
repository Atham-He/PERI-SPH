"""Independent exact/95-digit actual PD reference and signed MLS geometry."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from functools import lru_cache
import math
from peri_normal_projection_audit import solve as rational_solve
from peri_dense_contact_audit import determinant as rational_det
from peri_many_contact_audit import exact_rank

def dec(x):return D.from_float(x) if isinstance(x,float) else D(str(x))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),D(0))
def flatten(values):return [x for row in values for x in row]
def norm(values):return dot(values,values).sqrt()
def close_tree(actual,wanted,scale=D(1),gate=D('1e-10')):
    if isinstance(wanted,(list,tuple)):
        assert len(actual)==len(wanted),('shape',len(actual),len(wanted))
        return max((close_tree(x,y,scale,gate) for x,y in zip(actual,wanted)),default=0.)
    error=abs(dec(actual)-dec(wanted))/max(D(1),dec(scale))
    assert error<=gate,('reference_cache',str(error),str(actual),str(wanted))
    return float(error)

def decimal_inverse(matrix):
    n=len(matrix);a=[list(row)+[D(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for col in range(n):
        pivot=max(range(col,n),key=lambda i:abs(a[i][col]))
        assert a[pivot][col]!=0,'reference shape singular'
        a[col],a[pivot]=a[pivot],a[col];p=a[col][col];a[col]=[x/p for x in a[col]]
        for i in range(n):
            if i==col:continue
            m=a[i][col];a[i]=[x-m*y for x,y in zip(a[i],a[col])]
    return [row[n:] for row in a]

def jacobi_eigenvalues(matrix):
    """Independent cyclic Jacobi rotations, not native Eigen factorization."""
    a=[list(map(float,row)) for row in matrix];n=len(a);scale=max(abs(x) for x in flatten(a))
    if scale==0:return [0.]*n
    for sweep in range(40):
        maximum=max((abs(a[p][q]) for p in range(n) for q in range(p+1,n)),default=0.)
        if maximum<=1e-14*scale:return sorted(a[i][i] for i in range(n))
        for p in range(n):
            for q in range(p+1,n):
                off=a[p][q]
                if abs(off)<=1e-16*scale:continue
                tau=(a[q][q]-a[p][p])/(2*off)
                t=math.copysign(1.,tau)/(abs(tau)+math.hypot(1.,tau))
                c=1/math.hypot(1.,t);s=t*c
                a[p][p]-=t*off;a[q][q]+=t*off;a[p][q]=a[q][p]=0.
                for k in range(n):
                    if k==p or k==q:continue
                    x,y=a[k][p],a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y
    raise AssertionError('independent Jacobi spectrum did not converge')

def base_sphere(rotation):
    points=[]
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-2,3):
                if i*i+j*j+k*k<=4:
                    x,y,z=Q(i,2),Q(j,2),Q(k,2)
                    points.append([-y,x,z] if rotation else [x,y,z])
    return points

def qdecimal(x):return D(x.numerator)/D(x.denominator)

@lru_cache(None)
def expected_reference(rotation,volumes,horizon):
    with localcontext() as ctx:
        ctx.prec=95
        x=base_sphere(rotation);n=len(x);v=[Q.from_float(z) for z in volumes];delta=Q.from_float(horizon)
        offsets=[0];targets=[];weights=[];inverse=[];gradients=[];measure=[];ratios=[];local=[];rigidity=[]
        for i in range(n):
            shape=[[Q(0)]*3 for _ in range(3)];neighbors=[];local_weights=[]
            for j in range(n):
                if i==j:continue
                xi=[x[j][a]-x[i][a] for a in range(3)];radius2=sum(t*t for t in xi)
                if radius2>=delta*delta:continue
                w=(1-radius2/(delta*delta))**3;neighbors.append(j);local_weights.append(w)
                for a in range(3):
                    for b in range(3):shape[a][b]+=w*v[j]*xi[a]*xi[b]
                if j>i:
                    row=[Q(0)]*(3*n)
                    for a in range(3):row[3*j+a]=xi[a];row[3*i+a]=-xi[a]
                    rigidity.append(row)
            assert all(rational_det([row[:k] for row in shape[:k]])>0 for k in range(1,4))
            inv_columns=[rational_solve(shape,[Q(int(k==c)) for k in range(3)]) for c in range(3)]
            inv=[[inv_columns[b][a] for b in range(3)] for a in range(3)]
            gi=[[Q(0)]*3 for _ in range(n)]
            for j,w in zip(neighbors,local_weights):
                xi=[x[j][a]-x[i][a] for a in range(3)]
                for a in range(3):
                    gi[j][a]=w*v[j]*sum(inv[a][b]*xi[b] for b in range(3));gi[i][a]-=gi[j][a]
            assert all(sum(gi[j][a] for j in range(n))==0 for a in range(3))
            for a in range(3):
                for b in range(3):assert sum(x[j][a]*gi[j][b] for j in range(n))==int(a==b)
            eig=jacobi_eigenvalues([[float(z) for z in row] for row in shape]);assert eig[0]/eig[-1]>1e-12
            ratios.append(D(str(eig[0]/eig[-1])))
            targets.extend(neighbors);weights.extend(local_weights);offsets.append(len(targets));inverse.append(flatten(inv));gradients.extend(gi)
            measure.append(sum(w*v[j] for j,w in zip(neighbors,local_weights)));local.append([i]+neighbors)
        assert exact_rank(rigidity)==3*n-6,'intact bond rigidity requires exactly six rigid modes'
        g=[[qdecimal(t) for t in row] for row in gradients];l=[[D(0)]*n for _ in range(n)]
        for i in range(n):
            for p in range(offsets[i],offsets[i+1]):
                j=targets[p];xi=[qdecimal(x[j][a]-x[i][a]) for a in range(3)]
                h=[D(int(k==j)-int(k==i))-dot(g[i*n+k],xi) for k in range(n)]
                coefficient=qdecimal(v[i]*weights[p]*v[j]/(delta*delta*measure[i]))
                nonzero=[k for k in range(n) if h[k]]
                for a in nonzero:
                    for b in nonzero:l[a][b]+=coefficient*h[a]*h[b]
        return dict(position=[[qdecimal(t) for t in row] for row in x],volume=list(map(qdecimal,v)),
            offsets=offsets,target=targets,weight=list(map(qdecimal,weights)),active=[1]*len(targets),
            inverse_shape=[[qdecimal(t) for t in row] for row in inverse],shape_gradient=g,neighbor_measure=list(map(qdecimal,measure)),
            shape_eigenvalue_ratio=ratios,nonaffine_unit=flatten(l),L=l,local=local,rigidity_rank=3*n-6)

def expected_surface_points(kind,rotation,azimuths,cap):
    points=[];areas=[];pi=math.pi
    def rotate(p):return [-p[1],p[0],p[2]] if rotation else p
    if kind=='surface':
        mu=[0.,(1-cap)/2,1-cap,1-cap/2,1.]
        weights=[(1-cap)/6,2*(1-cap)/3,1/6,2*cap/3,cap/6]
        for k in range(5):
            for side in [-1,1]:
                if k==0 and side==1:continue
                w=2*weights[k] if k==0 else weights[k]
                if k==4:points.append(rotate([0.,side*1.,0.]));areas.append(2*pi*w);continue
                radial=math.sqrt((1-mu[k])*(1+mu[k]))
                for a in range(azimuths):
                    angle=2*pi*(a+.5)/azimuths;points.append(rotate([radial*math.cos(angle),side*mu[k],radial*math.sin(angle)]));areas.append(2*pi*w/azimuths)
    else:
        big=math.sqrt((3+2*math.sqrt(6/5))/7);small=math.sqrt((3-2*math.sqrt(6/5))/7)
        roots=[big,small,-small,-big];wg=[(18-math.sqrt(30))/72,(18+math.sqrt(30))/72,(18+math.sqrt(30))/72,(18-math.sqrt(30))/72]
        for panel in range(5):
            lower=0. if panel==0 else 2.**(panel-5);upper=2.**(panel-4)
            for root,w in zip(roots,wg):
                z=lower+(upper-lower)*(1+root)/2;mu=1-z*z;weight=(upper-lower)*w;radial=math.sqrt((1-mu)*(1+mu))
                for side in [-1,1]:
                    for a in range(16):
                        angle=2*pi*(a+.5)/16;points.append(rotate([radial*math.cos(angle),side*mu,radial*math.sin(angle)]));areas.append(2*z*weight*2*pi/16)
    return points,areas

@lru_cache(None)
def expected_mls(particles,samples,support):
    with localcontext() as ctx:
        ctx.prec=95
        x=[list(map(dec,p)) for p in particles];z=[list(map(dec,p)) for p in samples];delta=dec(support)
        offsets=[0];targets=[];weights=[]
        for point in z:
            moment=[[D(0)]*4 for _ in range(4)];basis=[];kernel=[];local=[]
            for i,particle in enumerate(x):
                relative=[(a-b)/delta for a,b in zip(particle,point)];distance=dot(relative,relative)
                if distance>=1:continue
                row=[D(1),*relative];w=(1-distance)**3;basis.append(row);kernel.append(w);local.append(i)
                for a in range(4):
                    for b in range(4):moment[a][b]+=w*row[a]*row[b]
            inv=decimal_inverse(moment);coeff=[row[0] for row in inv]
            phi=[w*dot(row,coeff) for w,row in zip(kernel,basis)]
            assert abs(sum(phi,D(0))-1)<D('1e-80')
            for a in range(3):assert abs(sum((v*(x[j][a]-point[a]) for j,v in zip(local,phi)),D(0)))<D('1e-80')
            targets.extend(local);weights.extend(phi);offsets.append(len(targets))
        return dict(offsets=offsets,targets=targets,weights=weights)

def check_geometry(case,header):
    with localcontext() as ctx:
        ctx.prec=95
        rotation=case['rotation'];r=header['reference'];n=33
        assert header['actual_volume_nodes']==n and r['generation']==1 and r['horizon']==1.25
        assert header['normal_only'] and not header['full_S3_complete'] and header['normalized_benchmark']
        expected=expected_reference(rotation,tuple(r['volume']),r['horizon']);maximum=0.
        def compare(actual,wanted):
            nonlocal maximum
            maximum=max(maximum,close_tree(actual,wanted))
        compare(r['position'],expected['position']);compare(r['volume'],[4*math.pi/(3*n)]*n)
        for field in ['offsets','target','active']:assert r[field]==expected[field]
        for field in ['weight','inverse_shape','shape_gradient','neighbor_measure','shape_eigenvalue_ratio','nonaffine_unit']:compare(r[field],expected[field])
        compare(header['physical_mass'],[dec(case['rho0'])*v for v in expected['volume']]);compare(header['rho0'],case['rho0']);compare(header['eta'],case['eta'])
        for kind in ['surface','dense_surface']:
            actual=header[kind];points,areas=expected_surface_points(kind,rotation,case['azimuths'],case['cap_width'])
            assert actual['generation']==1 and actual['particle_count']==n
            compare(actual['support'],case['support']);compare(actual['reference'],points);compare(actual['areas'],areas)
            assert all(x>0 for x in actual['areas'])
            key_x=tuple(tuple(float(x) for x in row) for row in r['position']);key_s=tuple(tuple(float(x) for x in row) for row in actual['reference'])
            mls=expected_mls(key_x,key_s,actual['support'])
            assert actual['offsets']==mls['offsets'] and actual['targets']==mls['targets'];compare(actual['weights'],mls['weights'])
        normals=[[1,0,0],[-1,0,0]] if rotation else [[0,1,0],[0,-1,0]]
        if rotation:normals=[[-1,0,0],[1,0,0]]
        compare(header['plate_normal'],normals)
        coefficients=[];ids=[];surface=header['surface']
        if case['mode']!=2:
            for s in range(len(surface['reference'])):
                for plate in range(2):
                    row=[D(0)]*(3*n)
                    for p in range(surface['offsets'][s],surface['offsets'][s+1]):
                        j=surface['targets'][p]
                        for a in range(3):row[3*j+a]=dec(surface['weights'][p])*normals[plate][a]
                    coefficients.append(row);ids.append(plate)
        compare(header['C'],coefficients);assert header['plate_ids']==ids
        assert header['dense_penetration_stroke_gate']==.01 and not header['area_union_is_integral']
        return expected,dict(maximum_normalized_cache_error=maximum,exact_rigidity_rank=expected['rigidity_rank'],
            nodes=33,physical_DOF=99,normal_rows=len(coefficients),dense_observation_points=len(header['surface']['reference'])+len(header['dense_surface']['reference']))

@lru_cache(None)
def stability_spectrum(rotation,volumes,horizon,bulk,mu,eta):
    assert bulk>0 and mu>0 and eta>0
    ref=expected_reference(rotation,volumes,horizon);n=33
    g=[list(map(float,row)) for row in ref['shape_gradient']]
    k=[[0.]*(3*n) for _ in range(3*n)]
    for i in range(n):
        for j in ref['local'][i]:
            gj=g[i*n+j]
            for p in ref['local'][i]:
                gp=g[i*n+p];inner=sum(x*y for x,y in zip(gj,gp))
                for a in range(3):
                    for c in range(3):
                        k[3*j+a][3*p+c]+=volumes[i]*(mu*(inner if a==c else 0)+mu*gj[c]*gp[a]+(bulk-2*mu/3)*gj[a]*gp[c])
    for i in range(n):
        for j in range(n):
            for a in range(3):k[3*i+a][3*j+a]+=eta*mu*float(ref['L'][i][j])
    eig=jacobi_eigenvalues(k);scale=eig[-1];assert scale>0
    assert max(map(abs,eig[:6]))<=1e-10*scale and eig[6]>1e-9*scale
    # Exact rigidity rank plus positive material/H^TH factors proves exactly
    # six ideal rigid modes. Independent Jacobi checks the numerical gap.
    return dict(exact_rigidity_rank=ref['rigidity_rank'],rigid_modes=6,
        seventh_eigenvalue=eig[6],largest_eigenvalue=scale,relative_positive_gap=eig[6]/scale)

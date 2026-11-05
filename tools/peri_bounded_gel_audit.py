"""Independent actual-PD normal sphere audit; no native solver is called."""
from decimal import Decimal as D, localcontext
import json
from peri_pd_reference_audit import dec,dot,norm,flatten,close_tree,check_geometry,stability_spectrum,expected_mls
from peri_pd_material_reference import material_step
from peri_gel_affine_audit import green,sub,pk1,multiply,transpose,determinant

def vectors(values):return [list(map(dec,row)) for row in values]
def flat_vectors(values):return [list(values[3*i:3*i+3]) for i in range(len(values)//3)]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def add3(a,b):return [x+y for x,y in zip(a,b)]

def invalid_message(flag):
    # The injected non-deviatoric Green history is rejected by the material
    # history guard before the later accepted-F/history cursor comparison.
    return {1:'bounded gel model: physical mass mismatch',2:'invalid gel material history',
        3:'fixed PD reference: invalid directed packet',4:'stored signed surface weight mismatch'}[flag]
def matrix_vector(m,v):return [dot(row,v) for row in m]

def cubic_minimum(coefficients):
    c0,c1,c2,c3=coefficients;candidates=[D(0),D(1)]
    a,b,c=3*c3,2*c2,c1
    if a==0:
        if b!=0:candidates.append(-c/b)
    else:
        discriminant=b*b-4*a*c
        if discriminant>=0:candidates.extend([(-b-discriminant.sqrt())/(2*a),(-b+discriminant.sqrt())/(2*a)])
    return min(c0+s*(c1+s*(c2+s*c3)) for s in candidates if 0<=s<=1)

def deform(ref,i,q):
    n=len(q);g=ref['shape_gradient']
    return [sum((q[k][a]*g[i*n+k][b] for k in range(n)),D(0)) for a in range(3) for b in range(3)]

def nonaffine_energy(ref,q,eta_mu,horizon):
    n=len(q);total=D(0)
    for i in range(n):
        f=deform(ref,i,q)
        for p in range(ref['offsets'][i],ref['offsets'][i+1]):
            j=ref['target'][p];xi=sub(ref['position'][j],ref['position'][i])
            residual=sub(sub(q[j],q[i]),[dot(f[3*a:3*a+3],xi) for a in range(3)])
            total+=eta_mu*ref['volume'][i]*ref['weight'][p]*ref['volume'][j]*dot(residual,residual)/(2*horizon*horizon*ref['neighbor_measure'][i])
    return total

def evaluate(ref,header,old,velocity,dt,zold,lossold,tangent,emitted=None):
    """Independent95digit point integrals and actual local transpose assembly."""
    n=33;q0=vectors(old['position']);v0=vectors(old['velocity']);v=flat_vectors(list(map(dec,velocity)))
    # Candidate positions are the actual binary64 endpoint arithmetic. Accepted
    # records additionally verify this same endpoint against the committed q.
    q1=[[dec(float(old['position'][i][a])+float(dt)*float(velocity[3*i+a])) for a in range(3)] for i in range(n)]
    p=list(map(dec,header['parameters']));eta_mu=dec(header['eta'])*p[1];horizon=dec(header['reference']['horizon'])
    mass=list(map(dec,header['physical_mass']));point=[];material_force=[[D(0)]*3 for _ in range(n)]
    jac=[[D(0)]*(3*n) for _ in range(3*n)]
    for i in range(n):
        f0=deform(ref,i,q0);f1=deform(ref,i,q1)
        response=material_step(f0,f1,zold[i],lossold[i],dt,p,tangent)
        derivative=response['endpoint_tangent'] if tangent else [[dec(x) for x in row] for row in emitted['points'][i]['endpoint_tangent']]
        point.append(response);volume=ref['volume'][i]
        for k in ref['local'][i]:
            gk=ref['shape_gradient'][i*n+k]
            for a in range(3):material_force[k][a]-=volume*dot(gk,response['average_pk1'][3*a:3*a+3])
        # Compute T*B once per local column, then apply the actual B^T.
        tb={}
        for k in ref['local'][i]:
            gk=ref['shape_gradient'][i*n+k]
            for c in range(3):tb[k,c]=[sum((derivative[3*c+d][r]*gk[d] for d in range(3)),D(0)) for r in range(9)]
        for k in ref['local'][i]:
            gk=ref['shape_gradient'][i*n+k]
            for j in ref['local'][i]:
                for a in range(3):
                    for c in range(3):jac[3*k+a][3*j+c]+=dt*dt*volume*dot(gk,tb[j,c][3*a:3*a+3])
    average=[[(x+y)/2 for x,y in zip(a,b)] for a,b in zip(q0,q1)]
    stabilization=[[-eta_mu*sum((ref['L'][i][j]*average[j][a] for j in range(n)),D(0)) for a in range(3)] for i in range(n)]
    force=[[a+b for a,b in zip(x,y)] for x,y in zip(material_force,stabilization)]
    momentum=[]
    for i in range(n):
        for a in range(3):
            row=3*i+a;momentum.append(mass[i]*(v[i][a]-v0[i][a])-dt*force[i][a]);jac[row][row]+=mass[i]
            for j in range(n):jac[row][3*j+a]+=dt*dt*eta_mu*ref['L'][i][j]/2
    stored0=sum((ref['volume'][i]*point[i]['stored_old'] for i in range(n)),D(0))
    stored1=sum((ref['volume'][i]*point[i]['stored_new'] for i in range(n)),D(0))
    loss=sum((ref['volume'][i]*point[i]['physical_dissipation'] for i in range(n)),D(0))
    work=sum((ref['volume'][i]*point[i]['work'] for i in range(n)),D(0))
    hg0=nonaffine_energy(ref,q0,eta_mu,horizon);hg1=nonaffine_energy(ref,q1,eta_mu,horizon)
    return dict(points=point,q1=q1,momentum=momentum,jacobian=jac,material_force=material_force,stabilization_force=stabilization,internal_force=force,
        material_stored_old=stored0,material_stored_new=stored1,physical_material_loss=loss,material_work=work,
        stabilization_stored_old=hg0,stabilization_stored_new=hg1,
        transpose_work_residual=-dt*sum((dot(f,u) for f,u in zip(material_force,v)),D(0))-work,
        stabilization_work_residual=-dt*sum((dot(f,u) for f,u in zip(stabilization,v)),D(0))-(hg1-hg0))

class Inspector:
    def __init__(self,case,header,full,expect_normal_only=True):
        self.case=case;self.h=header;self.full=full;self.maximum=0.;self.trial_points=0
        self.ref,self.geometry=check_geometry(case,header,expect_normal_only)
        self.mass=list(map(dec,header['physical_mass']));self.I=sum(self.mass,D(0));self.E=max(self.I,sum(self.ref['volume'],D(0))*max(map(dec,header['parameters'][:2])))
        self.C=[list(map(dec,row)) for row in header['C']]
        self.ids=header['plate_ids'];self.cache={}
    def close(self,actual,wanted,scale=D(1)):
        self.maximum=max(self.maximum,close_tree(actual,wanted,scale))
    def expected_end_time(self):return dec(self.case['dt'])*self.case['steps']
    def initial(self):
        c=self.case;h=self.h;old=h['initial_state'];ref=self.ref;n=33
        self.close(h['parameters'],c['parameters']);assert h['mode']==c['mode'] and h['steps_requested']==c['steps']
        self.close(h['nominal_dt'],c['dt']);self.close(h['end_time'],self.expected_end_time());self.close(h['stroke'],c['stroke'])
        q=[row.copy() for row in ref['position']]
        if c['mode']==0:
            f=list(map(dec,c['initial_affine']));dq=list(map(dec,c['local_perturbation']))
            for i in range(n):
                x=ref['position'][i];local=[x[1],-x[0],x[2]] if c['rotation'] else x
                y=[dot(f[3*a:3*a+3],local)+(dq[a] if i==c['perturb_node'] else 0) for a in range(3)]
                q[i]=[-y[1],y[0],y[2]] if c['rotation'] else y
        self.close(old['position'],q);assert old['time']==0 and old['step']==0
        expected_v=[D(0)]*3
        if c['mode']==2:
            v=list(map(dec,c['translation']));expected_v=[-v[1],v[0],v[2]] if c['rotation'] else v
        self.close(old['velocity'],[expected_v]*n)
        z=[];loss=[D(0)]*n;p=list(map(dec,h['parameters']))
        for i in range(n):
            f=deform(ref,i,q);e=green(f);zi=[dec(c['relaxed_fraction'])*x for x in e] if c['mode']==0 else [D(0)]*9;z.append(zi)
            self.close(old['green'][i],e);self.close(old['viscous'][i],zi);self.close(old['dissipation'][i],0)
            stress=[x/determinant(f) for x in multiply(pk1(f,sub(e,zi),p),transpose(f))]
            self.close(old['accepted_stress'][i],stress,max(p[:3]));self.close(old['density'][i],dec(h['rho0'])/determinant(f))
        return old,z,loss
    def get_evaluation(self,old,v,dt,z,loss,emitted=None):
        key=(old['step'],tuple(v),str(dt),self.full)
        if key not in self.cache:self.cache[key]=evaluate(self.ref,self.h,old,v,dt,z,loss,self.full,emitted)
        return self.cache[key]
    def evaluation(self,native,reference,old,v,dt):
        pscale=max(map(dec,self.h['parameters'][:3]));n=33
        for i,point in enumerate(reference['points']):
            got=native['points'][i]
            self.close(got['F0'],deform(self.ref,i,vectors(old['position'])));self.close(got['F1'],deform(self.ref,i,reference['q1']))
            for field in ['green','viscous','minimum_determinant']:self.close(got[field],point[field])
            for field in ['average_pk1','total_dissipation','stored_old','stored_new','physical_dissipation','work']:self.close(got[field],point[field],pscale)
            self.close(got['energy_residual'],0,pscale)
            if self.full:self.close(got['endpoint_tangent'],point['endpoint_tangent'],pscale)
        for field in ['material_stored_old','material_stored_new','physical_material_loss','material_work','stabilization_stored_old','stabilization_stored_new','transpose_work_residual','stabilization_work_residual']:
            self.close(native[field],reference[field],self.E)
        for field in ['material_force','stabilization_force','internal_force']:self.close(native[field],reference[field],self.I/dt)
        self.close(native['momentum'],reference['momentum'],self.I);self.close(native['jacobian'],reference['jacobian'],self.I)
        jac=reference['jacobian'];square=sum((x*x for row in jac for x in row),D(0))
        asym=sum(((jac[i][j]-jac[j][i])**2 for i in range(99) for j in range(99)),D(0))
        self.close(native['jacobian_asymmetry'],(asym/square).sqrt())
    def bounds(self,old,dt,offsets):
        q=flatten(vectors(old['position']))
        return [(dec(d)-dot(c,q))/dt for c,d in zip(self.C,offsets)]
    def trace(self,record,old,dt,z,loss):
        solve=record['solve'];trace=solve['trace'];assert trace and solve['converged'];n=99;m=len(self.C)
        bound=self.bounds(old,dt,record['new_offsets'])
        for index,t in enumerate(trace):
            assert t['iteration']==index
            v=list(map(dec,t['velocity']));lam=list(map(dec,t['normal_impulse']));f=list(map(dec,t['normalized_residual']))
            matrix=[list(map(dec,row)) for row in t['generalized_jacobian']]
            assert len(v)==n and len(lam)==m and len(f)==len(matrix)==n+m
            assert all(len(row)==n+m for row in matrix)
            for i in range(n):
                for s in range(m):self.close(matrix[i][n+s],-self.C[s][i])
            if self.full:
                ref=self.get_evaluation(old,t['velocity'],dt,z,loss)
                expected_f=[(ref['momentum'][i]-sum((self.C[s][i]*lam[s] for s in range(m)),D(0)))/self.I for i in range(n)]
                self.close(f[:n],expected_f)
                self.close([row[:n] for row in matrix[:n]],[[x/self.I for x in row] for row in ref['jacobian']])
                self.trial_points+=33
            for s in range(m):
                gap=(dot(self.C[s],v)-bound[s]);b=dec(t['normalized_gaps'][s]);a=lam[s]/self.I;radius=(a*a+b*b).sqrt()
                self.close(b,gap);self.close(f[n+s],radius-a-b)
                da=a/radius-1 if radius else D('.5').sqrt()-1;db=b/radius-1 if radius else D('.5').sqrt()-1
                self.close(matrix[n+s][:n],[db*x for x in self.C[s]])
                self.close(matrix[n+s][n:],[da if j==s else D(0) for j in range(m)])
            self.close(t['squared_norm'],dot(f,f));self.close(t['merit'],max(map(abs,f)))
            assert len(t['accumulation_roundoff_bounds'])==99 and max(t['accumulation_roundoff_bounds'])<=1e-12
            direction=list(map(dec,t['direction']))
            if direction:
                assert len(direction)==n+m and 0<t['linear_rank']<=n+m
                error=max(abs(dot(row,direction)+b) for row,b in zip(matrix,f))
                denominator=max(sum(map(abs,row),D(0)) for row in matrix)*max(map(abs,direction))+max(map(abs,f))
                assert error/denominator<=D('1e-12') and t['linear_residual']<=1e-12
            if index+1<len(trace):
                next_t=trace[index+1];alpha=dec(t['alpha']);assert direction and 0<alpha<=1
                self.close(next_t['velocity'],[x+alpha*d for x,d in zip(v,direction[:n])])
                self.close(next_t['normal_impulse'],[x+alpha*self.I*d for x,d in zip(lam,direction[n:])],self.I)
                nf=list(map(dec,next_t['normalized_residual']))
                assert dot(nf,nf)<=(1-D('1e-4')*alpha)*dot(f,f)+D('1e-25')
        assert trace[-1]['merit']<=1e-12 and trace[-1]['KKT']<=1e-12 and solve['normalized_KKT']<=1e-12
        self.close(solve['velocity'],trace[-1]['velocity']);self.close(solve['normal_impulse'],trace[-1]['normal_impulse'],self.I)
        assert not solve['precision_rejection']['present']
    def dense(self,q,time):
        h=self.h;n=33;minimum=[D('Infinity')]*2;normals=vectors(h['plate_normal']);u=self.displacement(time)
        particles=tuple(tuple(x) for x in h['reference']['position'])
        for name in ['surface','dense_surface']:
            map_=h[name];samples=tuple(tuple(x) for x in map_['reference'])
            expected=expected_mls(particles,samples,map_['support'])
            for s in range(len(samples)):
                x=[sum((expected['weights'][p]*q[expected['targets'][p]][a] for p in range(expected['offsets'][s],expected['offsets'][s+1])),D(0)) for a in range(3)]
                for plate in range(2):minimum[plate]=min(minimum[plate],dot(normals[plate],x)+1-(u if plate else 0))
        return minimum
    def dense_interval(self,q0,q1,t0,t1):
        h=self.h;minimum=[D('Infinity')]*2;normals=vectors(h['plate_normal'])
        theta=t0/dec(h['end_time']);delta=(t1-t0)/dec(h['end_time']);stroke=dec(h['stroke'])
        particles=tuple(tuple(x) for x in h['reference']['position'])
        for name in ['surface','dense_surface']:
            map_=h[name];samples=tuple(tuple(x) for x in map_['reference']);mls=expected_mls(particles,samples,map_['support'])
            for s in range(len(samples)):
                x=[]
                for q in [q0,q1]:
                    x.append([sum((mls['weights'][p]*q[mls['targets'][p]][a] for p in range(mls['offsets'][s],mls['offsets'][s+1])),D(0)) for a in range(3)])
                for plate in range(2):
                    coefficients=[dot(normals[plate],x[0])+1-(self.displacement(t0) if plate else 0),dot(normals[plate],sub(x[1],x[0])),D(0),D(0)]
                    if plate:
                        coefficients[1]-=stroke*(6*theta-6*theta*theta)*delta
                        coefficients[2]=-stroke*(3-6*theta)*delta*delta
                        coefficients[3]=2*stroke*delta**3
                    minimum[plate]=min(minimum[plate],cubic_minimum(coefficients))
        return minimum
    def displacement(self,time):
        s=time/dec(self.h['end_time']);return dec(self.h['stroke'])*s*s*(3-2*s)
    def step(self,row,old,z,loss):
        assert row['type']=='step' and row['committed'] and row['local_root_accepted']
        dt=dec(row['actual_dt']);new=row['state'];assert dt>0 and new['step']==old['step']+1
        self.close(row['old_time'],old['time']);self.close(new['time'],dec(old['time'])+dt);self.close(row['candidate_time'],new['time'])
        if new['step']==self.case['steps']:self.close(new['time'],self.h['end_time'])
        else:self.close(dt,self.case['dt'])
        self.close(row['upper_old'],self.displacement(dec(old['time'])));self.close(row['upper_new'],self.displacement(dec(new['time'])))
        for s,plate in enumerate(self.ids):
            self.close(row['old_offsets'][s],-1+(self.displacement(dec(old['time'])) if plate else 0))
            self.close(row['new_offsets'][s],-1+(self.displacement(dec(new['time'])) if plate else 0))
        v=row['solve']['velocity'];ref=self.get_evaluation(old,v,dt,z,loss,row['evaluation'])
        self.evaluation(row['evaluation'],ref,old,v,dt);self.trace(row,old,dt,z,loss)
        self.close(new['position'],ref['q1']);self.close(new['velocity'],flat_vectors(list(map(dec,v))))
        self.close(new['last_dt'],dt);self.close(new['last_normal_impulse'],row['solve']['normal_impulse'],self.I)
        self.close(new['last_KKT'],row['solve']['normalized_KKT']);assert new['last_newton_iterations']==row['solve']['trace'][-1]['iteration']
        p=list(map(dec,self.h['parameters']));q1=vectors(new['position']);v1=vectors(new['velocity']);q0=vectors(old['position']);v0=vectors(old['velocity'])
        for i,point in enumerate(ref['points']):
            self.close(new['green'][i],point['green']);self.close(new['viscous'][i],point['viscous']);self.close(new['dissipation'][i],point['total_dissipation'],max(p[:3]))
            f=deform(self.ref,i,q1);stress=[x/determinant(f) for x in multiply(pk1(f,sub(green(f),point['viscous']),p),transpose(f))]
            self.close(new['accepted_stress'][i],stress,max(p[:3]));self.close(new['density'][i],dec(self.h['rho0'])/determinant(f))
        lam=list(map(dec,row['solve']['normal_impulse']));bound=self.bounds(old,dt,row['new_offsets']);fv=flatten(v1);flatq0=flatten(q0)
        gaps=[dot(c,fv)-b for c,b in zip(self.C,bound)]
        assert min(lam,default=D(0))>=-D('1e-10')*self.I and min(gaps,default=D(0))>=-D('1e-10')
        assert sum((abs(l*g) for l,g in zip(lam,gaps)),D(0))<=D('1e-10')*self.E
        impulse=[sum((self.C[s][i]*lam[s] for s in range(len(lam))),D(0)) for i in range(99)]
        self.close(row['contact_impulse'],flat_vectors(impulse),self.I)
        normalized=[(r-i)/self.I for r,i in zip(ref['momentum'],impulse)];assert max(map(abs,normalized))<=D('1e-10')
        self.close(row['normalized_momentum_residual'],flat_vectors(normalized))
        plate_impulse=[[sum((sum((self.C[s][3*i+a] for i in range(33)),D(0))*lam[s] for s in range(len(lam)) if self.ids[s]==plate),D(0)) for a in range(3)] for plate in range(2)]
        self.close(row['particle_plate_impulse'],plate_impulse,self.I);self.close(row['plate_reaction_interval'],[[-x/dt for x in v] for v in plate_impulse],self.I/dt)
        assert row['force_cursor']=='impulse_over_actual_interval_dt_not_instantaneous_terminal_force'
        work=[sum((lam[s]*(dec(row['new_offsets'][s])-dec(row['old_offsets'][s]))/dt for s in range(len(lam)) if self.ids[s]==plate),D(0)) for plate in range(2)]
        normal_loss=sum((lam[s]*(dot(self.C[s],flatq0)-dec(row['old_offsets'][s]))/dt for s in range(len(lam))),D(0))
        k0=sum((self.mass[i]*dot(v0[i],v0[i])/2 for i in range(33)),D(0));k1=sum((self.mass[i]*dot(v1[i],v1[i])/2 for i in range(33)),D(0))
        algorithmic=sum((self.mass[i]*dot(sub(v1[i],v0[i]),sub(v1[i],v0[i]))/2 for i in range(33)),D(0))
        for name,wanted in [('plate_work',work),('normal_contact_loss',normal_loss),('kinetic_old',k0),('kinetic_new',k1),('algorithmic_loss',algorithmic)]:self.close(row[name],wanted,self.E)
        energy=sum(work,D(0))-(k1-k0+ref['material_stored_new']-ref['material_stored_old']+ref['stabilization_stored_new']-ref['stabilization_stored_old']+ref['physical_material_loss']+algorithmic+normal_loss)
        assert abs(energy)<=D('1e-10')*self.E;self.close(row['energy_residual'],energy,self.E);assert not row['friction_implemented']
        net=[D(0)]*3;force_sum=[D(0)]*3;angular=[D(0)]*3;torque=[D(0)]*3
        for i in range(33):
            for a in range(3):
                net[a]+=self.mass[i]*(v1[i][a]-v0[i][a])-impulse[3*i+a];force_sum[a]+=ref['internal_force'][i][a]
            external=cross(q0[i],impulse[3*i:3*i+3]);delta=sub(cross(q1[i],v1[i]),cross(q0[i],v0[i]))
            angular=add3(angular,[self.mass[i]*x-y for x,y in zip(delta,external)]);torque=add3(torque,cross(q0[i],ref['internal_force'][i]))
        self.close(row['total_momentum_residual'],net,self.I);self.close(row['internal_force_sum'],force_sum,self.I/dt)
        self.close(row['angular_drift_minus_external'],angular,self.I);self.close(row['internal_old_torque'],torque,self.I/dt)
        assert norm(net)<=D('1e-10')*self.I and norm(force_sum)<=D('1e-10')*self.I/dt
        assert norm(sub(angular,[dt*x for x in torque]))<=D('1e-10')*self.I*D('1.25')
        if self.case['mode']==1:
            assert row['geometry_gate_applicable'];dense=self.dense_interval(q0,q1,dec(old['time']),dec(new['time']))
            endpoint=self.dense(q1,dec(new['time']));self.close(row['dense_endpoint_minimum_plate_gap'],endpoint);self.close(row['dense_minimum_plate_gap'],dense)
            assert min(dense)>=-D('.01')*dec(self.h['stroke'])
        else:assert not row['geometry_gate_applicable'] and row['dense_minimum_plate_gap'] is None and row['dense_endpoint_minimum_plate_gap'] is None
        znew=[x['viscous'] for x in ref['points']];lossnew=[x['total_dissipation'] for x in ref['points']]
        return new,znew,lossnew

def inspect(case,records,code,full=False):
    with localcontext() as ctx:
        ctx.prec=95
        assert code==case.get('expected_code',0),(case['name'],'return_code',code)
        assert records and records[0]['type']=='identity'
        if code:
            assert code==1 and case['invalid_injection'] in [1,2,3,4]
            assert records[-1]['type']=='error' and not any(r['type']=='step' for r in records)
            assert records[-1]['message']==invalid_message(case['invalid_injection']),records[-1]['message']
            if case['invalid_injection']==2:
                assert records[1]['type']=='header';check_geometry(case,records[1])
                damaged=records[1]['initial_state']['green'][0]
                assert abs(sum(dec(damaged[k]) for k in [0,4,8])-D('.01'))<=D('1e-12')
            return dict(name=case['name'],status='declared_invalid_physical_state_rejected_without_step',reason=case['invalid_injection'])
        assert records[1]['type']=='header' and records[-1]['type']=='end' and records[-1]['completed']
        h=records[1];audit=Inspector(case,h,full);old,z,loss=audit.initial()
        assert h['geometry_cursor']=='minimum_over_entire_linear_particle_path_and_true_cubic_wall_interval'
        steps=[r for r in records if r['type']=='step']
        assert not any(r['type'] in ['rejected_step','error'] for r in records)
        if case['mode']==0:
            assert len(records)==4 and records[2]['type']=='operator_sample' and not records[2]['accepted_physical_step']
            row=records[2];dt=dec(row['dt']);ref=audit.get_evaluation(old,row['velocity'],dt,z,loss,row['evaluation'])
            audit.evaluation(row['evaluation'],ref,old,row['velocity'],dt);assert records[-1]['accepted_steps']==0
            velocity=flat_vectors(list(map(dec,row['velocity'])))
            audit.close(row['trial_kinetic_actual_mass'],sum((audit.mass[i]*dot(velocity[i],velocity[i])/2 for i in range(33)),D(0)),audit.E)
            assert row['evaluation']['jacobian_asymmetry']>1e-9
        else:
            assert len(steps)==case['steps']
            for row in steps:old,z,loss=audit.step(row,old,z,loss)
            assert records[-1]['accepted_steps']==case['steps'];audit.close(records[-1]['last_accepted_time'],h['end_time'])
        stability=None
        if full:stability=stability_spectrum(case['rotation'],tuple(h['reference']['volume']),h['reference']['horizon'],float(h['parameters'][0]),float(h['parameters'][1]),float(h['eta']))
        return dict(name=case['name'],status='independent_actual_PD_reference_material_global_normal_consistency_passed' if full else 'runtime_original_PD_gates_passed_pending_full_reference',
            accepted_steps=len(steps),maximum_normalized_error=max(audit.maximum,audit.geometry['maximum_normalized_cache_error']),
            geometry=audit.geometry,independent_stability=stability,all_trial_material_tangents_checked=full,
            independent_trial_material_points=audit.trial_points,normal_only=True,trajectory_convergence_qualified=False,angular_conservation_qualified=False,
            full_restart_qualified=False,S3_complete=False,calibration_allowed=False)

def read_jsonl(path):
    with open(path,encoding='utf8') as stream:return [json.loads(line) for line in stream if line.strip()]

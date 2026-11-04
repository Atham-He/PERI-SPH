"""Streaming original 3D gel/friction reference; never calls a native solver.

The default normal-only reader remains strict. This reader explicitly admits
normal_only=False, uses exact original contact coefficients/normalization, and
independently reconstructs point material R/J at every Newton observation.
"""
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import json,math
from peri_coulomb_contact_audit import frac,decimal,fdot,coefficients,check_trace,dec,close_tree
from peri_bounded_gel_audit import (
    Inspector,vectors,flat_vectors,flatten,deform,green,sub,multiply,pk1,transpose,determinant,
    cross,add3,dot,norm,cubic_minimum,expected_mls,stability_spectrum)

CURSOR='actual_interval_mean_prescribed_displacement_difference_over_actual_dt'


def protocol(case,step,time):
    """95-digit frozen piecewise motion, original binary endpoint clock."""
    dt=float(case['dt']);ends=[16*dt,24*dt,40*dt]
    phase=0 if step<16 else (1 if step<24 else 2)
    assert 0<=step<(16 if case['protocol_kind']==0 else 40)
    begin=0. if phase==0 else ends[phase-1];end=ends[phase]
    actual_dt=end-time if step in [15,23,39] else dt;next_time=time+actual_dt
    assert begin<=time<end and actual_dt>0 and next_time<=end
    duration=dec(end-begin);theta=(dec(time)-dec(begin))/duration;delta=dec(actual_dt)/duration
    def shape(t):
        if phase==1:return D(1)
        x=(dec(t)-dec(begin))/duration;s=x*x*(3-2*x)
        return s if phase==0 else 1-s
    a,b=shape(time),shape(next_time);stroke=dec(case['stroke'])
    amplitude=stroke if phase==0 else -stroke
    cubic=[-1+stroke*a,D(0),D(0),D(0)]
    if phase!=1:
        cubic[1]=amplitude*(6*theta-6*theta*theta)*delta
        cubic[2]=amplitude*(3-6*theta)*delta*delta
        cubic[3]=-2*amplitude*delta**3
    tangent=list(map(dec,case['tangent_stroke']))
    return dict(phase=phase,step=step,old_time=time,new_time=next_time,actual_dt=actual_dt,
        old_offset=[D(-1),-1+stroke*a],new_offset=[D(-1),-1+stroke*b],
        normal_cubic=[[D(-1),D(0),D(0),D(0)],cubic],
        old_tangent_displacement=[[D(0),D(0)],[x*a for x in tangent]],
        new_tangent_displacement=[[D(0),D(0)],[x*b for x in tangent]])


def neumaier(value,total,correction):
    nxt=total+value
    correction+=((total-nxt)+value) if abs(total)>=abs(value) else ((value-nxt)+total)
    return nxt,correction


def expected_contacts(header,interval):
    surface=header['surface'];rows=[]
    for s in range(len(surface['reference'])):
        weights=[0.]*33
        for k in range(surface['offsets'][s],surface['offsets'][s+1]):
            weights[surface['targets'][k]]=surface['weights'][k]
        for plate in range(2):
            rows.append(dict(weights=weights.copy(),normal=header['plate_normal'][plate],
                tangent=header['plate_tangent'],old_offset=interval['old_offset'][plate],
                new_offset=interval['new_offset'][plate],mu=header['mu'][plate],
                wall_tangent_velocity=interval['wall_tangent_velocity'][plate]))
    return rows


class FrictionInspector(Inspector):
    def __init__(self,case,header,full):
        super().__init__(case,header,full,expect_normal_only=False)
        mass=0.
        for value in header['physical_mass']:mass+=value
        assert header['normalization_mass_scale']==mass
        assert header['normalization_impulse_scale']==mass*header['solver_velocity_scale']
        self.I=dec(header['normalization_impulse_scale'])
        self.E=max(self.I*dec(header['solver_velocity_scale']),sum(self.ref['volume'],D(0))*max(map(dec,header['parameters'][:2])))
        assert header['format']==1 and header['actual_PD_body']
        assert header['physical_mass']==[header['rho0']*v for v in header['reference']['volume']]
        assert header['primal_dofs']==99 and header['contact_count']==len(header['plate_ids'])
        assert header['unknowns']==99+3*len(header['plate_ids'])
        assert header['solver_velocity_scale']==1
        expected_frame=[0,-1,0,1,0,0,0,0,1] if case['rotation'] else [1,0,0,0,1,0,0,0,1]
        expected_tangent=[[0,1,0],[0,0,1]] if case['rotation'] else [[1,0,0],[0,0,1]]
        self.close(header['frame'],expected_frame);self.close(header['plate_tangent'],expected_tangent)
        assert header['protocol_kind']==case['protocol_kind'] and header['phase_step_counts']==[16,8,16]
        self.close(header['phase_end_times'],[16*case['dt'],24*case['dt'],40*case['dt']])
        self.close(header['mu'],case['mu']);self.close(header['tangent_stroke'],case['tangent_stroke'])
        assert header['wall_velocity_cursor']==CURSOR
        assert header['geometry_cursor']=='minimum_over_each_linear_particle_and_true_prescribed_cubic_wall_interval'
        self.G=coefficients(header['initial_contacts'],99)
        assert header['C_full_three_components']==self.G and header['C']==self.G[::3]
        self.ever_active=set();self.slip_signs=[set(),set()]
        self.statistics=dict(sliding_rows=0,strict_stick_rows=0,previously_active_separation_rows=0,
            Maxwell_loss=D(0),normal_loss=D(0),friction_loss=D(0),algorithmic_loss=D(0),
            positive_input_work=D(0),minimum_observed_gap=D('Infinity'))

    def initial(self):
        old,z,loss=super().initial();history=self.h['initial_friction_history'];m=len(self.ids)
        for field in ['last_tangent_impulse','relative_tangent_travel_sum','relative_tangent_travel_correction']:
            assert history[field]==[[0,0] for _ in range(m)]
        assert history['ledger_order']=='lower_Ixyz,upper_Ixyz,lower_work,upper_work,Maxwell,normal,friction,algorithmic'
        assert history['ledger_sum']==history['ledger_correction']==[0]*12
        assert old['last_normal_impulse']==[0]*m and old['last_dt']==old['last_KKT']==old['last_newton_iterations']==0
        assert history['relative_travel_includes_open_rows']
        return old,z,loss,history

    def interval(self,interval,old):
        k=old['step'];expected=protocol(self.case,k,old['time'])
        for field in ['phase','step','old_time','new_time','actual_dt']:
            assert interval[field]==expected[field],('original_clock',field)
        for field in ['old_offset','new_offset','normal_cubic','old_tangent_displacement','new_tangent_displacement']:
            self.close(interval[field],expected[field])
        assert interval['wall_velocity_cursor']==CURSOR
        for plate in range(2):
            for a in range(2):
                want=(interval['new_tangent_displacement'][plate][a]-interval['old_tangent_displacement'][plate][a])/interval['actual_dt']
                assert interval['wall_tangent_velocity'][plate][a]==want
        return expected

    def geometry_interval(self,q0,q1,interval,expected):
        minimum=[D('Infinity')]*2;endpoint=[D('Infinity')]*2
        particles=tuple(tuple(x) for x in self.h['reference']['position'])
        normals=vectors(self.h['plate_normal'])
        for name in ['surface','dense_surface']:
            map_=self.h[name];samples=tuple(tuple(x) for x in map_['reference'])
            mls=expected_mls(particles,samples,map_['support'])
            for s in range(len(samples)):
                x=[[sum((mls['weights'][k]*q[mls['targets'][k]][a]
                    for k in range(mls['offsets'][s],mls['offsets'][s+1])),D(0)) for a in range(3)] for q in [q0,q1]]
                for plate in range(2):
                    c=expected['normal_cubic'][plate]
                    gap=[dot(normals[plate],x[0])-c[0],dot(normals[plate],sub(x[1],x[0]))-c[1],-c[2],-c[3]]
                    minimum[plate]=min(minimum[plate],cubic_minimum(gap))
                    endpoint[plate]=min(endpoint[plate],dot(normals[plate],x[1])-expected['new_offset'][plate])
        return minimum,endpoint

    def step(self,row,old,z,loss,history):
        assert row['type']=='step' and row['committed'] and row['local_step_gates_passed'] and row['candidate_ready']
        assert row['reason']=='all_applicable_step_gates_passed'
        interval=row['interval'];expected=self.interval(interval,old);dt=dec(interval['actual_dt']);new=row['state']
        assert new['time']==interval['new_time'] and new['step']==old['step']+1
        rows=expected_contacts(self.h,interval);assert rows==row['contacts']
        if old['step']==0:assert rows==self.h['initial_contacts']
        assert coefficients(rows,99)==self.G
        def material(v,emitted):
            if not self.full:return emitted['supplied_momentum'],emitted['supplied_jacobian']
            reference=self.get_evaluation(old,v,dt,z,loss)
            self.trial_points+=33
            return reference['momentum'],reference['jacobian']
        self.maximum=max(self.maximum,check_trace(flatten(old['position']),rows,interval['actual_dt'],
            self.h['normalization_impulse_scale'],self.h['solver_velocity_scale'],row['solve'],material,self.G))
        v=row['solve']['velocity'];p=row['solve']['impulses']
        reference=self.get_evaluation(old,v,dt,z,loss,row['evaluation'])
        self.evaluation(row['evaluation'],reference,old,v,dt)
        last=row['solve']['trace'][-1]['original']
        assert row['evaluation']['momentum']==last['supplied_momentum']
        assert row['evaluation']['jacobian']==last['supplied_jacobian']
        self.close(new['position'],reference['q1']);self.close(new['velocity'],flat_vectors(list(map(dec,v))))
        assert new['last_dt']==interval['actual_dt'] and new['last_KKT']==row['solve']['normalized_KKT']
        assert new['last_newton_iterations']==row['solve']['trace'][-1]['iteration']
        assert new['last_normal_impulse']==p[::3]
        hnew=row['friction_history'];assert hnew['last_tangent_impulse']==[p[3*s+1:3*s+3] for s in range(len(rows))]
        for field in ['ledger_order','relative_travel_includes_open_rows']:assert hnew[field]==history[field]
        parameters=list(map(dec,self.h['parameters']));q0=vectors(old['position']);q1=vectors(new['position'])
        v0=vectors(old['velocity']);v1=vectors(new['velocity'])
        for i,point in enumerate(reference['points']):
            for field,wanted in [('green',point['green']),('viscous',point['viscous']),('dissipation',point['total_dissipation'])]:
                self.close(new[field][i],wanted,max(parameters[:3]) if field=='dissipation' else D(1))
            f=deform(self.ref,i,q1);stress=[x/determinant(f) for x in multiply(pk1(f,sub(green(f),point['viscous']),parameters),transpose(f))]
            self.close(new['accepted_stress'][i],stress,max(parameters[:3]));self.close(new['density'][i],dec(self.h['rho0'])/determinant(f))
        fp=list(map(frac,p));fv=list(map(frac,v));flatq=flatten(old['position']);qdt=frac(interval['actual_dt'])
        impulse=[fdot([r[i] for r in self.G],p) for i in range(99)]
        assert flatten(row['contact_impulse'])==[float(x) for x in impulse]
        original_balance=[(frac(a)-b)/frac(self.h['normalization_impulse_scale']) for a,b in zip(row['evaluation']['momentum'],impulse)]
        assert flatten(row['normalized_momentum_residual'])==[float(x) for x in original_balance]
        true_balance=[(r-decimal(p))/self.I for r,p in zip(reference['momentum'],impulse)]
        assert max(map(abs,true_balance))<=D('1e-10')
        plate_n=[[F(0)]*3 for _ in range(2)];plate_t=[[F(0)]*3 for _ in range(2)]
        Wn=[F(0)]*2;Wt=[F(0)]*2;Dn=F(0);Df=F(0)
        travel=[row.copy() for row in history['relative_tangent_travel_sum']]
        travel_correction=[row.copy() for row in history['relative_tangent_travel_correction']]
        for s,contact in enumerate(rows):
            plate=self.ids[s];lam=fp[3*s]
            for a in range(3):
                plate_n[plate][a]+=sum((frac(self.G[3*s][3*i+a])*lam for i in range(33)),F(0))
                plate_t[plate][a]+=sum((frac(self.G[3*s+b][3*i+a])*fp[3*s+b] for i in range(33) for b in [1,2]),F(0))
            Wn[plate]+=lam*(frac(contact['new_offset'])-frac(contact['old_offset']))/qdt
            Dn+=lam*(fdot(self.G[3*s],flatq)-frac(contact['old_offset']))/qdt
            relative=[fdot(self.G[3*s+b],fv)-frac(contact['wall_tangent_velocity'][b-1]) for b in [1,2]]
            Wt[plate]+=fdot(fp[3*s+1:3*s+3],contact['wall_tangent_velocity']);Df-=fdot(fp[3*s+1:3*s+3],relative)
            for a in range(2):
                travel[s][a],travel_correction[s][a]=neumaier(float(qdt*relative[a]),travel[s][a],travel_correction[s][a])
            active=lam>frac(self.I)*frac(1e-8)
            if active:
                self.ever_active.add(s)
                radius=dec(contact['mu'])*decimal(lam);tn=sum(decimal(x)**2 for x in fp[3*s+1:3*s+3]).sqrt()
                speed=sum(decimal(x)**2 for x in relative).sqrt()
                if speed>D('1e-8'):
                    self.statistics['sliding_rows']+=1
                    for a in range(2):
                        if abs(relative[a])>frac(1e-8):self.slip_signs[a].add(1 if relative[a]>0 else -1)
                elif radius-tn>D('1e-8')*self.I:self.statistics['strict_stick_rows']+=1
            else:
                endpoint_gap=fdot(self.G[3*s],flatten(new['position']))-frac(contact['new_offset'])
                if s in self.ever_active and endpoint_gap>frac(1e-8):self.statistics['previously_active_separation_rows']+=1
        assert hnew['relative_tangent_travel_sum']==travel and hnew['relative_tangent_travel_correction']==travel_correction
        total=[[a+b for a,b in zip(x,y)] for x,y in zip(plate_n,plate_t)]
        for field,wanted in [('normal_plate_impulse',plate_n),('tangent_plate_impulse',plate_t),('particle_plate_impulse',total)]:
            assert row[field]==[[float(x) for x in v] for v in wanted]
        self.close(row['plate_reaction_interval'],[[-decimal(x)/dt for x in v] for v in total],self.I/dt)
        assert row['force_cursor']=='impulse_over_actual_interval_dt_not_instantaneous_terminal_force'
        k0=sum((frac(self.h['physical_mass'][i])*fdot(old['velocity'][i],old['velocity'][i])/2 for i in range(33)),F(0))
        k1=sum((frac(self.h['physical_mass'][i])*fdot(new['velocity'][i],new['velocity'][i])/2 for i in range(33)),F(0))
        da=sum((frac(self.h['physical_mass'][i])*sum(((frac(a)-frac(b))**2 for a,b in zip(new['velocity'][i],old['velocity'][i])),F(0))/2 for i in range(33)),F(0))
        for field,wanted in [('kinetic_old',k0),('kinetic_new',k1),('algorithmic_loss',da),('normal_contact_loss',Dn),('friction_loss',Df)]:
            assert row[field]==float(wanted)
        assert row['normal_plate_work']==list(map(float,Wn)) and row['tangent_plate_work']==list(map(float,Wt))
        energy=decimal(sum(Wn)+sum(Wt))-(decimal(k1-k0)+reference['material_stored_new']-reference['material_stored_old']+
            reference['stabilization_stored_new']-reference['stabilization_stored_old']+reference['physical_material_loss']+decimal(da+Dn+Df))
        assert abs(energy)<=D('1e-10')*self.E and min(decimal(Dn),decimal(Df),decimal(da))>=-D('1e-10')*self.E
        self.close(row['energy_residual'],energy,self.E)
        net=[D(0)]*3;force_sum=[D(0)]*3;angular=[D(0)]*3;torque=[D(0)]*3
        for i in range(33):
            for a in range(3):
                net[a]+=self.mass[i]*(v1[i][a]-v0[i][a])-decimal(impulse[3*i+a])
                force_sum[a]+=reference['internal_force'][i][a]
            angular=add3(angular,sub([self.mass[i]*x for x in sub(cross(q1[i],v1[i]),cross(q0[i],v0[i]))],
                                   cross(q0[i],list(map(decimal,impulse[3*i:3*i+3])))))
            torque=add3(torque,cross(q0[i],reference['internal_force'][i]))
        for field,wanted in [('total_momentum_residual',net),('internal_force_sum',force_sum),
                             ('angular_drift_minus_external',angular),('internal_old_torque',torque)]:
            self.close(row[field],wanted,self.I/dt if field in ['internal_force_sum','internal_old_torque'] else self.I)
        assert norm(net)<=D('1e-10')*self.I and norm(force_sum)<=D('1e-10')*self.I/dt
        L=max(D('1.25'),max(map(norm,q0)))
        assert norm(sub(angular,[dt*x for x in torque]))<=D('1e-10')*self.I*L
        dense,endpoint=self.geometry_interval(q0,q1,interval,expected)
        assert row['geometry_gate_applicable']
        self.close(row['dense_minimum_plate_gap'],dense);self.close(row['dense_endpoint_minimum_plate_gap'],endpoint)
        assert min(dense)>=-D('.01')*dec(self.h['stroke'])
        self.statistics['minimum_observed_gap']=min(self.statistics['minimum_observed_gap'],*dense)
        sums=history['ledger_sum'].copy();corrections=history['ledger_correction'].copy()
        increments=flatten(row['particle_plate_impulse'])+[row['normal_plate_work'][k]+row['tangent_plate_work'][k] for k in range(2)]
        increments+=[row['evaluation']['physical_material_loss'],row['normal_contact_loss'],row['friction_loss'],row['algorithmic_loss']]
        for i,value in enumerate(increments):sums[i],corrections[i]=neumaier(value,sums[i],corrections[i])
        assert hnew['ledger_sum']==sums and hnew['ledger_correction']==corrections
        for field,value in [('Maxwell_loss',reference['physical_material_loss']),('normal_loss',decimal(Dn)),('friction_loss',decimal(Df)),('algorithmic_loss',decimal(da))]:
            self.statistics[field]+=value
        self.statistics['positive_input_work']+=sum((max(D(0),decimal(x+y)) for x,y in zip(Wn,Wt)),D(0))
        self.cache.clear()  # retain at most one accepted step's Newton material packets
        return new,[x['viscous'] for x in reference['points']],[x['total_dissipation'] for x in reference['points']],hnew


def inspect_stream(case,path,code,full=False):
    with localcontext() as ctx:
        ctx.prec=95
        assert code==case.get('expected_code',0)==0,(case['name'],'return_code',code)
        with open(path,encoding='utf8') as stream:
            meta=json.loads(next(stream));h=json.loads(next(stream))
            assert meta['type']=='identity' and h['type']=='header'
            audit=FrictionInspector(case,h,full);old,z,loss,history=audit.initial();count=0;end=None
            for line in stream:
                row=json.loads(line)
                if row['type']=='end':
                    assert end is None;end=row;break
                assert row['type']=='step'
                old,z,loss,history=audit.step(row,old,z,loss,history);count+=1
            assert end and end['completed'] and count==case['steps']==end['accepted_steps']
            assert end['last_accepted_time']==old['time']==h['end_time']
            assert not any(line.strip() for line in stream)
        stat=audit.statistics
        for field in case.get('require_statistics',[]):assert stat[field]>0,(case['name'],'missing_observation',field)
        if case.get('require_two_way_slip'):
            assert audit.slip_signs==[{1,-1},{1,-1}],(case['name'],'two_way_actual_relative_slip',audit.slip_signs)
        spectrum=stability_spectrum(case['rotation'],tuple(h['reference']['volume']),h['reference']['horizon'],
            *[float(h['parameters'][k]) for k in [0,1]],float(h['eta'])) if full else None
        return dict(name=case['name'],status='independent_actual_PD_original_friction_consistency_passed' if full else 'runtime_pending_full_point_tangent_review',
            accepted_steps=count,maximum_normalized_error=max(audit.maximum,audit.geometry['maximum_normalized_cache_error']),
            geometry=audit.geometry,independent_stability=spectrum,all_trial_point_tangents_checked=full,
            independent_trial_material_points=audit.trial_points,
            observations={k:float(v) if isinstance(v,D) else v for k,v in stat.items()},
            algorithmic_loss_positive_input=float(stat['algorithmic_loss']/stat['positive_input_work']) if stat['positive_input_work'] else None,
            normal_only=False,trajectory_convergence_qualified=False,angular_conservation_qualified=False,
            full_restart_qualified=False,S3_complete=False,calibration_allowed=False)

"""Frozen inputs, durable-cursor inspection and new cross-case comparisons."""
import json,math
from decimal import localcontext,Decimal as D
from pathlib import Path
from peri_elastic_normal import read,record,verify
from peri_coulomb_contact_audit import close_tree


def controls():
    def base(name,kind=0):
        return dict(name=name,mode=1,steps=16 if kind==0 else 40,protocol_kind=kind,dt=.01,stroke=.01,
            rotation=0,azimuths=8,rho0=1.,parameters=[5.,1.,.5,.3],eta=.1,support=1.25,cap_width=1/64,
            tangent_stroke=[.003,.002],mu=[.2,.2],expected_code=0,
            native_wall_seconds=300 if kind==0 else 900,output_GiB=4,
            reader_CPU_seconds=7200 if kind==0 else 14400,
            require_statistics=['nonzero_tangent0_rows','nonzero_tangent1_rows'])
    press=base('actual_PD_two_tangent_press')
    rotated=base('actual_PD_two_tangent_rotated');rotated['rotation']=1
    reverse=base('actual_PD_hold_release_reverse',1)
    reverse.update(require_two_way_slip=True,require_statistics=['nonzero_tangent0_rows','nonzero_tangent1_rows','sliding_rows','previously_active_separation_rows'])
    stick=base('actual_PD_strict_stick_hold_release',1)
    stick.update(mu=[2.,2.],tangent_stroke=[.0001,.00005],require_statistics=['nonzero_tangent0_rows','nonzero_tangent1_rows','strict_stick_rows'])
    zero=base('same114sample_zero_mu')
    zero.update(azimuths=16,mu=[0.,0.],tangent_stroke=[0.,0.],native_wall_seconds=1500,require_statistics=[])
    return [press,rotated,reverse,stick,zero]


def encode(c):
    return ' '.join(map(str,[c['protocol_kind'],c['rotation'],c['azimuths'],c['dt'],c['stroke'],*c['tangent_stroke'],
        *c['mu'],c['rho0'],*c['parameters'],c['eta'],c['support'],c['cap_width']]))+'\n'


def complete_cursor(stdout,stderr,return_code=0):
    """A regular-file record is not durable if the native fsync failed."""
    complete_steps=0;complete_time=0.;last_type=None;error_cursor=None;meta=None;malformed=False
    with Path(stdout).open(encoding='utf8') as stream:
        for line in stream:
            try:r=json.loads(line)
            except (ValueError,TypeError):malformed=True;break
            last_type=r.get('type')
            if meta is None:meta=r
            if r.get('type')=='step' and r.get('committed'):
                complete_steps=r['state']['step'];complete_time=r['state']['time']
            if 'last_durable_accepted_steps' in r:
                error_cursor=(r['last_durable_accepted_steps'],r['last_durable_accepted_time'])
    # stderr is written first on exception, even if stdout cannot flush anymore.
    with Path(stderr).open(encoding='utf8') as stream:
        for line in stream:
            try:r=json.loads(line)
            except (ValueError,TypeError):continue
            if 'last_durable_accepted_steps' in r:error_cursor=(r['last_durable_accepted_steps'],r['last_durable_accepted_time'])
    if error_cursor:
        assert 0<=error_cursor[0]<=complete_steps and 0<=error_cursor[1]<=complete_time
    # A hard signal may land after the newline but before fsync. Only a normal
    # completed native exit, a known explicit gate exit, or the native exception
    # witness proves the durable cursor. Preserve complete-looking rows
    # separately when this evidence is absent.
    confirmed=error_cursor is not None or (return_code in [0,2,3] and last_type=='end' and not malformed)
    steps,time=(error_cursor if error_cursor is not None else (complete_steps,complete_time)) if confirmed else (None,None)
    return dict(complete_looking_steps=complete_steps,complete_looking_time=complete_time,
        last_durable_accepted_steps=steps,last_durable_accepted_time=time,native_error_cursor_used=error_cursor is not None,
        durability_confirmed=confirmed,last_complete_record_type=last_type,partial_record_present=malformed),meta


def anchored_normal(root):
    review=read(root/'independent_review.json');terminal=read(root/'terminal.json');prep=read(root/'preparation.json')
    assert review['status']=='actual_bounded_PD_normal_reference_material_operator_trajectory_independently_passed'
    assert review['normal_only'] and not review['S3_complete']
    a=next(a for a in terminal['attempts'] if a['name']=='normal_compression')
    assert a['return_code']==0
    artifacts=[record(root/name) for name in ['independent_review.json','state.json','terminal.json','preparation.json','plan.json']]
    artifacts.extend([a['stdout'],a['stderr']])
    # Only immutable historical snapshots, never demand that live current
    # source still equals an old qualified source version.
    artifacts.extend(e for e in prep['artifacts'] if e['path'].startswith('/mnt/f/PERI-SPH-verification/'))
    verify(artifacts)
    return dict(review=artifacts[0],stdout=a['stdout'],stderr=a['stderr'],artifacts=artifacts,
        old_qualified_native_command=a['command'],old_reviewer_not_rerun=True)


def step_stream(path):
    with Path(path).open(encoding='utf8') as stream:
        for line in stream:
            r=json.loads(line)
            if r['type']=='header':yield r
            elif r['type']=='step':yield r


def rotate_vector(v):return [-v[1],v[0],v[2]]
def rotate_tensor(t):
    # R T R^T for the same90degree z rotation used by the frozen scene.
    order=[1,0,2];sign=[-1,1,1]
    return [sign[i]*sign[j]*t[3*order[i]+order[j]] for i in range(3) for j in range(3)]


def compare_rotated(base_path,rotated_path):
    with localcontext() as ctx:
        ctx.prec=95;maximum=0.;a=step_stream(base_path);b=step_stream(rotated_path)
        ha,hb=next(a),next(b);assert ha['steps_requested']==hb['steps_requested']==16
        I=D.from_float(ha['normalization_impulse_scale'])
        E=max(I,sum(map(D.from_float,ha['reference']['volume']))*max(map(D.from_float,ha['parameters'][:2])))
        P=max(map(D.from_float,ha['parameters'][:3]))
        count=0
        for x,y in zip(a,b,strict=True):
            assert x['state']['time']==y['state']['time'] and x['interval']==y['interval']
            for field in ['position','velocity']:
                maximum=max(maximum,close_tree(y['state'][field],[rotate_vector(v) for v in x['state'][field]]))
            for field in ['green','viscous','accepted_stress']:
                maximum=max(maximum,close_tree(y['state'][field],[rotate_tensor(v) for v in x['state'][field]],P if field=='accepted_stress' else D(1)))
            for field in ['density','dissipation']:maximum=max(maximum,close_tree(y['state'][field],x['state'][field],P if field=='dissipation' else D.from_float(ha['rho0'])))
            for field in ['contact_impulse','particle_plate_impulse','normal_plate_impulse','tangent_plate_impulse','plate_reaction_interval']:
                maximum=max(maximum,close_tree(y[field],[rotate_vector(v) for v in x[field]],I/D.from_float(x['interval']['actual_dt']) if field=='plate_reaction_interval' else I))
            for field in ['normal_plate_work','tangent_plate_work','kinetic_old','kinetic_new','algorithmic_loss','normal_contact_loss','friction_loss','energy_residual']:
                maximum=max(maximum,close_tree(y[field],x[field],E))
            count+=1
        assert count==16
        return dict(status='new_actual_PD_frame_covariance_passed',steps=count,maximum_normalized_error=maximum,individual_dual_equality_required=False)


def compare_normal_anchor(new_path,anchor):
    verify([anchor['review'],anchor['stdout'],anchor['stderr']])
    with localcontext() as ctx:
        ctx.prec=95;maximum=0.;a=step_stream(anchor['stdout']['path']);b=step_stream(new_path)
        ha,hb=next(a),next(b)
        for field in ['reference','surface','dense_surface','C','plate_ids','parameters','rho0','eta','physical_mass','stroke','nominal_dt','end_time']:
            assert ha[field]==hb[field],('same_discrete_normal_model',field)
        assert len(hb['surface']['reference'])==114 and len(hb['plate_ids'])==228 and hb['unknowns']==783
        assert hb['mu']==[0,0] and hb['tangent_stroke']==[0,0]
        I=D.from_float(hb['normalization_impulse_scale'])
        E=max(I,sum(map(D.from_float,hb['reference']['volume']))*max(map(D.from_float,hb['parameters'][:2])))
        P=max(map(D.from_float,hb['parameters'][:3]))
        count=0
        for old,new in zip(a,b,strict=True):
            assert old['actual_dt']==new['interval']['actual_dt']
            assert old['state']['time']==new['state']['time']
            for field in ['position','velocity','green','viscous','accepted_stress','density','dissipation']:
                scale=P if field in ['accepted_stress','dissipation'] else (D.from_float(hb['rho0']) if field=='density' else D(1))
                maximum=max(maximum,close_tree(new['state'][field],old['state'][field],scale))
            for field in ['contact_impulse','particle_plate_impulse','plate_reaction_interval','kinetic_new','kinetic_old','algorithmic_loss','normal_contact_loss']:
                scale=I/D.from_float(old['actual_dt']) if field=='plate_reaction_interval' else (I if field in ['contact_impulse','particle_plate_impulse'] else E)
                maximum=max(maximum,close_tree(new[field],old[field],scale))
            for field in ['material_stored_old','material_stored_new','physical_material_loss','material_work','stabilization_stored_old','stabilization_stored_new']:
                maximum=max(maximum,close_tree(new['evaluation'][field],old['evaluation'][field],E))
            assert max(abs(D.from_float(p)) for p in new['solve']['impulses'][1::3]+new['solve']['impulses'][2::3])<=D('1e-10')*I
            assert abs(D.from_float(new['friction_loss']))<=D('1e-10')*E
            count+=1
        assert count==16
        return dict(status='new_zero_mu_same114sample_normal_limit_passed',steps=count,
            maximum_normalized_error=maximum,old_qualified_review_imported=anchor['review'],old_solver_dispatched=False,old_review_repeated=False)

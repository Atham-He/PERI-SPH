"""Explicit v2 complete physical bytes/clock/whole-trial comparisons; no solver."""
import hashlib,json,math,re,struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import body_fields,mechanical_fields,loads,packed
from peri_physical_time_friction_audit import schedule_grid
from peri_actual_checkpoint_support import record,file_tree,assert_preserved,restored_noadvance_metadata

def header_body_fields(header):
    final={'state':header['initial_state']}
    friction=header['initial_friction_history']
    fields = {}

    def add(name, value, tag=3, shape=()):
        assert name not in fields
        def flatten(x):
            if isinstance(x, list):
                for y in x:
                    yield from flatten(y)
            else:
                yield x
        values = list(flatten(value))
        assert len(values) == math.prod(shape)
        if tag == 3:
            assert all(math.isfinite(x) for x in values)
        else:
            assert all(isinstance(x, int) and 0 <= x < (256 if tag == 1 else 2**64) for x in values)
        fields[name] = dict(tag=tag, shape=list(shape), raw=struct.pack('<'+{1:'B', 2:'Q', 3:'d'}[tag]*len(values), *values))

    def array(name, value, tag=3, element=()):
        add(name+'/count', len(value), 2)
        add(name, value, tag, (len(value),)+element)

    ref = header['reference']
    add('model/PDReferenceGeneration', ref['generation'], 2)
    add('model/PDReferenceHorizon', ref['horizon'])
    for key, label, tag, element in [
        ('position','PDReferencePosition',3,(3,)), ('volume','PDReferenceVolume',3,()),
        ('offsets','PDReferenceOffsets',2,()), ('target','PDReferenceTarget',2,()),
        ('weight','PDReferenceWeight',3,()), ('active','PDReferenceActive',1,()),
        ('inverse_shape','PDShapeInverse',3,(3,3)), ('shape_gradient','PDShapeGradient',3,(3,)),
        ('neighbor_measure','PDNeighborMeasure',3,()), ('shape_eigenvalue_ratio','PDShapeEigenvalueRatio',3,()),
        ('nonaffine_unit','PDNonaffineUnitCache',3,())]:
        array('model/'+label, ref[key], tag, element)
    array('model/PhysicalNodeMass', header['physical_mass'])
    add('model/ReferenceDensity', header['rho0'])
    add('model/NonaffineEta', header['eta'])
    for label, value in zip(['GelBulk','GelEquilibriumShear','GelBranchShear','GelRelaxationTime'], header['parameters'], strict=True):
        add('model/'+label, value)
    for prefix in ['surface','dense_surface']:
        surface = header[prefix]
        add(prefix+'/SurfaceGeneration', surface['generation'], 2)
        add(prefix+'/SurfaceParticleCount', surface['particle_count'], 2)
        add(prefix+'/SurfaceReferenceSupport', surface['support'])
        for key, label, tag, element in [('offsets','SurfaceOffsets',2,()), ('targets','SurfaceTargets',2,()),
                ('weights','SurfaceWeights',3,()), ('reference','SurfaceReference',3,(3,)), ('areas','SurfaceAreas',3,())]:
            array(prefix+'/'+label, surface[key], tag, element)
    state = final['state']
    for key, label, element in [('position','GelPosition',(3,)), ('velocity','GelVelocity',(3,)),
            ('green','GelGreen',(3,3)), ('viscous','GelViscous',(3,3)),
            ('dissipation','GelPhysicalDissipationDensity',()), ('accepted_stress','GelAcceptedCauchyStress',(3,3)),
            ('density','GelAcceptedDensity',()), ('last_normal_impulse','LastNormalImpulse',())]:
        array('accepted/'+label, state[key], 3, element)
    for key, label, tag in [('time','AcceptedTime',3), ('step','AcceptedStep',2),
            ('last_dt','LastActualDt',3), ('last_KKT','LastKKT',3), ('last_newton_iterations','LastNewtonIterations',2)]:
        add('accepted/'+label, state[key], tag)
    for key, label, element in [('last_tangent_impulse','LastTwoTangentImpulse',(2,)),
            ('relative_tangent_travel_sum','TangentSlipSum',(2,)),
            ('relative_tangent_travel_correction','TangentSlipCorrection',(2,)),
            ('ledger_sum','ContactCumulativeLedgerSum',()), ('ledger_correction','ContactCumulativeLedgerCorrection',())]:
        array('accepted/'+label, friction[key], 3, element)
    add('contact/Frame', header['frame'], 3, (3,3))
    for i, label in enumerate(['LowerNormal','UpperNormal']):
        add('contact/'+label, header['plate_normal'][i], 3, (3,))
    for i in range(2):
        add('contact/Tangent'+str(i), header['plate_tangent'][i], 3, (3,))
    add('contact/LowerMu', header['mu'][0])
    add('contact/UpperMu', header['mu'][1])
    array('contact/PlateIDs', header['plate_ids'], 2)
    return fields

def physical_fields(body,denominator,row=None):
    fields=mechanical_fields(body,row)
    for k,v in [('RecordSchema',2),('TimeDenominator',denominator),
                ('OldTick',row['interval']['old_tick'] if row else 0),
                ('AcceptedTick',row['interval']['new_tick'] if row else 0)]:fields[k]=packed(2,(),v)
    assert len(fields)==142
    return fields

def initial_contract(root,case,expected_trusted):
    root=Path(root);generation,fields=parse((root/'initial_candidate.bin').read_bytes())
    assert generation==1 and len(fields)==149
    contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    assert contract['Schema']==packed(2,(),2) and contract['PhysicalStreamSchema']==packed(2,(),2)
    p=case['physical_time'];grid=schedule_grid(case)
    for name,v in [('TimeDenominator',p['denominator']),('NominalTicks',p['nominal_ticks'])]:assert contract[name]==packed(2,(),v)
    for k in range(3):
        assert contract['NormalEvent'+str(k)]==packed(2,(),p['normal_event_ticks'][k])
        assert contract['TangentEvent'+str(k)]==packed(2,(),p['tangent_event_ticks'][k])
    assert contract['CheckpointSteps']==packed(2,(len(case['checkpoint_steps']),),case['checkpoint_steps'])
    trusted=hashlib.sha256(encode(contract,1)).digest()
    assert trusted==expected_trusted,'externally anchored contract required, never trust self-derived digest alone'
    body={k.removeprefix('body/'):v for k,v in fields.items() if k.startswith('body/')};assert len(body)==98
    header=b'PERIPHY2'+struct.pack('<Q',2)+trusted
    initial_record=encode(physical_fields(body,p['denominator']),1)
    physical=header+struct.pack('<Q',len(initial_record))+initial_record
    assert (root/'initial.phy').read_bytes()==physical
    assert fields['output/HeaderSHA256']==packed(1,(32,),list(hashlib.sha256(header).digest()))
    assert fields['output/PrefixSHA256']==packed(1,(32,),list(hashlib.sha256(physical).digest()))
    assert fields['accepted_algorithm/AcceptedTick']==packed(2,(),0)
    return dict(fields=fields,body=body,trusted=trusted,header=header,physical=physical,grid=grid)

def trace_packets(path,case,first,last,*,resume_generation=None):
    grid=schedule_grid(case);assert 1<=first<=last<=len(grid)-1==case['steps'];den=case['physical_time']['denominator']
    with Path(path).open('rb') as stream:
        def line():
            raw=next(stream);assert raw.endswith(b'\n'),'partial native packet';return loads(raw),raw
        meta,raw=line();assert meta['type']=='identity';yield meta,raw
        h,raw=line()
        if first==1:
            assert resume_generation is None and h['type']=='header' and h['format']==2
            assert h['initial_state']['time']==0 and h['time_denominator']==den
        else:
            assert resume_generation is not None and h['type']=='resume_header' and h['new_initial_state'] is False
            assert h['selected_generation']==resume_generation and h['accepted_step']==first-1
            assert struct.pack('<d',h['accepted_time'])==struct.pack('<d',grid[first-1]/den)
        yield h,raw
        for step in range(first,last+1):
            row,raw=line();assert row['type']=='step' and row['committed'] and row['local_step_gates_passed'] and row['candidate_ready']
            assert row['reason']=='all_applicable_step_gates_passed'
            i=row['interval'];state=row['state'];old_time,new_time=grid[step-1]/den,grid[step]/den
            assert i['old_tick']==grid[step-1] and i['new_tick']==grid[step] and i['time_denominator']==den
            assert i['step']==step-1 and state['step']==step
            for actual,want in [(i['old_time'],old_time),(i['new_time'],new_time),(state['time'],new_time),
                                (i['actual_dt'],new_time-old_time),(state['last_dt'],new_time-old_time)]:
                assert struct.pack('<d',actual)==struct.pack('<d',want),'exact rational/binary native clock'
            yield row,raw
        end,raw=line();assert end['type']=='end' and end['completed']
        assert end['accepted_steps']==last and end['last_accepted_time']==grid[last]/den and end['actual_solver_calls']==last-first+1
        yield end,raw
        assert stream.read()==b'','extra native tail records'

def trace_index(path,case,first,last,*,resume_generation=None):
    rows=[];header=None;identity=None;header_sha256=None;header_bytes=0
    for row,raw in trace_packets(path,case,first,last,resume_generation=resume_generation):
        if row['type']=='identity':identity=row
        elif row['type'] in ['header','resume_header']:
            header=row;header_sha256=hashlib.sha256(raw).hexdigest();header_bytes=len(raw)
        elif row['type']=='step':rows.append(dict(step=row['state']['step'],time=row['state']['time'],bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return dict(raw=record(path),identity=identity,header=header,header_sha256=header_sha256,header_bytes=header_bytes,
                steps=rows,first_step=first,last_step=last,scope='ENTIRE native line incl ALL Newton/trial/pointT;clock explicitly v2')

def compare_new_branch(qualified_continuous,new_branch,case):
    assert qualified_continuous['first_step']==1 and qualified_continuous['last_step']==case['steps']
    first,last=new_branch['first_step'],new_branch['last_step'];assert 1<=first<=last<=case['steps']
    assert len(new_branch['steps'])==last-first+1
    if first==1:
        assert (new_branch['header_sha256'],new_branch['header_bytes'])==(qualified_continuous['header_sha256'],qualified_continuous['header_bytes'])
    else:
        previous=qualified_continuous['steps'][first-2]
        assert new_branch['header']['accepted_step']==previous['step']
        assert struct.pack('<d',new_branch['header']['accepted_time'])==struct.pack('<d',previous['time'])
    old={r['step']:r for r in qualified_continuous['steps']};different=[]
    for row in new_branch['steps']:
        before=old[row['step']];assert struct.pack('<d',row['time'])==struct.pack('<d',before['time'])
        if (row['bytes'],row['sha256'])!=(before['bytes'],before['sha256']):different.append(row['step'])
    return dict(first_step=first,last_step=last,compared_complete_steps=last-first+1,different_steps=different,
                original_full_reference_not_repeated=not different,full_restart_qualified=False)

def candidate_fields(initial,body,prefix,case,step,time):
    grid=initial['grid'];p=case['physical_time'];tick=grid[step]
    assert struct.pack('<d',time)==struct.pack('<d',tick/p['denominator'])
    fields=dict(initial['fields']);fields.update({'body/'+k:v for k,v in body.items()})
    def put(name,value,tag=2,shape=()):
        assert fields[name]['tag']==tag and fields[name]['shape']==list(shape);fields[name]=packed(tag,shape,value)
    phase=3 if step==case['steps'] else next(k for k,e in enumerate(p['normal_event_ticks']) if tick<e)
    put('accepted_algorithm/NextPhase',phase);put('accepted_algorithm/AcceptedTick',tick)
    for name,value in [('Records',step+1),('Bytes',len(prefix)),('LastStep',step)]:put('output/'+name,value)
    put('output/LastTime',time,3);put('output/PrefixSHA256',list(hashlib.sha256(prefix).digest()),1,(32,))
    assert len(fields)==149
    return fields

def derive_continuous(root,path,case,expected_trusted):
    initial=initial_contract(root,case,expected_trusted);prefix=initial['physical'];checkpoints={};index=[];header=None;header_sha256=None;header_bytes=0
    for row,raw in trace_packets(path,case,1,case['steps']):
        if row['type']=='header':
            header=row;header_sha256=hashlib.sha256(raw).hexdigest();header_bytes=len(raw)
            assert header_body_fields(header)==initial['body'],'all98 actual header body fields must match authenticated initial CP bytes'
        elif row['type']=='step':
            body=body_fields(initial['body'],row['state'],row['friction_history']);step=row['state']['step']
            physical=encode(physical_fields(body,case['physical_time']['denominator'],row),step+1)
            prefix+=struct.pack('<Q',len(physical))+physical
            index.append(dict(step=step,time=row['state']['time'],bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
            if step in case['checkpoint_steps']:
                fields=candidate_fields(initial,body,prefix,case,step,row['state']['time'])
                checkpoints[step]=dict(archive=encode(fields,step),physical=prefix,body=body,
                                       full_native_state=row['state'],full_native_friction_history=row['friction_history'])
    assert len(index)==case['steps'] and set(checkpoints)==set(case['checkpoint_steps'])
    return dict(initial=initial,checkpoints=checkpoints,index=dict(raw=record(path),header=header,header_sha256=header_sha256,header_bytes=header_bytes,
                                                                steps=index,first_step=1,last_step=case['steps']))

def selected(store,leaf,generation,expected,trusted):
    store=Path(store);path=store/leaf;assert path.parent==store and not path.is_symlink()
    gen,fields=parse(path.read_bytes());assert gen==generation and set(fields)=={'SelectionSchema','Branch','ArchiveSHA256','PhysicalPrefixSHA256','ContractSHA256'}
    assert fields['SelectionSchema']==packed(2,(),1) and fields['ContractSHA256']==packed(1,(32,),list(trusted))
    branch=fields['Branch']['raw'].decode('ascii');assert re.fullmatch(r'g_[A-Za-z0-9_-]{1,94}',branch)
    a,p=store/branch/'state.bin',store/branch/'prefix.phy';assert not a.is_symlink() and not p.is_symlink()
    archive,prefix=a.read_bytes(),p.read_bytes();g,candidate=parse(archive);assert g==generation and len(candidate)==149
    assert candidate['contract/Schema']==packed(2,(),2) and prefix[:48]==b'PERIPHY2'+struct.pack('<Q',2)+trusted
    assert fields['ArchiveSHA256']==packed(1,(32,),list(hashlib.sha256(archive).digest()))
    assert fields['PhysicalPrefixSHA256']==packed(1,(32,),list(hashlib.sha256(prefix).digest()))
    assert archive==expected['archive'],'complete v2 candidate bytes differ'
    assert prefix==expected['physical'],'complete v2 physical prefix bytes differ'
    return dict(generation=generation,branch=branch,selection=record(path),archive=record(a),physical=record(p),
                full_candidate_byte_equal=True,full_physical_prefix_byte_equal=True)

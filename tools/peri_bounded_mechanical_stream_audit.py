"""Independent conversion/byte oracle for preserved actual-PD observations.

No solver dispatch. This tests a NEW physical-stream format and reader, never
changes the verdict of the original physical study or independent reviews.
"""
import hashlib
import json
import math
import struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode, make_fields

RAW_SHA='7c9f8097bded6bcc079a07dcc18e582ca460d909edf4d12986321094fd26c4ec'
STATE={'GelPosition':'position','GelVelocity':'velocity','GelGreen':'green','GelViscous':'viscous',
       'GelPhysicalDissipationDensity':'dissipation','GelAcceptedCauchyStress':'accepted_stress',
       'GelAcceptedDensity':'density','LastNormalImpulse':'last_normal_impulse','AcceptedTime':'time',
       'AcceptedStep':'step','LastActualDt':'last_dt','LastKKT':'last_KKT','LastNewtonIterations':'last_newton_iterations'}
FRICTION={'LastTwoTangentImpulse':'last_tangent_impulse','TangentSlipSum':'relative_tangent_travel_sum',
          'TangentSlipCorrection':'relative_tangent_travel_correction','ContactCumulativeLedgerSum':'ledger_sum',
          'ContactCumulativeLedgerCorrection':'ledger_correction'}
ENERGY={'KineticOld':'kinetic_old','KineticNew':'kinetic_new','AlgorithmicLoss':'algorithmic_loss',
        'NormalContactLoss':'normal_contact_loss','FrictionLoss':'friction_loss','EnergyResidual':'energy_residual'}
EVALUATION={'MaterialStoredOld':'material_stored_old','MaterialStoredNew':'material_stored_new',
            'NonaffineStoredOld':'stabilization_stored_old','NonaffineStoredNew':'stabilization_stored_new',
            'MaxwellLoss':'physical_material_loss','MaterialTransposeWorkResidual':'transpose_work_residual',
            'NonaffineWorkResidual':'stabilization_work_residual'}
VECTORS={'TotalMomentumResidual':'total_momentum_residual','InternalForceSum':'internal_force_sum',
         'AngularDriftMinusExternal':'angular_drift_minus_external','OldInternalTorque':'internal_old_torque'}


def loads(line):
    # C++17-digit JSON spells an IEEE negative zero as integer token -0.
    return json.loads(line,parse_int=lambda text:-0.0 if text=='-0' else int(text))


def packed(tag,shape,value):
    def flatten(x):
        if isinstance(x,list):
            for y in x:yield from flatten(y)
        else:yield x
    x=list(flatten(value));assert len(x)==math.prod(shape)
    return dict(tag=tag,shape=list(shape),raw=struct.pack('<'+{1:'B',2:'Q',3:'d'}[tag]*len(x),*x))


def body_fields(template,state,history):
    fields=dict(template)
    assert history['relative_travel_includes_open_rows'] is True
    for mapping,source in [(STATE,state),(FRICTION,history)]:
        for label,key in mapping.items():
            name='accepted/'+label;old=fields[name]
            fields[name]=packed(old['tag'],old['shape'],source[key])
    return fields


def mechanical_fields(body,row=None):
    fields={'body/'+k:v for k,v in body.items()}
    def add(name,value,shape=(),tag=3):fields[name]=packed(tag,shape,value)
    def array(name,value):
        add(name+'/count',len(value),tag=2);add(name,value,(len(value),3))
    step=struct.unpack('<Q',body['accepted/AcceptedStep']['raw'])[0]
    add('RecordKind',int(step>0),tag=2);add('IntervalPhase',row['interval']['phase'] if row else 0,tag=2)
    add('OldTime',row['interval']['old_time'] if row else 0.)
    array('observation/ContactImpulse',row['contact_impulse'] if row else [])
    array('observation/NormalizedMomentumResidual',row['normalized_momentum_residual'] if row else [])
    plate_map={'Impulse':'particle_plate_impulse','NormalImpulse':'normal_plate_impulse',
               'TwoTangentImpulse':'tangent_plate_impulse','IntervalMeanReaction':'plate_reaction_interval',
               'NormalWork':'normal_plate_work','TangentWork':'tangent_plate_work',
               'WholeCubicMinimumGap':'dense_minimum_plate_gap','EndpointMinimumGap':'dense_endpoint_minimum_plate_gap'}
    for p,name in enumerate(['lower','upper']):
        for label,key in plate_map.items():
            shape=(3,) if label in ['Impulse','NormalImpulse','TwoTangentImpulse','IntervalMeanReaction'] else ()
            add('observation/'+name+'/'+label,row[key][p] if row else ([0.,0.,0.] if shape else 0.),shape)
    for label,key in ENERGY.items():add('observation/'+label,row[key] if row else 0.)
    for label,key in EVALUATION.items():add('observation/'+label,row['evaluation'][key] if row else 0.)
    for label,key in VECTORS.items():add('observation/'+label,row[key] if row else [0.,0.,0.],(3,))
    assert len(fields)==138
    return fields


def source_records(raw,template):
    with Path(raw).open('rb') as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==RAW_SHA
    assert struct.pack('<d',loads('-0'))==struct.pack('<Q',1<<63)
    count=0
    with Path(raw).open(encoding='utf8') as stream:
        for line in stream:
            row=loads(line)
            if 'reference' in row:
                yield mechanical_fields(body_fields(template,row['initial_state'],row['initial_friction_history']))
            elif 'state' in row:
                count+=1;assert row['committed'] and row['state']['step']==count
                yield mechanical_fields(body_fields(template,row['state'],row['friction_history']),row)
    assert count==40


def write(path,data):
    with Path(path).open('xb') as stream:stream.write(data)


def create(root,body_fixture,original_input,identities,raw):
    root=Path(root);body_data=Path(body_fixture).read_bytes();_,template=parse(body_data)
    candidate_fields,trusted=make_fields(body_data,original_input,identities)
    records=list(source_records(raw,template))
    assert len(records)==41
    chunks=[encode(fields,step+1) for step,fields in enumerate(records)]
    header=b'PERIPHY1'+struct.pack('<Q',1)+bytes.fromhex(trusted)
    def stream(chunks):return header+b''.join(struct.pack('<Q',len(x))+x for x in chunks)
    baseline=stream(chunks)
    def cursor(data,body_change=None):
        fields=dict(candidate_fields)
        fields['output/Bytes']=packed(2,(),len(data))
        fields['output/HeaderSHA256']=packed(1,(32,),list(hashlib.sha256(data[:48]).digest()))
        fields['output/PrefixSHA256']=packed(1,(32,),list(hashlib.sha256(data).digest()))
        if body_change:body_change(fields)
        return encode(fields,37)
    write(root/'physical.phy',baseline);write(root/'session_candidate.bin',cursor(baseline))
    tail=baseline+b'preserved_later_unselected_tail_NOT_A_RECORD\n'
    write(root/'later_tail.phy',tail)
    evidence=root/'evidence';evidence.mkdir();cases=[];evidence_inputs=[]
    def add(label,reason,data,change=None):
        c=label+'.candidate.bin';s=label+'.phy'
        write(evidence/c,cursor(data,change));write(evidence/s,data)
        cases.append((label,c,s,reason));evidence_inputs.extend([evidence/c,evidence/s])
    def mutate(fields,key,index,value):
        old=fields[key];fmt={1:'B',2:'Q',3:'d'}[old['tag']]
        values=list(struct.unpack('<'+fmt*(len(old['raw'])//struct.calcsize(fmt)),old['raw']))
        values[index]=value(values[index]);fields[key]=packed(old['tag'],old['shape'],values)
    def bad_record(label,reason,step,change):
        selected=dict(records[step]);change(selected);changed=list(chunks);changed[step]=encode(selected,step+1)
        add(label,reason,stream(changed))
    bad_header=bytearray(baseline);bad_header[0]^=1
    add('trusted_header','trusted physical header',bytes(bad_header))
    bad=bytearray(chunks[1]);bad[-40]^=1;bad_chunks=list(chunks);bad_chunks[1]=bytes(bad)
    add('record_checksum','SHA256 mismatch',stream(bad_chunks))
    bad_length=bytearray(baseline);bad_length[48:56]=struct.pack('<Q',16*1024*1024+1)
    add('record_length','record length',bytes(bad_length))
    bad_record('unknown_record_schema','exact schema names',1,lambda f:f.update(extra_unknown=packed(2,(),1)))
    bad_record('accepted_time_ulp','adjacent exact clock',1,lambda f:mutate(f,'body/accepted/AcceptedTime',0,lambda x:math.nextafter(x,1)))
    bad_record('actual_dt_ulp','adjacent exact clock',1,lambda f:mutate(f,'body/accepted/LastActualDt',0,lambda x:math.nextafter(x,0)))
    bad_record('KKT_cursor_ulp','original KKT cursor',1,lambda f:mutate(f,'body/accepted/LastKKT',0,lambda x:math.nextafter(x,1)))
    for label,key in [('reaction','observation/lower/IntervalMeanReaction'),('plate_impulse','observation/upper/Impulse'),
                      ('plate_work','observation/upper/NormalWork'),('friction_loss','observation/FrictionLoss'),
                      ('algorithmic_loss','observation/AlgorithmicLoss')]:
        bad_record(label,'physical observation byte mismatch',1,lambda f,k=key:mutate(f,k,0,lambda x:x+.01))
    bad_record('ledger_compensation','accepted history byte mismatch',1,
        lambda f:mutate(f,'body/accepted/ContactCumulativeLedgerCorrection',11,lambda x:x+.001))
    bad_record('travel_compensation','accepted history byte mismatch',1,
        lambda f:mutate(f,'body/accepted/TangentSlipCorrection',0,lambda x:x+.001))
    def dual_change(f):
        key='body/accepted/LastTwoTangentImpulse';values=struct.unpack('<232d',f[key]['raw'])
        row=max(range(116),key=lambda i:values[2*i]**2+values[2*i+1]**2)
        assert values[2*row]**2+values[2*row+1]**2>1e-20
        mutate(f,key,2*row,lambda x:x*.5);mutate(f,key,2*row+1,lambda x:x*.5)
    bad_record('dual_original_equation','original nonsymmetric equation/KKT or precision gate',1,dual_change)
    bad_record('material_position','accepted material cursor',1,
        lambda f:mutate(f,'body/accepted/GelPosition',0,lambda x:x+1e-5))
    bad_record('accepted_density','accepted material cursor',1,
        lambda f:mutate(f,'body/accepted/GelAcceptedDensity',0,lambda x:x*1.01))
    bad_record('static_friction','static body identity',0,
        lambda f:mutate(f,'body/contact/LowerMu',0,lambda x:x*2))
    bad_record('signed_MLS','stored signed surface weight mismatch',0,
        lambda f:mutate(f,'body/surface/SurfaceWeights',0,lambda x:x+.01))
    bad_record('initial_velocity','initial kinematics',0,
        lambda f:mutate(f,'body/accepted/GelVelocity',0,lambda x:x+.001))
    bad_record('initial_ledger','initial cumulative ledger',0,
        lambda f:mutate(f,'body/accepted/ContactCumulativeLedgerSum',6,lambda x:x+.001))
    add('final_archive_body','final archive/physical output body mismatch',baseline,
        lambda f:mutate(f,'body/accepted/ContactCumulativeLedgerCorrection',11,lambda x:x+.001))
    assert len(cases)==22
    (root/'cases.tsv').write_text(''.join('\t'.join(row)+'\n' for row in cases))
    return dict(trusted=trusted,records=records,bytes=len(baseline),cases=cases,
                input_files=[root/'physical.phy',root/'session_candidate.bin',root/'later_tail.phy',root/'cases.tsv']+evidence_inputs,
                original_raw_SHA256=RAW_SHA)


def parse_stream(data):
    assert data[:8]==b'PERIPHY1' and struct.unpack('<Q',data[8:16])[0]==1
    offset=48;result=[]
    while offset<len(data):
        length=struct.unpack('<Q',data[offset:offset+8])[0];offset+=8
        assert 72<=length<=16*1024*1024 and offset+length<=len(data)
        generation,fields=parse(data[offset:offset+length]);offset+=length
        assert generation==len(result)+1
        result.append(fields)
    assert offset==len(data)
    return result


def assess(root,prepared,native):
    root=Path(root);data=(root/'physical.phy').read_bytes()
    assert parse_stream(data)==prepared['records']
    assert (root/'staged/physical.phy').read_bytes()==data==(root/'staged/from_later_prefix.phy').read_bytes()
    assert (root/'later_tail.phy').read_bytes().startswith(data) and len((root/'later_tail.phy').read_bytes())>len(data)
    assert not (root/'staged/must_not_exist.phy').exists()
    expected={row[0] for row in prepared['cases']}|{'exact_snapshot_later_tail','changed_candidate_after_verification','semantic_read_bound'}
    labels=[label for label,reason in native['refusals']]
    assert len(labels)==len(set(labels))==25 and set(labels)==expected
    for label,candidate_name,stream_name,reason in prepared['cases']:
        _,candidate=parse((root/'evidence'/candidate_name).read_bytes());stream=(root/'evidence'/stream_name).read_bytes()
        assert candidate['output/PrefixSHA256']['raw']==hashlib.sha256(stream).digest()
        if label not in ['trusted_header','record_checksum','record_length']:
            parsed=parse_stream(stream)
            if label=='final_archive_body':assert parsed==prepared['records']
            elif label=='unknown_record_schema':assert set(parsed[1])-set(prepared['records'][1])=={'extra_unknown'}
            else:
                differences=[i for i in range(41) if parsed[i]!=prepared['records'][i]]
                assert len(differences)==1
                i=differences[0]
                assert set(parsed[i])==set(prepared['records'][i])
                assert len([k for k in parsed[i] if parsed[i][k]!=prepared['records'][i][k]])==1
    assert native['records']==41 and native['accepted_steps']==40 and native['bytes']==len(data) and native['positive_prefixes']==2
    assert native['nonlinear_solver_calls']==native['physical_advances']==native['live_commits']==0
    assert native['filesystem_transaction_qualified'] is False and native['full_restart_qualified'] is False
    return dict(status='converted_existing_stream_all_records_original_proposals_and_prefix_bytes_independently_matched_NOT_new_physical_restart',
        records=41,accepted_observations=40,fields_per_record=138,stream_bytes=len(data),positive_prefixes=2,
        correct_hash_semantic_refusals=20,integrity_refusals=2,external_bound_or_staging_refusals=3,
        preserved_later_tail=True,staged_full_prefix_bytes_identical=True,original_raw_SHA256=RAW_SHA,
        negative_zero_token_preserved=True,nonlinear_solver_calls=0,physical_advances=0,live_commits=0,
        filesystem_transaction_qualified=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,
        limits='Converted old raw under new codec; no new native physical input, no old physical verdict change, no transaction publish or full restart trajectory qualification.')

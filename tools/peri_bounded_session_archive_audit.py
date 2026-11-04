"""Independent byte and accepted-clock oracle, NOT physical-prefix verification."""
import hashlib
import math
import struct
from pathlib import Path
from peri_canonical_archive_audit import parse


def encode(fields, generation):
    pieces=[b'PERIAR01',struct.pack('<QQQ',1,generation,len(fields))]
    for name,field in sorted(fields.items()):
        key=name.encode('ascii');shape=field['shape'];raw=field['raw']
        pieces += [struct.pack('<Q',len(key)),key,bytes([field['tag'],len(shape)]),
                   struct.pack('<'+'Q'*len(shape),*shape),struct.pack('<Q',len(raw)),raw]
    prefix=b''.join(pieces)
    return prefix+hashlib.sha256(prefix).digest()+b'PERIEND1'


def make_fields(body_data, original_input, identities):
    generation,body=parse(body_data)
    assert generation==37 and len(body)==98
    values=Path(original_input).read_text().split()
    assert len(values)==17
    kind,rotation,azimuths=map(int,values[:3]);dt,stroke,t0,t1=map(float,values[3:7])
    assert kind==2 and rotation==0 and azimuths==8
    fields={'body/'+k:dict(v) for k,v in body.items()}
    def add(name,value,tag=2,shape=()):
        assert name not in fields
        x=value if isinstance(value,list) else [value]
        assert len(x)==math.prod(shape)
        fields[name]=dict(tag=tag,shape=list(shape),raw=struct.pack('<'+{1:'B',2:'Q',3:'d'}[tag]*len(x),*x))
    def digest(name,value):
        raw=bytes.fromhex(value)
        assert len(raw)==32 and any(raw)
        add(name+'/count',32)
        add(name,list(raw),1,(32,))
    for name in ['Schema','OriginalEquationVersion','AffineFMaterialPathVersion','SignedTransferVersion',
                 'CompensationVersion','FloatingPointProfile','PhysicalStreamSchema','SceneMode']:
        add('contract/'+name,1)
    for name,value in [('NormalBaseSteps',16),('Rotation',rotation),('Azimuths',azimuths),('ProtocolKind',kind),
                       ('MaxNewton',48),('MaxLineSearch',32)]:
        add('contract/'+name,value)
    for name,value in [('CapWidth',float(values[16])),('PenetrationStrokeGate',.01),('NominalDt',dt),
                       ('NormalStroke',stroke),('SolverVelocityScale',1.),('SolverRelativeTolerance',1e-12)]:
        add('contract/'+name,value,3)
    add('contract/TwoTangentStrokes',[t0,t1],3,(2,))
    static={k:v for k,v in body.items() if not k.startswith('accepted/')}
    fingerprints=dict(identities,Input=hashlib.sha256(Path(original_input).read_bytes()).hexdigest(),
                      StaticBody=hashlib.sha256(encode(static,1)).hexdigest())
    assert set(fingerprints)=={'SourceManifest','ELF','BuildRuntimeManifest','Input','StaticBody'}
    for key,value in fingerprints.items():digest('contract/'+key+'SHA256',value)
    add('accepted_algorithm/NoInflightTrial',1)
    add('accepted_algorithm/NextPhase',3)
    add('output/Schema',1);add('output/Records',41);add('output/Bytes',len(body_data));add('output/LastStep',40)
    add('output/LastTime',.4,3)
    # Actual engineering files supply these hashes. They are deliberately NOT
    # a physical stream: the candidate remains unusable for a live restore.
    digest('output/HeaderSHA256',hashlib.sha256(Path(original_input).read_bytes()).hexdigest())
    digest('output/PrefixSHA256',hashlib.sha256(body_data).hexdigest())
    contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    assert len(fields)==140
    return fields,hashlib.sha256(encode(contract,1)).hexdigest()


CHANGES={
    'source_identity':'contract/SourceManifestSHA256','ELF_identity':'contract/ELFSHA256',
    'build_runtime_identity':'contract/BuildRuntimeManifestSHA256','input_identity':'contract/InputSHA256',
    'protocol_tangent':'contract/TwoTangentStrokes','cap_width':'contract/CapWidth',
    'solver_gate':'contract/SolverRelativeTolerance','geometry_gate':'contract/PenetrationStrokeGate',
    'algorithm_version':'contract/OriginalEquationVersion','static_friction':'body/contact/LowerMu',
    'static_surface_area':'body/surface/SurfaceAreas','accepted_time_ulp':'body/accepted/AcceptedTime',
    'actual_dt_ulp':'body/accepted/LastActualDt','accepted_boundary':'accepted_algorithm/NoInflightTrial',
    'next_phase':'accepted_algorithm/NextPhase','output_records':'output/Records',
    'output_step':'output/LastStep','output_time_ulp':'output/LastTime','output_bytes':'output/Bytes'}
EXTERNAL={'wrong_external_trust','wrong_explicit_generation','different_rounding_mode','flush_subnormal_mode'}


def expected_clocks(dt):
    clocks=[]
    bits=lambda x:struct.unpack('<Q',struct.pack('<d',x))[0]
    for kind in range(3):
        time=0.
        for step in range(16 if kind==0 else 40):
            end=(16 if step<16 else (24 if step<24 else 40))*dt
            actual=end-time if step in [15,23,39] else dt
            time=time+actual
            clocks.append([bits(time),bits(actual)])
        assert time==(16 if kind==0 else 40)*dt
    return clocks


def assess(fields,fixture,evidence,native):
    data=Path(fixture).read_bytes();generation,decoded=parse(data)
    assert generation==37 and decoded==fields and data==encode(fields,37)
    assert Path(evidence,'roundtrip.bin').read_bytes()==data
    assert native['status']=='detached_session_contract_component_checks_passed_NOT_output_or_restart_qualification'
    assert native['physical_advances']==native['live_commits']==0
    assert native['output_prefix_verified'] is False and native['full_restart_qualified'] is False
    labels=[label for label,reason in native['refusals']]
    assert len(labels)==len(set(labels))==24
    assert set(labels)==set(CHANGES)|EXTERNAL|{'unknown_schema'}
    for label in set(CHANGES)|{'unknown_schema'}:
        cg,changed=parse(Path(evidence,label+'.bin').read_bytes())
        assert cg==37
        if label=='unknown_schema':
            assert set(changed)==set(fields)|{'unregistered'}
            assert all(changed[k]==v for k,v in fields.items())
        else:
            assert set(changed)==set(fields)
            assert {k for k in fields if fields[k]!=changed[k]}=={CHANGES[label]}
            key=CHANGES[label]
            assert fields[key]['shape']==changed[key]['shape'] and fields[key]['tag']==changed[key]['tag']
    dt=struct.unpack('<d',fields['contract/NominalDt']['raw'])[0]
    assert native['clock_checks']==96 and native['clock_bits']==expected_clocks(dt)
    return dict(status='detached_contract_bytes_external_identity_and_exact_clocks_independently_matched_NOT_full_CP',
        fields=140,archive_bytes=len(data),valid_hash_semantic_refusals=20,external_generation_or_FP_refusals=4,
        exact_protocol_time_dt_pairs=96,unchanged_baseline_bytes=True,physical_advances=0,live_commits=0,
        output_prefix_verified=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,
        limits='Engineering cursor hashes intentionally refer to real non-stream fixtures. This candidate has no live commit API; whole physical output semantics, transaction and restart remain unqualified.')

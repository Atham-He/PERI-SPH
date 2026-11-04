"""Independent struct/hashlib body-codec oracle using preserved nonzero raw.

No solver, material update, reference constructor or live restore is called.
This component audit cannot qualify a complete session restart.
"""
import hashlib
import json
import math
import struct
from pathlib import Path

from peri_canonical_archive_audit import parse

RAW_SHA = '7c9f8097bded6bcc079a07dcc18e582ca460d909edf4d12986321094fd26c4ec'
GENERATION = 37


def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(path), size=path.stat().st_size, sha256=digest)


def source_fields(path):
    assert record(path)['sha256'] == RAW_SHA
    header = final = None
    steps = 0
    with Path(path).open(encoding='utf8') as stream:
        for line in stream:
            row = json.loads(line)
            if 'reference' in row:
                assert header is None
                header = row
            if 'state' in row:
                steps += 1
                assert row['committed'] and row['local_step_gates_passed']
                assert row['state']['step'] == steps
                final = row
    assert steps == 40 and final['state']['time'] == .4
    friction = final['friction_history']
    assert friction['relative_travel_includes_open_rows'] is True
    assert friction['ledger_order'] == 'lower_Ixyz,upper_Ixyz,lower_work,upper_work,Maxwell,normal,friction,algorithmic'
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


def build(fields):
    pieces = [b'PERIAR01', struct.pack('<QQQ', 1, GENERATION, len(fields))]
    for name, field in sorted(fields.items()):
        encoded = name.encode('ascii')
        shape = field['shape']
        pieces += [struct.pack('<Q',len(encoded)), encoded, bytes([field['tag'],len(shape)]),
                   struct.pack('<'+'Q'*len(shape),*shape), struct.pack('<Q',len(field['raw'])), field['raw']]
    prefix = b''.join(pieces)
    return prefix+hashlib.sha256(prefix).digest()+b'PERIEND1'


CHANGES = {
    'active_mask': 'model/PDReferenceActive', 'kernel_weight': 'model/PDReferenceWeight',
    'shape_inverse': 'model/PDShapeInverse', 'shape_gradient_null_affine': 'model/PDShapeGradient',
    'nonaffine_cache': 'model/PDNonaffineUnitCache', 'signed_MLS': 'surface/SurfaceWeights',
    'physical_mass': 'model/PhysicalNodeMass', 'material_green': 'accepted/GelGreen',
    'material_viscous': 'accepted/GelViscous', 'material_dissipation': 'accepted/GelPhysicalDissipationDensity',
    'accepted_stress': 'accepted/GelAcceptedCauchyStress', 'accepted_density': 'accepted/GelAcceptedDensity',
    'normal_dual': 'accepted/LastNormalImpulse', 'tangent_dual': 'accepted/LastTwoTangentImpulse',
    'Maxwell_compensation': 'accepted/ContactCumulativeLedgerCorrection',
    'contact_order': 'contact/PlateIDs', 'contact_frame': 'contact/Tangent0'}


def assess(raw, fixture, output, native):
    expected = source_fields(raw)
    data = Path(fixture).read_bytes()
    generation, fields = parse(data)
    assert generation == GENERATION and fields == expected and data == build(expected)
    assert Path(output, 'roundtrip.bin').read_bytes() == data
    assert native['status'] == 'nonzero_body_component_codec_checks_passed_not_full_session_CP'
    assert native['physical_advances'] == 0 and native['reference_rebuilt'] is False
    assert native['full_restart_qualified'] is False and native['bytes'] == len(data)
    labels = [item[0] for item in native['refusals']]
    assert len(labels) == len(set(labels)) == 18 and set(labels) == set(CHANGES)|{'unknown_schema'}
    evidence = []
    for label in labels:
        path = Path(output, label+'.bin')
        cg, changed = parse(path.read_bytes())  # Valid independent SHA/envelope despite corrupt semantics.
        assert cg == GENERATION
        if label == 'unknown_schema':
            assert set(changed) == set(expected)|{'unknown_field'}
            assert all(changed[k] == v for k,v in expected.items())
        else:
            assert set(changed) == set(expected)
            differences = {k for k in expected if changed[k] != expected[k]}
            assert differences == {CHANGES[label]}
            assert changed[CHANGES[label]]['tag'] == expected[CHANGES[label]]['tag']
            assert changed[CHANGES[label]]['shape'] == expected[CHANGES[label]]['shape']
        evidence.append(record(path))
    assert all(native[k] > 0 for k in ['nonzero_viscous_nodes','nonzero_tangent_rows','nonzero_travel_corrections'])
    return dict(status='existing_nonzero_body_all_fields_bits_and_valid_hash_semantic_refusals_independently_matched',
        field_count=len(fields), payload_bytes=sum(len(x['raw']) for x in fields.values()), archive_bytes=len(data),
        array_count_companions=sum(k.endswith('/count') for k in fields), valid_hash_semantic_refusals=18,
        corruption_evidence=evidence, raw=record(raw), fixture=record(fixture),
        nonzero_viscous_nodes=native['nonzero_viscous_nodes'],nonzero_tangent_rows=native['nonzero_tangent_rows'],
        nonzero_travel_corrections=native['nonzero_travel_corrections'],
        physical_advances=0, reference_rebuilt=False, full_restart_qualified=False, S3_complete=False,
        limits='Stored body only; protocol/algorithm/output transaction and physical restart trajectories remain unqualified.')

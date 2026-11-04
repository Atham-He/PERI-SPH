"""Zero-advance driver build/initialization only; no physical campaign release."""
import argparse, hashlib, json, shutil, struct, subprocess, traceback
from pathlib import Path
from peri_archive_component_build import prepare, compile_and_freeze
from peri_run_body_archive_engineering import run, reserves, save, stamp, identity
from peri_bounded_body_archive_audit import record
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import mechanical_fields, packed


def assess(root, baseline, build, native):
    candidate=(root/'initial_candidate.bin').read_bytes(); generation,fields=parse(candidate)
    assert generation==1 and len(fields)==140
    data=(root/'initial.phy').read_bytes(); trusted=native['contract_sha256']
    assert data[:48]==b'PERIPHY1'+struct.pack('<Q',1)+bytes.fromhex(trusted)
    size=struct.unpack('<Q',data[48:56])[0]; assert len(data)==56+size
    _,initial=parse(data[56:]); assert len(initial)==138
    with baseline.open('rb') as stream:
        assert stream.read(8)==b'PERIPHY1';stream.read(40)
        length=struct.unpack('<Q',stream.read(8))[0];_,old=parse(stream.read(length))
    # Independent t0 regression only. Same geometry/material/mass/reference;
    # new input has different prescribed tangent amplitudes and friction.
    body={k.removeprefix('body/'):v for k,v in old.items() if k.startswith('body/')}
    for side in ['LowerMu','UpperMu']:body['contact/'+side]=packed(3,(),.22)
    expected=mechanical_fields(body)
    assert initial==expected and data[56:]==encode(expected,1)
    assert {k.removeprefix('body/'):v for k,v in fields.items() if k.startswith('body/')}==body
    contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    static={k:v for k,v in body.items() if not k.startswith('accepted/')}
    fingerprints=dict(SourceManifest=build['source_manifest']['sha256'],ELF=build['ELF']['sha256'],
        BuildRuntimeManifest=build['build_runtime_manifest']['sha256'],Input=record(root/'original.input')['sha256'],
        StaticBody=hashlib.sha256(encode(static,1)).hexdigest())
    expected_contract={}
    def add(name,value,tag=2,shape=()):expected_contract[name]=packed(tag,shape,value)
    for name in ['Schema','OriginalEquationVersion','AffineFMaterialPathVersion','SignedTransferVersion',
                 'CompensationVersion','FloatingPointProfile','PhysicalStreamSchema','SceneMode']:add(name,1)
    for name,value in [('NormalBaseSteps',16),('Rotation',0),('Azimuths',8),('ProtocolKind',2),('MaxNewton',48),('MaxLineSearch',32)]:add(name,value)
    for name,value in [('CapWidth',.015625),('PenetrationStrokeGate',.01),('NominalDt',.01),('NormalStroke',.01),('SolverVelocityScale',1.),('SolverRelativeTolerance',1e-12)]:add(name,value,3)
    add('TwoTangentStrokes',[.0006,-.0004],3,(2,))
    for key,value in fingerprints.items():add(key+'SHA256/count',32);add(key+'SHA256',list(bytes.fromhex(value)),1,(32,))
    assert contract==expected_contract
    encoded=encode(contract,1);assert (root/'contract.bin').read_bytes()==encoded
    assert hashlib.sha256(encoded).hexdigest()==trusted
    def field(name,value,tag=2,shape=()):assert fields[name]==packed(tag,shape,value)
    for name,value in [('accepted_algorithm/NoInflightTrial',1),('accepted_algorithm/NextPhase',0),
                       ('output/Schema',1),('output/Records',1),('output/Bytes',len(data)),('output/LastStep',0)]:field(name,value)
    field('output/LastTime',0.,3)
    for label,prefix in [('HeaderSHA256',data[:48]),('PrefixSHA256',data)]:
        field('output/'+label+'/count',32);field('output/'+label,list(hashlib.sha256(prefix).digest()),1,(32,))
    assert candidate==encode(fields,1)
    assert native['new_physical_steps']==0 and native['full_restart_qualified'] is False
    return dict(status='actual_driver_zero_advance_initialization_bytes_independently_matched_NOT_physical_restart',
        candidate_fields=140,physical_record_fields=138,candidate_bytes=len(candidate),physical_bytes=len(data),
        initial_body_matches_preserved_t0_except_preregistered_mu=True,new_physical_steps=0,contract_sha256=trusted,
        full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);args=parser.parse_args()
    root=Path(args.root).resolve();repo=Path(__file__).resolve().parents[1]
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    reserves(root);assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)=='';root.mkdir(mode=0o700)
    try:
        stage='tests/3d_examples/extra_3d_PD_framework_tests/'
        files=[stage+x for x in ['bounded_gel_checkpoint_probe.cpp','bounded_gel_friction_observer.h','coulomb_contact_probe_io.h','bounded_gel_probe_io.h','bounded_gel_friction_probe.cpp']]
        files+=['tools/peri_actual_checkpoint_preflight.py','tools/peri_archive_component_build.py','tools/peri_run_body_archive_engineering.py',
            'tools/peri_bounded_body_archive_audit.py','tools/peri_canonical_archive_audit.py','tools/peri_bounded_session_archive_audit.py',
            'tools/peri_bounded_mechanical_stream_audit.py','docs/peri/ACTUAL_CHECKPOINT_DRIVER_20261105.md']
        prepared=prepare(root,repo,files,stage+'bounded_gel_checkpoint_probe.cpp','native_checkpoint')
        baseline=root.parent/'engineering-bounded-mechanical-stream/physical.phy'
        protected=json.loads((root.parent/'cross-runtime-v19-raw/plan.json').read_text())['artifacts']
        protected=[x for x in protected if x['path'].startswith('/home/')]
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        anchors=protected+[record(production),record(baseline),record(root.parent/'cross-runtime-v19-raw/plan.json')]
        assert all(record(x['path'])==x for x in anchors)
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        save(root/'plan.json',dict(prepared_at=stamp(),source_commit=source,source_snapshots=prepared['source_snapshots'],anchors=anchors,
            compile_command=prepared['compile_command'],link_command=prepared['link_command'],build_cwd=str(prepared['builddir']),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),native_budget=dict(CPU_seconds=30,wall_seconds=45,AS_GiB=1,threads=1),
            declared_modes=['--capabilities','--initialize','--diagnostic-fixture'],new_physical_steps=0,physical_dispatch_allowed=False,
            gates=['All original/actual includes/source/runtime/production/reference hashes unchanged',
                'Compile actual original nonsymmetric Newton entry, but never invoke physical modes',
                'Explicit initialized140field candidate/138field record and old t0 body independently byte-matched except new mu',
                'Actual self/source/build-runtime/newinput/external contract/header/fullprefix SHA independently reconstructed',
                'No claim of native trajectory, noadvance/split/tail physical restart, S3/S4 or calibration qualification'],
            failure_policy='Stop on unexpected compile/execution/integrity/resource/accuracy failure; preserve raw/returncode; do not redispatch this state.',
            total_bound_GiB=1,file_bound_MiB=256,host_reserve_GiB=120,output_reserve_GiB=150))
        save(root/'state.json',dict(status='zero_advance_driver_preflight_attempted',identity=identity(__import__('os').getpid()),started_at=stamp()))
        build=compile_and_freeze(prepared)
        assert 'SPH::peridynamics::solveCoulombContact(' in subprocess.check_output(['nm','-C',build['ELF']['path']],text=True)
        with (root/'original.input').open('x') as stream:stream.write('2 0 8 0.01 0.01 0.0006 -0.0004 0.22 0.22 1.0 5.0 1.0 0.5 0.3 0.1 1.25 0.015625\n')
        inputs=[record(root/'original.input')]; frozen=dict(plan=record(root/'plan.json'),build=build,inputs=inputs,anchors=anchors)
        save(root/'frozen_execution.json',frozen)
        phases=[]
        for label,arguments in [('capabilities',['--capabilities']),('initialize',['--initialize',str(root)]),('diagnostic_fixture',['--diagnostic-fixture'])]:
            native=run(root,label,[build['ELF']['path']]+arguments,root,30,45,1)
            rows=[json.loads(x) for x in (root/(label+'.stdout')).read_text().splitlines()]
            assert native['peak_parent_threads']<=1
            for x in anchors+build['support']+build['source_snapshots']+inputs+[build['ELF'],build['source_manifest'],build['build_runtime_manifest']]:assert record(x['path'])==x
            phases.append(native)
            if label=='capabilities':assert rows[-1]['actual_PD_checkpoint_driver'] and not rows[-1]['full_restart_qualified']
            elif label=='initialize':assessment=assess(root,baseline,build,rows[-1])
            else:
                assert rows[-1]==dict(type='diagnostic_fixture',controlled_formatter_only=True,new_physical_steps=0)
                diagnostic=json.loads((root/(label+'.stderr')).read_text())
                assert diagnostic==dict(type='checkpoint_driver_failure',mechanical_durable_step=7,mechanical_durable_time=.07,
                    mechanical_durable_bytes=123456,last_full_debug_step=6,last_full_debug_time=.06,actual_solver_calls=7,
                    publication_error_number=5,publication_stage=3,message='controlled formatter fixture: quote" backslash\\ newline\n',append_witness=None)
                assessment['allocation_free_failure_formatter_fixture_independently_matched']=True
        assessment.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),frozen_execution=record(root/'frozen_execution.json'),
            actual_compiler_dependency_count=build['actual_compiler_dependency_count'],phases=phases)
        save(root/'independent_preflight.json',assessment);save(root/'terminal.json',dict(status=assessment['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_preflight.json')))
        print(json.dumps(assessment))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),new_physical_steps=0,full_restart_qualified=False));raise


if __name__=='__main__':main()

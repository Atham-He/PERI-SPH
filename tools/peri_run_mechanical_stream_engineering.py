"""Single converted-data prefix/readonly proposal test, no physical solver run."""
import argparse,json,os,subprocess,traceback
from pathlib import Path
from peri_run_body_archive_engineering import run,reserves,save,identity,stamp
from peri_bounded_body_archive_audit import record
from peri_archive_component_build import prepare,compile_and_freeze
from peri_bounded_mechanical_stream_audit import create,assess


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);args=parser.parse_args()
    root=Path(args.root).resolve();repo=Path(__file__).resolve().parents[1]
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    reserves(root);assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    root.mkdir(mode=0o700)
    try:
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        body_root=root.parent/'engineering-nonzero-body-archive';body=body_root/'nonzero_body.bin'
        body_receipt=json.loads((body_root/'engineering_receipt.json').read_text());assert record(body)==body_receipt['fixture']
        original_input=root.parent/'study-v19-early-friction-cycle/actual_PD_early_biaxial_slip_cycle.input'
        raw=root.parent/'study-v19-early-friction-cycle/actual_PD_early_biaxial_slip_cycle.jsonl'
        active=root.parent/'cross-runtime-v19-raw/plan.json';active_plan=json.loads(active.read_text())
        protected=[x for x in active_plan['artifacts'] if x['path'].startswith('/home/')]
        assert all(record(x['path'])==x for x in protected)
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        anchors=[record(body),record(body_root/'engineering_receipt.json'),record(original_input),record(raw),record(active),record(production)]+protected
        files=['tools/bounded_mechanical_stream_checks.cpp','src/shared/peridynamics/pd_bounded_mechanical_stream.h',
            'src/shared/peridynamics/pd_bounded_gel_friction.h','src/shared/peridynamics/pd_bounded_friction_protocol.h',
            'src/shared/peridynamics/pd_output_prefix.h','tools/peri_bounded_mechanical_stream_audit.py',
            'tools/peri_run_mechanical_stream_engineering.py','tools/peri_archive_component_build.py',
            'tools/peri_run_body_archive_engineering.py','tools/peri_bounded_body_archive_audit.py',
            'tools/peri_bounded_session_archive_audit.py','tools/peri_canonical_archive_audit.py',
            'docs/peri/BOUNDED_MECHANICAL_STREAM_20261105.md']
        prepared=prepare(root,repo,files,'tools/bounded_mechanical_stream_checks.cpp','mechanical')
        plan=dict(prepared_at=stamp(),source_commit=source,anchors=anchors,source_snapshots=prepared['source_snapshots'],
            compile_command=prepared['compile_command'],link_command=prepared['link_command'],build_cwd=str(prepared['builddir']),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),
            native_budget=dict(CPU_seconds=900,wall_seconds=1200,AS_GiB=1,threads=1),
            native_cost_basis='Preserved original40step native205.775s includes all Newton solves. Three full readonly prefix walks plus22 early-cause controls, each scratch accepted proposal only, no Newton. Bound900CPU/1200wall prospectively exceeds3x full-native cost plus early-control allowance; no old failed budget extended.',
            fixture_recipe='Convert existing raw header+40accepted records to new138field canonical physical stream, preserve negative-zero JSON tokens. Actual new checker ELF/manifests are engineering identities, not original raw producer identities. Original physical verdict unchanged.',
            file_bound_MiB=256,total_bound_GiB=1,host_reserve_GiB=120,output_reserve_GiB=150,
            gates=['41 records/all138fields independently encoded and parsed; all accepted proposal body and physical-observation bytes exactly match original saved values',
                'Original final nonsymmetric R/KKT/precision, trueM ledgers, material history, geometry and cumulative sums/corrections re-evaluated in scratch without nonlinear solve',
                'Full selected prefix and exclusive staged bytes identical, later unselected tail preserved',
                '20 correct-hash semantic+2 integrity+3 external-bound/staging refusals with expected causes; no live commit',
                'All old source/raw/production and active reference artifacts unchanged'],
            failure_policy='Preserve exact error/cursor/returncode/raw and all frozen evidence; stop this state. No old physical rerun/review, no nonlinear solver dispatch/reset/fallback/gate relaxation.',
            nonlinear_solver_calls=0,physical_advances=0,live_commits=0,filesystem_transaction_qualified=False,
            full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan);save(root/'state.json',dict(status='single_converted_stream_engineering_attempted',identity=identity(os.getpid()),started_at=stamp()))
        build=compile_and_freeze(prepared)
        symbols=subprocess.check_output(['nm','-C',build['ELF']['path']],text=True)
        assert 'SPH::peridynamics::solveCoulombContact(' not in symbols,'unexpected nonlinear solver symbol in readonly checker'
        save(root/'no_nonlinear_solver_symbol.json',dict(proven_scope='source call graph plus actual linked symbol inspection; no nonlinear solver entrypoint in this reader',
            absent_symbol='SPH::peridynamics::solveCoulombContact(',ELF=build['ELF']))
        identities=dict(SourceManifest=build['source_manifest']['sha256'],BuildRuntimeManifest=build['build_runtime_manifest']['sha256'],ELF=build['ELF']['sha256'])
        fixture=create(root,body,original_input,identities,raw)
        save(root/'external_trust.json',dict(contract_SHA256=fixture['trusted'],identities=identities,raw=record(raw),converted_existing_data_only=True))
        (root/'staged').mkdir()
        input_records=[record(p) for p in fixture['input_files']]
        execution=dict(recorded_at=stamp(),plan=record(root/'plan.json'),build=build,anchors=anchors,
            inputs=input_records,external_trust=record(root/'external_trust.json'),no_nonlinear_solver_symbol=record(root/'no_nonlinear_solver_symbol.json'))
        for item in anchors+build['support']+build['source_snapshots']:assert record(item['path'])==item
        save(root/'frozen_execution.json',execution)
        command=[build['ELF']['path'],str(root),fixture['trusted']]
        native=run(root,'native',command,root,900,1200,1)
        for item in anchors+build['support']+build['source_snapshots']+input_records+[execution[k] for k in ['plan','external_trust','no_nonlinear_solver_symbol']]+[build['ELF'],build['source_manifest'],build['build_runtime_manifest']]:
            assert record(item['path'])==item
        report=json.loads((root/'native.stdout').read_text())
        if native['identity']:assert report['native_identity']==native['identity']
        else:
            assert report['native_identity']['argv']==command and report['native_identity']['ppid']==os.getpid()
            assert report['native_identity']['boot_id']==identity(os.getpid())['boot_id']
        result=assess(root,fixture,report)
        result.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),frozen_execution=record(root/'frozen_execution.json'),
            native=native,staged=[record(p) for p in (root/'staged').iterdir()],actual_compiler_dependency_count=build['actual_compiler_dependency_count'],
            support_original_and_copy_records=build['support_original_and_copy_records'],no_nonlinear_solver_symbol=record(root/'no_nonlinear_solver_symbol.json'))
        save(root/'engineering_receipt.json',result)
        save(root/'terminal.json',dict(status=result['status'],finished_at=stamp(),return_code=0,receipt=record(root/'engineering_receipt.json'),
            physical_advances=0,live_commits=0,filesystem_transaction_qualified=False,full_restart_qualified=False))
        print(json.dumps(dict(status=result['status'],records=41,refusals=25)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),
            nonlinear_solver_calls=0,physical_advances=0,live_commits=0,full_restart_qualified=False,S3_complete=False))
        raise


if __name__=='__main__':main()

"""One new prospective two-process generation-store test, no physical solver."""
import argparse,json,os,shutil,subprocess,traceback
from pathlib import Path
from peri_run_body_archive_engineering import run,reserves,save,identity,stamp
from peri_bounded_body_archive_audit import record
from peri_archive_component_build import prepare,compile_and_freeze
from peri_bounded_generation_store_audit import assess_publish,assess_restore


def verify_native(report,execution,command):
    if execution['identity']:assert report['native_identity']==execution['identity']
    else:
        assert report['native_identity']['argv']==command and report['native_identity']['ppid']==os.getpid()
        assert report['native_identity']['boot_id']==identity(os.getpid())['boot_id']
    assert execution['peak_parent_threads']<=1


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);args=parser.parse_args()
    root=Path(args.root).resolve();repo=Path(__file__).resolve().parents[1]
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    reserves(root);assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    root.mkdir(mode=0o700)
    try:
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        old=root.parent/'engineering-bounded-mechanical-stream';old_receipt=json.loads((old/'engineering_receipt.json').read_text())
        assert old_receipt['status']=='converted_existing_stream_all_records_original_proposals_and_prefix_bytes_independently_matched_NOT_new_physical_restart'
        assert (old/'physical.phy').stat().st_size==12394445
        inputs=[]
        for name in ['physical.phy','session_candidate.bin']:
            shutil.copy2(old/name,root/name);inputs += [record(old/name),record(root/name)]
        trusted=json.loads((old/'external_trust.json').read_text())['contract_SHA256']
        active=root.parent/'cross-runtime-v19-raw/plan.json';active_plan=json.loads(active.read_text())
        protected=[x for x in active_plan['artifacts'] if x['path'].startswith('/home/')]
        assert all(record(x['path'])==x for x in protected)
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        anchors=inputs+[record(old/'engineering_receipt.json'),record(old/'external_trust.json'),record(production),record(active)]+protected
        files=['tools/bounded_generation_store_checks.cpp','src/shared/peridynamics/pd_bounded_generation_store.h',
            'src/shared/peridynamics/pd_output_prefix.h','src/shared/peridynamics/pd_bounded_mechanical_stream.h',
            'src/shared/peridynamics/pd_canonical_archive.h','tools/peri_bounded_generation_store_audit.py',
            'tools/peri_run_generation_store_engineering.py','tools/peri_archive_component_build.py',
            'tools/peri_run_body_archive_engineering.py','tools/peri_bounded_body_archive_audit.py',
            'tools/peri_bounded_session_archive_audit.py','tools/peri_canonical_archive_audit.py',
            'docs/peri/BOUNDED_GENERATION_STORE_20261105.md']
        prepared=prepare(root,repo,files,'tools/bounded_generation_store_checks.cpp','generation')
        plan=dict(prepared_at=stamp(),source_commit=source,anchors=anchors,source_snapshots=prepared['source_snapshots'],
            compile_command=prepared['compile_command'],link_command=prepared['link_command'],build_cwd=str(prepared['builddir']),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),
            native_budget_per_phase=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=1,threads=1),
            cost_basis='Previous converted-prefix checker3full walks+25controls took22.487s. This is two declared processes with2full walks in publish and1in explicit restore,copy/fsync/rename/lock controls. Bound180CPU/240wall per phase prospectively includes ample filesystem margin, not an old failure retry.',
            phases=['--publish: actual F publication,19expected refusal controls,2durable publishes and1explicit injected post-rename uncertainty; preserve source and synthetic abandoned-working tail',
                    '--restore: only after phase1 independent byte oracle passes and outputs freeze; separate new process explicitly selects older generation37 while current HEAD38, publishes new branch with identical archive/prefix bytes'],
            external_trusted_contract=trusted,file_bound_MiB=256,total_bound_GiB=1,host_reserve_GiB=120,output_reserve_GiB=150,
            gates=['Current F file+directory fsync/exclusive create/flock/atomic rename succeed on actual volume; unsupported call fails this state without durability fallback',
                'All selected generations decoded from explicit external generation,not inferred from HEAD or later tail',
                'All candidate/prefix semantics and proof binding before private staging; all allocation/open/hash/fsync before HEAD rename',
                'Pre-publication negatives leave live owned bytes/owner identity/HEAD unchanged; input-invalid cases create no private branch',
                'Controlled post-rename failure returns published-but-uncertain,new owner known,previous pointer/branches retained,blocks further same-process publish',
                'Independent struct/hashlib bytes verify every archive/prefix/pointer,retained newer branch and synthetic tail; separate-process restore fully identical same-generation archive',
                'No actual power-loss claim,new physical solver,accepted-step writer or complete physical restart qualification'],
            failure_policy='Any unexpected execution/integrity/resource/accuracy/FS capability failure preserves precise exception/returncode/artifacts and stops later phase. No retry,old solver run/fallback/reset/gate relaxation.',
            physical_advances=0,running_solver_ELF_identity_bound=False,accepted_step_writer_qualified=False,
            full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan);save(root/'state.json',dict(status='two_process_generation_engineering_attempted',identity=identity(os.getpid()),started_at=stamp()))
        build=compile_and_freeze(prepared)
        assert 'SPH::peridynamics::solveCoulombContact(' not in subprocess.check_output(['nm','-C',build['ELF']['path']],text=True)
        (root/'store').mkdir()
        execution=dict(recorded_at=stamp(),plan=record(root/'plan.json'),build=build,anchors=anchors,inputs=inputs,external_trusted_contract=trusted)
        for item in anchors+build['support']+build['source_snapshots']:assert record(item['path'])==item
        save(root/'frozen_execution.json',execution)
        command=[build['ELF']['path'],'--publish',str(root),trusted]
        first=run(root,'publish',command,root,180,240,1);first_report=json.loads((root/'publish.stdout').read_text())
        verify_native(first_report,first,command);publish=assess_publish(root,first_report)
        save(root/'publish_independent_assessment.json',publish)
        # Freeze phase1 evidence BEFORE the separate restore process; HEAD is
        # the only permitted overwrite,all other phase1 evidence is retained.
        preserved=[record(p) for p in root.rglob('*') if p.is_file() and not p.is_symlink() and p!=root/'store/HEAD']
        symlinks=[dict(path=str(p),target=os.readlink(p)) for p in root.rglob('*') if p.is_symlink()]
        save(root/'restore_inputs.json',dict(recorded_at=stamp(),explicit_generation=37,
            selected=record(root/'selection_initial.bin'),expected_current=record(root/'selection_uncertain.bin'),
            prior_independent_assessment=record(root/'publish_independent_assessment.json'),preserved=preserved,symlinks=symlinks))
        second_command=[build['ELF']['path'],'--restore',str(root),trusted]
        second=run(root,'restore',second_command,root,180,240,1);second_report=json.loads((root/'restore.stdout').read_text())
        verify_native(second_report,second,second_command)
        for item in preserved+anchors+build['support']+build['source_snapshots']+[execution['plan'],build['ELF'],build['source_manifest'],build['build_runtime_manifest']]:
            assert record(item['path'])==item
        assert all(Path(x['path']).is_symlink() and os.readlink(x['path'])==x['target'] for x in symlinks)
        result=assess_restore(root,publish,second_report)
        result.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),frozen_execution=record(root/'frozen_execution.json'),
            publish_phase=first,restore_phase=second,publish_assessment=record(root/'publish_independent_assessment.json'),
            restore_inputs=record(root/'restore_inputs.json'),actual_compiler_dependency_count=build['actual_compiler_dependency_count'],
            support_original_and_copy_records=build['support_original_and_copy_records'])
        save(root/'engineering_receipt.json',result)
        save(root/'terminal.json',dict(status=result['status'],finished_at=stamp(),return_code=0,receipt=record(root/'engineering_receipt.json'),
            physical_advances=0,accepted_step_writer_qualified=False,full_restart_qualified=False))
        print(json.dumps(dict(status=result['status'],durable_publishes=3,controlled_uncertain_publish=1,expected_refusals=19)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),
            physical_advances=0,accepted_step_writer_qualified=False,full_restart_qualified=False,S3_complete=False))
        raise


if __name__=='__main__':main()

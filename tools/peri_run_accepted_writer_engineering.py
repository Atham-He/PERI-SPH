"""New bounded five-phase saved-record writer engineering,zero nonlinear solves."""
import argparse,json,shutil,subprocess,traceback
from pathlib import Path
from peri_run_body_archive_engineering import run,reserves,save,identity,stamp
from peri_bounded_body_archive_audit import record
from peri_archive_component_build import prepare,compile_and_freeze
from peri_run_generation_store_engineering import verify_native
from peri_bounded_accepted_writer_audit import create,assess_phase,final_assessment


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);args=parser.parse_args()
    root=Path(args.root).resolve();repo=Path(__file__).resolve().parents[1]
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    reserves(root);assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)=='';root.mkdir(mode=0o700)
    try:
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        old=root.parent/'engineering-bounded-mechanical-stream'
        old_receipt=json.loads((old/'engineering_receipt.json').read_text())
        assert old_receipt['status']=='converted_existing_stream_all_records_original_proposals_and_prefix_bytes_independently_matched_NOT_new_physical_restart'
        original=root.parent/'study-v19-early-friction-cycle/actual_PD_early_biaxial_slip_cycle.input'
        shutil.copy2(original,root/'original.input')
        active=root.parent/'cross-runtime-v19-raw/plan.json';active_plan=json.loads(active.read_text())
        protected=[x for x in active_plan['artifacts'] if x['path'].startswith('/home/')]
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        anchors=[record(old/name) for name in ['engineering_receipt.json','physical.phy','session_candidate.bin']]+[record(original),record(production),record(active)]+protected
        assert all(record(x['path'])==x for x in anchors)
        files=['tools/bounded_accepted_writer_checks.cpp','src/shared/peridynamics/pd_bounded_accepted_writer.h',
            'src/shared/peridynamics/pd_bounded_generation_store.h','tools/peri_bounded_accepted_writer_audit.py',
            'tools/peri_run_accepted_writer_engineering.py','tools/peri_archive_component_build.py','tools/peri_run_generation_store_engineering.py',
            'tools/peri_run_body_archive_engineering.py','tools/peri_bounded_body_archive_audit.py','tools/peri_bounded_generation_store_audit.py',
            'tools/peri_bounded_session_archive_audit.py','tools/peri_bounded_mechanical_stream_audit.py','tools/peri_canonical_archive_audit.py',
            'docs/peri/BOUNDED_ACCEPTED_WRITER_20261105.md']
        prepared=prepare(root,repo,files,'tools/bounded_accepted_writer_checks.cpp','accepted')
        modes=['--faults','--continuous','--split-first','--resume','--rollback']
        plan=dict(prepared_at=stamp(),source_commit=source,source_snapshots=prepared['source_snapshots'],anchors=anchors,
            compile_command=prepared['compile_command'],link_command=prepared['link_command'],build_cwd=str(prepared['builddir']),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),native_budget_per_phase=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=1,threads=1),
            modes=modes,cost_basis='Prior converted prefix3walks+25controls22.487s and two transaction phases15.604/7.916s. New phases each <=~80saved endpoint proposal checks plusprefix IO,without Newton;180CPU/240wall prospectively bounds each phase,not an old failed run extension.',
            record_scope='104saved-record COMMIT replays (40continuous+16split-first+24resume+24rollback),0new nonlinear physical solves. Old producer identity preserved separately;new checker actual ELF/source/build-runtime/input contract independently frozen.',
            file_bound_MiB=256,total_bound_GiB=1,host_reserve_GiB=120,output_reserve_GiB=150,
            gates=['Actual /proc/self/exe and source/build-runtime/original-input hashes match bound contract,current FP profile,unchanged descriptor identities before append',
                'All next arrays/frame/cursors/hashes ready beforewrite;direct append readback and fsync success thenonly noexcept candidate ownership+primitivecursor commit',
                '13primaryfault controls and18blocked continuations;partial/fulluncommitted tails preserved,oldaccepted archive/HEADunchanged,no truncate/reset/fallback',
                'Independent struct/hashlib expected records/candidateprefixes;continuous/reopened selected step16/split/tailrollback exactsame full G100 archive andphysical stream bytes',
                'Each concrete phase independent byte assessment passes andnew selection input hashes freeze before dependent next process',
                'All prior physical raw,production andactive independent reference artifacts unchanged;new physical CP/S3/S4 not released'],
            failure_policy='Any unexpected execution/integrity/resource/accuracy failure preserves exact exception/returncode/tail/lastoldcursor and stops later modes. No redispatch or physical solver replay,reset/fallback/gate relaxation.',
            new_physical_solver_calls=0,new_physical_steps=0,full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan);save(root/'state.json',dict(status='five_phase_saved_record_writer_attempted',identity=identity(__import__('os').getpid()),started_at=stamp()))
        build=compile_and_freeze(prepared)
        assert 'SPH::peridynamics::solveCoulombContact(' not in subprocess.check_output(['nm','-C',build['ELF']['path']],text=True)
        identities=dict(SourceManifest=build['source_manifest']['sha256'],BuildRuntimeManifest=build['build_runtime_manifest']['sha256'],ELF=build['ELF']['sha256'])
        fixture=create(root,old,root/'original.input',identities)
        for name in ['continuous','split','faults']:(root/name).mkdir()
        input_records=[record(x) for x in fixture['inputs']]+[record(root/'original.input')]
        save(root/'external_trust.json',dict(contract_SHA256=fixture['trusted'],identities=identities,
            existing_raw_producer_scope='Original v19 old frozen ELF/raw;new replay consumer identities do not replace original provenance.',new_physical_steps=0))
        execution=dict(recorded_at=stamp(),plan=record(root/'plan.json'),build=build,anchors=anchors,inputs=input_records,external_trust=record(root/'external_trust.json'))
        for x in anchors+build['support']+build['source_snapshots']+input_records:assert record(x['path'])==x
        save(root/'frozen_execution.json',execution)
        results=[];assessments=[]
        for i,mode in enumerate(modes):
            label=mode[2:].replace('-','_')
            frozen_inputs=[]
            for name in ['split_selection_16.bin','split_selection_100.bin','split_first_final.bin']:
                if (root/name).exists():frozen_inputs.append(record(root/name))
            save(root/(label+'.inputs.json'),dict(recorded_at=stamp(),mode=mode,inputs=frozen_inputs,prior_assessments=assessments))
            command=[build['ELF']['path'],mode,str(root),fixture['trusted']]
            native=run(root,label,command,root,180,240,1);report=json.loads((root/(label+'.stdout')).read_text())
            verify_native(report,native,command)
            for item in anchors+build['support']+build['source_snapshots']+input_records+frozen_inputs+[execution['plan'],execution['external_trust'],build['ELF'],build['source_manifest'],build['build_runtime_manifest']]:
                assert record(item['path'])==item
            assessed=assess_phase(root,fixture,mode,report);assessed.update(recorded_at=stamp(),native=native)
            save(root/(label+'.independent_assessment.json'),assessed);assessments.append(record(root/(label+'.independent_assessment.json')));results.append(assessed)
        result=final_assessment(root,fixture,results)
        result.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),frozen_execution=record(root/'frozen_execution.json'),
            phase_assessments=assessments,actual_compiler_dependency_count=build['actual_compiler_dependency_count'],support_original_and_copy_records=build['support_original_and_copy_records'])
        save(root/'engineering_receipt.json',result)
        save(root/'terminal.json',dict(status=result['status'],finished_at=stamp(),return_code=0,receipt=record(root/'engineering_receipt.json'),
            saved_record_commits=104,new_physical_steps=0,full_restart_qualified=False))
        print(json.dumps(dict(status=result['status'],saved_record_commits=104,primary_fault_controls=13,blocked_continuation_controls=18)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),new_physical_steps=0,full_restart_qualified=False,S3_complete=False))
        raise


if __name__=='__main__':main()

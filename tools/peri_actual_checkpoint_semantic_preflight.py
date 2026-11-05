"""NEW zero-step semantic validator controls;not nonzero physical restart."""
import argparse,json,shutil,subprocess,sys,traceback
from pathlib import Path
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_actual_checkpoint_negative import specification,run_controls
from peri_actual_checkpoint_support import initial_contract
from peri_bounded_session_archive_audit import encode
from peri_elastic_normal import read,record,verify,save,reserves,stamp


def main(root):
    repo=Path(__file__).resolve().parents[1];producer=root.parent/'engineering-actual-checkpoint-driver-preflight-v5'
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique engineering state required'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    component=verified_execution_component(root.parent)
    build=read(producer/'frozen_execution.json')['build'];proof=read(producer/'independent_preflight.json')
    assert proof['owned_restore_noadvance_control']['actual_solver_calls']==0
    verify(build['support']+build['source_snapshots']+[build['ELF']]);reserves(root.parent)
    root.mkdir(mode=0o700);frozen=root/'frozen';frozen.mkdir()
    try:
        artifacts=[]
        names=['peri_actual_checkpoint_campaign.py','peri_actual_checkpoint_negative.py','peri_actual_checkpoint_semantic_preflight.py',
               'peri_actual_checkpoint_support.py','peri_canonical_archive_audit.py','peri_bounded_session_archive_audit.py',
               'peri_bounded_mechanical_stream_audit.py','peri_elastic_normal.py']
        for name in names:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        for name in ['source_manifest.json','build_runtime_manifest.json','original.input','initial_candidate.bin','initial.phy','contract.bin']:
            source=producer/name;target=root/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        executable=frozen/'actual_CP_validator';shutil.copy2(build['ELF']['path'],executable);artifacts += [build['ELF'],record(executable)]
        initial=initial_contract(root);assert initial['trusted'].hex()==proof['contract_sha256']
        archive=encode(initial['fields'],7);derived=dict(checkpoints={0:dict(archive=archive,physical=initial['physical'])})
        critical=[record(executable)]+[record(root/name) for name in ['source_manifest.json','build_runtime_manifest.json','original.input']]
        critical += [x for x in build['support'] if '.so' in Path(x['path']).name]
        artifacts += [record(producer/'independent_preflight.json'),record(producer/'frozen_execution.json')]+critical
        artifacts += [component['receipt'],component['frozen_controller']]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            executable=str(executable),trusted_contract=initial['trusted'].hex(),artifacts=artifacts,critical_runtime=critical,
            negative_specification=specification(),logical_generation=7,physical_step=0,
            execution_component=component,
            scope='NEW readonly semantic validator controls from already zero-advance-qualified initial fixture. Correct SHA corruption,0Newton solves. Not physical7step history or nonzero restart.',
            controls=17,native_per_control=dict(CPU_seconds=900,wall_seconds=900,AS_GiB=1,threads=1),
            file_bound_GiB=4,total_output_GiB=4,host_F_reserve_GiB=120,output_reserve_GiB=150,
            expected_negative_codes=1,positive_binding_code=0,
            failure_policy='Any unexpected execution/identity/integrity/resource/accuracy error preserves exactraw/cause/returncode andstops remainingcontrols. Neverredispatch this state.',
            full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        verify(artifacts);save(root/'plan.json',plan)
        result=run_controls(root,plan,Execution(root,plan),derived,step=0,generation=7)
        reviewed=read(result['path']);assert reviewed['semantic_rejections']==16 and not reviewed['qualified_physical_nonzero_seed']
        verify(artifacts);verify(build['support']+build['source_snapshots'])
        receipt=dict(status='new_zero_step_actual_CP_validator_correct_SHA_semantics_and_binding_independently_passed_NOT_physical_restart',
            recorded_at=stamp(),source_commit=plan['source_commit'],review=result,plan=record(root/'plan.json'),
            positive_binding_controls=1,semantic_rejections=16,logical_generation=7,physical_step=0,
            native_solver_calls=0,full_restart_qualified=False,S3_complete=False)
        save(root/'engineering_receipt.json',receipt);save(root/'terminal.json',dict(status=receipt['status'],return_code=0,finished_at=stamp(),receipt=record(root/'engineering_receipt.json')))
        print(json.dumps(dict(status=receipt['status'],semantic_rejections=16,native_solver_calls=0)))
    except Exception as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),native_solver_calls=0,S3_complete=False));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True,type=Path);args=parser.parse_args();main(args.root)

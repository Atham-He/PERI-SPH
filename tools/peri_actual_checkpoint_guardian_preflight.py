"""NEW execution-guardian controls only;zero material/native solver calls."""
import argparse,hashlib,json,os,shutil,subprocess,sys,traceback
from pathlib import Path
from peri_actual_checkpoint_campaign import Execution
from peri_elastic_normal import record,save,reserves,stamp


def main(root):
    repo=Path(__file__).resolve().parents[1]
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(root.parent);root.mkdir(mode=0o700);frozen=root/'frozen';frozen.mkdir()
    try:
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip();artifacts=[]
        for name in ['peri_actual_checkpoint_campaign.py','peri_actual_checkpoint_reference.py','peri_actual_checkpoint_guardian_preflight.py','peri_elastic_normal.py']:
            original=repo/'tools'/name;copy=frozen/name;shutil.copy2(original,copy);artifacts += [record(original),record(copy)]
        original=repo/'docs/peri/ACTUAL_CHECKPOINT_CAMPAIGN_20261105.md';copy=frozen/original.name
        shutil.copy2(original,copy);artifacts += [record(original),record(copy)]
        cross=root.parent/'cross-runtime-v19-raw/frozen'
        for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py']:
            original=cross/name;copy=frozen/name;shutil.copy2(original,copy);artifacts += [record(original),record(copy)]
        artifacts += [record(p) for p in sorted((cross/'reference').glob('*.py'))]
        artifacts.append(record(Path(sys.executable).resolve()))
        programs={
            'fast':"""import pathlib,os,sys,json
p=pathlib.Path('/proc/self');f=p.joinpath('stat').read_text().rsplit(')',1)[1].split()
print(json.dumps(dict(type='synthetic_guardian_identity',pid=os.getpid(),ppid=os.getppid(),start_ticks=f[19],boot_id=pathlib.Path('/proc/sys/kernel/random/boot_id').read_text().strip(),argv=[sys.executable]+sys.argv)))
""",
            'failure':"""import sys,json
print(json.dumps(dict(type='checkpoint_driver_failure',synthetic_control_only=True,mechanical_durable_step=7,last_full_debug_step=6,actual_solver_calls=0,message='controlled synthetic witness')),file=sys.stderr)
sys.exit(1)
""",
            'wall':"""import time
time.sleep(10)
"""}
        for name,text in programs.items():
            p=frozen/(name+'.py');p.write_text(text);artifacts.append(record(p))
        save(root/'plan.json',dict(prepared_at=stamp(),source_commit=source,artifacts=artifacts,
            controls=['fast process birth+full emitted identity','fixed pure backend/95digit context binding','controlled exit1/errorwitness','duplicate dispatch refused before launch','controlled1s wall stop of sleeping synthetic process'],
            file_bound_GiB=4,total_output_GiB=4,host_reserve_GiB=120,output_reserve_GiB=150,
            control_CPU_seconds=5,control_AS_GiB=1,control_wall_seconds=5,controlled_wall_negative_seconds=1,
            scope='New guardian only;no actual CPP/material/reference/CP trajectory or old review calls. Deliberate exit1/wall stop are predeclared negatives,not real study failures.',
            failure_policy='Unexpected control/identity/integrity/resource failure stops this unique engineeringstate andpreserves raw/returncode;no redispatch.',
            physical_solver_calls=0,S3_complete=False,full_restart_qualified=False))
        assert all(record(x['path'])==x for x in artifacts)
        execution=Execution(root,{})
        fast=execution.run('fast',['/usr/bin/python3',str(frozen/'fast.py')],wall=5,cpu=5,gib=1)
        identity=json.loads((root/'fast.stdout').read_text())
        assert all(identity[k]==v for k,v in fast['birth_identity'].items()) and identity['argv']==fast['command']
        if fast['child_identity']:assert all(identity[k]==v for k,v in fast['child_identity'].items())
        binding=execution.run('reference_binding',['/usr/bin/python3',str(frozen/'peri_actual_checkpoint_reference.py'),'--selftest-binding',str(cross/'reference')],wall=5,cpu=5,gib=1)
        bound=json.loads((root/'reference_binding.stdout').read_text())
        assert bound==dict(status='new_actual_CP_reference_binding_and_95digit_context_preflight_passed',backend=str(frozen/'_pydecimal.py'),
            precision=95,C_decimal_loaded=False,native_solver_calls=0,material_calculations=0,old_raw_reviews_repeated=0)
        failed=execution.run('controlled_failure',['/usr/bin/python3',str(frozen/'failure.py')],wall=5,cpu=5,gib=1,expected_code=1)
        assert failed['native_failure_witnesses']==[dict(type='checkpoint_driver_failure',synthetic_control_only=True,
            mechanical_durable_step=7,last_full_debug_step=6,actual_solver_calls=0,message='controlled synthetic witness')]
        before={p.name:record(p) for p in root.iterdir() if p.is_file()}
        try:execution.run('fast',fast['command'],wall=5,cpu=5,gib=1)
        except AssertionError as error:assert 'never repeat a phase' in str(error)
        else:raise AssertionError('duplicate dispatched')
        assert before=={p.name:record(p) for p in root.iterdir() if p.is_file()}
        try:execution.run('controlled_wall',['/usr/bin/python3',str(frozen/'wall.py')],wall=1,cpu=5,gib=1)
        except AssertionError:pass
        else:raise AssertionError('wall control did not stop')
        wall=json.loads((root/'controlled_wall.attempt.json').read_text())
        assert wall['return_code']==-9 and 'wall budget' in wall['stop_reason']
        assert all(record(x['path'])==x for x in artifacts)
        result=dict(status='new_CP_execution_guardian_birth_identity_failure_duplicate_and_wall_controls_passed_NOT_physical_restart',
            recorded_at=stamp(),source_commit=source,control_count=5,fast=record(root/'fast.attempt.json'),reference_binding=record(root/'reference_binding.attempt.json'),controlled_failure=record(root/'controlled_failure.attempt.json'),
            controlled_wall=record(root/'controlled_wall.attempt.json'),plan=record(root/'plan.json'),
            native_solver_calls=0,material_calculations=0,old_physical_reviews_repeated=0,full_restart_qualified=False,S3_complete=False)
        save(root/'engineering_receipt.json',result);save(root/'terminal.json',dict(status=result['status'],return_code=0,finished_at=stamp(),receipt=record(root/'engineering_receipt.json')))
        print(json.dumps(dict(status=result['status'],control_count=5,native_solver_calls=0)))
    except Exception as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),native_solver_calls=0,S3_complete=False));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True,type=Path);args=parser.parse_args();main(args.root)

"""Prospective new streamed-index engine structural admission,zero old raw."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;component=base/'engineering-s4-incremental-typed-bytes-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    end=read(component/'terminal.json');assert end['return_code']==0 and end['receipt']==record(component/'independent_review.json')
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_streamed_raw_index.py','test_streamed_raw_index.py','run_streamed_index_controls.py','peri_streamed_index_preflight.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']+[record(component/n) for n in ['plan.json','independent_review.json','terminal.json']]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,methods=4,NEW_synthetic_nominal_ticks=[1,3,5,20],
            CPU_seconds=60,wall_seconds=120,AS_GiB=1,threads=1,
            scope='Distinct streamed index structural engine,exact full JSON bytes equal legacy oracle on NEW synthetic nondivisor grids.Full T payload mutation must change row/prefix/wholeSHA.No old8/13/18/14/stress rerun.',
            old_index_crash_root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0,
            actual_index_or_material_dispatch_admitted=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Retain exact evidence,stop dependencies,no same-state retry or gate/precision/payload relaxation.')
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_streamed_index_structural_attempted',started_at=stamp(),identity=identity()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_streamed_index_controls.py')],wall=120,cpu=60,gib=1)
        q=read(root/'controls.stdout');assert q['tests_run']==4 and q['errors']==q['failures']==0
        assert q['precision']==95 and not q['C_decimal_loaded'] and q['physical_solver_calls']==q['material_calculations']==q['old_raw_reads']==0
        verify(artifacts);save(root/'independent_review.json',dict(status=q['status'],plan=record(root/'plan.json'),attempt=record(root/'controls.attempt.json'),
            review=a['stdout'],methods=4,complete_index_JSON_bytes_equal_on_four_NEW_grids=True,full_T_mutation_detected=True,
            old_index_crash_root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0,
            actual_index_or_material_dispatch_admitted=False,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=4,root_cause_resolved=False)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))

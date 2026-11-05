"""Prospective component admission before actual existing-artifact assessment."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;component=base/'engineering-s4-streamed-index-engine-v2'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen);artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        names=['peri_buffered_metadata.py','test_buffered_metadata.py','peri_existing_index_assessment.py','test_existing_index_assessment.py',
            'run_existing_index_assessment_controls.py','peri_existing_index_controls_preflight.py']
        for name in names:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']+[record(component/n) for n in ['plan.json','independent_review.json','terminal.json']]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,methods=11,buffered_and_pin_controls=6,existing_artifact_controls=5,
            CPU_seconds=60,wall_seconds=120,AS_GiB=1,threads=1,buffer_bytes=2*2**20,
            scope='NEW synthetic pin-dedup/buffered UTF8 line and existing-artifact metadata assessment controls;no rebuilding old index/no old tests or original raw reading.',
            physical_solver_calls=0,material_calculations=0,old_raw_reads=0,actual_metadata_review_admitted=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Keep exact failure and stop dependencies,no retry old root or relax checks/resources.')
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_existing_artifact_component_attempted',started_at=stamp(),identity=identity()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_existing_index_assessment_controls.py')],wall=120,cpu=60,gib=1)
        q=read(root/'controls.stdout');assert q['tests_run']==11 and q['errors']==q['failures']==0 and q['precision']==95 and not q['C_decimal_loaded']
        assert q['physical_solver_calls']==q['material_calculations']==q['old_raw_reads']==0
        verify(artifacts);save(root/'independent_review.json',dict(status=q['status'],plan=record(root/'plan.json'),attempt=record(root/'controls.attempt.json'),review=a['stdout'],
            methods=11,all_unique_pins_preserved_conflicts_rejected=True,buffered_all_bytes_preserved=True,complete_existing_index_payload_controls=True,
            physical_solver_calls=0,material_calculations=0,old_raw_reads=0,actual_metadata_review_admitted=False,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=11)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))

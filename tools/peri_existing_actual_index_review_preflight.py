"""New complete existing-artifact assessment,not producer restart/extension."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_buffered_metadata import verify_unique


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;component=base/'engineering-s4-existing-index-assessment-v1'
    producer=base/'analysis-s4-actual-streamed-index-v2';old=base/'study-v23-physical-time-complete-history'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    end=read(component/'terminal.json');qa=read(component/'independent_review.json')
    assert end['return_code']==0 and end['receipt']==record(component/'independent_review.json') and qa['methods']==11
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_review_existing_actual_index.py','peri_existing_actual_index_review_preflight.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']+[record(component/n) for n in ['plan.json','independent_review.json','terminal.json']]
        artifacts+=read(producer/'plan.json')['artifacts']
        native=read(old/'continuous.attempt.json');producer_attempt=read(producer/'index.attempt.json')
        assert producer_attempt['return_code']==-9 and 'wall budget' in producer_attempt['stop_reason']
        index=record(producer/'accepted_binary_index.json');assert index['size']==1725554 and index['sha256']=='ca0281850f9aef4528765d8728cc33c445c06596536fb72d42ecc67c7ee1f981'
        fronts=[record(producer/'frontiers'/f'{k}.json') for k in range(81)]
        input_records=[native['stdout'],native['stderr'],record(old/'continuous.attempt.json'),record(old/'plan.json'),
            record(producer/'plan.json'),record(producer/'index.attempt.json'),record(producer/'failure.json'),record(producer/'stage_witness.jsonl'),index]+fronts
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,input_records=input_records,execution_component=execution,raw=native['stdout'],native_attempt=record(old/'continuous.attempt.json'),
            original_plan=record(old/'plan.json'),producer_attempt=record(producer/'index.attempt.json'),index_artifact=index,frontiers=fronts,
            CPU_seconds=180,wall_seconds=240,AS_GiB=8,threads=1,buffer_bytes=2*2**20,
            prior_metadata_attempts=[dict(root=str(base/'analysis-s4-actual-durable-index-v1'),return_code=-11,CPU_reservation_seconds=180,CPU_final_unknown=True,wall_cap_seconds=240),
                dict(root=str(producer),return_code=-9,CPU_reservation_seconds=180,CPU_final_unknown=True,wall_cap_seconds=240,wall_actual_seconds=240.069058337)],
            prospective_metadata_CPU_reservation_total_seconds=540,
            scope='NEW independent complete assessment of already existing index artifact.All81frontier/fullnative row typedSHA/line/prefix/offset/time/birth/wholeend/source/support;no builder/native/material calls.Preserves both producers failed and budgets unchanged.',
            semantic_checks_not_removed=True,all_unique_pins_before_and_after=True,
            actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Retain exact failed attempt/last durable verified frontier/stage/code/raw SHA and stop dependents;no same-root replay,producer rebuild,timeout extension or partialpass.')
        verify_unique(artifacts+input_records);save(root/'plan.json',plan)
        save(root/'state.json',dict(status='new_existing_index_independent_assessment_attempted',started_at=stamp(),identity=identity()))
        Execution(root,plan).run('review',[sys.executable,'-B',str(frozen/'peri_review_existing_actual_index.py'),str(root)],wall=240,cpu=180,gib=8)
        result=read(root/'independent_review.json')
        assert result['assessment']['accepted_steps']==80 and result['assessment']['frontiers']==81 and result['assessment']['point_metadata']==11847
        assert result['index_builder_called'] is False and result['physical_solver_calls']==result['material_calculations']==0
        verify_unique(artifacts+input_records)
        save(root/'terminal.json',dict(status=result['status'],return_code=0,finished_at=stamp(),attempt=record(root/'review.attempt.json'),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=result['status'],CPU_seconds=result['CPU_seconds'],assessment=result['assessment'],producer_status_changed=False)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))

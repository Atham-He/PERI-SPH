"""Distinct prospectively frozen metadata admission after byte-equivalence QA."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;component=base/'engineering-s4-streamed-index-engine-v2';old=base/'study-v23-physical-time-complete-history'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    end=read(component/'terminal.json');q=read(component/'independent_review.json')
    assert end['return_code']==0 and end['receipt']==record(component/'independent_review.json') and q['methods']==4 and q['full_T_mutation_detected']
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_actual_streamed_index_v2.py','peri_actual_streamed_index_preflight_v2.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']+[record(component/n) for n in ['plan.json','independent_review.json','terminal.json']]
        I=Path('/mnt/e/PERI-SPH/PER_sph_casule/framework_work/2026-11-05_S4_convergence_stage/incidents/20261105T113813Z-actual-index-sigsegv')
        diagnostic=base/'engineering-s4-synthetic-json-heap-diagnostic-v1';forensic=base/'forensics-s4-core-gdb-v1'
        artifacts += [record(I/n) for n in ['incident_receipt.json','failed_index_root_size_SHA.json']]
        artifacts += [record(diagnostic/n) for n in ['plan.json','independent_diagnostic.json','terminal.json']]
        artifacts += [record(forensic/n) for n in ['core_unwind_plan.json','core_unwind.attempt.json','core_unwind.stdout','debuginfo_fetch.attempt.json']]
        native=read(old/'continuous.attempt.json')
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,raw=native['stdout'],native_attempt=record(old/'continuous.attempt.json'),original_plan=record(old/'plan.json'),
            CPU_seconds=180,wall_seconds=240,AS_GiB=8,threads=1,
            prior_failed_metadata_reservation=dict(root=str(base/'analysis-s4-actual-durable-index-v1'),CPU_reserved_seconds=180,CPU_final_unknown=True,
                wall_reserved_seconds=240,wall_actual_seconds=110.019292049,return_code=-11,old_attempt_immutable=True),
            prospective_total_metadata_CPU_reserved_seconds=360,
            admission_reason='Distinct deterministic streamed typed-byte engine after8type-byte and4NEWgrid equivalence controls plus scoped synthetic diagnostics;not an exception-triggered backend fallback or a retry of old plan/source/state.',
            unchanged_old_attempt_budget=True,unchanged_current_CPU_wall_limits=True,old_SIGSEGV_root_cause_resolved=False,
            scope='Readonly complete native80 raw binary index;full legacy typed-byte/hash format,all original fields/trials retained;faulthandler+durable stage/all81frontiers.No material equations/native run.',
            physical_solver_calls=0,material_calculations=0,actual_material_equation_reviews=0,actual_full_material_audit_admitted=False,
            actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='On any execution/integrity/resource/metadata failure preserve exact code/stage/frontier/core and stop dependents.No same-root retry,partialpass,old budget change,lower precision or omitted T.')
        verify(artifacts+[plan['raw'],plan['native_attempt'],plan['original_plan']]);save(root/'plan.json',plan)
        save(root/'state.json',dict(status='new_streamed_actual_index_attempted',identity=identity(),started_at=stamp()))
        Execution(root,plan).run('index',[sys.executable,'-B',str(frozen/'peri_actual_streamed_index_v2.py'),str(root)],wall=240,cpu=180,gib=8)
        result=read(root/'independent_review.json');assert result['accepted_steps']==80 and result['material_point_metadata']==11847 and result['durable_frontiers']==81
        assert result['physical_solver_calls']==result['material_calculations']==result['actual_material_equation_reviews']==0
        verify(artifacts+[plan['raw'],plan['native_attempt'],plan['original_plan']])
        save(root/'terminal.json',dict(status=result['status'],return_code=0,finished_at=stamp(),attempt=record(root/'index.attempt.json'),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=result['status'],steps=80,CPU_seconds=result['CPU_seconds'],index=result['index'],old_SIGSEGV_root_cause_resolved=False)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))

"""One new frozen byte-equivalence component,zero original raw/math replay."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;component=base/'engineering-s4-durable-coordinator-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_streaming_typed_json.py','test_streaming_typed_json.py','run_streaming_typed_controls.py','peri_streaming_typed_preflight.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']
        diagnostic=base/'engineering-s4-synthetic-json-heap-diagnostic-v1'
        artifacts += [record(diagnostic/n) for n in ['plan.json','independent_diagnostic.json','terminal.json']]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,methods=8,large_matrix_unknowns=447,CPU_seconds=60,wall_seconds=120,AS_GiB=1,threads=1,
            scope='New incremental exact canonical(pack(value)) bytes/hash component only;legacy oracle on NEW synthetic values/alltypes/nonsymmetric447matrix/context.No old13/18/14/8stress repeated.',
            known_index_crash_root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0,
            actual_index_or_material_audit_admitted=False,original_frozen_encoders_unchanged=True,S4_complete=False,calibration_allowed=False,
            failure_policy='Preserve exact failure/returncode/input/output and stop dependent use,no same-root retry/format change/drop of payload or numerical gates.')
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_streaming_typed_component_attempted',identity=identity(),started_at=stamp()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_streaming_typed_controls.py')],wall=120,cpu=60,gib=1)
        q=read(root/'controls.stdout');assert q['tests_run']==8 and q['errors']==q['failures']==0
        assert q['precision']==95 and not q['C_decimal_loaded'] and q['physical_solver_calls']==q['material_calculations']==q['old_raw_reads']==0
        verify(artifacts);save(root/'independent_review.json',dict(status=q['status'],plan=record(root/'plan.json'),attempt=record(root/'controls.attempt.json'),
            review=a['stdout'],all_declared_bytes_and_hashes_equal=True,methods=8,known_index_crash_root_cause_resolved=False,
            actual_index_or_material_audit_admitted=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=8,CPU_seconds=q['CPU_seconds'],root_cause_resolved=False)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))

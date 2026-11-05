"""One-shot new structural controls;no original math or native replay."""
from pathlib import Path
import argparse,ast,json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent
    old=base/'study-v23-physical-time-complete-history';state=base/'engineering-s4-durable-reference-v1';budget=base/'engineering-s4-reference-budget-owner-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique unused root only'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True);artifacts=[]
        for source in sorted((old/'frozen/reference').glob('*.py')):
            target=reference/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py','peri_physical_checkpoint_reference.py']:
            source=old/'frozen'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        for source,target in [(state/'frozen/reference/peri_reference_checkpoint.py',reference/'peri_reference_checkpoint.py'),
                              (state/'frozen/reference_checkpoint_fixture.py',frozen/'reference_checkpoint_fixture.py'),
                              (budget/'frozen/peri_reference_budget.py',frozen/'peri_reference_budget.py')]:
            shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        names=['peri_durable_reference_stream.py','peri_durable_reference_launch.py','peri_durable_reference_worker.py',
               'peri_durable_reference_summary.py','durable_reference_fixture.py','durable_reference_launch_child.py',
               'test_durable_reference_coordinator.py','run_durable_reference_coordinator_controls.py','peri_durable_reference_coordinator_preflight.py']
        for name in names:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        tests=[n.name for n in ast.walk(ast.parse((frozen/'test_durable_reference_coordinator.py').read_text())) if isinstance(n,ast.FunctionDef) and n.name.startswith('test_')]
        assert len(tests)==14
        artifacts+=read(base/'engineering-s4-physical-header-v1/plan.json')['software']
        prior=[record(p/n) for p in [state,budget] for n in ['plan.json','independent_review.json','terminal.json']]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,prior_component_receipts=prior,controls=tests,declared_methods=14,
            stream_methods=9,launch_methods=5,synthetic_steps=8,correctSHA_clock_identity_cases=4,
            CPU_seconds=60,wall_seconds=180,AS_GiB=1,threads=1,
            scope='NEW binary-index/accepted-boundary loop/explicit store+resourceHEAD/launch handshake structure only. Synthetic memory steps,real bounded metadata children.No real trajectory raw or equations.',
            worker_scope='Worker implemented with original Inspector methods and AST-identical original finalization body;actual full-material execution NOT admitted by this preflight.',
            comparison_policy='Preserve all original snapshot bytes and actual cpu_completed;comparison projection explicitly excludes ONLY cpu_completed and retains original/derivedSHA plus excluded fields/measurement.',
            previous_13_and_18_controls_reexecuted=False,physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
            actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Retain exact attempt/returncode/raw/fixtures and stop dependent use.No same-root retry,oldsource reset or gate relaxation.')
        verify(artifacts+prior);save(root/'plan.json',plan)
        save(root/'state.json',dict(status='new_durable_coordinator_structural_attempted',started_at=stamp(),identity=identity()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_durable_reference_coordinator_controls.py')],wall=180,cpu=60,gib=1)
        q=read(root/'controls.stdout');assert q['tests_run']==14 and q['failures']==q['errors']==0
        assert q['precision']==95 and not q['C_decimal_loaded']
        assert q['physical_solver_calls']==q['material_calculations']==q['real_trajectory_raw_reads']==0
        sem=read(root/'work/test_correctSHA_wrong_step_clock_and_native_birth_rejected/correctSHA_controls.json')
        assert len(sem)==4 and all(x['rejected'] is True for x in sem)
        good=read(root/'work/test_launch_is_reserved_born_released_acknowledged_then_wait4_settled/attempt/terminal.json')
        bad=read(root/'work/test_nonzero_complete_child_records_exact_exit_and_halts/attempt/terminal.json')
        early=read(root/'work/test_exit_before_handshake_retains_reservation_and_known_code/attempt/failure.json')
        assert (good['return_code'],bad['return_code'],early['return_code'])==(0,9,7)
        evidence=[record(p) for p in sorted((root/'work').rglob('*')) if p.is_file() and not p.is_symlink()]
        verify(artifacts+prior)
        save(root/'independent_review.json',dict(status=q['status'],finished_at=stamp(),plan=record(root/'plan.json'),
            attempt=record(root/'controls.attempt.json'),review=a['stdout'],controls=tests,correctSHA_semantic_cases=4,
            real_launch_control_exit_codes=[0,9,7],control_evidence=evidence,original_finalization_body_AST_identical=True,
            explicit_CP_selection_and_resource_HEAD_structural_qualified=True,previous_13_and_18_controls_reexecuted=False,
            actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,
            physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(status=q['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=14,real_launch_exit_codes=[0,9,7],physical_solver_calls=0)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),physical_solver_calls=0));raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root)

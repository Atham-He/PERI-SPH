"""One new frozen engineering attempt;no trajectory/raw re-audit or solver."""
from pathlib import Path
import argparse
import ast
import json
import shutil
import subprocess
import sys
import traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent
    old=base/'study-v23-physical-time-complete-history'
    component=base/'engineering-s4-durable-reference-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique new root only'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True);artifacts=[]
        for source in sorted((old/'frozen/reference').glob('*.py')):
            target=reference/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py']:
            source=old/'frozen'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        # Reuse qualified state component bytes;never execute its13 terminal tests.
        for name,where in [('peri_reference_checkpoint.py',reference),('reference_checkpoint_fixture.py',frozen)]:
            source=component/'frozen'/('reference' if name=='peri_reference_checkpoint.py' else '')/name
            target=where/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        names=['peri_reference_budget.py','test_reference_budget_owner.py','reference_budget_child.py',
               'run_reference_budget_owner_controls.py','peri_reference_budget_owner_preflight.py']
        for name in names:
            source=repo/'tools'/name;target=frozen/name
            shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        tests=[n.name for n in ast.walk(ast.parse((frozen/'test_reference_budget_owner.py').read_text()))
               if isinstance(n,ast.FunctionDef) and n.name.startswith('test_')]
        assert len(tests)==18
        artifacts+=read(base/'engineering-s4-physical-header-v1/plan.json')['software']
        prior=[record(component/'plan.json'),record(component/'independent_review.json'),record(component/'terminal.json')]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,prior_component=prior,
            controls=tests,declared_methods=18,ledger_methods=15,real_class_memory_methods=3,
            correctSHA_semantic_cases=12,declared_child_kills=['before_HEAD','after_HEAD'],real_wait4_witnesses=1,
            CPU_seconds=60,wall_seconds=180,AS_GiB=1,threads=1,
            resource_policy='Separate append-only CAS ledger;explicit reserve/intent before launch;complete wait4 identity-bound settlement only;unknown fully charged;old mathematical selection cannot rewind cost;orphan fails closed;nonzero/over-reservation halts.',
            owner_scope='Real PhysicalTimeFrictionInspector and OneStepCapture classes populated with synthetic33node116contact memory,19member ABI,static identities and inherited close checked.No computed material trajectory or full reader restart qualification.',
            scope='Engineering controls only;synthetic resource/owner memory plus controlled wait4/SIGKILL children.No old13tests rerun/no actual trajectory raw parse.',
            physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
            actual_computed_trajectory_owner_qualified=False,actual_review_restart_qualified=False,
            hardware_power_loss_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Retain exact returncode/stdout/stderr/fixtures and frozen attempt;stop dependent use,no same-root redispatch or gate relaxation.')
        verify(artifacts+prior);save(root/'plan.json',plan)
        save(root/'state.json',dict(status='new_budget_owner_preflight_attempted',identity=identity(),started_at=stamp()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_reference_budget_owner_controls.py')],wall=180,cpu=60,gib=1)
        q=read(root/'controls.stdout')
        assert q['tests_run']==18 and q['failures']==q['errors']==0
        assert q['precision']==95 and q['C_decimal_loaded'] is False
        assert q['physical_solver_calls']==q['material_calculations']==q['real_trajectory_raw_reads']==0
        semantic=read(root/'work/test_correct_SHA_semantic_corruptions_leave_ledger_unchanged/semantic_controls.json')
        assert len(semantic)==12 and all(x['rejected'] is True for x in semantic)
        for name in ['test_real_child_kill_before_HEAD_fails_closed_with_full_charge','test_real_child_kill_after_HEAD_keeps_unsettled_charge']:
            assert read(root/'work'/name/'child.attempt.json')['return_code']==-9
        measured=read(root/'work/test_real_wait4_witness_complete_process_accounting/measured.wait4.json')
        assert measured['complete'] is True and measured['return_code']==0 and measured['CPU_ns']>0 and measured['wall_ns']>0
        evidence=[record(p) for p in sorted((root/'work').rglob('*')) if p.is_file() and not p.is_symlink()]
        aliases=[str(p) for p in (root/'work').rglob('*') if p.is_symlink()]
        verify(artifacts+prior)
        save(root/'independent_review.json',dict(status=q['status'],finished_at=stamp(),plan=record(root/'plan.json'),
            attempt=record(root/'controls.attempt.json'),review=a['stdout'],controls=tests,correctSHA_semantic_cases=12,
            actual_child_kills=2,real_wait4_witness=measured,control_evidence=evidence,explicit_invalid_fixture_aliases=aliases,
            cumulative_budget_ledger_component_qualified=True,real_Inspector_class_synthetic_memory_owner_qualified=True,
            actual_computed_trajectory_owner_qualified=False,actual_review_restart_qualified=False,
            physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
            hardware_power_loss_qualified=False,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(status=q['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=18,semantic_rejects=12,actual_child_kills=2,real_wait4=1,physical_solver_calls=0)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),physical_solver_calls=0));raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();main(a.root)

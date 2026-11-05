"""One-shot bounded synthetic reference checkpoint preflight,zero physics."""
from pathlib import Path
import argparse,ast,json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component

def inventory(path,cls):
    tree=ast.parse(Path(path).read_text());node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
    init=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    return {n.attr for n in ast.walk(init) if isinstance(n,ast.Attribute) and isinstance(n.ctx,ast.Store) and isinstance(n.value,ast.Name) and n.value.id=='self'}

def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;old=base/'study-v23-physical-time-complete-history'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'new unused preflight only'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True);artifacts=[]
        for source in sorted((old/'frozen/reference').glob('*.py')):
            target=reference/source.name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py']:
            source=old/'frozen'/name;target=frozen/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        for name in ['peri_reference_checkpoint.py','reference_checkpoint_fixture.py','reference_checkpoint_child.py',
                     'test_reference_checkpoint.py','run_reference_checkpoint_controls.py','peri_reference_checkpoint_preflight.py']:
            source=repo/'tools'/name;target=(reference if name=='peri_reference_checkpoint.py' else frozen)/name
            shutil.copy2(source,target);artifacts += [record(source),record(target)]
        observed=inventory(reference/'peri_bounded_gel_audit.py','Inspector')|inventory(reference/'peri_bounded_friction_audit.py','FrictionInspector')
        expected={'case','h','full','maximum','trial_points','ref','geometry','mass','I','E','C','ids','cache','G','ever_active','slip_signs','row_slip_signs','capture','statistics'}
        assert observed==expected and len(expected)==19,'actual frozen inspector member inventory changed'
        artifacts+=read(base/'engineering-s4-physical-header-v1/plan.json')['software']
        controls=[n.name for n in ast.walk(ast.parse((frozen/'test_reference_checkpoint.py').read_text())) if isinstance(n,ast.FunctionDef) and n.name.startswith('test_')]
        assert len(controls)==13
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,actual_inspector_member_inventory=sorted(observed),
            controls=controls,declared_methods=13,correctSHA_semantic_cases=12,declared_child_kills=['before_HEAD','after_HEAD'],
            synthetic_steps=8,synthetic_state_shape=dict(nodes=33,contacts=116),CPU_seconds=60,wall_seconds=180,AS_GiB=1,threads=1,
            scope='Exact typed pure95 memory/context/negativecapture and committed ancestor transaction COMPONENT only,with real childSIGKILL.No original audit/state/native/data retry.',
            missing_before_actual_use=['actual Inspector snapshot/install integration and no-advance state qualification','independent cumulative resource highwater ledger','new raw-index/source contract and fresh audit admission'],
            physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,old_math_reviews_repeated=False,
            hardware_power_loss_qualified=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Retain exact failed attempt/raw/fixtures and stop dependent use.No redispatch same state/fallback/precision or gate relaxation.')
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_reference_state_transaction_preflight_attempted',identity=identity(),started_at=stamp()))
        a=Execution(root,plan).run('controls',[sys.executable,'-B',str(frozen/'run_reference_checkpoint_controls.py')],wall=180,cpu=60,gib=1)
        q=read(root/'controls.stdout');assert q['tests_run']==13 and q['failures']==q['errors']==0 and q['precision']==95 and q['C_decimal_loaded'] is False
        assert q['physical_solver_calls']==q['material_calculations']==q['real_trajectory_raw_reads']==0
        evidence=[record(p) for p in sorted((root/'work').rglob('*')) if p.is_file()]
        for phase in ['test_actual_child_kill_before_HEAD_keeps_old_frontier','test_actual_child_kill_after_HEAD_exposes_only_committed_frontier']:
            attempt=read(root/'work'/phase/'child.attempt.json');assert attempt['return_code']==-9
        verify(artifacts);save(root/'independent_review.json',dict(status=q['status'],finished_at=stamp(),plan=record(root/'plan.json'),
            attempt=record(root/'controls.attempt.json'),review=a['stdout'],controls=controls,correctSHA_semantic_cases=12,actual_child_kills=2,
            complete_inspector_inventory=sorted(observed),control_evidence=evidence,physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
            actual_inspector_integration_qualified=False,cumulative_budget_ledger_qualified=False,actual_review_restart_qualified=False,
            hardware_power_loss_qualified=False,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(status=q['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],methods=13,correctSHA_semantic_cases=12,actual_child_kills=2,physical_solver_calls=0)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),physical_solver_calls=0));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();main(a.root)

"""One new metadata-only actual-raw admission,never native/material dispatch."""
from pathlib import Path
import argparse,json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;old=base/'study-v23-physical-time-complete-history'
    component=base/'engineering-s4-durable-coordinator-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique unused readonly metadata root'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    terminal=read(component/'terminal.json');review=read(component/'independent_review.json')
    assert terminal['return_code']==0 and terminal['receipt']==record(component/'independent_review.json')
    assert len(review['controls'])==14 and not review['actual_full_material_audit_admitted']
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_durable_real_index.py','peri_durable_real_index_preflight.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']
        artifacts += [record(component/'plan.json'),record(component/'independent_review.json'),record(component/'terminal.json')]
        incident=Path('/mnt/e/PERI-SPH/PER_sph_casule/framework_work/2026-11-05_S4_convergence_stage/incidents/20261105T1033Z-v23-host-reboot')
        artifacts += [record(incident/'incident_receipt.json'),record(incident/'root_size_SHA_index.json')]
        native=read(old/'continuous.attempt.json')
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,raw=native['stdout'],native_attempt=record(old/'continuous.attempt.json'),
            original_plan=record(old/'plan.json'),CPU_seconds=180,wall_seconds=240,AS_GiB=8,threads=1,
            scope='Readonly complete471838271B native80 binary index and separately coded81frontier prefix/line/offset/pointcount check only.No old raw equations assessed,no S3math or native replay.',
            physical_solver_calls=0,material_calculations=0,actual_material_equation_reviews=0,
            actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Preserve exact failure/returncode/raw/last complete metadata and stop dependents,no same-root retry or budget extension.')
        verify(artifacts+[plan['raw'],plan['native_attempt'],plan['original_plan']]);save(root/'plan.json',plan)
        save(root/'state.json',dict(status='new_actual_raw_index_attempted',identity=identity(),started_at=stamp()))
        a=Execution(root,plan).run('index',[sys.executable,'-B',str(frozen/'peri_durable_real_index.py'),'--root',str(root)],wall=240,cpu=180,gib=8)
        q=read(root/'independent_review.json')
        assert q['accepted_steps']==80 and q['material_point_metadata']==11847 and q['independent_all81offset_prefix_lineSHA_and_pointcounter_passed']
        assert q['physical_solver_calls']==q['material_calculations']==q['actual_material_equation_reviews']==0
        verify(artifacts+[plan['raw'],plan['native_attempt'],plan['original_plan']])
        save(root/'terminal.json',dict(status=q['status'],return_code=0,finished_at=stamp(),attempt=record(root/'index.attempt.json'),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=q['status'],steps=80,point_metadata=11847,CPU_seconds=q['CPU_seconds'],index=q['index'])))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),physical_solver_calls=0));raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root)

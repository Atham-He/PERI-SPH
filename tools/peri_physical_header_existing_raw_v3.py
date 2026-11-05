"""Unique new fixed-reader analysis of existing header raw; ZERO native calls."""
import argparse,json,shutil,subprocess,sys,traceback
from pathlib import Path
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_physical_session_preflight import assess
from peri_physical_header_analysis_inputs import copy_admission
from peri_elastic_normal import read,save,verify,record,reserves,stamp,identity

def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;native=base/'engineering-s4-physical-header-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique read-only state;never redispatch'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    old=read(native/'header_independent_review.json')
    assert old['status']=='new_v2_runtime_header_admission_failed_preserve_evidence'
    assert old['failure']=="TypeError(\"string indices must be integers, not 'str'\")" and old['native_solver_calls']==0
    for label in ['initialize','describe_initial']:
        a=read(native/(label+'.attempt.json'));assert a['return_code']==0 and a['stop_reason'] is None;verify([a['stdout'],a['stderr']])
    assert read(native/'header_reference.attempt.json')['return_code']==1
    reserves(base);component=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True)
        artifacts=[]
        for p in sorted((native/'frozen/reference').glob('*.py')):
            target=reference/p.name;shutil.copy2(p,target);artifacts += [record(p),record(target)]
        for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py']:
            p=native/'frozen'/name;target=frozen/name;shutil.copy2(p,target);artifacts += [record(p),record(target)]
        for name in ['peri_physical_header_readonly_v2.py','peri_physical_header_existing_raw_v3.py','test_peri_header_readonly_metadata.py']:
            p=repo/'tools'/name;target=frozen/name;shutil.copy2(p,target);artifacts += [record(p),record(target)]
        manifest=read(native/'header_reference_admission.json');artifacts += manifest['artifacts']
        artifacts += copy_admission(native,root)
        for name in ['peri_physical_header_analysis_inputs.py','test_peri_physical_header_analysis_inputs.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        prior_failed=base/'analysis-s4-header-v2-existing-raw'
        for name in ['plan.json','failure.json','header_reference.attempt.json','header_reference.stderr']:
            artifacts.append(record(prior_failed/name))

        artifacts += [record(native/n) for n in ['plan.json','failure.json','frozen_execution.json','header_independent_review.json',
            'header_reference.attempt.json','header_reference.stdout','header_reference.stderr','header_reference_admission.json',
            'initialize.attempt.json','initialize.stdout','describe_initial.attempt.json','describe_initial.stdout']]
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69';artifacts.append(record(production))
        plan=dict(prepared_at=stamp(),native_root=str(native),native_producer_source=read(native/'plan.json')['source_commit'],
                  source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                  auditor=record(frozen/'peri_physical_header_readonly_v2.py'),artifacts=artifacts,execution_component=component,
                  CPU_seconds=180,wall_seconds=240,AS_GiB=2,threads=1,normalized_gate=1e-10,
                  scope='NEW admission-path transport and pre-execution selfcheck; numeric reader/data/gates unchanged. Both original failures preserved.',
                  negative_metadata_QA=dict(methods=3,passed=True,tool_chunk='3ae5a1',seconds=.001),
                  path_metadata_QA=dict(methods=4,passed=True,tool_chunk='e54535',seconds=.003),
                  physical_solver_calls=0,native_dispatch_allowed=False,old_successful_math_reviews_repeated=0,
                  failure_policy='Preserve exact failure and stop dependents;do not re-run any native or same attempted audit state.',S4_complete=False,calibration_allowed=False)
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_existing_header_readonly_v3_path_admission_attempted',identity=identity(),started_at=stamp()))
        build=read(native/'frozen_execution.json')['build']
        rows=[json.loads(line) for line in (native/'initialize.stdout').read_text().splitlines()]
        initial=assess(native,base/'engineering-actual-checkpoint-driver-preflight-v5',build,rows[-1])
        executor=Execution(root,plan);a=executor.run('header_reference',[sys.executable,str(frozen/'peri_physical_header_readonly_v2.py'),'--root',str(root)],wall=240,cpu=180,gib=2)
        q=read(root/'header_independent_review.json');assert 'failure' not in q and q['maximum_normalized_error']<=1e-10
        assert q['new_PD_steps']==q['independent_trial_material_points']==0 and len(q['clock_header_negative_controls'])==6
        verify(artifacts)
        receipt=dict(status='existing_actual_v2_header_independently_qualified_by_NEW_path_admitted_reader_ZERO_steps',
                     finished_at=stamp(),plan=record(root/'plan.json'),review=record(root/'header_independent_review.json'),
                     attempt=record(root/'header_reference.attempt.json'),native_producer=record(native/'frozen_execution.json'),
                     old_failed_review=record(native/'header_independent_review.json'),old_failure_preserved=True,prior_readonly_path_failure=record(prior_failed/'failure.json'),
                     initial_complete_bytes=initial,physical_solver_calls=0,old_successful_math_reviews_repeated=0,
                     full_nonzero_restart_qualified=False,time_convergence_qualified=False,S4_complete=False,calibration_allowed=False)
        save(root/'independent_review.json',receipt);save(root/'terminal.json',dict(status=receipt['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json'),physical_solver_calls=0))
        print(json.dumps(dict(status=receipt['status'],maximum_header_error=q['maximum_normalized_error'],physical_solver_calls=0)))
    except BaseException as e:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(e),traceback=traceback.format_exc(),physical_solver_calls=0));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();main(a.root)

"""A NEW fixed pure-Python decimal reference of existing native raw.

No solver launch, exception-triggered backend switch, gate change or edit of
the original failed review. The full original reader and corruption controls
run in a fresh process with the declared backend selected before imports.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback


def select_backend(reference):
    assert 'decimal' not in sys.modules and 'fractions' not in sys.modules,'fresh fixed-backend process required'
    import _pydecimal as backend
    sys.modules['decimal']=backend
    sys.path.insert(0,str(reference))
    import fractions
    import peri_bounded_friction_audit as friction
    import peri_bounded_gel_audit as gel
    import peri_pd_material_reference as material
    import peri_gel_affine_audit as affine
    import peri_coulomb_contact_audit as contact
    assert fractions.Decimal is backend.Decimal
    assert all(m.D is backend.Decimal for m in [friction,gel,material,affine,contact])
    return backend


def selftest(reference):
    backend=select_backend(reference)
    from fractions import Fraction
    with backend.localcontext() as context:
        context.prec=95
        assert Fraction(backend.Decimal('.1'))==Fraction(1,10)
        assert backend.Decimal(4).sqrt()==2
        assert context.prec==95
    print(json.dumps(dict(status='fixed_pure_backend_binding_passed',module=backend.__file__,
        C_decimal_loaded='_decimal' in sys.modules,physical_solver_calls=0)))


def run(root):
    assert not (root/'audit_state.json').exists(),'never repeat a fixed-reference audit state'
    backend=select_backend(root/'frozen/reference')
    from peri_elastic_normal import read,record,verify,save,reserves,identity,stamp
    from peri_bounded_friction_audit import FrictionInspector,inspect_stream
    from peri_bounded_friction_negative import OneStepCapture,negatives
    plan=read(root/'plan.json');verify(plan['artifacts'])
    assert plan['reference_backend']=='fixed_pure_Python_decimal_before_import'
    assert record(Path(__file__).resolve())==plan['auditor']
    original=Path(plan['original_root']);case_name=plan['case']
    old_failure=read(original/(case_name+'.independent_review.json'))
    assert old_failure['failure']=="RuntimeError('internal error in flags_as_exception')"
    assert read(original/'terminal.json')['status']=='stopped_for_review'
    save(root/'audit_state.json',dict(status='attempted_new_fixed_reference',identity=identity(),started_at=stamp()))
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30))
    resource.setrlimit(resource.RLIMIT_CPU,(21600,21600))
    resource.setrlimit(resource.RLIMIT_FSIZE,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    start=time.monotonic();cursor={};original_step=FrictionInspector.step
    receipt=dict(status='fixed_pure_reference_audit_failed',identity=identity(),started_at=stamp(),
        physical_solver_calls=0,backend=backend.__file__,backend_selected_before_reference_imports=True,
        exception_triggered_backend_switch=False,original_failed_review_preserved=True,
        original_runtime_cause_resolved=False,original_v19_study_passed=False,
        S3_complete=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    def observed_step(self,row,old,z,loss,history):
        cursor.update(old_step=old['step'],recorded_step=row['state']['step'],time=row['state']['time'])
        result=original_step(self,row,old,z,loss,history)
        with (root/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(recorded_at=stamp(),cursor=cursor.copy(),trial_points=self.trial_points,
                maximum_normalized_error=self.maximum,CPU_seconds=time.process_time(),
                RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),allow_nan=False)+'\n')
            stream.flush();os.fsync(stream.fileno())
        print(json.dumps(dict(checked_steps=row['state']['step'],trial_points=self.trial_points)),flush=True)
        return result
    FrictionInspector.step=observed_step
    try:
        reserves(root)
        case=next(c for c in read(original/'plan.json')['controls'] if c['name']==case_name)
        attempt=read(original/(case_name+'.attempt.json'));verify([attempt['stdout'],attempt['stderr']])
        assert attempt['return_code']==0 and attempt['stop_reason'] is None
        assert attempt['cursor']['durability_confirmed'] and attempt['cursor']['last_durable_accepted_steps']==40
        with Path(attempt['stdout']['path']).open() as stream:meta=json.loads(next(stream))
        assert meta==attempt['native_identity'] and meta['argv']==attempt['command']
        assert meta['pid']==attempt['child_pid'] and meta['ppid']==attempt['identity']['pid']
        capture=OneStepCapture()
        receipt['assessment']=inspect_stream(case,attempt['stdout']['path'],0,True,capture)
        receipt['negative_controls']=negatives(capture,root/'single_step_negatives',attempt['stdout'])
        assert len(receipt['negative_controls'])==8 and all(x['rejected'] for x in receipt['negative_controls'])
        assert receipt['assessment']['all_trial_point_tangents_checked']
        assert all(n>0 for n in receipt['assessment']['temporally_reversed_contacts_per_tangent'])
        verify(plan['artifacts']);verify([attempt['stdout'],attempt['stderr']])
        receipt.update(status='existing_v19_raw_independently_qualified_by_new_fixed_pure_Python_reference',
            scope='This complete existing native40step raw and original accuracy/coverage/8corruption gates. Original C-runtime cause unresolved; old failed study unchanged. Not CP,convergence or tail qualification.',
            raw=attempt['stdout'],original_failed_review=record(original/(case_name+'.independent_review.json')))
    except Exception as error:
        receipt.update(failure=repr(error),failure_traceback=traceback.format_exc(),last_cursor=cursor.copy())
    finally:FrictionInspector.step=original_step
    receipt.update(finished_at=stamp(),wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time(),
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        plan=record(root/'plan.json'),audit_state=record(root/'audit_state.json'),
        progress=record(root/'progress.jsonl') if (root/'progress.jsonl').exists() else None)
    save(root/'independent_review.json',receipt)
    print(json.dumps(dict(status=receipt['status'],failure=receipt.get('failure'))),flush=True)
    return int('failure' in receipt)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path);parser.add_argument('--selftest-reference',type=Path)
    args=parser.parse_args()
    if args.selftest_reference:selftest(args.selftest_reference)
    else:sys.exit(run(args.root))

"""Fresh fixed95 full-equation worker;explicit admitted request only.

This entry never dispatches a native solver or selects a recovery frontier.
Unchanged original material/global/contact/geometry methods do the math.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback


def read_record(item):
    path=Path(item['path'])
    if path.is_symlink() or not path.is_file() or path.stat().st_size!=item['size']:
        raise ValueError('regular frozen input/size required')
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=item['sha256']:raise ValueError('frozen input SHA')
    return json.loads(data)


def run(request_path,request_sha):
    data=Path(request_path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=request_sha:raise ValueError('explicit admission request SHA')
    request=json.loads(data);plan=read_record(request['plan'])
    if request['material_audit_admitted'] is not True or plan['material_audit_admitted'] is not True:
        raise ValueError('full material audit not prospectively admitted')
    resource_contract=read_record(request['resource_contract'])
    from peri_durable_reference_launch import grant,atomic_new
    launch=Path(request['launch_directory'])
    grant(launch,request['resource_ledger'],resource_contract,request['attempt_id'])
    # The backend is selected once,before Fraction or any numerical reader.
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).resolve().parent/'reference');backend.getcontext().prec=95
    from peri_elastic_normal import verify,record,stamp
    from peri_reference_checkpoint import Store
    from peri_durable_reference_stream import DurableSession
    from peri_durable_reference_summary import original_summary
    from peri_physical_friction_negative import negatives
    from peri_physical_checkpoint_reference import require_full_reference
    verify(plan['artifacts']);index=read_record(request['index']);contract=read_record(request['state_contract'])
    if contract['raw_sha256']!=index['raw']['sha256'] or index['case']!=plan['case']:
        raise ValueError('raw index/state contract/case admission')
    verify([index['raw']]);target=request['end_step']
    if type(target) is not int or not 0<=target<=contract['total_steps']:raise ValueError('explicit segment endpoint')
    receipt=dict(status='new_durable_reference_failed',started_at=stamp(),native_solver_calls=0,
        backend_selected_before_reference_imports=True,exception_triggered_backend_switch=False,
        precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
        request=record(request_path),plan=request['plan'],raw=index['raw'],
        full_restart_qualified=False,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False)
    session=None;commits=[];cursor=None
    try:
        with backend.localcontext() as context:
            context.prec=95
            session=DurableSession(index['raw']['path'],index,contract,Store(request['store'],contract),
                request['expected_checkpoint_HEAD_sha256'],request['selected_leaf'])
            if target<session.step:raise ValueError('endpoint precedes explicit selected frontier')
            prefix=request['attempt_id']
            if target==session.step and request['selected_leaf'] is not None:
                commits.append(session.commit(prefix+'_noadvance',cpu=session.restored_CPU_diagnostic))
            else:
                if request['selected_leaf'] is None:commits.append(session.commit(prefix+'_s0'))
                while session.step<target:
                    session.advance();cursor=dict(step=session.step,time=session.old['time'])
                    commit=session.commit(prefix+'_s'+str(session.step));commits.append(commit)
                    atomic_new(launch/('accepted_'+str(session.step)+'.json'),commit)
            if target==contract['total_steps']:
                session.finish_raw()
                receipt['assessment']=original_summary(session.audit,session.case,session.header,session.step,True)
                receipt['negative_controls']=negatives(session.capture,launch/'negative_controls',index['raw'])
                receipt['status']='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart'
                require_full_reference(receipt,session.case)
            else:
                receipt['status']='new_durable_reference_segment_committed_NOT_full_equation_review'
            receipt.update(accepted_step=session.step,checkpoint_HEAD_sha256=session.expected_head,commits=commits,
                complete_accepted_state=True,actual_review_restart_qualified=False,
                comparison_scope='explicit mathematical_state_sha256 excludes only preserved measured CPU diagnostic;original full state SHA/CPU also retained')
        verify(plan['artifacts']+[index['raw'],request['index'],request['state_contract'],request['resource_contract'],request['plan']])
    except BaseException as error:
        receipt.update(status='new_durable_reference_failed',failure=repr(error),failure_traceback=traceback.format_exc(),
                       last_complete_commit=commits[-1] if commits else None,last_computed_cursor=cursor)
    finally:
        if session is not None:session.close()
    receipt.update(finished_at=stamp(),CPU_seconds=time.process_time())
    atomic_new(launch/'reference_review.json',receipt)
    print(json.dumps(dict(status=receipt['status'],failure=receipt.get('failure'))),flush=True)
    return int('failure' in receipt)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--request',required=True,type=Path);p.add_argument('--request-sha',required=True)
    a=p.parse_args();sys.exit(run(a.request,a.request_sha))

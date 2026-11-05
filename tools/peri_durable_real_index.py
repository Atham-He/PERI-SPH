"""New readonly actual-v23 byte index,zero independent material equations."""
from pathlib import Path
import argparse,json,sys,time


def run(root):
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
    from peri_elastic_normal import read,record,verify,save,stamp
    from peri_durable_reference_stream import build_index
    plan=read(root/'plan.json');verify(plan['artifacts']+[plan['raw'],plan['native_attempt'],plan['original_plan']])
    native=read(plan['native_attempt']['path']);old=read(plan['original_plan']['path'])
    assert native['return_code']==0 and native['stop_reason'] is None and native['stderr']['size']==0
    expected=dict(native['birth_identity'],argv=native['command'])
    index=build_index(plan['raw']['path'],old['case'],plan['raw'],expected)
    assert index['trial_points']==11847 and len(index['index'])==81
    assert index['index'][-1]['time_hex']==float(.4).hex()
    save(root/'accepted_binary_index.json',index)
    # A separately coded byte pass checks every offset,prefix,line,complete
    # ordered counter and Newton count. This is not an equation reassessment.
    import hashlib
    prefix=hashlib.sha256();points=0;count=0
    with Path(plan['raw']['path']).open('rb') as stream:
        for line in stream:
            prefix.update(line);row=json.loads(line)
            if row['type'] not in ('header','step'):continue
            step=0 if row['type']=='header' else row['state']['step']
            if step:
                assert step==count+1;count=step;points+=33*len(row['solve']['trace'])
            cursor=index['index'][step]['raw_cursor']
            assert stream.tell()==cursor['offset'] and prefix.hexdigest()==cursor['prefix_sha256']
            assert hashlib.sha256(line).hexdigest()==cursor['last_line_sha256']
            assert index['index'][step]['trial_points']==points
    assert count==80 and points==11847 and prefix.hexdigest()==plan['raw']['sha256']
    verify(plan['artifacts']+[plan['raw'],plan['native_attempt'],plan['original_plan']])
    result=dict(status='new_actual80_binary_frontier_index_qualified_NOT_equation_review',finished_at=stamp(),
        plan=record(root/'plan.json'),index=record(root/'accepted_binary_index.json'),raw=plan['raw'],
        native_attempt=plan['native_attempt'],accepted_steps=count,material_point_metadata=points,
        final_time=.4,independent_all81offset_prefix_lineSHA_and_pointcounter_passed=True,
        precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
        physical_solver_calls=0,material_calculations=0,actual_material_equation_reviews=0,
        actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,
        S4_complete=False,calibration_allowed=False,CPU_seconds=time.process_time())
    save(root/'independent_review.json',result);print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();run(a.root)

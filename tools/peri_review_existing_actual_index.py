"""NEW readonly independent review of an existing producer artifact."""
from pathlib import Path
import faulthandler,json,os,resource,sys,time


def main(root):
    fault=(root/'faulthandler.txt').open('x');faulthandler.enable(file=fault,all_threads=True)
    def stage(name,step=-1):
        with (root/'stage_witness.jsonl').open('a') as f:
            f.write(json.dumps(dict(stage=name,step=step,CPU_seconds=time.process_time(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))+'\n');f.flush();os.fsync(f.fileno())
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
    from peri_elastic_normal import read,save,stamp
    from peri_buffered_metadata import verify_unique,buffered_record
    from peri_existing_index_assessment import assess
    plan=read(root/'plan.json');pins=plan['artifacts']+plan['input_records']
    stage('before_initial_unique_pin_verification');initial=verify_unique(pins);stage('after_initial_unique_pin_verification')
    native=read(plan['native_attempt']['path']);producer=read(plan['producer_attempt']['path']);original=read(plan['original_plan']['path'])
    assert native['return_code']==0 and native['stop_reason'] is None and native['stderr']['size']==0
    assert producer['return_code']==-9 and 'wall budget' in producer['stop_reason']
    index=read(plan['index_artifact']['path']);frontiers=[read(p['path']) for p in plan['frontiers']]
    assessment=assess(index,plan['raw']['path'],original['case'],dict(native['birth_identity'],argv=native['command']),plan['raw'],frontiers,stage)
    assert assessment['accepted_steps']==80 and assessment['frontiers']==81 and assessment['point_metadata']==11847 and assessment['final_time']==.4
    stage('before_final_unique_pin_verification');final=verify_unique(pins);stage('after_final_unique_pin_verification')
    assert initial==final
    q=dict(status='new_existing_actual80_index_artifact_independently_qualified_producer_failure_preserved_NOT_equation_review',finished_at=stamp(),
        plan=buffered_record(str(root/'plan.json')),raw=plan['raw'],index_artifact=plan['index_artifact'],native_attempt=plan['native_attempt'],
        producer_attempt=plan['producer_attempt'],producer_return_code_preserved=-9,producer_status_changed=False,
        assessment=assessment,unique_pin_verification=initial,stage_witness=buffered_record(str(root/'stage_witness.jsonl')),
        precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
        CPU_seconds=time.process_time(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        index_builder_called=False,physical_solver_calls=0,material_calculations=0,old_math_reviews_repeated=False,
        old_SIGSEGV_root_cause_resolved=False,actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    save(root/'independent_review.json',q);fault.close();print(json.dumps(q))


if __name__=='__main__':main(Path(sys.argv[1]))

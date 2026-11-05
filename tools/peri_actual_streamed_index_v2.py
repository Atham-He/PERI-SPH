"""New actual raw indexing admission using qualified incremental typed bytes.

No solver/material equations or old state recovery. Old SIGSEGV root cause
remains unresolved even if this distinct implementation completes.
"""
from pathlib import Path
import faulthandler,hashlib,json,os,resource,sys,time


def main(root):
    fault=(root/'faulthandler.txt').open('x');faulthandler.enable(file=fault,all_threads=True)
    def stage(name,step):
        item=dict(stage=name,step=step,CPU_seconds=time.process_time(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        with (root/'stage_witness.jsonl').open('a') as stream:stream.write(json.dumps(item)+'\n');stream.flush();os.fsync(stream.fileno())
    stage('before_fixed_backend',0)
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
    from peri_elastic_normal import read,record,verify,save,stamp
    from peri_streamed_raw_index import build_streamed_index
    plan=read(root/'plan.json');verify(plan['artifacts']+[plan['raw'],plan['native_attempt'],plan['original_plan']])
    native=read(plan['native_attempt']['path']);old=read(plan['original_plan']['path'])
    assert native['return_code']==0 and native['stop_reason'] is None and native['stderr']['size']==0
    expected=dict(native['birth_identity'],argv=native['command']);frontiers=root/'frontiers';frontiers.mkdir()
    def accepted(step,value):
        target=frontiers/(str(step)+'.json');assert not target.exists();save(target,value)
        stage('complete_index_frontier_durable',step)
    stage('before_streamed_index',0)
    index=build_streamed_index(plan['raw']['path'],old['case'],plan['raw'],expected,stage,accepted)
    assert index['trial_points']==11847 and len(index['index'])==81 and index['index'][-1]['time_hex']==float(.4).hex()
    save(root/'accepted_binary_index.json',index);stage('complete_index_saved',80)
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
            assert read(frontiers/(str(step)+'.json'))==index['index'][step]
    assert count==80 and points==11847 and prefix.hexdigest()==plan['raw']['sha256'];stage('independent_metadata_pass_complete',80)
    verify(plan['artifacts']+[plan['raw'],plan['native_attempt'],plan['original_plan']])
    q=dict(status='new_actual80_streamed_binary_index_qualified_NOT_equation_review',finished_at=stamp(),plan=record(root/'plan.json'),
        index=record(root/'accepted_binary_index.json'),raw=plan['raw'],native_attempt=plan['native_attempt'],
        accepted_steps=count,material_point_metadata=points,final_time=.4,durable_frontiers=81,
        independent_all81offset_prefix_lineSHA_and_pointcounter_passed=True,
        stage_witness=record(root/'stage_witness.jsonl'),faulthandler=record(root/'faulthandler.txt'),
        CPU_seconds=time.process_time(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
        originalfailedindex_retried=False,old_SIGSEGV_root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,
        actual_material_equation_reviews=0,actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    save(root/'independent_review.json',q);fault.close();print(json.dumps(q))


if __name__=='__main__':main(Path(sys.argv[1]))

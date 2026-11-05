"""One prospectively frozen fixed-backend reference process per NEW CP raw.

No native dispatch, old raw review rerun or exception-triggered backend choice.
Identical complete trial packets import authenticated continuous evidence;
different NEW packets use the qualified selected checkpoint as initial data.
"""
import argparse,json,os,resource,sys,time,traceback
from pathlib import Path


def configure_backend(reference):
    # Packaged script and _pydecimal reside in frozen/, before any reference
    # or Fraction import. Both actual paths are frozen by the future controller.
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(reference)
    backend.getcontext().prec=95
    assert backend.getcontext().prec==95
    return backend


def selftest_binding(reference):
    backend=configure_backend(reference)
    from fractions import Fraction
    assert Fraction(backend.Decimal('.1'))==Fraction(1,10)
    with backend.localcontext() as context:assert context.prec==95
    import peri_physical_time_friction_audit as physical
    import peri_physical_friction_negative as negatives
    assert physical.D is backend.Decimal and negatives.D is backend.Decimal
    print(json.dumps(dict(status='new_physical_time_CP_reference_binding_and_95digit_context_preflight_passed',
        backend=backend.__file__,precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
        native_solver_calls=0,material_calculations=0,old_raw_reviews_repeated=0)))


def require_full_reference(review,case):
    import math
    assert review['status']=='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart'
    assert review['backend_selected_before_reference_imports'] is True and review['exception_triggered_backend_switch'] is False
    assert review['native_solver_calls']==0
    q=review['assessment'];assert q['name']==case['name'] and q['accepted_steps']==case['steps']
    assert q['all_trial_point_tangents_checked'] is True and q['independent_trial_material_points']>0
    assert math.isfinite(q['maximum_normalized_error']) and 0<=q['maximum_normalized_error']<=1e-10
    assert q['geometry']['nodes']==33 and q['geometry']['physical_DOF']==99 and q['geometry']['normal_rows']==116
    expected={'wrong_physical_mass','altered_tangent_dual_history','kinetic_using_wrong_mass','swapped_material_algorithmic_losses',
              'wrong_friction_loss','wrong_wall_velocity_cursor','corrupt_compensation_history','endpoint_substituted_for_whole_interval',
              'wrong_physical_tick_cursor','wrong_protocol_denominator'}
    assert len(review['negative_controls'])==10 and {c['name'] for c in review['negative_controls']}==expected
    assert all(c['rejected'] is True for c in review['negative_controls'])


def run(root,label):
    backend=configure_backend(root/'frozen/reference')
    from peri_elastic_normal import read,record,verify,save,identity,stamp
    from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector as FrictionInspector,inspect_physical_stream
    from peri_physical_friction_negative import OneStepCapture,negatives
    from peri_physical_checkpoint_support import derive_continuous,trace_index,trace_packets,compare_new_branch,loads
    from peri_coulomb_contact_audit import dec,frac
    plan=read(root/'plan.json');attempt=read(root/(label+'.attempt.json'));budget=read(root/(label+'.reference_admission.json'))
    assert not (root/(label+'.reference_state.json')).exists(),'never repeat an attempted reference state'
    plan_record=read(root/'state.json')['plan']
    verify([plan_record]+plan['artifacts']);verify([attempt['stdout'],attempt['stderr'],budget['attempt']])
    assert attempt['return_code']==0 and attempt['stop_reason'] is None
    save(root/(label+'.reference_state.json'),dict(status='new_physical_time_CP_reference_attempted',identity=identity(),started_at=stamp(),admission=record(root/(label+'.reference_admission.json'))))
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30));resource.setrlimit(resource.RLIMIT_CPU,(budget['CPU_seconds'],budget['CPU_seconds']))
    resource.setrlimit(resource.RLIMIT_FSIZE,(256*2**20,256*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    start=time.monotonic();case=plan['case'];cursor={};old_step=FrictionInspector.step
    receipt=dict(status='physical_time_CP_new_reference_failed',recorded_at=stamp(),identity=identity(),backend=backend.__file__,
        backend_selected_before_reference_imports=True,exception_triggered_backend_switch=False,native_solver_calls=0,
        full_restart_qualified=False,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False)
    def observed(self,row,old,z,loss,history):
        cursor.update(old_step=old['step'],recorded_step=row['state']['step'],time=row['state']['time'])
        result=old_step(self,row,old,z,loss,history)
        with (root/(label+'.reference_progress.jsonl')).open('a') as stream:
            stream.write(json.dumps(dict(recorded_at=stamp(),cursor=cursor,trial_points=self.trial_points,
                maximum_normalized_error=self.maximum,CPU_seconds=time.process_time()),allow_nan=False)+'\n')
            stream.flush();os.fsync(stream.fileno())
        return result
    FrictionInspector.step=observed
    try:
        with Path(attempt['stdout']['path']).open() as stream:meta=json.loads(next(stream))
        assert all(meta[k]==v for k,v in attempt['birth_identity'].items()) and meta['argv']==attempt['command']
        if attempt['child_identity']:assert all(meta[k]==v for k,v in attempt['child_identity'].items())
        if label=='continuous':
            capture=OneStepCapture();receipt['assessment']=inspect_physical_stream(case,attempt['stdout']['path'],full=True,capture=capture)
            receipt['negative_controls']=negatives(capture,root/'continuous_negatives',attempt['stdout'])
            assert len(receipt['negative_controls'])==10 and all(x['rejected'] for x in receipt['negative_controls'])
            assert receipt['assessment']['accepted_steps']==case['steps'] and receipt['assessment']['all_trial_point_tangents_checked']
            receipt['status']='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart'
            require_full_reference(receipt,case)
        else:
            previous=read(root/'continuous.independent_review.json')
            assert previous['status']=='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart'
            require_full_reference(previous,case)
            assert previous['assessment']['all_trial_point_tangents_checked'] and len(previous['negative_controls'])==10 and all(x['rejected'] for x in previous['negative_controls'])
            verify([previous['attempt'],previous['raw']])
            byteproof=read(root/(label+'.bytes_review.json'))
            assert byteproof['full_candidate_byte_equal'] and byteproof['full_physical_prefix_byte_equal']
            verify([byteproof['native_attempt']]);receipt['byte_assessment']=record(root/(label+'.bytes_review.json'))
            derived=derive_continuous(root,previous['raw']['path'],case,bytes.fromhex(plan['trusted_contract']))
            first,last=plan['segments'][label]
            generation=None if first==1 else first-1
            branch=trace_index(attempt['stdout']['path'],case,first,last,resume_generation=generation)
            comparison=compare_new_branch(derived['index'],branch,case);receipt['complete_trial_comparison']=comparison
            receipt['imported_continuous_review']=record(root/'continuous.independent_review.json')
            if not comparison['different_steps']:
                receipt.update(status='new_physical_time_CP_branch_full_trial_packets_identical_math_evidence_imported',new_reference_trial_points=0,
                    old_full_reference_repeated=False,scope='Entire new accepted-step line incl ALL Newton/trial/pointT bytes matches authenticated continuous;full140field/physicalprefix comparison separately required.')
            else:
                # The native accepted state/prefix has ALREADY matched all bytes.
                # A distinct new trial context is evaluated under original gates.
                audit=FrictionInspector(case,derived['index']['header'],True)
                old,z,loss,history=audit.initial()
                if first>1:
                    checkpoint=derived['checkpoints'][first-1]
                    old=checkpoint['full_native_state'];history=checkpoint['full_native_friction_history']
                    z=[list(map(dec,row)) for row in old['viscous']];loss=list(map(dec,old['dissipation']))
                    # Diagnostic active-history comes from already qualified
                    # prefix packets; no old material equation is re-assessed.
                    with Path(previous['raw']['path']).open('rb') as stream:
                        for line in stream:
                            row=loads(line)
                            if row.get('type')!='step':continue
                            if row['state']['step']>=first:break
                            for k,normal in enumerate(row['state']['last_normal_impulse']):
                                if frac(normal)>frac(audit.I)*frac(1e-8):audit.ever_active.add(k)
                for row,_ in trace_packets(attempt['stdout']['path'],case,first,last,resume_generation=generation):
                    if row['type']=='step':old,z,loss,history=audit.step(row,old,z,loss,history)
                assert old['step']==last
                receipt.update(status='new_physical_time_CP_branch_distinct_trial_context_independently_passed',
                    new_reference_trial_points=audit.trial_points,maximum_normalized_error=max(audit.maximum,audit.geometry['maximum_normalized_cache_error']),
                    all_new_trial_point_tangents_checked=True,reference_seed='Authenticated selected CP exactbinary viscous/dissipation/kinematics/travel+correction+ledger;no native history reset or replay of passed prefix math.',
                    independently_assessed_new_segment=[first,last],old_full_reference_repeated=False)
        verify([plan_record]+plan['artifacts']);verify([attempt['stdout'],attempt['stderr'],budget['attempt']])
        receipt.update(raw=attempt['stdout'],attempt=record(root/(label+'.attempt.json')))
    except Exception as error:receipt.update(failure=repr(error),failure_traceback=traceback.format_exc(),last_cursor=cursor)
    finally:FrictionInspector.step=old_step
    receipt.update(finished_at=stamp(),wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time(),
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,plan=record(root/'plan.json'),
        admission=record(root/(label+'.reference_admission.json')))
    save(root/(label+'.independent_review.json'),receipt)
    print(json.dumps(dict(status=receipt['status'],failure=receipt.get('failure'))),flush=True)
    return int('failure' in receipt)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path);parser.add_argument('--label');parser.add_argument('--selftest-binding',type=Path)
    args=parser.parse_args()
    if args.selftest_binding:selftest_binding(args.selftest_binding)
    else:
        assert args.root and args.label,'explicit root/label required'
        sys.exit(run(args.root,args.label))

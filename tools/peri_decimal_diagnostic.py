"""Instrument a failed reader once; never call or redispatch a solver.

Original Decimal arithmetic, precision, physical gates and raw input stay
unchanged. Extra wrappers observe cursors/operands and an exception only.
Even a clean diagnostic replay does not overwrite the failed study verdict.
"""
import argparse
import ctypes
import decimal
import dis
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import traceback


def context_snapshot(context=None):
    context=decimal.getcontext() if context is None else context
    # Read ONLY this live owned object's advertised fixed allocation. Do not
    # infer/write C fields or traverse arbitrary addresses. The raw observation
    # supplements the public context and has no mechanical/checkpoint meaning.
    raw=ctypes.string_at(id(context),type(context).__basicsize__)
    return dict(precision=context.prec,Emin=context.Emin,Emax=context.Emax,
        rounding=context.rounding,clamp=context.clamp,
        flags={k.__name__:v for k,v in context.flags.items()},
        traps={k.__name__:v for k,v in context.traps.items()},
        context_type_basicsize=type(context).__basicsize__,owned_context_object_hex=raw.hex())


def simple(value):
    if isinstance(value,decimal.Decimal):return dict(decimal=str(value))
    if isinstance(value,(int,float,str,bool)) or value is None:return value
    if isinstance(value,(list,tuple)) and len(value)<=128:return [simple(x) for x in value]
    return dict(type=type(value).__name__)


def exception_operands(error):
    selected=[];tb=error.__traceback__
    wanted={'inverse_transpose':['f','j','a','b','c','d','e','z','g','h','i'],
        'pk1':['f','x','p','j'],
        'average':['f0','f1','x0','dt','p','h','s','w','decay','df','x'],
        'material_step':['f0','f1','z0','loss0','dt','p','epsilon','column','plus','minus','min_j'],
        'evaluate':['i','v','dt','p','f0','f1']}
    while tb is not None:
        frame=tb.tb_frame;name=frame.f_code.co_name
        if name in wanted:
            instructions=[i for i in dis.get_instructions(frame.f_code) if i.offset==tb.tb_lasti]
            instruction=instructions[0] if instructions else None
            selected.append(dict(function=name,file=frame.f_code.co_filename,line=tb.tb_lineno,
                bytecode_offset=tb.tb_lasti,opcode=instruction.opname if instruction else None,
                operation=instruction.argrepr if instruction else None,
                source_columns=[instruction.positions.col_offset,instruction.positions.end_col_offset] if instruction else None,
                locals={k:simple(frame.f_locals[k]) for k in wanted[name] if k in frame.f_locals}))
        tb=tb.tb_next
    return selected


def failed_context(error):
    tb=error.__traceback__
    while tb is not None:
        context=tb.tb_frame.f_locals.get('ctx')
        if isinstance(context,decimal.Context):return context_snapshot(context)
        tb=tb.tb_next
    return context_snapshot()


def run(root):
    assert not (root/'state.json').exists(),'never repeat an attempted diagnostic'
    plan=json.loads((root/'plan.json').read_text())
    original=Path(plan['original_root']);sys.path.insert(0,str(original/'frozen'))
    from peri_elastic_normal import identity,read,record,verify,save,stamp,reserves
    from peri_bounded_friction_audit import FrictionInspector,inspect_stream
    from peri_bounded_friction_negative import OneStepCapture
    import peri_gel_affine_audit as affine
    import _decimal
    verify(plan['artifacts'])
    assert record(Path(__file__).resolve())==plan['diagnostic_source']
    assert read(original/'terminal.json')['status']=='stopped_for_review'
    prior=read(original/(plan['case']+'.independent_review.json'))
    assert prior['failure']=="RuntimeError('internal error in flags_as_exception')"
    save(root/'state.json',dict(status='attempted_read_only_runtime_diagnostic',identity=identity(),started_at=stamp()))
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30))
    resource.setrlimit(resource.RLIMIT_CPU,(14400,14400))
    resource.setrlimit(resource.RLIMIT_FSIZE,(64*2**20,64*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    reserves(root)
    start=time.monotonic();cursor={};packet_written=False
    result=dict(status='instrumented_reader_diagnostic_failed',recorded_at=stamp(),identity=identity(),
        solver_calls=0,original_failure_preserved=True,study_passed=False,S3_complete=False,
        full_restart_qualified=False,calibration_allowed=False,arithmetic_replaced=False,
        decimal_backend=dict(module=_decimal.__file__,version=_decimal.__version__,libmpdec=_decimal.__libmpdec_version__))
    original_inverse=affine.inverse_transpose
    original_step=FrictionInspector.step
    def observed_inverse(f):
        nonlocal packet_written
        try:return original_inverse(f)
        except Exception as error:
            if not packet_written:
                packet_written=True
                save(root/'exception_packet.json',dict(cursor=cursor.copy(),failure=repr(error),
                    context=context_snapshot(),operands=exception_operands(error),input_f=simple(f),
                    resource_usage=dict(max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime),
                    statement='Observation only: no retried operation, trap clearing, rounding or precision change.'))
            raise
    def observed_step(self,row,old,z,loss,history):
        cursor.update(step_before=old['step'],recorded_step=row['state']['step'],
                      old_time=old['time'],new_time=row['state']['time'],phase=row['interval']['phase'])
        answer=original_step(self,row,old,z,loss,history)
        event=dict(recorded_at=stamp(),cursor=cursor.copy(),checked_trial_points=self.trial_points,
            maximum_normalized_error=self.maximum,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        with (root/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(event,allow_nan=False)+'\n');stream.flush()
            import os
            os.fsync(stream.fileno())
        print(json.dumps(dict(checked_steps=row['state']['step'],points=self.trial_points)),flush=True)
        return answer
    affine.inverse_transpose=observed_inverse
    FrictionInspector.step=observed_step
    try:
        prepare=read(original/'preparation.json');verify(prepare['artifacts']+prepare['inputs'])
        case=next(c for c in read(original/'plan.json')['controls'] if c['name']==plan['case'])
        attempt=read(original/(plan['case']+'.attempt.json'))
        verify([attempt['stdout'],attempt['stderr']]);assert attempt['return_code']==0
        result['assessment']=inspect_stream(case,attempt['stdout']['path'],0,True,OneStepCapture())
        # This is deliberately not a replacement positive case review. The
        # original unexplained failure and missing negative-control qualification
        # remain until separately adjudicated from explicit evidence.
        result['status']='instrumented_original_reader_completed_without_reproducing_decimal_exception_NOT_STUDY_PASS'
        verify(plan['artifacts']);verify(prepare['artifacts']+prepare['inputs'])
    except Exception as error:
        result.update(failure=repr(error),failure_traceback=traceback.format_exc(),
            last_cursor=cursor.copy(),decimal_context_at_failure=failed_context(error),
            exception_operands=exception_operands(error))
        if isinstance(error,RuntimeError) and str(error)=='internal error in flags_as_exception':
            result['status']='original_decimal_exception_reproduced_operands_preserved_NOT_STUDY_PASS'
    finally:
        affine.inverse_transpose=original_inverse;FrictionInspector.step=original_step
    result.update(finished_at=stamp(),wall_seconds=time.monotonic()-start,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        exception_packet=record(root/'exception_packet.json') if packet_written else None,
        progress=record(root/'progress.jsonl') if (root/'progress.jsonl').exists() else None,
        plan=record(root/'plan.json'),state=record(root/'state.json'))
    save(root/'diagnostic_receipt.json',result)
    print(json.dumps({k:result.get(k) for k in ['status','failure','wall_seconds','peak_RSS_KiB']}),flush=True)
    return int('failure' in result and result['status']!='original_decimal_exception_reproduced_operands_preserved_NOT_STUDY_PASS')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    sys.exit(run(parser.parse_args().root))

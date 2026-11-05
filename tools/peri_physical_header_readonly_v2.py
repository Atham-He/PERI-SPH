"""Fresh fixed-pure95 actual v2 header/initial-state admission, ZERO PD steps."""
import argparse,json,os,sys,time,traceback
from pathlib import Path

def verify_one(receipt):
    from peri_elastic_normal import verify
    verify([receipt])

def main(root):
    from peri_fixed_python_raw_audit import select_backend
    native_root=Path(json.loads((root/'plan.json').read_text())['native_root'])
    backend=select_backend(root/'frozen/reference')
    backend.getcontext().prec=95
    from fractions import Fraction
    assert Fraction(backend.Decimal('.1'))==Fraction(1,10)
    from peri_elastic_normal import read,save,verify,record,identity,stamp
    from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector
    from peri_canonical_archive_audit import parse
    from peri_bounded_mechanical_stream_audit import loads
    from peri_physical_checkpoint_support import header_body_fields
    manifest=read(native_root/'header_reference_admission.json');verify(manifest['artifacts'])
    assert not (root/'header_reference_state.json').exists(),'never repeat this attempted reference'
    save(root/'header_reference_state.json',dict(started_at=stamp(),identity=identity(),admission=record(root/'header_reference_admission.json')))
    start=time.monotonic();cpu=time.process_time()
    result=dict(S4_complete=False,calibration_allowed=False,new_PD_steps=0,native_solver_calls=0,
                all_trial_point_tangents_checked=False,scope='actual header/initial geometry and algebra ONLY,not a trajectory or nonzero CP',
                backend=backend.__file__,backend_selected_before_reference_imports=True,exception_triggered_backend_switch=False)
    try:
        raw=Path(manifest['raw']['path']);verify_one(manifest['raw'])
        rows=[loads(line) for line in raw.read_bytes().splitlines()]
        assert len(rows)==3 and rows[0]['type']=='identity' and rows[1]['type']=='header'
        assert rows[2]==dict(type='initial_header_only',physical_solver_calls=0,new_physical_steps=0,owned_state_installed=False)
        attempt=read(native_root/'describe_initial.attempt.json')
        assert all(rows[0][k]==v for k,v in attempt['birth_identity'].items()) and rows[0]['argv']==attempt['command']
        case=manifest['case'];header=rows[1]
        audit=PhysicalTimeFrictionInspector(case,header,True);old,z,loss,history=audit.initial()
        assert old['step']==0 and old['time']==0 and audit.trial_points==0
        _,fields=parse((native_root/'initial_candidate.bin').read_bytes())
        body={k.removeprefix('body/'):v for k,v in fields.items() if k.startswith('body/')}
        assert header_body_fields(header)==body and len(body)==98
        negatives=[]
        changes=[('old_format','format',1),('wrong_denominator','time_denominator',2000000),
                 ('wrong_nominal_ticks','nominal_ticks',2500),('changed_physical_end','phase_end_times',[.08,.12,.2]),
                 ('wrong_total_steps','steps_requested',79),('wrong_checkpoint_selection','checkpoint_steps',[16,80])]
        for name,key,value in changes:
            altered=dict(header);altered[key]=value
            try:audit.validate_protocol_header(case,altered)
            except AssertionError:negatives.append(dict(name=name,rejected=True))
            else:raise AssertionError(('new runtime header corruption accepted',name))
        verify(manifest['artifacts'])
        result.update(status='existing_actual_v2_header_independently_qualified_by_new_readonly_reader_ZERO_steps',
                      maximum_normalized_error=max(audit.maximum,audit.geometry['maximum_normalized_cache_error']),
                      geometry=audit.geometry,accepted_initial_body_all_bytes_match=True,
                      planned_steps=header['steps_requested'],fixed_end_time=header['end_time'],phase_step_counts=header['phase_step_counts'],
                      clock_header_negative_controls=negatives,independent_trial_material_points=0,
                      reference=record(native_root/'header_reference_admission.json'),raw=manifest['raw'])
    except Exception as e:result.update(status='new_v2_runtime_header_admission_failed_preserve_evidence',failure=repr(e),traceback=traceback.format_exc())
    result.update(finished_at=stamp(),wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,identity=identity())
    save(root/'header_independent_review.json',result)
    print(json.dumps(dict(status=result['status'],failure=result.get('failure'),new_PD_steps=0)))
    return int('failure' in result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();sys.exit(main(a.root))

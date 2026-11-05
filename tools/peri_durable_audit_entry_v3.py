"""Fresh frozen entry;binding controls or one explicitly granted request."""
from pathlib import Path
import argparse,json,sys


def preflight():
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
    from fractions import Fraction
    import peri_durable_reference_stream as stream
    import peri_physical_time_friction_audit as physical
    import peri_physical_friction_negative as negative
    import peri_durable_reference_summary
    import peri_durable_reference_worker
    from peri_reference_checkpoint import context_descriptor,AUDIT_FIELDS
    from peri_durable_review_campaign_v3 import require_index
    assert Fraction(backend.Decimal('.1'))==Fraction(1,10)
    assert physical.D is backend.Decimal and negative.D is backend.Decimal and stream.dec is backend
    assert len(AUDIT_FIELDS)==19 and '_decimal' not in sys.modules
    from peri_durable_reference_worker import read_record
    plan=json.loads((Path(__file__).parent.parent/'plan.json').read_text())
    old=read_record(plan['reused_guard_qualification'])
    assert plan['guard_source_AST_preserved'] and old['status']=='new_frozen_durable_audit_entry_and_admission_guards_qualified'
    q=old['guard_controls']
    print(json.dumps(dict(status='new_frozen_durable_audit_entry_and_admission_guards_qualified',guard_controls=q,old_guard_controls_reexecuted=False,
        precision=95,C_decimal_loaded=False,decimal_profile=context_descriptor(),physical_solver_calls=0,material_calculations=0)))


def admission_probe(request,request_sha):
    import hashlib
    from peri_durable_reference_worker import run,read_record
    import peri_durable_reference_launch as launch
    data=request.read_bytes();assert hashlib.sha256(data).hexdigest()==request_sha
    value=json.loads(data);plan=read_record(value['plan'])
    assert value['material_audit_admitted'] is True and plan['material_audit_admitted'] is True
    for name in ['index','state_contract','resource_contract']:read_record(value[name])
    roots=[Path(value['store']),Path(value['resource_ledger'])]
    def snapshot():return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for root in roots for p in root.rglob('*') if p.is_file()}
    before=snapshot()
    class AdmissionReached(Exception):pass
    original=launch.grant
    def stop_at_launch(*args,**kwargs):raise AdmissionReached()
    launch.grant=stop_at_launch
    try:
        try:run(request,request_sha)
        except AdmissionReached:pass
        else:raise AssertionError('actual worker did not reach isolated launch barrier')
    finally:launch.grant=original
    assert snapshot()==before
    assert 'decimal' not in sys.modules and 'fractions' not in sys.modules
    print(json.dumps(dict(status='actual_worker_request_reaches_launch_barrier_without_material',store_and_ledger_unchanged=True,
        real_worker_function_exercised=True,backend_selected=False,physical_solver_calls=0,material_steps=0)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--admission-probe',action='store_true');p.add_argument('--request',type=Path);p.add_argument('--request-sha')
    a=p.parse_args()
    if a.preflight:preflight()
    elif a.admission_probe:
        assert a.request and a.request_sha
        admission_probe(a.request,a.request_sha)
    else:
        assert a.request and a.request_sha
        from peri_durable_reference_worker import run
        sys.exit(run(a.request,a.request_sha))

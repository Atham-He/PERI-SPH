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
    from peri_durable_review_campaign import guard_controls
    assert Fraction(backend.Decimal('.1'))==Fraction(1,10)
    assert physical.D is backend.Decimal and negative.D is backend.Decimal and stream.dec is backend
    assert len(AUDIT_FIELDS)==19 and '_decimal' not in sys.modules
    q=guard_controls()
    print(json.dumps(dict(status='new_frozen_durable_audit_entry_and_admission_guards_qualified',guard_controls=q,
        precision=95,C_decimal_loaded=False,decimal_profile=context_descriptor(),physical_solver_calls=0,material_calculations=0)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--request',type=Path);p.add_argument('--request-sha')
    a=p.parse_args()
    if a.preflight:preflight()
    else:
        assert a.request and a.request_sha
        from peri_durable_reference_worker import run
        sys.exit(run(a.request,a.request_sha))

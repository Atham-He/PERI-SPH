"""Prospective new structural controls in a fixed pure95 fresh process."""
from pathlib import Path
import io,json,os,sys,unittest
frozen=Path(__file__).parent
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(frozen/'reference');backend.getcontext().prec=95
os.environ['DURABLE_COORDINATOR_WORK']=str(frozen.parent/'work');Path(os.environ['DURABLE_COORDINATOR_WORK']).mkdir()
import test_durable_reference_coordinator
stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test_durable_reference_coordinator))
print(json.dumps(dict(status='new_durable_reference_stream_and_coordinator_structural_controls_passed' if result.wasSuccessful() else 'new_durable_coordinator_controls_failed',
    tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),log=stream.getvalue(),precision=backend.getcontext().prec,
    C_decimal_loaded='_decimal' in sys.modules,physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
    actual_full_material_audit_admitted=False,actual_review_restart_qualified=False,S4_complete=False,calibration_allowed=False),indent=2))
sys.exit(0 if result.wasSuccessful() else 1)

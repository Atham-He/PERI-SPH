"""Fresh fixed-pure95 process for declared synthetic audit state controls."""
from pathlib import Path
import io,json,os,sys,unittest
frozen=Path(__file__).resolve().parent;sys.path.insert(0,str(frozen))
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(frozen/'reference');backend.getcontext().prec=95
os.environ['REFERENCE_CP_CONTROL_WORK']=str(frozen.parent/'work')
import test_reference_checkpoint
stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test_reference_checkpoint))
receipt=dict(status='new_lossless_pure95_reference_state_transaction_controls_passed' if result.wasSuccessful() else 'new_reference_state_transaction_controls_failed',
    tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),log=stream.getvalue(),backend=backend.__file__,precision=backend.getcontext().prec,
    C_decimal_loaded='_decimal' in sys.modules,physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
    actual_inspector_integration_qualified=False,cumulative_budget_ledger_qualified=False,actual_review_restart_qualified=False,
    hardware_power_loss_qualified=False,S4_complete=False,calibration_allowed=False)
print(json.dumps(receipt,indent=2));sys.exit(0 if result.wasSuccessful() else 1)

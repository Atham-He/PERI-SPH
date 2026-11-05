"""Fresh pure95 process;new ledger plus actual-class memory controls only."""
from pathlib import Path
import io
import json
import os
import sys
import unittest

frozen=Path(__file__).resolve().parent
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(frozen/'reference');backend.getcontext().prec=95
os.environ['REFERENCE_BUDGET_WORK']=str(frozen.parent/'work')
Path(os.environ['REFERENCE_BUDGET_WORK']).mkdir()
import test_reference_budget_owner
stream=io.StringIO()
result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test_reference_budget_owner))
print(json.dumps(dict(status='new_cumulative_resource_ledger_and_actual_class_memory_controls_passed' if result.wasSuccessful() else 'new_budget_owner_controls_failed',
    tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),log=stream.getvalue(),
    precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
    physical_solver_calls=0,material_calculations=0,real_trajectory_raw_reads=0,
    actual_computed_trajectory_owner_qualified=False,actual_review_restart_qualified=False,
    hardware_power_loss_qualified=False,S4_complete=False,calibration_allowed=False),indent=2))
sys.exit(0 if result.wasSuccessful() else 1)

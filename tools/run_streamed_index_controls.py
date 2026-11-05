from pathlib import Path
import io,json,sys,unittest
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
import test_streamed_raw_index
log=io.StringIO();r=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test_streamed_raw_index))
print(json.dumps(dict(status='new_streamed_index_full_bytes_equivalent_on_declared_synthetic_grids' if r.wasSuccessful() else 'new_streamed_index_controls_failed',
    tests_run=r.testsRun,failures=len(r.failures),errors=len(r.errors),log=log.getvalue(),precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
    root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0),indent=2))
sys.exit(0 if r.wasSuccessful() else 1)

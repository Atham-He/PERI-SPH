from pathlib import Path
import io,json,sys,unittest
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
import test_buffered_metadata,test_existing_index_assessment
suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromModule(m) for m in [test_buffered_metadata,test_existing_index_assessment]])
log=io.StringIO();r=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
print(json.dumps(dict(status='new_existing_index_buffered_independent_assessment_controls_passed' if r.wasSuccessful() else 'new_existing_index_controls_failed',
    tests_run=r.testsRun,errors=len(r.errors),failures=len(r.failures),log=log.getvalue(),precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules,
    physical_solver_calls=0,material_calculations=0,old_raw_reads=0),indent=2))
sys.exit(0 if r.wasSuccessful() else 1)

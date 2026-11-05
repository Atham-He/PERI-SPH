from pathlib import Path
import json,sys
from peri_fixed_python_raw_audit import select_backend
backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
from test_streaming_typed_json import run
q=run();q.update(precision=backend.getcontext().prec,C_decimal_loaded='_decimal' in sys.modules)
print(json.dumps(q,indent=2));sys.exit(0 if not q['errors'] and not q['failures'] else 1)

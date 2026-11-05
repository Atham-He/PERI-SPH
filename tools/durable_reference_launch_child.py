"""Only launch handshake control actions;no reference or solver imports."""
import json
from pathlib import Path
import sys
from peri_durable_reference_launch import grant,atomic_new

root=Path(sys.argv[1]);mode=sys.argv[2]
if mode=='early_exit':sys.exit(7)
contract=json.loads(Path(sys.argv[4]).read_bytes())
grant(root,sys.argv[3],contract,sys.argv[5])
atomic_new(root/'controlled_work.json',dict(only_after_committed_grant=True,material_calls=0))
sys.exit(9 if mode=='controlled_failure' else 0)

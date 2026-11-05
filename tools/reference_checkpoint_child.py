"""Declared process-interruption control;not a physical/native executable."""
from pathlib import Path
import json,os,resource,signal,sys
frozen=Path(__file__).resolve().parent
sys.path.insert(0,str(frozen));from peri_fixed_python_raw_audit import select_backend
backend=select_backend(frozen/'reference');backend.getcontext().prec=95
from peri_reference_checkpoint import Store,decode,exclusive_file
root=Path(sys.argv[1]);stage=sys.argv[2];expected=sys.argv[3]
contract=json.loads((root/'contract.json').read_text());state=decode((root/'candidate.bin').read_bytes(),contract)
resource.setrlimit(resource.RLIMIT_CPU,(10,10));resource.setrlimit(resource.RLIMIT_AS,(512*2**20,512*2**20))
fields=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()
print(json.dumps(dict(pid=os.getpid(),ppid=os.getppid(),start_ticks=fields[19],boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),argv=Path('/proc/self/cmdline').read_bytes().decode().strip(chr(0)).split(chr(0)))),flush=True)
def hook(point):
    if point==stage:
        exclusive_file(root/(point+'.witness'),point.encode());os.kill(os.getpid(),signal.SIGKILL)
Store(root/'store',contract).publish(state,'child_'+stage,expected,hook)
raise AssertionError('declared interruption not reached')

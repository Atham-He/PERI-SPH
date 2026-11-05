"""Bounded resource-transaction child controls;never a physical solver."""
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
from peri_reference_budget import Ledger, canonical, write_new

root = Path(sys.argv[1]); mode = sys.argv[2]
resource.setrlimit(resource.RLIMIT_CPU,(5,5)); resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
fields = Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()
identity = dict(pid=os.getpid(),ppid=os.getppid(),start_ticks=int(fields[19]),
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
    argv=Path('/proc/self/cmdline').read_bytes().decode().rstrip('\x00').split('\x00'))
if mode == 'usage':
    temp=root/'usage.identity.pending';write_new(temp,canonical(identity));os.replace(temp,root/'usage.identity.json')
    deadline=time.monotonic()+10
    while not (root/'usage.release').exists():
        if time.monotonic()>deadline:raise RuntimeError('bounded usage child release timeout')
        time.sleep(.005)
    sum(i*i for i in range(10000))
else:
    contract=json.loads((root/'contract.json').read_text())
    def hook(stage):
        if stage==mode:
            write_new(root/(stage+'.witness'),canonical(dict(stage=stage,identity=identity)))
            os.kill(os.getpid(),signal.SIGKILL)
    Ledger(root/'ledger',contract).append(dict(kind='authorize',id='attempt_1'),sys.argv[3],hook)
    raise AssertionError('declared interruption not reached')

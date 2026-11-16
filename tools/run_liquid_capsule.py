#!/usr/bin/env python3
"""Run one fixed-input liquid capsule history without modifying source files."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, default=ROOT/'build/examples/liquid_capsule_240um/peri_sph_liquid_capsule_240um')
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--threads', type=int, default=min(12, os.cpu_count() or 1))
parser.add_argument('--memory-gib', type=float, default=32.)
parser.add_argument('--max-steps', type=int, default=0)
parser.add_argument('--resume', action='store_true')
args = parser.parse_args()
if args.threads < 1 or args.memory_gib <= 0 or args.max_steps < 0:
    parser.error('threads and memory must be positive; max-steps must be nonnegative')
binary = args.binary.resolve()
if not binary.is_file():
    parser.error('Build the executable first: '+str(binary))
out = args.output.resolve()
checkpoint = out/'checkpoint.bin'
if args.resume and not checkpoint.is_file():
    parser.error('--resume requires an existing complete checkpoint')
if not args.resume and out.exists() and any(out.iterdir()):
    parser.error('Output directory is not empty; choose a new directory or reviewed --resume')
out.mkdir(parents=True, exist_ok=True)
binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
identity = out/'run_identity.json'
if args.resume:
    if not identity.is_file():
        parser.error('Missing run identity; raw checkpoint migration requires explicit review')
    if json.loads(identity.read_text())['binary_sha256'] != binary_hash:
        parser.error('Executable changed; checkpoint compatibility must be reviewed before continuation')
else:
    identity.write_text(json.dumps(dict(binary_sha256=binary_hash,
        source_provenance=json.loads((ROOT/'docs/source_provenance.json').read_text())), indent=2)+'\n')
runtime = out/'runtime'
runtime.mkdir(exist_ok=True)
cmd = [str(binary), '--output', str(out), '--threads', str(args.threads)]
if args.resume:
    cmd.append('--resume')
if args.max_steps:
    cmd += ['--max-steps', str(args.max_steps)]
def budget():
    limit = int(args.memory_gib*1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
env = dict(os.environ, MALLOC_ARENA_MAX='2', MALLOC_TRIM_THRESHOLD_='131072')
receipt = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    command=cmd, binary_sha256=binary_hash, physical_inputs_changed=False,
    memory_budget_gib=args.memory_gib)
with (out/'run.log').open('a' if args.resume else 'w') as log:
    process = subprocess.Popen(cmd, cwd=runtime, env=env, preexec_fn=budget,
                               stdout=log, stderr=subprocess.STDOUT)
    receipt['native_pid'] = process.pid
    (out/'launch.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('Running PID', process.pid, '| diagnostics:', out/'diagnostics.csv', flush=True)
    try:
        code = process.wait()
    except KeyboardInterrupt:
        process.send_signal(signal.SIGTERM)
        code = process.wait()
receipt.update(returncode=code, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(out/'launch.json').write_text(json.dumps(receipt, indent=2)+'\n')
raise SystemExit(code)

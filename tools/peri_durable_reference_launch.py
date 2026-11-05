"""One-shot reserved launch/ready/birth/release/wait4 coordinator.

No numerical state retry,implicit CP selection or resource refund on unknown
exit. Caller must authenticate the frozen plan/command/limits and inputs.
"""
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import time
import json
from peri_reference_budget import Ledger,canonical,digest,write_new,sync,identity


def current_identity(pid=None):
    pid=os.getpid() if pid is None else pid;p=Path('/proc')/str(pid)
    fields=(p/'stat').read_text().rsplit(')',1)[1].split()
    return dict(pid=pid,ppid=int(fields[1]),start_ticks=int(fields[19]),
        boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
        argv=(p/'cmdline').read_bytes().decode().rstrip('\x00').split('\x00'))


def atomic_new(path,value):
    path=Path(path);temp=path.with_name(path.name+'.pending')
    if path.exists() or path.is_symlink():raise ValueError('never overwrite launch witness')
    write_new(temp,canonical(value));os.replace(temp,path);sync(path.parent)


def grant(root,ledger_root,contract,attempt_id,timeout=30):
    """Child must call before backend imports/material work,no automatic retry."""
    root=Path(root);me=current_identity();identity(me)
    atomic_new(root/'ready.json',dict(identity=me,attempt_id=attempt_id,plan_sha256=contract['plan_sha256']))
    deadline=time.monotonic()+timeout
    while not (root/'release.json').exists():
        if time.monotonic()>deadline:raise TimeoutError('committed launch release absent')
        time.sleep(.01)
    release=json.loads((root/'release.json').read_bytes());ledger=Ledger(ledger_root,contract)
    state=ledger.read(release['resource_HEAD_sha256']);attempt=state['attempts'][attempt_id]
    if release['attempt_id']!=attempt_id or release['plan_sha256']!=contract['plan_sha256']:
        raise ValueError('launch grant binding')
    if attempt['phase']!='born' or attempt['process']!=me or attempt['birth_sha256']!=release['birth_sha256']:
        raise ValueError('launch grant process/committed birth')
    atomic_new(root/'acknowledged.json',release)
    return release


def ceil_ns(seconds):
    from fractions import Fraction
    q=Fraction(seconds)*10**9
    return -(-q.numerator//q.denominator)


def run(root,ledger_root,contract,attempt_id,command,expected_head,CPU_ns,wall_ns,limits):
    """Exactly one attempt;resource store and numerical store remain separate."""
    root=Path(root)
    if root.exists():raise ValueError('new unused launch directory required')
    if not command or not all(type(x) is str and x for x in command):raise ValueError('full frozen argv required')
    if set(limits)!={'RSS_bytes','AS_bytes','threads','output_bytes','output_roots','reserves'}:raise ValueError('explicit runtime limits')
    if any(type(limits[k]) is not int or limits[k]<=0 for k in ['RSS_bytes','AS_bytes','threads','output_bytes']):raise ValueError('positive runtime limits')
    def reserves():
        for check in limits['reserves']:
            free=os.statvfs(check['path']);available=free.f_bavail*free.f_frsize
            if available<check['minimum_bytes']:raise RuntimeError('disk reserve before/during audit')
    reserves();root.mkdir();supervisor=current_identity();ledger=Ledger(ledger_root,contract)
    frozen=dict(attempt_id=attempt_id,command=command,CPU_ns=CPU_ns,wall_ns=wall_ns,limits=limits,
                plan_sha256=contract['plan_sha256'],supervisor=supervisor,expected_resource_HEAD=expected_head)
    atomic_new(root/'dispatch.json',frozen)
    reservation=ledger.append(dict(kind='reserve',id=attempt_id,CPU_ns=CPU_ns,wall_ns=wall_ns,
        plan_sha256=contract['plan_sha256'],command=command,supervisor=supervisor),expected_head)
    authorized=ledger.append(dict(kind='authorize',id=attempt_id),reservation['HEAD_sha256'])
    head=authorized['HEAD_sha256'];start=time.monotonic_ns();child=None;process=None;birth=None;usage=None;code=None;stop=None
    peak_RSS=peak_threads=0
    def child_limits():
        cpu=max(1,math.ceil(CPU_ns/10**9))
        resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu));resource.setrlimit(resource.RLIMIT_AS,(limits['AS_bytes'],limits['AS_bytes']))
        resource.setrlimit(resource.RLIMIT_FSIZE,(limits['output_bytes'],limits['output_bytes']));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    try:
        with (root/'stdout').open('xb') as out,(root/'stderr').open('xb') as err:
            child=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True,preexec_fn=child_limits,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_MAX_THREADS='1'))
            ready=None
            while True:
                pid,status,found=os.wait4(child.pid,os.WNOHANG)
                if pid:
                    usage=found;code=os.waitstatus_to_exitcode(status);child.returncode=code;break
                now=time.monotonic_ns()
                try:
                    current=current_identity(child.pid)
                    if process is not None and current!=process:raise RuntimeError('audit child identity changed')
                    if ready is None and (root/'ready.json').exists():
                        ready=json.loads((root/'ready.json').read_bytes())
                        if ready!={'identity':current,'attempt_id':attempt_id,'plan_sha256':contract['plan_sha256']} or current['argv']!=command:
                            raise ValueError('actual child ready/birth/argv binding')
                        process=current
                        birth=ledger.append(dict(kind='born',id=attempt_id,process=process),head);head=birth['HEAD_sha256']
                        atomic_new(root/'release.json',dict(attempt_id=attempt_id,plan_sha256=contract['plan_sha256'],
                            resource_HEAD_sha256=head,birth_sha256=birth['HEAD_sha256']))
                    status_fields={k:v.strip() for line in (Path('/proc')/str(child.pid)/'status').read_text().splitlines() if ':' in line for k,v in [line.split(':',1)]}
                    rss=max(int(status_fields.get(k,'0 kB').split()[0]) for k in ['VmRSS','VmHWM'])*1024;threads=int(status_fields.get('Threads','0'))
                    peak_RSS=max(peak_RSS,rss);peak_threads=max(peak_threads,threads)
                    if rss>limits['RSS_bytes']:raise RuntimeError('RSS budget')
                    if threads>limits['threads']:raise RuntimeError('thread budget')
                    outputs={p for r in limits['output_roots'] for p in Path(r).rglob('*')}
                    if any(p.is_symlink() for p in outputs):raise RuntimeError('output alias rejected')
                    total=sum(p.stat().st_size for p in outputs if p.is_file())
                    if total>limits['output_bytes']:raise RuntimeError('launch output budget')
                    reserves()
                    if now-start>wall_ns:raise TimeoutError('audit wall budget')
                    if ready is None and now-start>30*10**9:raise TimeoutError('audit birth handshake absent')
                except (ProcessLookupError,FileNotFoundError):
                    continue  # Reap the same known child,never infer an exit code.
                time.sleep(.02)
        end=time.monotonic_ns()
        if process is None or birth is None:raise RuntimeError('child exited without complete born/release witness')
        if not (root/'acknowledged.json').exists():raise RuntimeError('child did not acknowledge committed grant')
        proof=dict(schema=1,complete=True,attempt_id=attempt_id,plan_sha256=contract['plan_sha256'],
            intent_sha256=authorized['HEAD_sha256'],birth_sha256=birth['HEAD_sha256'],process=process,supervisor=supervisor,
            measurement=contract['measurement'],CPU_ns=ceil_ns(usage.ru_utime)+ceil_ns(usage.ru_stime),wall_ns=end-start,
            start_monotonic_ns=start,end_monotonic_ns=end,return_code=code)
        atomic_new(root/'wait4.json',proof)
        settled=ledger.append(dict(kind='settle',id=attempt_id,evidence=proof,evidence_sha256=digest(canonical(proof))),head)
        head=settled['HEAD_sha256']
        result=dict(status='complete_one_shot_resource_witness' if code==0 and not settled['state']['halted'] else 'execution_or_resource_failure',
            return_code=code,stop_reason=None,process=process,peak_RSS_bytes=peak_RSS,peak_threads=peak_threads,
            resource_HEAD_sha256=head,wait4_sha256=digest(canonical(proof)),CPU_ns=proof['CPU_ns'],wall_ns=proof['wall_ns'])
        atomic_new(root/'terminal.json',result);return result
    except BaseException as error:
        stop=repr(error)
        if child is not None and code is None:
            try:os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            pid,status,usage=os.wait4(child.pid,0);code=os.waitstatus_to_exitcode(status);child.returncode=code
        # Exact measured partial usage is retained but cannot refund a failed
        # handshake/monitor attempt;the reservation stays fully charged.
        failure=dict(status='new_launch_failed_reservation_retained_no_dependent_dispatch',failure=stop,return_code=code,
            resource_HEAD_sha256=head,process=process,CPU_ns=None if usage is None else ceil_ns(usage.ru_utime)+ceil_ns(usage.ru_stime),
            wall_ns=time.monotonic_ns()-start,peak_RSS_bytes=peak_RSS,peak_threads=peak_threads)
        atomic_new(root/'failure.json',failure);raise

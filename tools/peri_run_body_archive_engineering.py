"""Single isolated zero-advance component engineering state, not a solver run."""
import argparse
import datetime
import hashlib
import json
import os
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

from peri_bounded_body_archive_audit import source_fields, build, assess, record


def save(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def identity(pid):
    p = Path('/proc')/str(pid)
    fields = p.joinpath('stat').read_text().rsplit(')',1)[1].split()
    return dict(pid=pid, ppid=int(fields[1]), start_ticks=fields[19],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=p.joinpath('cmdline').read_bytes().decode().strip('\0').split('\0'))


def reserves(root):
    assert shutil.disk_usage('/mnt/f').free >= 150*2**30
    memory = dict(line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())
    assert int(memory['MemAvailable'].split()[0])*1024 >= 8*2**30


def run(root, label, command, cwd, cpu, wall, gib):
    def bounds():
        resource.setrlimit(resource.RLIMIT_CPU, (cpu,cpu))
        resource.setrlimit(resource.RLIMIT_AS, (gib*2**30,gib*2**30))
        resource.setrlimit(resource.RLIMIT_FSIZE, (256*2**20,256*2**20))
        os.umask(0o077)
    start = time.monotonic()
    peak = count = 0
    stop = None
    with (root/(label+'.stdout')).open('xb') as out, (root/(label+'.stderr')).open('xb') as err:
        child = subprocess.Popen(command, cwd=cwd, stdout=out, stderr=err, start_new_session=True,
            preexec_fn=bounds, env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                                      MKL_NUM_THREADS='1', NUMEXPR_MAX_THREADS='1'))
        witness = None
        # Popen may return before exec. Never bind a launch witness to the
        # transient pre-exec Python argv, then mistake exec for PID reuse.
        while child.poll() is None:
            try:
                candidate = identity(child.pid)
                if candidate['argv'] == command:
                    witness = candidate
                    break
            except FileNotFoundError:
                pass
            if time.monotonic()-start > 2:
                break
            time.sleep(.001)
        save(root/(label+'.launch.json'), dict(command=command, cwd=str(cwd), identity=witness, started_at=stamp()))
        while child.poll() is None:
            try:
                if witness:
                    assert identity(child.pid) == witness
                status = Path('/proc',str(child.pid),'status').read_text().splitlines()
                rss = next((int(x.split()[1]) for x in status if x.startswith('VmRSS:')),0)
                threads = next((int(x.split()[1]) for x in status if x.startswith('Threads:')),0)
                peak = max(peak,rss)
                count = max(count,threads)
                reserves(root)
                assert rss <= gib*1024**2 and threads <= 8
                assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) <= 2**30
                assert time.monotonic()-start <= wall
            except FileNotFoundError:
                continue
            except Exception as error:
                stop = repr(error)
                os.killpg(child.pid, signal.SIGKILL)
                break
            time.sleep(.1)
        code = child.wait()
        out.flush(); err.flush(); os.fsync(out.fileno()); os.fsync(err.fileno())
    result = dict(return_code=code, stop_reason=stop, wall_seconds=time.monotonic()-start,
        peak_parent_RSS_KiB=peak, peak_parent_threads=count, identity=witness,
        stdout=record(root/(label+'.stdout')), stderr=record(root/(label+'.stderr')))
    save(root/(label+'.terminal.json'), result)
    assert code == 0 and stop is None, result
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    assert root.parent == Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material')
    assert not root.exists(), 'never redispatch a component engineering state'
    reserves(root)
    repo = Path(__file__).resolve().parents[1]
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True) == ''
    root.mkdir(mode=0o700)
    try:
        source = subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        raw = root.parent/'study-v19-early-friction-cycle/actual_PD_early_biaxial_slip_cycle.jsonl'
        files = ['tools/peri_bounded_body_archive_audit.py','tools/peri_canonical_archive_audit.py',
                 'tools/peri_run_body_archive_engineering.py','tools/bounded_archive_body_checks.cpp',
                 'src/shared/peridynamics/pd_bounded_archive_body.h','docs/peri/BOUNDED_BODY_ARCHIVE_20261105.md']
        frozen = root/'frozen'
        for name in files:
            target = frozen/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(repo/name,target)
        # Hash all local compile inputs before compiling; MD selects the actual
        # set afterward. Imported system headers are also frozen before native.
        before = {str(p.resolve()):record(p) for p in (repo/'src').rglob('*') if p.suffix in ['.h','.hpp','.cpp']}
        snapshot = []
        for name in files:
            snapshot += [record(repo/name),record(frozen/name)]
        fixture = root/'nonzero_body.bin'
        fields = source_fields(raw)
        with fixture.open('xb') as stream:
            stream.write(build(fields));stream.flush();os.fsync(stream.fileno())
        builddir = repo/'build/peri-latest-20261103'
        original_object = 'extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_bounded_gel_friction_probe.dir/bounded_gel_friction_probe.cpp.o'
        original_binary = 'extra_3d_PD_framework_tests/bin/extra_3d_PD_bounded_gel_friction_probe'
        original = subprocess.check_output(['ninja','-t','commands','-s',original_object],cwd=builddir,text=True).strip()
        command = shlex.split(original)
        for flag, value in [('-o',root/'body.o'),('-MF',root/'body.d'),('-MT',root/'body.o'),('-c',frozen/'tools/bounded_archive_body_checks.cpp')]:
            command[command.index(flag)+1] = str(value)
        command += ['-Wall','-Wextra','-Wno-misleading-indentation','-Wno-unused-parameter','-Werror']
        link_text = subprocess.check_output(['ninja','-t','commands','-s',original_binary],cwd=builddir,text=True).strip()
        assert link_text.startswith(': && ') and link_text.endswith(' && :')
        link = shlex.split(link_text[5:-5])
        assert link.count(original_object) == 1
        link[link.index(original_object)] = str(root/'body.o')
        link[link.index('-o')+1] = str(root/'body_checks')
        plan = dict(status='prepared_component_only', prepared_at=stamp(), source_commit=source,
            physical_advances=0, raw=record(raw), fixture=record(fixture), source_snapshots=snapshot,
            compile_command=command, link_command=link, build_cwd=str(builddir),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),
            native_budget=dict(CPU_seconds=30,wall_seconds=60,AS_GiB=1,threads=1),
            file_bound_MiB=256,total_bound_GiB=1,host_reserve_GiB=120,output_reserve_GiB=150,
            gates=['all source-derived field IEEE bits and array counts exactly match',
                   'full archive byte equality','18 specific valid-hash semantic refusals',
                   'baseline unchanged after every refused decode','nonzero viscous, tangent, travel correction',
                   'zero physical advance and zero fresh reference construction'],
            failure_policy='Preserve exact failure/returncode/inputs and stop this state; no solver invocation/retry/fallback/reset/gate relaxation.',
            full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan)
        save(root/'state.json',dict(status='single_component_engineering_attempted',identity=identity(os.getpid()),started_at=stamp()))
        compile_result = run(root,'compile',command,builddir,180,240,2)
        link_result = run(root,'link',link,builddir,30,60,2)
        dependency_text = (root/'body.d').read_text().replace('\\\n',' ')
        dependencies = [Path(x).resolve() for x in shlex.split(dependency_text.split(':',1)[1])]
        # Freeze actual original headers/libraries + their copy records; every
        # working-tree input must still match its pre-compilation identity.
        linked = [builddir/x for x in link if x.endswith(('.a','.so')) or '.so.' in x]
        dynamic_text = subprocess.check_output(['ldd',str(root/'body_checks')],text=True)
        with (root/'dynamic_libraries.txt').open('x') as stream:
            stream.write(dynamic_text)
        for line in dynamic_text.splitlines():
            tokens = line.strip().split()
            paths = [word for word in tokens if word.startswith('/')]
            assert 'not found' not in line
            linked += [Path(word) for word in paths]
        inputs = sorted(set(dependencies+linked+[Path(command[0]),Path(link[0])]),key=str)
        support = []
        for i,p in enumerate(inputs):
            original_record = record(p)
            if str(p) in before:
                assert original_record == before[str(p)], ('compile input changed',p)
            target = root/'support'/str(i)/p.name
            target.parent.mkdir(parents=True,exist_ok=False)
            shutil.copy2(p,target)
            support += [original_record,record(target)]
        for item in snapshot:
            assert record(item['path']) == item
        frozen_execution = dict(recorded_at=stamp(),plan=record(root/'plan.json'),ELF=record(root/'body_checks'),
            fixture=record(fixture),support=support,source_snapshots=snapshot,compile=compile_result,link=link_result)
        save(root/'frozen_execution.json',frozen_execution)
        output = root/'evidence'
        output.mkdir()
        native = run(root,'native',[str(root/'body_checks'),str(fixture),str(output)],root,30,60,1)
        for item in support+snapshot+[frozen_execution['ELF'],frozen_execution['fixture'],frozen_execution['plan']]:
            assert record(item['path']) == item
        assert record(raw) == plan['raw']
        native_report = json.loads((root/'native.stdout').read_text())
        if native['identity']:
            assert native_report['native_identity'] == native['identity']
        else:
            assert native_report['native_identity']['argv'] == [str(root/'body_checks'),str(fixture),str(output)]
            assert native_report['native_identity']['ppid'] == os.getpid()
            assert native_report['native_identity']['boot_id'] == identity(os.getpid())['boot_id']
        result = assess(raw,fixture,output,native_report)
        result.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),
                      frozen_execution=record(root/'frozen_execution.json'),native=native,
                      actual_compiler_dependency_count=len(dependencies),support_original_and_copy_records=len(support))
        save(root/'engineering_receipt.json',result)
        save(root/'terminal.json',dict(status=result['status'],finished_at=stamp(),return_code=0,
            receipt=record(root/'engineering_receipt.json'),physical_advances=0,full_restart_qualified=False))
        print(json.dumps(dict(status=result['status'],field_count=result['field_count'],refusals=18)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),
            physical_advances=0,full_restart_qualified=False,S3_complete=False))
        raise


if __name__ == '__main__':
    main()

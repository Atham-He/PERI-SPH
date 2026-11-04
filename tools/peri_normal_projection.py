#!/usr/bin/env python3
"""Freeze/run/review bounded normal-constraint mathematics, never a gel study."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess


def read(path): return json.loads(path.read_text())
def digest(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream, 'sha256').hexdigest()
def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n'); temporary.replace(path)
def identity():
    tail = Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()
    return dict(pid=os.getpid(), ppid=os.getppid(), start_ticks=tail[19],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=Path('/proc/self/cmdline').read_bytes().decode().rstrip('\0').split('\0'))
def resources(root):
    if not root.resolve().is_relative_to(Path('/mnt/f').resolve()) or shutil.disk_usage('/mnt/f').free < 150 * 2**30:
        raise ValueError('Prototype output requires the real F host and >=150GiB reserve (also covers120GiB host reserve)')


def row(targets=(0,), weights=(1,), normal=1, old=0, new=0):
    return dict(targets=list(targets), weights=list(weights), normal=normal, old_offset=old, new_offset=new)
def probe(mass, q, v, dt, rows, **extra):
    value = dict(mass=mass, position=q, velocity=v, dt=dt, rows=rows, velocity_scale=1, tolerance=1e-12, max_sweeps=10000, expected_code=0)
    value.update(extra)
    return value


def prepare(args):
    root, executable = args.root.resolve(), args.executable.resolve(strict=True)
    if root.exists(): raise FileExistsError('Prototype states cannot be reused')
    resources(root.parent)
    repo = Path(__file__).resolve().parents[1]
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo).strip(): raise ValueError('Commit the prototype before freezing it')
    cases = dict(
        free=probe([1, 2], [1, 1], [-.1, .2], .1, []),
        fixed=probe([2], [.1], [-1], .2, [row()]),
        moving=probe([2], [0], [0], .2, [row(new=.1)]),
        signed_mls=probe([2, 1], [.1, .1], [-1, -1], .2, [row((0, 1), (1.25, -.25))]),
        two_walls=probe([1, 1], [0, 1], [-1, 1], .1, [row(), row((1,), normal=-1, old=-1, new=-1)]),
        redundant=probe([1], [0], [-2], .1, [row(), row()]),
        coupled=probe([2, 3], [0, 0], [-1, -.5], .1, [row((0, 1), (.7, .3)), row((0, 1), (.2, .8))]),
        si_units=probe([1e-9], [1e-6], [-1e-3], .002, [row()]),
        infeasible=probe([1], [.5], [0], .1, [row(new=.7), row(normal=-1, old=-1, new=-.3)]),
        invalid=probe([1], [-.1], [0], .1, [row()]))
    cases['si_units']['velocity_scale'] = .001
    cases['infeasible'].update(expected_code=2, max_sweeps=64)
    cases['invalid']['expected_code'] = 1
    root.mkdir(); (root / 'frozen_tools').mkdir(); (root / 'frozen_source').mkdir()
    shutil.copy2(executable, root / 'native_probe')
    names = ['src/shared/peridynamics/pd_normal_projection.h',
             'tests/3d_examples/extra_3d_PD_framework_tests/normal_contact_probe.cpp',
             'tests/3d_examples/extra_3d_PD_framework_tests/normal_projection_tests.cpp',
             'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
             'tools/peri_normal_projection.py', 'tools/peri_normal_projection_audit.py',
             'docs/peri/NORMAL_CONSTRAINT_PROTOTYPE_20261104.md']
    sources = {}
    for name in names:
        target = root / 'frozen_source' / name; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / name, target); sources[name] = digest(target)
    for name in ('peri_normal_projection.py', 'peri_normal_projection_audit.py'):
        shutil.copy2(repo / 'tools' / name, root / 'frozen_tools' / name)
    libraries = {}
    for line in subprocess.check_output(['ldd', str(root / 'native_probe')]).decode().splitlines():
        for token in line.split():
            if token.startswith('/'):
                path = Path(token).resolve(strict=True); libraries[str(path)] = digest(path); break
    plan = dict(scope='Low-dimensional normal-only perfectly inelastic contact mathematics; not integrated into PD/gel and not calibration release',
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo).decode().strip(),
        source_sha256=sources, executable_sha256=digest(root / 'native_probe'), runtime_library_sha256=libraries,
        support_sha256={path.name: digest(path) for path in (root / 'frozen_tools').iterdir()}, cases=cases,
        impact_trajectories=[dict(name=f'impact_{i}', dt=dt, end_time=1.7, mass=1, initial_position=1, initial_velocity=-1)
                             for i, dt in enumerate([.36, .18, .09, .045])],
        resources=dict(max_native_calls=256, native_timeout_s=5, native_RSS_MiB=256, total_output_MiB=32, controller_RSS_MiB=1024),
        independent_tolerance=1e-10, policy='No existing state rerun, fallback or initial penetration repair; unexpected failures stop subsequent probes',
        production_integrated=False, material_calibration_allowed=False)
    save(root / 'plan.json', plan); save(root / 'state.json', dict(status='prepared', plan_sha256=digest(root / 'plan.json'), calls=[]))
    print(json.dumps(dict(status='prepared', root=str(root), source_commit=plan['source_commit'])))


def verify(root, plan, state):
    if digest(root / 'plan.json') != state['plan_sha256'] or digest(root / 'native_probe') != plan['executable_sha256']:
        raise ValueError('Prototype plan/ELF changed')
    for name, expected in plan['source_sha256'].items():
        if digest(root / 'frozen_source' / name) != expected: raise ValueError('Frozen prototype source changed')
    for name, expected in plan['support_sha256'].items():
        if digest(root / 'frozen_tools' / name) != expected: raise ValueError('Frozen independent support changed')
    for name, expected in plan['runtime_library_sha256'].items():
        if digest(Path(name)) != expected: raise ValueError('Prototype runtime library changed')


def encode(case):
    def numbers(values): return ' '.join(format(x, '.17g') for x in values)
    lines = [numbers([len(case['mass']), len(case['rows']), case['dt'], case['velocity_scale'], case['tolerance'], case['max_sweeps']])]
    lines += [numbers(case[key]) for key in ('mass', 'position', 'velocity')]
    for c in case['rows']:
        lines.append(numbers([c['normal'], c['old_offset'], c['new_offset'], len(c['targets'])] +
                             [value for pair in zip(c['targets'], c['weights']) for value in pair]))
    return '\n'.join(lines) + '\n'


def run(root):
    root = root.resolve(strict=True); plan, state = read(root / 'plan.json'), read(root / 'state.json')
    if state['status'] != 'prepared': raise ValueError('Existing prototype state cannot be redispatched')
    verify(root, plan, state); resources(root)
    state.update(status='running', controller=identity()); save(root / 'state.json', state)
    (root / 'raw').mkdir()
    def call(name, case):
        resources(root)
        if len(state['calls']) >= plan['resources']['max_native_calls']: raise ValueError('native call budget exceeded')
        if digest(root / 'native_probe') != plan['executable_sha256']: raise ValueError('ELF changed before probe')
        directory = root / 'raw' / name; directory.mkdir(exist_ok=False)
        save(directory / 'input.json', case); (directory / 'input.txt').write_text(encode(case))
        command = [str(root / 'native_probe'), str(directory / 'input.txt')]
        state['active_case'] = name; save(root / 'state.json', state)
        with (directory / 'stdout.json').open('wb') as stdout, (directory / 'stderr.txt').open('wb') as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=plan['resources']['native_timeout_s'])
        save(directory / 'return.json', dict(return_code=process.returncode, command=command, recorded_at_utc=datetime.now(timezone.utc).isoformat()))
        output = read(directory / 'stdout.json'); actual = output['identity']
        if actual['argv'] != command or actual['executable'] != str(root / 'native_probe') or actual['ppid'] != os.getpid() or \
                actual['boot_id'] != state['controller']['boot_id'] or not actual['start_ticks'].isdigit():
            raise ValueError('Actual native OS identity/argv mismatch')
        state['calls'].append(dict(name=name, expected_code=case['expected_code'], return_code=process.returncode, identity=actual))
        save(root / 'state.json', state)
        if process.returncode != case['expected_code']: raise ValueError('Unexpected native failure or false negative-control pass')
        if process.returncode == 0 and not output['converged']: raise ValueError('Native convergence gate failed')
        if resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss > 256*1024 or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 1024*1024:
            raise ValueError('Prototype memory budget exceeded')
        return output
    try:
        for name, case in plan['cases'].items(): call(name, case)
        for trajectory in plan['impact_trajectories']:
            t, q, v, step = 0., trajectory['initial_position'], trajectory['initial_velocity'], 0
            while t < trajectory['end_time']:
                dt = min(trajectory['dt'], trajectory['end_time']-t)
                case = probe([trajectory['mass']], [q], [v], dt, [row()])
                output = call(trajectory['name'] + f'_step_{step:03d}', case)
                q, v = output['position'][0], output['velocity'][0]
                t = trajectory['end_time'] if dt == trajectory['end_time']-t else t+dt; step += 1
        state['status'] = 'normal_prototype_runtime_controls_completed'
    except Exception as error:
        state.update(status='prototype_stopped_for_review', failure=str(error)); raise
    finally:
        state.update(active_case=None, finished_at_utc=datetime.now(timezone.utc).isoformat())
        files = [path for path in (root / 'raw').rglob('*') if path.is_file()]
        state['artifact_index'] = {path.relative_to(root).as_posix(): dict(size=path.stat().st_size, sha256=digest(path)) for path in files}
        if sum(path.stat().st_size for path in files) > 32*2**20: state.update(status='prototype_stopped_for_review', failure='output budget exceeded')
        save(root / 'state.json', state)
    print(json.dumps(dict(status=state['status'], native_calls=len(state['calls']))))


def review(root):
    from peri_normal_projection_audit import check_native
    root = root.resolve(strict=True); plan, state = read(root / 'plan.json'), read(root / 'state.json')
    if state['status'] in ('prepared', 'running'): print(json.dumps(dict(status='waiting'))); return
    output = root / 'independent_review.json'; state_sha = digest(root / 'state.json')
    if output.exists():
        receipt = read(output)
        if receipt['state_sha256'] != state_sha: raise ValueError('Reviewed prototype state changed')
        print(json.dumps(dict(status=receipt['status'], idempotent=True))); return
    receipt = dict(state_sha256=state_sha, plan_sha256=digest(root / 'plan.json'), checks={}, impact_trajectories={},
                   production_integrated=False, material_calibration_allowed=False)
    try:
        verify(root, plan, state)
        for name, expected in state['artifact_index'].items():
            path = root / name
            if path.stat().st_size != expected['size'] or digest(path) != expected['sha256']: raise ValueError('Prototype raw size/SHA changed')
        for call in state['calls']:
            directory = root / 'raw' / call['name']; case = read(directory / 'input.json')
            if call['name'] in plan['cases'] and case != plan['cases'][call['name']]: raise ValueError('Probe inputs differ from frozen scalar case')
            if (directory / 'input.txt').read_text() != encode(case): raise ValueError('Actual native input differs from captured numerical input')
            native = read(directory / 'stdout.json'); returned = read(directory / 'return.json'); code = returned['return_code']
            actual = native['identity']
            if not native['prototype_only'] or actual['argv'] != returned['command'] or actual['executable'] != str(root / 'native_probe') or \
                    actual['ppid'] != state['controller']['pid'] or actual['boot_id'] != state['controller']['boot_id'] or not actual['start_ticks'].isdigit():
                raise ValueError('Independent actual native identity/argv check failed')
            receipt['checks'][call['name']] = check_native(case, native, code)
        for trajectory in plan['impact_trajectories']:
            calls = [call for call in state['calls'] if call['name'].startswith(trajectory['name'] + '_step_')]
            q, v, t, impulses, losses, peak_force = trajectory['initial_position'], trajectory['initial_velocity'], 0., [], [], 0.
            for call in calls:
                directory = root / 'raw' / call['name']; case, native = read(directory / 'input.json'), read(directory / 'stdout.json')
                dt = min(trajectory['dt'], trajectory['end_time']-t)
                if case != probe([trajectory['mass']], [q], [v], dt, [row()]):
                    raise ValueError('Trajectory was not advanced from its own actual previous output')
                impulses.append(native['normal_impulse'][0]); losses.append(native['contact_step_dissipation'])
                peak_force = max(peak_force, native['normal_impulse'][0]/dt)
                q, v = native['position'][0], native['velocity'][0]
                if q < -1e-10: raise ValueError('Accepted prototype trajectory penetrated')
                t = trajectory['end_time'] if dt == trajectory['end_time']-t else t+dt
            if t != trajectory['end_time'] or abs(math.fsum(impulses)-1) > 1e-10 or abs(math.fsum(losses)-.5) > 1e-10 or abs(q) > 1e-10 or abs(v) > 1e-10:
                raise ValueError('Analytic perfectly inelastic impact impulse/energy limit failed')
            receipt['impact_trajectories'][trajectory['name']] = dict(dt=trajectory['dt'], steps=len(calls), final_position=q,
                final_velocity=v, total_normal_impulse=math.fsum(impulses), contact_step_dissipation=math.fsum(losses), peak_step_force=peak_force)
        receipt.update(status='independent_normal_projection_prototype_passed' if state['status'] == 'normal_prototype_runtime_controls_completed' else
                       'independent_prototype_review_completed_remains_unqualified', runtime_failure=state.get('failure'), scope=plan['scope'])
    except Exception as error:
        receipt.update(status='independent_prototype_check_failed', failure=str(error)); raise
    finally: save(output, receipt)
    print(json.dumps(dict(status=receipt['status'], checked_native_calls=len(receipt['checks']), impact_trajectories=receipt['impact_trajectories']), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare'); prep.add_argument('--root', type=Path, required=True); prep.add_argument('--executable', type=Path, required=True)
    for action in ('run', 'review'): sub.add_parser(action).add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare': prepare(args)
    elif args.action == 'run': run(args.root)
    else: review(args.root)
if __name__ == '__main__': main()

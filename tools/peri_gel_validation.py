#!/usr/bin/env python3
"""Versioned gel time/contact/space and full-history restart verification.

New numerical controls only. Failed cases stop the sequence; existing states
cannot be dispatched again. No material calibration is performed.
"""
import argparse
import csv
import hashlib
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

from peri_run import digest, parameters, source_snapshot_files


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def resources(root, host):
    if not root.resolve().is_relative_to(host.resolve()):
        raise ValueError('Verification output must reside on the declared actual host volume')
    if shutil.disk_usage(host).free < 120 * 2**30 or shutil.disk_usage(root).free < 150 * 2**30:
        raise ValueError('disk_reserve: host 120 GiB / output filesystem 150 GiB required')


def identity(pid):
    proc = Path('/proc') / str(pid)
    argv = (proc / 'cmdline').read_bytes().decode().rstrip('\0').split('\0')
    fields = (proc / 'stat').read_text().rsplit(')', 1)[1].split()
    return dict(pid=pid, ppid=int(fields[1]), start_ticks=fields[19],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=argv, executable=str((proc / 'exe').resolve(strict=True)),
                argv_observed=bool(argv and argv[0]))


def rows(path):
    with path.open(newline='') as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def force_error(a, b):
    """Weighted L2 error within actual shared support; never extrapolate."""
    for curve in (a, b):
        if len(curve) < 2 or any(v['stroke_m'] <= u['stroke_m'] for u, v in zip(curve, curve[1:])):
            raise ValueError('force comparison needs strictly increasing supported stroke')
    if a[0]['stroke_m'] != b[0]['stroke_m'] or a[-1]['stroke_m'] != b[-1]['stroke_m']:
        raise ValueError('force curves do not cover the same actual stroke')
    errors, reference, grid = [], [], []
    j = 0
    for row in b:
        x = row['stroke_m']
        while j + 1 < len(a) - 1 and a[j + 1]['stroke_m'] < x:
            j += 1
        left, right = a[j], a[j + 1]
        t = (x - left['stroke_m']) / (right['stroke_m'] - left['stroke_m'])
        value = left['upper_observed_Fy_N'] + t * (right['upper_observed_Fy_N'] - left['upper_observed_Fy_N'])
        errors.append((value - row['upper_observed_Fy_N'])**2)
        reference.append(row['upper_observed_Fy_N']**2)
        grid.append(x)
    def integrate(values):
        return math.fsum((u + v) * (y - x) / 2 for u, v, x, y in zip(values, values[1:], grid, grid[1:]))
    denominator = integrate(reference)
    if denominator <= 1e-30:
        raise ValueError('UNRESOLVED: force denominator is at the noise floor')
    return math.sqrt(integrate(errors) / denominator)


def prepare(args):
    root, host = args.root.resolve(), args.host_storage.resolve(strict=True)
    if root.exists():
        raise FileExistsError('A validation directory cannot be reused')
    resources(root.parent, host)
    values = parameters(args.config)
    if not values.get('gel_viscoelastic') or values.get('viscoplastic'):
        raise ValueError('This verification requires the recoverable gel material')
    repo = Path(__file__).resolve().parents[1]
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo).strip():
        raise ValueError('Commit the reviewed source before freezing validation inputs')
    root.mkdir()
    source_files = source_snapshot_files(repo)
    source_hashes = {str(path.relative_to(repo)): digest(path) for path in source_files}
    with tarfile.open(root / 'source_snapshot.tar.gz', 'w:gz', dereference=True) as archive:
        for path in source_files:
            archive.add(path, arcname=str(path.relative_to(repo)), recursive=False)
    shutil.copy2(args.executable, root / 'solver')
    (root / 'frozen_tools').mkdir()
    for name in ['peri_gel_validation.py', 'peri_run.py', 'peri_gel_audit.py']:
        shutil.copy2(repo / 'tools' / name, root / 'frozen_tools' / name)
    # The physical nonlocal horizon is 40um at every spatial level.
    values.update(horizon_ratio=4, pd_cfl=0.1, output_frames=40)
    if args.surface_contact:
        values.update(reference_surface_contact=1, surface_support=4e-5, contact_penalty_length=4e-5,
                      surface_rings=128, surface_azimuths=16)
    cases = []
    specifications = [
        ('baseline', {}), ('dt_half', {'pd_cfl': 0.05}), ('dt_quarter', {'pd_cfl': 0.025}),
        ('contact_half', {'pd_cfl': 0.025, 'friction_velocity': values['friction_velocity'] / 2}),
        ('contact_quarter', {'pd_cfl': 0.025, 'friction_velocity': values['friction_velocity'] / 4}),
    ]
    if args.surface_contact:
        specifications += [
            ('quadrature_32', {'pd_cfl': 0.025, 'surface_rings': 32}),
            ('quadrature_64', {'pd_cfl': 0.025, 'surface_rings': 64}),
            ('space_medium', {'pd_cfl': 0.025, 'spacing': 7.5e-6, 'horizon_ratio': 4e-5 / 7.5e-6}),
            ('space_fine', {'pd_cfl': 0.025, 'spacing': 5e-6, 'horizon_ratio': 8}),
        ]
    else:
        specifications += [
            ('space_coarse', {'pd_cfl': 0.025, 'spacing': 2e-5, 'horizon_ratio': 2}),
            ('space_medium', {'pd_cfl': 0.025, 'spacing': 1.5e-5, 'horizon_ratio': 4e-5 / 1.5e-5})]
    for name, edits in specifications:
        configuration = dict(values, **edits)
        file = root / (name + '.cfg')
        file.write_text(''.join(f'{key}={value:.17g}\n' for key, value in sorted(configuration.items())))
        cases.append(dict(name=name, configuration=file.name, configuration_sha256=digest(file)))
    # Restart uses the exact baseline inputs and executable, not a new reference.
    plan = dict(schema_version=1, scope='small effective-gel numerical control; not original-size or experimental calibration',
                source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo).decode().strip(),
                source_sha256=source_hashes, source_snapshot_sha256=digest(root / 'source_snapshot.tar.gz'),
                executable_sha256=digest(root / 'solver'), host_storage=str(host), cases=cases,
                tolerances=dict(force_relative=0.02, mass_relative=1e-10, max_penetration_over_spacing=0.1),
                resource_limits=dict(host_free_GiB=120, output_free_GiB=150, memory_GiB=8, steps=98304,
                                     threads=8, particles=1000, seconds_per_case=3600, output_bytes_per_case=2 * 2**30),
                policy='Stop on any execution/integrity/resource/accuracy failure; no retry, fallback or reference refresh.',
                reference_surface_contact=args.surface_contact, material_calibration_allowed=False)
    if args.surface_contact:
        plan['resource_limits'].update(particles=8192, output_bytes_per_case=4 * 2**30)
        plan['prospective_cost_basis'] = '912 nodes/920 steps: 2.764s, 9.94MiB CP; fixed-horizon dx/2 predicts ~128x work/~64x neighbor storage; old R8 guard unchanged.'
    plan['support_sha256'] = {p.name: digest(p) for p in (root / 'frozen_tools').iterdir()}
    libraries = []
    for line in subprocess.check_output(['ldd', str(root / 'solver')]).decode().splitlines():
        for token in line.split():
            if token.startswith('/'):
                library = Path(token).resolve(strict=True)
                libraries.append(library)
                break
    plan['runtime_library_sha256'] = {str(path): digest(path) for path in libraries}
    save(root / 'plan.json', plan)
    save(root / 'state.json', dict(status='prepared', plan_sha256=digest(root / 'plan.json'), results=[]))
    print(json.dumps(dict(root=str(root), status='prepared', source_commit=plan['source_commit'])))


def validate_case(directory, configuration, tolerance):
    data, trace = rows(directory / 'output/diagnostics.csv'), rows(directory / 'output/step_trace.csv')
    status = (directory / 'output/run_status.txt').read_text().splitlines()
    if not status or status[0] != 'completed':
        raise ValueError('case did not report an accepted complete trajectory')
    duration = configuration['compression_stroke'] / configuration['compression_speed']
    effective = parameters(directory / 'output/effective_parameters.cfg')
    if any(effective.get(key) != value for key, value in configuration.items()):
        raise ValueError('actual effective inputs differ from the frozen configuration')
    if not data or data[-1]['time_s'] != duration or data[-1]['stroke_m'] != configuration['compression_stroke']:
        raise ValueError('incomplete actual time/stroke support')
    if not trace or len(trace) != int(data[-1]['pd_steps']):
        raise ValueError('missing accepted-step trace')
    for row in data + trace:
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError('non-finite independent diagnostic value')
    if max(abs(row['mass_relative_error']) for row in data) > tolerance['mass_relative'] or \
       max(row['particle_mass_error'] for row in data) > tolerance['mass_relative']:
        raise ValueError('mass integrity failed')
    if max(row['upper_penetration_over_spacing'] for row in data) > tolerance['max_penetration_over_spacing']:
        raise ValueError('contact penetration accuracy failed')
    if any(row['material_energy_available'] != 1 for row in data):
        raise ValueError('gel free energy/dissipation unavailable')
    if any(v['gel_dissipation_J'] < u['gel_dissipation_J'] for u, v in zip(data, data[1:])):
        raise ValueError('physical gel dissipation decreased')
    prior = 0
    for number, row in enumerate(trace, 1):
        if row['pd_steps'] != number or row['dt_s'] <= 0 or row['time_s'] <= prior:
            raise ValueError('invalid accepted-step cursor')
        if not math.isclose(row['time_s'] - prior, row['dt_s'], rel_tol=1e-10, abs_tol=1e-20):
            raise ValueError('actual step length does not match the time cursor')
        if not math.isclose(row['material_history_time_s'], row['time_s'], rel_tol=1e-12, abs_tol=1e-20):
            raise ValueError('accepted material memory has a different time cursor')
        prior = row['time_s']
    return dict(accepted_steps=len(trace), initial_actual_dt_s=trace[0]['dt_s'],
                max_penetration_over_spacing=max(row['upper_penetration_over_spacing'] for row in data),
                final_gel_dissipation_J=data[-1]['gel_dissipation_J'],
                final_mechanical_residual_J=data[-1]['mechanical_residual_J'])


def run(root):
    root = root.resolve(strict=True)
    plan, state = read(root / 'plan.json'), read(root / 'state.json')
    if state['status'] != 'prepared':
        raise ValueError('Existing validation state cannot be dispatched again')
    if digest(root / 'plan.json') != state['plan_sha256'] or digest(root / 'solver') != plan['executable_sha256']:
        raise ValueError('Frozen plan/executable changed')
    for name, expected in plan['support_sha256'].items():
        if digest(root / 'frozen_tools' / name) != expected:
            raise ValueError('Frozen verification support changed')
    for name, expected in plan['runtime_library_sha256'].items():
        if digest(Path(name)) != expected: raise ValueError('Frozen dynamic dependency changed')
    host = Path(plan['host_storage'])
    environment = dict(os.environ, PERI_HOST_STORAGE=str(host))
    state.update(status='running', controller=identity(os.getpid()))
    save(root / 'state.json', state)
    results = {}

    def execute(name, configuration, *, steps=0, stop=False, resume=False, expected_code=0):
        resources(root, host)
        frozen_case = next(case for case in plan['cases'] if root / case['configuration'] == configuration)
        if digest(configuration) != frozen_case['configuration_sha256'] or digest(root / 'solver') != plan['executable_sha256']:
            raise ValueError('frozen solver/case changed before dispatch')
        if digest(root / 'source_snapshot.tar.gz') != plan['source_snapshot_sha256']:
            raise ValueError('frozen source archive changed before dispatch')
        for library, expected in plan['runtime_library_sha256'].items():
            if digest(Path(library)) != expected: raise ValueError('runtime dependency changed before dispatch')
        directory = root / 'runs' / name
        if not resume:
            directory.mkdir(parents=True, exist_ok=False)
        command = [str(root / 'solver'), '--case-config', str(configuration), '--state_recording=true', '--pd-threads=8',
                   '--pd-checkpoint-seconds=0', f'--pd-checkpoint-steps={steps}']
        if stop: command.append('--pd-stop-after-checkpoint')
        if resume: command.append('--pd-resume=checkpoints')
        state.update(active_case=name, active_command=command)
        with (directory / 'run.log').open('a') as log:
            process = subprocess.Popen(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT)
            try:
                state['active_process'] = identity(process.pid)
            except FileNotFoundError:
                state['active_process'] = None
            save(root / 'state.json', state)
            try:
                code = process.wait(timeout=3605)
            except subprocess.TimeoutExpired:
                process.terminate()
                try: process.wait(timeout=10)
                except subprocess.TimeoutExpired: process.kill(); process.wait()
                code = 124
        if code != expected_code:
            state['last_return_code'] = code
            raise ValueError(f'{name}: expected return {expected_code}, actual {code}; original run.log/status preserved')
        if sum(path.stat().st_size for path in directory.rglob('*') if path.is_file()) > plan['resource_limits']['output_bytes_per_case']:
            raise ValueError(name + ': output cost budget exceeded')
        return directory

    def complete(case):
        file = root / case['configuration']
        if digest(file) != case['configuration_sha256']:
            raise ValueError('frozen case inputs changed')
        directory = execute(case['name'], file)
        metrics = validate_case(directory, parameters(file), plan['tolerances'])
        state['results'].append(dict(name=case['name'], status='checked', **metrics))
        results[case['name']] = rows(directory / 'output/diagnostics.csv')
        save(root / 'state.json', state)
        print(json.dumps(dict(case=case['name'], **metrics)), flush=True)

    def compare(group, coarse, fine):
        value = force_error(results[coarse], results[fine])
        state.setdefault('comparisons', []).append(dict(group=group, coarse=coarse, fine=fine, relative_force_error=value))
        save(root / 'state.json', state)
        if value > plan['tolerances']['force_relative']:
            raise ValueError(f'{group}: force difference {value:.12g} exceeds predefined 0.02')

    def latest(directory):
        return directory / 'checkpoints' / (directory / 'checkpoints/LATEST').read_text().strip()

    try:
        for case in plan['cases'][:3]: complete(case)
        initial = [row['initial_actual_dt_s'] for row in state['results']]
        if not math.isclose(initial[1] / initial[0], 0.5, rel_tol=1e-12) or not math.isclose(initial[2] / initial[1], 0.5, rel_tol=1e-12):
            raise ValueError('actual initial steps were not halved')
        compare('time', 'baseline', 'dt_quarter'); compare('time', 'dt_half', 'dt_quarter')

        baseline = root / 'runs/baseline'
        split = execute('restart_split', root / 'baseline.cfg', steps=max(1, state['results'][0]['accepted_steps'] // 2), stop=True, expected_code=75)
        saved = read(latest(split) / 'manifest.json')
        evidence = root / 'restart_evidence'
        shutil.copytree(latest(split), evidence / 'before_noadvance')
        if saved['time_s'] <= 0:
            raise ValueError('restart fixture did not preserve nonzero history')
        execute('restart_split', root / 'baseline.cfg', stop=True, resume=True, expected_code=75)
        restored = read(latest(split) / 'manifest.json')
        shutil.copytree(latest(split), evidence / 'after_noadvance')
        if any(saved[key] != restored[key] for key in ['physical_state_fnv1a64', 'time_s', 'pd_steps', 'completed_frame']):
            raise ValueError('no-advance restore changed physical history/cursor')
        for name in ['diagnostics.csv', 'step_trace.csv']:
            with (split / 'output' / name).open('a') as stream: stream.write('999,uncommitted-tail\n')
        sentinel = split / 'output/GelParticle_ite_1000001234.vtp'
        sentinel.write_text('uncommitted snapshot tail; preserve outside active output')
        (evidence / 'uncommitted_output').mkdir()
        for file in [split / 'output/diagnostics.csv', split / 'output/step_trace.csv', sentinel]:
            shutil.copy2(file, evidence / 'uncommitted_output' / file.name)
        execute('restart_split', root / 'baseline.cfg', resume=True)
        flags = {name: (baseline / 'output' / name).read_bytes() == (split / 'output' / name).read_bytes()
                 for name in ['diagnostics.csv', 'step_trace.csv']}
        flags['full_physical_state'] = read(latest(baseline) / 'manifest.json')['physical_state_fnv1a64'] == read(latest(split) / 'manifest.json')['physical_state_fnv1a64']
        flags['snapshot_tail_preserved'] = not sentinel.exists() and bool(list((split / 'output').glob('orphaned-after-resume-*/' + sentinel.name)))
        state['restart_checks'] = flags
        save(root / 'state.json', state)
        if not all(flags.values()): raise ValueError('split-continuous or tail rollback failed')
        for case in plan['cases'][3:5]: complete(case)
        compare('contact_regularization', 'dt_quarter', 'contact_quarter')
        compare('contact_regularization', 'contact_half', 'contact_quarter')
        for case in plan['cases']:
            if case['name'].startswith('quadrature_'): complete(case)
        if plan.get('reference_surface_contact'):
            compare('surface_quadrature', 'quadrature_64', 'dt_quarter')
        for case in plan['cases']:
            if case['name'].startswith('space_'): complete(case)
        fine = 'space_fine' if plan.get('reference_surface_contact') else 'dt_quarter'
        coarse = 'dt_quarter' if plan.get('reference_surface_contact') else 'space_coarse'
        coarse_error = force_error(results[coarse], results[fine])
        medium_error = force_error(results['space_medium'], results[fine])
        state['spatial_refinement'] = dict(coarse_to_fine=coarse_error, medium_to_fine=medium_error,
                                          physical_horizon_m=4e-5)
        if medium_error > coarse_error and medium_error > 1e-10:
            raise ValueError('spatial force error did not decrease under refinement')
        compare('space_fixed_horizon', 'space_medium', fine)
        state.update(status='small_control_checks_passed', material_calibration_allowed=False)
    except Exception as error:
        state.update(status='validation_stopped_for_review', failure=str(error), material_calibration_allowed=False)
        raise
    finally:
        state['active_case'] = None
        state['active_process'] = None
        state['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        state['artifact_index'] = {str(path.relative_to(root)): dict(size=path.stat().st_size, sha256=digest(path))
                                   for folder in [root / 'runs', root / 'restart_evidence']
                                   for path in folder.rglob('*') if path.is_file()}
        save(root / 'state.json', state)
        save(root / 'runtime_checks.json', dict(status=state['status'], source_commit=plan['source_commit'],
             plan_sha256=digest(root / 'plan.json'), state_sha256=digest(root / 'state.json'),
             comparisons=state.get('comparisons', []), restart_checks=state.get('restart_checks', {}),
             failure=state.get('failure'), scope=plan['scope'], material_calibration_allowed=False))


def review(root):
    """Separate process: re-read bytes and independently reconstruct raw fields."""
    from peri_gel_audit import audit_frames, checkpoint_physical_digest
    root = root.resolve(strict=True)
    state, plan = read(root / 'state.json'), read(root / 'plan.json')
    if state['status'] in ('prepared', 'running'):
        print(json.dumps(dict(status='waiting_for_terminal_state')))
        return
    state_sha = digest(root / 'state.json')
    output = root / 'independent_review.json'
    if output.exists():
        receipt = read(output)
        if receipt['state_sha256'] != state_sha: raise ValueError('Reviewed state changed; preserve the prior receipt')
        print(json.dumps(dict(status=receipt['status'], idempotent=True)))
        return
    receipt = dict(state_sha256=state_sha, plan_sha256=digest(root / 'plan.json'), source_commit=plan['source_commit'],
                   recorded_at_utc=datetime.now(timezone.utc).isoformat(), raw_frame_checks={}, material_calibration_allowed=False)
    try:
        if receipt['plan_sha256'] != state['plan_sha256'] or digest(root / 'solver') != plan['executable_sha256']:
            raise ValueError('plan or executable integrity failed')
        if digest(root / 'source_snapshot.tar.gz') != plan['source_snapshot_sha256']:
            raise ValueError('source snapshot integrity failed')
        with tarfile.open(root / 'source_snapshot.tar.gz', 'r:gz') as archive:
            for name, expected in plan['source_sha256'].items():
                stream = archive.extractfile(name)
                if stream is None or hashlib.sha256(stream.read()).hexdigest() != expected:
                    raise ValueError('archived source differs: ' + name)
        for name, expected in plan['support_sha256'].items():
            if digest(root / 'frozen_tools' / name) != expected: raise ValueError('review dependency changed')
        for name, expected in plan['runtime_library_sha256'].items():
            if digest(Path(name)) != expected: raise ValueError('dynamic dependency changed')
        for name, expected in state['artifact_index'].items():
            path = root / name
            if path.stat().st_size != expected['size'] or digest(path) != expected['sha256']:
                raise ValueError('raw size/SHA failed: ' + name)
        fresh = {}
        effective_inputs = {}
        for result in state['results']:
            name = result['name']
            case = next(case for case in plan['cases'] if case['name'] == name)
            config_path = root / case['configuration']
            if digest(config_path) != case['configuration_sha256']: raise ValueError('frozen case changed')
            directory = root / 'runs' / name
            validate_case(directory, parameters(config_path), plan['tolerances'])
            effective_inputs[name] = parameters(directory / 'output/effective_parameters.cfg')
            receipt['raw_frame_checks'][name] = audit_frames(directory, parameters(config_path))
            fresh[name] = rows(directory / 'output/diagnostics.csv')
        baseline_inputs = effective_inputs.get('baseline')
        if baseline_inputs:
            for name, effective in effective_inputs.items():
                allowed = ({'pd_cfl', 'spacing', 'horizon_ratio'} if name.startswith('space_') else
                           {'pd_cfl', 'surface_rings'} if name.startswith('quadrature_') else
                           {'pd_cfl', 'friction_velocity'} if name.startswith('contact_') else {'pd_cfl'})
                if effective.keys() != baseline_inputs.keys() or any(
                        effective[key] != value for key, value in baseline_inputs.items() if key not in allowed):
                    raise ValueError('uncontrolled effective parameter change: ' + name)
                if name.startswith('space_') and not math.isclose(effective['spacing'] * effective['horizon_ratio'],
                        baseline_inputs['spacing'] * baseline_inputs['horizon_ratio'], rel_tol=1e-14):
                    raise ValueError('physical horizon changed during spatial refinement')
            receipt['effective_parameter_control_checked'] = True
        receipt['force_comparisons'] = [dict(item, independently_recomputed_error=force_error(fresh[item['coarse']], fresh[item['fine']]))
                                        for item in state.get('comparisons', [])]
        if state.get('restart_checks'):
            baseline, split = root / 'runs/baseline', root / 'runs/restart_split'
            for name in ['diagnostics.csv', 'step_trace.csv']:
                if (baseline / 'output' / name).read_bytes() != (split / 'output' / name).read_bytes():
                    raise ValueError('independent restart byte comparison failed')
            receipt['raw_frame_checks']['restart_split'] = audit_frames(split, parameters(root / 'baseline.cfg'))
            def latest(directory):
                return directory / 'checkpoints' / (directory / 'checkpoints/LATEST').read_text().strip()
            digests = {}
            for name, checkpoint, output_folder in [
                ('continuous', latest(baseline), baseline / 'output'), ('split', latest(split), split / 'output'),
                ('before_noadvance', root / 'restart_evidence/before_noadvance', split / 'output'),
                ('after_noadvance', root / 'restart_evidence/after_noadvance', split / 'output')]:
                digests[name] = checkpoint_physical_digest(checkpoint / 'state.bin', output_folder)
                if digests[name]['physical_fnv1a64'] != read(checkpoint / 'manifest.json')['physical_state_fnv1a64']:
                    raise ValueError('independent physical-region FNV disagrees with native receipt')
            if digests['continuous'] != digests['split'] or digests['before_noadvance'] != digests['after_noadvance']:
                raise ValueError('independent full-state SHA/cursor restart comparison failed')
            receipt['restart_checks'] = dict(state['restart_checks'], independent_full_physical_bytes=digests)
        receipt['status'] = ('independent_small_control_checks_passed' if state['status'] == 'small_control_checks_passed'
                             else 'independent_review_completed_case_remains_unqualified')
        receipt['runtime_failure'] = state.get('failure')
        receipt['scope'] = plan['scope']
    except Exception as error:
        receipt.update(status='independent_check_failed', failure=str(error))
        raise
    finally:
        save(output, receipt)
    print(json.dumps(dict(status=receipt['status'], checked_cases=list(receipt['raw_frame_checks']),
                          runtime_failure=receipt.get('runtime_failure')), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--executable', type=Path, required=True)
    prep.add_argument('--config', type=Path, required=True)
    prep.add_argument('--host-storage', type=Path, required=True)
    prep.add_argument('--surface-contact', action='store_true', help='New independent surface quadrature and dx=10/7.5/5um spatial study')
    execute = sub.add_parser('run'); execute.add_argument('--root', type=Path, required=True)
    inspect = sub.add_parser('review'); inspect.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare': prepare(args)
    elif args.action == 'run': run(args.root)
    else: review(args.root)


if __name__ == '__main__': main()

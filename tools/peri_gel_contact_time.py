#!/usr/bin/env python3
"""Anchored adaptive-time diagnosis of an already reviewed contact trajectory.

The anchor is read, never dispatched or re-reviewed. Its exact ELF and physical
source archive are reused for new smaller-step independent initial states.
This is diagnosis only: it cannot release a failed normal-penalty study.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

from peri_gel_validation import (digest, force_error, identity, parameters, read,
                                 resources, rows, save, validate_case, verify_prerequisites)


def require_time_only_inputs(anchor, values, factor):
    if values.keys() != anchor.keys() or any(values[key] != value for key, value in anchor.items() if key not in ('pd_cfl', 'maximum_dt')) or \
            values['pd_cfl'] != anchor['pd_cfl'] * factor or values['maximum_dt'] != anchor['maximum_dt'] * factor:
        raise ValueError('Diagnosis permits only the declared CFL/cap halving; physical and contact inputs are frozen')


def prepare(args):
    anchor, root = args.anchor.resolve(strict=True), args.root.resolve()
    if root.exists(): raise FileExistsError('Diagnostic state cannot be reused')
    old_plan, old_state, old_review = [read(anchor / name) for name in ('plan.json', 'state.json', 'independent_review.json')]
    if old_review['status'] not in ('independent_small_control_checks_passed', 'independent_review_completed_case_remains_unqualified') or \
            old_review['state_sha256'] != digest(anchor / 'state.json') or old_review['plan_sha256'] != digest(anchor / 'plan.json'):
        raise ValueError('Anchor requires a matching completed mathematical/integrity review')
    if 'baseline' not in old_review['raw_frame_checks']: raise ValueError('Anchor baseline raw fields were not checked')
    host = Path(old_plan['host_storage'])
    resources(root.parent, host)
    repo = Path(__file__).resolve().parents[1]
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo).strip():
        raise ValueError('Commit the diagnostic controller before freezing it')
    dependencies = [anchor / name for name in ('plan.json', 'state.json', 'independent_review.json')]
    dependencies += [anchor / 'runs/baseline/output' / name for name in ('diagnostics.csv', 'step_trace.csv', 'effective_parameters.cfg')]
    for path in dependencies[3:]:
        expected = old_state['artifact_index'][str(path.relative_to(anchor))]
        if path.stat().st_size != expected['size'] or digest(path) != expected['sha256']:
            raise ValueError('Anchor output differs from its terminal size/SHA receipt')
    if digest(anchor / 'solver') != old_plan['executable_sha256'] or digest(anchor / 'source_snapshot.tar.gz') != old_plan['source_snapshot_sha256']:
        raise ValueError('Anchor ELF/source archive changed')
    root.mkdir()
    shutil.copy2(anchor / 'solver', root / 'solver')
    shutil.copy2(anchor / 'source_snapshot.tar.gz', root / 'source_snapshot.tar.gz')
    (root / 'frozen_tools').mkdir()
    for name in ('peri_gel_contact_time.py', 'peri_gel_validation.py', 'peri_gel_audit.py', 'peri_run.py'):
        shutil.copy2(repo / 'tools' / name, root / 'frozen_tools' / name)
    values = parameters(anchor / 'runs/baseline/output/effective_parameters.cfg')
    if not values['maximum_dt'] or not values['reference_surface_contact']:
        raise ValueError('This diagnosis requires a capped reference-surface anchor')
    cases = []
    for name, factor in [('dt_half', .5), ('dt_quarter', .25)]:
        config = dict(values, pd_cfl=values['pd_cfl'] * factor, maximum_dt=values['maximum_dt'] * factor)
        path = root / (name + '.cfg')
        path.write_text(''.join(f'{key}={value:.17g}\n' for key, value in sorted(config.items())))
        cases.append(dict(name=name, config=path.name, sha256=digest(path)))
    plan = dict(scope='Adaptive CFL/cap diagnosis at unchanged physical loading, not rigid-contact or calibration release',
        source_commit=old_plan['source_commit'], verification_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo).decode().strip(),
        source_sha256=old_plan['source_sha256'], source_snapshot_sha256=old_plan['source_snapshot_sha256'],
        executable_sha256=old_plan['executable_sha256'], runtime_library_sha256=old_plan['runtime_library_sha256'],
        prerequisite_sha256={str(path): digest(path) for path in dependencies},
        host_storage=str(host), resource_limits=old_plan['resource_limits'], cases=cases,
        support_sha256={path.name: digest(path) for path in (root / 'frozen_tools').iterdir()},
        anchor_curve=str(dependencies[3]), anchor_effective_inputs=str(dependencies[5]),
        anchor_metrics=next(item for item in old_state['results'] if item['name'] == 'baseline'),
        tolerances=dict(force_relative=.02, mass_relative=1e-10, max_penetration_over_spacing=.1),
        policy='Stop on any execution/integrity/resource/accuracy failure; no existing state redispatch or anchor rerun.',
        median_halving_claim=False, material_calibration_allowed=False, rigid_contact_qualified=False)
    save(root / 'plan.json', plan)
    save(root / 'state.json', dict(status='prepared', plan_sha256=digest(root / 'plan.json'), results=[]))
    print(json.dumps(dict(status='prepared', root=str(root), physical_source_commit=plan['source_commit'], verification_commit=plan['verification_commit'])))


def verify(root, plan, state):
    if digest(root / 'plan.json') != state['plan_sha256'] or digest(root / 'solver') != plan['executable_sha256'] or \
            digest(root / 'source_snapshot.tar.gz') != plan['source_snapshot_sha256']:
        raise ValueError('Frozen diagnostic plan/ELF/source changed')
    verify_prerequisites(plan)
    for name, expected in plan['runtime_library_sha256'].items():
        if digest(Path(name)) != expected: raise ValueError('Dynamic dependency changed')
    for name, expected in plan['support_sha256'].items():
        if digest(root / 'frozen_tools' / name) != expected: raise ValueError('Diagnostic support changed')


def run(root):
    root = root.resolve(strict=True)
    plan, state = read(root / 'plan.json'), read(root / 'state.json')
    if state['status'] != 'prepared': raise ValueError('Existing diagnostic cannot be redispatched')
    verify(root, plan, state)
    state.update(status='running', controller=identity(os.getpid()))
    save(root / 'state.json', state)
    previous, metrics = rows(Path(plan['anchor_curve'])), plan['anchor_metrics']
    try:
        for case in plan['cases']:
            verify(root, plan, state); resources(root, Path(plan['host_storage']))
            config = root / case['config']
            if digest(config) != case['sha256']: raise ValueError('Frozen diagnostic config changed')
            directory = root / 'runs' / case['name']; directory.mkdir(parents=True, exist_ok=False)
            command = [str(root / 'solver'), '--case-config', str(config), '--state_recording=true', '--pd-threads=8',
                       '--pd-checkpoint-seconds=0', '--pd-checkpoint-steps=0']
            state.update(active_case=case['name'], active_command=command)
            with (directory / 'run.log').open('w') as log:
                process = subprocess.Popen(command, cwd=directory, env=dict(os.environ, PERI_HOST_STORAGE=plan['host_storage']), stdout=log, stderr=subprocess.STDOUT)
                try:
                    snapshot = identity(process.pid)
                    for _ in range(10):
                        if snapshot['argv_observed'] or process.poll() is not None: break
                        time.sleep(.005); snapshot = identity(process.pid)
                    state['active_process'] = snapshot
                    if process.poll() is None and (snapshot['argv'] != command or snapshot['ppid'] != os.getpid() or
                            snapshot['executable'] != str(root / 'solver') or snapshot['boot_id'] != state['controller']['boot_id']):
                        raise ValueError('Actual full child identity does not match the planned diagnostic')
                    save(root / 'state.json', state)
                    try: code = process.wait(timeout=3605)
                    except subprocess.TimeoutExpired:
                        state['last_return_code'] = 124
                        raise ValueError('wall_time_budget: diagnostic exceeded 3605s')
                    state['last_return_code'] = code
                    if code != 0: raise ValueError(f"{case['name']}: solver return {code}; original logs preserved")
                finally:
                    if process.poll() is None:
                        process.terminate()
                        try: process.wait(timeout=10)
                        except subprocess.TimeoutExpired: process.kill(); process.wait()
            if sum(path.stat().st_size for path in directory.rglob('*') if path.is_file()) > plan['resource_limits']['output_bytes_per_case']:
                raise ValueError('Diagnostic output budget exceeded')
            current = validate_case(directory, parameters(config), plan['tolerances'])
            curve = rows(directory / 'output/diagnostics.csv')
            error = force_error(previous, curve)
            comparison = dict(initial_dt_ratio=current['initial_actual_dt_s'] / metrics['initial_actual_dt_s'],
                              median_dt_ratio=current['median_actual_dt_s'] / metrics['median_actual_dt_s'],
                              relative_stored_force_curve_error=error)
            state['results'].append(dict(name=case['name'], metrics=current, comparison=comparison))
            save(root / 'state.json', state)
            if not math.isclose(comparison['initial_dt_ratio'], .5, rel_tol=1e-12) or comparison['median_dt_ratio'] >= 1:
                raise ValueError('Actual initial step did not halve or actual median did not decrease')
            if error > plan['tolerances']['force_relative']:
                raise ValueError(f"{case['name']}: stored force curve error {error:.12g} exceeds predefined 0.02")
            previous, metrics = curve, current
            print(json.dumps(state['results'][-1]), flush=True)
        state['status'] = 'adaptive_time_diagnostic_checks_passed'
    except Exception as error:
        state.update(status='diagnostic_stopped_for_review', failure=str(error))
        raise
    finally:
        state.update(active_case=None, active_process=None, finished_at_utc=datetime.now(timezone.utc).isoformat())
        state['artifact_index'] = {str(path.relative_to(root)): dict(size=path.stat().st_size, sha256=digest(path))
                                   for path in (root / 'runs').rglob('*') if path.is_file()}
        save(root / 'state.json', state)


def review(root):
    from peri_gel_audit import audit_frames
    root = root.resolve(strict=True)
    plan, state = read(root / 'plan.json'), read(root / 'state.json')
    if state['status'] in ('prepared', 'running'): print(json.dumps(dict(status='waiting'))); return
    output, state_sha = root / 'independent_review.json', digest(root / 'state.json')
    if output.exists():
        receipt = read(output)
        if receipt['state_sha256'] != state_sha: raise ValueError('Previously reviewed state changed')
        print(json.dumps(dict(status=receipt['status'], idempotent=True))); return
    receipt = dict(recorded_at_utc=datetime.now(timezone.utc).isoformat(), state_sha256=state_sha,
                   plan_sha256=digest(root / 'plan.json'), raw_frame_checks={}, material_calibration_allowed=False, rigid_contact_qualified=False)
    try:
        verify(root, plan, state)
        with tarfile.open(root / 'source_snapshot.tar.gz') as archive:
            for name, expected in plan['source_sha256'].items():
                stream = archive.extractfile(name)
                if stream is None or hashlib.sha256(stream.read()).hexdigest() != expected: raise ValueError('Physical source archive mismatch')
        for name, expected in state['artifact_index'].items():
            path = root / name
            if path.stat().st_size != expected['size'] or digest(path) != expected['sha256']: raise ValueError('Raw diagnostic size/SHA mismatch')
        previous = rows(Path(plan['anchor_curve']))
        previous_metrics = plan['anchor_metrics']
        anchor = parameters(Path(plan['anchor_effective_inputs']))
        receipt['comparisons'] = []
        for item in state['results']:
            case = next(case for case in plan['cases'] if case['name'] == item['name'])
            config = root / case['config']
            if digest(config) != case['sha256']: raise ValueError('Diagnostic case changed')
            values = parameters(config)
            require_time_only_inputs(anchor, values, .5 if item['name'] == 'dt_half' else .25)
            directory = root / 'runs' / item['name']
            metrics = validate_case(directory, values, plan['tolerances'])
            receipt['raw_frame_checks'][item['name']] = audit_frames(directory, values)
            curve = rows(directory / 'output/diagnostics.csv')
            error = force_error(previous, curve)
            ratios = dict(initial=metrics['initial_actual_dt_s'] / previous_metrics['initial_actual_dt_s'],
                          median=metrics['median_actual_dt_s'] / previous_metrics['median_actual_dt_s'])
            if not math.isclose(ratios['initial'], .5, rel_tol=1e-12) or ratios['median'] >= 1:
                raise ValueError('Independent actual-step refinement check failed')
            if state['status'] == 'adaptive_time_diagnostic_checks_passed' and error > plan['tolerances']['force_relative']:
                raise ValueError('Independent force accuracy gate disagrees with runtime pass')
            receipt['comparisons'].append(dict(name=item['name'], independently_recomputed_force_error=error, actual_dt_ratios=ratios, metrics=metrics))
            previous, previous_metrics = curve, metrics
        receipt.update(status='independent_adaptive_time_diagnostic_passed' if state['status'] == 'adaptive_time_diagnostic_checks_passed' else
                       'independent_diagnostic_review_completed_remains_unqualified', runtime_failure=state.get('failure'),
                       scope=plan['scope'], sparse_raw_not_full_step_energy=True, anchor_redispatched=False, anchor_re_reviewed=False)
    except Exception as error:
        receipt.update(status='independent_diagnostic_check_failed', failure=str(error)); raise
    finally: save(output, receipt)
    print(json.dumps(dict(status=receipt['status'], comparisons=receipt.get('comparisons'), runtime_failure=receipt.get('runtime_failure')), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare'); prep.add_argument('--root', type=Path, required=True); prep.add_argument('--anchor', type=Path, required=True)
    for name in ('run', 'review'): sub.add_parser(name).add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare': prepare(args)
    elif args.action == 'run': run(args.root)
    else: review(args.root)


if __name__ == '__main__': main()

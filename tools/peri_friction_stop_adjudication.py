"""Read-only adjudication of a stopped actual-PD campaign, never redispatch.

Equation consistency and required coverage have separate verdicts. A native
return zero cannot turn missing same-contact temporal reversal into a pass.
The original terminal, case configuration and raw bytes remain immutable.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path
import resource
import sys
import time
import traceback


def contact_events(header, row):
    """Independent binary-input rational gather with the frozen thresholds."""
    f = Fraction
    velocity = list(map(f, row['solve']['velocity']))
    matrix = header['C_full_three_components']
    threshold = f(header['normalization_impulse_scale']) * f(1e-8)
    result = []
    for s, contact in enumerate(row['contacts']):
        normal = f(row['solve']['impulses'][3*s])
        slip = [sum((f(a)*b for a, b in zip(matrix[3*s+k+1], velocity)), f(0))
                - f(contact['wall_tangent_velocity'][k]) for k in range(2)]
        active = normal > threshold
        moving = sum((x*x for x in slip), f(0)) > f(Decimal('1e-8'))**2
        signs = [1 if x > f(1e-8) else (-1 if x < -f(1e-8) else 0) for x in slip]
        result.append(dict(row=s, plate=header['plate_ids'][s], active=active,
                           sliding=active and moving, signs=signs,
                           relative_velocity=list(map(float, slip)), normal_impulse=float(normal)))
    return result


def summarize(events):
    signs = {}
    phases = {}
    for step in events:
        phase = phases.setdefault(step['phase'], dict(steps=0, active=0, sliding=0,
             upper_wall_velocity=[], active_rows=[], sliding_rows=[]))
        phase['steps'] += 1
        phase['upper_wall_velocity'].append(step['upper_wall_velocity'])
        for event in step['contacts']:
            if event['active']:
                phase['active'] += 1
                phase['active_rows'].append([step['step'], event['row']])
            if event['sliding']:
                phase['sliding'] += 1
                phase['sliding_rows'].append([step['step'], event['row'], event['signs']])
                for a, sign in enumerate(event['signs']):
                    if sign:
                        signs.setdefault((event['row'], a), set()).add(sign)
    reverse = [sum(value == {-1, 1} for (s, k), value in signs.items() if k == a) for a in range(2)]
    return dict(phases=phases, same_contact_temporal_reversals=reverse,
                per_row_signs=[dict(row=s, tangent=a, signs=sorted(value))
                               for (s, a), value in sorted(signs.items())])


def run(root, output):
    # Import only the prior frozen reference, not a solver or modified reader.
    sys.path.insert(0, str(root/'frozen'))
    from peri_bounded_friction_audit import FrictionInspector, reversed_contact_counts
    from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
    assert not output.exists(), 'a read-only adjudication is attempted only once'
    plan_path = output.parent/'adjudication_plan.json'
    plan = read(plan_path)
    verify(plan['artifacts'])
    assert record(Path(__file__).resolve()) == plan['adjudicator']
    output.parent.mkdir(exist_ok=True)
    state_path = output.parent/'adjudication_state.json'
    save(state_path, dict(status='attempted_read_only', identity=identity(), started_at=stamp()))
    resource.setrlimit(resource.RLIMIT_AS, (8*2**30, 8*2**30))
    resource.setrlimit(resource.RLIMIT_CPU, (14400, 14400))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    start = time.monotonic()
    receipt = dict(status='read_only_stop_adjudication_failed', identity=identity(),
                   solver_calls=0, old_reviews_repeated=False, study_passed=False,
                   S3_complete=False, full_restart_qualified=False, calibration_allowed=False)
    try:
        prep = read(root/'preparation.json')
        verify(prep['artifacts']+prep['inputs'])
        terminal = read(root/'terminal.json')
        assert terminal['status'] == 'stopped_for_review'
        assert terminal['failure'] == "AssertionError(('actual_PD_hold_release_reverse', 'same_contact_temporal_slip_reversal'))"
        assert len(terminal['attempts']) == 3
        case = next(c for c in read(root/'plan.json')['controls'] if c['name'] == 'actual_PD_hold_release_reverse')
        attempt = terminal['attempts'][-1]
        verify([attempt['stdout'], attempt['stderr']])
        assert attempt['return_code'] == 0 and attempt['stop_reason'] is None
        assert not (root/(case['name']+'.independent_review.json')).exists()
        assert attempt['cursor']['durability_confirmed'] and attempt['cursor']['last_durable_accepted_steps'] == 40
        reserves(output.parent)
        events = []
        with localcontext() as context:
            context.prec = 95
            with open(attempt['stdout']['path'], encoding='utf8') as stream:
                meta = json.loads(next(stream)); header = json.loads(next(stream))
                assert meta == attempt['native_identity'] and meta['argv'] == attempt['command']
                assert meta['pid'] == attempt['child_pid'] and meta['ppid'] == attempt['identity']['pid']
                assert meta['boot_id'] == attempt['identity']['boot_id']
                audit = FrictionInspector(case, header, True)
                old, z, loss, history = audit.initial()
                count = 0
                for line in stream:
                    row = json.loads(line)
                    if row['type'] == 'end':
                        end = row
                        break
                    events.append(dict(step=row['state']['step'], time=row['state']['time'],
                         phase=row['interval']['phase'],
                         upper_wall_velocity=row['interval']['wall_tangent_velocity'][1],
                         contacts=contact_events(header, row)))
                    old, z, loss, history = audit.step(row, old, z, loss, history)
                    count += 1
                    print(json.dumps(dict(checked_steps=count, points=audit.trial_points)), flush=True)
                assert not any(line.strip() for line in stream)
                assert end['completed'] and end['accepted_steps'] == count == case['steps'] == 40
                assert end['last_accepted_time'] == old['time'] == header['end_time']
            coverage = summarize(events)
            assert coverage['same_contact_temporal_reversals'] == reversed_contact_counts(audit.row_slip_signs)
            for field in case['require_statistics']:
                assert audit.statistics[field] > 0, ('other_required_coverage_failed', field)
            # Preserve, rather than bypass, the original failing postcondition.
            assert case['require_two_way_slip'] and not all(n > 0 for n in coverage['same_contact_temporal_reversals'])
            receipt.update(status='original_equations_independently_consistent_required_reversal_coverage_FAILED',
                original_failure=terminal['failure'], required_reversal_coverage='FAILED',
                accepted_steps=count, all_trial_point_tangents_checked=True,
                independent_trial_material_points=audit.trial_points,
                maximum_normalized_error=max(audit.maximum, audit.geometry['maximum_normalized_cache_error']),
                coverage=coverage, observations={k:float(v) if isinstance(v,Decimal) else v for k,v in audit.statistics.items()},
                raw=attempt['stdout'], native_attempt=record(root/(case['name']+'.attempt.json')))
        save(output.parent/'coverage_events.json', events)
        receipt['coverage_events'] = record(output.parent/'coverage_events.json')
        verify(plan['artifacts']); verify(prep['artifacts']+prep['inputs'])
        verify([attempt['stdout'], attempt['stderr']])
    except Exception as error:
        receipt.update(failure=repr(error), traceback=traceback.format_exc())
    receipt.update(finished_at=stamp(), wall_seconds=time.monotonic()-start,
                   max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   plan=record(plan_path), state=record(state_path))
    save(output, receipt)
    print(json.dumps(dict(status=receipt['status'], failure=receipt.get('failure'))), flush=True)
    return int('failure' in receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.exit(run(args.root, args.output))

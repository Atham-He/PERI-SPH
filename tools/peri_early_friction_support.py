"""A new coverage input plus two never-attempted tails, not an old retry."""
from pathlib import Path
from peri_bounded_friction_support import controls, encode
from peri_elastic_normal import read, record, verify


def continuation_controls(attempted_names):
    prior = controls()
    assert attempted_names == [c['name'] for c in prior[:3]], 'unexpected previous attempts'
    cycle = dict(prior[2], name='actual_PD_early_biaxial_slip_cycle', protocol_kind=2,
                 tangent_stroke=[.0005, .0003333333333333333])
    return [cycle, prior[3], prior[4]]


def qualification_prerequisite(prior_root, adjudication_root):
    prior_root, adjudication_root = Path(prior_root), Path(adjudication_root)
    stop = read(prior_root/'terminal.json')
    assert stop['status'] == 'stopped_for_review'
    assert stop['failure'] == "AssertionError(('actual_PD_hold_release_reverse', 'same_contact_temporal_slip_reversal'))"
    assert len(stop['attempts']) == 3
    case_list = continuation_controls([a['name'] for a in stop['attempts']])
    adjudication = read(adjudication_root/'independent_stop_adjudication.json')
    adjudication_plan = read(adjudication_root/'adjudication_plan.json')
    verify(adjudication_plan['artifacts'])
    verify([adjudication['plan'],adjudication['state'],adjudication['coverage_events']])
    assert adjudication['status'] == 'original_equations_independently_consistent_required_reversal_coverage_FAILED'
    assert adjudication['required_reversal_coverage'] == 'FAILED' and not adjudication['study_passed']
    assert adjudication['all_trial_point_tangents_checked'] and adjudication['accepted_steps'] == 40
    assert adjudication['solver_calls'] == 0 and not adjudication['old_reviews_repeated']
    artifacts = [record(prior_root/name) for name in ['preparation.json','state.json','terminal.json','plan.json']]
    state = read(prior_root/'state.json')
    verify([state['preparation']])
    prep = read(prior_root/'preparation.json')
    # Historical source has legitimately advanced. Its immutable frozen files,
    # rather than the current working source, are the old qualification basis.
    snapshots = [e for e in prep['artifacts'] if e['path'].startswith(str(prior_root/'frozen')+'/')]
    assert snapshots
    verify(snapshots)
    artifacts.extend(snapshots)
    for a in stop['attempts']:
        verify([a['stdout'],a['stderr']]);artifacts.extend([a['stdout'],a['stderr']])
    for a in stop['attempts'][:2]:
        verify([a['independent_review']])
        p = prior_root/(a['name']+'.independent_review.json');receipt = read(p)
        assert receipt['status'] == 'actual_PD_friction_case_independently_passed'
        artifacts.append(record(p))
    assert adjudication['raw'] == stop['attempts'][-1]['stdout']
    artifacts.extend([record(adjudication_root/name) for name in
                      ['adjudication_plan.json','adjudication_state.json','independent_stop_adjudication.json','coverage_events.json']])
    for case in case_list[1:]:
        name = case['name']
        assert not any((prior_root/(name+suffix)).exists() for suffix in
                       ['.attempt.json','.jsonl','.review_state.json','.independent_review.json'])
        assert (prior_root/(name+'.input')).read_text() == encode(case)
        artifacts.append(record(prior_root/(name+'.input')))
    assert encode(case_list[0]) not in [(prior_root/(a['name']+'.input')).read_text() for a in stop['attempts']]
    verify(artifacts)
    return dict(artifacts=artifacts,qualified_previous_cases=2,old_study_passed=False,
                adjudicated_missing_coverage=True,old_solver_calls_repeated=0,
                old_case_reviews_repeated=0,continuation_controls=case_list)

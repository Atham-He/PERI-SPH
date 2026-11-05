"""Original frozen inspect_stream finalization body,AST identical."""
from decimal import Decimal as D
from peri_bounded_friction_audit import reversed_contact_counts,stability_spectrum

def original_summary(audit, case, h, count, full):
    stat = audit.statistics
    for field in case.get('require_statistics', []):
        assert stat[field] > 0, (case['name'], 'missing_observation', field)
    if case.get('require_two_way_slip'):
        assert all((n > 0 for n in reversed_contact_counts(audit.row_slip_signs))), (case['name'], 'same_contact_temporal_slip_reversal')
    spectrum = stability_spectrum(case['rotation'], tuple(h['reference']['volume']), h['reference']['horizon'], *[float(h['parameters'][k]) for k in [0, 1]], float(h['eta'])) if full else None
    return dict(name=case['name'], status='independent_actual_PD_original_friction_consistency_passed' if full else 'runtime_pending_full_point_tangent_review', accepted_steps=count, maximum_normalized_error=max(audit.maximum, audit.geometry['maximum_normalized_cache_error']), geometry=audit.geometry, independent_stability=spectrum, all_trial_point_tangents_checked=full, independent_trial_material_points=audit.trial_points, observations={k: float(v) if isinstance(v, D) else v for k, v in stat.items()}, temporally_reversed_contacts_per_tangent=reversed_contact_counts(audit.row_slip_signs), algorithmic_loss_positive_input=float(stat['algorithmic_loss'] / stat['positive_input_work']) if stat['positive_input_work'] else None, normal_only=False, trajectory_convergence_qualified=False, angular_conservation_qualified=False, full_restart_qualified=False, S3_complete=False, calibration_allowed=False)

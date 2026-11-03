#!/usr/bin/env python3
"""Prepare a traceable calibration sequence; never dispatch or fit a simulation.

Contact and material parameters are separate blocks. The experimental row order
does not provide the missing time axis or individual specimen geometry.
"""
import argparse
import json
from pathlib import Path

from peri_run import digest


PARAMETER_BLOCKS = {
    'contact': frozenset({'friction_mu', 'friction_velocity'}),
    'material': frozenset({'youngs_modulus', 'poisson_ratio', 'yield_stress',
                           'hardening_modulus', 'plastic_viscosity'}),
}


def parameter_block(varying_parameters):
    """Reject simultaneous contact/material compensation and untracked knobs."""
    varying = set(varying_parameters)
    if not varying:
        raise ValueError('At least one explicit varying parameter is required')
    for name, allowed in PARAMETER_BLOCKS.items():
        if varying <= allowed:
            return name
    raise ValueError('Vary one contact or material block at a time; force scaling and hidden knobs are forbidden')


def prepare_plan(manifest_path, target='Micro-H_0TO-3SA-H'):
    manifest_path = Path(manifest_path).resolve(strict=True)
    source = json.loads(manifest_path.read_text(encoding='utf-8'))
    curves = source['curves']
    selected = [curve for curve in curves if curve['sheet'] + '_' + curve['sample_label'] == target]
    if len(selected) != 1:
        raise ValueError('Target must identify exactly one exported curve')
    curve = selected[0]
    if curve['experiment_kind'] != 'compression' or not curve['source_order_preserved']:
        raise ValueError('Baseline requires a compression curve in original source order')
    if curve.get('nominal_tung_oil_loading_percent') != 0 or curve.get('hydration_state') != 'hydrated':
        raise ValueError('This sequence starts with a zero-oil hydrated baseline; oil-bearing fits come later')
    csv_path = (manifest_path.parent / curve['csv']).resolve(strict=True)
    if not csv_path.is_relative_to(manifest_path.parent) or digest(csv_path) != curve['csv_sha256']:
        raise ValueError('Experimental CSV path or SHA does not match its source manifest')
    duplicate_groups = [group for group in source['duplicate_groups'] if target in group]
    return {
        'schema_version': 1,
        'status': 'development_sequence_only_no_fit_or_dispatch',
        'preparation_script_sha256': digest(Path(__file__)),
        'target': target,
        'experimental_source': {
            'manifest': str(manifest_path), 'manifest_sha256': digest(manifest_path),
            'workbook_sha256': source['source_sha256'], 'csv': str(csv_path),
            'csv_sha256': curve['csv_sha256'], 'rows': curve['rows'],
            'duplicate_groups': duplicate_groups,
            'independent_replicate_identity': curve['replicate_identity'],
        },
        'unknowns': {
            'specimen_diameter_m': curve['specimen_diameter_m'],
            'micro_loading_rate_m_s': curve['loading_rate_m_s'],
            'physical_time_available': curve['physical_time_available'],
            'contact_zero_verified': curve['contact_zero_verified'],
        },
        'stages': [
            {'name': 'code_and_numerical_verification',
             'requires': ['new-version small-control regression', 'contact force/potential/reaction checks',
                          'read-only observation checks', 'restart no-advance/rollback/split verification'],
             'released': False},
            {'name': 'zero_oil_gel_geometry_and_contact',
             'requires': ['solid effective-gel geometry justified', 'individual geometry or explicit geometry uncertainty',
                          'fixed material during contact sensitivity', 'spacing and friction regularization sensitivity'],
             'released': False},
            {'name': 'zero_oil_material_calibration',
             'requires': ['numerically qualified baseline', 'fixed qualified contact settings',
                          'validated gel constitutive model', 'metadata sufficient for the fitted parameters',
                          'repeat contact penetration/convergence checks over candidate E/nu because normal penalty depends on material moduli'],
             'rate_parameter_fit_allowed': bool(curve['physical_time_available'] and curve['loading_rate_m_s'] is not None),
             'released': False},
            {'name': 'contact_refinement_with_frozen_baseline_material',
             'requires': ['freeze baseline material before changing contact',
                          'repeat held-out numerical/contact checks before another material fit'],
             'released': False},
            {'name': 'oil_bearing_particle_calibration',
             'requires': ['freeze qualified zero-oil baseline', 'justify chamber geometry or homogenization',
                          'independent oil amount/transport assumptions', 'coupled mass/impulse/work validation',
                          'qualified material and full-state restart for the added histories'],
             'released': False},
        ],
        'parameter_blocks': {name: sorted(values) for name, values in PARAMETER_BLOCKS.items()},
        'comparison_policy': {
            'force_unit': 'N', 'displacement_unit': 'm', 'force_scale_fit_allowed': False,
            'common_supported_stroke_only': True, 'extrapolation_allowed': False,
            'material_contact_joint_fit_allowed': False,
            'group_D50_is_individual_diameter': False,
            'source_rows_are_physical_timestamps': False,
            'equal_labelled_curves_count_as_independent_validation': False,
            'normal_penalty_is_independent_of_material': False,
            'friction_velocity_is_numerical_regularization': True,
        },
        'solver_dispatched': False, 'material_parameters_fitted': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--target', default='Micro-H_0TO-3SA-H')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = prepare_plan(args.manifest, args.target)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(plan, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': plan['status'], 'target': plan['target'], 'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()

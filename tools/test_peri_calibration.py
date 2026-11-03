import json
from pathlib import Path
import tempfile
import unittest

from peri_calibration import parameter_block, prepare_plan
from peri_run import digest


class CalibrationProtocolTests(unittest.TestCase):
    def test_contact_material_compensation_and_force_rescaling_are_rejected(self):
        self.assertEqual(parameter_block(['friction_mu', 'friction_velocity']), 'contact')
        self.assertEqual(parameter_block(['youngs_modulus', 'poisson_ratio']), 'material')
        for values in [[], ['friction_mu', 'youngs_modulus'], ['force_scale'], ['contact_stiffness_factor']]:
            with self.assertRaises(ValueError):
                parameter_block(values)

    def test_preparation_preserves_unknowns_and_cannot_release_fitting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            curve_file = root / 'curve.csv'
            curve_file.write_text('source_row,displacement_m,force_N,time_s\n4,0,0,\n5,1e-6,1e-3,\n')
            curve = dict(sheet='Micro-H', sample_label='0TO-3SA-H', csv=curve_file.name,
                         csv_sha256=digest(curve_file), experiment_kind='compression',
                         source_order_preserved=True, replicate_identity=None, rows=2,
                         specimen_diameter_m=None, loading_rate_m_s=None,
                         nominal_tung_oil_loading_percent=0, hydration_state='hydrated',
                         physical_time_available=False, contact_zero_verified=False)
            manifest = root / 'manifest.json'
            manifest.write_text(json.dumps(dict(curves=[curve], duplicate_groups=[], source_sha256='workbook-hash')))
            plan = prepare_plan(manifest)
            self.assertIsNone(plan['unknowns']['specimen_diameter_m'])
            self.assertFalse(plan['unknowns']['physical_time_available'])
            self.assertFalse(any(stage['released'] for stage in plan['stages']))
            self.assertFalse(plan['stages'][2]['rate_parameter_fit_allowed'])
            self.assertFalse(plan['solver_dispatched'])
            self.assertFalse(plan['material_parameters_fitted'])
            curve_file.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'SHA'):
                prepare_plan(manifest)

    def test_nanoindentation_and_duplicate_target_records_are_not_silent_micro_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / 'manifest.json'
            target = dict(sheet='Micro-H', sample_label='0TO-3SA-H',
                          experiment_kind='nanoindentation', source_order_preserved=True)
            manifest.write_text(json.dumps(dict(curves=[target], duplicate_groups=[])))
            with self.assertRaisesRegex(ValueError, 'compression curve'):
                prepare_plan(manifest)
            manifest.write_text(json.dumps(dict(curves=[target, target], duplicate_groups=[])))
            with self.assertRaisesRegex(ValueError, 'exactly one'):
                prepare_plan(manifest)

    def test_oil_bearing_data_cannot_skip_the_zero_oil_baseline(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / 'manifest.json'
            curve = dict(sheet='Micro-H', sample_label='20TO-3SA-H', experiment_kind='compression',
                         source_order_preserved=True, nominal_tung_oil_loading_percent=20,
                         hydration_state='hydrated')
            manifest.write_text(json.dumps(dict(curves=[curve], duplicate_groups=[])))
            with self.assertRaisesRegex(ValueError, 'zero-oil hydrated baseline'):
                prepare_plan(manifest, 'Micro-H_20TO-3SA-H')


if __name__ == '__main__':
    unittest.main()

import csv
import math
from pathlib import Path
import struct
import tempfile
import unittest
from peri_gel_validation import force_error, validate_case
from peri_gel_audit import material_response, matrix, checkpoint_physical_digest, surface_reference, surface_contact_fields
from peri_gel_contact_time import require_time_only_inputs


class GelValidationTests(unittest.TestCase):
    def test_contact_diagnosis_cannot_change_speed_mass_material_or_penalty(self):
        anchor = dict(pd_cfl=.1, maximum_dt=4e-9, compression_speed=.1, density=1000, youngs_modulus=1e6, contact_penalty_length=1e-9)
        fine = dict(anchor, pd_cfl=.05, maximum_dt=2e-9)
        require_time_only_inputs(anchor, fine, .5)
        for key in ('compression_speed', 'density', 'youngs_modulus', 'contact_penalty_length'):
            changed = dict(fine); changed[key] *= .5
            with self.assertRaises(ValueError): require_time_only_inputs(anchor, changed, .5)
        with self.assertRaises(ValueError): require_time_only_inputs(anchor, dict(fine, maximum_dt=4e-9), .5)

    def test_forged_step_cap_and_two_wall_penetration_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); output = root / 'output'; output.mkdir()
            config = dict(compression_stroke=1, compression_speed=1, relaxation_time=100, maximum_dt=1, spacing=1)
            (output / 'effective_parameters.cfg').write_text(''.join(f'{key}={value}\n' for key, value in config.items()))
            (output / 'run_status.txt').write_text('completed\n')
            data = dict(time_s=1, stroke_m=1, pd_steps=1, mass_relative_error=0, particle_mass_error=0,
                upper_penetration_over_spacing=0, upper_max_penetration_m=0, lower_max_penetration_m=0,
                total_penetration_over_total_stroke=0, material_energy_available=1, gel_dissipation_J=0, mechanical_residual_J=0)
            step = dict(time_s=1, dt_s=1, pd_steps=1, last_contact_sample_time_s=.5, material_history_time_s=1,
                pd_wave_dt_s=2, configured_candidate_dt_s=1, output_limited_dt_s=1, contact_limited_dt_s=1)
            def write():
                for name, row in [('diagnostics.csv', data), ('step_trace.csv', step)]:
                    with (output/name).open('w', newline='') as stream:
                        writer = csv.DictWriter(stream, row.keys()); writer.writeheader(); writer.writerow(row)
            tolerance = dict(mass_relative=1e-10, max_penetration_over_spacing=.1)
            write(); self.assertEqual(validate_case(root, config, tolerance)['accepted_steps'], 1)
            step['configured_candidate_dt_s'] = 2; write()
            with self.assertRaisesRegex(ValueError, 'bounds'): validate_case(root, config, tolerance)
            step['configured_candidate_dt_s'] = 1
            data.update(upper_max_penetration_m=.05, lower_max_penetration_m=.05); write()
            with self.assertRaisesRegex(ValueError, 'two-wall'): validate_case(root, config, tolerance)

    def test_surface_cache_rejects_clipped_partition_and_duplicate_targets(self):
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / 'surface.csv'
            header = 'sample,x,y,z,area,particle,weight\n'
            file.write_text(header + '0,0,0,0,1,0,1.2\n0,0,0,0,1,1,-0.2\n')
            self.assertEqual(len(surface_reference(file)), 1)
            file.write_text(header + '0,0,0,0,1,0,1.2\n')
            with self.assertRaises(ValueError): surface_reference(file)
            file.write_text(header + '0,0,0,0,1,0,0.5\n0,0,0,0,1,0,0.5\n')
            with self.assertRaises(ValueError): surface_reference(file)

    def test_missing_full_surface_rule_is_not_a_zero_force_pass(self):
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            surface_contact_fields({'Mass': [1]}, {}, dict(compression_stroke=1, compression_speed=1,
                youngs_modulus=1, poisson_ratio=.3, relaxing_shear_ratio=1, contact_penalty_length=1,
                surface_rings=4, surface_azimuths=8, radius=1), [dict(area=1)])

    def test_checkpoint_hash_excludes_wall_time_but_includes_internal_history_and_checks_prefix(self):
        def number(value): return struct.pack('<Q', value)
        def label(value): return number(len(value)) + value.encode()
        def checkpoint(wall, history):
            header = label('PERI-SPH-CPU-CHECKPOINT-1') + b''.join(number(v) for v in [0x0102030405060708, 8, 8, 3])
            header += label('immutable-identity') + struct.pack('<ddQQQ', 1., wall, 5, 5, 2)
            physical = number(1) + b'GelViscousStrain' + struct.pack('<d', history)
            attachments = number(3)
            for name in ['effective_parameters.cfg', 'diagnostics.csv', 'step_trace.csv']:
                payload = name.encode()
                attachments += label(name) + number(len(payload)) + payload
            payload = header + physical + attachments
            checksum = 14695981039346656037
            for byte in payload: checksum = ((checksum ^ byte) * 1099511628211) & ((1 << 64) - 1)
            return payload + struct.pack('<QQ', checksum, 0x50455249444F4E45)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['effective_parameters.cfg', 'diagnostics.csv', 'step_trace.csv']:
                (root/name).write_bytes(name.encode())
            file = root/'state.bin'
            file.write_bytes(checkpoint(17, 1.2)); a = checkpoint_physical_digest(file, root)
            file.write_bytes(checkpoint(230, 1.2)); b = checkpoint_physical_digest(file, root)
            self.assertEqual(a, b)
            file.write_bytes(checkpoint(230, 1.3)); c = checkpoint_physical_digest(file, root)
            self.assertNotEqual(a['physical_sha256'], c['physical_sha256'])
            (root/'diagnostics.csv').write_bytes(b'changed-prefix')
            with self.assertRaises(ValueError): checkpoint_physical_digest(file, root)

    def test_known_force_error_and_no_extrapolation_or_zero_force_pass(self):
        baseline = [dict(stroke_m=x, upper_observed_Fy_N=2*x) for x in [0., .2, .7, 1.]]
        scaled = [dict(stroke_m=x, upper_observed_Fy_N=2.02*x) for x in [0., .5, 1.]]
        self.assertAlmostEqual(force_error(scaled, baseline), .01)
        with self.assertRaises(ValueError): force_error(scaled[:-1], baseline)
        zeros = [dict(stroke_m=x, upper_observed_Fy_N=0) for x in [0., 1.]]
        with self.assertRaisesRegex(ValueError, 'UNRESOLVED'): force_error(zeros, zeros)

    def test_independent_raw_matrix_order_preserves_nonsymmetric_shear(self):
        self.assertEqual(matrix([1, 0, 0, .2, 1, 0, 0, 0, 1]), [[1, .2, 0], [0, 1, 0], [0, 0, 1]])

    def test_independent_material_energy_in_simple_shear_and_volumetric_compression(self):
        config = dict(youngs_modulus=1e6, poisson_ratio=.35, relaxing_shear_ratio=1)
        shear = config['youngs_modulus'] / (2 * (1 + config['poisson_ratio']))
        gamma = .2
        F = [[1, gamma, 0], [0, 1, 0], [0, 0, 1]]
        zero = [[0.]*3 for _ in range(3)]
        energy, _, green = material_response(F, zero, config)
        expected = shear / 2 * gamma**2 + shear * math.fsum(v*v for row in green for v in row)
        self.assertAlmostEqual(energy, expected, places=8)
        scale = .9
        F = [[scale if i == j else 0 for j in range(3)] for i in range(3)]
        J = scale**3; bulk = 1e6 / (3 * (1 - .7))
        energy, _, _ = material_response(F, zero, config)
        self.assertAlmostEqual(energy, bulk / 4 * (J*J - 1 - 2*math.log(J)), places=8)


if __name__ == '__main__': unittest.main()

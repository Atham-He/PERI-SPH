from pathlib import Path
import tempfile
import csv
import unittest
from unittest.mock import patch
from xml.parsers import expat

from peri_audit import rows, scan_vtp, inspect, acoustic_work_assessment
from peri_run import require_completion_audit


VTP = '''<VTKFile><PolyData><FieldData><DataArray Name="TimeValue" type="Float64" NumberOfTuples="1" format="ascii">1.234567890123456e-11</DataArray></FieldData>
<Piece NumberOfPoints="1"><Points><DataArray Name="Position" type="Float64" NumberOfComponents="3" format="ascii">1e-9 2e-9 3e-9</DataArray></Points><PointData>
<DataArray Name="Mass" type="Float64" format="ascii">9.350000000000002e-13</DataArray>
<DataArray Name="Density" type="Float64" format="ascii">935</DataArray>
<DataArray Name="VolumetricMeasure" type="Float64" format="ascii">1.0000000000000003e-15</DataArray>
<DataArray Name="Velocity" type="Float64" NumberOfComponents="3" format="ascii">0.01 0.02 0.03</DataArray>
</PointData></Piece></PolyData></VTKFile>'''


class SnapshotAuditTests(unittest.TestCase):
    def scan(self, xml):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"LiquidCore_ite_0000000000.vtp"
            path.write_text(xml)
            return scan_vtp(path)

    def test_completion_gate_rejects_partial_warned_and_failed_audits(self):
        scenarios = [
            {"status": "sampled_checks_passed", "simulation_completed": False},
            {"status": "needs_attention", "simulation_completed": True,
             "warnings": ["negative viscous loss"]},
            {"status": "failed", "simulation_completed": True, "reason": "bad snapshot"},
        ]
        with tempfile.TemporaryDirectory() as directory:
            for audit in scenarios:
                with self.subTest(audit=audit), patch("peri_audit.inspect", return_value=audit):
                    with self.assertRaisesRegex(ValueError, "completion audit rejected"):
                        require_completion_audit(Path(directory))

    def test_completion_gate_accepts_and_records_full_healthy_audit(self):
        audit = {"status": "sampled_checks_passed", "simulation_completed": True, "warnings": []}
        with tempfile.TemporaryDirectory() as directory:
            with patch("peri_audit.inspect", return_value=audit) as inspect:
                self.assertEqual(require_completion_audit(Path(directory), snapshots=True), audit)
                inspect.assert_called_once_with(Path(directory), snapshots=True)
            self.assertTrue((Path(directory)/"completion_audit.json").exists())

    def test_small_mass_survives_and_density_volume_reconciles(self):
        result = self.scan(VTP)
        self.assertEqual(result["points"], 1)
        self.assertEqual(result["mass_min_kg"], 9.350000000000002e-13)
        self.assertLess(result["density_volume_mass_relative_error"], 1e-14)

    def test_rounded_to_zero_particle_mass_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-positive particle Mass"):
            self.scan(VTP.replace("9.350000000000002e-13", "0.000000000"))

    def test_nonfinite_velocity_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-finite Velocity"):
            self.scan(VTP.replace("0.01 0.02 0.03", "0.01 NaN 0.03"))

    def test_declared_float32_is_rejected_in_precision_mode(self):
        with self.assertRaisesRegex(ValueError, "not Float64"):
            self.scan(VTP.replace('Name="Mass" type="Float64"', 'Name="Mass" type="Float32"'))

    def test_wrong_array_length_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "wrong tuple count"):
            self.scan(VTP.replace("0.01 0.02 0.03", "0.01 0.02"))

    def test_committed_nonfinite_csv_row_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/"progress.csv"
            p.write_text("time_s,value\n0,1\n1,NaN\n")
            with self.assertRaisesRegex(ValueError,"non-finite committed row"):
                rows(p)
            p.write_text("time_s,value\n0,1\n1,Na")
            self.assertEqual(rows(p),[{"time_s":0.,"value":1.}])


class AcousticBudgetTests(unittest.TestCase):
    def audit(self, **overrides):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory); output=run/"output"; output.mkdir()
            (output/"effective_parameters.cfg").write_text("compression_distance=1\ncompression_speed=1\nmax_density_error=0.1\nmax_fluid_mach=0.3\nc_f=10\nwith_liquid=1\n")
            (output/"run_status.txt").write_text("completed\n")
            energy={"time_s":1.,"kinetic_J":1.,"gravity_potential_change_J":0.,
                    "plate_contact_work_J":0.,"plate_fluid_pressure_work_J":0.,"plate_fluid_viscous_work_J":0.,
                    "pd_damping_signed_loss_J":0.,"acoustic_inner_signed_loss_J":0.,"acoustic_wall_signed_loss_J":0.,
                    "viscous_signed_loss_J":0.,"friction_signed_loss_J":0.,"hourglass_signed_transfer_J":0.,
                    "mechanical_residual_J":0.,"mechanical_residual_relative":0.,
                    "shell_kinetic_identity_residual_J":0.,"fluid_kinetic_identity_residual_J":0.}
            energy.update(overrides)
            diagnostic={"time_s":1.,"numerical_failures":0.,"damage_mean":0.,"damage_max":0.,
                        "fluid_density_error_max":0.,"fluid_speed_max_m_s":0.1}
            mass={"time_s":1.,"shell_count":1.,"liquid_count":1.,"shell_mass_kg":1.,"liquid_mass_kg":1.,
                  "shell_mass_relative_error":0.,"liquid_mass_relative_error":0.,"max_particle_mass_relative_error":0.,
                  "liquid_mass_outside_domain_kg":0.,"density_volume_mass_residual_kg":0.}
            for name,row in [("energy_budget.csv",energy),("diagnostics.csv",diagnostic),("mass_budget.csv",mass)]:
                with (output/name).open("w",newline="") as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader();writer.writerow(row)
            return inspect(run)

    def test_component_signs_are_not_separate_grad_div_passivity_tests(self):
        result=self.audit(acoustic_inner_signed_loss_J=-0.02,acoustic_wall_signed_loss_J=0.03,acoustic_boundary_work_J=0.)
        self.assertEqual(result["status"],"sampled_checks_passed")
        acoustic=result["acoustic_work_assessment"]
        self.assertAlmostEqual(acoustic["combined_dissipation_last_J"],0.01)
        self.assertTrue(acoustic["dissipation_sign_checked"])
        self.assertEqual(len(acoustic["negative_component_observations"]),1)

    def test_moving_wall_counterexample_has_negative_raw_work_and_positive_dissipation(self):
        # One fluid mass and a prescribed wall, using the production grad-div
        # pair formula with unit volume/density and G=dW*e=1.
        dt,coefficient,fluid_velocity,wall_velocity=0.01,0.1,0.,-0.1
        divergence=-2*(fluid_velocity-wall_velocity)
        force=2*coefficient*divergence
        fluid_work=dt*force*(fluid_velocity+0.5*dt*force)
        boundary_work=dt*force*wall_velocity
        raw_transfer=-fluid_work
        self.assertLess(raw_transfer,0)
        self.assertGreater(boundary_work-fluid_work,0)
        result=self.audit(kinetic_J=fluid_work,acoustic_wall_signed_loss_J=raw_transfer)
        self.assertEqual(result["status"],"sampled_checks_passed")
        acoustic=result["acoustic_work_assessment"]
        self.assertFalse(acoustic["dissipation_sign_checked"])
        self.assertIsNone(acoustic["combined_dissipation_last_J"])
        self.assertTrue(acoustic["negative_component_observations"])
        measured=self.audit(kinetic_J=fluid_work,acoustic_wall_signed_loss_J=raw_transfer,acoustic_boundary_work_J=boundary_work)
        self.assertGreater(measured["acoustic_work_assessment"]["combined_dissipation_last_J"],0)

    def test_measured_negative_combined_acoustic_dissipation_still_warns(self):
        result=self.audit(acoustic_inner_signed_loss_J=-0.2,acoustic_wall_signed_loss_J=0.1,acoustic_boundary_work_J=0.)
        self.assertEqual(result["status"],"needs_attention")
        self.assertEqual(result["negative_loss_flags"][0]["term"],"acoustic_combined_dissipation_J")

    def test_unambiguous_negative_dissipation_terms_still_warn(self):
        for term in ("pd_damping_signed_loss_J","viscous_signed_loss_J","friction_signed_loss_J"):
            with self.subTest(term=term):
                result=self.audit(**{term:-0.1})
                self.assertEqual(result["status"],"needs_attention")
                self.assertEqual(result["negative_loss_flags"][0]["term"],term)

    def test_mechanical_residual_gate_is_unchanged(self):
        result=self.audit(acoustic_wall_signed_loss_J=-0.1,mechanical_residual_relative=0.06)
        self.assertEqual(result["status"],"needs_attention")
        self.assertIn("5%",result["warnings"][0])


if __name__=="__main__":
    unittest.main()

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from peri_study import compare_group, relative_curve_error, run as run_study
from peri_run import parameters, digest


class StudyAnalysisTests(unittest.TestCase):
    def test_warned_completed_case_stops_queue_before_next_case(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/"frozen_source").mkdir()
            (root/"runs/baseline").mkdir(parents=True)
            config = root/"baseline.cfg"; config.write_text("shell_layers=3\n")
            data = {"source_sha256": {}, "cases": [
                {"name": "baseline", "group": "baseline", "config": "baseline.cfg", "config_sha256": digest(config)},
                {"name": "dt_half", "group": "time", "config": "not-yet-needed.cfg"}]}
            (root/"study.json").write_text(json.dumps(data))
            args = SimpleNamespace(root=root, memory_limit_gib=15, baseline_memory_gib=13, no_states=True)
            with patch("peri_study.verify"), patch("peri_study.require_completion_audit", side_effect=ValueError("negative viscous loss")), patch("peri_study.subprocess.run") as launch:
                self.assertEqual(run_study(args), 2)
                launch.assert_not_called()
            result = json.loads((root/"status.json").read_text())
            self.assertEqual(result["status"], "audit_failed")
            self.assertEqual(result["failed_case"], "baseline")
            self.assertIn("negative viscous", result["reason"])

    def test_full_case_inputs_accept_boolean_tokens(self):
        config = parameters(Path(__file__).resolve().parents[1] / "tests/3d_examples/extra_3d_PD_capsule_compression/configs/full_capsule.cfg")
        self.assertEqual(config["with_liquid"], 1)
        self.assertEqual(config["reinitialize_density"], 0)
        self.assertEqual(config["capsule_radius_out"], 0.0005)

    def test_curve_norm_matches_known_scale_error(self):
        reference = [{"stroke_m": x, "total_Fy_N": 2*x} for x in (0, 0.25, 0.5, 1)]
        coarse = [{"stroke_m": x, "total_Fy_N": 3*x} for x in (0, 0.5, 1)]
        self.assertAlmostEqual(relative_curve_error(coarse, reference), 0.5)

    def test_partial_stroke_cannot_be_called_convergence(self):
        a = [{"stroke_m": 0., "total_Fy_N": 0.}, {"stroke_m": 0.5, "total_Fy_N": 1.}]
        b = [{"stroke_m": 0., "total_Fy_N": 0.}, {"stroke_m": 1., "total_Fy_N": 2.}]
        with self.assertRaisesRegex(ValueError, "same full stroke"):
            relative_curve_error(a, b)

    def test_material_change_is_rejected_even_if_runs_completed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = ["baseline", "dt_half", "dt_quarter"]
            data = {"cases": [{"name": name} for name in names],
                    "executable": "frozen_source/bin/solver", "source_sha256": {"bin/solver": "same"},
                    "force_relative_tolerance": 0.05, "energy_monitoring_tolerance": 0.05}
            for i, name in enumerate(names):
                run = root / "runs" / name
                output = run / "output"
                output.mkdir(parents=True)
                (run / "run_manifest.json").write_text(json.dumps({"executable_sha256": "same"}))
                (output / "effective_parameters.cfg").write_text(f"youngs_shell={1e7 if i != 1 else 2e7}\npd_cfl={0.1/2**i}\nacoustic_cfl={0.3/2**i}\n")
                (output / "diagnostics.csv").write_text("stroke_m,total_Fy_N,damage_mean\n0,0,0\n1,2,0\n")
                (output / "energy_budget.csv").write_text("plate_contact_work_J,plate_fluid_pressure_work_J,plate_fluid_viscous_work_J\n1,0,0\n")
            with patch("peri_study.verify", return_value={"mechanical_energy_residual_max_relative": 0.001}), patch("peri_study.require_completion_audit", return_value={"simulation_completed": True}):
                with self.assertRaisesRegex(ValueError, "uncontrolled parameter"):
                    compare_group(root, data, names, "time")

    def test_accidental_coarse_fine_agreement_is_not_convergence(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); names=["baseline","dt_half","dt_quarter"]
            data={"executable":"frozen_source/bin/solver","source_sha256":{"bin/solver":"same"},
                  "force_relative_tolerance":0.5,"energy_monitoring_tolerance":0.05}
            for i,(name,factor) in enumerate(zip(names,(1,0.8,1))):
                output=root/"runs"/name/"output";output.mkdir(parents=True)
                (output.parent/"run_manifest.json").write_text(json.dumps({"executable_sha256":"same"}))
                (output/"effective_parameters.cfg").write_text(f"youngs_shell=10000000\npd_cfl={0.1/2**i}\nacoustic_cfl={0.3/2**i}\n")
                (output/"diagnostics.csv").write_text(f"stroke_m,total_Fy_N,damage_mean\n0,0,0\n1,{2*factor},0\n")
                (output/"energy_budget.csv").write_text(f"plate_contact_work_J,plate_fluid_pressure_work_J,plate_fluid_viscous_work_J\n{factor},0,0\n")
            with patch("peri_study.verify",return_value={"mechanical_energy_residual_max_relative":0.001}), patch("peri_study.require_completion_audit", return_value={"simulation_completed": True}):
                result=compare_group(root,data,names,"time")
            self.assertTrue(result["adjacent_difference_within_tolerance"])
            self.assertFalse(result["refinement_difference_trend"]["reaction_curve_relative_L2"])
            self.assertFalse(result["accepted"])


if __name__ == "__main__":
    unittest.main()

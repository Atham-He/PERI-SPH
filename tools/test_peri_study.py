import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from peri_study import compare_group, relative_curve_error


class StudyAnalysisTests(unittest.TestCase):
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
            with patch("peri_study.verify", return_value={"mechanical_energy_residual_max_relative": 0.001}):
                with self.assertRaisesRegex(ValueError, "uncontrolled parameter"):
                    compare_group(root, data, names, "time")


if __name__ == "__main__":
    unittest.main()

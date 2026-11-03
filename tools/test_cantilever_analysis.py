"""Scientific-analysis regressions: known damped modes and misleading completions."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from analyze_cantilever import analyze, frequency


class CantileverAnalysisTests(unittest.TestCase):
    def test_constant_or_strongly_multimodal_signal_is_not_a_valid_mode(self):
        t = np.linspace(0, 2, 321)
        self.assertFalse(frequency(t, np.ones_like(t))["fit_valid"])
        y = np.sin(2*np.pi*0.91*t) + 0.5*np.sin(2*np.pi*5.7*t)
        self.assertFalse(frequency(t, y)["fit_valid"])

    def test_frequency_and_decay_from_known_signal(self):
        t = np.linspace(0, 2, 321)
        y = 0.0007 * np.exp(-0.08*t) * np.sin(2*np.pi*1.073*t + 0.13) + 0.00001
        result = frequency(t, y)
        self.assertTrue(result["fit_valid"])
        self.assertAlmostEqual(result["frequency_ratio"], 1.073, places=6)
        self.assertAlmostEqual(result["decay_per_EB_period"], 0.08, places=6)
        self.assertLess(result["fit_residual_relative"], 1e-7)

    def test_known_undamped_mode_with_small_second_mode(self):
        t = np.linspace(0, 2, 321)
        y = np.sin(2*np.pi*0.91*t) + 0.001*np.sin(2*np.pi*5.7*t)
        result = frequency(t, y)
        self.assertLess(abs(result["frequency_ratio"] - 0.91), 1e-4)
        self.assertLess(result["fit_residual_relative"], 0.002)

    def test_failed_run_is_not_fitted_as_success(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/"input.json").write_text(json.dumps({"parameters": {"load_mode": 0}}))
            (p/"result.json").write_text(json.dumps({"status": "failed", "elapsed_seconds": 1, "warnings": ["failed"]}))
            (p/"failure.txt").write_text("small strain envelope exceeded")
            result = analyze(p)
            self.assertEqual(result["status"], "failed")
            self.assertNotIn("frequency_ratio", result)

    def test_completed_run_with_spurious_kinetic_growth_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/"input.json").write_text(json.dumps({"parameters": {
                "load_mode": 0, "cycles": 2, "samples_per_cycle": 1, "amplitude_ratio": 0.001}}))
            (p/"result.json").write_text(json.dumps({"status": "completed", "elapsed_seconds": 1}))
            (p/"completion.json").write_text(json.dumps({"time": 2}))
            (p/"metadata.json").write_text(json.dumps({"period_eb": 1, "particles": 10, "directed_bonds": 20}))
            names = ["time", "t_over_T", "kinetic", "max_F_minus_I", "clamp_drift_over_H",
                     "modal_over_H", "hg_net_relative", "hg_torque_relative", "internal_net_relative",
                     "live_count", "damaged_count", "failure_count", "steps"]
            values = np.zeros((3, len(names)))
            values[:, 0] = values[:, 1] = [0, 1, 2]
            values[:, 2] = [1, 10, 100]
            values[:, 3] = 0.001  # Still small strain; completion alone is misleading.
            values[:, 5] = 0.001
            values[:, 9] = 10
            values[:, 12] = [0, 100, 200]
            np.savetxt(p/"history.csv", values, delimiter=",", header=",".join(names), comments="")
            (p/"profiles.csv").write_text("t_over_T,x_over_L\n0,0.5\n2,0.5\n")
            result = analyze(p)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["quality_status"], "rejected_energy_growth")
            self.assertFalse(result["fit_valid"])
            self.assertNotIn("frequency_ratio", result)


if __name__ == "__main__":
    unittest.main()

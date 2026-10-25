"""Synthetic signals exercise scientific exclusions without model integrations."""
import contextlib
import csv
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/summarize_experiments.py"
SPEC = importlib.util.spec_from_file_location("experiment_summary", str(SCRIPT))
summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary)


class ExperimentSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.options = summary.parser().parse_args([str(self.root), "--output", str(self.root / "out")])
        self.config = {"initial": {"jet_width_m": 1e6, "jet_radius_m": 15e6},
                       "planet": {"rotation_rate_s": 0.000165},
                       "output": {"diagnostic_interval_s": 2000.0, "m_max": 32}}

    def tearDown(self):
        self.temp.cleanup()

    def columns(self, growth=2e-6, angular=1.7e-6, mode=6):
        t = np.arange(0, 202000.0, 2000.0)
        amp = 1500 * np.exp(growth * t)
        coefficient = amp / 2 * np.exp(1j * (-2.7 - mode * angular * t))
        columns = {"time_s": t, "mean_peak_speed_m_s": np.full(len(t), 100.0),
                   "jet_coverage": np.ones(len(t))}
        for m in range(1, 33):
            for suffix in ("amplitude_m", "real_m", "imag_m"):
                columns["m%d_%s" % (m, suffix)] = np.zeros(len(t))
        columns.update({"m%d_amplitude_m" % mode: amp, "m%d_real_m" % mode: coefficient.real,
                        "m%d_imag_m" % mode: coefficient.imag})
        return columns

    def test_recovers_growing_rotating_wave_and_frame(self):
        result = summary.measure_mode(self.columns(), self.config, 6, self.options)
        self.assertAlmostEqual(result["growth_rate_s"], 2e-6, places=17)
        self.assertAlmostEqual(result["angular_pattern_speed_model_s"], 1.7e-6, places=17)
        self.assertAlmostEqual(result["angular_pattern_speed_system_iii_s"], 2.9150098226496324e-6, places=17)
        self.assertAlmostEqual(result["log_amplitude_r2"], 1)
        self.assertAlmostEqual(result["phase_r2"], 1)
        self.assertEqual(result["fit_start_s"], 20000)
        self.assertEqual(result["fit_end_s"], 200000)

    def test_never_unwraps_across_invalid_gap(self):
        columns = self.columns()
        bad = (columns["time_s"] >= 80000) & (columns["time_s"] <= 120000)
        columns["jet_coverage"][bad] = np.nan
        # An arbitrary independent orientation after the gap must not bias the fit.
        late = columns["time_s"] > 120000
        z = columns["m6_real_m"] + 1j * columns["m6_imag_m"]
        z[late] *= np.exp(2.4j)
        columns["m6_real_m"], columns["m6_imag_m"] = z.real, z.imag
        result = summary.measure_mode(columns, self.config, 6, self.options)
        self.assertEqual(result["qualifying_segments"], 2)
        self.assertEqual(result["fit_start_s"], 122000)
        self.assertAlmostEqual(result["angular_pattern_speed_model_s"], 1.7e-6, places=17)

    def test_missing_rows_break_segments(self):
        columns = self.columns()
        keep = (columns["time_s"] < 80000) | (columns["time_s"] > 120000)
        columns = {key: value[keep] for key, value in columns.items()}
        result = summary.measure_mode(columns, self.config, 6, self.options)
        self.assertEqual(result["qualifying_segments"], 2)
        self.assertEqual(result["fit_start_s"], 122000)

    def test_weak_nonfinite_and_short_segments_are_unavailable(self):
        for damage in ("weak", "nonfinite", "short"):
            columns = self.columns()
            if damage == "weak":
                for key in ("m6_amplitude_m", "m6_real_m", "m6_imag_m"):
                    columns[key] /= 100
            elif damage == "nonfinite":
                columns["m6_imag_m"][:] = np.nan
            else:
                columns["jet_coverage"][:] = 0
                columns["jet_coverage"][10:19] = 1
            result = summary.measure_mode(columns, self.config, 6, self.options)
            self.assertIsNone(result["growth_rate_s"], damage)
            self.assertIsNone(result["angular_pattern_speed_model_s"], damage)

    def test_linear_ceiling_and_mean_wind_gate(self):
        columns = self.columns(growth=3e-5)
        columns["mean_peak_speed_m_s"][columns["time_s"] > 100000] = 106
        result = summary.measure_mode(columns, self.config, 6, self.options)
        self.assertEqual(result["fit_end_s"], 100000)
        self.assertAlmostEqual(result["growth_rate_s"], 3e-5, places=17)
        columns["mean_peak_speed_m_s"][:] = 100
        result = summary.measure_mode(columns, self.config, 6, self.options)
        self.assertLess(result["fit_end_s"], 120000)
        self.assertLessEqual(columns["m6_amplitude_m"][int(result["fit_end_s"] / 2000)], 50000)

    def test_ambiguous_sampling_suppresses_phase_only(self):
        angular = 0.9 * math.pi / (6 * 2000)
        result = summary.measure_mode(self.columns(angular=angular), self.config, 6, self.options)
        self.assertEqual(result["phase_status"], "ambiguous_phase_increment")
        self.assertIsNone(result["angular_pattern_speed_model_s"])
        self.assertAlmostEqual(result["growth_rate_s"], 2e-6, places=17)

    def test_fast_competing_mode_vetoes_late_linear_fit(self):
        columns = self.columns()
        competitor = 1000 * np.exp(5e-5 * columns["time_s"])
        columns["m14_amplitude_m"] = competitor
        columns["m14_real_m"] = competitor / 2
        result = summary.measure_mode(columns, self.config, 6, self.options)
        self.assertEqual(result["fit_end_s"], 84000)
        self.assertLessEqual(result["maximum_resolved_centreline_rms_m"], 50000)
        self.assertAlmostEqual(result["growth_rate_s"], 2e-6, places=17)
        # Mode six alone stays small throughout; the competing mode sets the cutoff.
        self.assertLess(np.max(columns["m6_amplitude_m"]), 50000)
        self.assertGreater(result["exclusion_counts_nonexclusive"]["combined_resolved_amplitude_outside_linear_interval"], 0)
        columns["m14_amplitude_m"][:] = np.nan
        self.assertIsNone(summary.measure_mode(columns, self.config, 6, self.options)["growth_rate_s"])
        del columns["m14_amplitude_m"]
        self.assertEqual(summary.measure_mode(columns, self.config, 6, self.options)["fit_status"],
                         "missing_configured_mode_amplitudes")

    def write_run(self, complete=True):
        path = self.root / "run"
        path.mkdir()
        (path / "config.toml").write_text('''[experiment]
name = "synthetic"
[initial]
jet_width_m = 1000000.0
jet_radius_m = 15000000.0
[planet]
rotation_rate_s = 0.000165
[numerics]
duration_s = 200000.0
[output]
diagnostic_interval_s = 2000.0
m_max = 32
''')
        (path / "metadata.toml").write_text('device = "synthetic"\n')
        columns = self.columns()
        with (path / "diagnostics.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerows(zip(*columns.values()))
        if complete:
            (path / "COMPLETED").write_text("200000.0\n")
        return path

    def test_failed_and_incomplete_are_not_completed(self):
        path = self.write_run(False)
        self.assertEqual(summary.summarize_run(path, self.options)["execution_status"], "incomplete")
        (path / "COMPLETED").write_text("200000\n")
        self.assertEqual(summary.summarize_run(path, self.options)["execution_status"], "completed")
        (path / "COMPLETED").write_text("100000\n")
        self.assertEqual(summary.summarize_run(path, self.options)["execution_status"], "incomplete")
        (path / "FAILED.txt").write_text("CFL failure")
        self.assertEqual(summary.summarize_run(path, self.options)["execution_status"], "failed")

    def test_cli_strict_json_and_no_overwrite(self):
        path = self.write_run()
        output = self.root / "summary"
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(summary.main([str(path), "--output", str(output)]), 0)
        result = json.loads((output / "summary.json").read_text(), parse_constant=lambda s: self.fail(s))
        self.assertEqual(result["runs"][0]["mode_six"]["fit_status"], "estimated_provisional")
        self.assertIsNone(result["runs"][0]["modes"][0]["growth_rate_s"])
        with (output / "modes.csv").open() as stream:
            self.assertEqual(len(list(csv.DictReader(stream))), 32)
        old = (output / "summary.json").read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            summary.main([str(path), "--output", str(output)])
        self.assertEqual(old, (output / "summary.json").read_bytes())

    def test_invalid_chronology_is_not_silently_sorted(self):
        path = self.write_run()
        csvpath = path / "diagnostics.csv"
        lines = csvpath.read_text().splitlines()
        lines[2], lines[3] = lines[3], lines[2]
        csvpath.write_text("\n".join(lines) + "\n")
        result = summary.summarize_run(path, self.options)
        self.assertEqual(result["input_status"], "invalid_or_missing_input")
        self.assertEqual(result["modes"], [])

    def test_optional_regime_extrema_preserve_missing_and_nonfinite_status(self):
        finite_ld = {"physics": {"deformation_radius_m": 2e6}}
        columns = {"max_abs_hprime_over_H": np.array([0.1, 0.3, 0.2]),
                   "min_h_over_H": np.array([0.9, 0.7, 0.8]),
                   "max_abs_zeta_over_f": np.array([0.2, np.nan, 0.5])}
        result = summary.regime_summary(columns, finite_ld)
        self.assertEqual(result["maximum_abs_hprime_over_H"], 0.3)
        self.assertEqual(result["minimum_h_over_H"], 0.7)
        self.assertEqual(result["maximum_abs_zeta_over_f"], 0.5)
        self.assertEqual(result["columns"]["max_abs_zeta_over_f"]["status"], "partly_nonfinite")
        self.assertEqual(result["columns"]["max_abs_zeta_over_f"]["finite_samples"], 2)
        old = summary.regime_summary({}, finite_ld)
        self.assertIsNone(old["maximum_abs_hprime_over_H"])
        self.assertEqual(old["columns"]["max_abs_hprime_over_H"]["status"], "absent_legacy_or_unrecorded_column")

    def test_barotropic_thickness_is_undefined_not_zero(self):
        columns = {"max_abs_hprime_over_H": np.array([np.nan, np.nan]),
                   "min_h_over_H": np.array([np.nan, np.nan]),
                   "max_abs_zeta_over_f": np.array([0.2, 0.4])}
        result = summary.regime_summary(columns, {"physics": {"deformation_radius_m": math.inf}})
        self.assertIsNone(result["maximum_abs_hprime_over_H"])
        self.assertIsNone(result["minimum_h_over_H"])
        self.assertEqual(result["columns"]["min_h_over_H"]["status"], "not_applicable_barotropic")
        self.assertEqual(result["maximum_abs_zeta_over_f"], 0.4)
        invalid = summary.regime_summary(columns, {"physics": {"deformation_radius_m": 2e6}})
        self.assertEqual(invalid["columns"]["min_h_over_H"]["status"], "no_finite_samples")
        columns["min_h_over_H"][:] = 1
        unexpected = summary.regime_summary(columns, {"physics": {"deformation_radius_m": math.inf}})
        self.assertEqual(unexpected["columns"]["min_h_over_H"]["status"], "unexpected_finite_barotropic_depth_values")
        self.assertIsNone(unexpected["minimum_h_over_H"])

    def test_regime_columns_are_optional_and_exported_per_run(self):
        path = self.write_run()
        old = summary.summarize_run(path, self.options)
        self.assertEqual(old["regime_summary"]["columns"]["max_abs_zeta_over_f"]["status"],
                         "absent_legacy_or_unrecorded_column")
        csvpath = path / "diagnostics.csv"
        with csvpath.open(newline="") as stream:
            rows = list(csv.reader(stream))
        rows[0].extend(["max_abs_hprime_over_H", "min_h_over_H", "max_abs_zeta_over_f"])
        for row in rows[1:]:
            row.extend(["nan", "nan", "0.27"])
        with csvpath.open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)
        with (path / "config.toml").open("a") as stream:
            stream.write("\n[physics]\ndeformation_radius_m = inf\n")
        output = self.root / "regime_output"
        with contextlib.redirect_stdout(io.StringIO()):
            summary.main([str(path), "--output", str(output)])
        result = json.loads((output / "summary.json").read_text(), parse_constant=lambda s: self.fail(s))
        regime = result["runs"][0]["regime_summary"]
        self.assertEqual(regime["columns"]["min_h_over_H"]["status"], "not_applicable_barotropic")
        self.assertEqual(regime["maximum_abs_zeta_over_f"], 0.27)
        self.assertIn("max |zeta/f|=0.27", (output / "report.txt").read_text())


if __name__ == "__main__":
    unittest.main()

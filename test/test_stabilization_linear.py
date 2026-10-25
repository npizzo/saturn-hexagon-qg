"""Small analytic/manifest tests only; no large eigenproblems or cluster jobs."""
import contextlib
import copy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import toml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stabilization_linear", str(ROOT/"scripts/stabilization_linear.py"))
linear = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(linear)
radial = linear.radial


class StabilizationLinearTests(unittest.TestCase):
    def config(self, matched=True):
        result = {"initial": {"jet_speed_m_s": 2.0, "jet_radius_m": 1.0, "jet_width_m": 0.2,
                            "cyclone_speed_m_s": 0.4, "cyclone_radius_m": 0.3},
                "physics": {"deformation_radius_m": 0.7, "deep_jet_speed_m_s": 0.0 if matched else 0.8,
                            "deep_jet_radius_m": 1.0, "deep_jet_width_m": 0.2,
                            "deep_profile": "matched_upper" if matched else "annular"},
                "planet": {"rotation_rate_s": 0.1, "radius_m": 4.0},
                "grid": {"nx": 128, "length_m": 8.0},
                "boundary": {"inner_radius_m": 2.0}}
        if matched:
            result["physics"]["deep_match_fraction"] = 0.75
        return result

    def fixture(self, directory, configs):
        project = Path(directory)/"public checkout"; (project/"scripts").mkdir(parents=True)
        (project/"src").mkdir(); (project/"src/example.jl").write_text("# immutable source fixture\n")
        (project/"Manifest.toml").write_text('julia_version="1.10.12"\n')
        paths = ["Manifest.toml", "src/example.jl", "scripts/radial_stability.py", "scripts/stabilization_linear.py"]
        for name in paths[2:]:
            (project/name).write_bytes((ROOT/name).read_bytes())
        campaign = project/"configs/stabilization_01"; stage = campaign/"screen"; stage.mkdir(parents=True)
        (campaign/"provenance.tsv").write_text("".join(hashlib.sha256((project/name).read_bytes()).hexdigest()+"\t"+name+"\n" for name in paths))
        members = []
        for index, config in enumerate(configs):
            name = "case%d.toml" % index
            data = toml.dumps(config).encode(); (stage/name).write_bytes(data)
            members.append(hashlib.sha256(data).hexdigest()+"\t"+name+"\n")
        (stage/"members.tsv").write_text("".join(members))
        return project, stage

    def test_fixed_config_controls_and_complete_spectra(self):
        c = self.config(); before = copy.deepcopy(c)
        rows, spectra, functions = linear.screen_config(c, [16], 1, 3, 2.0, 1e-12)
        self.assertEqual(c, before)
        self.assertEqual(len(rows), 3); self.assertEqual(len(spectra), 48)
        self.assertTrue(all("single_radial_grid" in row["flags"] for row in rows))
        self.assertTrue(all("weak_radial_resolution" in row["flags"] for row in rows))
        self.assertIn("N16_r_m", functions)
        r, op = radial.radial_operator(16, 2.0, 2, 0.7)
        b = radial.circular_base(c, r, 0.7, match_fraction=0.75)
        expected = radial.solve_modes(r, op, b["omega"], b["q_r_over_r"], 2.0, 2)
        self.assertAlmostEqual(rows[1]["growth_rate_s"], expected["growth_rate_s"], places=14)

    def test_neutral_spectrum_has_no_arbitrary_drift(self):
        c = self.config(); c["initial"]["jet_speed_m_s"] = 0; c["initial"]["cyclone_speed_m_s"] = 0
        rows, spectra, functions = linear.screen_config(c, [16], 1, 3, 2.0, 1e-12)
        self.assertTrue(all(row["growing_angular_pattern_speed_s"] is None for row in rows))
        self.assertTrue(all(row["growing_frequency_real_s"] is None for row in rows))
        self.assertFalse(any(key.endswith("_phi") for key in functions))
        self.assertEqual(len(spectra), 48)

    def test_high_mode_edge_and_resolution_flags_are_explicit(self):
        flags, cells = linear.resolution_flags(self.config(), 2.0, 16, 1, 3, [0, 1e-3, 2e-3], 1e-12, [16])
        self.assertIn("fastest_growth_at_high_m_edge_extend_search", flags)
        self.assertIn("unexamined_modes_above_3", flags)
        self.assertLess(cells["jet_width"], 8)

    def test_stage_hashes_and_outputs_in_space_path(self):
        with tempfile.TemporaryDirectory() as folder:
            project, stage = self.fixture(folder, [self.config(), self.config(False)])
            out = Path(folder)/"linear output"
            with contextlib.redirect_stdout(io.StringIO()):
                result = linear.run(stage, out, nodes=(16,), modes=(1, 3), project=project)
            self.assertEqual(len(result), 2)
            self.assertTrue((out/"COMPLETED").is_file())
            meta = json.loads((out/"case0/metadata.json").read_text())
            self.assertEqual(meta["deep_match_fraction"], 0.75)
            self.assertEqual(meta["deformation_radius_m"], 0.7)
            self.assertEqual((out/"case0/config.toml").read_bytes(), (stage/"case0.toml").read_bytes())
            self.assertEqual(len((out/"case0/eigenvalues.csv").read_text().splitlines()), 49)
            with self.assertRaises(ValueError):
                linear.run(stage, out, nodes=(16,), modes=(1, 3), project=project)
            (stage/"case0.toml").write_text("# corrupted\n")
            with self.assertRaisesRegex(ValueError, "Hash mismatch"):
                linear.load_stage(stage, project)

    def test_changed_source_and_escape_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            project, stage = self.fixture(folder, [self.config()])
            (project/"src/example.jl").write_text("# changed\n")
            with self.assertRaisesRegex(ValueError, "Hash mismatch"):
                linear.load_stage(stage, project)
            (stage/"members.tsv").write_text("0"*64+"\t../escape.toml\n")
            with self.assertRaises(ValueError):
                linear.hashed_files(stage/"members.tsv", stage, basenames=True)

    def test_cli_rejects_ignored_deep_control_and_old_annular_cli_works(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"config.toml"
            path.write_text(toml.dumps(self.config()))
            common = [str(path), "--output", str(Path(folder)/"out"), "--m-range", "1:1",
                      "--radial-nodes", "8", "--outer-radius", "2"]
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                radial.main(common+["--deep-ratios", "0.5"])
            with contextlib.redirect_stdout(io.StringIO()):
                matched = radial.main(common+["--match-fractions", "0.75"])
            self.assertEqual(matched[0]["deep_match_fraction"], 0.75)
            self.assertIsNone(matched[0]["deep_speed_m_s"])
            c = self.config(False); c["physics"].pop("deep_profile")
            path.write_text(toml.dumps(c))
            common[2] = str(Path(folder)/"old_cli")
            with contextlib.redirect_stdout(io.StringIO()):
                annular = radial.main(common+["--deep-ratios", "0.5"])
            self.assertEqual(annular[0]["deep_speed_m_s"], 1.0)
            self.assertEqual(annular[0]["deep_profile"], "annular")


if __name__ == "__main__":
    unittest.main(verbosity=2)

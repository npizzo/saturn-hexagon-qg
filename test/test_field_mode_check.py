"""Analytic field modes and Julia/HDF5 orientation; no model integrations."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

SPEC = importlib.util.spec_from_file_location("field_mode_check", str(Path(__file__).resolve().parents[1]/"scripts/field_mode_check.py"))
check = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(check)


class FieldModeTests(unittest.TestCase):
    def signal(self, phase=0.37, mode=6, nx=256, ny=288):
        x, y = np.linspace(-4, 4, nx, endpoint=False), np.linspace(-4, 4, ny, endpoint=False)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        r = np.hypot(xx, yy); theta = np.arctan2(yy, xx)
        safe = np.where(r == 0, 1, r)
        g = r**mode*np.exp(-mode*(r*r-1)/2)
        gr = mode*(1/safe-r)*g
        psi = 4*g*np.cos(mode*theta+phase)
        dr = 4*gr*np.cos(mode*theta+phase)
        angular = -4*mode*g*np.sin(mode*theta+phase)/safe
        u = -np.sin(theta)*dr-np.cos(theta)*angular
        v = np.cos(theta)*dr-np.sin(theta)*angular
        return (psi, u, v), x, y

    def measure(self, fields, x, y):
        return check.field_modes(*fields, x, y, 1, 0.1, 2, 2)

    def test_analytic_norm_amplitude_velocity_and_phase(self):
        fields, x, y = self.signal()
        first = self.measure(fields, x, y)
        moved = self.measure(self.signal(phase=0.91)[0], x, y)
        r = first["radii"]
        g = r**6*np.exp(-3*(r*r-1))
        expected = 2*g*np.exp(0.37j)
        self.assertLess(np.max(np.abs(first["psi"][:, 5]-expected)), 5e-4)
        row = check.summarize_modes(first, first)[5]
        self.assertAlmostEqual(row["psi_mode_radial_rms_m2_s"], np.sqrt(check.radial_inner(expected, expected, r).real), delta=2e-4)
        self.assertAlmostEqual(row["equivalent_displacement_at_Rj_m"], 2, delta=2e-4)
        velocity = 12*(1/r-r)*g
        self.assertAlmostEqual(row["utheta_mode_radial_rms_m_s"], np.sqrt(check.radial_inner(velocity, velocity, r).real), delta=4e-4)
        self.assertAlmostEqual(check.summarize_modes(moved, first)[5]["m6_projection_phase_rad"], 0.54, delta=1e-5)
        self.assertIsNone(check.summarize_modes(moved, first)[4]["m6_projection_phase_rad"])
        other = np.delete(first["psi"], 5, axis=1)
        self.assertLess(np.max(np.abs(other)), 2e-4)

    def test_phase_is_wrapped_and_zero_reference_is_null(self):
        fields, x, y = self.signal(phase=0)
        first = self.measure(fields, x, y)
        later = self.measure(self.signal(phase=3.5)[0], x, y)
        self.assertAlmostEqual(check.summarize_modes(later, first)[5]["m6_projection_phase_rad"], 3.5-2*np.pi, delta=1e-5)
        zero = self.measure([np.zeros_like(f) for f in fields], x, y)
        self.assertIsNone(check.summarize_modes(later, zero)[5]["m6_projection_phase_rad"])
        self.assertIsNone(check.summarize_modes(zero, first)[5]["m6_projection_phase_rad"])

    def test_julia_orientation_nonsquare_roundtrip(self):
        fields, x, y = self.signal()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"snapshot_000000001.jld2"
            with h5py.File(str(path), "w") as f:
                f["time_s"] = 123.0
                for name, a in zip(("psi", "u", "v"), fields):
                    f[name] = a.T
            time, restored = check.read_snapshot(path)
            self.assertEqual(time, 123)
            for a, b in zip(fields, restored):
                np.testing.assert_array_equal(a, b)
            self.assertAlmostEqual(check.summarize_modes(self.measure(restored, x, y), self.measure(fields, x, y))[5]["m6_projection_phase_rad"], 0, places=14)

    def test_rejects_bad_geometry_and_nonfinite_fields(self):
        fields, x, y = self.signal()
        with self.assertRaises(ValueError):
            self.measure(fields, x[::-1], y)
        with self.assertRaises(ValueError):
            self.measure([f.T for f in fields], x, y)
        fields[0][0, 0] = np.nan
        with self.assertRaises(ValueError):
            self.measure(fields, x, y)

    def test_saved_run_units_ordering_and_new_output_contract(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)/"results"; run.mkdir()
            (run/"config.toml").write_text('[grid]\nnx=128\nlength_m=8.0\n[initial]\njet_radius_m=1.0\njet_width_m=0.1\njet_speed_m_s=2.0\n[boundary]\ninner_radius_m=2.0\n')
            (run/"metadata.toml").write_text('units="SI; planar distance in m; velocity m/s; psi m2/s"\naxes="arrays indexed [x,y]"\n')
            for name, time, phase in (("snapshot_000000009.jld2", 0, 0), ("snapshot_000000001.jld2", 50000, 0.7)):
                fields, x, y = self.signal(phase=phase, nx=128, ny=128)
                with h5py.File(str(run/name), "w") as f:
                    f["time_s"] = time
                    for key, field in zip(("psi", "u", "v"), fields):
                        f[key] = field.T
            with h5py.File(str(run/"domain.jld2"), "w") as f:
                f["x"] = x; f["y"] = y
            output = Path(folder)/"audit"
            rows = check.analyze(run, output)
            self.assertEqual(len(rows), 64)
            self.assertEqual(rows[0]["time_s"], 0)
            self.assertEqual(rows[-1]["time_s"], 50000)
            self.assertAlmostEqual(rows[37]["m6_projection_phase_rad"], 0.7, delta=1e-4)
            self.assertTrue((output/"fields.csv").is_file())
            self.assertEqual(json.loads((output/"metadata.json").read_text())["sampling"]["ntheta"], 720)
            with self.assertRaises(ValueError):
                check.analyze(run, output)
            (run/"metadata.toml").write_text('units="unknown"\naxes="arrays indexed [x,y]"\n')
            with self.assertRaises(ValueError):
                check.analyze(run, Path(folder)/"bad")
            self.assertFalse((Path(folder)/"bad").exists())


if __name__ == "__main__":
    unittest.main()

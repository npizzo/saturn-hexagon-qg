"""Analytic and convergence checks; no Julia evolution or cluster submission."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import jn_zeros, jv

path = Path(__file__).resolve().parents[1] / "scripts/radial_stability.py"
spec = importlib.util.spec_from_file_location("radial_stability", str(path))
radial = importlib.util.module_from_spec(spec)
spec.loader.exec_module(radial)


class RadialStabilityTests(unittest.TestCase):
    def config(self):
        return {"initial": {"jet_speed_m_s": 2.0, "jet_radius_m": 1.0,
                            "jet_width_m": 0.2, "cyclone_speed_m_s": 0.4,
                            "cyclone_radius_m": 0.3},
                "physics": {"deep_jet_radius_m": 1.0, "deep_jet_width_m": 0.2},
                "planet": {"rotation_rate_s": 0.1, "radius_m": 4.0}}

    def test_analytic_profile_derivatives(self):
        r = np.linspace(0.05, 1.4, 37)
        for profile, args in [(radial.ring_profile, (2.0, 1.0, 0.2)),
                              (radial.cyclone_profile, (0.4, 0.3))]:
            u, up, upp = profile(r, *args)
            complex_u, complex_up, _ = profile(r + 1e-24j, *args)
            np.testing.assert_allclose(up, complex_u.imag / 1e-24, rtol=3e-12, atol=1e-12)
            np.testing.assert_allclose(upp, complex_up.imag / 1e-24, rtol=3e-12, atol=1e-12)

    def test_pv_gradient_and_pole_regularity(self):
        c = self.config()
        r = np.linspace(0.01, 1.8, 100)
        b = radial.circular_base(c, r, 0.7, 0.8)
        direct = (b["u_second"] + b["u_prime"] / r - b["u"] / r**2
                  - (b["u"] - b["u_deep"]) / 0.7**2 - b["gamma_code"] * r) / r
        np.testing.assert_allclose(b["q_r_over_r"], direct, rtol=3e-12, atol=2e-11)
        pole = radial.circular_base(c, np.array([0.0, 1e-12, 1e-9]), 0.7, 0.8)
        self.assertTrue(all(np.all(np.isfinite(v)) for v in pole.values()))
        np.testing.assert_allclose(pole["omega"], pole["omega"][0], rtol=1e-14)
        np.testing.assert_allclose(pole["q_r_over_r"], pole["q_r_over_r"][0], rtol=1e-14)
        self.assertEqual(pole["u"][0], 0.0)
        self.assertEqual(pole["q_r"][0], 0.0)

    def test_matching_deep_velocity_cancels_stretching_pv_gradient(self):
        c = self.config()
        c["initial"]["cyclone_speed_m_s"] = 0.0
        r = np.linspace(0.0, 2.0, 80)
        finite = radial.circular_base(c, r, 0.4, c["initial"]["jet_speed_m_s"])
        barotropic = radial.circular_base(c, r, np.inf, 0.0)
        np.testing.assert_allclose(finite["q_r_over_r"], barotropic["q_r_over_r"], rtol=0, atol=1e-14)

    def test_matched_upper_includes_cyclone_and_rejects_ignored_overrides(self):
        c = self.config()
        c["physics"].update(deep_profile="matched_upper", deep_match_fraction=1.0)
        r = np.linspace(0.0, 2.0, 80)
        finite = radial.circular_base(c, r, 0.4)
        barotropic = radial.circular_base(c, r, np.inf)
        np.testing.assert_allclose(finite["u_deep"], finite["u"], rtol=0, atol=0)
        np.testing.assert_allclose(finite["q_r_over_r"], barotropic["q_r_over_r"], rtol=0, atol=1e-14)
        c["physics"].pop("deep_match_fraction")
        default = radial.circular_base(c, r, 0.4)
        np.testing.assert_allclose(default["u_deep"], finite["u_deep"], rtol=0, atol=0)
        c["physics"]["deep_match_fraction"] = 1.0
        half = radial.circular_base(c, r, 0.4, match_fraction=0.5)
        np.testing.assert_allclose(half["u_deep"], 0.5*finite["u"], rtol=0, atol=0)
        np.testing.assert_allclose(half["q_r_over_r"]-finite["q_r_over_r"], -0.5*finite["omega"]/0.4**2, atol=2e-14)
        with self.assertRaises(ValueError):
            radial.circular_base(c, r, 0.4, deep_speed=0.0)
        for fraction in (-0.1, np.nan, np.inf, True):
            with self.assertRaises(ValueError):
                radial.circular_base(c, r, 0.4, match_fraction=fraction)
        c["physics"]["deep_jet_speed_m_s"] = 0.2
        with self.assertRaisesRegex(ValueError, "deep_jet_speed_m_s=0"):
            radial.circular_base(c, r, 0.4)
        c["physics"]["deep_jet_speed_m_s"] = 0.0
        c["physics"]["deep_profile"] = "annular"
        with self.assertRaisesRegex(ValueError, "deep_match_fraction applies"):
            radial.circular_base(c, r, 0.4)
        c["physics"].pop("deep_match_fraction")
        with self.assertRaises(ValueError):
            radial.circular_base(c, r, 0.4, 0.2, match_fraction=0.5)

    def test_weighted_operator_is_symmetric_negative_definite(self):
        for mode in [1, 6, 16]:
            r, operator = radial.radial_operator(48, 2.0, mode, 0.7)
            np.testing.assert_array_equal(operator, operator.T)
            self.assertLess(eigvalsh(operator)[-1], 0.0)
            self.assertEqual(r[0], 2.0 / 49)
            self.assertLess(r[-1], 2.0)

    def test_gamma_plane_bessel_dispersion_and_second_order_convergence(self):
        radius, mode, gamma, ld = 2.0, 3, 0.07, 0.8
        zeros = jn_zeros(mode, 3)
        expected = -gamma / ((zeros / radius)**2 + ld**-2)
        errors = []
        for nodes in [48, 96, 192]:
            r, operator = radial.radial_operator(nodes, radius, mode, ld)
            result = radial.solve_modes(r, operator, 0.0, -gamma, radius, mode)
            actual = np.sort(result["eigenvalues"].real)[:3]
            self.assertLess(np.max(np.abs(result["eigenvalues"].imag)), 1e-12)
            errors.append(abs(actual[0] - expected[0]))
            np.testing.assert_allclose(actual, expected, rtol=0.01)
        self.assertGreater(errors[0] / errors[1], 3.5)
        self.assertGreater(errors[1] / errors[2], 3.5)

    def test_solid_body_rotation_dispersion(self):
        radius, mode, omega, omega_deep, gamma, ld = 1.5, 2, 0.6, 0.2, 0.1, 0.8
        r, operator = radial.radial_operator(128, radius, mode, ld)
        qr_over_r = -ld**-2 * (omega - omega_deep) - gamma
        result = radial.solve_modes(r, operator, omega, qr_over_r, radius, mode)
        expected = omega + qr_over_r / ((jn_zeros(mode, 3) / radius)**2 + ld**-2)
        np.testing.assert_allclose(np.sort(result["eigenvalues"].real)[:3], expected, rtol=5e-4)
        self.assertLess(abs(result["growth_rate_s"]), 1e-12)

    def test_uniform_pv_gradient_has_real_spectrum_even_with_shear(self):
        r, operator = radial.radial_operator(64, 2.0, 4, np.inf)
        omega = 0.2 + 0.3 * np.exp(-r * r)
        for qr_over_r in [-0.08, 0.0, 0.08]:
            result = radial.solve_modes(r, operator, omega, qr_over_r, 2.0, 4)
            self.assertLess(np.max(np.abs(result["eigenvalues"].imag)), 1e-11)
        zero_gradient = radial.solve_modes(r, operator, omega, 0.0, 2.0, 4)
        np.testing.assert_allclose(np.sort(zero_gradient["eigenvalues"].real), np.sort(omega), rtol=1e-10)

    def test_growth_sign_and_azimuthal_mode_factor(self):
        # A two-dimensional algebra fixture has B^-1 A = [[0,.5],[-.5,0]],
        # hence c=+/-0.5i. This checks conventions, not a physical radial base.
        operator = np.array([[-2.0, 1.0], [1.0, -2.0]])
        result = radial.solve_modes(np.array([0.3, 0.6]), operator,
                                    np.array([-1.0, 1.0]), np.array([2.5, -2.5]), 1.0, 6)
        self.assertAlmostEqual(result["leading_c"].imag, 0.5, places=13)
        self.assertAlmostEqual(result["growth_rate_s"], 3.0, places=13)
        self.assertAlmostEqual(result["angular_pattern_speed_s"], 0.0, places=13)
        self.assertLess(result["relative_residual"], 1e-13)

    def test_outer_boundary_scaling_and_regular_bessel_shape(self):
        mode, nodes = 1, 96
        for radius in [1.0, 2.0]:
            r, operator = radial.radial_operator(nodes, radius, mode, np.inf)
            eigenvalues, eigenvectors = np.linalg.eigh(operator)
            phi = eigenvectors[:, -1] / np.sqrt(r)
            expected = jv(mode, jn_zeros(mode, 1)[0] * r / radius)
            phi *= np.dot(phi, expected) / np.dot(phi, phi)
            np.testing.assert_allclose(phi, expected, atol=2e-4)
            with_boundaries = np.r_[0.0, phi, 0.0]
            self.assertEqual(with_boundaries[0], 0.0)
            self.assertEqual(with_boundaries[-1], 0.0)
            self.assertAlmostEqual(eigenvalues[-1], -jn_zeros(mode, 1)[0]**2, delta=0.005)
        with self.assertRaises(ValueError):
            radial.radial_operator(32, 1.0, 0, np.inf)
        with self.assertRaises(ValueError):
            radial.radial_operator(32, 0.0, 1, np.inf)


if __name__ == "__main__":
    unittest.main(verbosity=2)

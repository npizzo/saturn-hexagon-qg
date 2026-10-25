"""Independent analytic checks of the regularized 2017 QG radial base."""
import contextlib
import copy
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np
import toml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("rostami_radial", str(ROOT/"scripts/radial_stability.py"))
radial = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(radial)


class RostamiLinearTests(unittest.TestCase):
    def config(self):
        return {"initial": {"profile": "rostami2017", "jet_speed_m_s": 81.92,
                "jet_radius_m": 10784000., "jet_width_m": float(3200000/np.sqrt(3)),
                "cyclone_speed_m_s": 153.6,
                "cyclone_radius_m": 3200000*(.42/1.3)**(1/1.3),
                "core_regularization_m": 500000.},
                "physics": {"deformation_radius_m": 3200000., "deep_profile": "matched_components",
                    "deep_jet_fraction": 0., "deep_cyclone_fraction": 0.,
                    "deep_jet_speed_m_s": 0., "deep_jet_radius_m": 10784000.,
                    "deep_jet_width_m": float(3200000/np.sqrt(3))},
                "planet": {"radius_m": 55000000., "rotation_rate_s": .00016},
                "boundary": {"inner_radius_m": 28000000.},
                "grid": {"length_m": 80000000., "nx": 1024}}

    def test_physical_scaling_and_component_peak_normalization(self):
        c = self.config(); p = radial.rostami_parameters(c["initial"])
        self.assertAlmostEqual(p["cyclone_length_m"], 3200000., places=8)
        for function, key, speed in [(radial.rostami_jet_profile, "jet_peak_radius_m", 81.92),
                (radial.rostami_cyclone_profile, "cyclone_peak_radius_m", 153.6)]:
            peak = p[key]
            u, up, upp = function(peak, p)
            self.assertAlmostEqual(float(u), speed, places=10)
            self.assertLess(abs(float(up)), 1e-16)
            self.assertLess(float(upp), 0)
            flank, _, _ = function(np.array([.5,1.5])*peak,p)
            self.assertTrue(np.all(flank < speed))
        self.assertGreater(p["cyclone_peak_radius_m"], c["initial"]["cyclone_radius_m"])

    def test_derivatives_against_complex_step(self):
        c = self.config(); p = radial.rostami_parameters(c["initial"])
        r = np.array([0., 1., 1000., 1e5, 5e5, 1e6, 3e6, 7e6, 10e6, 11e6, 18e6])
        for function in (radial.rostami_jet_profile, radial.rostami_cyclone_profile):
            u, up, upp = function(r,p)
            cs_u, cs_up, _ = function(r+1e-16j,p)
            np.testing.assert_allclose(up,cs_u.imag/1e-16,rtol=8e-14,atol=1e-23)
            np.testing.assert_allclose(upp,cs_up.imag/1e-16,rtol=8e-13,atol=1e-28)

    def test_origin_finite_pv_gradient_and_direct_expression(self):
        c = self.config(); r = np.array([0.,1e-12,1e-7,1.])
        b = radial.circular_base(c,r,3.2e6)
        self.assertTrue(all(np.all(np.isfinite(v)) for v in b.values()))
        self.assertEqual(b["u"][0],0.); self.assertEqual(b["q_r"][0],0.)
        np.testing.assert_allclose(b["omega"],b["omega"][0],rtol=1e-11)
        np.testing.assert_allclose(b["q_r_over_r"],b["q_r_over_r"][0],rtol=1e-11)
        r = np.geomspace(1e3,2.5e7,101); b = radial.circular_base(c,r,3.2e6)
        direct = (b["u_second"]+b["u_prime"]/r-b["u"]/r**2)/r
        np.testing.assert_allclose(b["vorticity_gradient_over_r"],direct,rtol=1e-9,atol=5e-27)
        self.assertTrue(np.all(np.isfinite(radial.circular_base(c,0.,3.2e6)["q_r_over_r"])))

    def test_component_matching_and_ld_hold_upper_velocity_fixed(self):
        c = self.config(); r = np.linspace(0,2.5e7,301)
        rest = radial.circular_base(c,r,3.2e6)
        c["physics"].update(deep_jet_fraction=.7,deep_cyclone_fraction=.2)
        partial = radial.circular_base(c,r,3.2e6)
        np.testing.assert_array_equal(rest["u"],partial["u"])
        np.testing.assert_allclose(partial["u_deep"],.7*rest["u_jet"]+.2*rest["u_cyclone"])
        expected = (.7*rest["omega_jet"]+.2*rest["omega_cyclone"])/(3.2e6)**2
        np.testing.assert_allclose(partial["q_r_over_r"]-rest["q_r_over_r"],expected,rtol=2e-13,atol=1e-30)
        c["physics"].update(deep_jet_fraction=1.,deep_cyclone_fraction=1.)
        finite = radial.circular_base(c,r,2.4e6)
        bt = radial.circular_base(c,r,np.inf)
        np.testing.assert_array_equal(finite["u"],rest["u"])
        np.testing.assert_array_equal(bt["u"],rest["u"])
        np.testing.assert_array_equal(finite["q_r_over_r"],bt["q_r_over_r"])
        # Cancellation in the base PV gradient does not remove LD in inversion.
        _, finite_op = radial.radial_operator(16,2.5e7,6,2.4e6)
        _, bt_op = radial.radial_operator(16,2.5e7,6,np.inf)
        self.assertFalse(np.allclose(finite_op,bt_op))

    def test_polar_smoothing_is_explicit_and_approaches_paper_off_axis(self):
        c = self.config(); c["initial"]["core_regularization_m"] = 1.
        p = radial.rostami_parameters(c["initial"])
        r = np.linspace(2e5,1.5e7,100); l = 3.2e6
        a,b = .42,1.3; rp = l*(a/b)**(1/b)
        norm = (rp/l)**a*np.exp(-(rp/l)**b)
        paper_core = 153.6*(r/l)**a*np.exp(-(r/l)**b)/norm
        actual, _, _ = radial.rostami_cyclone_profile(r,p)
        np.testing.assert_allclose(actual,paper_core,rtol=1e-11)
        near_jet = np.linspace(8e6,14e6,100)
        paper_jet = 81.92*np.exp(-((near_jet-10784000)/(l/np.sqrt(3)))**2)
        actual, _, _ = radial.rostami_jet_profile(near_jet,p)
        np.testing.assert_allclose(actual,paper_jet,rtol=1e-14)

    def test_unsupported_profiles_and_ignored_parameters_rejected(self):
        base = self.config()
        for group,key,value in [("initial","profile","unknown"),("initial","core_regularization_m",0),
                ("initial","core_regularization_m",True),("physics","deep_jet_fraction",-.1),
                ("physics","deep_cyclone_fraction",np.inf),("physics","deep_match_fraction",1),
                ("physics","deep_jet_speed_m_s",1.)]:
            c = copy.deepcopy(base); c[group][key] = value
            with self.assertRaises(ValueError):
                radial.circular_base(c,np.array([0.,1e6]),3.2e6)
        for kwargs in ({"deep_speed":0.},{"match_fraction":1.}):
            with self.assertRaises(ValueError):
                radial.circular_base(base,np.array([0.]),3.2e6,**kwargs)

    def test_profile_audit_pressure_reference_and_deep_matching(self):
        import sys
        sys.path.insert(0,str(ROOT/"scripts"))
        import rostami_audit
        c = self.config()
        metrics,columns = rostami_audit.audit(c,samples=4097)
        self.assertAlmostEqual(metrics["paper_reference_1000_over_f0_s"],3125000.,places=8)
        self.assertGreater(metrics["max_abs_zeta_over_local_f"],.3)
        self.assertLess(metrics["centrifugal_pressure_center_minus_outer_over_f0sq_LDsq"],0)
        self.assertIn("NOT full-square",metrics["radial_pressure_reference"])
        c["physics"].update(deep_jet_fraction=1.,deep_cyclone_fraction=1.)
        matched,_ = rostami_audit.audit(c,samples=4097)
        self.assertEqual(matched["radial_max_abs_chi_relative_to_outer"],0.)
        self.assertEqual(metrics["max_abs_zeta_over_local_f"],matched["max_abs_zeta_over_local_f"])
        # These flags describe the model approximation, not eigenvalue signs.
        self.assertIn("velocity_preserving_qg_analogue_not_rsw_reproduction",metrics["flags"])

    def test_cli_keeps_explicit_component_fractions(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"case.toml"; path.write_text(toml.dumps(self.config()))
            common = [str(path),"--output",str(Path(folder)/"out"),"--radial-nodes","8",
                "--m-range","6:6","--outer-radius","25000000"]
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                radial.main(common+["--match-fractions",".5"])
            with contextlib.redirect_stdout(io.StringIO()):
                rows = radial.main(common)
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]["deep_jet_fraction"],0.)
            self.assertEqual(rows[0]["deep_cyclone_fraction"],0.)
            self.assertEqual(rows[0]["initial_profile"],"rostami2017")


if __name__ == "__main__":
    unittest.main(verbosity=2)

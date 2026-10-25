#!/usr/bin/env python3
"""Audit the printed Ptolemaic map, without constructing a model initializer.

Python >= 3.8; standard library only. Example:
  python scripts/audit_ptolemaic_map.py --output validation/ptolemaic-map-audit.json

Both signs of the printed Omega are evaluated deliberately. Agreement of one
sign with the conventions below does not identify the original simulation's
conventions or authorize changing production model equations.
"""

import argparse
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import sys


@dataclass(frozen=True)
class Parameters:
    n: int = 5
    epsilon: float = 0.05
    k: float = 1.0
    amplitude: float = 1.0
    gamma: float = 10.0
    f0: float = 10.0
    omega: float = 0.004


def quantities(b, p, sign):
    radial = math.exp(p.k * b)
    harmonic = p.epsilon * math.exp(p.n * p.k * b)
    determinant = radial**2 - p.n**2 * harmonic**2
    bracket = radial**2 - p.n * harmonic**2
    frequency = p.omega + sign * p.gamma * p.amplitude**2 * bracket / (2 * p.n)
    derivative = sign * p.gamma * p.amplitude**2 * p.k * determinant / p.n
    return radial, harmonic, determinant, bracket, frequency, derivative


def position(a, b, p):
    """The printed map at t=0, with the physical amplitude retained."""
    radial, harmonic, _, _, _, _ = quantities(b, p, 1)
    theta = p.k * a
    return (p.amplitude * (radial * math.sin(theta) - harmonic * math.sin(p.n * theta)),
            p.amplitude * (radial * math.cos(theta) + harmonic * math.cos(p.n * theta)))


def map_derivatives(a, b, p):
    radial, harmonic, _, _, _, _ = quantities(b, p, 1)
    theta = p.k * a
    factor = p.amplitude * p.k
    xa = factor * (radial * math.cos(theta) - p.n * harmonic * math.cos(p.n * theta))
    ya = factor * (-radial * math.sin(theta) - p.n * harmonic * math.sin(p.n * theta))
    xb = factor * (radial * math.sin(theta) - p.n * harmonic * math.sin(p.n * theta))
    yb = factor * (radial * math.cos(theta) + p.n * harmonic * math.cos(p.n * theta))
    return xa, ya, xb, yb


def inverse_position(x, y, a, b, p):
    """Newton inversion near a supplied label; this uses geometry, never q."""
    for _ in range(30):
        xx, yy = position(a, b, p)
        xa, ya, xb, yb = map_derivatives(a, b, p)
        jac = xa * yb - xb * ya
        if abs(jac) < 1e-14:
            raise ArithmeticError("Map inversion approached a singular Jacobian")
        da = ((x - xx) * yb - xb * (y - yy)) / jac
        db = (xa * (y - yy) - (x - xx) * ya) / jac
        a, b = a + da, b + db
        if abs(da) + abs(db) < 1e-13:
            xx, yy = position(a, b, p)
            if math.hypot(x - xx, y - yy) < 1e-12 * max(1.0, p.amplitude):
                return a, b
    raise ArithmeticError("Cartesian-to-label Newton inversion did not converge")


def velocity(a, b, p, sign):
    # Differentiate theta=ka-W(b)t at fixed material labels, evaluated at t=0.
    frequency = quantities(b, p, sign)[4]
    xa, ya, _, _ = map_derivatives(a, b, p)
    return -frequency * xa / p.k, -frequency * ya / p.k


def analytic_q(a, b, p, sign):
    """Direct map differentiation, not the paper's printed qJ formula."""
    radial, harmonic, determinant, _, frequency, derivative = quantities(b, p, sign)
    cosine = math.cos((p.n + 1) * p.k * a)
    zeta = (2 * frequency * (radial**2 + p.n**3 * harmonic**2)
            + derivative / p.k * (radial**2 + p.n**2 * harmonic**2
                                   - 2 * p.n * radial * harmonic * cosine)) / determinant
    x, y = position(a, b, p)
    return zeta + p.f0 - p.gamma * (x * x + y * y)


def cartesian_fd_q(a, b, p, sign, delta):
    """Independent Cartesian curl from nearby numerically inverted positions."""
    x, y = position(a, b, p)
    vplus = velocity(*inverse_position(x + delta, y, a, b, p), p, sign)[1]
    vminus = velocity(*inverse_position(x - delta, y, a, b, p), p, sign)[1]
    uplus = velocity(*inverse_position(x, y + delta, a, b, p), p, sign)[0]
    uminus = velocity(*inverse_position(x, y - delta, a, b, p), p, sign)[0]
    return (vplus - vminus - uplus + uminus) / (2 * delta) + p.f0 - p.gamma * (x*x + y*y)


def printed_q_unit_amplitude(b, p):
    """Literal equation (3.11), restricted to A0=1 to avoid rescaling ambiguity."""
    if p.amplitude != 1.0:
        return None
    radial, harmonic, determinant, _, _, _ = quantities(b, p, 1)
    return (p.gamma * (p.n-2) * radial**4
            + p.gamma * (1+2*p.n) * p.n**3 * harmonic**4
            - (p.f0+2*p.omega) * p.n * math.exp(2*p.k*p.n*b)
            - (p.f0-2*p.n*p.omega) * p.n**3 * harmonic**2) / (p.n*determinant)


def streamfunction(b, p, sign):
    bracket = quantities(b, p, sign)[3]
    return (p.amplitude**2 * p.omega * bracket / 2
            + sign * p.gamma * p.amplitude**4 * bracket**2 / (8 * p.n))


def audit_case(p, b, deltas):
    labels = [0.0, math.pi / ((p.n + 1) * p.k),
              math.pi / (2 * (p.n + 1) * p.k)]
    records = []
    for sign in (1, -1):
        exact = [analytic_q(a, b, p, sign) for a in labels]
        previous_range = None
        for delta in deltas:
            measured = [cartesian_fd_q(a, b, p, sign, delta) for a in labels]
            spread = max(measured) - min(measured)
            records.append({
                "omega_sign": sign,
                "sign_description": "printed_plus" if sign == 1 else "opposite_sign_audit_only",
                "cartesian_fd_step": delta,
                "q_at_same_b": measured,
                "same_b_pv_range": spread,
                "analytic_q_at_same_b": exact,
                "analytic_same_b_pv_range": max(exact) - min(exact),
                "maximum_analytic_vs_cartesian_fd_error": max(abs(x-y) for x, y in zip(exact, measured)),
                "range_reduction_from_previous_step": previous_range / spread if previous_range is not None and spread else None,
            })
            previous_range = spread

    circular = replace(p, epsilon=0.0)
    circular_records, psi_records = [], []
    for sign in (1, -1):
        radial, _, _, _, frequency, derivative = quantities(b, circular, sign)
        expected = 2 * frequency + derivative / p.k + p.f0 - p.gamma * p.amplitude**2 * radial**2
        circular_records.append({
            "omega_sign": sign,
            "expected_q_2W_plus_Wprime_over_k_plus_f": expected,
            "maximum_analytic_q_error": max(abs(analytic_q(a, b, circular, sign)-expected) for a in labels),
            "maximum_cartesian_fd_q_error": max(abs(cartesian_fd_q(a, b, circular, sign, deltas[-1])-expected) for a in labels),
            "true_particle_angular_rate": frequency,
            "printed_equation_3_13_angular_rate": (1-p.n) * frequency / 2,
            "printed_equation_3_11_q_at_unit_amplitude": printed_q_unit_amplitude(b, circular),
            "printed_q_minus_circular_map_q": (printed_q_unit_amplitude(b, circular)-expected) if circular.amplitude == 1.0 else None,
        })
        _, _, determinant, _, frequency, _ = quantities(b, p, sign)
        expected_psib = p.amplitude**2 * p.k * frequency * determinant
        db = 1e-5
        measured_psib = (streamfunction(b+db, p, sign)-streamfunction(b-db, p, sign)) / (2*db)
        velocity_errors = []
        for a in labels:
            x, y = position(a, b, p)
            delta = deltas[-1]
            bp = inverse_position(x, y+delta, a, b, p)[1]
            bm = inverse_position(x, y-delta, a, b, p)[1]
            u = -(streamfunction(bp, p, sign)-streamfunction(bm, p, sign)) / (2*delta)
            bp = inverse_position(x+delta, y, a, b, p)[1]
            bm = inverse_position(x-delta, y, a, b, p)[1]
            v = (streamfunction(bp, p, sign)-streamfunction(bm, p, sign)) / (2*delta)
            u0, v0 = velocity(a, b, p, sign)
            velocity_errors.append(max(abs(u-u0), abs(v-v0)))
        psi_records.append({"omega_sign": sign, "analytic_psi_b": expected_psib,
                            "finite_difference_psi_b": measured_psib,
                            "psi_b_absolute_error": abs(expected_psib-measured_psib),
                            "maximum_velocity_vs_cartesian_psi_gradient_error": max(velocity_errors)})
    return {"parameters": asdict(p), "b": b, "a_samples": labels,
            "finite_difference_comparisons": records,
            "epsilon_zero_circular_check": circular_records,
            "streamfunction_checks": psi_records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Save JSON; otherwise write JSON to stdout")
    args = parser.parse_args()
    p = Parameters()
    cases = [audit_case(p, -0.1, [1e-3, 1e-4, 1e-5]),
             audit_case(replace(p, amplitude=0.2, k=1.3), -0.1, [2e-4, 2e-5, 2e-6])]
    checks = {}
    for i, case in enumerate(cases):
        plus = case["finite_difference_comparisons"][2]
        minus = case["finite_difference_comparisons"][5]
        checks["case_%d_independent_cartesian_curl_matches_analytic" % i] = max(plus["maximum_analytic_vs_cartesian_fd_error"], minus["maximum_analytic_vs_cartesian_fd_error"]) < 1e-7
        checks["case_%d_printed_sign_has_nonconstant_contour_pv" % i] = plus["same_b_pv_range"] > 1e-3
        checks["case_%d_opposite_sign_contour_pv_converges_to_constant" % i] = minus["same_b_pv_range"] < 1e-7 and case["finite_difference_comparisons"][4]["range_reduction_from_previous_step"] > 50
        checks["case_%d_circular_and_streamfunction_limits" % i] = all(x["maximum_cartesian_fd_q_error"] < 1e-7 for x in case["epsilon_zero_circular_check"]) and all(x["maximum_velocity_vs_cartesian_psi_gradient_error"] < 1e-7 and x["psi_b_absolute_error"] < 1e-7 for x in case["streamfunction_checks"])
    result = {
        "schema_version": 1,
        "status": "convention_reconciliation_pending_original_code",
        "production_initializer_created_or_modified": False,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_version": sys.version.split()[0],
        "source": {"title": "Exact planetary waves and jet streams", "doi": "10.1017/jfm.2024.932",
                   "url": "https://doi.org/10.1017/jfm.2024.932", "map_equations": "3.6-3.10",
                   "suspect_printed_expressions": "3.9 (relative sign), 3.11 (qJ), 3.13 (mean angular rate)",
                   "local_arxiv_tex_lines": "264-307; the published Figure 7 uses epsilon=0.03, the local preprint 0.05"},
        "conventions": {"velocity": "u=-psi_y, v=psi_x", "relative_vorticity": "zeta=v_x-u_y",
                        "planetary_vorticity": "f=f0-gamma*(x*x+y*y)",
                        "phase": "theta=k*a-W(b)*t", "evaluation_time": 0.0,
                        "gamma_scaling": "gamma is the physical coefficient of r^2; retaining A0 introduces gamma*A0^2 in W; project gamma_code=2*gamma",
                        "parameters_are": "synthetic nondimensional audit values, not a reconstruction of a published simulation"},
        "identities": {"definitions": "p=exp(k*b), h=epsilon*exp(n*k*b), D=p^2-n^2*h^2, B=p^2-n*h^2, C=cos((n+1)*theta)",
                       "physical_jacobian": "A0^2*k^2*D",
                       "frequency": "W=omega+sign*gamma*A0^2*B/(2*n)",
                       "vorticity": "zeta=[2*W*(p^2+n^3*h^2)+(Wprime/k)*(p^2+n^2*h^2-2*n*p*h*C)]/D",
                       "pv_cosine_coefficient": "-2*p*h*(n*Wprime/(k*D)+gamma*A0^2)",
                       "contourwise_pv_requires": "Wprime=-gamma*A0^2*k*D/n",
                       "streamfunction": "psi=A0^2*omega*B/2+sign*gamma*A0^4*B^2/(8*n), up to an additive constant",
                       "streamfunction_b_derivative": "psi_b=A0^2*k*W*D",
                       "stationarity": "At fixed b, time only relabels a around a fixed contour; particle angular circulation is not pattern drift."},
        "interpretation": "With the stated map, vorticity and Coriolis conventions, the printed plus sign leaves nonconstant PV along a material contour. The opposite-sign diagnostic removes it. Reconcile conventions with the original author code before creating an initializer; this audit does not establish stability.",
        "cases": cases, "checks": checks, "all_checks_passed": all(checks.values()),
    }
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
        print("Saved %s; %d/%d checks passed; %s" % (args.output, sum(checks.values()), len(checks), result["status"]))
    else:
        sys.stdout.write(encoded)
    if not result["all_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Inviscid circular-base QG stability screen; not a nonlinear hexagon test.

For exp[i*m*(theta-c*t)], solve
    (Omega-c) L_m phi - (q_r/r) phi = 0,
    L_m = d_rr + d_r/r - m**2/r**2 - Ld**-2.
The growth rate is m*Im(c), and Re(c) is angular pattern speed. The upper
base is the configured analytic ring plus cyclone or regularized 2017 fit; the prescribed deep flow
is an annular jet, a fraction of the full circular upper flow, or independently
matched jet and cyclone components. Tapers, convection, sponge and hyperviscosity are absent
from this linear screen. The circular Dirichlet outer boundary is artificial
and differs from the evolution model's periodic continuation.

Example (explicitly small screen):
  python scripts/radial_stability.py configs/unity_followup_01/maintenance_unforced.toml \
    --output analysis/radial-pilot --ld-values 2000000,inf --deep-ratios 0,0.5,1 \
    --m-range 1:16 --radial-nodes 128,256 --outer-radius 24000000
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.linalg import eig, norm
from scipy.optimize import brentq


def ring_profile(r, speed, radius, width):
    """U_theta, U_theta', U_theta'' of the exact untapered Julia ring."""
    r = np.asarray(r)
    a = 1.0 / (radius * width)
    z = (r * r - radius * radius) * a / 2.0
    t = np.tanh(z)
    s = 1.0 - t * t
    k = speed / radius
    u = k * r * s
    up = k * s * (1.0 - 2.0 * a * r * r * t)
    upp = k * ((4.0 * t * t * s - 2.0 * s * s) * a * a * r**3
               - 6.0 * a * r * t * s)
    return u, up, upp


def cyclone_profile(r, speed, radius):
    """U_theta, U_theta', U_theta'' of the exact untapered Julia cyclone."""
    r = np.asarray(r)
    inverse_r2 = 1.0 / radius**2
    factor = (speed / radius) * np.exp(0.5 - r * r * inverse_r2 / 2.0)
    return (factor * r, factor * (1.0 - r * r * inverse_r2),
            factor * (r**3 * inverse_r2**2 - 3.0 * r * inverse_r2))


def _positive_parameter(value, name):
    if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
        raise ValueError(name + " must be positive and finite")
    return float(value)


def _fraction(value, name):
    if isinstance(value, bool) or not np.isfinite(value) or value < 0:
        raise ValueError(name + " must be finite and nonnegative")
    return float(value)


def rostami_parameters(initial):
    """Physical profile parameters; no length or speed is rescaled with LD.

    The Gaussian width is its 1/e half-width. cyclone_radius_m is the peak
    radius of the *unregularized* paper fit. Both speeds specify the isolated
    component's actual peak after polar regularization, not their sum's peak.
    """
    keys = ("jet_radius_m", "jet_width_m", "cyclone_radius_m", "core_regularization_m")
    values = {key: _positive_parameter(initial[key], key) for key in keys}
    for key in ("jet_speed_m_s", "cyclone_speed_m_s"):
        value = initial[key]
        if isinstance(value, bool) or not np.isfinite(value):
            raise ValueError(key + " must be finite")
        values[key] = float(value)
    a, b = 0.42, 1.3
    radius, width = values["jet_radius_m"], values["jet_width_m"]
    length = values["cyclone_radius_m"] / (a/b)**(1/b)
    epsilon = values["core_regularization_m"]
    # Scale roots to O(1); this avoids the absolute-tolerance sensitivity of
    # finding a dimensional radius directly and independently checks Julia.
    def jet_derivative_scaled(x):
        # Equivalent to U_r=0 after factoring positive exponentials.
        return radius/width-x*np.tanh(2*x*radius/width)
    upper = max(radius/width, 1.0)
    while jet_derivative_scaled(upper) > 0:
        upper *= 2
    jet_peak = width*brentq(jet_derivative_scaled, max(0.0, radius/width-1.0), upper, xtol=5e-15)
    shape_peak = (np.exp(-((jet_peak-radius)/width)**2)
                  * -np.expm1(-4*jet_peak*radius/width**2))
    def core_stationary(x):
        r = x*length; s2 = r*r + epsilon*epsilon
        return 1+(a-1)*r*r/s2-b*(r*r/s2)*(np.sqrt(s2)/length)**b
    upper = max(1.0, epsilon/length)
    while core_stationary(upper) > 0:
        upper *= 2
    core_peak = length*brentq(core_stationary, 0.0, upper, xtol=5e-15)
    sc = np.hypot(core_peak, epsilon)
    core_shape_peak = (core_peak/sc)*(sc/length)**a*np.exp(-(sc/length)**b)
    if not (shape_peak > 0 and core_shape_peak > 0):
        raise ValueError("Profile normalization underflow; incompatible profile scales")
    values.update(alpha=a, beta=b, cyclone_length_m=length,
                  jet_peak_radius_m=jet_peak, cyclone_peak_radius_m=core_peak,
                  jet_coefficient=values["jet_speed_m_s"]/shape_peak,
                  cyclone_coefficient=values["cyclone_speed_m_s"]/core_shape_peak)
    return values


def rostami_jet_profile(r, parameters):
    """Odd Gaussian, peak normalized; returns U, U_r, U_rr.

    The reflected Gaussian is a declared smooth-pole correction to the paper
    profile. It is exponentially small near the reference jet, but necessary
    to avoid a nonzero azimuthal velocity at r=0 on the Cartesian grid.
    """
    r = np.asarray(r)
    radius, width = parameters["jet_radius_m"], parameters["jet_width_m"]
    coefficient = parameters["jet_coefficient"]
    minus, plus = (r-radius)/width, (r+radius)/width
    em, ep = np.exp(-minus*minus), np.exp(-plus*plus)
    u = coefficient*em*(-np.expm1(-4*r*radius/width**2))
    up = coefficient*(-2*minus*em+2*plus*ep)/width
    upp = coefficient*((4*minus**2-2)*em-(4*plus**2-2)*ep)/width**2
    return u, up, upp


def rostami_cyclone_profile(r, parameters):
    """Smoothed published power-exponential core, peak normalized."""
    r = np.asarray(r)
    a, b = parameters["alpha"], parameters["beta"]
    length, epsilon = parameters["cyclone_length_m"], parameters["core_regularization_m"]
    s2 = r*r+epsilon*epsilon
    sr = np.sqrt(s2)
    angular = parameters["cyclone_coefficient"]*(sr/length)**a/sr*np.exp(-(sr/length)**b)
    g = (a-1)/s2-b*sr**(b-2)/length**b
    h = -2*(a-1)/s2**2-b*(b-2)*sr**(b-4)/length**b
    return angular*r, angular*(1+r*r*g), angular*r*(3*g+r*r*(g*g+h))


def rostami_angular_profiles(r, parameters):
    """Cancellation-free (U/r, zeta_r/r) for jet and cyclone, including r=0."""
    r = np.asarray(r, dtype=float)
    radius, width = parameters["jet_radius_m"], parameters["jet_width_m"]
    k = 2*radius/width**2
    # A power series for sinh(z)/z and its derivatives in z^2 avoids
    # subtracting singular terms in the vorticity gradient near the pole.
    z = k*r
    small = np.abs(z) < 0.1
    omega_j, gradient_j = np.empty_like(r), np.empty_like(r)
    rs = r[small]; t = z[small]**2
    coeff = [1/math.factorial(2*n+1) for n in range(8)]
    series = np.polynomial.polynomial.polyval(t, coeff)
    first = np.polynomial.polynomial.polyval(t, [n*coeff[n] for n in range(1,8)])
    second = np.polynomial.polynomial.polyval(t, [n*(n-1)*coeff[n] for n in range(2,8)])
    factor = 2*parameters["jet_coefficient"]*k*np.exp(-(rs*rs+radius*radius)/width**2)
    omega_j[small] = factor*series
    gradient_j[small] = factor*(8*(-series/width**2+k*k*first)
        +4*rs*rs*(series/width**4-2*k*k*first/width**2+k**4*second))
    rl = r[~small]
    u, up, upp = rostami_jet_profile(rl, parameters)
    omega_j[~small] = u/rl
    gradient_j[~small] = (upp+up/rl-u/rl**2)/rl
    a, b = parameters["alpha"], parameters["beta"]
    length, epsilon = parameters["cyclone_length_m"], parameters["core_regularization_m"]
    s2 = r*r+epsilon*epsilon; sr = np.sqrt(s2)
    omega_c = parameters["cyclone_coefficient"]*(sr/length)**a/sr*np.exp(-(sr/length)**b)
    g = (a-1)/s2-b*sr**(b-2)/length**b
    h = -2*(a-1)/s2**2-b*(b-2)*sr**(b-4)/length**b
    gradient_c = omega_c*(4*g+r*r*(g*g+h))
    return omega_j, gradient_j, omega_c, gradient_c


def upper_profiles(initial, r):
    """Component U,U_r,U_rr,Omega,zeta_r/r for the configured upper profile."""
    profile = initial.get("profile", "annular")
    if profile == "rostami2017":
        parameters = rostami_parameters(initial)
        jet = rostami_jet_profile(r, parameters)
        core = rostami_cyclone_profile(r, parameters)
        oj, gj, oc, gc = rostami_angular_profiles(r, parameters)
        return (*jet, oj, gj), (*core, oc, gc)
    if profile != "annular":
        raise ValueError("initial.profile must be annular or rostami2017")
    if "core_regularization_m" in initial:
        raise ValueError("core_regularization_m applies only to rostami2017")
    jet = ring_profile(r, initial["jet_speed_m_s"], initial["jet_radius_m"], initial["jet_width_m"])
    core = cyclone_profile(r, initial["cyclone_speed_m_s"], initial["cyclone_radius_m"])
    rj, wj = initial["jet_radius_m"], initial["jet_width_m"]
    alpha = 1/(rj*wj); t = np.tanh((r*r-rj*rj)*alpha/2); s = 1-t*t
    oj = (initial["jet_speed_m_s"]/rj)*s
    gj = (initial["jet_speed_m_s"]/rj)*((4*t*t*s-2*s*s)*alpha**2*r*r-8*alpha*t*s)
    rc = initial["cyclone_radius_m"]
    oc = (initial["cyclone_speed_m_s"]/rc)*np.exp(0.5-r*r/(2*rc**2))
    gc = oc*(r*r/rc**4-4/rc**2)
    return (*jet, oj, gj), (*core, oc, gc)


def circular_base(config, r, deformation_radius, deep_speed=None, match_fraction=None):
    """Analytic base, including cancellation-free q_r/r at the pole.

    Component matching holds the upper wind fixed and changes only prescribed
    deep velocity. Neither LD nor matching fractions rescale the upper fit.
    """
    a, d = config["initial"], config["physics"]
    r = np.asarray(r, dtype=float)
    if np.any(r < 0) or not np.all(np.isfinite(r)):
        raise ValueError("r must contain finite nonnegative radii")
    if isinstance(deformation_radius, bool) or not (deformation_radius > 0):
        raise ValueError("deformation radius must be positive or inf")
    profile = d.get("deep_profile", "annular")
    if profile not in ("annular", "matched_upper", "matched_components"):
        raise ValueError("deep_profile must be annular, matched_upper or matched_components")
    if profile != "matched_components" and any(key in d for key in ("deep_jet_fraction", "deep_cyclone_fraction")):
        raise ValueError("deep component fractions apply only to matched_components")
    if profile == "matched_upper":
        configured_fraction = _fraction(d.get("deep_match_fraction", 1.0), "deep_match_fraction")
        fraction = configured_fraction if match_fraction is None else _fraction(match_fraction, "deep_match_fraction")
    elif "deep_match_fraction" in d or match_fraction is not None:
        raise ValueError("deep_match_fraction applies only to matched_upper")
    if profile in ("matched_upper", "matched_components"):
        if d.get("deep_jet_speed_m_s", 0.0) != 0:
            raise ValueError(profile+" requires deep_jet_speed_m_s=0; use matching fractions")
        if deep_speed is not None:
            raise ValueError("deep_speed cannot override matched deep flow; use matching fractions")
    else:
        deep_speed = d.get("deep_jet_speed_m_s", 0.0) if deep_speed is None else deep_speed
        if isinstance(deep_speed, bool) or not np.isfinite(deep_speed):
            raise ValueError("deep speed must be finite")
    f = 0.0 if math.isinf(deformation_radius) else deformation_radius**-2
    (u, up, upp, omega_ring, grad_ring), (uc, ucp, ucpp, omega_cyclone, grad_cyclone) = upper_profiles(a, r)
    omega = omega_ring+omega_cyclone
    vort_grad_over_r = grad_ring+grad_cyclone
    if profile == "matched_upper":
        ud, omega_deep = fraction*(u+uc), fraction*omega
    elif profile == "matched_components":
        jet_fraction = _fraction(d["deep_jet_fraction"], "deep_jet_fraction")
        core_fraction = _fraction(d["deep_cyclone_fraction"], "deep_cyclone_fraction")
        ud = jet_fraction*u+core_fraction*uc
        omega_deep = jet_fraction*omega_ring+core_fraction*omega_cyclone
    else:
        ud, _, _ = ring_profile(r, deep_speed, d["deep_jet_radius_m"], d["deep_jet_width_m"])
        td = np.tanh((r*r-d["deep_jet_radius_m"]**2)/(2*d["deep_jet_radius_m"]*d["deep_jet_width_m"]))
        omega_deep = (deep_speed/d["deep_jet_radius_m"])*(1-td*td)
    gamma = 2.0*config["planet"]["rotation_rate_s"]/config["planet"]["radius_m"]**2
    qr_over_r = vort_grad_over_r-f*(omega-omega_deep)-gamma
    return {"r": r, "u": u+uc, "u_prime": up+ucp, "u_second": upp+ucpp,
            "u_jet": u, "u_cyclone": uc, "omega_jet": omega_ring, "omega_cyclone": omega_cyclone,
            "u_deep": ud, "omega": omega, "q_r": r*qr_over_r,
            "q_r_over_r": qr_over_r, "vorticity_gradient_over_r": vort_grad_over_r,
            "gamma_code": gamma}


def radial_operator(nodes, outer_radius, mode, deformation_radius):
    """Return r and the symmetric sqrt(r)-transformed R^2 L_m.

    Interior nodes r_j=j*R/(N+1), j=1..N. Fluxes use face radii r_j+/-h/2.
    phi(0)=phi(R)=0; the origin condition selects regular modes for m>=1.
    Multiplication by R^2 conditions the matrix without changing c.
    """
    if not isinstance(nodes, (int, np.integer)) or nodes < 8:
        raise ValueError("nodes must be an integer >= 8")
    if not isinstance(mode, (int, np.integer)) or mode < 1:
        raise ValueError("only regular nonaxisymmetric modes m>=1 are supported")
    if not np.isfinite(outer_radius) or outer_radius <= 0:
        raise ValueError("outer radius must be positive and finite")
    if not (deformation_radius > 0):
        raise ValueError("deformation radius must be positive or inf")
    h = 1.0 / (nodes + 1)
    rho = h * np.arange(1, nodes + 1, dtype=float)
    stretching = 0.0 if math.isinf(deformation_radius) else (outer_radius / deformation_radius)**2
    diagonal = -2.0 / h**2 - mode**2 / rho**2 - stretching
    off_diagonal = (rho[:-1] + h / 2.0) / (h**2 * np.sqrt(rho[:-1] * rho[1:]))
    matrix = np.diag(diagonal) + np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    return rho * outer_radius, matrix


def solve_modes(r, operator, omega, qr_over_r, outer_radius, mode):
    """Generalized eigensolve in y=sqrt(r)*phi; retain every eigenvalue."""
    r = np.asarray(r, dtype=float)
    omega, qr_over_r = np.broadcast_arrays(omega, qr_over_r, r)[:2]
    if not (np.all(np.isfinite(operator)) and np.all(np.isfinite(omega))
            and np.all(np.isfinite(qr_over_r))):
        raise ValueError("the eigenproblem must be finite")
    a = omega[:, None] * operator - np.diag(qr_over_r * outer_radius**2)
    values, vectors = eig(a, operator, check_finite=True)
    if not np.all(np.isfinite(values)):
        raise ArithmeticError("nonfinite generalized eigenvalues")
    leading = int(np.argmax(values.imag))
    value, vector = values[leading], vectors[:, leading]
    av, bv = a.dot(vector), operator.dot(vector)
    residual = norm(av - value * bv) / max(norm(av) + abs(value) * norm(bv), np.finfo(float).tiny)
    phi = vector / np.sqrt(r)
    phi /= max(np.max(np.abs(phi)), np.finfo(float).tiny)
    return {"eigenvalues": values, "leading_index": leading, "leading_c": value,
            "growth_rate_s": mode * value.imag, "angular_pattern_speed_s": value.real,
            "frequency_real_s": mode * value.real, "relative_residual": residual,
            "leading_phi": phi}


def load_config(path):
    try:
        import tomllib
        with open(path, "rb") as stream:
            return tomllib.load(stream)
    except ImportError:
        import toml
        return toml.load(str(path))


def number_list(value, integer=False):
    parser = int if integer else float
    return [parser(part) for part in value.split(",")]


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", required=True, type=Path, help="new analysis directory")
    parser.add_argument("--ld-values", help="comma-separated deformation radii in metres; inf allowed; default from config")
    parser.add_argument("--deep-ratios", help="annular only: deep ring / upper ring speed; default 0,0.5,1")
    parser.add_argument("--match-fractions", help="matched_upper only: comma-separated fractions; default from config")
    parser.add_argument("--m-range", default="1:16", help="inclusive azimuthal mode range FIRST:LAST")
    parser.add_argument("--radial-nodes", default="128,256", help="comma-separated counts of interior radial nodes")
    parser.add_argument("--outer-radius", type=float, default=24000000.0, help="metres; must lie inside or on untapered interior edge")
    parser.add_argument("--growth-tolerance", type=float, default=1e-12, help="classification threshold in s^-1; raw growth is retained")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    lds = number_list(args.ld_values) if args.ld_values else [float(config["physics"]["deformation_radius_m"])]
    profile = config["physics"].get("deep_profile", "annular")
    if profile not in ("annular", "matched_upper", "matched_components"):
        parser.error("unknown deep_profile")
    if profile == "matched_upper":
        if args.deep_ratios is not None:
            parser.error("--deep-ratios is incompatible with matched_upper; use --match-fractions")
        fractions = number_list(args.match_fractions) if args.match_fractions else [config["physics"].get("deep_match_fraction", 1.0)]
        if any(v is None or isinstance(v, bool) or not np.isfinite(v) or v < 0 for v in fractions):
            parser.error("matched fractions must be finite and nonnegative")
        ratios, controls = [], fractions
    elif profile == "matched_components":
        if args.deep_ratios is not None or args.match_fractions is not None:
            parser.error("matched_components uses the two configured fractions exactly; no CLI fraction overrides")
        ratios, fractions, controls = [], [], [None]
    else:
        if args.match_fractions is not None:
            parser.error("--match-fractions requires deep_profile=matched_upper")
        ratios, fractions = number_list(args.deep_ratios or "0,0.5,1"), []
        controls = ratios
    resolutions = number_list(args.radial_nodes, integer=True)
    first, last = (int(x) for x in args.m_range.split(":"))
    if first < 1 or last < first:
        parser.error("m-range must satisfy 1 <= FIRST <= LAST")
    if not (0 < args.outer_radius <= config["boundary"]["inner_radius_m"]):
        parser.error("outer radius must be positive and no larger than boundary.inner_radius_m; this tool omits the taper")
    if not np.isfinite(args.growth_tolerance) or args.growth_tolerance < 0:
        parser.error("growth tolerance must be finite and nonnegative")
    if any(not np.isfinite(ratio) for ratio in ratios) or any(not (ld > 0) for ld in lds):
        parser.error("deep ratios must be finite; deformation radii must be positive or inf")
    if config["initial"]["jet_speed_m_s"] == 0 and any(ratio != 0 for ratio in ratios):
        parser.error("nonzero deep ratios require a nonzero upper ring speed; ratios are defined relative to that coefficient")
    for n in resolutions:
        radial_operator(n, args.outer_radius, first, lds[0])
    args.output.mkdir(parents=True, exist_ok=False)
    source_root = Path(__file__).resolve().parents[1] / "src"
    source_hash = hashlib.sha256("\n".join(p.name + "\n" + p.read_text() for p in sorted(source_root.glob("*.jl"))).encode()).hexdigest()
    metadata = {
        "status": "circular_base_linear_screen_only", "config_path": str(args.config.resolve()),
        "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "julia_source_sha256": source_hash,
        "analysis_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "config": config, "outer_radius_m": args.outer_radius, "radial_nodes": resolutions,
        "deformation_radii_m": lds, "deep_profile": profile, "deep_ratios": ratios,
        "deep_match_fractions": fractions,
        "deep_jet_fraction": config["physics"].get("deep_jet_fraction"),
        "deep_cyclone_fraction": config["physics"].get("deep_cyclone_fraction"),
        "initial_profile": config["initial"].get("profile", "annular"), "mode_range": [first, last],
        "growth_tolerance_s": args.growth_tolerance,
        "method": "Second-order radial flux difference; sqrt(r)-weighted symmetric R^2 L_m; scipy generalized eigenproblem; exp[i*m*(theta-c*t)].",
        "assumptions": [
            "Exact configured untapered upper radial profile; nonaxisymmetric initial disturbance is omitted. rostami2017 includes the declared origin and core regularizations.",
            "Lower layer is prescribed and steady: annular uses a deep/upper ring coefficient, matched_upper multiplies the circular upper flow, and matched_components multiplies jet and cyclone separately.",
            "Regular m>=1 origin condition phi(0)=0 and artificial circular outer Dirichlet condition phi(R)=0.",
            "Convection, sponge and hyperviscosity are absent from this inviscid linear calculation, regardless of configuration values.",
            "This circular outer boundary is not the Julia model's periodic outer continuation. Check outer-radius and radial-grid sensitivity.",
            "Growth of a circular base is not stability of a nonlinear polygon, a matched rotating exact wave, or a claim about Saturn.",
            "Non-growing spectra do not rule out transient growth or finite-amplitude instability. Leading drift and eigenfunction are omitted when growth does not exceed tolerance; all raw eigenvalues are retained.",
        ]}
    (args.output / "metadata.json").write_text(json.dumps(clean_json(metadata), indent=2, allow_nan=False) + "\n")
    (args.output / "config.toml").write_bytes(args.config.read_bytes())
    records, spectra = [], {}
    case_index = 0
    for ld in lds:
        # Deep flow vanishes from unforced barotropic dynamics; do not duplicate it.
        active_controls = [0.0] if math.isinf(ld) else controls
        for control in active_controls:
            speed = control * config["initial"]["jet_speed_m_s"] if profile == "annular" else None
            fraction = control if profile == "matched_upper" else None
            for n in resolutions:
                for mode in range(first, last + 1):
                    r, operator = radial_operator(n, args.outer_radius, mode, ld)
                    base = circular_base(config, r, ld, speed, fraction)
                    result = solve_modes(r, operator, base["omega"], base["q_r_over_r"], args.outer_radius, mode)
                    if result["relative_residual"] > 1e-8:
                        raise ArithmeticError("leading eigenpair residual exceeds 1e-8")
                    key = "case_{:05d}".format(case_index)
                    spectra[key] = result["eigenvalues"]
                    growing = result["growth_rate_s"] > args.growth_tolerance
                    if growing:
                        spectra[key + "_leading_phi"] = result["leading_phi"]
                    row = {"spectrum_key": key, "deformation_radius_m": ld, "deep_profile": profile,
                           "deep_ratio": control if profile == "annular" else None, "deep_match_fraction": fraction,
                           "deep_speed_m_s": speed,
                           "deep_jet_fraction": config["physics"].get("deep_jet_fraction"),
                           "deep_cyclone_fraction": config["physics"].get("deep_cyclone_fraction"),
                           "initial_profile": config["initial"].get("profile", "annular"), "radial_nodes": n, "outer_radius_m": args.outer_radius,
                           "m": mode, "leading_growth_rate_s": result["growth_rate_s"],
                           "leading_angular_pattern_speed_s": result["angular_pattern_speed_s"] if growing else None,
                           "leading_frequency_real_s": result["frequency_real_s"] if growing else None,
                           "leading_relative_residual": result["relative_residual"],
                           "q_r_over_r_min": float(np.min(base["q_r_over_r"])),
                           "q_r_over_r_max": float(np.max(base["q_r_over_r"])),
                           "q_gradient_changes_sign": bool(np.min(base["q_r_over_r"]) < 0 < np.max(base["q_r_over_r"])),
                           "growing_above_tolerance": bool(growing)}
                    records.append(row)
                    case_index += 1
                print("screened Ld={} m, {} control={}, N={}, modes={}:{}".format(ld, profile, control, n, first, last), flush=True)
    with (args.output / "leading_modes.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    np.savez_compressed(str(args.output / "spectra.npz"), **spectra)
    (args.output / "COMPLETED").write_text("{} eigenproblems; circular-base linear screen only\n".format(len(records)))
    return records


if __name__ == "__main__":
    main()

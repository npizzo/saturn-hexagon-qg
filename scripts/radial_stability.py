#!/usr/bin/env python3
"""Inviscid circular-base QG stability screen; not a nonlinear hexagon test.

For exp[i*m*(theta-c*t)], solve
    (Omega-c) L_m phi - (q_r/r) phi = 0,
    L_m = d_rr + d_r/r - m**2/r**2 - Ld**-2.
The growth rate is m*Im(c), and Re(c) is angular pattern speed. The upper
base is the configured analytic ring plus cyclone; the prescribed deep flow
is an annular jet or a fraction of the full circular upper flow, including
the core cyclone. Tapers, convection, sponge and hyperviscosity are absent
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


def circular_base(config, r, deformation_radius, deep_speed=None, match_fraction=None):
    """Analytic base, including a cancellation-free q_r/r at the pole."""
    a, d = config["initial"], config["physics"]
    r = np.asarray(r, dtype=float)
    if np.any(r < 0) or not np.all(np.isfinite(r)):
        raise ValueError("r must contain finite nonnegative radii")
    if not (deformation_radius > 0):
        raise ValueError("deformation radius must be positive or inf")
    profile = d.get("deep_profile", "annular")
    if profile not in ("annular", "matched_upper"):
        raise ValueError("deep_profile must be annular or matched_upper")
    configured_fraction = d.get("deep_match_fraction", 1.0)
    if isinstance(configured_fraction, bool) or not np.isfinite(configured_fraction) or configured_fraction < 0:
        raise ValueError("deep_match_fraction must be finite and nonnegative")
    if profile == "matched_upper":
        if d.get("deep_jet_speed_m_s", 0.0) != 0:
            raise ValueError("matched_upper requires deep_jet_speed_m_s=0; use deep_match_fraction")
        if deep_speed is not None:
            raise ValueError("deep_speed cannot override matched_upper; use match_fraction")
        fraction = configured_fraction if match_fraction is None else match_fraction
        if isinstance(fraction, bool) or fraction is None or not np.isfinite(fraction) or fraction < 0:
            raise ValueError("matched_upper requires a finite nonnegative deep_match_fraction")
    else:
        if "deep_match_fraction" in d:
            raise ValueError("deep_match_fraction applies only to matched_upper")
        if match_fraction is not None:
            raise ValueError("match_fraction applies only to matched_upper")
        deep_speed = d["deep_jet_speed_m_s"] if deep_speed is None else deep_speed
        if not np.isfinite(deep_speed):
            raise ValueError("deep speed must be finite")
    f = 0.0 if math.isinf(deformation_radius) else deformation_radius**-2
    u, up, upp = ring_profile(r, a["jet_speed_m_s"], a["jet_radius_m"], a["jet_width_m"])
    uc, ucp, ucpp = cyclone_profile(r, a["cyclone_speed_m_s"], a["cyclone_radius_m"])
    rj, wj = a["jet_radius_m"], a["jet_width_m"]
    alpha = 1.0 / (rj * wj)
    t = np.tanh((r * r - rj * rj) * alpha / 2.0)
    s = 1.0 - t * t
    omega_ring = (a["jet_speed_m_s"] / rj) * s
    # Analytic (U'' + U'/r - U/r^2)/r; no subtraction of singular terms.
    vort_grad_over_r = (a["jet_speed_m_s"] / rj) * (
        (4.0 * t * t * s - 2.0 * s * s) * alpha**2 * r * r - 8.0 * alpha * t * s)
    rc = a["cyclone_radius_m"]
    omega_cyclone = (a["cyclone_speed_m_s"] / rc) * np.exp(0.5 - r * r / (2.0 * rc**2))
    vort_grad_over_r += omega_cyclone * (r * r / rc**4 - 4.0 / rc**2)
    omega = omega_ring + omega_cyclone
    if profile == "matched_upper":
        ud, omega_deep = fraction * (u + uc), fraction * omega
    else:
        ud, _, _ = ring_profile(r, deep_speed, d["deep_jet_radius_m"], d["deep_jet_width_m"])
        td = np.tanh((r * r - d["deep_jet_radius_m"]**2)
                     / (2.0 * d["deep_jet_radius_m"] * d["deep_jet_width_m"]))
        omega_deep = (deep_speed / d["deep_jet_radius_m"]) * (1.0 - td * td)
    gamma = 2.0 * config["planet"]["rotation_rate_s"] / config["planet"]["radius_m"]**2
    qr_over_r = vort_grad_over_r - f * (omega - omega_deep) - gamma
    return {"r": r, "u": u + uc, "u_prime": up + ucp, "u_second": upp + ucpp,
            "u_deep": ud, "omega": omega, "q_r": r * qr_over_r,
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
    if profile not in ("annular", "matched_upper"):
        parser.error("unknown deep_profile")
    if profile == "matched_upper":
        if args.deep_ratios is not None:
            parser.error("--deep-ratios is incompatible with matched_upper; use --match-fractions")
        fractions = number_list(args.match_fractions) if args.match_fractions else [config["physics"].get("deep_match_fraction", 1.0)]
        if any(v is None or isinstance(v, bool) or not np.isfinite(v) or v < 0 for v in fractions):
            parser.error("matched fractions must be finite and nonnegative")
        ratios, controls = [], fractions
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
        "deep_match_fractions": fractions, "mode_range": [first, last],
        "growth_tolerance_s": args.growth_tolerance,
        "method": "Second-order radial flux difference; sqrt(r)-weighted symmetric R^2 L_m; scipy generalized eigenproblem; exp[i*m*(theta-c*t)].",
        "assumptions": [
            "Exact analytic untapered upper ring plus central cyclone; any configured nonaxisymmetric initial disturbance is omitted.",
            "Lower layer is prescribed and steady: annular uses deep/upper ring coefficient; matched_upper multiplies the entire circular upper flow, including its core cyclone.",
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
                           "deep_speed_m_s": speed, "radial_nodes": n, "outer_radius_m": args.outer_radius,
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

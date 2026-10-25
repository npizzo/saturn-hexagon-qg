#!/usr/bin/env python3
"""Independent saved-field modes, not a jet-shape or stability diagnostic.

Usage: python scripts/field_mode_check.py RESULTS NEW_OUTPUT_DIRECTORY
Requires numpy, scipy, h5py, and toml on Python <3.11. No model is advanced.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import scipy
from scipy.ndimage import map_coordinates
try:
    import tomllib as toml
except ImportError:
    import toml


def radial_inner(a, b, radii):
    """Area-weighted, normalized radial inner product (no angular factor)."""
    return np.trapz(radii * np.conj(a) * b, radii) / np.trapz(radii, radii)


def field_modes(psi, u, v, x, y, radius, width, speed, cap, nr=65, ntheta=720, mmax=32):
    """Accept physical arrays indexed [x,y]; return complex FFT coefficients."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not (np.isfinite([radius, width, speed, cap]).all() and
            width > 0 and speed > 0 and 0 < radius-2*width < radius+2*width < cap):
        raise ValueError("Positive-speed annulus Rj +/- 2w must lie inside the usable cap")
    if nr < 3 or nr % 2 != 1 or mmax < 6 or ntheta <= 2*mmax:
        raise ValueError("Require odd nr >=3, mmax >=6, and ntheta >2*mmax")
    spacing = []
    for a in (x, y):
        if a.ndim != 1 or len(a) < 8 or not np.isfinite(a).all():
            raise ValueError("Coordinates must be finite one-dimensional vectors")
        d = a[1]-a[0]
        if d <= 0 or not np.allclose(np.diff(a), d, rtol=1e-10, atol=abs(d)*1e-10):
            raise ValueError("Coordinates must be strictly increasing and uniform")
        if not (a[2] < -cap and a[-3] > cap):
            raise ValueError("Usable cap must lie inside the Cartesian interpolation domain")
        spacing.append(d)
    radii = np.linspace(radius-2*width, radius+2*width, nr)
    theta = 2*np.pi*np.arange(ntheta)/ntheta  # No duplicated periodic endpoint.
    coords = np.array([(radii[:, None]*np.cos(theta)-x[0])/spacing[0],
                       (radii[:, None]*np.sin(theta)-y[0])/spacing[1]])
    sampled = []
    for f in (psi, u, v):
        f = np.asarray(f)
        if f.shape != (len(x), len(y)) or not np.isrealobj(f) or not np.isfinite(f).all():
            raise ValueError("Fields must be finite real arrays with physical [x,y] shape")
        sampled.append(map_coordinates(f, coords, order=3, mode="wrap", prefilter=True))
    psi_ring = sampled[0]
    ut = -sampled[1]*np.sin(theta)+sampled[2]*np.cos(theta)
    def coeff(a):
        return np.fft.rfft(a-a.mean(axis=1, keepdims=True), axis=1)[:, 1:mmax+1]/ntheta
    return {"radii": radii, "psi": coeff(psi_ring), "utheta": coeff(ut),
            "jet_speed": speed, "psi_scale": float(np.max(np.abs(psi_ring)))}


def summarize_modes(measured, reference):
    """Projection phase is a wrapped field-coefficient phase, never a speed."""
    r = measured["radii"]
    if not np.array_equal(r, reference["radii"]):
        raise ValueError("Reference and current radial grids differ")
    phase = None
    a, b = reference["psi"][:, 5], measured["psi"][:, 5]
    na = np.sqrt(max(0.0, radial_inner(a, a, r).real))
    nb = np.sqrt(max(0.0, radial_inner(b, b, r).real))
    projection = radial_inner(a, b, r)
    eps = 64*np.finfo(float).eps
    if (na > eps*reference["psi_scale"] and nb > eps*measured["psi_scale"] and
            abs(projection) > eps*na*nb):
        phase = float(np.angle(projection))
    rows = []
    for j in range(measured["psi"].shape[1]):
        p, velocity = measured["psi"][:, j], measured["utheta"][:, j]
        rows.append({"mode": j+1,
            "psi_mode_radial_rms_m2_s": float(np.sqrt(max(0.0, radial_inner(p, p, r).real))),
            "equivalent_displacement_at_Rj_m": float(2*abs(p[len(r)//2])/measured["jet_speed"]),
            "utheta_mode_radial_rms_m_s": float(np.sqrt(max(0.0, radial_inner(velocity, velocity, r).real))),
            "m6_projection_phase_rad": phase if j == 5 else None})
    return rows


def read_snapshot(path):
    """h5py reverses Julia matrix dimensions: transpose to physical [x,y]."""
    with h5py.File(str(path), "r") as f:
        t = float(f["time_s"][()])
        fields = [np.array(f[key]).T for key in ("psi", "u", "v")]
    if not np.isfinite(t):
        raise ValueError("Nonfinite saved time: " + str(path))
    return t, fields


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(results, output):
    results, output = Path(results).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("Output must be a NEW directory")
    config = toml.loads((results/"config.toml").read_text())
    metadata = toml.loads((results/"metadata.toml").read_text())
    units = metadata.get("units", "")
    if not all(s in units for s in ("planar distance in m", "velocity m/s", "psi m2/s")):
        raise ValueError("Saved metadata does not identify expected SI coordinate/field units")
    if "arrays indexed [x,y]" not in metadata.get("axes", ""):
        raise ValueError("Saved metadata does not identify Julia [x,y] orientation")
    with h5py.File(str(results/"domain.jld2"), "r") as f:
        x, y = np.array(f["x"]), np.array(f["y"])
    n, length = config["grid"]["nx"], config["grid"]["length_m"]
    expected = -length/2+np.arange(n)*length/n
    for a in (x, y):
        if a.shape != (n,) or not np.allclose(a, expected, rtol=1e-10, atol=length*1e-12):
            raise ValueError("Coordinates disagree with configured SI periodic grid (endpoint excluded)")
    files = list(results.glob("snapshot_[0-9]*.jld2"))
    records = []
    for path in files:
        with h5py.File(str(path), "r") as f:
            records.append((float(f["time_s"][()]), path))
    records.sort(key=lambda item: item[0])
    times = np.array([t for t, _ in records])
    if not len(times) or not np.isfinite(times).all() or not np.all(np.diff(times) > 0):
        raise ValueError("Require snapshots with distinct, finite increasing times")
    initial, boundary = config["initial"], config["boundary"]
    rows, reference, inputs = [], None, []
    for _, path in records:
        time, fields = read_snapshot(path)
        modes = field_modes(*fields, x, y, initial["jet_radius_m"], initial["jet_width_m"],
                            initial["jet_speed_m_s"], boundary["inner_radius_m"])
        if reference is None:
            reference = modes
        rows.extend(dict(snapshot=path.name, time_s=time, **row) for row in summarize_modes(modes, reference))
        stat = path.stat()
        inputs.append({"file": path.name, "time_s": time, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    info = {"schema_version": 1, "results": str(results), "snapshots": inputs,
        "sha256": {name: digest(results/name) for name in ("config.toml", "metadata.toml", "domain.jld2")},
        "script_sha256": digest(Path(__file__)), "reported_source_sha256": metadata.get("source_sha256"),
        "reported_manifest_sha256": metadata.get("manifest_sha256"),
        "libraries": {"numpy": np.__version__, "scipy": scipy.__version__, "h5py": h5py.__version__},
        "sampling": {"nr": 65, "ntheta": 720, "m_max": 32, "annulus_m": [float(reference["radii"][0]), float(reference["radii"][-1])],
                     "interpolation": "scipy order-3 spline, mode=wrap, prefilter=True; annulus inside usable cap"},
        "definition": "c_m(r)=mean_theta[(field-azimuthal_mean) exp(-i m theta)]; normalized radial L2 = sqrt(integral r|c_m|^2 dr / integral r dr); coefficient norms are not doubled",
        "displacement": "2|psi_m(Rj)|/Uj is a field-amplitude proxy, not a measured jet displacement",
        "phase": "m=6 only: arg integral r conj(psi_6(first)) psi_6(t) dr, in [-pi,pi]; empty if reference/current/projection numerically zero; no unwrapping or speed inference",
        "limitations": "No jet connectivity, stability, modal-growth significance, or interpolation convergence established. First saved snapshot is the reference, including resumed runs. Snapshot files identified by path/time/size/mtime, not content hashes."}
    output.mkdir(parents=True, exist_ok=False)
    with (output/"fields.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (output/"metadata.json").write_text(json.dumps(info, indent=2, allow_nan=False)+"\n")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results"); parser.add_argument("output")
    args = parser.parse_args()
    rows = analyze(args.results, args.output)
    print("Wrote %d field-mode records; phases are wrapped and no speeds are inferred." % len(rows))

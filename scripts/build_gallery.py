#!/usr/bin/env python3
"""Build a portable static gallery from completed SaturnPolarQG result directories.

This reads saved fields only: it never advances the model. Python 3.8+, numpy,
h5py, matplotlib, and toml (only on Python <3.11) are required. Optional MP4
uses ffmpeg and repeats actual saved frames, without temporal interpolation.
"""
import argparse
import csv
import hashlib
import html
import io
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import tomllib
    parse_toml = tomllib.loads
except ImportError:
    import toml
    parse_toml = toml.loads


SCHEMA = {
    "experiment": {"name", "purpose", "seed"},
    "grid": {"nx", "length_m"},
    "planet": {"radius_m", "rotation_rate_s"},
    "physics": {"deformation_radius_m", "deep_jet_speed_m_s", "deep_jet_radius_m", "deep_jet_width_m", "deep_profile", "deep_match_fraction"},
    "initial": {"jet_speed_m_s", "jet_radius_m", "jet_width_m", "cyclone_speed_m_s", "cyclone_radius_m", "perturbation", "displacement_m", "mode", "m_max"},
    "boundary": {"inner_radius_m", "outer_radius_m", "sponge_rate_s"},
    "numerics": {"dt_s", "duration_s", "hyper_order", "cutoff_damping_time_s", "cfl_limit"},
    "convection": {"enabled", "r_min_m", "r_max_m", "vortex_radius_m", "q_max_s", "bias_max", "psi_limit_m2_s", "interval_min_s", "interval_max_s"},
    "output": {"diagnostic_interval_s", "snapshot_interval_s", "checkpoint_interval_s", "m_max", "ntheta", "nr", "jet_search_halfwidth_m"},
}
CORE_COLUMNS = ["time_s", "step", "energy_m2_s2", "total_pv_enstrophy_s2", "mean_s_s",
                "mean_jet_radius_m", "mean_peak_speed_m_s", "jet_coverage"]
MODE_COLUMN = re.compile(r"m([1-9][0-9]*)_(amplitude|real|imag)_m$")
HEX_HASH = re.compile(r"[a-fA-F0-9]{64}$")
SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9 _.-]{0,99}$")
FIELD_KEYS = ("s", "psi", "u", "v", "totalpv")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def public_name(value, fallback):
    """Do not publish paths or arbitrary configuration strings as names."""
    return value if isinstance(value, str) and SAFE_NAME.fullmatch(value) and ".." not in value else fallback


def safe_config(config):
    result = {}
    for section, keys in SCHEMA.items():
        values = {}
        for key in sorted(keys):
            if key not in config.get(section, {}):
                continue
            value = config[section][key]
            if key == "name":
                value = public_name(value, "unnamed-run")
            elif key == "deep_profile":
                if value not in {"annular", "matched_upper"}:
                    raise ValueError("Unknown deep-flow profile")
            elif key in ("purpose", "perturbation"):
                if value not in {"maintenance", "selection", "emergence", "smoke", "validation", "screen", "benchmark", "none", "mode", "broadband"}:
                    continue
            elif not isinstance(value, (int, float, bool)):
                raise ValueError("Expected numeric/boolean configuration field " + section + "." + key)
            elif isinstance(value, float) and math.isnan(value):
                raise ValueError("NaN configuration field " + section + "." + key)
            values[key] = value
        result[section] = values
    return result


def toml_text(config):
    """Write the restricted scalar schema without copying source comments."""
    lines = ["# Public configuration: recognized model scalars only; source comments omitted."]
    for section, values in config.items():
        lines.extend(["", "[" + section + "]"])
        for key, value in values.items():
            if isinstance(value, bool):
                encoded = "true" if value else "false"
            elif isinstance(value, str):
                encoded = json.dumps(value)
            elif isinstance(value, float) and math.isinf(value):
                encoded = "inf" if value > 0 else "-inf"
            else:
                encoded = repr(value)
            lines.append(key + " = " + encoded)
    return "\n".join(lines) + "\n"


def increasing_times(values, label):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all() or np.any(values < 0) or np.any(np.diff(values) <= 0):
        raise ValueError(label + " must be finite, nonnegative, strictly increasing times")
    return values


def read_diagnostics(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Missing or duplicate diagnostic columns")
        # Do not copy arbitrary extra CSV columns or string-valued provenance.
        columns = [k for k in reader.fieldnames if k in CORE_COLUMNS or MODE_COLUMN.fullmatch(k)]
        if any(k not in columns for k in CORE_COLUMNS):
            raise ValueError("Missing required diagnostic columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Empty diagnostics")
    data = {key: np.array([float(row[key]) for row in rows]) for key in columns}
    increasing_times(data["time_s"], "Diagnostics")
    for key in CORE_COLUMNS:
        if not np.isfinite(data[key]).all():
            raise ValueError("Nonfinite core diagnostic: " + key)
    if np.any((data["jet_coverage"] < 0) | (data["jet_coverage"] > 1)):
        raise ValueError("Jet coverage outside [0,1]")
    for key in columns:
        if MODE_COLUMN.fullmatch(key):
            if np.isinf(data[key]).any():
                raise ValueError("Infinite mode diagnostic: " + key)
            if np.any(np.isnan(data[key]) & (data["jet_coverage"] >= 1)):
                raise ValueError("Missing mode despite full tracker coverage: " + key)
            if "amplitude" in key and np.any(data[key][np.isfinite(data[key])] < 0):
                raise ValueError("Negative mode amplitude")
    return data


def read_snapshot(path, shape):
    with h5py.File(str(path), "r") as f:
        t = float(f["time_s"][()])
        # JLD2/HDF5 reverses the Julia dimensions as seen by h5py.
        fields = {key: np.asarray(f[key]).T for key in FIELD_KEYS}
    if not math.isfinite(t) or t < 0:
        raise ValueError("Nonfinite or negative snapshot time")
    for key, value in fields.items():
        if value.shape != shape or not np.isrealobj(value) or not np.isfinite(value).all():
            raise ValueError("Invalid snapshot field " + key)
    return t, fields


def vorticity(fields, deformation_radius):
    if deformation_radius <= 0 or math.isnan(deformation_radius):
        raise ValueError("Deformation radius must be positive")
    return fields["s"] + fields["psi"] / deformation_radius**2


def display_array(field, maximum_side=512):
    """Block-average [x,y] cells for display only, then put y on image rows."""
    n, m = field.shape
    factor = max(1, int(math.ceil(max(n, m)/maximum_side)))
    # Avoid discarding edge cells or shifting the Cartesian origin.
    while factor > 1 and (n % factor or m % factor):
        factor += 1
        if factor > min(n, m):
            raise ValueError("Grid cannot be evenly block-averaged for display")
    reduced = field.reshape(n//factor, factor, m//factor, factor).mean(axis=(1, 3))
    return reduced.T, factor


def source_digest(source):
    files = sorted(Path(source).glob("*.jl"))
    if not files:
        return None
    text = "\n".join(p.name + "\n" + p.read_text() for p in files)
    return hashlib.sha256(text.encode()).hexdigest()


def validation_report(path, config, integrity, config_hash, diagnostics):
    """Verify a short CPU/GPU report and expose only recognized numerical fields.

    Smoke integration completion alone cannot establish that its subsequent
    CPU/GPU comparison passed. The saved validator, inputs and report must agree.
    """
    marker = path / "CPU_GPU_VALIDATION_PASSED"
    report_path = path / "cpu_gpu_validation.toml"
    validator = path / "validate_stabilization.jl"
    if not marker.is_file() or marker.read_text() != "Short-time CPU/GPU validation passed.\n" or not report_path.is_file() or not validator.is_file():
        raise ValueError("Smoke gallery requires a completed CPU/GPU validation report and saved validator")
    report = parse_toml(report_path.read_text())
    if type(report.get("schema_version")) is not int or report["schema_version"] != 1 or report.get("passed") is not True:
        raise ValueError("CPU/GPU validation did not pass or has an unsupported schema")
    if not integrity["source_verified"] or not integrity["manifest_verified"]:
        raise ValueError("CPU/GPU validation requires saved source and dependency manifest")
    hashes = {"source_sha256": integrity["source_sha256"],
              "manifest_sha256": integrity["manifest_sha256"],
              "config_sha256": config_hash, "validator_sha256": sha256(validator)}
    for key, expected in hashes.items():
        if report.get(key) != expected:
            raise ValueError("CPU/GPU validation " + key + " does not match saved inputs")

    def finite_number(value, label, nonnegative=True):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or (nonnegative and value < 0):
            raise ValueError("Invalid CPU/GPU validation " + label)
        return float(value)

    def table(parent, key):
        value = parent.get(key)
        if not isinstance(value, dict):
            raise ValueError("Missing or invalid CPU/GPU validation table " + key)
        return value

    tolerance = finite_number(report.get("relative_l2_tolerance"), "relative L2 tolerance")
    if not 0 < tolerance <= 1e-8:
        raise ValueError("CPU/GPU validation relative L2 tolerance exceeds the supported bound")
    errors = {}
    for phase in ("initial_errors", "final_errors"):
        supplied = table(report, phase)
        errors[phase] = {}
        for field in FIELD_KEYS:
            values = table(supplied, field)
            relative = finite_number(values.get("relative_l2"), phase + " " + field)
            maximum = finite_number(values.get("maximum_absolute"), phase + " " + field + " maximum absolute error")
            if relative > tolerance:
                raise ValueError("CPU/GPU validation field error exceeds tolerance")
            errors[phase][field] = {"relative_l2": relative, "maximum_absolute": maximum}

    duration = float(config["numerics"]["duration_s"])
    clock = table(report, "final_clock")
    time = finite_number(clock.get("time_s"), "final time")
    if not 0 < duration <= 20000 or not math.isclose(time, duration, rel_tol=1e-12, abs_tol=1e-7):
        raise ValueError("CPU/GPU validation final time differs from the short saved run")
    steps, events = clock.get("step"), clock.get("events")
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in (steps, events)) or steps != diagnostics["step"][-1]:
        raise ValueError("CPU/GPU validation step/event count is invalid")
    convection = config["convection"]["enabled"]
    if report.get("convection_enabled") is not convection or (not convection and events != 0):
        raise ValueError("CPU/GPU validation convection state differs from the configuration")
    physics = config["physics"]
    deep_profile = physics.get("deep_profile", "annular")
    if report.get("deep_profile") != deep_profile:
        raise ValueError("CPU/GPU validation deep profile differs from the configuration")
    applies = deep_profile == "matched_upper" and physics.get("deep_match_fraction", 1.0) == 1 and math.isfinite(physics["deformation_radius_m"])
    base = table(report, "full_match_base")
    if base.get("applicable") is not applies:
        raise ValueError("CPU/GPU full-match base check applicability differs from the configuration")
    public_base = {"applicable": applies}
    if applies:
        base_tolerance = finite_number(base.get("absolute_thickness_tolerance"), "base thickness tolerance")
        if not 0 < base_tolerance <= 1e-11:
            raise ValueError("CPU/GPU full-match base tolerance exceeds the supported bound")
        public_base["absolute_thickness_tolerance"] = base_tolerance
        for device in ("cpu", "gpu"):
            metrics = table(base, device)
            maximum = finite_number(metrics.get("max_abs_hprime_over_H"), "base thickness magnitude")
            minimum = finite_number(metrics.get("min_h_over_H"), "base minimum thickness", nonnegative=False)
            if maximum >= base_tolerance or abs(minimum-1) >= base_tolerance:
                raise ValueError("CPU/GPU full-match base thickness cancellation failed")
            public_base[device] = {"max_abs_hprime_over_H": maximum, "min_h_over_H": minimum}
    return dict(schema_version=1, passed=True, relative_l2_tolerance=tolerance,
                input_report_sha256=sha256(report_path), **hashes, **errors,
                final_time_s=time, steps=steps, convection_events=events,
                full_match_base=public_base,
                scope="Short-time CPU/GPU agreement and applicable full-match base thickness cancellation. This does not test hexagon persistence, spatial convergence, or long-time stability.")


def load_run(path):
    path = Path(path).resolve()
    if not (path / "COMPLETED").is_file() or any((path / name).exists() for name in ("FAILED", "FAILED.txt")):
        raise ValueError("Gallery requires a completed result directory")
    c = parse_toml((path / "config.toml").read_text())
    meta = parse_toml((path / "metadata.toml").read_text())
    public = safe_config(c)
    n = int(c["grid"]["nx"])
    length = float(c["grid"]["length_m"])
    cap = float(c["boundary"]["inner_radius_m"])
    ld = float(c["physics"]["deformation_radius_m"])
    if n < 8 or not math.isfinite(length) or not 0 < cap < length/2 or ld <= 0 or math.isnan(ld):
        raise ValueError("Invalid domain or deformation radius")
    data = read_diagnostics(path / "diagnostics.csv")
    completed = float((path / "COMPLETED").read_text().strip())
    duration = float(c["numerics"]["duration_s"])
    if not math.isfinite(completed) or not math.isclose(completed, duration, abs_tol=1e-7, rel_tol=1e-12) or not math.isclose(data["time_s"][-1], completed, abs_tol=1e-7, rel_tol=1e-12):
        raise ValueError("Completion marker, requested duration and diagnostics disagree")
    with h5py.File(str(path / "domain.jld2"), "r") as f:
        x, y = np.asarray(f["x"]), np.asarray(f["y"])
    for axis in (x, y):
        if axis.shape != (n,) or not np.isfinite(axis).all() or np.any(np.diff(axis) <= 0) or not np.allclose(np.diff(axis), length/n, rtol=1e-10, atol=1e-6):
            raise ValueError("Invalid Cartesian coordinates")
        if not axis[0] <= -cap < cap <= axis[-1]:
            raise ValueError("Science disk lies outside saved coordinates")
    snapshots = sorted(path.glob("snapshot_*.jld2"))
    if not snapshots or any(not re.fullmatch(r"snapshot_[0-9]+\.jld2", s.name) for s in snapshots):
        raise ValueError("Missing or unrecognized snapshot files")
    mask = np.hypot(x[:, None], y[None, :]) <= cap
    frames, zmax, vmax = [], 0.0, 0.0
    for snap in snapshots:
        time, fields = read_snapshot(snap, (n, n))
        z = vorticity(fields, ld)
        speed = np.hypot(fields["u"], fields["v"])
        zmax = max(zmax, float(np.max(np.abs(z[mask]))))
        vmax = max(vmax, float(np.max(speed[mask])))
        frames.append({"file": snap.name, "time_s": time, "sha256": sha256(snap)})
    increasing_times([s["time_s"] for s in frames], "Snapshots")
    if not math.isclose(frames[0]["time_s"], data["time_s"][0], abs_tol=1e-7, rel_tol=1e-12) or not math.isclose(frames[-1]["time_s"], completed, abs_tol=1e-7, rel_tol=1e-12):
        raise ValueError("Snapshots do not span the diagnostic interval")
    integrity = {}
    for key, actual in (("source_sha256", source_digest(path / "src")),
                        ("manifest_sha256", sha256(path / "Manifest.toml") if (path / "Manifest.toml").exists() else None)):
        expected = meta.get(key)
        if not isinstance(expected, str) or not HEX_HASH.fullmatch(expected):
            raise ValueError("Missing or invalid " + key)
        if actual is not None and actual != expected:
            raise ValueError("Saved " + key + " does not match metadata")
        integrity[key] = expected
        integrity[key.replace("sha256", "verified")] = actual is not None
    config_hash = sha256(path / "config.toml")
    validation = validation_report(path, c, integrity, config_hash, data) if c.get("experiment", {}).get("purpose") == "smoke" else None
    name = public_name(c.get("experiment", {}).get("name"), "run-" + config_hash[:8])
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-") + "-" + config_hash[:8]
    return dict(path=path, config=c, public_config=public, name=name, slug=slug,
                data=data, x=x, y=y, cap=cap, ld=ld, frames=frames,
                zmax=zmax, vmax=vmax, integrity=integrity,
                validation=validation,
                input_config_sha256=config_hash, input_diagnostics_sha256=sha256(path / "diagnostics.csv"))


def render_frame(run, frame, output, zlim, vlim):
    time, fields = read_snapshot(run["path"] / frame["file"], (len(run["x"]), len(run["y"])))
    cap = run["cap"]/1e6
    x, y = run["x"], run["y"]
    dx, dy = x[1]-x[0], y[1]-y[0]
    extent = [(x[0]-dx/2)/1e6, (x[-1]+dx/2)/1e6, (y[0]-dy/2)/1e6, (y[-1]+dy/2)/1e6]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.65), dpi=100)
    fig.subplots_adjust(left=.065, right=.975, bottom=.14, top=.82, wspace=.32)
    fig.suptitle("%s  |  %.3f days (%s s)" % (run["name"], time/86400, format(time, ".8g")), fontsize=12, fontweight="bold", y=.98)
    arrays = (vorticity(fields, run["ld"])*1e5, np.hypot(fields["u"], fields["v"]))
    styles = (("Relative vorticity", "RdBu_r", -zlim*1e5, zlim*1e5, "10⁻⁵ s⁻¹"),
              ("Horizontal speed", "viridis", 0, vlim, "m s⁻¹"))
    for ax, value, (label, cmap, low, high, unit) in zip(axes, arrays, styles):
        shown, factor = display_array(value)
        xx = x.reshape(-1, factor).mean(axis=1)
        yy = y.reshape(-1, factor).mean(axis=1)
        shown = np.ma.masked_where(np.hypot(xx[None, :], yy[:, None]) > run["cap"], shown)
        im = ax.imshow(shown, origin="lower", extent=extent, cmap=cmap, vmin=low, vmax=high, interpolation="nearest")
        ax.set(xlim=(-cap, cap), ylim=(-cap, cap), xlabel="x (Mm)", ylabel="y (Mm)", title=label)
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, shrink=.82, pad=.025, label=unit)
    fig.text(.5, .025, "Saved frame; no temporal interpolation · Common color limits · %d×%d display cells" % (len(x)//factor, len(y)//factor), ha="center", fontsize=8, color="#53606c")
    fig.savefig(str(output), dpi=100)
    plt.close(fig)


def render_diagnostics(run, output, mode_max_m=None):
    d = run["data"]
    days = d["time_s"]/86400
    modes = sorted(int(MODE_COLUMN.fullmatch(k).group(1)) for k in d if MODE_COLUMN.fullmatch(k) and "amplitude" in k)
    if not modes:
        raise ValueError("No mode amplitudes in diagnostics")
    amps = np.array([d["m%d_amplitude_m" % m] for m in modes])/1000
    display_floor_km = .01  # Readability limit, not a scientific detection threshold.
    finite_amplitudes = amps[np.isfinite(amps)]
    common_max_km = max(display_floor_km*10, mode_max_m/1000 if mode_max_m is not None else float(finite_amplitudes.max()) if finite_amplitudes.size else display_floor_km*10)
    def visible(a):
        return np.where(np.nan_to_num(a) >= display_floor_km, a, np.nan)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7), dpi=100)
    fig.subplots_adjust(left=.095, right=.98, bottom=.09, top=.91, hspace=.44, wspace=.32)
    fig.suptitle("Measured jet structure and budgets", fontsize=14, fontweight="bold")
    ax = axes[0, 0]
    competitors = sorted((float(np.nanmax(a)) if np.isfinite(a).any() else -1, m) for a, m in zip(amps, modes) if m != 6)[-3:]
    selected = {m for _, m in competitors} | {6}
    for m, a in zip(modes, amps):
        if m not in selected:
            ax.plot(days, visible(a), color="#b5c0c8", lw=.6, alpha=.6)
    for m in sorted(selected):
        if m in modes:
            a = amps[modes.index(m)]
            ax.plot(days, visible(a), lw=2 if m == 6 else 1.1,
                    color="#bd4c35" if m == 6 else None, label="m = %d" % m)
    amplitude_title = "Sixfold mode and strongest competitors" if 6 in modes else "Measured mode amplitudes"
    ax.set(yscale="log", ylim=(display_floor_km,common_max_km*1.3), xlabel="Simulated time (days)", ylabel="Centreline amplitude (km)", title=amplitude_title)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax = axes[0, 1]
    if len(days) > 1:
        edges = np.r_[days[0]-(days[1]-days[0])/2, (days[:-1]+days[1:])/2, days[-1]+(days[-1]-days[-2])/2]
    else:
        edges = np.array([days[0]-.5, days[0]+.5])
    # Preserve missing modes/measurements as masked cells, never numerical zero.
    full = np.full((max(modes), len(days)), np.nan)
    for m, a in zip(modes, amps):
        good = np.nan_to_num(a) >= display_floor_km
        full[m-1, good] = np.log10(a[good])
    im = ax.pcolormesh(edges, np.arange(.5, max(modes)+1.5), np.ma.masked_invalid(full), cmap="viridis", vmin=math.log10(display_floor_km), vmax=math.log10(common_max_km), rasterized=True)
    ax.set(xlabel="Simulated time (days)", ylabel="Azimuthal mode m", title="Resolved displacement spectrum", ylim=(.5,max(modes)+.5))
    fig.colorbar(im, ax=ax, label="log₁₀ amplitude (km)", pad=.02)
    ax = axes[1, 0]
    ax.plot(days, d["jet_coverage"], color="#216b83", lw=1.8)
    ax.set(xlabel="Simulated time (days)", ylabel="Fraction of valid angles", title="Jet tracker coverage", ylim=(-.03,1.05))
    ax.text(.02,.05,"Full coverage alone does not rule out ridge switching.", transform=ax.transAxes, fontsize=7, color="#53606c")
    ax = axes[1, 1]
    for key, label, color in (("energy_m2_s2", "QG energy", "#216b83"), ("total_pv_enstrophy_s2", "Total-PV enstrophy", "#bd4c35")):
        if d[key][0] > 0:
            ax.plot(days,100*(d[key]/d[key][0]-1), label=label, color=color)
    ax.set(xlabel="Simulated time (days)", ylabel="Change from initial (%)", title="Total changes, including forcing and dissipation")
    ax.set_title("Total changes, including forcing and dissipation", pad=18)
    ax.legend(fontsize=8, frameon=False)
    for ax in axes.flat:
        ax.grid(alpha=.13)
        if len(days)>1:
            ax.set_xlim(days[0],days[-1])
    fig.text(.5,.012,"Mode display floor: 10 m for readability; not a detection threshold. Full values remain in the CSV.",ha="center",fontsize=7,color="#53606c")
    fig.savefig(str(output), dpi=100)
    plt.close(fig)


def format_number(value, unit=""):
    return ("∞" if math.isinf(float(value)) else format(float(value), ".5g")) + (" " + unit if unit else "")


def deep_flow_label(physics):
    if physics.get("deep_profile", "annular") == "matched_upper":
        return "Whole upper base × " + format_number(physics.get("deep_match_fraction", 1.0)) + " (jet + cyclone)"
    return "Annular jet " + format_number(physics["deep_jet_speed_m_s"], "m/s")


def metadata_rows(run):
    c = run["public_config"]
    p, initial, numerical = c["physics"], c["initial"], c["numerics"]
    if initial["perturbation"] == "none":
        disturbance = "None"
    elif initial["perturbation"] == "broadband":
        count = initial["m_max"]
        disturbance = "Broadband; normalization d = %s km; per-mode parameter d/√M = %s km (M = %s); seed %s" % (
            format_number(initial["displacement_m"]/1000),
            format_number(initial["displacement_m"]/(1000*math.sqrt(count))), count, c["experiment"]["seed"])
    else:
        disturbance = "%s; %s km; seed %s" % (initial["perturbation"], format_number(initial["displacement_m"]/1000), c["experiment"]["seed"])
    rows = [
        ("Convection", "On" if c["convection"]["enabled"] else "Off"),
        ("Deformation radius", format_number(p["deformation_radius_m"]/1e6,"Mm")),
        ("Upper annular jet", format_number(initial["jet_speed_m_s"], "m/s")),
        ("Prescribed deep flow", deep_flow_label(p)),
        ("Jet radius / width", "%s / %s Mm" % (format_number(initial["jet_radius_m"]/1e6), format_number(initial["jet_width_m"]/1e6))),
        ("Initial disturbance", disturbance),
        ("Initial azimuthal mode(s)", str(initial["mode"]) if initial["perturbation"] == "mode" else "1–%s" % initial["m_max"] if initial["perturbation"] == "broadband" else "None"),
        ("Grid / timestep", "%s² / %s s" % (c["grid"]["nx"], format_number(numerical["dt_s"]))),
        ("Saved time span", "%s–%s days" % (format_number(run["frames"][0]["time_s"]/86400), format_number(run["frames"][-1]["time_s"]/86400))),
        ("Actual snapshots", str(len(run["frames"]))),
        ("Rotation rate", format_number(c["planet"]["rotation_rate_s"],"s⁻¹ (model frame)")),
    ]
    if run.get("validation"):
        rows.insert(0, ("Run purpose", "Short CPU/GPU software validation"))
        rows.insert(1, ("Compared duration", format_number(run["validation"]["final_time_s"], "s")))
    return rows


def page(title, body, prefix):
    return '<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><link rel="stylesheet" href="'+prefix+'assets/style.css"></head><body>'+body+'<footer>Saturn polar QG · Saved simulation output · No model is run in this browser.</footer><script src="'+prefix+'assets/gallery.js" defer></script></body></html>\n'


CSS = """
:root{--ink:#182d39;--muted:#586d78;--line:#d8e2e5;--paper:#f5f8f8;--accent:#236c83;--warm:#b74e37}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.6 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}a{color:var(--accent);text-underline-offset:3px}a:hover{color:var(--warm)}header,main,footer{max-width:1160px;margin:auto;padding:28px}header{padding-top:55px}.eyebrow{text-transform:uppercase;letter-spacing:.14em;font-size:.75rem;font-weight:700;color:var(--accent);margin-bottom:12px}h1{font-size:clamp(2rem,4vw,3.4rem);line-height:1.12;letter-spacing:-.035em;margin:0 0 22px;max-width:900px}h2{font-size:1.55rem;line-height:1.3;margin:0 0 16px}h3{font-size:1.15rem;margin:0 0 12px}.intro{font-size:1.12rem;color:var(--muted);max-width:850px}.notice{border-left:4px solid var(--warm);padding:14px 20px;background:#fff8f2;margin:25px 0}.notice p{margin:0}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:24px}.card,section.panel{background:white;border:1px solid var(--line);border-radius:12px;overflow:hidden;box-shadow:0 3px 14px #18303a05}.card img{display:block;width:100%;height:auto}.card-body{padding:24px}.card h2,h1{overflow-wrap:anywhere}.tags{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0}.tag{font-size:.76rem;padding:4px 9px;background:#eaf2f3;border-radius:5px;color:#295261}.tag.off{background:#edf0f2;color:#51606a}.card p{color:var(--muted);font-size:.94rem}.button{display:inline-block;border:1px solid var(--accent);border-radius:6px;padding:8px 14px;font-weight:650;text-decoration:none;background:white;color:var(--accent)}.button:hover{background:#eef5f6}section.panel{padding:26px;margin-bottom:24px}.frame{width:100%;max-width:960px;display:block;margin:auto}.controls{display:flex;align-items:center;gap:12px;margin:16px auto;max-width:900px}.controls input{flex:1;min-width:50px;accent-color:var(--accent)}button{font:inherit;cursor:pointer}.frame-label{font-variant-numeric:tabular-nums;font-weight:600;text-align:center}.fine{font-size:.87rem;color:var(--muted)}.meta{display:grid;grid-template-columns:1fr 1fr;gap:22px}table{border-collapse:collapse;width:100%;font-size:.91rem}th,td{border-bottom:1px solid var(--line);padding:9px 0;text-align:left;vertical-align:top}th{font-weight:500;color:var(--muted);width:46%;padding-right:20px}code{font-size:.78rem;overflow-wrap:anywhere}.downloads{display:flex;flex-wrap:wrap;gap:10px;margin:20px 0}.back{font-size:.88rem;display:inline-block;margin-bottom:28px}footer{font-size:.8rem;color:var(--muted);padding-top:12px;padding-bottom:35px}video{display:block;width:100%;max-width:960px;margin:auto}.comparison{display:grid;grid-template-columns:1fr 1fr;gap:18px}.comparison img{width:100%}.small-list{padding-left:20px;color:var(--muted)}.figure{width:100%;max-width:960px;display:block;margin:auto}@media(max-width:680px){header,main,footer{padding:20px}header{padding-top:32px}.cards,.meta,.comparison{grid-template-columns:1fr}.card-body,section.panel{padding:18px}.controls{gap:7px}.button{padding:6px 9px;font-size:.85rem}th{width:45%}}
"""

JAVASCRIPT = """
document.querySelectorAll('[data-frame-viewer]').forEach(function(viewer){
  const data=JSON.parse(viewer.querySelector('script[type="application/json"]').textContent);
  const slider=viewer.querySelector('input[type="range"]');
  const image=viewer.querySelector('img');
  const label=viewer.querySelector('[data-frame-label]');
  function update(){const i=Number(slider.value), f=data[i];image.src=f.image;image.alt='Relative vorticity and speed at '+f.time_s+' model seconds';label.textContent='Saved frame '+(i+1)+' / '+data.length+' · '+(f.time_s/86400).toFixed(3)+' days · '+f.time_s.toLocaleString()+' s';}
  slider.addEventListener('input',update);
  viewer.querySelector('[data-previous]').addEventListener('click',function(){slider.value=Math.max(0,Number(slider.value)-1);update();});
  viewer.querySelector('[data-next]').addEventListener('click',function(){slider.value=Math.min(data.length-1,Number(slider.value)+1);update();});
  update();
});
"""


def run_summary(run, frames, limits):
    c = run["public_config"]
    d = run["data"]
    invalid = np.flatnonzero(d["jet_coverage"] < 1)
    summary = dict(schema_version=1, name=run["name"], run_id=run["slug"], status="completed",
                start_s=float(d["time_s"][0]), end_s=float(d["time_s"][-1]),
                source_sha256=run["integrity"]["source_sha256"], manifest_sha256=run["integrity"]["manifest_sha256"],
                source_verified=run["integrity"]["source_verified"], manifest_verified=run["integrity"]["manifest_verified"],
                input_config_sha256=run["input_config_sha256"], input_diagnostics_sha256=run["input_diagnostics_sha256"],
                first_incomplete_tracker_time_s=float(d["time_s"][invalid[0]]) if len(invalid) else None,
                mode_numbers=sorted(int(MODE_COLUMN.fullmatch(k).group(1)) for k in d if MODE_COLUMN.fullmatch(k) and "amplitude" in k),
                snapshot_count=len(frames), frames=frames, shared_color_limits=limits,
                field_definition="Relative vorticity = s + psi / LD^2; speed = sqrt(u^2+v^2).",
                orientation="x = r cos(theta), y = r sin(theta); origin lower; positive theta counterclockwise in this north-polar view.",
                display="Fields block-averaged to at most 512 cells per side for rendering only; masked outside science radius.",
                phase_policy="No drift fit or temporal phase unwrapping in this gallery.",
                interpretation="Completion is not a stability certificate. An imposed mode is not evidence of spontaneous selection. Tracker gaps do not measure a physical lifetime.")
    if run.get("validation"):
        summary["purpose"] = "software_validation"
        summary["cpu_gpu_validation"] = run["validation"]
        summary["interpretation"] = run["validation"]["scope"]
    return summary


def write_run_page(run, directory, summary, movie):
    rows = "".join("<tr><th>"+html.escape(k)+"</th><td>"+html.escape(v)+"</td></tr>" for k,v in metadata_rows(run))
    frames = [{"image":f["image"],"time_s":f["time_s"]} for f in summary["frames"]]
    missing = summary["first_incomplete_tracker_time_s"]
    note = ("The original jet tracker first loses complete angular coverage at %.3f days. Later missing mode values remain gaps; this is not a measured physical lifetime." % (missing/86400)) if missing is not None else "The original tracker has full angular coverage, which alone does not establish a single continuous physical jet."
    validation = run.get("validation")
    eyebrow = "Software validation · CPU/GPU comparison passed" if validation else "Completed simulation · exploratory result"
    introduction = ("This %s-second run compares the CPU and GPU implementations. The saved fields document a short software check; they do not test whether a hexagon persists or becomes stable." % format_number(validation["final_time_s"])) if validation else "Polar flow in a one-active-layer quasigeostrophic model. Inspect the flow alongside its measured structure; a completed integration is not a stable hexagon."
    body = '<header><a class="back" href="../../index.html">← All runs</a><div class="eyebrow">'+eyebrow+'</div><h1>'+html.escape(run["name"])+ '</h1><p class="intro">'+introduction+'</p></header><main>'
    if validation:
        error_rows = ''.join('<tr><th>'+field+'</th><td>'+format(validation["initial_errors"][field]["relative_l2"], '.4e')+'</td><td>'+format(validation["final_errors"][field]["relative_l2"], '.4e')+'</td></tr>' for field in FIELD_KEYS)
        body += '<section class="panel"><h2>CPU/GPU agreement</h2><p>Each field passes a relative L2 tolerance of <code>'+format(validation["relative_l2_tolerance"], '.1e')+'</code> at the initial and final saved times.</p><table><thead><tr><th>Field</th><th>Initial relative L2 difference</th><th>Final relative L2 difference</th></tr></thead><tbody>'+error_rows+'</tbody></table><p class="fine">The report matches the saved configuration, source, dependency manifest and validator hashes. CPU/GPU agreement tests these implementations over this short interval; it does not establish spatial convergence or physical stability.</p>'
        if validation["full_match_base"]["applicable"]:
            body += '<p class="fine">An additional unperturbed full-match check passed on both CPU and GPU: the prescribed deep flow cancels the base interface displacement within '+format(validation["full_match_base"]["absolute_thickness_tolerance"], '.1e')+' H. This is a base-state identity, not a stability test.</p>'
        body += '</section>'
    body += '<section class="panel"><h2>Saved flow evolution</h2><div data-frame-viewer><img class="frame" src="frame_000.png" alt="Initial relative vorticity and horizontal speed"><div class="controls"><button class="button" data-previous aria-label="Previous saved frame">Previous</button><input aria-label="Saved simulation frame" type="range" min="0" max="'+str(len(frames)-1)+'" value="0" step="1"><button class="button" data-next aria-label="Next saved frame">Next</button></div><p class="frame-label" data-frame-label>Saved frame 1</p><script type="application/json">'+json.dumps(frames,separators=(',',':'))+'</script></div><p class="fine">Every frame is an actual saved state. Color limits are fixed across every frame and run in this gallery. The maps show the science disk, not a global sphere. No temporal interpolation or optical flow is used.</p></section>'
    if movie:
        body += '<section class="panel"><h2>Discrete frame playback</h2><video controls preload="none" poster="frame_000.png"><source src="movie.mp4" type="video/mp4"></video><p class="fine">Two saved frames per playback second, with equal display duration per frame. This is a sequence of sparse snapshots, not motion resolved between them. Actual model time is printed in each frame.</p></section>'
    initial=run["public_config"]["initial"]
    selection_note = 'Mode six was imposed in this initial condition; its dominance is not spontaneous selection.' if initial["perturbation"] == "mode" and initial["mode"] == 6 else 'Interpret mode growth in the context of the stated initial disturbance and the competing modes.'
    body += '<section class="panel"><h2>Jet measurements and total budgets</h2><img class="figure" src="diagnostics.png" alt="Mode amplitudes, spectrum, tracker coverage, energy and total-PV enstrophy"><div class="notice"><p>'+html.escape(note)+'</p></div><p class="fine">Amplitudes come from the runtime jet-centreline tracker and are provisional without field validation. Weak modes depend on interpolation. '+selection_note+' No late pattern drift is fitted. Energy and enstrophy curves include convection, dissipation and numerical effects; they are not conservation-error estimates.</p></section>'
    body += '<section class="panel"><h2>Configuration and provenance</h2><div class="meta"><table>'+rows+'</table><div><h3>Reproducible inputs</h3><p class="fine">Source SHA-256<br><code>'+summary["source_sha256"]+'</code></p><p class="fine">Dependency manifest SHA-256<br><code>'+summary["manifest_sha256"]+'</code></p><p class="fine">Original configuration SHA-256<br><code>'+summary["input_config_sha256"]+'</code></p><p class="fine">The downloadable configuration retains recognized model fields and omits source comments and arbitrary metadata. Its public-file hash is recorded separately in the summary. Absolute paths, machine names, event logs and raw model files are not included.</p></div></div><div class="downloads"><a class="button" href="config.toml" download>Configuration</a><a class="button" href="diagnostics.csv" download>Diagnostic CSV</a><a class="button" href="summary.json" download>Provenance JSON</a></div></section></main>'
    (directory/"index.html").write_text(page(run["name"], body, "../../"),encoding="utf-8")


def build_gallery(paths, output, title="Saturn polar QG experiments", description="Exploratory saved simulations; completion does not establish a stable hexagon.", mp4=False, ffmpeg="ffmpeg"):
    output = Path(output)
    if output.exists():
        raise FileExistsError("Output exists; choose a new gallery directory")
    # An explicit public title is allowed HTML-like text, escaped on output, but not local paths.
    if any(re.search(r"(?:/Users/|/home/|[A-Za-z]:\\)", value) for value in (title, description)):
        raise ValueError("Public title/description must not contain local filesystem paths")
    runs = [load_run(path) for path in paths]
    if not runs or len({r["slug"] for r in runs}) != len(runs):
        raise ValueError("Provide one or more distinct result directories")
    executable = shutil.which(ffmpeg) if mp4 else None
    if mp4 and executable is None:
        raise ValueError("Requested MP4 generation but ffmpeg is unavailable")
    mode_max = max([float(np.max(a[np.isfinite(a)])) for run in runs for key,a in run["data"].items() if MODE_COLUMN.fullmatch(key) and "amplitude" in key and np.isfinite(a).any()] or [0.0])
    limits = {"relative_vorticity_s_inverse": [-max(r["zmax"] for r in runs), max(r["zmax"] for r in runs)],
              "speed_m_s": [0.0, max(r["vmax"] for r in runs)], "mode_display_floor_m":10.0,"mode_amplitude_max_m":mode_max,
              "scope": "all supplied runs; field limits use all saved science-disk frames, mode limit uses all diagnostics"}
    zlim = max(limits["relative_vorticity_s_inverse"][1], np.finfo(float).tiny)
    vlim = max(limits["speed_m_s"][1], np.finfo(float).tiny)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white","savefig.facecolor":"white"})
    output.mkdir(parents=True, exist_ok=False)
    (output/"assets").mkdir()
    (output/"assets"/"style.css").write_text(CSS)
    (output/"assets"/"gallery.js").write_text(JAVASCRIPT)
    cards, summaries = [], []
    for run in runs:
        print("Rendering " + run["name"], flush=True)
        directory=output/"runs"/run["slug"]
        directory.mkdir(parents=True)
        frames=[]
        for index, original in enumerate(run["frames"]):
            name="frame_%03d.png" % index
            render_frame(run, original, directory/name, zlim, vlim)
            frames.append(dict(original, image=name))
        render_diagnostics(run, directory/"diagnostics.png", mode_max)
        (directory/"config.toml").write_text(toml_text(run["public_config"]))
        with (directory/"diagnostics.csv").open("w",newline="") as stream:
            writer=csv.writer(stream); keys=list(run["data"])
            writer.writerow(keys)
            writer.writerows(zip(*(run["data"][key] for key in keys)))
        if mp4:
            subprocess.run([executable,"-hide_banner","-loglevel","error","-framerate","2","-i",str(directory/"frame_%03d.png"),"-c:v","libx264","-pix_fmt","yuv420p","-vf","pad=ceil(iw/2)*2:ceil(ih/2)*2","-movflags","+faststart",str(directory/"movie.mp4")],check=True)
        summary=run_summary(run,frames,limits)
        summary["public_config_sha256"]=sha256(directory/"config.toml")
        summary["public_diagnostics_sha256"]=sha256(directory/"diagnostics.csv")
        summary["movie"]= "movie.mp4" if mp4 else None
        (directory/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
        write_run_page(run,directory,summary,mp4)
        prefix="runs/"+run["slug"]+"/"
        forced=run["config"]["convection"]["enabled"]
        tag='Convection on' if forced else 'Convection off'
        if run.get("validation"):
            tag = 'CPU/GPU validation passed'
        cards.append('<article class="card"><a href="'+prefix+'index.html"><img src="'+prefix+frames[-1]["image"]+'" alt="Final vorticity and speed for '+html.escape(run["name"],quote=True)+'" loading="lazy"></a><div class="card-body"><h2>'+html.escape(run["name"])+ '</h2><div class="tags"><span class="tag'+('' if forced else ' off')+'">'+tag+'</span><span class="tag">'+str(run["config"]["grid"]["nx"])+ '² grid</span><span class="tag">'+str(len(frames))+' saved frames</span></div><p>LD '+format_number(run["ld"]/1e6,'Mm')+' · '+html.escape(deep_flow_label(run["config"]["physics"]))+' · '+format_number(summary["end_s"]/86400,'model days')+'</p><a class="button" href="'+prefix+'index.html">Explore this run →</a></div></article>')
        summaries.append({"name":run["name"],"run_id":run["slug"],"page":prefix+"index.html","summary":prefix+"summary.json","source_sha256":summary["source_sha256"],"input_config_sha256":summary["input_config_sha256"],"snapshot_count":len(frames),"end_s":summary["end_s"]})
    notice='These are exploratory simulations, not established Saturn hexagons. Read the initial condition and configuration differences before comparing runs.'
    if all(r.get("validation") for r in runs):
        notice='These are short software-validation runs. CPU/GPU agreement does not establish hexagon persistence, spatial convergence, or physical stability. The nearly unchanged flow over this interval is not a stabilization result.'
    if all(r["config"]["initial"]["perturbation"] == "mode" and r["config"]["initial"]["mode"] == 6 for r in runs):
        notice += ' A sixfold initial disturbance tests its evolution; it does not test spontaneous selection of six sides.'
    body='<header><div class="eyebrow">Numerical experiments · North-polar flow</div><h1>'+html.escape(title)+'</h1><p class="intro">'+html.escape(description)+'</p><div class="notice"><p>'+notice+'</p></div></header><main><div class="cards">'+''.join(cards)+'</div><section class="panel" style="margin-top:28px"><h2>How to read the gallery</h2><ul class="small-list"><li>Every image is an actual saved time. Sparse snapshots can show changing structure but do not resolve all motion between frames.</li><li>Vorticity and speed color scales are shared across this gallery; each run uses its configured science disk.</li><li>Mode curves describe the measured jet centreline. Tracker gaps remain missing, and complete coverage alone does not guarantee a continuous physical jet.</li><li>Pattern speeds need a declared rotation frame and a valid coherent signal. This gallery does not fit drift from late rolled-up fields.</li><li>Comparisons should use matched physical times and account for all configuration differences. Completion is not timestep, spatial or physical validation.</li></ul><a href="summary.json" download>Download gallery manifest</a></section></main>'
    (output/"index.html").write_text(page(title,body,""),encoding="utf-8")
    manifest=dict(schema_version=1,title=title,description=description,runs=summaries,shared_color_limits=limits,media_policy="Actual saved frames only; optional MP4 uses equal display durations with no temporal interpolation.")
    (output/"summary.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    (output/"COMPLETED").write_text("Static gallery build completed.\n")
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results",nargs="+",type=Path,help="Completed local result directories; never published directly")
    parser.add_argument("--output",required=True,type=Path,help="New gallery directory; refuses overwrite")
    parser.add_argument("--title",default="Saturn polar QG experiments")
    parser.add_argument("--description",default="Exploratory saved simulations; completion does not establish a stable hexagon.")
    parser.add_argument("--mp4",action="store_true",help="Also encode actual saved frames with ffmpeg, without interpolation")
    parser.add_argument("--ffmpeg",default="ffmpeg",help="ffmpeg executable name or path (not published)")
    args=parser.parse_args()
    try:
        manifest=build_gallery(args.results,args.output,args.title,args.description,args.mp4,args.ffmpeg)
    except (ValueError,FileExistsError,OSError,KeyError,subprocess.CalledProcessError) as error:
        parser.exit(1,"Gallery build failed: "+str(error)+"\n")
    print("Built %d run pages with %d actual saved frames." % (len(manifest["runs"]),sum(r["snapshot_count"] for r in manifest["runs"])))


if __name__ == "__main__":
    main()

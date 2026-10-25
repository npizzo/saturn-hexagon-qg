#!/usr/bin/env python3
"""Provisional CSV growth/phase measurements; never a stability certificate.

Requires NumPy and, on Python <3.11, toml. No snapshot or JLD2 reads.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

SYSTEM_III_PERIOD_S = 38362.4
OBSERVED_2008_2014_PERIOD_S = 38363.01
LIMITS = [
    "CSV fits are provisional and require full-field and numerical validation; no run is classified stable.",
    "A quiet circular jet is not a maintained hexagon. A seeded mode-six pattern is not evidence of mode selection.",
    "Fits use one longest contiguous admissible segment (earliest on ties), not the segment with the largest growth.",
    "Early nonmodal adjustment and changing eigenfunction content can bias fitted growth and phase.",
    "Fits require both the selected-mode amplitude and sqrt(sum(A_m^2)/2) over every configured output mode to remain below the linear ceiling; missing/nonfinite configured amplitudes are rejected. Unresolved modes and full-field nonlinear changes still require validation.",
    "Forced-run slopes include the forcing response and cannot establish an unforced instability growth rate.",
    "The initial mean detected peak speed is a conditional contour diagnostic, not a complete mean-flow constraint.",
    "Complete angular coverage does not ensure the tracer follows one continuous physical jet; inspect fields for ridge switching.",
    "The wrapped-phase increment guard cannot exclude hidden integer-turn aliasing; cadence refinement is still required.",
    "System III conversion assumes configured rotation_rate_s is the frame rate and positive model theta is eastward.",
    "The 38363.01 s target is the 2008-2014 measurement, not a timeless exact target or an acceptance tolerance.",
    "Different source/configuration hashes are recorded for comparison; this tool does not establish that runs differ in only one parameter.",
    "Regime extrema cover saved diagnostic times only and do not certify QG validity; absent columns and undefined barotropic thickness are unavailable, not zero.",
]


def load_toml(path):
    try:
        import tomllib
    except ImportError:
        import toml
        return toml.load(str(path))
    with path.open("rb") as stream:
        return tomllib.load(stream)


def finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_sha256(directory):
    """Same filename/newline concatenation as Julia's source_digest()."""
    files = sorted(directory.glob("*.jl"))
    if not files:
        return None
    return hashlib.sha256(b"\n".join(p.name.encode() + b"\n" + p.read_bytes()
                                     for p in files)).hexdigest()


def linear_fit(t, values):
    x = t - np.mean(t)
    y = values - np.mean(values)
    slope = float(np.dot(x, y) / np.dot(x, x))
    residual = y - slope * x
    total = float(np.dot(y, y))
    r2 = None if total == 0 else 1.0 - float(np.dot(residual, residual)) / total
    return slope, r2


def segments(mask, times, maximum_gap):
    result, current = [], []
    for index, valid in enumerate(mask):
        if current and (not valid or times[index] - times[current[-1]] > maximum_gap):
            result.append(current)
            current = []
        if valid:
            current.append(index)
    if current:
        result.append(current)
    return result


def measure_mode(columns, config, mode, options):
    result = dict(mode=mode, fit_status="unavailable", phase_status="unavailable",
                  fit_start_s=None, fit_end_s=None, sample_count=0,
                  admissible_sample_count=0, qualifying_segments=0,
                  growth_rate_s=None, log_amplitude_r2=None, fitted_log_amplitude_change=None,
                  phase_r2=None, angular_pattern_speed_model_s=None,
                  angular_pattern_speed_system_iii_s=None,
                  model_speed_at_nominal_jet_m_s=None, system_iii_speed_at_nominal_jet_m_s=None,
                  angular_offset_from_observed_s=None, maximum_wrapped_phase_increment_rad=None,
                  initial_amplitude_m=None, final_amplitude_m=None,
                  maximum_resolved_centreline_rms_m=None)
    names = ["m%d_%s" % (mode, suffix) for suffix in ("amplitude_m", "real_m", "imag_m")]
    if not all(name in columns for name in names):
        result["fit_status"] = "missing_mode_columns"
        return result
    t = columns["time_s"]
    amp, real, imag = (columns[name] for name in names)
    result.update(initial_amplitude_m=finite(amp[0]), final_amplitude_m=finite(amp[-1]))
    wind = columns["mean_peak_speed_m_s"]
    if t[0] != 0 or not np.isfinite(wind[0]) or wind[0] <= 0:
        result["fit_status"] = "missing_original_initial_wind"
        return result
    width = float(config["initial"]["jet_width_m"])
    ceiling = options.amplitude_ceiling_fraction * width
    configured_names = ["m%d_amplitude_m" % m for m in range(1, config["output"]["m_max"] + 1)]
    if not all(name in columns for name in configured_names):
        result["fit_status"] = "missing_configured_mode_amplitudes"
        return result
    configured_amplitudes = np.array([columns[name] for name in configured_names])
    coeff_amp = 2 * np.hypot(real, imag)
    with np.errstate(invalid="ignore", over="ignore"):
        resolved_rms = np.sqrt(np.sum(configured_amplitudes**2, axis=0) / 2)
        rms_valid = (np.all(np.isfinite(configured_amplitudes) & (configured_amplitudes >= 0), axis=0)
                     & np.isfinite(resolved_rms) & (resolved_rms <= ceiling))
        tests = {
            "outside_requested_window": (t >= options.start_s) & (t <= options.end_s),
            "amplitude_outside_linear_interval": np.isfinite(amp) & (amp >= options.amplitude_floor_m) & (amp <= ceiling),
            "invalid_complex_coefficient": np.isfinite(real) & np.isfinite(imag) & (coeff_amp > 0),
            "amplitude_coefficient_mismatch": np.isclose(amp, coeff_amp, rtol=1e-6, atol=1e-6),
            "incomplete_angular_coverage": np.isfinite(columns["jet_coverage"]) & (np.abs(columns["jet_coverage"] - 1) <= 1e-12),
            "mean_peak_speed_changed": np.isfinite(wind) & (np.abs(wind / wind[0] - 1) <= options.max_mean_speed_drift_fraction),
            "combined_resolved_amplitude_outside_linear_interval": rms_valid,
        }
    valid = np.logical_and.reduce(list(tests.values()))
    result["exclusion_counts_nonexclusive"] = {key: int(np.count_nonzero(~value)) for key, value in tests.items()}
    result["admissible_sample_count"] = int(np.count_nonzero(valid))
    maximum_gap = options.max_gap_factor * float(config["output"]["diagnostic_interval_s"])
    groups = [g for g in segments(valid, t, maximum_gap) if len(g) >= options.min_samples]
    result["qualifying_segments"] = len(groups)
    if not groups:
        result["fit_status"] = "insufficient_contiguous_samples"
        return result
    selected = max(groups, key=lambda g: (t[g[-1]] - t[g[0]], -t[g[0]]))
    selected_t = t[selected]
    growth, growth_r2 = linear_fit(selected_t, np.log(amp[selected]))
    result.update(fit_status="estimated_provisional", fit_start_s=float(selected_t[0]),
                  fit_end_s=float(selected_t[-1]), sample_count=len(selected),
                  growth_rate_s=growth, log_amplitude_r2=growth_r2,
                  fitted_log_amplitude_change=growth * (selected_t[-1] - selected_t[0]),
                  maximum_resolved_centreline_rms_m=float(np.max(resolved_rms[selected])))
    angles = np.arctan2(imag[selected], real[selected])
    increments = np.angle(np.exp(1j * np.diff(angles)))
    largest = float(np.max(np.abs(increments)))
    result["maximum_wrapped_phase_increment_rad"] = largest
    if largest >= options.max_phase_step_rad:
        result["phase_status"] = "ambiguous_phase_increment"
        return result
    slope, phase_r2 = linear_fit(selected_t, np.unwrap(angles))
    angular = -slope / mode
    frame = float(config["planet"]["rotation_rate_s"])
    system_iii = angular + frame - 2 * math.pi / SYSTEM_III_PERIOD_S
    target = 2 * math.pi / OBSERVED_2008_2014_PERIOD_S - frame
    radius = float(config["initial"]["jet_radius_m"])
    result.update(phase_status="estimated_provisional_aliasing_not_excluded", phase_r2=phase_r2,
                  angular_pattern_speed_model_s=angular,
                  angular_pattern_speed_system_iii_s=system_iii,
                  model_speed_at_nominal_jet_m_s=radius * angular,
                  system_iii_speed_at_nominal_jet_m_s=radius * system_iii,
                  angular_offset_from_observed_s=angular - target)
    return result


def regime_summary(columns, config):
    """Optional sampled extrema, preserving missing/undefined/nonfinite states."""
    ld = config.get("physics", {}).get("deformation_radius_m")
    barotropic = isinstance(ld, (int, float)) and math.isinf(ld) and ld > 0
    result = {"scope": "all_saved_diagnostic_rows_not_continuous_time_extrema", "columns": {}}
    metrics = [("max_abs_hprime_over_H", "maximum_abs_hprime_over_H", np.max),
               ("min_h_over_H", "minimum_h_over_H", np.min),
               ("max_abs_zeta_over_f", "maximum_abs_zeta_over_f", np.max)]
    for source, target, operation in metrics:
        result[target] = None
        if source not in columns:
            result["columns"][source] = {"status": "absent_legacy_or_unrecorded_column"}
            continue
        values = np.asarray(columns[source], dtype=float)
        good = np.isfinite(values)
        finite_count = int(np.count_nonzero(good))
        info = {"finite_samples": finite_count, "nonfinite_samples": int(values.size - finite_count),
                "total_samples": int(values.size)}
        if barotropic and source != "max_abs_zeta_over_f":
            info["status"] = "unexpected_finite_barotropic_depth_values" if finite_count else "not_applicable_barotropic"
        elif finite_count:
            result[target] = float(operation(values[good]))
            info["status"] = "all_samples_finite" if finite_count == values.size else "partly_nonfinite"
        else:
            info["status"] = "no_finite_samples"
        result["columns"][source] = info
    return result


def summarize_run(directory, options):
    directory = Path(directory).resolve()
    run = dict(result_directory=str(directory), execution_status="incomplete",
               input_status="readable", issues=[], hashes={}, modes=[])
    config, metadata = {}, {}
    for filename, target in [("config.toml", config), ("metadata.toml", metadata)]:
        path = directory / filename
        try:
            target.update(load_toml(path))
            run["hashes"][filename + "_sha256"] = sha256(path)
        except Exception as error:
            run["issues"].append("Cannot read %s: %s" % (filename, error))
            run["input_status"] = "invalid_or_missing_input"
    run["name"] = config.get("experiment", {}).get("name", directory.name)
    run["forcing_enabled"] = config.get("convection", {}).get("enabled")
    run["initial_perturbation"] = config.get("initial", {}).get("perturbation")
    run["config"] = clean(config)
    ld = config.get("physics", {}).get("deformation_radius_m")
    if isinstance(ld, (int, float)) and math.isinf(ld):
        run["config"]["physics"]["deformation_radius_m"] = "inf"
    run["declared_metadata"] = clean(metadata)
    actual_source = source_sha256(directory / "src")
    actual_manifest = sha256(directory / "Manifest.toml") if (directory / "Manifest.toml").is_file() else None
    run["hashes"].update(snapshot_source_sha256=actual_source, snapshot_manifest_sha256=actual_manifest)
    run["integrity_checks"] = {}
    for label, actual, declared in [("source", actual_source, metadata.get("source_sha256")),
                                     ("manifest", actual_manifest, metadata.get("manifest_sha256"))]:
        status = "unverified" if actual is None or not declared else "matches" if actual == declared else "mismatch"
        run["integrity_checks"][label] = status
        if status != "matches":
            run["issues"].append("%s digest %s" % (label, status))
    completion = None
    marker = directory / "COMPLETED"
    if marker.is_file():
        completion = finite(marker.read_text().strip())
        if completion is None:
            run["issues"].append("COMPLETED marker has no finite final time")
    run["completed_marker_time_s"] = completion
    failed = (directory / "FAILED.txt").exists()
    if failed:
        run["execution_status"] = "failed"
        run["issues"].append("FAILED.txt is present")
    try:
        path = directory / "diagnostics.csv"
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            headers = reader.fieldnames
            if not headers or len(headers) != len(set(headers)):
                raise ValueError("missing or duplicate diagnostic headers")
            rows = list(reader)
        if not rows:
            raise ValueError("empty diagnostics")
        columns = {key: np.array([finite(row.get(key)) if finite(row.get(key)) is not None else np.nan
                                  for row in rows], dtype=float) for key in headers}
        for key in ("time_s", "jet_coverage", "mean_peak_speed_m_s"):
            if key not in columns:
                raise ValueError("missing diagnostics column " + key)
        t = columns["time_s"]
        if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0) or t[0] < 0:
            raise ValueError("diagnostic times must be finite, nonnegative and strictly increasing; no sorting or deduplication is performed")
        run["hashes"]["diagnostics.csv_sha256"] = sha256(path)
        run["regime_summary"] = regime_summary(columns, config)
        run.update(diagnostic_rows=len(rows), initial_time_s=float(t[0]), final_time_s=float(t[-1]))
        duration = float(config["numerics"]["duration_s"])
        tolerance = max(1e-6, duration * 1e-10)
        if not failed and completion is not None and abs(completion - duration) <= tolerance and abs(t[-1] - completion) <= tolerance:
            run["execution_status"] = "completed"
        elif not failed:
            run["issues"].append("No consistent completion marker and final diagnostic at configured duration")
        for value in (config["initial"]["jet_width_m"], config["initial"]["jet_radius_m"],
                      config["planet"]["rotation_rate_s"], config["output"]["diagnostic_interval_s"], duration):
            if finite(value) is None or float(value) <= 0:
                raise ValueError("required physical scales and time intervals must be finite and positive")
        run["initial_mean_peak_speed_m_s"] = finite(columns["mean_peak_speed_m_s"][0])
        run["minimum_jet_coverage"] = finite(np.min(columns["jet_coverage"]))
        run["amplitude_ceiling_m"] = options.amplitude_ceiling_fraction * float(config["initial"]["jet_width_m"])
        frame = float(config["planet"]["rotation_rate_s"])
        run["conditional_observed_angular_speed_model_s"] = 2 * math.pi / OBSERVED_2008_2014_PERIOD_S - frame
        run["modes"] = [measure_mode(columns, config, mode, options) for mode in range(1, 33)]
        run["mode_six"] = run["modes"][5]
    except Exception as error:
        run["input_status"] = "invalid_or_missing_input"
        run["issues"].append("Cannot estimate diagnostics: " + str(error))
    if run["execution_status"] != "completed":
        run["issues"].append("Partial-run fits, if any, must not be presented as a completed experiment")
    return clean(run)


def parser():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("results", nargs="+", type=Path, help="Result directories containing config, metadata, diagnostics and completion marker")
    cli.add_argument("--output", required=True, type=Path, help="NEW directory; existing paths are refused")
    cli.add_argument("--start-s", type=float, default=20000.0)
    cli.add_argument("--end-s", type=float, default=200000.0)
    cli.add_argument("--amplitude-floor-m", type=float, default=500.0)
    cli.add_argument("--amplitude-ceiling-fraction", type=float, default=0.05, help="Fraction of configured jet_width_m")
    cli.add_argument("--min-samples", type=int, default=10)
    cli.add_argument("--max-mean-speed-drift-fraction", type=float, default=0.05)
    cli.add_argument("--max-gap-factor", type=float, default=1.5, help="Break segments across gaps longer than this times configured diagnostic cadence")
    cli.add_argument("--max-phase-step-rad", type=float, default=math.pi / 2, help="Reject phase fits at or above this wrapped increment; cannot detect all aliases")
    return cli


def main(argv=None):
    cli = parser()
    options = cli.parse_args(argv)
    numeric = [options.start_s, options.end_s, options.amplitude_floor_m,
               options.amplitude_ceiling_fraction, options.max_mean_speed_drift_fraction,
               options.max_gap_factor, options.max_phase_step_rad]
    if not all(math.isfinite(x) for x in numeric) or not (0 <= options.start_s < options.end_s):
        cli.error("finite values and 0 <= start-s < end-s are required")
    if min(numeric[2:]) <= 0 or options.min_samples < 3 or options.max_gap_factor < 1 or not (0 < options.max_phase_step_rad < math.pi):
        cli.error("positive thresholds, min-samples >= 3, max-gap-factor >= 1 and 0 < max-phase-step-rad < pi are required")
    output = options.output.resolve()
    if output.exists() or output.is_symlink():
        cli.error("output already exists; choose a NEW directory")
    runs = [summarize_run(path, options) for path in options.results]
    settings = {key: value for key, value in vars(options).items() if key not in ("results", "output")}
    summary = dict(schema_version=1, interpretation="provisional_csv_measurements_not_stability_classification",
                   analysis_script_sha256=sha256(Path(__file__)), settings=settings,
                   system_iii_period_s=SYSTEM_III_PERIOD_S,
                   observed_2008_2014_period_s=OBSERVED_2008_2014_PERIOD_S,
                   observation_source="https://doi.org/10.1002/2013GL059078", limitations=LIMITS, runs=runs)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "summary.json").open("x") as stream:
        json.dump(clean(summary), stream, indent=2, allow_nan=False)
        stream.write("\n")
    csv_rows = []
    for run in runs:
        for mode in run["modes"]:
            row = {key: value for key, value in mode.items() if key != "exclusion_counts_nonexclusive"}
            row.update(result_directory=run["result_directory"], name=run["name"],
                       execution_status=run["execution_status"], input_status=run["input_status"])
            csv_rows.append(row)
    fields = list(csv_rows[0]) if csv_rows else ["result_directory", "mode", "fit_status"]
    with (output / "modes.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(csv_rows)
    lines = ["Provisional experiment summary", "Requested window: %.9g to %.9g s" % (options.start_s, options.end_s),
             "Blank CSV cells / JSON null mean unavailable, not zero.", ""]
    for run in runs:
        six = run.get("mode_six", {})
        lines.extend(["%s: %s; %s" % (run["name"], run["execution_status"], run["input_status"]),
                      "  " + run["result_directory"],
                      "  m6 growth=%s s^-1; c_model=%s rad/s; c_SystemIII=%s rad/s; status=%s / %s" %
                      (six.get("growth_rate_s"), six.get("angular_pattern_speed_model_s"),
                       six.get("angular_pattern_speed_system_iii_s"), six.get("fit_status"), six.get("phase_status"))])
        regime = run.get("regime_summary")
        if regime:
            lines.append("  Sampled regime extrema: max |h'/H|=%s; min h/H=%s; max |zeta/f|=%s" %
                         (regime["maximum_abs_hprime_over_H"], regime["minimum_h_over_H"],
                          regime["maximum_abs_zeta_over_f"]))
            lines.append("  Regime columns: " + "; ".join(key + "=" + info["status"]
                                                          for key, info in regime["columns"].items()))
        lines.extend("  Issue: " + issue for issue in run["issues"])
    lines.extend(["", "Interpretation limits:"] + ["- " + limit for limit in LIMITS])
    with (output / "report.txt").open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    print("Wrote provisional measurements for %d result directories to %s" % (len(runs), output))
    return 0


if __name__ == "__main__":
    sys.exit(main())

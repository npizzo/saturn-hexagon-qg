#!/usr/bin/env python3
"""Screen exactly one hashed stage of circular QG configurations on a CPU.

Usage: python scripts/stabilization_linear.py STAGE_DIRECTORY --output NEW_DIRECTORY
Dependencies: numpy, scipy, and toml for Python <3.11; install before the job.
No nonlinear integrations, installation, parameter sweeps, or submissions occur.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import radial_stability as radial

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def hashed_files(manifest, base, basenames=False):
    """Strict portable SHA256<TAB>relative-path manifest, with no escaping links."""
    result, seen = [], set()
    for line in manifest.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError("Malformed hash manifest: " + str(manifest))
        expected, name = parts
        relative = Path(name)
        if (len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected)
                or relative.is_absolute() or ".." in relative.parts or not name
                or name in seen or (basenames and relative.name != name)):
            raise ValueError("Invalid or duplicate manifest entry: " + name)
        path = (base/relative).resolve()
        if base.resolve() not in path.parents or not path.is_file():
            raise ValueError("Manifest file is missing or escapes its root: " + name)
        data = path.read_bytes()
        if sha(data) != expected:
            raise ValueError("Hash mismatch: " + str(path))
        seen.add(name); result.append((name, expected, data))
    if not result:
        raise ValueError("Empty hash manifest: " + str(manifest))
    return result


def load_stage(stage, project=ROOT):
    stage, project = Path(stage).resolve(), Path(project).resolve()
    provenance = stage.parent/"provenance.tsv"
    sources = hashed_files(provenance, project)
    names = {name for name, _, _ in sources}
    required = {"Manifest.toml", "scripts/radial_stability.py", "scripts/stabilization_linear.py"}
    required.update("src/"+p.name for p in (project/"src").glob("*.jl"))
    if not any(name.startswith("src/") for name in required) or not required <= names:
        raise ValueError("Provenance must cover Manifest, both analysis scripts, and all Julia source files")
    for local in (Path(__file__), Path(radial.__file__)):
        if sha(local.read_bytes()) != sha((project/"scripts"/local.name).read_bytes()):
            raise ValueError("Executing analysis script differs from the provenance checkout")
    members = hashed_files(stage/"members.tsv", stage, basenames=True)
    configs = []
    for name, expected, data in members:
        if not name.endswith(".toml"):
            raise ValueError("Stage members must be TOML configurations")
        try:
            import tomllib as toml
        except ImportError:
            import toml
        config = toml.loads(data.decode("utf-8"))
        configs.append((name, expected, data, config))
    return configs, {name: expected for name, expected, _ in sources}


def validate_config(c, outer):
    for group, keys in (("initial", ("jet_radius_m", "jet_width_m", "cyclone_radius_m")),
                        ("physics", ("deep_jet_radius_m", "deep_jet_width_m")),
                        ("planet", ("radius_m", "rotation_rate_s")),
                        ("grid", ("nx", "length_m"))):
        if any(not np.isfinite(c[group][k]) or c[group][k] <= 0 for k in keys):
            raise ValueError("Invalid positive physical parameters in " + group)
    if any(not np.isfinite(c["initial"][k]) for k in ("jet_speed_m_s", "cyclone_speed_m_s")):
        raise ValueError("Upper speeds must be finite")
    if not np.isfinite(outer) or not 0 < outer <= c["boundary"]["inner_radius_m"]:
        raise ValueError("Outer radius must be inside or on the untapered cap")
    if outer <= c["initial"]["jet_radius_m"] + 2*c["initial"]["jet_width_m"]:
        raise ValueError("Outer boundary must enclose the jet and its two-width annulus")
    # Uses configured LD and deep flow exactly; no ratio or fraction override.
    radial.circular_base(c, np.array([0.0, outer/2]), c["physics"]["deformation_radius_m"])


def resolution_flags(c, outer, nodes, first, last, growth, tolerance, grids):
    dr = outer/(nodes+1)
    scales = {"jet_width": c["initial"]["jet_width_m"]}
    if c["initial"]["cyclone_speed_m_s"] != 0:
        scales["cyclone_radius"] = c["initial"]["cyclone_radius_m"]
    ld = c["physics"]["deformation_radius_m"]
    if math.isfinite(ld):
        scales["deformation_radius"] = ld
        if c["physics"].get("deep_profile", "annular") == "annular" and c["physics"]["deep_jet_speed_m_s"] != 0:
            scales["deep_width"] = c["physics"]["deep_jet_width_m"]
    flags = ["unexamined_modes_above_%d" % last, "outer_boundary_sensitivity_not_tested"]
    if first > 1:
        flags.append("omitted_low_azimuthal_modes")
    if len(grids) == 1:
        flags.append("single_radial_grid_no_convergence_check")
    if any(scale/dr < 8 for scale in scales.values()):
        flags.append("weak_radial_resolution_fewer_than_8_cells_on_physical_scale")
    if growth[-1] > tolerance and growth[-1] >= max(growth):
        flags.append("fastest_growth_at_high_m_edge_extend_search")
    return flags, {key: float(scale/dr) for key, scale in scales.items()}


def screen_config(c, nodes, first, last, outer, tolerance):
    """All spectra retained; no neutral eigenvalue is presented as a pattern speed."""
    validate_config(c, outer)
    ld = c["physics"]["deformation_radius_m"]
    rows, spectra, eigenfunctions = [], [], {}
    for n in nodes:
        growth = []
        for m in range(first, last+1):
            r, operator = radial.radial_operator(n, outer, m, ld)
            base = radial.circular_base(c, r, ld)
            result = radial.solve_modes(r, operator, base["omega"], base["q_r_over_r"], outer, m)
            if result["relative_residual"] > 1e-8:
                raise ArithmeticError("Leading eigenpair residual exceeds 1e-8")
            rate = float(result["growth_rate_s"]); growth.append(rate)
            growing = rate > tolerance
            # Ties between neutral modes carry no unique modal drift information.
            angular = float(result["angular_pattern_speed_s"]) if growing else None
            rows.append({"radial_nodes": n, "m": m, "growth_rate_s": rate,
                         "growing_angular_pattern_speed_s": angular,
                         "growing_frequency_real_s": m*angular if growing else None,
                         "relative_residual": float(result["relative_residual"]),
                         "mode_status": "growth_above_tolerance" if growing else "no_growth_above_tolerance_in_this_discretization",
                         "flags": ""})
            for index, value in enumerate(result["eigenvalues"]):
                spectra.append({"radial_nodes": n, "m": m, "eigenvalue_index": index,
                                "c_real_s": float(value.real), "c_imag_s": float(value.imag),
                                "growth_rate_s": float(m*value.imag)})
            if growing:
                eigenfunctions["N%d_m%d_phi" % (n, m)] = result["leading_phi"]
            eigenfunctions["N%d_r_m" % n] = r
        flags, cells = resolution_flags(c, outer, n, first, last, growth, tolerance, nodes)
        for row in rows[-(last-first+1):]:
            row["flags"] = ";".join(flags)
            row["minimum_cells_per_physical_scale"] = min(cells.values())
    return rows, spectra, eigenfunctions


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def run(stage, output, nodes=(128,), modes=(1, 32), outer_radius=None, tolerance=1e-12, project=ROOT):
    stage, output = Path(stage).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("Output must be a NEW directory")
    first, last = modes
    if first < 1 or last < first or not nodes or len(set(nodes)) != len(nodes):
        raise ValueError("Invalid mode interval or duplicate/empty radial grids")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Growth tolerance must be finite and nonnegative")
    configs, provenance = load_stage(stage, project)
    for _, _, _, c in configs:
        outer = c["boundary"]["inner_radius_m"] if outer_radius is None else outer_radius
        validate_config(c, outer)
        for n in nodes:
            radial.radial_operator(n, outer, first, c["physics"]["deformation_radius_m"])
    output.mkdir(parents=True, exist_ok=False)
    summaries = []
    for name, digest, data, c in configs:
        outer = c["boundary"]["inner_radius_m"] if outer_radius is None else outer_radius
        rows, spectra, functions = screen_config(c, nodes, first, last, outer, tolerance)
        dest = output/Path(name).stem; dest.mkdir()
        (dest/"config.toml").write_bytes(data)
        write_csv(dest/"modes.csv", rows); write_csv(dest/"eigenvalues.csv", spectra)
        np.savez_compressed(str(dest/"growing_eigenfunctions.npz"), **functions)
        meta = {"status": "inviscid_circular_base_screen_not_a_stability_certificate",
                "config": name, "config_sha256": digest, "provenance_sha256": provenance,
                "deep_profile": c["physics"].get("deep_profile", "annular"),
                "deep_match_fraction": c["physics"].get("deep_match_fraction", 1.0) if c["physics"].get("deep_profile", "annular") == "matched_upper" else None,
                "deformation_radius_m": c["physics"]["deformation_radius_m"],
                "outer_radius_m": outer, "radial_nodes": list(nodes), "mode_range": [first, last],
                "growth_tolerance_s": tolerance, "python": platform.python_version(),
                "numpy": np.__version__, "scipy": scipy.__version__,
                "method": "Second-order radial flux differences, phi(0)=phi(R)=0; exp[i*m*(theta-c*t)], growth=m*Im(c), angular pattern speed=Re(c). SI units and configured coordinate frame; no observational-frame fit.",
                "assumptions": ["Untapered circular upper ring plus cyclone; configured nonaxisymmetric seed omitted.",
                    "Matched deep flow is a fraction of the entire circular upper flow; it contains no polygonal seed.",
                    "No convection, sponge, or hyperviscosity; this is not the nonlinear periodic model.",
                    "Artificial circular outer Dirichlet boundary and second-order radial differences require independent sensitivity tests.",
                    "No growing eigenvalue in a finite mode band and grid is not a stability certificate. Eight cells per scale is a warning heuristic, not convergence.",
                    "No leading drift is selected from a neutral spectrum. All raw eigenvalues remain in eigenvalues.csv.",
                    "Growing eigenfunctions have arbitrary complex phase and unit maximum magnitude; no eigenbranch tracking is performed."]}
        (dest/"metadata.json").write_text(json.dumps(radial.clean_json(meta), indent=2, allow_nan=False)+"\n")
        for n in nodes:
            group = [row for row in rows if row["radial_nodes"] == n]
            best = max(group, key=lambda row: row["growth_rate_s"])
            summaries.append(dict(config=name, config_sha256=digest,
                dominant_growing_m=best["m"] if best["growth_rate_s"] > tolerance else None,
                **{key: value for key, value in best.items() if key != "m"}))
        (dest/"COMPLETED").write_text("%d eigenproblems; no stability certificate\n" % len(rows))
        print("Screened %s: %d eigenproblems" % (name, len(rows)), flush=True)
    write_csv(output/"stage_summary.csv", summaries)
    (output/"COMPLETED").write_text("%d explicit stage configurations\n" % len(configs))
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", type=Path); parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--m-range", default="1:32")
    parser.add_argument("--radial-nodes", default="128")
    parser.add_argument("--outer-radius", type=float, help="default: each config's untapered-cap edge")
    parser.add_argument("--growth-tolerance", type=float, default=1e-12)
    args = parser.parse_args(argv)
    modes = tuple(int(v) for v in args.m_range.split(":"))
    if len(modes) != 2:
        parser.error("--m-range must be FIRST:LAST")
    return run(args.stage, args.output, tuple(radial.number_list(args.radial_nodes, integer=True)),
               modes, args.outer_radius, args.growth_tolerance)


if __name__ == "__main__":
    main()

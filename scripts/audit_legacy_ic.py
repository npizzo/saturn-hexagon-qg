#!/usr/bin/env python3
"""Read a legacy MATLAB initializer as text and audit its q0 formula.

No MATLAB evaluation, solver execution, directory change, or model initialization
is performed. Python >= 3.8, standard library only.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

from audit_ptolemaic_map import Parameters, analytic_q


def first_number(active, name):
    match = re.search(r"(?<![A-Za-z_0-9])" + re.escape(name)
                      + r"\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*;", active)
    if not match:
        raise ValueError("Cannot identify the first literal assignment to " + name)
    return float(match.group(1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("attachment", type=Path)
    parser.add_argument("--output", type=Path, help="JSON output path; default: stdout")
    args = parser.parse_args()
    raw = args.attachment.read_bytes()
    source = raw.decode("utf-8")
    initializer = source.split("%% now look at the data", 1)[0]
    # This is text normalization for a narrowly specified source audit, not a
    # MATLAB parser or evaluator. Unknown q0 formulas are rejected explicitly.
    statements = "\n".join(line.split("%", 1)[0] for line in initializer.splitlines())
    active = re.sub(r"\s+", "", statements.replace("...", ""))
    expected = ("0+gamma*(-exp(2*b(j))*n+exp(2*n*b(j))*n^3*eps^2+"
                "a0^2*(exp(4*b(j))*(2+n)-exp(4*n*b(j))*n^3*(1+2*n)*eps^4))/"
                "(-exp(2*b(j))*n+exp(2*n*b(j))*n^3*eps^2)")
    formula = re.search(r"q0\(i,j\)=(.*?);", active)
    if formula is None or formula.group(1) != expected or "gamma=2*omega;" not in active:
        raise ValueError("Unrecognized source formula; inspect it before extending this audit")
    n_value = first_number(statements, "n")
    if n_value != int(n_value) or n_value == 0:
        raise ValueError("This audit requires a nonzero integer n")
    n, epsilon = int(n_value), first_number(statements, "eps")
    amplitude, planetary_omega = first_number(statements, "a0"), first_number(statements, "omega")
    gamma = 2 * planetary_omega
    p = Parameters(n=n, epsilon=epsilon, amplitude=amplitude, k=1.0,
                   gamma=gamma, f0=gamma, omega=0.0)

    def attachment_formula(b):
        # Direct translation of the recognized MATLAB numerator/denominator,
        # independently compared below with the map-derived analytic curl.
        return gamma * (-math.exp(2*b)*n + math.exp(2*n*b)*n**3*epsilon**2
                        + amplitude**2*(math.exp(4*b)*(2+n)
                                        - math.exp(4*n*b)*n**3*(1+2*n)*epsilon**4)) / (
                            -math.exp(2*b)*n + math.exp(2*n*b)*n**3*epsilon**2)

    points = []
    for b in (-2.0, -1.0, -0.1, 0.0):
        for a in (0.0, 0.25, 1.0, 3.0):
            actual = attachment_formula(b)
            comparison = analytic_q(a, b, p, -1)
            points.append({"a": a, "b": b, "attachment_q0": actual,
                           "map_derived_total_pv": comparison,
                           "absolute_difference": abs(actual-comparison)})
    max_error = max(x["absolute_difference"] for x in points)
    required_pipeline = ["q2=griddata(x,y,q0,Xi,Eta,'nearest');",
                         "q2=q2-min(min(q2));", "q2(isnan(q2))=0;",
                         "q2(isinf(q2))=0;", "Q=q2;", "saveq.txtQ-ascii"]
    checks = {"recognized_attachment_formula": True,
              "recognized_active_saved_field_pipeline": all(x in active for x in required_pipeline),
              "recognized_sampling_and_grid": all(x in active for x in
                  ("b=-2:0.0025:0;", "a=0:.01:1.1*n*2*pi;", "ng=1*100;")),
              "c_is_assigned_without_an_active_reference": len(re.findall(r"\bc\b", statements)) == 1,
              "axisymmetric_subtraction_is_inactive": "q1(i,j)=" not in active and "q0=q0-q1;" not in active,
              "map_is_locally_admissible_at_b_zero": abs(n*epsilon) < 1,
              "sixteen_point_total_pv_formula_match": len(points) == 16 and max_error < 1e-12}
    radius_values = [amplitude*(1-abs(epsilon)), amplitude, amplitude*(1+abs(epsilon))]
    result = {
        "schema_version": 1,
        "status": "formula_verified_saved_field_solver_semantics_unresolved",
        "attachment": {"path": str(args.attachment.resolve()),
                       "sha256": hashlib.sha256(raw).hexdigest()},
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "map_audit_script_sha256": hashlib.sha256(Path(__file__).with_name("audit_ptolemaic_map.py").read_bytes()).hexdigest(),
        "python_version": sys.version.split()[0],
        "execution": {"matlab_executed": False, "legacy_solver_executed": False,
                      "model_initialized_or_integrated": False,
                      "production_files_modified": False,
                      "source_read_as_text_only": True},
        "active_parameters": {"n": n, "polygon_mode_abs_n_plus_one": abs(n+1),
                              "epsilon_used_to_construct_map_and_q0": epsilon,
                              "a0": amplitude, "planetary_omega": planetary_omega,
                              "gamma": gamma, "c_assigned_but_unused": first_number(statements, "c"),
                              "b_samples": "-2:0.0025:0", "a_samples": "0:0.01:1.1*n*2*pi",
                              "ng": 100, "grid_shape_expected": [201, 201],
                              "later_epsilon_reset": "eps=0 occurs after geometry/q0 computation; q1 stays zero"},
        "formula": {"definitions": "p=exp(b), h=epsilon*exp(n*b), D=p^2-n^2*h^2",
                    "attachment_q0": "gamma*[1-A0^2*((n+2)*p^4-n^3*(1+2*n)*h^4)/(n*D)]",
                    "matched_conventions": "u=-psi_y, v=psi_x; qBT=laplacian(psi)+f; f=gamma-gamma*r^2",
                    "matched_frequency": "W=-gamma*A0^2*(p^2-n*h^2)/(2*n)",
                    "matched_free_wave_frequency_constant": 0.0,
                    "interpretation": "Unshifted q0 matches total planar barotropic PV under these explicit conventions. This does not identify the legacy solver's q.txt variable."},
        "verification": {"number_of_points": len(points), "maximum_absolute_error": max_error,
                         "points": points, "boundary_q0_b_zero": attachment_formula(0.0),
                         "innermost_sample_q0_b_minus_two": attachment_formula(-2.0)},
        "saved_field_pipeline": ["Nearest-neighbor scattered interpolation of q0 to loaded Xi,Eta",
                                 "Subtract gridded finite minimum", "Replace NaN/Inf with zero",
                                 "Save Q=q2 to q.txt; computed spherical f and q3 mask are unused"],
        "interpolation_caution": "MATLAB griddata nearest returns NaN outside the sample convex hull; SciPy nearest extrapolation alone does not reproduce the zero-filled hull exterior.",
        "spherical_mapping": {
            "northern_stereographic_z": "z=(1-rho^2)/(1+rho^2)",
            "computed_but_unused_f": "2*planetary_omega*(1-rho^2)/(1+rho^2)",
            "spherical_f_quadratic_coefficient": 4*planetary_omega,
            "matched_planar_f_quadratic_coefficient": gamma,
            "coordinate_scaling_status": "Direct planar x=Xi,y=Eta identification has different quadratic Coriolis coefficients; chart/length/time conversion remains to be established.",
            "radius_and_implied_north_latitude_degrees": [
                {"stereographic_radius": r, "latitude_degrees": 90-2*math.degrees(math.atan(r))}
                for r in radius_values],
            "latitude_interpretation": "Geometry implied by the script's plotting map; not a calibrated Saturn hexagon position."},
        "unresolved": {
            "solver_folder_named_in_attachment": "Rick.Sphere.Code/vortex.waves",
            "solver_reader_inspected": False,
            "prior_targeted_read_result": "macOS Operation not permitted; this script does not retry access",
            "q_txt_semantics": "Unknown: relative vorticity, total/absolute vorticity, or a solver-specific anomaly cannot be assigned from filenames/comments alone.",
            "remaining_checks": ["Reader/inversion semantics and hemisphere signs", "Handling of global mean and minimum subtraction",
                                 "Coordinate/time units, map orientation and output ordering", "Complete outer patch, velocity and initial equilibrium residual"]},
        "campaign_scope": "This legacy IC verification does not block the existing annular-jet campaign. It does not certify a complete exact, stable, or correctly drifting Ptolemaic QG initial state.",
        "source_locations": {"parameters_and_map": "attachment lines 16-31", "q0": "lines 33-44",
                             "inactive_axisymmetric_subtraction": "lines 45-56", "interpolation_and_shift": "lines 68-78",
                             "unused_f_and_active_Q": "lines 81-95", "save": "line 146",
                             "stereographic_plotting_map": "lines 191-196, 270-272"},
        "documentation_sources": ["https://www.mathworks.com/help/matlab/ref/griddata.html",
                                  "https://www.mathworks.com/help/matlab/ref/double.min.html"],
        "checks": checks, "all_checks_passed": all(checks.values()),
    }
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print("Saved %s; %d/%d checks passed; maximum formula error %.3g" % (
            args.output, sum(checks.values()), len(checks), max_error))
    else:
        sys.stdout.write(text)
    if not result["all_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

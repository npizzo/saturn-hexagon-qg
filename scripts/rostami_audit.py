#!/usr/bin/env python3
"""Audit a regularized Rostami2017 circular QG base without integrating it.

Usage: python scripts/rostami_audit.py CONFIG --output NEW_DIRECTORY
The radial pressure calculations use psi(R)=psi_d(R)=0 at the audit's outer
radius. They do NOT use the evolution model's full-square mean gauge and are
not its diagnosed layer thickness. Numpy/scipy and toml (Python<3.11) suffice.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import radial_stability as radial


def audit(config, samples=16385, outer_radius=None):
    if config["initial"].get("profile") != "rostami2017":
        raise ValueError("This audit requires initial.profile=rostami2017")
    if not isinstance(samples, int) or isinstance(samples, bool) or samples < 257:
        raise ValueError("samples must be an integer >=257")
    outer = config["boundary"]["inner_radius_m"] if outer_radius is None else outer_radius
    if not np.isfinite(outer) or not 0 < outer <= config["boundary"]["inner_radius_m"]:
        raise ValueError("audit radius must lie inside the untapered science cap")
    r = np.linspace(0, outer, samples)
    ld = config["physics"]["deformation_radius_m"]
    base = radial.circular_base(config,r,ld)
    p = radial.rostami_parameters(config["initial"])
    f0 = 2*config["planet"]["rotation_rate_s"]
    f = f0-base["gamma_code"]*r*r/2
    if np.min(f) <= 0:
        raise ValueError("local gamma-plane Coriolis parameter must remain positive")
    zeta = base["u_prime"]+base["omega"]
    curvature_ratio = np.abs(base["omega"])/f
    relative_vorticity = np.abs(zeta)/f
    def outer_integral(values):
        cumulative = np.r_[0., np.cumsum((values[:-1]+values[1:])*np.diff(r)/2)]
        return cumulative-cumulative[-1]
    # Unit [m^2/s] streamfunction relative to the circle r=outer. No mean shift.
    psi_shear = outer_integral(base["u"]-base["u_deep"])
    if np.isfinite(ld):
        radial_chi = psi_shear/(f0*ld*ld)
        qg_pressure = outer_integral(f0*base["u"])/(f0*ld)**2
        localf_pressure = outer_integral(f*base["u"])/(f0*ld)**2
        centrifugal_pressure = outer_integral(base["u"]*base["omega"])/(f0*ld)**2
    else:
        radial_chi = qg_pressure = localf_pressure = centrifugal_pressure = None
    dx = config["grid"]["length_m"]/config["grid"]["nx"]
    flags = ["velocity_preserving_qg_analogue_not_rsw_reproduction",
             "paper_core_regularized_origin_and_peak_normalized",
             "gauge_of_radial_pressure_audit_differs_from_evolution"]
    if np.max(relative_vorticity) >= 1:
        flags.append("max_abs_zeta_over_local_f_at_least_1_qg_not_asymptotically_valid_in_core")
    if np.max(relative_vorticity) >= .3:
        flags.append("max_abs_zeta_over_local_f_at_least_0p3_qg_marginal")
    if p["core_regularization_m"]/dx < 8:
        flags.append("fewer_than_8_grid_cells_across_core_regularization")
    if np.max(curvature_ratio) >= .1:
        flags.append("centrifugal_acceleration_at_least_10percent_of_coriolis_somewhere")
    if radial_chi is not None and np.max(np.abs(radial_chi)) >= .1:
        flags.append("radial_reference_interface_displacement_at_least_0p1_check_saved_gauge")
    metrics = {"status":"initial_circular_profile_audit_not_a_stability_test",
        "source_doi":"10.1016/j.icarus.2017.06.006", "samples":samples,
        "outer_radius_m":outer,"rotation_rate_s":config["planet"]["rotation_rate_s"],
        "planet_radius_m":config["planet"]["radius_m"],"f0_s":f0,
        "gamma_s_m2":base["gamma_code"],"deformation_radius_m":ld,
        "paper_reference_1000_over_f0_s":1000/f0,
        "core_regularization_m":p["core_regularization_m"],
        "regularized_core_peak_radius_m":p["cyclone_peak_radius_m"],
        "jet_peak_radius_m":p["jet_peak_radius_m"],
        "component_peak_jet_m_s":p["jet_speed_m_s"],
        "component_peak_cyclone_m_s":p["cyclone_speed_m_s"],
        "grid_spacing_m":dx,"core_regularization_cells":p["core_regularization_m"]/dx,
        "jet_width_cells":p["jet_width_m"]/dx,
        "max_abs_zeta_over_local_f":float(np.max(relative_vorticity)),
        "max_abs_U_over_r_local_f":float(np.max(curvature_ratio)),
        "max_abs_u_m_s":float(np.max(np.abs(base["u"]))),
        "q_r_over_r_min_s_m2":float(np.min(base["q_r_over_r"])),
        "q_r_over_r_max_s_m2":float(np.max(base["q_r_over_r"])),
        "q_r_changes_sign":bool(np.min(base["q_r_over_r"]) < 0 < np.max(base["q_r_over_r"])),
        "radial_pressure_reference":"psi(R)=psi_d(R)=0; untapered disk, NOT full-square zero-mean evolution gauge",
        "radial_max_abs_chi_relative_to_outer":float(np.max(np.abs(radial_chi))) if radial_chi is not None else None,
        "qg_f0_pressure_center_minus_outer_over_f0sq_LDsq":float(qg_pressure[0]) if qg_pressure is not None else None,
        "variable_f_pressure_center_minus_outer_over_f0sq_LDsq":float(localf_pressure[0]) if localf_pressure is not None else None,
        "centrifugal_pressure_center_minus_outer_over_f0sq_LDsq":float(centrifugal_pressure[0]) if centrifugal_pressure is not None else None,
        "flags":flags,
        "scope":"Analytic untapered initial profiles only; no spectral projection, nonaxisymmetric seed, dynamics, convergence, observed-frame correction, or saved-field thickness inferred."}
    columns={"r_m":r,"u_m_s":base["u"],"jet_u_m_s":base["u_jet"],
             "cyclone_u_m_s":base["u_cyclone"],"deep_u_m_s":base["u_deep"],
             "zeta_s":zeta,"local_f_s":f,"q_r_over_r_s_m2":base["q_r_over_r"],
             "abs_U_over_rf":curvature_ratio,"abs_zeta_over_f":relative_vorticity}
    if radial_chi is not None:
        columns["chi_radial_outer_reference"] = radial_chi
    return metrics,columns


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config",type=Path);parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--samples",type=int,default=16385);parser.add_argument("--outer-radius",type=float)
    args=parser.parse_args(argv)
    config=radial.load_config(args.config)
    metrics,columns=audit(config,args.samples,args.outer_radius)
    metrics["config_sha256"]=hashlib.sha256(args.config.read_bytes()).hexdigest()
    metrics["analysis_sha256"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__),Path(radial.__file__))}
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/"config.toml").write_bytes(args.config.read_bytes())
    (args.output/"audit.json").write_text(json.dumps(radial.clean_json(metrics),indent=2,allow_nan=False)+"\n")
    with (args.output/"profiles.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    (args.output/"COMPLETED").write_text("Initial radial profile audit only\n")
    print(json.dumps(radial.clean_json(metrics),indent=2,allow_nan=False))
    return metrics


if __name__=="__main__":
    main()

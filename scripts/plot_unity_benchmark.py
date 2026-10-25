"""Reproduce the review figures for Unity benchmark job 66041133 only.

This case-specific report never advances the simulation. It requires the
downloaded benchmark and the offline diagnostic_sampling.csv audit.

Requires numpy, h5py, matplotlib, toml (or Python 3.11+ tomllib), and the outputs
of inspect_saved_fields.jl. JLD2's native HDF5 matrix dimension order is reversed
when read with h5py: transpose to retain the model's [x,y] convention.
"""
import argparse
import base64
import csv
import hashlib
import json
from pathlib import Path
import copy
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import tomllib
    parse_toml = tomllib.loads
except ImportError:
    import toml
    parse_toml = toml.loads


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("analysis", type=Path)
    args = parser.parse_args()
    result, out = args.results.resolve(), args.analysis.resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = parse_toml((result / "config.toml").read_text())
    checks = parse_toml((out / "field_checks.toml").read_text())
    meta = parse_toml((result / "metadata.toml").read_text())
    assert result.parent.name == "job-66041133.AAVu7fFZ", "This report is specific to job 66041133"
    assert (result / "COMPLETED").read_text().strip() == "20000.0"
    sampling = read_csv(out / "diagnostic_sampling.csv")
    cubic_rows = [r for r in sampling if r["interpolation"] == "cubic_lagrange"]
    bilinear_rows = [r for r in sampling if r["interpolation"] == "bilinear"]
    cubic = next(r for r in cubic_rows if r["nr"] == "1024" and r["ntheta"] == "1440")
    cubic_initial = np.array([float(cubic["m%d_initial_m" % m]) for m in range(1,13)])
    cubic_final = np.array([float(cubic["m%d_final_m" % m]) for m in range(1,13)])
    cubic_growth = [float(r["percent_change"]) for r in cubic_rows]
    bilinear_growth = [float(r["percent_change"]) for r in bilinear_rows]
    diagnostics = read_csv(result / "diagnostics.csv")
    d = {k: np.array([float(row[k]) for row in diagnostics]) for k in diagnostics[0]}
    assert all(np.isfinite(v).all() for v in d.values())
    events = read_csv(result / "events.csv")
    evtime = np.array([float(e["time_s"]) for e in events])
    de = np.array([float(e["delta_energy"]) for e in events])
    dz = np.array([float(e["delta_total_enstrophy"]) for e in events])
    injected = np.array([e["injected"].lower() == "true" for e in events])
    times = d["time_s"]
    ce = np.array([sum(de[evtime <= t + 1e-9]) for t in times])
    cz = np.array([sum(dz[evtime <= t + 1e-9]) for t in times])
    energy, enstrophy = d["energy_m2_s2"], d["total_pv_enstrophy_s2"]
    er, zr = energy-energy[0]-ce, enstrophy-enstrophy[0]-cz
    phase = np.unwrap(np.angle(d["m6_real_m"] + 1j*d["m6_imag_m"]))
    omega = -np.polyfit(times,phase,1)[0]/6
    with h5py.File(str(result / "domain.jld2"), "r") as f:
        x, y = np.array(f["x"]), np.array(f["y"])
    snapshots = []
    for info in checks["snapshots"]:
        with h5py.File(str(result / info["file"]), "r") as f:
            snap = {k: np.array(f[k]).T for k in ("psi","s","u","v","totalpv")}
            snap["time_s"] = float(f["time_s"][()])
            snap["step"] = int(f["step"][()])
        assert all(np.isfinite(snap[k]).all() for k in ("psi","s","u","v","totalpv"))
        assert max(info["relative_inversion_errors"].values()) < 1e-9
        assert info["pv_identity_max_error_s"] < 1e-15
        j = np.flatnonzero(np.isclose(times, snap["time_s"], rtol=0, atol=1e-7))
        assert len(j) == 1
        calc_e = np.mean(snap["u"]**2+snap["v"]**2+
                         (snap["psi"]/config["physics"]["deformation_radius_m"])**2)/2
        assert np.isclose(calc_e, energy[j[0]], rtol=1e-12, atol=0)
        snap["zeta"] = snap["s"] + snap["psi"]/config["physics"]["deformation_radius_m"]**2
        snapshots.append(snap)
    assert checks["checkpoint"]["snapshot_s_max_abs_difference_s"] < 1e-15
    first, final = snapshots[0], snapshots[-1]
    traces = np.genfromtxt(str(out/"jet_traces.csv"),delimiter=",",names=True)
    profiles = np.genfromtxt(str(out/"radial_profiles.csv"),delimiter=",",names=True)
    theta = traces["theta_rad"]
    initial_radius = traces["radius_step_0_m"]
    final_radius = traces["radius_step_207_m"]
    mode_initial = np.array([d["m%d_amplitude_m" % m][0] for m in range(1,13)])
    mode_final = np.array([d["m%d_amplitude_m" % m][-1] for m in range(1,13)])
    rest = mode_final.copy(); rest[5] = 0
    metrics = {
        "job_id":66041133,"result_directory":str(result),"device":meta["device"],
        "duration_s":float(times[-1]),"steps":int(d["step"][-1]),
        "diagnostic_rows":len(times),"snapshots":len(snapshots),
        "events_attempted":len(events),"events_injected":int(injected.sum()),
        "events_suppressed":int((~injected).sum()),
        "energy_initial_m2_s2":float(energy[0]),"energy_final_m2_s2":float(energy[-1]),
        "energy_change_percent":float(100*(energy[-1]/energy[0]-1)),
        "net_impulse_energy_m2_s2":float(de.sum()),
        "energy_change_minus_impulses_m2_s2":float(er[-1]),
        "energy_residual_relative_to_initial":float(er[-1]/energy[0]),
        "enstrophy_change_percent":float(100*(enstrophy[-1]/enstrophy[0]-1)),
        "enstrophy_change_minus_impulses_s2":float(zr[-1]),
        "max_abs_mean_s_s":float(np.max(np.abs(d["mean_s_s"]))),
        "min_jet_coverage":float(d["jet_coverage"].min()),
        "initial_mean_peak_speed_m_s":float(d["mean_peak_speed_m_s"][0]),
        "final_mean_peak_speed_m_s":float(d["mean_peak_speed_m_s"][-1]),
        "mean_jet_radius_change_m":float(d["mean_jet_radius_m"][-1]-d["mean_jet_radius_m"][0]),
        "initial_m6_displacement_km":float(mode_initial[5]/1000),
        "final_m6_displacement_km":float(mode_final[5]/1000),
        "m6_amplitude_change_percent":float(100*(mode_final[5]/mode_initial[5]-1)),
        "runtime_largest_other_final_mode":int(np.argmax(rest)+1),
        "runtime_largest_other_final_amplitude_km":float(rest.max()/1000),
        "diagnostic_sampling_audit":{
            "settings_per_interpolant":len(cubic_rows),
            "cubic_selected_nr":1024,"cubic_selected_ntheta":1440,
            "cubic_initial_m6_km":float(cubic_initial[5]/1000),
            "cubic_final_m6_km":float(cubic_final[5]/1000),
            "cubic_growth_percent_range":[min(cubic_growth),max(cubic_growth)],
            "bilinear_growth_percent_range":[min(bilinear_growth),max(bilinear_growth)],
            "interpretation":"Remeasurements of identical fields, not a confidence interval or PDE convergence test. Weak competing modes are provisional."},
        "transient_mathematical_m6_angular_speed_rad_s":float(omega),
        "transient_mathematical_m6_speed_at_mean_jet_m_s":float(omega*d["mean_jet_radius_m"].mean()),
        "nominal_jet_circuits":float(times[-1]*config["initial"]["jet_speed_m_s"]/(2*np.pi*config["initial"]["jet_radius_m"])),
        "max_m6_change_between_bilinear_128x360_and_512x720_m":max(abs(s["m6_finer_sampling_m"]-s["m6_amplitude_m"]) for s in checks["snapshots"]),
        "field_check_summary":checks,
        "hashes":{name:hashlib.sha256((result/name).read_bytes()).hexdigest() for name in ["config.toml","Manifest.toml","diagnostics.csv","events.csv"]},
        "limits":["Seeded mode six; no emergence or selection test.",
                  "Weak modes depend on diagnostic interpolation; use the sampling audit before interpreting them.",
                  "Impulse residuals include hyperdiffusion and numerical error; no separate dissipation integral.",
                  "The phase fit uses actual output times and the model angular convention, not an observed System III drift.",
                  "This is not a timestep, spatial, boundary, or GPU restart convergence test."]}
    (out/"summary.json").write_text(json.dumps(metrics,indent=2,allow_nan=False)+"\n")

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
        "axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white",
        "savefig.facecolor":"white","pdf.fonttype":42})
    blue, orange = "#145E91", "#BD4936"
    fig, axes = plt.subplots(2,2,figsize=(12.2,9.4),constrained_layout=True)
    fig.suptitle("Unity benchmark: a weak, seeded sixfold jet",fontsize=17,fontweight="bold")
    radius = np.hypot(x[:,None],y[None,:])
    extent = [x[0]/1e6,x[-1]/1e6,y[0]/1e6,y[-1]/1e6]
    cmap=copy.copy(plt.get_cmap("RdBu_r")); cmap.set_bad("#eeeeee")
    maps=[(final["zeta"]*1e5,"A  Final relative vorticity",r"$\zeta$ ($10^{-5}$ s$^{-1}$)"),
          ((final["zeta"]-first["zeta"])*1e6,"B  Change from the initial flow",r"$\Delta\zeta$ ($10^{-6}$ s$^{-1}$)")]
    for ax,(field,title,label) in zip(axes[0],maps):
        mask=radius>24e6
        value=np.ma.masked_where(mask,field)
        limit=float(np.max(np.abs(value)))
        im=ax.imshow(value.T,origin="lower",extent=extent,cmap=cmap,vmin=-limit,vmax=limit,interpolation="nearest")
        ax.plot(final_radius*np.cos(theta)/1e6,final_radius*np.sin(theta)/1e6,color="#172732",linewidth=.8)
        ax.set(xlim=(-22,22),ylim=(-22,22),xlabel="x (1,000 km)",ylabel="y (1,000 km)",title=title)
        fig.colorbar(im,ax=ax,shrink=.82,pad=.02,label=label)
        ax.set_aspect("equal")
    ax=axes[1,0]
    ax.plot(np.degrees(theta),(initial_radius-initial_radius.mean())/1000,color=blue,label="Initial")
    ax.plot(np.degrees(theta),(final_radius-final_radius.mean())/1000,color=orange,label="20,000 s")
    ax.set(xlabel="Azimuth (degrees; model convention)",ylabel="Jet displacement from its mean radius (km)",
           title="C  The sixfold meander, measured from velocity",xlim=(0,360),xticks=np.arange(0,361,60))
    ax.axhline(0,color=".75",linewidth=.6); ax.grid(alpha=.18); ax.legend(frameon=False)
    ax=axes[1,1]
    modes=np.arange(1,13)
    for a,color,label in [(cubic_initial,blue,"Initial"),(cubic_final,orange,"20,000 s")]:
        values=np.where(a/1000>=.01,a/1000,np.nan)
        ax.plot(modes,values,"o-",color=color,label=label,markersize=4)
    ax.set(yscale="log",ylim=(.01,200),xticks=modes,xlabel="Azimuthal mode m",ylabel="Jet-displacement amplitude (km)",
           title="D  Mode six dominates (offline cubic measurement)")
    ax.grid(alpha=.18); ax.legend(frameon=False)
    ax.text(.03,.03,"Weak modes remain provisional; values <0.01 km omitted.",transform=ax.transAxes,fontsize=8,color=".35")
    for suffix in ("png","pdf"):
        fig.savefig(str(out/("flow_overview."+suffix)),dpi=180)
    plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(12.2,8.0),constrained_layout=True)
    fig.suptitle("What changed during 5.56 simulated hours?",fontsize=17,fontweight="bold")
    hours=times/3600
    ax=axes[0,0]
    ax.plot(hours,100*(energy/energy[0]-1),"o-",markersize=3,color=blue,label="Measured energy change")
    et=np.r_[0,evtime,times[-1]]/3600
    ec=np.r_[0,np.cumsum(de),de.sum()]/energy[0]*100
    ax.step(et,ec,where="post",color=orange,linestyle="--",label="Sum of impulse energy changes")
    ax.set(xlabel="Simulated time (hours)",ylabel="Change from initial energy (%)",title="A  Impulses explain almost all the energy change")
    ax.legend(frameon=False,fontsize=8); ax.grid(alpha=.18)
    ax=axes[0,1]
    ax.bar(np.arange(1,len(de)+1),de,color=[orange if z<0 else blue for z in de])
    ax.plot(np.flatnonzero(~injected)+1,de[~injected],"x",color="black",markersize=7,label="Suppressed attempt")
    ax.axhline(0,color=".5",linewidth=.7)
    ax.set(xlabel="Event attempt",ylabel=r"Impulse energy change (m$^2$ s$^{-2}$)",title="B  Nine injections remove energy; three add it",xticks=np.arange(1,len(de)+1))
    ax.legend(frameon=False,fontsize=8)
    ax=axes[1,0]
    ax.plot(hours,d["m6_amplitude_m"]/1000,"o-",color=blue,markersize=3,label="Runtime bilinear diagnostic")
    ax.plot(hours[[0,-1]], [cubic_initial[5]/1000,cubic_final[5]/1000],"x",color=orange,markersize=8,label="Offline cubic endpoints")
    ax.set(xlabel="Simulated time (hours)",ylabel="Mode-six displacement amplitude (km)",title="C  Seeded amplitude grows by about 3.3%")
    ax.grid(alpha=.18); ax.legend(frameon=False,fontsize=8)
    ax=axes[1,1]
    ax.plot(profiles["radius_m"]/1e6,profiles["utheta_step_0_m_s"],color=blue,label="Initial")
    ax.plot(profiles["radius_m"]/1e6,profiles["utheta_step_207_m_s"],"--",color=orange,label="20,000 s")
    ax.set(xlabel="Radius (1,000 km)",ylabel="Azimuthal-mean tangential speed (m/s)",title="D  The mean jet profile changes very little")
    ax.grid(alpha=.18); ax.legend(frameon=False)
    for suffix in ("png","pdf"):
        fig.savefig(str(out/("budgets_and_jet."+suffix)),dpi=180)
    plt.close(fig)

    rows=[("Completed integration","20,000 s; 207 steps; GPU"),
          ("Jet coverage","100% in every diagnostic record"),
          ("Mode-six amplitude, runtime diagnostic","99.851 → 103.098 km (+3.252%)"),
          ("Mode-six amplitude, offline cubic","%.3f → %.3f km; growth %.4f–%.4f%% across sampling settings"%(cubic_initial[5]/1000,cubic_final[5]/1000,min(cubic_growth),max(cubic_growth))),
          ("Mean peak jet speed","100.5683 → 100.5637 m/s"),
          ("Energy change","−0.0038768%; net energy removed by impulses"),
          ("Forcing","13 attempts, 12 injections, 1 suppressed"),
          ("Mean PV anomaly", "Maximum absolute mean: %.3e s⁻¹"%metrics["max_abs_mean_s_s"]),
          ("Continuous energy change after subtracting impulses","%.3e m²/s² (includes dissipation and numerical error)"%er[-1]),
          ("Diagnostic limitation","Weak competing modes depend on the interpolation and sampling settings; no physical interpretation yet")]
    html="""<!doctype html><html><meta charset="utf-8"><title>Saturn QG — Unity benchmark review</title>
<style>body{font:16px/1.55 system-ui,sans-serif;color:#203342;max-width:1100px;margin:40px auto;padding:0 25px}h1,h2{line-height:1.2}table{border-collapse:collapse;width:100%}td{border-bottom:1px solid #ddd;padding:10px}td:first-child{width:36%;font-weight:600}img{width:100%;margin:20px 0}.note{background:#edf4f8;padding:16px;border-left:4px solid #145e91}a{color:#145e91}</style>
<h1>Saturn polar QG: first 1024² Unity run</h1>
<p>Review of job 66041133. Saved fields, diagnostic/event CSVs, configuration, source and manifest were inspected locally. All three snapshots are finite, satisfy the saved PV identity and spectral inversion checks, and agree with recorded energy. The final checkpoint agrees with the final saved PV field. This checks the checkpoint contents, not resumed GPU evolution.</p>
<p class="note"><b>Result:</b> the short run is numerically consistent at the checks performed. It contains a weak, deliberately seeded sixfold meander. It does not establish stabilization by convection or spontaneous selection of mode six.</p><table>"""
    html+="".join("<tr><td>%s</td><td>%s</td></tr>"%row for row in rows)+"</table>"
    for imagefile in ("flow_overview.png","budgets_and_jet.png"):
        encoded=base64.b64encode((out/imagefile).read_bytes()).decode("ascii")
        html+='<img alt="Scientific benchmark figure" src="data:image/png;base64,'+encoded+'">'
    html+="""<h2>Interpretation and next steps</h2>
<p>The displacement is less than 0.7% of the 15,032-km jet radius, so the full polar image looks nearly circular. The displacement and modal panels reveal the imposed sixfold structure. The run covers only 0.021 nominal jet circuits.</p>
<p>The sixfold amplitude increase is supported by an independent offline remeasurement using the existing tensor-cubic interpolant. Across nine sampling settings its growth is 3.2546–3.2620%, with less than 10 m variation in either endpoint amplitude. Oversampling the runtime bilinear interpolant gives non-monotone results (2.40–3.26% growth), so simply increasing the sampling density does not establish convergence. These are measurements of the same saved fields, not a statistical confidence interval or a test of the simulation's spatial or temporal convergence. Weak competing modes must remain provisional. The initial m=12 harmonic persists under the cubic audit and is not evidence of a newly selected mode. The modal plot uses the 1024-radial, 1440-angular cubic measurement; the jet traces use bilinear sampling.</p>
<p>Nine injections remove energy and three add it. A PV impulse interacts with the existing flow, so its energy change can have either sign; negative changes are not by themselves an implementation failure. Over this window, however, the forcing is a net energy sink. Calling it sustained positive energy input would be inaccurate. Compare forced and unforced trajectories at comparable mean jet strength before attributing stabilization.</p>
<p>The next prepared run halves the timestep to 50 s over the same 20,000 s. Compare actual output times and final fields, mode amplitude/phase and injection records. GPU restart equivalence should also be checked at matched timestep boundaries. The saved intermediate snapshot is at 10,028.955 s, not exactly 10,000 s; terminating a separate run at a new time can change the numerical step sequence.</p>
<p>A matched forced/unforced pair for 1,000,000 s (about 1.06 nominal circuits) is prepared, subject to the numerical checks. It starts each case from the same initialized velocity field. A seeded maintenance comparison does not test emergence. Further controls on deformation radius and deep flow, additional seeds, and longer durations follow any promising result.</p>
<p>Spatial sensitivity needs careful design. The present 250-km forcing radius spans 4.27 cells. The forcing center separation is 2r+dx, so changing resolution also changes the pair geometry. Holding the cutoff damping time fixed changes the physical hyperviscosity. A simple grid doubling therefore is not a fixed-physics convergence test.</p>
<p>Reproducibility: see summary.json, field_checks.toml, diagnostic_sampling.csv, diagnostic_sampling_metadata.toml, jet_traces.csv and radial_profiles.csv in this directory. Timing was user-reported as 56 s including startup and output; host MaxRSS was about 1.54 GiB. No GPU-memory or steady-throughput measurement is inferred from those values.</p></html>"""
    (out/"report.html").write_text(html,encoding="utf-8")
    print(json.dumps({k:v for k,v in metrics.items() if k not in ("field_check_summary","hashes","limits")},indent=2))


if __name__ == "__main__":
    main()

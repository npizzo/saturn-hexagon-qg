"""Review saved forced/unforced maintenance runs; do not evolve the model.

First run analyze_maintenance_jet.jl and inspect_saved_fields.jl. Uses numpy,
h5py, matplotlib, and toml (or tomllib). All field axes follow Julia [x,y].
Missing jet-mode measurements remain missing; they are never replaced by zero.
"""
import argparse
import base64
import json
import re
from pathlib import Path
import numpy as np
import h5py
import matplotlib.pyplot as plt
from compare_timesteps import read_csv, digest, config_diff, snapshot, parse_toml


def load_run(path, label, out):
    c = parse_toml((path / "config.toml").read_text())
    m = parse_toml((path / "metadata.toml").read_text())
    rows = read_csv(path / "diagnostics.csv")
    d = {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}
    assert float((path / "COMPLETED").read_text()) == d["time_s"][-1] == c["numerics"]["duration_s"]
    for key, values in d.items():
        missing = ~np.isfinite(values)
        if missing.any():
            assert re.match(r"m\d+_", key) and np.all(d["jet_coverage"][missing] < 1), key
    assert digest(path / "Manifest.toml") == m["manifest_sha256"]
    from hashlib import sha256
    source = "\n".join(p.name + "\n" + p.read_text() for p in sorted((path / "src").glob("*.jl")))
    assert sha256(source.encode()).hexdigest() == m["source_sha256"]
    checks = parse_toml((out / ("field-checks-" + label) / "field_checks.toml").read_text())
    assert checks["source_matches"]
    assert checks["checkpoint"]["snapshot_s_max_abs_difference_s"] < 1e-15
    for item in checks["snapshots"]:
        assert item["all_fields_finite"]
        assert max(item["relative_inversion_errors"].values()) < 1e-9
        assert item["pv_identity_max_error_s"] < 1e-15
    events_path = path / "events.csv"
    events = read_csv(events_path) if events_path.exists() else []
    snapshots = sorted(path.glob("snapshot_*.jld2"))
    assert len(snapshots) == len(checks["snapshots"]) == 11
    with h5py.File(str(path / "domain.jld2"), "r") as f:
        domain = {key: np.array(f[key]).T for key in ("x", "y", "background_pv", "psi_deep", "reference_s")}
    return dict(path=path, config=c, metadata=m, d=d, events=events,
                snapshots=snapshots, checks=checks, domain=domain)


def field_statistics(f, c, x, y):
    radius = np.hypot(x[:, None], y[None, :])
    outer = radius >= c["boundary"]["inner_radius_m"]
    ke = (f["u"]**2 + f["v"]**2) / 2
    pe = (f["psi"] / c["physics"]["deformation_radius_m"])**2 / 2
    zeta = f["s"] + f["psi"] / c["physics"]["deformation_radius_m"]**2
    omega, planet_radius = c["planet"]["rotation_rate_s"], c["planet"]["radius_m"]
    coriolis = 2*omega - omega*(radius/planet_radius)**2
    return dict(kinetic_energy_m2_s2=float(ke.mean()), stretching_energy_m2_s2=float(pe.mean()),
                max_speed_m_s=float(np.hypot(f["u"], f["v"]).max()),
                kinetic_energy_fraction_outside_science_radius=float(ke[outer].sum()/ke.sum()),
                max_speed_outside_science_radius_m_s=float(np.hypot(f["u"], f["v"])[outer].max()),
                max_abs_zeta_over_local_f_inside_science_radius=float(np.max(np.abs(zeta[~outer])/coriolis[~outer])))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("forced", type=Path)
    p.add_argument("unforced", type=Path)
    p.add_argument("restart", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    out = args.output.resolve()
    a, b = [load_run(path.resolve(), label, out) for path, label in
            ((args.forced, "forced"), (args.unforced, "unforced"))]
    diff = config_diff(a["config"], b["config"])
    assert set(diff) == {"convection.enabled", "experiment.name"}, diff
    assert a["config"]["convection"]["enabled"] and not b["config"]["convection"]["enabled"]
    for key in ("source_sha256", "manifest_sha256", "julia_version", "device"):
        assert a["metadata"][key] == b["metadata"][key]
    assert all(np.array_equal(a["domain"][k], b["domain"][k]) for k in a["domain"])
    initial_a, _, _ = snapshot(a["snapshots"][0])
    initial_b, _, _ = snapshot(b["snapshots"][0])
    assert all(np.array_equal(initial_a[k], initial_b[k]) for k in initial_a)
    restart = parse_toml((args.restart / "report.toml").read_text())
    assert restart["device"] == "gpu" and restart["passed"]
    assert (args.restart / "RESTART_VALIDATION_PASSED").is_file()
    assert all(item["relative_l2"] <= restart["relative_l2_tolerance"] for item in restart["final_errors"].values())
    assert restart["source_sha256"] == a["metadata"]["source_sha256"]
    assert restart["manifest_sha256"] == a["metadata"]["manifest_sha256"]
    polar = read_csv(out / "polar_cubic.csv")
    profiles = read_csv(out / "profiles.csv")
    assert (out / "polar_cubic_metadata.toml").is_file()
    assert len(polar) == 48
    def polar_rows(run, scope="configured", nr="512"):
        return sorted([r for r in polar if r["run_id"] == run["path"].parent.name
                       and r["scope"] == scope and r["nr"] == nr], key=lambda r: float(r["time_s"]))
    def profile(run, time):
        rows = [r for r in profiles if r["run_id"] == run["path"].parent.name and
                r["scope"] == "wide_sensitivity" and r["nr"] == "512" and float(r["time_s"]) == time]
        assert len(rows) == 512
        return {key: np.array([float(r[key]) for r in rows]) for key in ("radius_m", "mean_utheta_m_s")}
    summary = dict(config_differences=diff, identical_initial_fields=True,
                   identical_domains=True, matching_source_dependencies=True,
                   gpu_restart_passed=True, gpu_restart_max_relative_difference=max(v["relative_l2"] for v in restart["final_errors"].values()),
                   runs={}, limitations=["One parameter choice and one forced realization, with m6 imposed initially.",
                       "20,000-second timestep and restart tests do not establish nonlinear million-second convergence.",
                       "Late global displacement modes are unavailable; wider searches switch between different ridges.",
                       "Only convection is toggled; this pair does not isolate deformation radius or prescribed deep flow.",
                       "Spatial, dissipation, forcing-radius and boundary sensitivity remain open."])
    maps = {}
    for label, run in (("forced", a), ("unforced", b)):
        d = run["d"]
        rows = polar_rows(run)
        wide = polar_rows(run, "wide_sensitivity")
        assert len(rows) == len(wide) == 11
        bad = np.flatnonzero(d["jet_coverage"] < 1)
        assert len(bad) and np.all(~np.isfinite(d["m6_amplitude_m"][bad]))
        valid = np.isfinite(d["m6_amplitude_m"])
        assert np.all(valid[:bad[0]]) and not np.any(valid[bad[0]:])
        energy = d["energy_m2_s2"]
        z = d["total_pv_enstrophy_s2"]
        injected = [r for r in run["events"] if r["injected"] == "true"]
        suppressed = [r for r in run["events"] if r["injected"] != "true"]
        event_e = sum(float(r["delta_energy"]) for r in run["events"])
        event_z = sum(float(r["delta_total_enstrophy"]) for r in run["events"])
        initial, _, _ = snapshot(run["snapshots"][0])
        final, final_time, _ = snapshot(run["snapshots"][-1])
        x,y = run["domain"]["x"],run["domain"]["y"]
        final_wide = wide[-1]
        summary["runs"][label] = dict(path=str(run["path"]),
            duration_s=float(d["time_s"][-1]), steps=int(d["step"][-1]),
            energy_change_percent=float(100*(energy[-1]/energy[0]-1)),
            enstrophy_change_percent=float(100*(z[-1]/z[0]-1)),
            mean_s_max_absolute_s=float(np.max(np.abs(d["mean_s_s"]))),
            first_incomplete_runtime_coverage_s=float(d["time_s"][bad[0]]),
            last_valid_runtime_m6_time_s=float(d["time_s"][bad[0]-1]),
            last_valid_runtime_m6_amplitude_m=float(d["m6_amplitude_m"][bad[0]-1]),
            final_runtime_jet_coverage=float(d["jet_coverage"][-1]),
            final_runtime_m6_amplitude_m=None,
            final_wide_mean_profile_peak_m_s=float(final_wide["mean_profile_peak_speed_m_s"]),
            final_wide_mean_profile_peak_radius_m=float(final_wide["mean_profile_peak_radius_m"]),
            final_wide_profile_envelope_fwhm_m=float(final_wide["mean_profile_fwhm_m"]),
            final_wide_trace_max_adjacent_jump_m=float(final_wide["max_adjacent_valid_radius_jump_m"]),
            attempted_events=len(run["events"]), injected_events=len(injected), suppressed_events=len(suppressed),
            negative_energy_injections=sum(float(r["delta_energy"])<0 for r in injected),
            positive_energy_injections=sum(float(r["delta_energy"])>0 for r in injected),
            logged_event_energy_change_m2_s2=event_e,
            logged_suppressed_event_energy_change_m2_s2=sum(float(r["delta_energy"]) for r in suppressed),
            energy_change_minus_logged_events_m2_s2=float(energy[-1]-energy[0]-event_e),
            enstrophy_change_minus_logged_events_s2=float(z[-1]-z[0]-event_z),
            initial_fields=field_statistics(initial,run["config"],x,y),
            final_fields=field_statistics(final,run["config"],x,y),
            hashes={name:digest(run["path"]/name) for name in ("config.toml","Manifest.toml","diagnostics.csv")})
        maps[label] = []
        for idx in (0,2,5,10):
            fields,time,_ = snapshot(run["snapshots"][idx])
            maps[label].append((time,fields["s"]+fields["psi"]/run["config"]["physics"]["deformation_radius_m"]**2))
    times=[100000,300000,500000,570000]
    shared=[]
    for time in times:
        vals=[]
        for run in (a,b):
            valid=np.isfinite(run["d"]["m6_amplitude_m"])
            t=run["d"]["time_s"][valid]
            assert t[0] <= time <= t[-1]
            vals.append(float(np.interp(time,t,run["d"]["m6_amplitude_m"][valid])))
        shared.append(dict(time_s=time,forced_m=vals[0],unforced_m=vals[1],difference_percent=100*(vals[0]/vals[1]-1)))
    summary["runtime_m6_comparison_linearly_interpolated_in_valid_interval"]=shared
    summary["conclusion"]="Both flows amplify the seeded sixfold pattern and roll up into lobes and vortices; this pair does not show a convection-maintained stable hexagonal jet. Late tracker coverage is not a physical lifetime measure."
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,
                         "axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
    blue,orange="#145E91","#BD4936"
    fig,axs=plt.subplots(2,4,figsize=(15.5,8.0),constrained_layout=True)
    fig.suptitle("Both flows amplify the seeded pattern and roll up",fontsize=17,fontweight="bold")
    x,y=a["domain"]["x"],a["domain"]["y"]
    mask=np.hypot(x[:,None],y[None,:])>24e6
    limit=max(np.max(np.abs(z[~mask])) for fields in maps.values() for time,z in fields)*1e5
    for i,label in enumerate(("forced","unforced")):
        for j,(time,z) in enumerate(maps[label]):
            im=axs[i,j].imshow(np.ma.masked_where(mask,z*1e5).T,origin="lower",
                 extent=[x[0]/1e6,x[-1]/1e6,y[0]/1e6,y[-1]/1e6],cmap="RdBu_r",vmin=-limit,vmax=limit,interpolation="nearest")
            axs[i,j].set(xlim=(-24,24),ylim=(-24,24),xticks=[-20,0,20],yticks=[-20,0,20],
                xlabel="x (1,000 km)",ylabel="y (1,000 km)",
                title=("Convection on" if i==0 else "Convection off")+" · %.2f days"%(time/86400))
    fig.colorbar(im,ax=axs,label=r"Relative vorticity ($10^{-5}$ s$^{-1}$)",shrink=.75,pad=.015)
    for suffix in ("png","pdf"):fig.savefig(str(out/("flow_evolution."+suffix)),dpi=180)
    plt.close(fig)

    fig,axs=plt.subplots(2,3,figsize=(15,8.6),constrained_layout=True)
    fig.suptitle("Convection changes the details; strong meander growth occurs in both",fontsize=16,fontweight="bold")
    for label,run,color in (("Convection on",a,orange),("Convection off",b,blue)):
        d=run["d"];hours=d["time_s"]/86400
        axs[0,0].plot(hours,d["m6_amplitude_m"]/1000,color=color,label=label)
        pr=polar_rows(run)
        axs[0,0].plot([float(r["time_s"])/86400 for r in pr],[float(r["m6_amplitude_m"])/1000 for r in pr],"o",color=color,ms=4)
        axs[0,1].plot(hours,d["jet_coverage"],color=color)
        pf=profile(run,d["time_s"][-1])
        axs[0,2].plot(pf["radius_m"]/1e6,pf["mean_utheta_m_s"],color=color,label=label+", final")
        axs[1,0].plot(hours,100*(d["energy_m2_s2"]/d["energy_m2_s2"][0]-1),color=color,label=label)
        valid=np.isfinite(d["m6_real_m"])
        phase=np.unwrap(np.angle(d["m6_real_m"][valid]+1j*d["m6_imag_m"][valid]))
        axs[1,1].plot(hours[valid],np.degrees(-(phase-phase[0])/6),color=color)
        mid=pr[5]
        amps=np.array([float(mid["m%d_amplitude_m"%m])/1000 for m in range(1,13)])
        axs[1,2].plot(range(1,13),np.where(amps>=.01,amps,np.nan),"o-",color=color,label=label,ms=4)
    initial=profile(a,0.0)
    axs[0,2].plot(initial["radius_m"]/1e6,initial["mean_utheta_m_s"],"--",color=".35",label="Common initial profile")
    axs[0,0].set(title="A  Sixfold displacement while measurable",xlabel="Simulated time (days)",ylabel="Displacement amplitude (km)",xlim=(0,1e6/86400))
    axs[0,0].legend(frameon=False,fontsize=8)
    axs[0,0].text(.03,.95,"Dots: independent cubic extraction",transform=axs[0,0].transAxes,va="top",fontsize=8)
    axs[0,1].set(title="B  Coverage of the original jet search",xlabel="Simulated time (days)",ylabel="Fraction of angles with valid peak",ylim=(0,1.05))
    axs[0,2].set(title="C  Mean flow spreads over a broad radial band",xlabel="Radius (1,000 km)",ylabel="Azimuthal-mean tangential speed (m/s)")
    axs[0,2].legend(frameon=False,fontsize=8)
    event_time=np.array([float(r["time_s"]) for r in a["events"]])
    cum=np.cumsum([float(r["delta_energy"]) for r in a["events"]])
    axs[1,0].plot(event_time/86400,100*cum/a["d"]["energy_m2_s2"][0],"--",color=".3",lw=1,label="Logged event updates")
    axs[1,0].set(title="D  Forcing removes energy overall",xlabel="Simulated time (days)",ylabel="Change from initial QG energy (%)")
    axs[1,0].legend(frameon=False,fontsize=8)
    axs[1,1].set(title="E  Pattern orientation before tracker loss",xlabel="Simulated time (days)",ylabel="Angle relative to initial pattern (degrees)",xlim=(0,1e6/86400))
    axs[1,2].set(title="F  Modes near 500,000 s (cubic extraction)",xlabel="Azimuthal mode m",ylabel="Displacement amplitude (km)",yscale="log",xticks=[1,2,4,6,8,10,12],ylim=(.01,5000))
    axs[1,2].text(.03,.03,"Values below 0.01 km omitted",transform=axs[1,2].transAxes,fontsize=8)
    for ax in axs.flat:ax.grid(alpha=.16)
    for suffix in ("png","pdf"):fig.savefig(str(out/("maintenance_diagnostics."+suffix)),dpi=180)
    plt.close(fig)

    f,u=summary["runs"]["forced"],summary["runs"]["unforced"]
    rows=[("Integration", "1,000,000 s; 10,332 steps", "1,000,000 s; 10,000 steps"),
          ("A6 near 500,000 s", "≈2,005 km", "≈2,030 km"),
          ("First incomplete runtime jet coverage", "%.3f days"%(f["first_incomplete_runtime_coverage_s"]/86400),"%.3f days"%(u["first_incomplete_runtime_coverage_s"]/86400)),
          ("Final runtime coverage", "%.1f%%"%(100*f["final_runtime_jet_coverage"]),"%.1f%%"%(100*u["final_runtime_jet_coverage"])),
          ("Final global A6", "Unavailable", "Unavailable"),
          ("Final azimuthal-mean speed peak within 7–23 Mm", "%.2f m/s"%f["final_wide_mean_profile_peak_m_s"],"%.2f m/s"%u["final_wide_mean_profile_peak_m_s"]),
          ("Total QG energy change", "%.6f%%"%f["energy_change_percent"],"%.6f%%"%u["energy_change_percent"]),
          ("Total-PV enstrophy change", "%.4f%%"%f["enstrophy_change_percent"],"%.4f%%"%u["enstrophy_change_percent"]),
          ("Convection", "673 attempts; 482 injections; 191 suppressed", "Disabled")]
    body="""<!doctype html><html><meta charset="utf-8"><title>Saturn QG: first maintenance pair</title>
<style>body{font:16px/1.55 system-ui,sans-serif;color:#203342;max-width:1250px;margin:40px auto;padding:0 25px}h1,h2{line-height:1.2}table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid #ddd;padding:10px;text-align:left}td:first-child{width:40%;font-weight:600}img{width:100%;margin:20px 0}.note{background:#edf4f8;padding:16px;border-left:4px solid #145e91}a{color:#145e91}</style>
<h1>First million-second Saturn QG maintenance pair</h1>
<p>Convection on: Unity 66043704. Convection off: 66043705. GPU restart check: 66043703.</p>
<p class="note"><b>Scientific result:</b> both flows amplify the imposed sixfold meander and develop rolled-up lobes and vortices. At this parameter choice, convection does not prevent that nonlinear evolution. This pair does not show a stable, convection-maintained hexagonal jet, and it does not rule out such a state elsewhere in parameter space.</p>
<table><tr><th>Quantity</th><th>Convection on</th><th>Convection off</th></tr>"""
    body+="".join("<tr>"+"".join("<td>"+v+"</td>" for v in row)+"</tr>" for row in rows)+"</table>"
    for name in ("flow_evolution.png","maintenance_diagnostics.png"):
        body+='<img alt="Scientific maintenance-pair comparison" src="data:image/png;base64,'+base64.b64encode((out/name).read_bytes()).decode()+ '">'
    body+="""<h2>Controlled comparison and numerical checks</h2>
<p>Initial physical fields and saved domains are exactly identical. The configuration differs only in the experiment label and convection toggle; source, pinned dependencies and Julia version match. Both runs retain finite physical fields in all 22 saved snapshots, satisfy the PV/inversion checks and agree with the recorded invariants. Final checkpoints reconstruct the saved PV to about 4×10⁻¹⁶ relative accuracy. The GPU restart check reproduced the reference fields, RNG and event schedule exactly through 20,000 s. It tested a new model restored on the same device/allocation; the million-second nonlinear evolution has not yet passed a comparable timestep or spatial sensitivity check. Scheduler completion was supplied by the user; raw batch logs were not downloaded.</p>
<h2>What happened to the jet?</h2>
<p>The initial 100-km sixfold displacement grows to roughly 2,000 km by 500,000 s in both runs. Using linear interpolation only within the valid diagnostic interval, the forced amplitude is about 1.2% smaller at 500,000 s. Independent cubic extraction of nearby saved fields supports the strong common growth. Convection introduces additional asymmetry and modifies the roll-up; an early amplitude difference is not evidence of a permanently stabilized state. Only convection is toggled here, so the pair does not isolate the effect of finite deformation radius or the prescribed deep jet.</p>
<p>The original 11–19 Mm search loses full angular coverage after roughly 6.7–6.9 days. Late mode values are intentionally missing rather than zero. Widening the search to 7–23 Mm restores numerical coverage but produces adjacent trace jumps of about 12,538 km forced and 10,249 km unforced. Those jumps persist at doubled sampling: the algorithm is switching between different velocity ridges, not identifying one continuous hexagonal jet. Consequently neither late wide-band A6 nor the late conditional average of detected local peak speeds is a reliable whole-jet metric. The earlier loss of coverage in the forced run is not a demonstrated decrease in physical lifetime.</p>
<p>The full-circle mean velocity profile in the jet region broadens and its peak decreases from about 100 m/s to 28.6 m/s forced and 32.0 m/s unforced. The unforced final peak lies near 10.35 Mm, inside the original search's inner bound. The broad mean-profile envelope spans about 8.8 Mm at half height, compared with 1.76 Mm initially; this is not the thickness of each locally meandering or rolled-up jet segment. Total kinetic energy also falls in both cases while the stretching-energy contribution increases. The total QG energy is nearly conserved in the unforced run, so mean-profile weakening primarily describes flow redistribution rather than loss of the entire flow's energy.</p>
<h2>Convection and budgets</h2>
<p>Of 482 actual injections, 287 remove QG energy and 195 add it. The logged event updates sum to −1.28953647 m²/s², compared with a total energy change of −1.29547372. The residual is −0.00593725, versus −0.00534988 in the unforced run. These residuals contain configured hyperdiffusion and numerical effects; no independent time integral of dissipation was recorded. Suppressed events still pass through the model's spectral-state update and contribute a tiny aggregate logged energy change (−1.73×10⁻⁷). Thus the event-update budget should not be attributed entirely to nonzero physical injections. The forcing is a net energy sink in this realization. Mean prognostic PV remains negligible.</p>
<p>Mode-six phase can be followed continuously only before coverage loss. Its propagation is time dependent, and its model angular convention has not been calibrated to an observed Saturn reference frame. Mode six dominates the measured displacement modes before loss of coverage; this is not mode-selection evidence because mode six was seeded. Higher harmonics and rolled-up vortices do not, by themselves, establish a maintained hexagonal jet.</p>
<h2>Next scientific checks</h2>
<p>The primary next question is whether the lower layer supports a persistent finite hexagon with the observed pattern drift and upper winds. Stability and phase speed must be tested together. The current rotation parameter corresponds to a 10 h 34 m 40 s frame, not System III; register wind and pattern-speed comparisons consistently. Hold the upper velocity fixed in the circular screen and vary deformation radius and the prescribed deep profile without convection. A separate circular-base eigenvalue calculation can screen growth rates cheaply; it cannot establish stability of a finite-amplitude hexagon. Use the Pizzo–Salmon (2024) Ptolemaic family as a reference after checking its formula conventions, Eulerian reconstruction and ambient patch. A verified stationary barotropic equilibrium remains an equilibrium at finite deformation radius over resting deep flow, but a general circular deep jet changes the equilibrium and requires continuation or a residual check. The prepared dt=50 s repeat remains useful for validating the old pilot's nonlinear evolution. Broadband selection and convection should follow the equilibrium/stability controls rather than replace them.</p>
<p>Spatial and dissipation sensitivity still need a deliberate design: fixed cutoff damping time changes the physical hyperviscosity when resolution changes, and the convective pair separation currently includes a grid-spacing contribution. Late flow also reaches the polar-background continuation region: outside 24 Mm, peak speed is about 12.3 m/s forced and 9.7 m/s unforced, although only 0.36% and 0.32% of kinetic energy lies there. This warrants a boundary-sensitivity check, not a claim of boundary contamination. The maximum |relative vorticity|/local f inside 24 Mm rises from about 0.26 initially to 0.35 forced and 0.28 unforced; the QG approximation and illustrative jet width should be revisited when the observed mean velocity profile is available.</p>
<p>Reproducibility: summary.json, checkpoint-artifact-audit.toml, field-checks-forced/, field-checks-unforced/, polar_cubic.csv, jettraces.csv, profiles.csv and polar_cubic_metadata.toml accompany this report. The montage uses one shared vorticity color scale, actual saved times rounded to days, and displays only radii ≤24 Mm. Analysis scripts do not modify the model or saved data.</p></html>"""
    (out/"report.html").write_text(body,encoding="utf-8")
    print(json.dumps({"conclusion":summary["conclusion"],"runs":summary["runs"]},indent=2))


if __name__ == "__main__":
    main()

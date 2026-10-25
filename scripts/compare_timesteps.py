"""Compare two saved convective runs differing only in timestep; never advance the model.

Dependencies: numpy, h5py, matplotlib and toml (or Python >=3.11 tomllib).
JLD2 matrix axes are transposed on HDF5 input to recover the Julia [x,y] order.
All comparisons use recorded times; diagnostic row indices need not align.
"""
import argparse
import base64
import csv
import hashlib
import html
import json
from pathlib import Path
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
    with path.open() as f:
        return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config_diff(a, b, prefix=""):
    result = {}
    for key in sorted(set(a) | set(b)):
        name = prefix + key
        if key not in a or key not in b:
            result[name] = [a.get(key), b.get(key)]
        elif isinstance(a[key], dict) and isinstance(b[key], dict):
            result.update(config_diff(a[key], b[key], name + "."))
        elif a[key] != b[key]:
            result[name] = [a[key], b[key]]
    return result


def snapshot(path):
    with h5py.File(str(path), "r") as f:
        fields = {k: np.array(f[k]).T for k in ("psi", "s", "u", "v", "totalpv")}
        time = float(f["time_s"][()])
        step = int(f["step"][()])
    assert all(np.isfinite(v).all() for v in fields.values()), path
    return fields, time, step


def load_run(path):
    config = parse_toml((path / "config.toml").read_text())
    metadata = parse_toml((path / "metadata.toml").read_text())
    assert digest(path / "Manifest.toml") == metadata["manifest_sha256"]
    names = sorted((path / "src").glob("*.jl"))
    source = "\n".join(p.name + "\n" + p.read_text() for p in names)
    assert hashlib.sha256(source.encode()).hexdigest() == metadata["source_sha256"]
    records = read_csv(path / "diagnostics.csv")
    diag = {k: np.array([float(r[k]) for r in records]) for k in records[0]}
    assert all(np.isfinite(v).all() for v in diag.values())
    snapshots = {}
    for name in sorted(path.glob("snapshot_*.jld2")):
        fields, time, step = snapshot(name)
        assert time not in snapshots
        snapshots[time] = {"fields": fields, "step": step, "file": name.name}
    final_time = max(snapshots)
    assert final_time == float((path / "COMPLETED").read_text())
    assert final_time == config["numerics"]["duration_s"] == diag["time_s"][-1]
    with h5py.File(str(path / "domain.jld2"), "r") as f:
        domain = {k: np.array(f[k]).T for k in ("x", "y", "background_pv", "psi_deep", "reference_s")}
    return {"path": path, "config": config, "metadata": metadata, "diag": diag,
            "events": read_csv(path / "events.csv"), "snapshots": snapshots, "domain": domain}


def compare_fields(a, b, initial):
    metrics = {}
    for key in a:
        delta = b[key] - a[key]
        evolved = np.linalg.norm(b[key] - initial[key])
        metrics[key] = {
            "relative_l2": float(np.linalg.norm(delta) / np.linalg.norm(b[key])),
            "max_abs": float(np.max(np.abs(delta))),
            "l2_difference_over_evolution": float(np.linalg.norm(delta) / evolved) if evolved else None}
    du, dv = b["u"] - a["u"], b["v"] - a["v"]
    metrics["velocity_vector"] = {
        "relative_l2": float(np.sqrt(np.sum(du**2 + dv**2) / np.sum(b["u"]**2 + b["v"]**2))),
        "max_pointwise_difference_m_s": float(np.max(np.hypot(du, dv)))}
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("refined", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    a, b = load_run(args.baseline.resolve()), load_run(args.refined.resolve())
    out = args.output.resolve()
    assert out != a["path"] and out != b["path"]
    out.mkdir(parents=True, exist_ok=True)
    polar_audit = read_csv(out / "jet_measurements.csv")
    cubic = next(row for row in polar_audit if row["phase"] == "final"
                 and row["interpolation"] == "cubic_lagrange"
                 and row["nr"] == "1024" and row["ntheta"] == "1440")
    cubic_metadata = parse_toml((out / "jet_measurement_metadata.toml").read_text())
    assert Path(cubic_metadata["base_results"]) == a["path"]
    assert Path(cubic_metadata["refined_results"]) == b["path"]
    differences = config_diff(a["config"], b["config"])
    assert set(differences) == {"numerics.dt_s"}, differences
    dt_a, dt_b = differences["numerics.dt_s"]
    assert 0 < dt_b < dt_a
    for key in ("source_sha256", "manifest_sha256", "julia_version", "device"):
        assert a["metadata"][key] == b["metadata"][key], key
    assert digest(a["path"] / "Project.toml") == digest(b["path"] / "Project.toml")
    assert all(np.array_equal(a["domain"][k], b["domain"][k]) for k in a["domain"])
    assert 0.0 in a["snapshots"] and 0.0 in b["snapshots"]
    initial_a = a["snapshots"][0.0]["fields"]
    initial_b = b["snapshots"][0.0]["fields"]
    assert all(np.array_equal(initial_a[k], initial_b[k]) for k in initial_a)
    final_time = max(a["snapshots"])
    assert final_time == max(b["snapshots"])
    fields_a, fields_b = a["snapshots"][final_time]["fields"], b["snapshots"][final_time]["fields"]
    fields = compare_fields(fields_a, fields_b, initial_b)
    matched_snapshots = sorted(set(a["snapshots"]) & set(b["snapshots"]))
    snapshot_metrics = {str(t): compare_fields(a["snapshots"][t]["fields"], b["snapshots"][t]["fields"], initial_b)
                        for t in matched_snapshots}
    da, db = a["diag"], b["diag"]
    matched_times, ia, ib = np.intersect1d(da["time_s"], db["time_s"], return_indices=True)
    ea, eb = a["events"], b["events"]
    assert len(ea) == len(eb)
    schedule_keys = ("attempted", "injected", "event", "time_s", "next_event_s", "x_m", "y_m",
                     "cyclone_x_m", "cyclone_y_m", "anticyclone_x_m", "anticyclone_y_m")
    assert all(x[k] == y[k] for x, y in zip(ea, eb) for k in schedule_keys)
    qmax = a["config"]["convection"]["q_max_s"]
    dq = max(abs(float(y[k])-float(x[k])) for x,y in zip(ea,eb)
             for k in ("q_cyclone_s", "q_anticyclone_s"))
    impulse_a = sum(float(e["delta_energy"]) for e in ea)
    impulse_b = sum(float(e["delta_energy"]) for e in eb)
    ca = complex(da["m6_real_m"][-1], da["m6_imag_m"][-1])
    cb = complex(db["m6_real_m"][-1], db["m6_imag_m"][-1])
    phase = float(np.angle(cb / ca))
    summary = {
        "baseline":str(a["path"]), "refined":str(b["path"]), "configuration_differences":differences,
        "identical_initial_fields":True, "identical_domains":True, "matching_source_dependencies_device":True,
        "duration_s":final_time, "steps":[int(da["step"][-1]), int(db["step"][-1])],
        "matched_snapshot_times_s":matched_snapshots, "matched_diagnostic_times_s":matched_times.tolist(),
        "diagnostic_rows":[len(da["time_s"]),len(db["time_s"])],
        "final_field_comparison":fields, "matched_snapshot_comparisons":snapshot_metrics,
        "final_runtime_m6_amplitude_m":[float(da["m6_amplitude_m"][-1]),float(db["m6_amplitude_m"][-1])],
        "final_runtime_m6_difference_m":float(db["m6_amplitude_m"][-1]-da["m6_amplitude_m"][-1]),
        "final_runtime_m6_coefficient_phase_difference_rad":phase,
        "final_runtime_pattern_angle_difference_rad":-phase/6,
        "final_cubic_m6_difference_m":float(cubic["A6_refined_minus_base_m"]),
        "final_cubic_m6_coefficient_phase_difference_rad":float(cubic["m6_wrapped_phase_difference_rad"]),
        "final_mean_radius_difference_m":float(db["mean_jet_radius_m"][-1]-da["mean_jet_radius_m"][-1]),
        "final_energy_relative_difference":float(db["energy_m2_s2"][-1]/da["energy_m2_s2"][-1]-1),
        "max_runtime_m6_difference_at_matched_times_m":float(np.max(np.abs(db["m6_amplitude_m"][ib]-da["m6_amplitude_m"][ia]))),
        "events":{"attempts_each":len(ea),"injected_each":sum(e["injected"]=="true" for e in ea),
            "schedule_locations_decisions_identical":True,"max_strength_difference_s":dq,
            "max_strength_difference_over_configured_qmax":dq/qmax,
            "net_impulse_energy_m2_s2":[impulse_a,impulse_b],
            "energy_change_minus_impulses_m2_s2":[float(da["energy_m2_s2"][-1]-da["energy_m2_s2"][0]-impulse_a),
                                                    float(db["energy_m2_s2"][-1]-db["energy_m2_s2"][0]-impulse_b)]},
        "provenance":{label:{name:digest(run["path"]/name) for name in ("config.toml","Manifest.toml","diagnostics.csv","events.csv")}
                      for label,run in (("baseline",a),("refined",b))},
        "interpretation":"Very small short-time timestep sensitivity supports dt=100 s for the next exploratory pair; two stepsizes do not establish convergence order, long-time robustness, spatial convergence, or hexagon stabilization."}
    (out/"timestep_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    with (out/"matched_diagnostics.csv").open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["time_s","m6_dt_coarse_m","m6_dt_fine_m","difference_m"])
        writer.writerows(zip(matched_times,da["m6_amplitude_m"][ia],db["m6_amplitude_m"][ib],
                            db["m6_amplitude_m"][ib]-da["m6_amplitude_m"][ia]))
    plt.rcParams.update({"font.size":10,"font.family":"DejaVu Sans","pdf.fonttype":42,
                         "axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    fig.suptitle("Halving the timestep barely changes the short run",fontsize=16,fontweight="bold")
    blue,orange="#145E91","#BD4936"
    axs[0,0].plot(da["time_s"]/3600,da["m6_amplitude_m"]/1000,color=blue,lw=2,label="dt = %g s"%dt_a)
    axs[0,0].plot(db["time_s"]/3600,db["m6_amplitude_m"]/1000,"--",color=orange,label="dt = %g s"%dt_b)
    axs[0,0].set(xlabel="Simulated time (hours)",ylabel="Sixfold displacement (km)",title="A  The two amplitude histories overlap")
    axs[0,0].legend(frameon=False)
    axs[0,1].plot(matched_times/3600,1000*(db["m6_amplitude_m"][ib]-da["m6_amplitude_m"][ia]),"o-",color=blue,markersize=4)
    axs[0,1].set(xlabel="Simulated time (hours; exactly matched records)",ylabel="Amplitude difference (mm)",
                 title="B  Difference: refined minus baseline")
    keys=["psi","s","u","v","velocity_vector"]
    axs[1,0].bar(range(len(keys)),[fields[k]["relative_l2"] for k in keys],color=blue,width=.6)
    axs[1,0].set(yscale="log",xticks=range(len(keys)),xticklabels=[r"$\psi$",r"$s$",r"$u$",r"$v$","velocity"],
                 ylabel="Relative L2 difference",title="C  Final-field differences at 20,000 s")
    x,y=b["domain"]["x"],b["domain"]["y"]
    ld=b["config"]["physics"]["deformation_radius_m"]
    dz=(fields_b["s"]-fields_a["s"])+(fields_b["psi"]-fields_a["psi"])/ld**2
    mask=np.hypot(x[:,None],y[None,:])>24e6
    z=np.ma.masked_where(mask,dz*1e10)
    limit=float(np.max(np.abs(z)))
    im=axs[1,1].imshow(z.T,origin="lower",extent=[x[0]/1e6,x[-1]/1e6,y[0]/1e6,y[-1]/1e6],
                       cmap="RdBu_r",vmin=-limit,vmax=limit,interpolation="nearest")
    axs[1,1].set(xlim=(-22,22),ylim=(-22,22),xlabel="x (1,000 km)",ylabel="y (1,000 km)",
                 title="D  Final vorticity difference (expanded scale)")
    fig.colorbar(im,ax=axs[1,1],label=r"$\zeta_{50}-\zeta_{100}$ ($10^{-10}$ s$^{-1}$)",shrink=.85)
    for ax in axs.flat:
        ax.grid(alpha=.15)
    for ext in ("png","pdf"):
        fig.savefig(str(out/("timestep_comparison."+ext)),dpi=180)
    plt.close(fig)
    rows=[("Runs", "%s and %s"%(a["path"].parent.name,b["path"].parent.name)),
          ("Nominal timesteps / steps", "%g s / %d; %g s / %d"%(dt_a,da["step"][-1],dt_b,db["step"][-1])),
          ("Initial state and environment", "Identical initial fields, domain, source, Julia version and dependencies"),
          ("Final velocity relative L2 difference", "%.3e"%fields["velocity_vector"]["relative_l2"]),
          ("Largest local velocity difference", "%.3e m/s"%fields["velocity_vector"]["max_pointwise_difference_m_s"]),
          ("Final dynamic PV relative L2 difference", "%.3e"%fields["s"]["relative_l2"]),
          ("Final sixfold displacement difference", "%.4f mm (runtime diagnostic)"%(1000*summary["final_runtime_m6_difference_m"])),
          ("Independent cubic displacement comparison", "%.4f mm (1024 radial / 1440 angular samples)"%(1000*summary["final_cubic_m6_difference_m"])),
          ("Mode-six coefficient phase difference", "%.3e radians; pattern-angle difference is −1/6 of this"%phase),
          ("Forcing sequence", "Same %d attempts, %d injections; identical event times, locations and suppression decisions"%(len(ea),summary["events"]["injected_each"])),
          ("Largest injection-strength difference", "%.3e s⁻¹ (%.3e of configured PV cap)"%(dq,dq/qmax))]
    body="""<!doctype html><html><meta charset="utf-8"><title>Saturn QG timestep comparison</title>
<style>body{font:16px/1.55 system-ui,sans-serif;color:#203342;max-width:1100px;margin:40px auto;padding:0 25px}h1,h2{line-height:1.2}table{border-collapse:collapse;width:100%}td{border-bottom:1px solid #ddd;padding:10px}td:first-child{width:38%;font-weight:600}img{width:100%;margin:20px 0}.note{background:#edf4f8;padding:16px;border-left:4px solid #145e91}a{color:#145e91}</style>
<h1>Saturn QG: 100-second versus 50-second timestep</h1>
<p class="note"><b>Assessment:</b> timestep sensitivity is very small over the common 20,000-second window. Retaining 100 seconds is reasonable for the next exploratory matched forced/unforced pair. This establishes neither an integration convergence order nor long-time physical stability.</p><table>"""
    body+="".join("<tr><td>%s</td><td>%s</td></tr>"%(html.escape(k),html.escape(v)) for k,v in rows)+"</table>"
    encoded=base64.b64encode((out/"timestep_comparison.png").read_bytes()).decode("ascii")
    body+='<img alt="Matched-time timestep comparison" src="data:image/png;base64,'+encoded+'">'
    body+="""<h2>What the comparison establishes</h2>
<p>The configurations differ only in nominal timestep. Source and environment hashes are verified against the downloaded files, and initial fields are identical. Both runs end at 20,000 s. All three saved snapshot times match, while only 15 of 21 diagnostic times match exactly. The amplitude-difference panel uses only exactly matched diagnostic times; row-by-row subtraction of all records would introduce a sampling-time error.</p>
<p>The runs used different Unity GPU nodes (uri-gpu003 and uri-gpu002) under the same software environment. Their tiny differences can therefore include hardware-dependent floating-point roundoff as well as timestep sensitivity.</p>
<p>The amplitude histories above use the unchanged runtime bilinear jet extraction. A separate cubic remeasurement of the same snapshots is recorded in the adjacent jet-measurement audit. Diagnostic interpolation sensitivity and timestep sensitivity are distinct: the tiny differences between these runs do not validate weak competing modes.</p>
<p>Convection has the same random event schedule and locations, but its state-dependent strengths can differ slightly. The net impulse energy remains negative in both runs. The remaining energy change after subtracting impulses includes configured hyperdiffusion and numerical error; it is not a pure error estimate.</p>
<p>The next prepared check compares uninterrupted and checkpoint-resumed GPU trajectories with a split at an existing integration boundary. Longer forced and unforced runs can be scheduled to start only after that check succeeds. Both longer cases must start fresh from the same initialized flow, and still represent a single seeded realization over about one nominal jet circuit. Spatial resolution, forcing-radius, boundary and longer-time sensitivity remain open.</p>
<p>Reproduction: first run scripts/compare_jet_measurements.jl BASELINE_RESULTS REFINED_RESULTS OUTPUT using the pinned Julia environment, then scripts/compare_timesteps.py with the same arguments. Machine-readable values and input hashes are in timestep_summary.json; exact matched-time records are in matched_diagnostics.csv. The final-field differences are normalized by the refined run. Comparison with the actual change from the initial state is also recorded in the JSON. No model source or saved simulation data were modified.</p></html>"""
    (out/"report.html").write_text(body,encoding="utf-8")
    print(json.dumps({"final_field_comparison":fields,"runtime_m6_difference_m":summary["final_runtime_m6_difference_m"],
                      "matched_diagnostic_rows":len(matched_times),"events":summary["events"]},indent=2))


if __name__ == "__main__":
    main()

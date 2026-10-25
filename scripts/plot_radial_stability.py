#!/usr/bin/env python3
"""Summarize completed circular-base screens without modifying simulations."""
import argparse
import base64
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import numpy as np

from radial_stability import circular_base


def read_screen(path):
    if not (path / "COMPLETED").is_file():
        raise ValueError("incomplete radial screen: " + str(path))
    meta = json.loads((path / "metadata.json").read_text())
    rows = list(csv.DictReader((path / "leading_modes.csv").open()))
    groups = {}
    for row in rows:
        key = (float(row["deformation_radius_m"]), float(row["deep_ratio"]), int(row["radial_nodes"]))
        groups.setdefault(key, {})[int(row["m"])] = row
    return meta, groups


def rate(row):
    return float(row["leading_growth_rate_s"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fine", type=Path)
    ap.add_argument("coarse", type=Path)
    ap.add_argument("boundary", type=Path)
    ap.add_argument("--small-ld-refinement", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    fm, fine = read_screen(a.fine)
    cm, coarse = read_screen(a.coarse)
    bm, boundary = read_screen(a.boundary)
    assert fm["config_sha256"] == cm["config_sha256"] == bm["config_sha256"]
    assert fm["julia_source_sha256"] == cm["julia_source_sha256"] == bm["julia_source_sha256"]
    refinement = json.loads(a.small_ld_refinement.read_text())
    assert refinement["metadata"]["config_sha256"] == fm["config_sha256"]
    assert refinement["metadata"]["julia_source_sha256"] == fm["julia_source_sha256"]
    cfg = fm["config"]
    radius = cfg["initial"]["jet_radius_m"]
    model_rotation = cfg["planet"]["rotation_rate_s"]
    system_iii_period = 10*3600 + 39*60 + 22.4
    observed_period = 10*3600 + 39*60 + 23.01
    observed_c = 2*math.pi/observed_period - model_rotation
    # A growth tie does not select a physical neutral branch. Omit its speed.
    phase_display_floor = 1e-8
    a.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    comparison = []
    for (ld, ratio, nodes), modes in fine.items():
        best = max(modes, key=lambda m: rate(modes[m]))
        g = rate(modes[best])
        summaries.append(dict(ld_m="inf" if math.isinf(ld) else ld, deep_ratio=ratio,
                              fastest_m=best if g > 1e-12 else None, fastest_growth_s=g,
                              fastest_efolding_days=1/g/86400 if g > 1e-12 else None,
                              m6_growth_s=rate(modes[6]),
                              m6_angular_speed_s=(float(modes[6]["leading_angular_pattern_speed_s"])
                                                  if rate(modes[6]) > phase_display_floor else None),
                              largest_tested_m=max(modes), radial_nodes=nodes,
                              fastest_at_upper_mode_limit=g > 1e-12 and best == max(modes)))
        for label, other in (("radial_256_vs_512", coarse), ("outer_22_vs_24_Mm", boundary)):
            candidates = [(k, v) for k, v in other.items() if k[:2] == (ld, ratio)]
            if not candidates:
                continue
            key, old = max(candidates, key=lambda kv: kv[0][2])
            shared = sorted(set(modes) & set(old))
            changes = [abs(rate(modes[m])-rate(old[m])) for m in shared]
            comparison.append(dict(check=label, ld_m="inf" if math.isinf(ld) else ld,
                                   deep_ratio=ratio, old_nodes=key[2], new_nodes=nodes,
                                   shared_mode_max=max(shared), max_absolute_growth_change_s=max(changes),
                                   m6_growth_change_s=rate(modes[6])-rate(old[6]),
                                   m6_angular_speed_change_s=(float(modes[6]["leading_angular_pattern_speed_s"])-float(old[6]["leading_angular_pattern_speed_s"])) if min(rate(modes[6]), rate(old[6])) > phase_display_floor else None,
                                   m6_relative_growth_change=(rate(modes[6])/rate(old[6])-1) if rate(old[6])>1e-12 else None))
    result = dict(scope="Inviscid circular analytic pilot with an artificial disk wall; not the 2024 hexagon or a Saturn stability threshold.",
                  cases=summaries, comparisons=comparison,
                  growing_cases=sum(s["fastest_growth_s"]>1e-12 for s in summaries),
                  total_cases=len(summaries), sources=[str(a.fine), str(a.coarse), str(a.boundary)])
    result["small_ld_refinement"] = dict(source=str(a.small_ld_refinement),
                                           summary=refinement["summary"],
                                           radial_nodes=refinement["metadata"]["interior_radial_nodes"])
    result["phase_reference"] = dict(
        assumption="Positive mathematical theta is eastward, and configured planet.rotation_rate_s is the coordinate-frame rotation. Neither is registered to goflow observations yet.",
        model_rotation_s=model_rotation, model_rotation_period_s=2*math.pi/model_rotation,
        system_iii_period_s=system_iii_period, observed_2008_2014_period_s=observed_period,
        observed_angular_speed_in_model_frame_s=observed_c,
        observed_linear_speed_at_model_ring_m_s=radius*observed_c,
        phase_display_growth_floor_s=phase_display_floor,
        source="https://doi.org/10.1002/2013GL059078",
        limitations="Only the fastest-growing radial branch at m=6 is displayed; no neutral branch is selected and nonlinear pattern drift is not predicted. Weak growth below the display floor is unresolved, not certified stable.")
    (a.output/"summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False})
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.2), constrained_layout=True)
    colors = {-0.5:"#8e44ad", 0.:"#6a737d", .5:"#087e8b", 1.:"#d87420"}
    for ratio, color in colors.items():
        cases = sorted([s for s in summaries if s["deep_ratio"]==ratio and s["ld_m"]!="inf"], key=lambda s:s["ld_m"])
        label = "Deep/upper ring speed = {:+g}".format(ratio)
        for ax, field in ((axs[0,0], "m6_growth_s"), (axs[0,1], "fastest_growth_s")):
            ax.plot([s["ld_m"]/1e6 for s in cases], [s[field]*86400 for s in cases], "o-", color=color, label=label)
    bt = next(s for s in summaries if s["ld_m"]=="inf")
    for ax, field in ((axs[0,0], "m6_growth_s"), (axs[0,1], "fastest_growth_s")):
        ax.axhline(bt[field]*86400, color="black", ls="--", lw=1.1, label="Barotropic reference")
        ax.set(xlabel="Deformation radius (million metres)", ylabel="Growth rate (day⁻¹)", xscale="log")
        ax.xaxis.set_major_locator(FixedLocator([.5,1,2,4]))
        ax.xaxis.set_major_formatter(FixedFormatter(["0.5","1","2","4"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.2)
    axs[0,0].set_title("A  Sixfold disturbance of the circular jet")
    axs[0,1].set_title("B  Fastest growing disturbance among m = 1–32")
    axs[0,0].legend(fontsize=8, loc="best")
    n = max(k[2] for k in fine)
    for ratio, color in colors.items():
        modes = fine[(2e6, ratio, n)]
        axs[1,0].plot(sorted(modes), [rate(modes[m])*86400 for m in sorted(modes)], color=color, lw=1.8,
                      label="Deep ratio {:+g}".format(ratio))
    axs[1,0].axvline(6, color="black", ls=":", lw=1)
    axs[1,0].set(title="C  Mode competition at LD = 2 million metres", xlabel="Azimuthal mode m", ylabel="Growth rate (day⁻¹)", xlim=(1,32), ylim=(0,None))
    axs[1,0].grid(alpha=.2)
    r = np.linspace(7e6,23e6,800)
    for ratio, color in colors.items():
        base = circular_base(cfg,r,2e6,ratio*cfg["initial"]["jet_speed_m_s"])
        axs[1,1].plot(r/1e6, base["q_r"]*1e11, color=color, lw=1.6)
    axs[1,1].axhline(0,color="black",lw=.8)
    axs[1,1].set(title="D  Radial PV gradient at LD = 2 million metres", xlabel="Radius (million metres)", ylabel="PV gradient (10⁻¹¹ s⁻¹ m⁻¹)")
    axs[1,1].grid(alpha=.2)
    fig.suptitle("Lower-layer screening for the current circular pilot\nFixed upper wind; no convection; finite-amplitude hexagon stability remains untested", fontsize=14)
    fig.savefig(str(a.output/"radial_stability.png"), dpi=160)
    fig.savefig(str(a.output/"radial_stability.svg"))
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10.5, 5), constrained_layout=True)
    for ratio, color in colors.items():
        cases = sorted([s for s in summaries if s["deep_ratio"]==ratio
                        and s["ld_m"]!="inf" and s["ld_m"] >= 1e6], key=lambda s:s["ld_m"])
        ax.plot([s["ld_m"]/1e6 for s in cases],
                [radius*s["m6_angular_speed_s"] if s["m6_angular_speed_s"] is not None else np.nan for s in cases],
                "o-", color=color, label="Deep/upper ring speed = {:+g}".format(ratio))
    ax.axhline(radius*observed_c, color="#b3261e", ls="--", lw=1.5,
               label="2008–14 target if Ωref = 1.65×10⁻⁴ s⁻¹ and +θ is eastward")
    ax.axhline(0, color="black", lw=.7)
    ax.set(xlabel="Deformation radius (million metres)",
           ylabel="Angular pattern speed × 15 Mm (m/s; model frame)", xticks=[1,2,4])
    ax.set_title("Drift of the fastest-growing m = 6 radial branch\nCircular-base guide; no stable finite-amplitude hexagon has been identified")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8, loc="best")
    fig.savefig(str(a.output/"phase_speed.png"), dpi=160)
    fig.savefig(str(a.output/"phase_speed.svg"))
    plt.close(fig)
    table = []
    for s in summaries:
        ld = "∞" if s["ld_m"]=="inf" else "{:g}".format(s["ld_m"]/1e6)
        efold = "none detected" if s["fastest_efolding_days"] is None else "{:.3f}".format(s["fastest_efolding_days"])
        table.append("<tr><td>{}</td><td>{:+g}</td><td>{}</td><td>{}</td><td>{:.4f}</td></tr>".format(ld,s["deep_ratio"],s["fastest_m"] if s["fastest_m"] is not None else "—",efold,s["m6_growth_s"]*86400))
    img = base64.b64encode((a.output/"radial_stability.png").read_bytes()).decode()
    phase_img = base64.b64encode((a.output/"phase_speed.png").read_bytes()).decode()
    html = """<!doctype html><html lang="en"><meta charset="utf-8"><title>Lower-layer stability screen</title>
<style>body{font:17px/1.55 system-ui;max-width:1100px;margin:45px auto;padding:0 25px;color:#24313b}h1,h2{line-height:1.2}img{max-width:100%}table{border-collapse:collapse;width:100%;font-size:15px}td,th{border-bottom:1px solid #ddd;text-align:left;padding:7px}.note{background:#edf4f7;border-left:4px solid #087e8b;padding:14px}code{font-size:.9em}</style>
<h1>Can the lower layer support a hexagon with the observed drift?</h1>
<p class="note">This is an initial screen of the existing circular pilot. It does not yet test the finite-amplitude Pizzo–Salmon hexagon, and it does not establish a Saturn stability boundary.</p>
<p>All comparisons retain the upper 100 m/s annular jet, 15,000-km nominal radius, 1,000-km width parameter, and 50 m/s central cyclone. Only deformation radius and the prescribed annular deep speed change. Deep ratio 1 matches the annular component, not the cyclone. The lower layer is prescribed and time independent. Convection and dissipation are absent.</p>
<h2>Screening results</h2>"""
    html += '<p>Growing discrete modes occur in {} of {} sampled cases. Absence of a growing mode at finite resolution is not a stability certificate.</p>'.format(result["growing_cases"],result["total_cases"])
    html += '<img alt="Circular-base growth and PV-gradient screen" src="data:image/png;base64,'+img+'">'
    html += '<table><tr><th>LD (million m)</th><th>Deep/upper ring speed</th><th>Fastest m</th><th>Fastest e-fold time (days)</th><th>m=6 growth (day⁻¹)</th></tr>'+''.join(table)+'</table>'
    html += '<h2>Stability and phase speed must be satisfied together</h2>'
    html += '<p>A stable circular jet is not a finite hexagon, and a persistent hexagon at the wrong drift is not an observational match. The target is a sixfold relative equilibrium with the observed shape and upper winds, the correct angular drift in a common frame, and robustness to disturbances.</p>'
    html += '<img alt="Growing mode-six branch drift and conditional observational reference" src="data:image/png;base64,'+phase_img+'">'
    html += '<p>The configured rotation rate is {:.6g} s⁻¹ (period {:.3f} s), whereas System III has period {:.1f} s. Assuming positive model θ is eastward and this rotation defines its frame, the <a href="https://doi.org/10.1002/2013GL059078">2008–2014 observed period {:.2f} s</a> corresponds to c={:.6g} s⁻¹, or {:.3f} m/s at the model ring. This is a conditional coordinate conversion, not an observation at the pilot’s exact radius. The longitude convention, projection and observation epoch still need registration to goflow. Both measured winds and pattern drift must be transformed together.</p>'.format(model_rotation, 2*math.pi/model_rotation, system_iii_period, observed_period, observed_c, radius*observed_c)
    html += '<p>The plotted c is Re(c) of the fastest-growing radial eigenmode at m=6, not the drift of an equilibrated nonlinear hexagon. Lines may switch radial branches and are not solution continuation. Neutral ties and weak rates below 10⁻⁸ s⁻¹ are omitted from this speed plot; omitting a point does not certify stability. A reference-frame change shifts model and observational velocities together and cannot eliminate their physical mismatch.</p>'
    html += """<h2>Interpretation</h2><p>For this analytic circular profile, adding a deep annular flow in the same direction does not universally suppress growing modes. A smaller deformation radius changes both the PV gradient and the perturbation inversion. A quiet m=6 alone would not establish stability if other azimuthal modes grow. These results concern the chosen profile and are not evidence against lower-layer stabilization for other equilibria.</p>
<p>The circular-base equation is <code>(U/r−c)[∂rr+r⁻¹∂r−m²/r²−LD⁻²]φ−(Q_r/r)φ=0</code>, with <code>Q_r=U″+U′/r−U/r²−γr+LD⁻²(Ud−U)</code>. Disturbances vary as exp[im(θ−ct)], so growth is m Im(c). Deep flow acts through PV, not through a direct damping term.</p>
<h2>Numerical scope</h2><p>The fine screen uses 512 radial interior points and m=1–32, with a circular impermeable boundary at 24 million metres. Lower-resolution screens use 128 and 256 points over m=1–16. A separate 22-million-metre boundary check uses 470 points over m=1–24 for LD=2 million metres and deep ratios 0, 0.5 and 1, plus the barotropic reference. Full per-mode comparisons and raw spectra are saved. Nine independent analytic and convergence tests passed, including the disk Bessel-wave dispersion relation.</p>
<p>The 0.5-million-metre cases were additionally refined to 768 radial points using an equivalent matrix formulation checked against the primary generalized solver. At this resolution the maximum growth rates for deep ratios −0.5, 0, 0.5 and 1 are 0, 1.67×10⁻¹⁰, 3.44×10⁻⁸ and 1.88×10⁻⁷ s⁻¹. The weak positive growth rates and leading branches change substantially with resolution. They do not define a converged stability boundary; the very small resting-deep rate is particularly inconclusive. The counterflow case stays numerically quiet across the tested modes, but this is still a circular-base result. All refined leading eigenpairs satisfy the original matrix equation to relative residual below 1.14×10⁻¹³.</p>
<p>The radial boundary differs from the nonlinear Julia model's periodic square and continued planetary PV. Therefore even converged radial growth rates need comparison with weak-perturbation two-dimensional runs. Neither disk eigenvalues nor an axisymmetric mean profile determine the stability of a nonlinear hexagonal equilibrium.</p>
<h2>Using the 2024 solutions</h2><p><a href="https://doi.org/10.1017/jfm.2024.932">Pizzo &amp; Salmon (2024)</a> provides the nonlinear reference family. Published Figure 7 has n=5 and ε=0.03; the polygon is sixfold. First verify the full Eulerian initial state and ambient patch. A steady barotropic equilibrium ψ0 remains an exact finite-LD equilibrium over resting deep flow because qF=qBT−LD⁻²ψ0 and J(ψ0,ψ0)=0. Stability may still change. Adding a general circular deep jet usually requires recomputing the equilibrium before testing stability.</p>
<p>The independent map audit records a phase/PV convention mismatch in the printed formulas. It tests both signs using an independent Cartesian curl and does not silently modify the model. The original initializer has not yet been located. The downloaded SaturnHexagonIC.m instead builds a Gaussian jet around a Ptolemaic-shaped boundary.</p>
<p>The next physics study should continue verified finite-amplitude hexagonal equilibria across LD and circular deep profiles while solving for pattern speed, then test stability. Seek the overlap of acceptable drift, acceptable winds/shape and stability. Fixing the pattern speed to the observed value while allowing appropriate deep-profile parameters to vary is an alternative constrained continuation. Convection and spontaneous mode selection then become separate tests. The current million-second pair changed convection alone and cannot answer the lower-layer question.</p>
<p>Reproducibility: summary.json records all cases and convergence comparisons. Source directories retain the input configuration, code hashes, complete eigenvalue spectra, leading radial eigenfunctions and completion markers. Production Julia dynamics and dependencies are unchanged.</p></html>"""
    (a.output/"report.html").write_text(html, encoding="utf-8")
    print(json.dumps({"growing_cases":result["growing_cases"],"total_cases":result["total_cases"]}))


if __name__ == "__main__":
    main()

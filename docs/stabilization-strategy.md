# Stabilization and phase-speed strategy

**This is a research proposal; it neither implements a new experiment suite nor
establishes a stable hexagon.** It addresses two complementary questions: can a prescribed
deep layer suppress destructive instability of a circular jet, and can it move
an already persistent hexagon toward the observed drift while retaining its
shape and upper winds?

The two completed million-second runs show rapid amplification and roll-up of
the current seeded annular profile, with and without convection. They neither
isolate the deep-layer effect nor rule out other equilibria. A separate
preliminary circular linear check with complete deep-flow matching is now
finished. At 256 radial points, sampled fastest modes remain strongly unstable
for `LD = 1000, 2000, 4000 km` (respectively m12, m13, m14; e-folding times about
1.47, 0.86, 0.76 days). At 500 km, refinement from 128 to 256 to 512 radial
points changes the apparent fastest growth from zero to 1.00e-7 to 2.61e-7 s^-1,
and the leading mode changes. At 512 points, m6 also grows (2.28e-7 s^-1).
A 22-Mm outer disk at comparable spacing gives a maximum of 2.84e-7 s^-1.
Estimated corotation layers remain smaller than a radial cell: this is a
weak-growth, unresolved regime, not an established stable candidate. The
calculations cover m1-24 and are saved locally as `full-deep-match-01` and
`full-deep-match-refine-01`. These are linear results, not nonlinear hexagon tests.

The archived pair supports the observation that convection made little
difference here. Initial fields are identical and the only physics change is
convection enabled/disabled. Forcing was active: 482 nonzero injections from
673 attempts. Provisional runtime m6-amplitude differences are about 0.21% at
300,000 s and 1.21% at 500,000 s, using interpolated forced diagnostics.
Measured convective energy jumps sum to +0.264306 and −1.553842 m²/s² for
positive and negative changes: summed absolute change 1.818148 (0.159%)
and net −1.289536 (−0.113%), relative to initial total QG energy 1140.095763 m²/s².
This is a small forcing budget alongside rapid jet instability; it does not
prove that convection never matters. Keep the next benchmarks unforced and
calibrate forcing separately before interpreting a maintenance effect.

## Priority 1: isolate stabilization of the circular flow

[Thomson & McIntyre (2016)](https://doi.org/10.1175/JAS-D-14-0370.1) motivate
testing finite deformation together with deep-flow structure: both disturbance
inversion and the background PV gradient matter. Their straight Jovian jets do
not establish the existence or drift of a polar hexagon.

Keep the upper axisymmetric streamfunction `psi0(r)` fixed, including the
annular jet **and central cyclone**. Add a prescribed circular control

\[
U_d(r)=\alpha U_{\rm total}(r),\qquad
\psi_d(r)=\alpha\psi_0(r),\qquad F=L_D^{-2}.
\]

Use the same reference constant, outer continuation and spectral projection.
The current deep initializer contains only an annular jet; it cannot implement
this complete match by changing its existing speed parameter. A production
deep-profile extension is required; the preliminary check above uses a separate
analysis implementation. The time-evolution solver can be retained.

Start with `alpha = {0, 0.5, 1}` and `LD = {2000, 3200, 4000} km`, plus one
barotropic control. These are mechanism brackets, not fitted Saturn parameters.
For complete matching, also investigate the 500-1000 km transition, with
refinement of the weak modes and explicit disturbed-layer thickness checks.
The old small-LD thickness failure from an unmatched central cyclone does not
apply to the exactly matched base; small LD still needs a separate physical
assessment for a finite meander and for partially matched deep profiles.
Retain the existing annulus-only cases as a separate comparison. Removing the
cyclone is another controlled experiment, not part of the fixed-upper-flow sweep.
Recompute PV from the fixed streamfunction for each case; holding PV fixed
during Helmholtz inversion would change the upper winds.
The illustrative 1-Mm width favors m≈12–15 at larger `LD`: vary `LD/w` and
observationally admissible jet width/cyclone shape, since broadening can change
mode choice and reduce Rossby number, while retaining goflow wind constraints
instead of tuning deep amplitude alone.

In the untapered interior,

\[
Q_r=U''+U'/r-U/r^2-\gamma r+F(U_d-U).
\]

At `alpha=1`, stretching drops out of this **base** PV gradient while finite
`LD` remains in perturbation inversion. Matching `psi_d=psi0` also removes
base interface displacement under the shared gauge. Neither cancellation
proves stability, fixes pattern speed, or reduces the upper flow's Rossby
number; disturbed states can still have substantial interface displacement.

First compute the unforced spectrum across all resolved azimuthal modes,
initially `m=1–32`, extending the range until the high-mode result is resolved.
Track eigenfunctions as parameters change; refine radial resolution, outer
radius, and then the two-dimensional mesh, timestep, boundary continuation
and dissipation. Confirm selected linear rates with small disturbances in the
evolution model. A circular disk calculation alone cannot certify the periodic
square implementation.

Search near a resolved **m6 neutral transition**, while checking competing modes,
for a nonlinear sixfold branch. Follow that branch with a steady/rotating
equilibrium calculation and independently perturbed integrations. A stable
circle, a quiet m6 alongside growing competitors, and a persistent finite
hexagon are different outcomes.

## Priority 2: continue a verified Pizzo–Salmon wave

[Pizzo & Salmon (2024)](https://doi.org/10.1017/jfm.2024.932) provide the
barotropic Ptolemaic family and examples of finite-time persistence. Reconstruct
a complete sixfold member (`n=5`), including its constant-PV ambient region and
far-field flow. Resolve the map inversion, dimensional scales, gamma convention,
hemisphere, velocity/PV signs and reference frame. The
[legacy initializer audit](legacy-ic-audit.md) verifies a local formula, not this
complete equilibrium. Particle circulation frequency must not be substituted
for Eulerian pattern drift.

Verify the reconstructed velocity, PV inversion and rotating-equilibrium
residual before adding finite deformation. For a circular prescribed deep
streamfunction, the target is

\[
J\!\left(\psi-c r^2/2,
        (\Delta-F)\psi+f+F\psi_d\right)=0.
\]

Solve for **shape and angular pattern speed `c` jointly** as `F` and the circular
deep profile change. Specify the branch's PV/invariant constraints, an amplitude
constraint and a phase condition so continuation selects a definite family.
Measure residual convergence in the physical disk and sensitivity to the
ambient/boundary treatment; the periodic square is not globally rotation-invariant.

A useful analytic check is that adding finite `F` to a verified barotropic
relative equilibrium gives the residual
`F (c - Ud/r) partial_theta(psi)` if its shape and speed are held unchanged.
Thus resting deep flow preserves a stationary member, and `Ud=c*r` preserves a
rotating member in the ideal interior. A generic annular deep jet requires
continuation. Do not prescribe a deep hexagon to hold the upper one in place:
the proposed physical tests keep the deep layer circular.

Then measure stability of the **finite-amplitude equilibrium**, using its
two-dimensional linearization or controlled independent perturbations. The
circular spectrum does not supply this result.

## Priority 3: reproduce the persistent shallow-water precedent

[Rostami, Zeitlin & Spiga (2017)](https://doi.org/10.1016/j.icarus.2017.06.006)
provide a closely related nonlinear precedent: a rotating-shallow-water
jet–cyclone system supports a long-lived hexagonal structure, with a reference
deformation radius near 3200 km; its drift remains an observational limitation.
First reproduce their specified jet, cyclone, balanced thickness, perturbation,
domain and numerical conditions in rotating shallow water, including the
jet-only comparison. Obtain the persistent state and drift before modifying it.

Construct a **balanced QG counterpart**, not an exact transfer of the shallow-water
state. Compare winds, PV structure and geometry while quantifying the effects of
discarded divergence, gravity waves and finite-Rossby-number balance. Assess
`|zeta|/f0`, shear/curvature Rossby measures, thickness variations and balance
residuals, especially in the central cyclone. If that cyclone requires dynamics
outside QG, retain shallow water for that route rather than interpreting QG
failure as failure of the hexagon mechanism. A credible QG counterpart can then
enter the same circular-deep-flow continuation and stability tests as route 2.

## Shared observation test; convection last

In every route constrain jet radius, width, upper winds and
cyclone against the eventual goflow estimates and their uncertainties. Use
one registered epoch, longitude sign and rotation frame for model and observed
`c`. Changing coordinates does not correct a physical drift mismatch; changing
planetary rotation changes the model dynamics.

Advance only candidates with converged residuals, credible QG ordering where
applicable, persistent connected jet geometry, and acceptable drift and winds.
Calibrate convection afterward, then compare matched unforced, unbiased and
biased cases, recording realized injections and budgets. It should test maintenance or
selection of an understood state, not conceal an initialization imbalance.

# Local checks for the 2017 QG benchmark

These checks concern the regularized QG initial conditions and software. They
do not establish a persistent nonlinear hexagon or reproduce the paper's RSW
result. The frozen experiment set is `configs/rostami_01` (32 configurations).

The complete Julia test suite passed **3,418 assertions**. Python discovery
passed **109 tests**, covering radial stability, campaign manifests, mocked
Slurm submission/collection safeguards, diagnostics and gallery handling.
All ten campaign stages passed their source/configuration provenance checks.
The dependencies remain the existing pinned Julia 1.10.12 environment.
An independent review also checked the radial derivatives and polar limits
against high-precision analytic differentiation (maximum relative discrepancy
4.3e-15 in the points checked).

## Actual-grid initialization checks

`scripts/check_rostami_initialization.jl` constructs the production CPU state,
measures `J(psi,q)` from the projected fields, compares velocity to the analytic
profile, and advances three steps through the production timestep guards. The
audit removes the nonaxisymmetric seed explicitly. These are short checks,
not substitutes for the Unity axisymmetric controls or perturbed integrations.

The residual normalization is `Tjet * RMS(J) / Qjet`, where
`Tjet=2*pi*Rjet/Ujet=827122.44 s` and `Qjet=Ujet/wjet=4.43405e-5 s^-1`.
It is an advective-form residual of the sampled fields, not a norm of a
long-time integration error. Three-step evolution retains hyperviscosity.

| Circular state | Grid | Core smoothing | Disk residual | Initial max abs(zeta/f) |
|---|---:|---:|---:|---:|
| Jet only | 1024² | 500 km | 4.92e-9 | 0.138 |
| Reference core | 1024² | 500 km | 9.79e-5 | 1.702 |
| Stronger core | 1024² | 500 km | 1.74e-4 | 2.269 |
| Reference core | 2048² | 500 km | 2.24e-9 | 1.702 |
| Reference core | 2048² | 250 km | 1.17e-4 | 2.568 |
| 10% reference core | 1024² | 500 km | 1.09e-6 | 0.170 |

The largest three-step relative streamfunction change was 5.63e-12. At 1024²
the reference core's RMS velocity projection error is 1.70e-6 m/s in the core
and 9.10e-8 m/s in the jet annulus. Its jet-annulus residual is 1.23e-5.
Increasing the resolution at fixed 500-km smoothing substantially reduces the
residual; decreasing smoothing at fixed grid makes the core harder to resolve.
The prescribed paired core/grid tests are therefore part of the experiment,
not optional evidence for a later convergence claim.

The reference state's saved-gauge maximum abs(hprime/H) is 0.293. Full matching
can remove base interface displacement, but cannot reduce its order-one
upper-layer Rossby number. The reference/stronger-core runs are mechanistic
QG analogues. The weaker-core controls probe a more defensible QG regime.

## Preliminary linear screen

The frozen three circular controls were screened for modes 1–16 at 128 and
256 radial nodes, with the circular boundary at 28 Mm. All three selected
mode 6 as the fastest growing mode **within this tested band** on both grids.
At 256 nodes the reference core gave growth 5.182e-6 s^-1 (e-folding time
about 2.23 Earth days) and angular drift 2.686e-6 rad/s in the model frame.

A follow-up for the reference core at 512 radial nodes extended the band to
modes 1–32. Mode 6 remained fastest, with growth 5.1843e-6 s^-1 and angular
drift 2.6870e-6 rad/s. The leading eigenpair residual was 3.72e-11; these
growth/drift values differ from the 256-node values by less than 0.06%.
The outer-boundary and nonlinear checks remain outstanding.

This is a circular-base, inviscid linear result. It says nothing yet about
nonlinear saturation or the stability of the resulting polygon. The coarse
core resolution, outer boundary, unexamined higher modes and observation-frame
registration require separate checks. No observed phase-speed agreement is
claimed.

## Unity validation

CPU/GPU agreement for this release must be checked on Unity using the
`validation` stage. The old campaign's successful GPU checks do not validate
this new initializer. GPU receipts and reports, once produced, record the exact
deployed source and configuration hashes. See [the Unity runbook](unity-rostami.md).

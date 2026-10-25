# A testable Saturn hexagon programme

The initial scientific question is whether finite deformation radius, prescribed deep circulation, and parameterized convective PV injection increase the persistence of a coherent hexagonal circumpolar jet. A stronger, separate question is whether the same dynamics preferentially produce sixfold structure from states without a seeded hexagon. These are experimental goals; the initial implementation and example configurations do not establish either result.

**Current priority (24 October 2026): lower-layer stabilization of the vortex,
using Pizzo & Salmon's 2024 solutions as reference states, with observed
pattern drift and upper winds as joint constraints.** The completed
million-second convection pair held the deep flow and deformation radius fixed.
Both runs amplified the imposed meander and developed rolled-up lobes, so that
pair cannot establish whether the lower layer stabilizes the vortex. The next
physical comparison must change the lower-layer coupling. Convection and
emergence follow after identifying and checking suitable equilibria. See the
[maintenance-pair analysis](../gallery/maintenance-01/index.html) and
the specific construction and stability tests below.

Thomson and McIntyre's result motivates the controls: finite deformation radius and deep jets can suppress long-wave meandering and support comparatively straight jets. For Saturn, suppressing all meanders would erase the target hexagon. Success therefore requires a persistent, finite-amplitude meander with the appropriate drift and jet structure, not merely a quieter or stronger zonal flow. Their beta-plane forcing and stability conclusions must be assessed again in polar geometry. [Thomson & McIntyre (2016)](https://doi.org/10.1175/JAS-D-14-0370.1)

## 1. Establish the numerical and physical reference

Begin with smooth axisymmetric states, weak disturbances, and no forcing. Verify the inversion sign, polar velocity convention, stationary axisymmetric solution, barotropic limit, and energy/PV budgets with negligible dissipation. Repeat with finite deformation radius and a prescribed deep profile. Check that deep-flow settings have no effect in the **unforced** `Ld = Inf` limit: the factor multiplying the prescribed deep streamfunction in PV then vanishes exactly. Enabled convection can still depend on deep streamfunction through its coldness proxy.

At each resolution record the **physical** box size, jet width, deformation radius, forcing radius, injection rate per unit time, dissipation coefficient, and resolved spectral cutoff. Increasing resolution must not shrink the injected vortices or move the boundary. A forcing radius expressed in grid cells tests a different physical model each time.

The outer polar PV continuation is part of the experiment. Its radius, transition width, and distance from the jet require convergence tests at fixed interior physics. A periodic Cartesian square also has fourfold grid symmetry. Rotate perturbations relative to the grid and compare resolutions and domains, especially when evaluating mode selection. A trap or continuation should not set the selected wavelength or pattern speed.

## 2. Three initial-condition families

| Family | Initial condition | What the result can establish |
|---|---|---|
| Maintenance | Jet plus a small or equilibrated `m = 6` displacement, with and without the central cyclone | Survival, decay, breakup, and the basin of attraction of an existing hexagon |
| Selection | Same axisymmetric jet, equal treatment of a broad range of azimuthal modes, no preferred phase or privileged `m = 6` amplitude | Which disturbances grow and whether nonlinear equilibration preferentially produces six |
| Emergence | Weak random flow or rest under convective forcing, initially without an upper jet or imposed sixfold component | Whether a jet and polygon organize from less structured initial states |

For the emergence family, distinguish runs with a prescribed deep annular jet from runs without one. A deep jet already supplies a preferred radius and mean structure. Emergence above a prescribed deep jet is a useful conditional result, but is a different claim from forming the whole polar jet system without that template. Fully unforced random-flow decay is another useful comparison.

Siegelman et al. show polar vortex crystals emerging from freely decaying, initially random barotropic turbulence on a gamma plane. That supplies a competing outcome and a geometry reference, rather than a prediction that convection will make a Saturn-like jet. Their trapping boundary must also be distinguished from whichever smooth periodic continuation is used here. [Siegelman, Young & Ingersoll (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170063/)

## 3. Paired controls before a broad parameter sweep

Use the same upper-layer velocity profile, physical perturbation spectrum, and random seed when changing one control. Recompute PV from that velocity consistently for each deformation radius; holding PV fixed instead would change the initial wind. Finite deformation radius changes potential energy even for a common velocity, so report both the initial kinetic and total QG energy.

| Control | Required comparisons | Interpretation |
|---|---|---|
| Deformation | Barotropic limit versus several finite `Ld / jet_width` values | Effect of PV stretching and altered inversion |
| Deep flow | Absent versus prescribed radial profile at finite `Ld` | Whether coupling to deep structure changes stability and drift |
| Convection | Unforced; unbiased injection; biased injection | Distinguishes stretching/deep-flow effects from forcing effects |
| Forcing dependence | Fixed-amplitude bias versus interface-dependent bias/amplitude, where implemented | Separates a sign preference from state-dependent feedback |
| Polar cyclone | Absent versus present, with a small range of strengths | Sensitivity to the central circulation |
| Initial disturbance | Seeded six, broadband, and no pre-existing jet | Maintenance versus selection versus emergence |

The unforced barotropic deep-flow comparison is a numerical identity test, not an additional independent physical factor. With the current convection prescription, deep streamfunction can affect the coldness proxy even in that limit; label such differences as forcing-only effects and do not attribute them to dynamical deep-layer coupling. A resting lower layer still has stretching when `Ld` is finite. A spatially constant streamfunction offset is also relevant to interface-dependent forcing and needs an explicit gauge/reference convention. The current zero-domain-mean convention and its domain-size sensitivity are described in [model.md](model.md).

Start with a modest paired set, then add seeds to quantify variability near transitions. Do not treat different snapshots of one chaotic run as independent realizations. Keep a planned stopping time and analysis window; extend clearly transitional runs rather than selecting only visually attractive times. The first examples are pilots, not a completed factorial campaign.

For biased versus unbiased forcing, equal event rate and peak PV do not necessarily give equal energy input. Record every realized impulse's energy and enstrophy increments and its spatial mean. Report nominal forcing parameters alongside cumulative realized injection. If using an equal-energy-input comparison, explicitly normalize it and retain a separate equal-event-statistics comparison. Mean PV compensation and any quasi-frictional effects belong in this accounting.

The first convection implementation should be described by its actual equations and simplifications. Do not label every stochastic vortex-pair forcing as a full reproduction of Thomson–McIntyre: state dependence, clipping, bias, event timing, and mean compensation each matter.

## 4. Classify the evolving structure quantitatively

The leading observable is the jet centreline `r_j(theta,t)`, measured from tangential velocity maxima within a specified annulus. Extract its displacement spectrum and phase. For

`C_m(t) = mean_theta[(r_j(theta,t) - mean_theta(r_j)) exp(-i m theta)]`,

the mathematical angular pattern speed is `Omega_m = -(1/m) d arg(C_m)/dt`. Unwrap phase only when amplitude is appreciable, and document the mapping from mathematical theta to the observed longitude convention. The System III pattern drift and the wind speed through the wave are distinct observables.

Track, or add as analysis capabilities before making each associated claim:

- Jet radius, peak speed, width, closed-contour continuity, and its displacement amplitudes over several `m`; the dominance of six should not depend on looking at only `m = 6`.
- Sixfold phase coherence and pattern drift over a prespecified window, together with amplitude and lifetime in jet-circulation periods.
- Mean/eddy kinetic energy and QG total energy, forcing increments, dissipation and residuals; radial Reynolds stress and PV flux where available.
- Central-cyclone strength and off-axis vortex positions, circulation and lifetime; vortex tracking is needed to separate a coherent jet from six isolated vortices.
- The full azimuthal-mean velocity/PV profiles and their evolution, not only a final image.

A ring of six vortices can have a strong `m = 6` Fourier signal. It is not sufficient evidence of a hexagonal jet. Conversely, a straightened circular jet is a possible stabilized state but fails the finite-amplitude hexagon target. Predefine classification thresholds from measurement uncertainty and numerical sensitivity instead of choosing them after seeing the simulation outcomes. The diagnostics above are an analysis specification; consult the actual output schema for what is currently implemented.

## 5. Use the upcoming goflow velocities as observations

Preserve observation time, longitude convention, latitude definition (planetographic versus planetocentric), projection, velocity units, spatial resolution, masks, and velocity uncertainties. Fit the mean jet and polar cyclone with uncertainty ranges rather than one hand-picked profile. Compare model velocities after applying the observation's smoothing, masks, sampling and projection; include pattern speed, displacement spectra, mean flow, and eddy statistics.

Differentiate observed velocities only after accounting for noise amplification. Relative vorticity may be estimated from wind gradients; the complete QG PV also depends on deformation radius, prescribed deep streamfunction, and the interface/streamfunction reference. Treat inferred PV as model-dependent. Do not present cloud-top velocities as direct measurements of the deep circulation or convective forcing.

Reserve some epochs or regions for evaluation rather than fitting every observable. Carry plausible wind uncertainties through initial conditions and compare predictive ranges. This can determine which mechanisms are distinguishable with the new measurements, without fitting one striking snapshot.

## 6. Decision points

1. **Numerical credibility:** pass conservation, restart and CPU/GPU checks; show resolution and boundary convergence for the short-time instability and drift.
2. **Maintenance:** find paired cases where finite `Ld`, deep flow or convection improves jet coherence while retaining nonzero sixfold amplitude; report cases that straighten or destroy it.
3. **Selection:** test whether six is robust across broadband seeds and neighboring parameters, including competition with other wavenumbers.
4. **Emergence:** attempt long forced runs without an initialized upper jet, after the forcing and energy balance are independently verified.
5. **Observational and spherical comparison:** confront the surviving regimes with goflow uncertainties and your collaborator's spherical model. Slow drift can be sensitive to geometric corrections even when the polar morphology looks similar.

The next model increment should be chosen by these results. The separate
radial stability screen described below makes circular-base mode selection
and phase-speed interpretation cheaper. A gamma-plane model is a focused
mechanism experiment, while comparison on the sphere tests how much that
mechanism depends on the polar approximation.

## 7. Lower-layer stability and the 2024 exact-wave family

The relevant source is [Pizzo & Salmon (2024), *Exact planetary waves and jet
streams*, JFM 998, A43](https://doi.org/10.1017/jfm.2024.932), especially sections
3.2, 3.3 and 4. In the positive-n family the polygon has `m = n+1` sides:
the hexagonal example is `n = 5`, not `n = 6`. The published Figure 7 uses
epsilon = 0.03, A0 = 0.2, gamma = 10; the preceding text specifies omega =
0.004. The local June draft and arXiv version use epsilon = 0.05 for this
example. These dimensionless parameters are not yet Saturn-calibrated inputs.
The smooth-map condition `|n epsilon| < 1` prevents cusps; it is not a
stability criterion. The paper demonstrates persistence over its integration,
not a complete stability boundary or spontaneous selection of six.

### Preserve the equilibrium before testing stability

Write F = LD^-2 and let psi0 be a **verified complete steady barotropic
solution**, including its ambient patch. With the code's velocity convention,

    qBT = Laplacian(psi0) + f,
    J(psi0, qBT) = 0.

Over a resting prescribed deep layer, the same upper streamfunction is an
exact inviscid 1.5-layer equilibrium for every F:

    qF = qBT - F psi0,
    J(psi0, qF) = J(psi0, qBT) = 0.

This is a useful controlled extension of the 2024 solution: equilibrium
existence is preserved while perturbation inversion and stability change.
It is **not** a proof that finite LD stabilizes it. Hold psi0/upper velocity
fixed and recompute PV for each LD; holding barotropic PV fixed in the
Helmholtz inversion would change the upper flow.

For a barotropic relative equilibrium rotating at angular speed c, define
Phi = psi0 - c r^2/2. A fixed circular deep flow Ud(r) introduces the exact
inviscid residual

    partial_t qF + J(psi0, qF)
        = F [c - Ud(r)/r] partial_theta psi0.

Thus a solid-body deep flow Ud = c r preserves the original rotating
equilibrium, and resting deep flow preserves a stationary one. A general
annular deep jet changes the equilibrium as well as its stability. Continue
the solution to that background and measure the steady/rotating residual
before attributing its subsequent evolution to instability. A stationary
patterned deep streamfunction H(psi0) also preserves a stationary equilibrium,
but it supplies the polygon's geometry below; use that only as an explicit
mechanism control, not evidence of sixfold emergence above a circular deep flow.

The paper's Lagrangian phase advances particles along a fixed-b contour.
That particle circulation frequency is not automatically the Eulerian
pattern speed. Infer pattern drift from the actual full field in a stated
planetary reference frame.

### Joint requirement: a stable hexagon with the observed phase speed

Stability alone is insufficient. A circular state can be stable without
supporting a finite hexagon; a persistent hexagon can have the wrong drift.
For circular prescribed deep flow, continue relative equilibria satisfying

    J(psi - c*r^2/2, (Laplacian-F)*psi + f + F*psi_d) = 0.

Here c is the angular pattern speed, not a material circulation frequency.
The finite-amplitude shape, c and equilibrium residual must be determined
together. Continue in LD and deep-profile parameters with a phase condition
and explicit normalization (for example, sixfold amplitude), then compute
perturbation growth in the moving frame. Alternatively fix c to the observed
value and solve for a suitable free deep-profile parameter. A changed
parameter must still satisfy the observed upper wind constraints. Seek the
intersection of persistent sixfold shape, appropriate drift and stability;
the circular-mode screen only helps locate candidate branches.

The production pilot has `planet.rotation_rate_s = 0.000165`, corresponding
to 10 h 34 m 39.911 s. This is **not System III**, whose period is
10 h 39 m 22.4 s. If this configured rate defines the model frame and positive
mathematical theta is eastward, angular speeds transform as

    c_SystemIII = c_model + Omega_model - Omega_SystemIII,
    c_observed_model = 2*pi/T_hexagon - Omega_model.

The [2008–2014 measurement of Sánchez-Lavega et al. (2014)](https://doi.org/10.1002/2013GL059078)
gives T_hexagon = 10 h 39 m 23.01 s. It corresponds to
`c_observed_model = -1.21761412414e-6 s^-1`, or **-18.2642 m/s at the model's
15-Mm radius**, under the stated frame assumption. This is a conversion of
an observed angular rate, not a new measurement at that radius. The paper's
positive west-longitude drift +0.0129 deg/day must have its sign reversed
when expressed as east-positive angular drift. Record the observation epoch,
longitude sign, radius/latitude projection and uncertainties explicitly;
do not treat a single historical period as a timeless exact target.

A coordinate change by delta Omega shifts both c and U/r by -delta Omega,
leaving the intrinsic speed c-U/r unchanged. Transform the measured winds,
deep velocities and pattern speed consistently. Simply replacing the
configured planetary rotation changes the gamma-plane PV gradient and is
not a postprocessing frame change. Existing runs stay in their original
frame. The reported [2021 alternative-frame speed near -22 m/s](https://doi.org/10.1029/2021GL092461)
illustrates this convention dependence; it is not an additional independent
measurement to fit.

For each continued finite hexagon, retain angular c, model-frame and
observation-frame drift, uncertainty, amplitude and shape, maximum
perturbation growth, and fit to the upper velocity field. Track the same
branch by eigenfunction/state overlap. A neutral spectrum has many real
eigenvalues: taking the first maximum-imaginary eigenvalue does not select
the relevant sixfold wave. The screening plot therefore suppresses neutral
ties and unresolved weak-growth speeds rather than assigning them a
hexagon drift. Neither the fastest circular m=6 branch nor a transient
jet-centreline phase measurement is automatically the nonlinear pattern
speed.

### Circular precursor and screening calculation

For an axisymmetric upper tangential velocity U(r), prescribed deep velocity
Ud(r), and the code's f = constant - gamma r^2/2, the radial PV gradient is

    Q_r = U_rr + U_r/r - U/r^2 - gamma r + F (Ud - U).

The upper flow stays fixed during this comparison. Changing the deep flow
changes Q_r through F(Ud-U), while changing LD also changes the disturbance
inversion. If Ud=U everywhere, the stretching contribution to the base PV
gradient cancels; the perturbation still feels finite LD. Matching only the
annular jet does not also match the present central cyclone.

For a disturbance phi(r) exp[i m(theta-c t)], Omega(r)=U(r)/r gives

    (Omega-c) [d_rr + r^-1 d_r - m^2/r^2 - F] phi - (Q_r/r) phi = 0,
    growth rate = m Im(c).

A one-signed Q_r rules out exponentially growing inviscid normal modes under
the usual regular-origin/no-flux boundary assumptions. A sign change is only
a necessary condition for instability. It does not prove instability, and
nonmonotonic PV can still be stable through the finite-LD wave-interaction
mechanism emphasized by Thomson & McIntyre. Their beta-plane threshold must
not be transplanted into polar geometry.

One useful polar analogue of their second-Arnold-theorem argument is a
**sufficient linear bound**, derived here for an impermeable circular disk
of radius Rb. For a constant angular speed c_star, define

    A(r) = -Q_r / [r (Omega - c_star)].

If, throughout the disk, `0 < A(r) < F + (j_m,1/Rb)^2` with strict uniform
bounds, then a sign-definite conserved linear quadratic excludes growth in
that azimuthal sector. Here j_m,1 is the first positive zero of Bessel J_m;
use m=1 for a bound covering all nonaxisymmetric sectors. The result follows
from E' <= ||q'||^2/[2(F+(j_m,1/Rb)^2)] and conservation of

    E' + (1/2) integral [r(Omega-c_star)/Q_r] |q'|^2 dA.

Any simultaneous zeros in the ratio require a regular limiting value; do
not certify a profile from sampled ratios with singularities. This is a
sufficient linear screen under the stated disk assumptions, not a necessary
threshold or a nonlinear stability proof for the finite hexagon. In
particular Q_r itself changes with F: reducing LD alone is not guaranteed to
improve this bound. For a deep solid-body flow Ud=c_star*r, let
G=U_rr+U_r/r-U/r^2-gamma*r and D=U-Ud. The bound becomes
`-(j_m,1/Rb)^2 < G/D < F`; increasing F cannot repair a failure of its lower
inequality. This explains why both lower-layer structure and LD must be
controlled rather than treating the lower layer as generic damping. These
analytic bounds are documented separately; the implemented tool computes
eigenvalues and sampled PV-gradient signs, not an Arnold certificate.

The separate Python tool `scripts/radial_stability.py` screens the circular
analytic jet/cyclone profile. Its circular impermeable outer boundary is
different from the nonlinear model's periodic square and PV continuation.
Check radial resolution and outer radius, and verify selected growth rates
in weak-perturbation 2D runs. This tool does not linearize about the finite
hexagon: that requires a two-dimensional linear operator in the pattern
frame, or carefully controlled perturbation integrations about a verified
nonlinear equilibrium.

Use a first dimensionless screen with LD/w = {1,2,4} and deep annular speed
Ud,max/Uj = {-0.5,0,0.5,1}, plus a barotropic reference. These values are
exploratory brackets, not observational estimates. Include LD/w = 0.5 only
with adequate resolution. Near any transition, refine LD, then vary deep
width wd/w = {0.5,1,2} and radial offset (Rd-Rj)/w = {-1,0,1}. Keep the upper
profile, cyclone and perturbation fixed. Repeat with the cyclone removed to
separate polar-core and annular-jet instabilities. Examine all resolved m,
including low modes, rather than declaring success when m=6 alone is quiet.

### Import and validation requirements

The current Julia initializer creates a perturbed tanh annular jet; it is not
an implementation of the exact 2024 solution. The located file
`SaturnHexagonIC.m` likewise uses the Ptolemaic boundary
shape to build a Gaussian material jet, rather than importing the exact
Lagrangian PV/velocity family. It is not established as the paper's initializer.

Before adding an exact-wave initializer:

1. Reconcile the published phase/PV formulas with the original calculation.
   Direct differentiation of the printed map exposes a sign/convention
   inconsistency between its phase law and the stated PV equation; the
   printed qJ expression also fails a straightforward zero-amplitude check.
   Reconstruct vorticity from the map's velocity and Jacobian, then verify
   material PV conservation independently. Do not silently patch the formulas
   or label a shape-only initializer an exact solution.
2. Restore dimensions consistently: gamma_paper multiplying r^2 corresponds
   to gamma_code/2, and rescaling x=A0*x_tilde introduces A0^2 in the planetary
   term. Keep particle labels, physical coordinates, angular pattern speed
   and velocity units distinct.
3. Invert the nonsingular Lagrangian map numerically, reconstruct psi and
   velocity, and test inverse-map accuracy, incompressibility, circulation,
   PV, and J(psi,q). Include the ambient-region interface and far-field
   transition in these residuals; a nominal constant-PV patch is not by itself
   evidence of a globally matched steady solution.
4. Reproduce the barotropic reference and its mesh convergence first. Add
   finite LD with resting deep flow as an equilibrium-preserving control.
   Then continue to circular deep profiles and study perturbation growth.
5. Only after identifying persistent finite-amplitude hexagons, add convection
   and test their basin of attraction. Separately ask whether broadband
   perturbations of circular equilibria select six and whether weak initial
   flow organizes into a jet. A stable circular vortex alone does not explain
   a finite-amplitude hexagon.

The goflow velocities will constrain Eulerian U(r,theta), jet width, core
circulation, sixfold deformation and drift. The paper's Lagrangian mean flow
must be converted to Eulerian sampling before comparison. A successful fit
may restrict combinations of LD and deep circulation; surface winds alone
need not uniquely determine either.

## 8. Pilot configuration and paired-suite generation

`configs/pilot.toml` is an illustrative 1024², 60,000-km-box setup, not a fit to
Saturn data. It initializes a 100 m/s jet at radius 15,000 km with a 1000-km
profile width, a central cyclone, a finite 2000-km deformation radius, and a
50 m/s prescribed deep jet. Its deliberately seeded mode-six displacement is
100 km. The five-million-second duration is about 58 Earth days or 5.3 nominal
jet-circulation periods; it is an initial pilot window, not a test of decades of
persistence. The nominal jet-width Rossby number is about 0.3, so the QG
approximation itself deserves assessment when interpreting this chosen profile.

The 250-km convective radius occupies approximately 4.27 cells at this grid
spacing. Passing the minimum-resolution check does not establish convergence.
Keep this physical radius fixed in subsequent resolution tests. The default
outer sponge is disabled; the specified periodic background continuation and
small-scale hyperdiffusion remain part of the numerical model.

Generate the paired controls without running them:

```sh
julia --project=. scripts/make_suite.jl configs/pilot.toml suites/pilot-v1 --seeds 42:44
```

The destination must be new. The default creates 54 configurations: three
physics settings (barotropic/no deep flow, finite deformation/no deep flow,
finite deformation/with deep flow), three forcing settings (off, unbiased,
biased), two initial-condition families (maintenance, broadband selection), and
three seeds. Maintenance keeps the seeded six; selection uses equal spectral
amplitudes and random phases over the initializer's configured mode range.
Unforced maintenance cases with different seeds are identical in the current
initializer and are not independent realizations. They remain in the complete
paired index for bookkeeping and may be omitted from an actual resource request.

Add `--include-emergence` to a fresh destination to create 72 configurations.
The extra 18 cases use only unbiased/biased forcing and set upper-jet speed,
central-cyclone speed and displacement to zero, with perturbation `none`.
The nominal jet radius remains only a diagnostic search location in those
quiescent initial states. Finite-deformation/deep-flow emergence still contains
an imposed deep jet and must be interpreted accordingly. The same pilot duration
is inherited; it is not a guarantee of equilibration from rest.

The generator uses the package's shared input validator and checks every TOML
after serialization without importing GPU dependencies or allocating model
grids. It writes `config_paths.txt` for inspection and a separate `manifest.tsv`
containing `SHA256<TAB>absolute path`, compatible with the scheduler's integrity
checks. The latter is read-only and its digest is printed. These checks validate
configuration consistency, not timestep stability, runtime, or physical results.
No jobs are launched by generation.

Generate the suite at its final location on the execution host. If it is copied
from another machine or moved, regenerate the scheduler manifest with new
absolute paths rather than using the old one. For example, from the project
directory on the cluster:

```sh
bash slurm/prepare_array.sh "$PWD/suites/pilot-v1/manifest-on-cluster.tsv" "$PWD"/suites/pilot-v1/*.toml
```

Use the manifest path and hash printed by that command with the documented
Slurm environment. `config_paths.txt` alone is not a valid scheduler manifest.
Benchmark one reviewed case before assigning an array concurrency or wall time;
the templates do not supply a measured runtime for these physical pilots.

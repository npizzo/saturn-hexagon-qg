# Rostami–Zeitlin–Spiga (2017): a controlled QG formation benchmark

This campaign reconstructs the **circular initial velocity profile** used by
[Rostami, Zeitlin & Spiga (2017), Icarus 297, 59–70](https://doi.org/10.1016/j.icarus.2017.06.006)
and asks whether its QG counterpart develops and maintains a hexagonal jet.
It does not import their mature nonlinear hexagon, reproduce their full
rotating-shallow-water (RSW) equations, or establish nonlinear stability.
Their analytic profiles are sufficient to construct this formation experiment;
a saved two-dimensional mature field is not required.

The paper's distinction is important: a jet-only initial condition develops
vortex streets and destructive meandering, whereas adding a central cyclone
produces a much more persistent polygon. This is a nonlinear saturation result,
not a claim that the initial circular flow is linearly stable. Their Fig. 15
also tests competition among equal-amplitude random-phase azimuthal modes 2–8.
The paper's persistent hexagon retains an excessive positive angular drift
(§5, Table 1). Both morphology and drift must be assessed here.

## Source parameters and dimensional conventions

The source is Eq. (2.10), Fig. 1 and §2.1 of the paper. Its reference units are
`L0 = 3.2e6 m`, `f0 = 3.2e-4 s^-1`, and `U0 = f0*L0 = 1024 m/s`.
The reference planet radius is `a = 55e6 m` and rotation rate
`Omega = f0/2 = 1.6e-4 s^-1`. These differ from earlier generic campaign
values; the difference must not be hidden in a phase-speed comparison.
The approximate `f0` quoted in the paper is used as this benchmark's model
convention. It does not independently establish the observational reference
frame of the fitted winds. Any System III conversion remains conditional until
both wind and longitude conventions are registered consistently.

| Parameter | Published nondimensional value | This benchmark's dimensional interpretation |
|---|---:|---:|
| Jet peak speed | 0.08 | 81.92 m/s |
| Jet radius | 3.37 | 10,784,000 m |
| Gaussian exponent coefficient | 3 | 1/e half-width `L0/sqrt(3) = 1,847,521 m` |
| Reference cyclone peak speed | 0.15 | 153.6 m/s |
| Stronger cyclone peak speed, Fig. 13 | 0.20 | 204.8 m/s |
| Core exponents | alpha=0.42, beta=1.3 | Held fixed |
| Core scale in exponent | 1 | `Lcore = L0` |
| Unregularized cyclone velocity-peak radius | `(alpha/beta)^(1/beta)` | about 1,342 km |
| Reference deformation radius | 1 | 3,200,000 m |
| Illustrated late time | `1000/f0` | 3,125,000 s, about 36.17 Earth days |

All speeds above are peaks of the **individual components**, not the peak of
the sum. `initial.cyclone_radius_m` specifies the unregularized peak radius;
regularization moves the actual peak outward and the audit reports it.
When `deformation_radius_m` changes, the dimensional velocity profile, jet
radius, width, core shape and planetary parameters stay fixed. The reference
scale `L0` does not follow the varied deformation radius. Otherwise a putative
LD experiment would also change the upper wind being tested.

The paper uses `f(r)=f0*(1-r^2/(2a^2))`. This solver writes it as
`f0-gamma*r^2/2`, so `gamma=f0/a^2`. In the prognostic QG PV the constant
`f0` is dynamically immaterial and omitted; it remains essential in the
balance and physical-validity diagnostics.

## Explicit regularization of the published profiles

The raw jet is `Uj(r)=81.92*exp[-((r-Rj)/wj)^2]`. Its value at the origin is
exponentially small but nonzero, which does not define a smooth Cartesian
azimuthal velocity. The initializer uses the odd reflected Gaussian

```text
J(r) = exp[-((r-Rj)/wj)^2] - exp[-((r+Rj)/wj)^2]
Uj(r) = Uj_peak * J(r)/max_{r>=0} J(r).
```

This enforces `Uj(0)=0` and smoothness of `Uj(r)/r`. At the reference jet the
correction is negligible; it is still a documented modification, not an exact
copy of the analytic formula at every radius. The peak is found by solving
`Rj-r*tanh(2*r*Rj/wj^2)=0`.

The paper's normalized core fit is proportional to
`(r/Lcore)^alpha*exp[-(r/Lcore)^beta]`. Because `alpha=0.42<1`, its relative
vorticity diverges at the origin. A Cartesian QG experiment cannot silently
sample that singularity. The selected core regularization is

```text
s(r) = sqrt(r^2 + epsilon^2)
C(r) = (r/s)*(s/Lcore)^alpha*exp[-(s/Lcore)^beta]
Uc(r) = Uc_peak*C(r)/max_{r>=0} C(r)
Lcore = Rc_unregularized/(alpha/beta)^(1/beta).
```

`epsilon = initial.core_regularization_m` is an explicit physical sensitivity
parameter, independent of the grid. `C(r)` is linear at the origin and tends
to the paper fit away from the pole as epsilon decreases. The regularized
peak satisfies
`1+(alpha-1)*r^2/s^2-beta*(r^2/s^2)*(s/Lcore)^beta=0`.
The default 500 km smoothing is not observationally measured. The companion
250 and 1,000 km cases, with grid checks, are therefore part of interpretation.
Smoothing changes both core vorticity and the outer profile slightly through
peak normalization; the saved profile audit quantifies these changes.

The total circular upper streamfunction is obtained by integrating
`psi_r=Uj+Uc`; the solver applies its existing cap taper, gauge, spectral
projection and periodic continuation. Thus the analytic untapered circle is
a continuum QG equilibrium, while the actual Cartesian representation still
requires a quiet, zero-perturbation control on the intended grid. A seed uses
the radial displacement construction of the jet; it is **not** the paper's
RSW eigenfunction, and equal streamfunction/displacement coefficients are
not a claim of equal velocity or energy in every mode.

## QG projection and prescribed deep-flow controls

The reconstruction preserves upper velocity, then computes consistent QG PV:

```text
F = LD^-2
q = Laplacian(psi) - F*(psi-psi_d) + f_gamma
u = -psi_y, v = psi_x, U_theta = psi_r.
```

`deep_profile="matched_components"` sets

```text
psi_d = alpha_j*psi_jet + alpha_c*psi_cyclone
U_d   = alpha_j*U_jet   + alpha_c*U_cyclone.
```

The fractions `deep_jet_fraction` and `deep_cyclone_fraction` are separately
recorded. They apply only to the axisymmetric base, not its polygonal seed.
A resting deep layer `(0,0)` is the finite-LD QG analogue closest to the
paper's one-layer experiment. Comparing `(1,0)`, `(0,1)` and `(1,1)` separates
the roles of deep flow beneath the jet and beneath the cyclone while holding
the upper winds fixed. The deep flow is prescribed: this is not a dynamically
evolving second layer and cannot model lower-layer feedback.

For circular flow the PV gradient is

```text
q_r/r = (U_rr + U_r/r - U/r^2)/r - gamma + F*(U_d-U)/r.
```

Complete matching cancels the stretching contribution to the **base** PV
gradient. It does not remove F from perturbation inversion and does not itself
prove stability. The radial linear tool solves

```text
(Omega_flow-c)*L_m(phi) - (q_r/r)*phi = 0
Omega_flow = U/r
L_m = d_rr + d_r/r - m^2/r^2 - F
perturbation = Re[phi(r)*exp(i*m*(theta-c*t))]
growth rate = m*Im(c), angular drift = Re(c).
```

The Python formulas independently differentiate both regularized profiles and
use cancellation-free expressions for `U/r` and `zeta_r/r` at the pole.
`stabilization_linear.py` reads the exact hashed configurations. The standalone
`radial_stability.py` also supports the profile and component matching, but
rejects global fraction overrides for `matched_components`; edit explicit
configs instead of allowing an ambiguous control to be silently ignored.
Its artificial circular Dirichlet boundary, second-order radial grid and
omission of dissipation differ from the Cartesian simulation. Inspect radial
resolution, outer-radius sensitivity, critical levels and the high-mode edge.
A non-growing finite spectrum is not a certificate of nonlinear stability.

## The central cyclone challenges the QG approximation

RSW requires the gradient-wind balance

```text
g*H_r = f(r)*U + U^2/r,
```

whereas the QG construction neglects centrifugal acceleration and linearizes
thickness/PV. Holding velocity fixed therefore changes pressure/thickness and
PV. The reference unregularized cyclone has a centrifugal/Coriolis ratio
about 0.36 at its velocity peak, already appreciable. More seriously, for the
500 km regularized reference core, the analytic audit gives approximately

- `max |zeta/f(r)| = 1.70`;
- `max |U/(r*f(r))| = 0.85`, attained at the pole as a regular limit;
- actual component peak radius about 1,520 km;
- only 6.4 cells across the smoothing radius on the 80 Mm/1024 grid.

This core is outside a controlled small-Rossby-number QG limit. A functioning
QG run remains a mechanistic experiment, not a faithful RSW reproduction.
Matching the lower-layer core can reduce interface displacement but **cannot
reduce upper vorticity or centrifugal acceleration**. Reduced-amplitude core
controls and core/grid sensitivity are needed to decide whether a result is
robust or relies on an unresolved, strongly non-QG center. The `qg_regime`
stage sets the core speed to 15.36, 38.4 and 76.8 m/s (0.1, 0.25 and 0.5 of
the paper reference) while keeping the jet, core geometry and 500 km smoothing
fixed. These deliberately depart from the paper's observed core fit; they test
how the mechanism changes as the center approaches the QG regime. They do not
replace the reference/strong-core experiments or justify a claim of agreement
with the observed cyclone. Saved-field diagnostics determine the actual
projected Rossby numbers and interface displacement.

Run the reproducible radial audit before interpreting a candidate:

```bash
python scripts/rostami_audit.py PATH_TO_CONFIG.toml --output NEW_AUDIT_DIRECTORY
```

It writes `audit.json` and `profiles.csv`, checks the exact upper and deep
profiles, reports physical scales and flags large Rossby numbers. Its pressure
integrals deliberately use `psi(R)=psi_d(R)=0` at the untapered disk edge.
**That gauge differs from the evolution model's full-square zero-mean gauge.**
The resulting radial pressure drop is useful for comparing omitted RSW balance
terms, but is not the saved-field value of `min_h_over_H` or an inferred
physical negative layer depth. The evolution diagnostics must supply actual
`max_abs_hprime_over_H`, `min_h_over_H` and `max_abs_zeta_over_f` after projection.
The initial conditions, perturbations and subsequent state can all affect them.

## Interpreting the experiment series

Start with software validation and circular quiet controls; then compare
jet-only, reference cyclone, stronger cyclone, and reduced core amplitudes.
For formation, compare a sixfold seed with modes 2–8 and independent phases;
a wider modal band tests whether restricting the seed selects the outcome.
No claim of spontaneous selection follows from a sixfold-only seed. The
initial source profile is unstable by design: distinguish growth of the circle,
nonlinear saturation into a polygon and persistence of that resulting polygon.

Introduce the deep component controls at fixed upper winds, initially without
convection, so deterministic changes can be attributed. Only after a candidate
survives the unforced checks should paired convection experiments be interpreted.
The small convection-on/off difference in earlier generic runs does not settle
the role of convection for this different PV distribution and forcing budget.

Track the jet's connected shape, all competing azimuthal modes, mean winds,
core strength, PV rearrangement, energy and numerical drift. At minimum test
core size, timestep, grid, physical dissipation and outer boundary for any
claimed persistent candidate. The paper's late snapshot is 3.125 Ms: 400 s
validation and short pilots do not probe that persistence horizon. All trials
retain their stated finite observation horizon; persistence is not a proof of
asymptotic stability.

Track pattern speed only while the relevant mode remains coherent, and state
its fit interval and uncertainty. The paper's Table 1 first jet+cyclone case
has `c/f0=7.83e-3`, or about `2.506e-6 rad/s` in its stated rotating frame.
Do not compare it directly with System III drift without registering longitude
sense, planetary rotation rate and observational convention. In an eastward
positive convention, changing from frame Omega_model to Omega_obs gives
`c_obs=c_model+Omega_model-Omega_obs`; the upper and prescribed deep winds must
be transformed consistently too. A frame adjustment cannot independently fix
pattern speed while leaving an incompatible upper wind unexplained. Future
goflow estimates should constrain the full wind profile and phase speed jointly.

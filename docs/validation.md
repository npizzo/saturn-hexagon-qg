# Initial numerical validation

These are CPU checks of the implemented equations and diagnostics, completed on
2026-10-24 with Julia 1.10.12, two Julia/FFTW threads, and the supplied dependency
Manifest. They establish useful numerical baselines. They do **not** demonstrate
a persistent hexagon, spatial convergence of a turbulent forced run, or GPU
execution.

The measured values and source/Manifest SHA-256 hashes are in
the locally retained `validation/numerical_checks.csv`. The
reproducible driver is [`scripts/validate_numerics.jl`](../scripts/validate_numerics.jl).
From the project directory, after instantiating its environment, run:

```sh
julia --project=. --threads=2 scripts/validate_numerics.jl --output validation/recheck.csv
```

The driver refuses to overwrite a result. Omitting `--output` selects
`validation/numerical_checks.csv`. A `.partial` file is published as the requested
CSV only after all four cases finish and the source and Manifest hashes are
confirmed unchanged.

## Analytic and dynamical unit checks

All **44 assertions** in [`test/dynamics_tests.jl`](../test/dynamics_tests.jl)
passed in the local CPU environment:

- Analytic Helmholtz inversion, velocity signs, and energy normalization for
  finite deformation radius and the barotropic limit: 10 assertions.
- Westward Rossby-wave phase and fourth-order RK4 timestep convergence in both
  limits: 14 assertions.
- Inviscid energy and **total-PV** enstrophy conservation with a nonuniform
  periodic background PV: 8 assertions. The example changes anomaly-only
  enstrophy, so that incorrect invariant cannot pass as a substitute.
- Prescribed deep flow changes background PV while preserving the common initial
  upper velocity; it disappears from the unforced barotropic-limit dynamics:
  12 assertions.

These 44 assertions are the dynamics subset; they are not a count of every test
in the package.

## Short inviscid polar integrations

The driver uses the physical profiles in [`configs/smoke.toml`](../configs/smoke.toml),
including its imposed weak mode-six disturbance, finite deformation radius
`L_D = 2,000 km`, and prescribed deep flow. It disables convection, the boundary
sponge, and hyperviscosity, keeps all physical lengths and amplitudes fixed, and
changes only grid resolution. Each run covers 4,000 s with `dt = 20 s`. A further
64² run uses `dt = 10 s`.

Energy is `E = <|u|² + psi²/L_D²>/2`; total-PV enstrophy is
`Z = <(s+B)²>/2`. Angle brackets denote the whole periodic square, including its
outer continuation. The table reports the maximum absolute relative departure
from the initial value, sampled after every step.

| Grid | Grid spacing (km) | Timestep (s) | Maximum `|ΔE/E₀|` | Maximum `|ΔZ/Z₀|` |
| --- | ---: | ---: | ---: | ---: |
| 64² | 937.500 | 20 | 3.963 × 10⁻¹¹ | 1.330 × 10⁻⁹ |
| 128² | 468.750 | 20 | 8.882 × 10⁻¹⁶ | 3.775 × 10⁻¹⁵ |
| 256² | 234.375 | 20 | 4.441 × 10⁻¹⁶ | 8.882 × 10⁻¹⁶ |
| 64² | 937.500 | 10 | 3.963 × 10⁻¹¹ | 1.330 × 10⁻⁹ |

The maximum absolute change in spatial-mean `s` is below `7 × 10⁻²¹ s⁻¹` in
these runs. The finer-grid invariant changes are near floating-point precision.
Halving the timestep does not materially reduce the coarse-grid drift, which
points to a spatial/discrete-projection contribution; this check does not isolate
its cause. The separate analytic Rossby-wave test supplies the evidence for
fourth-order temporal convergence.

This is a short integration: 4,000 s is about 0.13 of the nominal jet-width
advection time `3,000 km / (100 m s⁻¹)`, and about 0.0053 of the nominal azimuthal
circulation time. It is not a long-time instability, stationarity, or boundary
convergence test. Small invariant drift also does not establish that the physical
solution or its selected polygonal mode has converged.

## Axisymmetric baseline for the convection limiter

A separate initial state at each resolution sets `initial.perturbation = "none"`.
The driver samples its total PV at `t = 0`; it does not inject vortices or evolve
this state. An exactly axisymmetric continuum field has zero azimuthal anomaly.
The measured residual therefore includes Cartesian/spectral representation and
interpolation errors; it is not physical eddy PV.

An initial implementation bilinearly sampled an array produced by
`radial_anomaly`, which subtracts a radial-bin mean. This produced an unacceptable
baseline for the limiter and was replaced. The **production** limiter now calls
`local_pv_anomaly`: it uses the same fourth-order tensor cubic interpolant at the
event point and around a circle at that point's exact radius, then subtracts the
circle mean. The default circle quadrature has
`max(64, ceil(4πr/min(dx,dy)))` angular samples. `radial_anomaly` remains available
as an approximate grid diagnostic and is not used by the production limiter.

The comparison below evaluates both methods at identical points: 16 radii with
equally spaced `r²`, including both annulus endpoints, and 48 angles
`theta = 2πk/48`, `k = 0,…,47`. Equal sample weights approximate an area-weighted
annular statistic. The jet annulus is `R_jet ± width`, or 9,000–15,000 km. The
convection annulus is the configured **pair-midpoint** region,
3,000–10,000 km. All entries are ratios to `q_max = 2 × 10⁻⁶ s⁻¹`.

| Grid | Legacy jet RMS | Production jet RMS | Legacy convection RMS | Production convection RMS | Production convection maximum |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64² | 1.4521 | 0.04185 | 0.5043 | 0.01816 | 0.06165 |
| 128² | 0.4143 | 0.002437 | 0.1267 | 0.0005346 | 0.001725 |
| 256² | 0.1057 | 0.0001600 | 0.03022 | 0.00003224 | 0.0001386 |

The replacement reduces the 64² convection-region RMS from **50% to 1.8% of the
cap** on these samples. At 256² the production RMS is **0.0032% of the cap**, with
a sampled maximum of **0.014%**. The CSV also retains the production jet maxima,
legacy sampled maxima, and the original legacy grid-point RMS/maxima. Grid-point
statistics can differ from event-point statistics because bilinear interpolation
at an annulus endpoint touches cells beyond that annulus.

The legacy method uses `max(16,nx÷2)` radial bins: 32, 64, and 128 here. Its
original grid-point convection RMS values were 0.4286, 0.1261, and 0.03019 times
`q_max`; those early measurements motivated the correction. The production
method has no radial bins. Its remaining error decreases strongly with grid
refinement in this example. Separate sampler tests check known angular modes
and approximately fourth-order convergence for a smooth analytic radial field.

These are sampled initial-state errors, not upper bounds over arbitrary event
positions or a later turbulent flow. Before interpreting forced experiments,
repeat the check for the actual physical profile, grid, convection region, and
`q_max`, and check circle-quadrature and grid sensitivity. The pilot uses different
physical settings; this table neither certifies that pilot nor supplies an
extrapolated error estimate for it. The smoke case remains a software integration
check, not a converged experiment on convective feedback.

## Initial jet-shape diagnostic sensitivity

The same three resolutions also diagnose the *perturbed* initial state before
time integration. The angular and radial peak-search settings are held fixed at
the smoke values, `ntheta = 96` and `nr = 32`. All cases detect a jet at every
angle. The following are measured radial-displacement amplitudes in kilometres,
not Fourier amplitudes of PV or velocity:

| Grid | `m = 2` | `m = 4` | `m = 6` | `m = 10` | `m = 12` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64² | 65.339 | 0.743 | 293.777 | 41.321 | 41.348 |
| 128² | 2.546 | 3.013 | 356.162 | 5.280 | 19.046 |
| 256² | 0.318 | 1.417 | 353.337 | 0.738 | 17.877 |

The large change in the secondary `m = 2` and `m = 10` components shows why the
coarse smoke diagnostic must not be read as physical mode competition. The
reported `m = 6` amplitude also changes substantially between 64² and 128².
These results combine the representation of the flow with bilinear polar
sampling and peak fitting; they do not isolate one source of error.

Not every non-six component is an error: a finite perturbation and the nonlinear
map from streamfunction to a velocity-maximum curve can generate higher
harmonics, including `m = 12`. The configured streamfunction-displacement
parameter is not an exact promise for the measured peak-curve amplitude.
Convergence of the physical grid, angular sampling, and radial peak search should
all be established before estimating growth rates or small secondary modes.

The next scientific checks are long-time and spatial convergence of growth,
phase speed, and mode competition; sensitivity to the artificial outer region;
realized forcing and dissipation budgets; and CPU/GPU agreement. Follow
[`experiment-design.md`](experiment-design.md) before interpreting a retained or
selected hexagon.

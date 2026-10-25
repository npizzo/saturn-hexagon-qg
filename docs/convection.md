# Impulsive convection option

This optional forcing is a **Thomson–McIntyre-inspired polar adaptation**, not
an exact reproduction of the 2016 Jupiter experiments or a physical convection
model. The numerical dynamics remain one active QG layer above a fixed deep flow.
The reference is Thomson & McIntyre (2016), equations (3.1)–(3.10):
[published author PDF](https://www.damtp.cam.ac.uk/user/mem2/thomson-mci-jupiters-unearthly-jets-jas16.pdf).

Each event chooses a midpoint uniformly **by area** between `r_min_m` and
`r_max_m` about the polar origin. Its radius is the square root of a uniform
sample in squared radius. Two compact parabolic vortices of radius
`vortex_radius_m` are placed on the local tangent, separated by
`2*vortex_radius_m + max(dx,dy)`. Cyclone/anticyclone order alternates between
events. The positive vortex is cyclonic for the north-polar sign convention.
The midpoint is in the configured annulus; the finite vortex disks can extend
beyond it. Keep the entire pair safely inside the physical polar region and
away from the numerical buffer. The caller validates that support and its
resolution before running.

The midpoint's nondimensional coldness is

```
z = (psi_deep - psi) / psi_limit_m2_s
R = clamp(0.5 + z, 0, 1)
b = bias_max * R
B = 1 - b
```

The first two lines implement the three-piece ramp. The initial anticyclone
magnitude is `A = q_max_s*R`; the cyclone is `B*A`. The midpoint's instantaneous
PV anomaly `qprime` relative to its azimuthal mean imposes the cap

```
max(B*A + qprime, A - qprime) <= q_max_s.
```

Thus `A` is reduced to the minimum of its initial value, `q_max_s+qprime`, and
`(q_max_s-qprime)/B` when `B>0`. It is never negative. If the background already
has `abs(qprime)>=q_max_s`, this implementation suppresses the event entirely.
The cap concerns the common sampled **midpoint** anomaly and nominal vortex
peaks. It is not a pointwise bound on total PV over the resolved disks, where
the background can vary. The driver obtains `qprime` directly from total PV
using `local_pv_anomaly`: identical periodic tensor cubic Lagrange interpolation
at the midpoint and on a circle at exactly the midpoint's radius. The angular
mean defaults to `max(64,ceil(4π*r/min(dx,dy)))` samples; it does not use radial
bins. A smooth axisymmetric background still leaves a fourth-order Cartesian
interpolation residual, which must be checked against `q_max_s` at the chosen
resolution. `psi` and `psi_deep` are sampled by periodic bilinear interpolation
and must have consistent gauges; otherwise an
arbitrary constant offsets the convection sensitivity.

Unlike the paper's separate constant-bias anticyclone-only option, this driver
uses the continuous `b=bias_max*R` rule even at `bias_max=1`. That value gives
anticyclones alone only where the coldness ramp saturates. `bias_max=0` gives
equal nominal peaks. Annular sampling, tangent orientation, and an azimuthal
mean for `qprime` replace the paper's channel geometry and zonal mean.

The grid-sampled increment is compensated by subtracting its domain mean;
a final roundoff residual is removed. This preserves the preexisting mean of
the evolved anomaly `s`. The compensation is spatially uniform and generally
nonzero for biased pairs. It does not impose zero mean on the preexisting
flow. Returned vortex peaks refer to the compact profiles **before** this
compensation. The compensating term is part of the actual applied increment.

## Timing and reproducibility

`ConvectionState(cfg, seed)` contains an independent `MersenneTwister`, the
event count, copied configuration, and next event time. The first event is
uniformly scheduled in `[interval_min_s, interval_max_s]` after time zero.
Each due call to `kick!` attempts exactly one pair and samples the next
interval from the same range. A warm or capped event still counts as an
attempt and advances the schedule and alternating order. Disabled or early
calls change neither field nor random state. The integration driver must land
on the next event time rather than accumulate several overdue events.

Use no random draws in Runge–Kutta substeps. Apply the kick between complete
steps, then update the inversion and derived variables. Save and restore the
entire `ConvectionState`, including the RNG, together with the model state and
time; restarting from only the original seed would repeat the wrong sequence.
The convection unit tests check exact repeatability within one Julia/runtime
setup; the full solver restart is checked to numerical tolerance. Different
Julia versions are not promised identical random streams.

The forcing works on CPU `Matrix{Float64}` arrays in `(x,y)` order. A GPU driver
copies fields to the CPU at these sparse events, uses the direct common-radius
cubic PV-anomaly sampler, applies the kick, and copies the updated `s` back.
The standalone `kick!` API also accepts a supplied anomaly array sampled
bilinearly when no `qprime_sampler` callback is provided; the simulation driver
always supplies the callback above. The returned record includes location, coldness, sampled
anomaly, bias, nominal signed peaks, mean compensation, and actual increment
extrema/RMS. Compute the energy jump before/after inversion in the caller;
neither positive energy input nor statistical equilibration is guaranteed by
this heuristic.

## Configuration

| Key | Meaning |
|---|---|
| `enabled` | Boolean; false disables events |
| `r_min_m`, `r_max_m` | Midpoint sampling annulus radii in metres |
| `vortex_radius_m` | Radius of each compact parabolic profile in metres |
| `q_max_s` | Positive nominal amplitude/anomaly cap in inverse seconds |
| `bias_max` | Maximum fractional anticyclonic bias, between zero and one |
| `psi_limit_m2_s` | Positive coldness sensitivity scale in m²/s |
| `interval_min_s`, `interval_max_s` | Positive inter-event interval bounds in seconds |

Test the unforced model first. For forced experiments inspect event-resolved
energy changes, mean compensation, cap activations, and energy/enstrophy
budgets; establish sensitivity to the injection radius, time interval, and
spatial resolution rather than assuming the forcing replicates TM16's
statistical steady states.

# Comparing the simulation with winds and the hexagon

The working input contract is an axisymmetric, azimuthal wind profile. Full
gridded spacecraft winds, cloud images, masks, and uncertainty products need a
separate adapter after their actual file schema and metadata have been supplied.
This project does not assume a schema for a particular wind-retrieval package.

## Profile input

`load_wind_profile(path; required_radius_range=(rmin,rmax))` accepts a plain
two-column CSV:

```csv
# Units, coordinate conversion and reference frame must be recorded separately.
radius_m,u_theta_m_s
0,0
1000000,25
2000000,10
```

These numbers illustrate the format; they are not Saturn observations. Blank
lines and whole-line comments beginning with `#` are allowed. The header and
column units are literal. Radii must be finite, nonnegative and strictly
increasing, and speeds must be finite. The optional coverage interval is checked;
there is no silent extrapolation. Without that keyword the loader reports the
available interval and cannot decide whether it covers an experiment. A file
starting outside the polar vortex does not constrain the vortex or the pole.
Any continuation of missing regions is a named modeling assumption.

The return value contains `radius_m`, `u_theta_m_s`, `radius_min_m` and
`radius_max_m`. The loader makes no frame or coordinate conversions.

To compare a converted observational profile with a model profile in this same
format, run:

```sh
julia --project=. scripts/compare_profile.jl observed.csv model.csv summary.csv
```

Omit the final path to print the summary CSV to standard output. The script
linearly interpolates model speeds at observed radii inside the model's supplied
interval, with no extrapolation. It reports sample counts, compared radial range,
coverage, unweighted RMSE and bias (model minus observation). Partial coverage is
reported explicitly; zero overlap is an error. The RMSE describes those sampled
radii only, not an area-weighted polar-cap score. This profile comparison cannot
measure a hexagonal mode, phase drift, or the accuracy of a pixel-level wind
retrieval. Those need time-resolved, gridded fields and the adapter described above.

## Coordinate and frame contract

The simulation uses a Cartesian plane centered on the north pole. Arrays are
indexed `[x_index,y_index]`, with uniformly spaced, endpoint-excluded periodic
coordinates in meters. The convention is

\[
u=-\psi_y,\quad v=\psi_x,\quad
u_r=(xu+yv)/r,\quad u_\theta=(-yu+xv)/r.
\]

The angle is `atan(y,x)`, measured counterclockwise from positive x. Positive
`u_theta` is counterclockwise when viewed from above the north pole. Define which
longitude is on positive x and whether the observation map uses eastward or
westward longitude before comparing a phase or a drift. Rotating or reflecting
an image without transforming its vector components changes the inferred winds.

Record the exact rotation/reference frame, longitude convention, latitude type,
planet shape/radius, map projection, observation time, and wind uncertainty with
every observational input. Convert cloud velocities and pattern motion into the
same frame. A change of rotating frame changes both the jet speed and pattern
speed; in a reduced-gravity model the prescribed deep wind must be transformed
consistently as well. Do not subtract a pattern drift from just one of them.

For a spherical geodesic polar mapping, `r=a*(π/2-φ_c)`, where `φ_c` is
planetocentric latitude in radians and `a` is the chosen spherical radius. A
planetographic latitude cannot be substituted directly. For an oblate spheroid
the relation is `tan(φ_g)=(a_eq/a_p)^2*tan(φ_c)`; use the documented radii before
applying any further spherical idealization. A stereographic or orthographic
map has a different radial coordinate and metric. Reproject positions and rotate
and scale vectors with that projection's Jacobian; changing positions alone is
insufficient. The plane is an approximation to a spherical polar cap.

## Shared diagnostic definitions

`polar_diagnostics` returns azimuthal means of radial and azimuthal velocity and
total PV over an annulus, plus the radial/azimuthal eddy covariance
`mean((u_r-mean(u_r))*(u_theta-mean(u_theta)))` in m²/s². `jet_halfwidth_m` is the
halfwidth of the *search annulus*, not a fitted width. It must contain the whole
meandering jet and enough of its flanks. The returned `radii` cover that annulus.

At each angle the jet is the interior positive maximum of `u_theta(r,theta)`,
refined by a three-point parabola. Flat, missing, negative-only, or
boundary-peaked profiles are invalid. `coverage` is the fraction of angles with
a valid jet; `radial_coverage` separately gives the valid velocity fraction on
each ring. Complete angular coverage is required for the simple Fourier
coefficients. Otherwise coefficients are NaN, because zero-filling missing
sectors would create artificial modes. A nonzero partial `meanjet` is not a
successful full-circle detection; inspect `jet_detected` and `coverage`.
The default `min_peak_contrast_fraction=1e-3` also rejects radial speed ranges
below 0.1% of peak absolute speed, preventing interpolation ripples in a flat
profile from masquerading as a jet. Report and sensitivity-test this threshold
for weak observed jets; it is not a statistical significance test.

The Fourier convention is

\[
c_m=N_\theta^{-1}\sum_j[r_j-\bar r]e^{-im\theta_j},\qquad A_m=2|c_m|.
\]

Thus `r(theta)=r0+A*cos(m*theta+φ)` has `c_m=(A/2)*exp(iφ)`. For a pattern moving
counterclockwise at angular speed `Ω_pattern`, its coefficient phase decreases
at `m*Ω_pattern`; infer `Ω_pattern=-d(unwrapped angle(c_m))/dt/m`. Unwrap the
*coefficient* phase before dividing by m. Weak or intermittent mode amplitude
requires a quality mask; phase of a near-zero complex coefficient is not a
measurement. `meanjet` is a radius in meters, and `mean_peak_speed` is m/s.

The polar wind/jet diagnostic sampler is periodic bilinear interpolation. It propagates missing
corners rather than inferring a wind vector. This numerical interpolation does
not replace uncertainty propagation. Future gridded adapters must preserve
quality masks and sampling density, define minimum coverage, and estimate mode
and drift uncertainty, for example by resampling the supplied wind uncertainties.
Uneven angular coverage needs an explicitly weighted fit and an assessment of
mode identifiability, not the present complete-circle DFT.

Apply the same resolution, smoothing/filtering, annulus, interpolation, and
jet-tracing rule to model and observation fields before comparing amplitudes,
curvature, PV gradients or eddy fluxes. Record the effective filter width in
kilometers. Derivatives amplify noise; wind-derived observational PV also needs
a thickness/deep-flow assumption and should not be treated as directly observed
total PV. A cloud brightness boundary and a velocity maximum are different
observables; report which defines the measured hexagon.

`radial_anomaly` subtracts an interpolated radial-bin mean. This is an
approximation on the Cartesian grid, not an exact azimuthal projection. Its
residual for an axisymmetric field depends on binning and resolution and can
be significant beside a narrow jet. It is retained for grid diagnostics and
is **not used by the convection limiter**.

For that limiter, `local_pv_anomaly(field,x,y,xp,yp; ntheta=nothing)` instead
subtracts the field's azimuthal mean at precisely `r=hypot(xp,yp)` from its
local value. Both terms use the same periodic tensor cubic Lagrange
interpolant. The circle has `max(64,ceil(4π*r/min(dx,dy)))` equally spaced
samples by default, without a repeated endpoint. Explicit angular resolution
can be supplied for convergence checks. A required nonfinite interpolation
stencil makes the anomaly NaN; missing circle sectors are not zero-filled.
Keep the full circle inside the usable polar cap, despite periodic wrapping
being supported numerically.

This avoids radial-bin bias but is still a numerical azimuthal projection.
A smooth axisymmetric field retains a fourth-order Cartesian interpolation
residual. The test suite measures a unit-amplitude Gaussian PV ring of radius
2 and standard deviation 0.35 in an 8-wide Cartesian box, plus a constant and
quadratic radial background. Across 13 radii from 1.3 to 2.7 and 48 angular
positions, the residual RMS is approximately `1.53e-4`, `1.03e-5`, and
`6.30e-7` at grid sizes 64, 128, and 256, respectively. These are interpolation
checks for that analytic field, not error bounds on a Saturn experiment.
Measure the actual jet's residual relative to `q_max_s`, and refine both the
grid and angular sampling before interpreting limiter activity physically.

## Sources for the observational context

- [Antuñano et al. (2015), Dynamics of Saturn's polar regions](https://doi.org/10.1002/2014JE004709): polar wind fields and hexagon geometry.
- [Sánchez-Lavega et al. (2014), The long-term steady motion of Saturn's hexagon](https://doi.org/10.1002/2013GL059078): pattern motion and seasonal persistence.
- [Rostami et al. (2017), On the dynamical nature of Saturn's North Polar hexagon](https://doi.org/10.1016/j.icarus.2017.06.006): reduced-model jet/vortex comparison and pattern-speed limitations.

The conventions and quality controls above are this project's explicit analysis
contract, not a claim that these papers use the same data-file format.

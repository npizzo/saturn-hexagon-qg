# Readiness of the two published initializer routes

Neither route is currently a verified, complete published-state initializer
for this QG solver. Both are useful next steps, with different reconstruction
requirements. These requirements **do not block the existing annular-jet
campaign**, whose initializer and claims are documented separately in
[unity-stabilization.md](unity-stabilization.md).

No saved full state is available. The route forward is to reconstruct the
fields from the published families and the supplied formula. In particular,
the ambient PV/velocity continuation and Eulerian pattern rotation are
**unknown family parameters to specify or solve for**, subject to interface,
boundary and equilibrium compatibility. They are not fixed by the interior
PV formula alone, and neither zero ambient PV nor zero pattern rotation should
be silently assumed.

## What is established

| Route | Established basis | Remaining gap |
|---|---|---|
| Pizzo–Salmon (2024) | Exact barotropic reference family; the supplied legacy code's unshifted active `q0` agrees with independently reconstructed planar PV under explicit conventions. | The legacy script's output convention, ambient construction, pattern rotation, spherical conventions and solver input semantics have not been verified as a complete equilibrium. |
| Rostami–Zeitlin–Spiga (2017) | A persistent nonlinear hexagonal jet in one-layer rotating shallow water, with the central cyclone included. | A shallow-water state is not automatically a QG equilibrium; its complete fields and the changes introduced by QG projection must be established. |

For [Pizzo & Salmon (2024)](https://doi.org/10.1017/jfm.2024.932), the attached
legacy code actively uses `n=2`, `eps=0.05`, `a0=0.3`: its symmetry is
`m=n+1=3`, not six. The published hexagonal example uses `n=5`. Agreement of
the unshifted formula at the audited points does not validate the script's
subsequent minimum shift, interpolation and truncation as a complete state. The original reader's
meaning of `q`, its mean convention, and the coordinate conversion between
the planar formula and spherical chart remain unresolved. The
[legacy audit](legacy-ic-audit.md) records the verified formula, provenance,
sign convention and limitations; no further access to the restricted legacy
folder is required for this readiness assessment.

For [Rostami et al. (2017)](https://doi.org/10.1016/j.icarus.2017.06.006),
centrifugal acceleration enters the gradient-wind balance, particularly for
the central cyclone. Their nonlinear persistence result is not a proof of
stability against arbitrary perturbations or a ready-made QG solution. A
conversion that preserves velocity will generally change the pressure/thickness
balance and PV; preserving their shallow-water PV instead generally changes
velocity. The nominal deformation radius alone does not specify the state.

## Gates before treating either route as a QG reference state

1. **Construct a complete, reproducible field.** Specify the source version,
   parameters, epoch, phase, spatial domain and boundary conditions. Reconstruct
   velocity, streamfunction and PV consistently throughout the interior,
   ambient region and outer patch/transition, including interface matching.
   Explicitly parameterize the missing ambient PV and circulation, and solve
   or justify their matching to the interior together with pattern rotation.
   For the shallow-water route also retain thickness/pressure and the gravity
   and depth parameters defining its balance. A contour, interior PV formula
   or cropped image alone is insufficient. Record any field reconstructed
   from a published profile rather than supplied numerically.

2. **Fix coordinates, dimensions and the observation frame.** Record physical
   length/time scales, SI velocities, hemisphere, metric, latitude definition,
   longitude sign, rotation rate and streamfunction/PV gauges. Distinguish
   planetary rotation, particle circulation frequency and Eulerian angular
   pattern speed. The legacy code's unused `c=0` does not establish a stationary
   Eulerian pattern. On
   the sphere use the spherical metric and Coriolis field; a polar-plane map
   cannot simply inherit a spherical interpretation. Transform the pattern,
   upper winds and prescribed deep winds consistently into one observation
   frame, with epoch and uncertainty, following [observations.md](observations.md).
   Constrain the upper wind field as well as drift. If observed pattern speed
   is imposed, solve for compatible free family parameters and shape; a frame
   change alone cannot repair the physical drift mismatch.

3. **Declare and measure the QG projection.** State whether the construction
   preserves upper velocity, PV or another quantity. Recompute QG PV from
   the actual projected streamfunction, deformation radius and deep field.
   Archive the original and projected fields and quantify changes in velocity,
   divergence, PV, pressure/thickness balance, jet width, cyclone circulation
   and polygon amplitude. Record filtering, interpolation, mean removal,
   smoothing, tapering and boundary changes separately. Evaluate Rossby number
   and interface-displacement smallness; a numerical equilibrium can still lie
   outside a defensible QG regime.

4. **Verify the inviscid moving equilibrium on each grid.** With the current
   polar-plane convention and circular prescribed deep flow, form

   ```text
   F = LD^-2
   qF = (Laplacian - F)*psi + f + F*psi_d
   R = -c*partial_theta(qF) + J(psi, qF)
     = J(psi - c*r^2/2, qF).
   ```

   Here `c` is the Eulerian angular pattern speed in the stated model frame.
   On the sphere use `-c*partial_longitude(qF) + J_sphere(psi,qF)` with the
   model's own spherical inversion. Reconstruct from the same physical state
   at each resolution, then recompute inversion and residual using that grid's
   actual operators and spectral projection. Report dimensional RMS and maximum
   residuals, separately for the jet, ambient/interface region and full domain;
   include boundary compatibility. Also report
   `epsilon_R = T_ref * RMS(R) / Q_ref`, with a declared physical turnover time
   and nonzero PV scale fixed across the grid comparison. Do not normalize by
   a large irrelevant constant PV offset. Require convergence and a residual
   small relative to the intended perturbation signal, with tolerances set
   before interpreting growth. Forcing and dissipation must not conceal an
   inviscid imbalance.

5. **Test perturbed persistence after the equilibrium gate.** Compare an
   unperturbed trajectory with independently perturbed copies of the same
   state, including modes outside sixfold symmetry, several phases and
   amplitudes, and grid-relative rotations. Check timestep, resolution,
   physical dissipation and boundary sensitivity. Separate neutral phase drift
   and baseline adjustment from growing shape disturbances. Track connected
   jet geometry, competing modes, upper winds and pattern speed over a declared
   number of circulation periods. Finite-time persistence is evidence with a
   measured horizon, not asymptotic stability or spontaneous sixfold selection.

Finite deformation and deep circulation must be introduced consistently with
the equilibrium. A verified stationary barotropic state is preserved over a
resting deep layer when its upper streamfunction is held fixed and PV is
recomputed. For a rotating state, a solid-body deep flow `U_d=c*r` is an
equilibrium-preserving polar-plane control; an arbitrary circular deep jet
generally requires continuation of the shape and speed together. The derivation
and the joint drift/wind constraint are in
[experiment-design.md, section 7](experiment-design.md#7-lower-layer-stability-and-the-2024-exact-wave-family).

Passing these gates would make a reconstructed state ready for the proposed
stability and deep-flow experiments. This audit supplies the acceptance
criteria; it does not add an initializer, certify a published state as runnable,
or launch simulations.

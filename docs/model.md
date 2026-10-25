# Model definition and numerical limits

This is a one-active-layer quasigeostrophic model on a Cartesian north-polar gamma plane. It is intended for controlled experiments on a circumpolar jet, a central cyclone, prescribed deep flow, and parameterized convective PV injection. It does not evolve a second layer, resolve convection, represent a spherical metric, or establish that Saturn's observed hexagon is stable.

The physical inspiration combines the finite-deformation/deep-flow formulation of [Thomson & McIntyre (2016)](https://doi.org/10.1175/JAS-D-14-0370.1) with the polar geometry of [Siegelman, Young & Ingersoll (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170063/). The forcing and outer-boundary treatment are adaptations, not exact reproductions of those experiments.

[Pizzo & Salmon (2024)](https://doi.org/10.1017/jfm.2024.932) supplies the
intended nonlinear vortex-wave reference family. The current annular-jet
initializer does not yet implement that family. Its import and lower-layer
stability controls are specified in [experiment-design.md](experiment-design.md#7-lower-layer-stability-and-the-2024-exact-wave-family).

## Prognostic variable and sign conventions

All configuration quantities are SI. Cartesian coordinates are in metres, time is in seconds, streamfunction is in m²/s, velocity is in m/s, and PV is in s⁻¹. Array dimensions are `(x,y)`; the polar origin is `(0,0)`.

With `r² = x² + y²`, define

\[
\gamma=\frac{2\Omega}{a^2},\qquad
\kappa^2=L_D^{-2},\qquad
s=(\nabla^2-\kappa^2)\psi,
\]
\[
B=f_\gamma(r)+\kappa^2\psi_d(r),\qquad q=s+B.
\]

The velocity convention is

\[
u=-\psi_y,\qquad v=\psi_x,\qquad
u_\theta=\psi_r.
\]

Positive tangential velocity is counterclockwise when viewed in the model's Cartesian plane. Positive relative vorticity is cyclonic at the north pole. Observed east/west longitude conventions and image orientation require an explicit conversion when comparing with goflow.

The evolving equation is

\[
\partial_t s+J(\psi,s+B)
=-\nu(-\nabla^2)^p s+F_{\rm sponge}
+\sum_j\Delta s_j\,\delta(t-t_j),
\qquad
J(\psi,q)=\psi_xq_y-\psi_yq_x.
\]

The code evolves `s`, **not** total PV `q`. The prescribed field `B` remains fixed. Dissipation acts on `s`; it does not diffuse the prescribed planetary/deep PV. There is no additional imposed Cartesian mean velocity. The entire upper-layer circulation is part of the prognostic field. The deep flow affects `B` and, when enabled, the forcing coldness; it does not separately advect the upper layer.

The spectral inversion is

\[
\widehat\psi(\boldsymbol k)
=-\frac{\widehat s(\boldsymbol k)}{|\boldsymbol k|^2+\kappa^2}.
\]

For `deformation_radius_m = inf`, `kappa² = 0` and the zero streamfunction mode is set by the barotropic inversion convention. This recovers barotropic vorticity dynamics. A finite deformation radius with a zero prescribed deep profile remains a 1.5-layer model with a resting deep layer. For `matched_upper`, zero `deep_jet_speed_m_s` alone does not imply a resting deep layer; its amplitude is set by `deep_match_fraction`.

## Periodic outer continuation

The square domain is periodic, with endpoint-excluded coordinates from `-L/2` to `L/2-dx`. It contains an interior disk bounded by `inner_radius_m = r_i` and a transition ending at `outer_radius_m = r_o`, where `r_o < L/2`.

The radial window is

\[
W(r)=\begin{cases}
1,&r\le r_i,\\
1-10\xi^3+15\xi^4-6\xi^5,
&r_i<r<r_o,\quad\xi=(r-r_i)/(r_o-r_i),\\
0,&r\ge r_o.
\end{cases}
\]

Planetary PV is constructed as

\[
f_\gamma(r)=-\gamma\int_0^r\rho W(\rho)\,d\rho.
\]

Before spectral projection, this is exactly `-gamma*r²/2` in the interior and constant outside `r_o`. Thus the **PV gradient** goes smoothly to zero in the exterior; there is no boundary jump or rigid circular wall. This differs from the conservative vortex trap in Siegelman et al. It need not prevent disturbances from reaching the exterior or interacting with periodic copies.

The upper and deep initial streamfunctions are multiplied by `W`. Background PV, deep streamfunction, initial `s`, and the sponge reference are projected onto the retained Fourier modes. This ensures consistent products but introduces a small resolution-dependent departure from the analytic radial profiles, including inside the nominal interior disk. Boundary position, width, box size, and resolution must be varied before interpreting mode selection or drift.

An optional sponge restores `s` toward the unperturbed initial axisymmetric `s_ref` outside the interior. Its rate is `sponge_rate_s`, its mask is the spectral projection of `1-W`, its spatial mean tendency is removed, and its tendency is projected again. A truncated smooth mask can have small sign/overshoot ripples; it is not a rigorous nonnegative local damping function. Keep the sponge disabled for initial conservation experiments. When enabled, its energy/PV effects and sensitivity require separate accounting.

## Streamfunction reference and interface interpretation

Initial `psi`, `psi_deep`, and the unperturbed reference streamfunction each have zero area mean over the **entire periodic square**. A constant is also subtracted from `B`; that constant cannot affect advection or the radial PV anomaly used by the forcing. Mean-zero convective compensation and sponge tendencies preserve the initial mean of `s` up to numerical error.

Forcing uses `(psi_deep-psi)/psi_limit_m2_s` as its coldness proxy. Choosing both streamfunction means to be zero is therefore a physical modelling convention for the reference interface, not just an irrelevant plotting choice. Its value has not been constrained by Saturn observations. Changing box size while retaining zero-domain-mean references can shift this proxy in the physical disk, even if the local winds are unchanged. A forced domain-convergence study must inspect that change; independent constraints or an explicit reference-interface parameter would be needed to remove it.

At infinite deformation radius, the deep streamfunction disappears from the **PV dynamics**. However, the current convection prescription still samples it. Deep-present and deep-absent barotropic runs are identical only without convection, or when the prescribed forcing is held identical. A difference caused solely by this proxy is a forcing-control effect, not evidence of dynamical coupling to a deep layer.

## Initial jet, cyclone, and disturbances

Before windowing and zero-mean adjustment, the annular jet uses

\[
\psi_j(r)=Uw\left[\tanh\left(\frac{r^2-R^2}{2Rw}\right)-1\right],
\qquad
u_{\theta,j}(r)=U\frac rR\operatorname{sech}^2\left(\frac{r^2-R^2}{2Rw}\right).
\]

This is regular at the origin and has speed `U` at `r=R`. Its actual peak is slightly displaced from `R`; `w` is a profile parameter rather than a measured full width at half maximum. Fit and compare the realized profile. The default deep profile (`physics.deep_profile = "annular"`, also used when omitted) uses the same ring family with independent parameters.

The central cyclone uses

\[
\psi_c(r)=-U_cR_c\exp\left(\frac12-\frac{r^2}{2R_c^2}\right).
\]

Before windowing, its tangential speed has maximum `U_c` at `r=R_c`. A single-mode or broadband perturbation is added through a streamfunction disturbance localized around the jet. Each mode has a smooth polar envelope proportional to `(r/R)^m exp[-m((r/R)²-1)/2]`. This removes the coordinate singularity of `cos(m theta)` at the origin.

With `physics.deep_profile = "matched_upper"`, the deep streamfunction is
`deep_match_fraction` times the unperturbed **jet plus cyclone**, with the same
taper and mean reference. The fraction defaults to one, must be nonnegative,
and is written explicitly in the stabilization campaign. The independent
`deep_jet_speed_m_s` must be zero in this mode to prevent silently ignoring a
second prescribed flow. Its radius and width remain positive configuration
placeholders. No azimuthal perturbation is imposed in the lower layer.

Changing the match fraction or deformation radius leaves the initial upper
streamfunction and winds fixed; the initial PV is recalculated. At fraction
one, the axisymmetric stretching term cancels: `psi_base = psi_deep`. On the
stored grid the matching branch applies the same projection passes as the
upper PV initializer, avoiding a cutoff-edge mismatch in the pinned FFT
library. The legacy annular branch retains its previous projection behavior.
This matching removes the *base-state* interface anomaly but neither its
perturbations nor the upper flow's Rossby number. It does not establish
stability or small-amplitude QG validity.

New diagnostic columns measure the evolved state inside `inner_radius_m`:
`max_abs_hprime_over_H`, `min_h_over_H`, and `max_abs_zeta_over_f`.
The first two use `hprime/H = (psi-psi_deep)/(2*Omega*LD^2)` in the model's
zero-domain-mean convention; they are NaN in the barotropic limit.
The last uses local `f=2*Omega-gamma*r^2/2` and
`zeta=s+psi/LD^2`. These are regime indicators, not automatic acceptance
criteria; inspect them alongside numerical convergence and the actual fields.

`displacement_m` sets the intended small-displacement scale near the jet, not an exactly prescribed nonlinear contour. A single mode has fixed initial phase; broadband phases use the experiment seed. For `N` broadband modes, each nominal component amplitude is `displacement_m/sqrt(N)`, making the root-sum-square amplitude fixed and the angular RMS approximately `displacement_m/sqrt(2)`. Measure the realized centreline after projection. Use displacements small relative to the jet width for linear selection experiments, and independently check that all initialized azimuthal wavelengths are resolved.

Setting both upper jet and cyclone speeds to zero starts from rest. Perturbations scale with upper jet speed, so a nonzero `displacement_m` alone does not seed motion in that case; enabled convection supplies the initial disturbances. `experiment.purpose` records intent and is not a separate dynamical equation.

## Numerics and timestep guards

The solver uses GeophysicalFlows' `SingleLayerQG`, double precision, two-thirds Fourier truncation, and explicit fourth-order Runge–Kutta stepping. Convective impulses occur **between complete steps**. The driver shortens a step to land on the next event or requested end time. It does not adapt the physical forcing radius or silently switch device.

Hyperviscosity is set from the requested cutoff damping time:

\[
k_c=\frac{2\pi}{L}\lfloor n/3\rfloor,
\qquad \nu=\frac{1}{\tau_c k_c^{2p}}.
\]

This defines the decay time at the stated scalar wavenumber `k_c`; corner modes of the retained square spectral region have larger magnitude and damp faster. Changing resolution with fixed `tau_c` changes `nu`. Decide deliberately whether a convergence experiment holds the physical coefficient or cutoff decay time fixed, and record that choice.

Before a step, the code checks the conservative velocity estimate `dt*(max|u|+max|v|)/dx`, an explicit damping bound, and a background-PV estimate

\[
\omega_B\le\|\nabla B\|_\infty
\max_{\boldsymbol k\ne0}\frac{|\boldsymbol k|}{|\boldsymbol k|^2+\kappa^2},
\qquad dt\,\omega_B\le0.5.
\]

The latter guards against fast linear background-PV waves even when the initial velocity is zero. These checks catch unsafe parameter choices but are not a proof of combined nonlinear RK stability or convergence. Repeat with a smaller timestep and inspect budgets. A rejected step asks for a revised configuration rather than silently changing the model's nominal timestep.

Minimum grid checks prevent obviously unresolved upper jets, forcing disks, and boundary transitions. They are starting limits, not sufficient accuracy criteria. Deep-jet width, cyclone radius, initialized azimuthal wavelengths, and diagnostic interpolation also need convergence checks.

## Conserved quantities and interpretation

For the unforced inviscid periodic model with fixed `B`, the active-layer QG Hamiltonian is

\[
E=\frac12\left\langle |\nabla\psi|^2+\kappa^2\psi^2\right\rangle,
\]

and total-PV enstrophy is `Z = <(s+B)²>/2`. The reported means are over the whole square, including the exterior. `E` is the invariant of this prescribed-background model; it is not the energy budget of two freely evolving physical layers. A PV constant changes the numerical value of `Z` by a constant when mean PV is conserved, without changing its conservation test.

Convection, the sponge, and hyperviscosity alter these invariants. Record actual energy and total-enstrophy jumps across each projected impulse. A positive PV-event amplitude does not guarantee positive energy input, and a long-lived pattern does not establish statistical equilibrium. See [convection.md](convection.md) for the implemented forcing and [experiment-design.md](experiment-design.md) for the scientific controls.

# Scientific basis and scope

The hypothesis tested by this project is that a finite-deformation QG weather
layer with a prescribed deep circulation can admit, select, or maintain a
six-fold circumpolar *jet meander*. Those are separate claims. A mode-six initial
condition that survives briefly is evidence of neither spontaneous selection nor
agreement with Saturn. Compare broadband perturbations, unforced relaxation and
forced maintenance; measure pattern speed, jet shape, mean flow and budgets.

## Primary sources

**Thomson & McIntyre (2016), “Jupiter's Unearthly Jets: A New Turbulent Model
Exhibiting Statistical Steadiness without Large-Scale Dissipation.”**
[Journal paper](https://doi.org/10.1175/JAS-D-14-0370.1),
[author PDF](https://www.damtp.cam.ac.uk/user/mem2/thomson-mci-jupiters-unearthly-jets-jas16.pdf).
The motivating 1.5-layer framework includes prescribed deep jets and stochastic,
interface-sensitive cyclone/anticyclone injections. It explains comparatively
straight Jovian jets. The Saturn question is whether a different background
flow and polar geometry can support a meandering equilibrium. Our optional
convection is a documented polar adaptation, not a reproduction of their channel
experiments, and statistical steadiness must be demonstrated rather than assumed.
See [the forcing description](convection.md).

For the lower-layer stability question, sections 4b and 6f are especially
relevant: the finite-deformation inversion weakens interaction between Rossby
waves, but removing the deep jets also changes the base PV gradient. In their
particular sinusoidal setup, that compensating change removes the beneficial
LD dependence of the stability threshold. Their result is not a statement
that any deep coflow, or small LD by itself, stabilizes any vortex.

**Pizzo & Salmon (2024), “Exact planetary waves and jet streams.”**
[Published paper](https://doi.org/10.1017/jfm.2024.932).
This is the intended nonlinear reference family for the next stage. Its
Ptolemaic interior solution with n=5 has sixfold geometry; published Figure 7
uses epsilon=0.03 and shows approximate preservation over the simulated
interval. The barotropic theory assumes infinite deformation radius and
requires an ambient PV/velocity construction. The article does not map the
full stability boundary or demonstrate spontaneous sixfold selection.
The [experiment design](experiment-design.md#7-lower-layer-stability-and-the-2024-exact-wave-family)
derives equilibrium-preserving finite-LD controls and records the formula and
frame audit needed before importing the exact family. The current pilot
initializer is a different annular-jet profile.

**Dritschel, Hmidi & Renault (2019), “Imperfect bifurcation for the
quasi-geostrophic shallow-water equations.”**
[Paper](https://doi.org/10.1007/s00205-018-1312-7),
[accepted manuscript](https://research-repository.st-andrews.ac.uk/bitstream/10023/18652/1/Dritschel_2018_Imperfect_bifurcation_ARMA_AAM.pdf).
For QG with Helmholtz inversion, their near-circle analysis establishes
branches of m-fold rotating vortex patches at arbitrary deformation length.
Their numerical stability calculations distinguish portions of the twofold
and threefold branches. This supports constructing finite-LD shape equilibria
and assessing their stability separately. The constant-f, patch-PV model
does not supply a stability result for Saturn's sixfold gamma-plane jet or
for an arbitrary imposed deep circulation.

**Barbosa Aguiar et al. (2010), “A laboratory model of Saturn's North Polar
Hexagon.”** [Paper](https://doi.org/10.1016/j.icarus.2009.10.022),
[author PDF](https://people.seas.harvard.edu/~rwordsworth/papers/barbosa2010laboratory.pdf).
Laboratory polygons and an equivalent-barotropic linear analysis motivate
instability and deformation-radius controls. Their fastest mode changes from
about 13 in the barotropic limit to about 5–6 at their quoted 2500-km deformation
scale. Their perturbation inversion uses `π²/L_D²`, so their numerical radius
cannot be transplanted into a `1/L_D²` operator without checking definitions and
the basic-state PV gradient. Linear mode selection is not a nonlinear longevity
or pattern-speed result.

**Morales-Juberías et al. (2011), “Emergence of polar-jet polygons from jet
instabilities in a Saturn model.”** [Paper](https://doi.org/10.1016/j.icarus.2010.11.006).
Multilevel EPIC simulations produce stable polygons through jet instability, but
the vortex-street solution does not match the slow observed propagation under
the observed wind constraints. This illustrates why the geometry alone is an
insufficient success criterion.

**Morales-Juberías et al. (2015), “Meandering Shallow Atmospheric Jet as a Model
of Saturn's North-polar Hexagon.”** [Paper](https://doi.org/10.1088/2041-8205/806/1/L18).
Shallow-jet EPIC simulations reproduce meander morphology and phase speed, with
deep winds as important as cloud-level winds. A prescribed lower circulation is
therefore a useful control, but a 1.5-layer model does not reproduce all the
vertical dynamics of a multilevel atmosphere.
This is an exception to a blanket claim that all persistent numerical
hexagons drift too fast. It does not establish that the stable members of
the 2024 exact family have the observed speed. A further constraint is that
[Antuñano et al. (2019)](https://doi.org/10.1029/2018JE005764) infer vertical
wind shears over 1–500 mbar about four times smaller than those used in the
2015 simulations; matching the drift alone does not validate the vertical
structure.

**Sánchez-Lavega et al. (2014), “The long-term steady motion of Saturn's
hexagon and the stability of its enclosed jet stream under seasonal
changes.”** [Paper](https://doi.org/10.1002/2013GL059078).
For August 2008–January 2014 the measured drift is +0.0129 ± 0.0020 degrees
per Earth day in west-positive System III longitude, equivalent to a small
westward east-positive velocity. The reported period is 10 h 39 m 23.01 s.
Use the angular rate with an explicit epoch and frame. The reported period
uncertainty of 0.01 s does not match straightforward propagation of the
printed drift uncertainty (approximately 0.095 s); reconcile their error
definitions before using it as a strict fitting tolerance.

**Sánchez-Lavega et al. (2021), “Interaction of Saturn's Hexagon With
Convective Storms.”** [Paper](https://doi.org/10.1029/2021GL092461).
The 2018–2020 drift is consistent with System III, and section 6 explicitly
explains how adopting a faster reference rotation near 10 h 33 min changes
the apparent pattern speed to roughly -22 m/s while shifting the winds
consistently. Reference-frame conversion is essential; it cannot physically
repair a model-observation mismatch.

**Rostami, Zeitlin & Spiga (2017), “On the dynamical nature of Saturn's North
Polar hexagon.”** [Paper](https://doi.org/10.1016/j.icarus.2017.06.006),
[author PDF](https://web.lmd.jussieu.fr/~aslmd/pub/REF/2017Icar..297...59R.pdf).
One-layer rotating shallow water on a polar gamma plane supports a persistent
hexagonal jet when the observed central cyclone is included; the jet-only case
eventually deforms. Their reference deformation radius is approximately
3200 km. Excessive pattern drift and low-latitude zonal velocities remain
limitations, and they propose full-sphere calculations to test the latter.
This is the nearest reduced nonlinear hexagon precedent. Our QG approximation
adds a further limitation for a strong central cyclone.

**Farrell & Ioannou (2017), “Statistical state dynamics based theory for the
formation and equilibration of Saturn's north polar jet.”**
[Paper](https://doi.org/10.1103/PhysRevFluids.2.073801),
[full text](https://arxiv.org/abs/1609.06170).
A two-active-layer QG beta-plane statistical closure produces a nearly neutral
wave-six perturbation coupled to a turbulence-maintained jet. The interpretation
uses enhanced effective beta from a sloping deep stable layer. This is a
different equilibration mechanism and a different vertical model; it does not
establish that an inactive lower layer has identical dynamics.

**Siegelman, Young & Ingersoll (2022), “Polar vortex crystals: Emergence and
structure.”** [Paper](https://doi.org/10.1073/pnas.2120486119),
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170063/).
This uses barotropic QG (`L_D=∞`) on a Cartesian gamma plane with a doubly
periodic pseudospectral solver. An outer PV jump provides a conservative vortex
trap; initial-vorticity tapering is a separate operation. It motivates the
polar numerical setting and long-time diagnostics, but Jupiter's vortex crystal
and Saturn's meandering jet are different structures. Our finite deformation,
prescribed deep flow and exterior buffer must each be tested independently.

**Shi & Kang (2026), “Polar vortex dynamics on gas giants: Insights from 2D
energy cascades.”** [Paper](https://doi.org/10.1073/pnas.2500791123),
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12849701/),
[supporting code](https://doi.org/10.5281/zenodo.17392274).
This is a direct computational precedent for forced 1.5-layer QG on a gamma
plane with GeophysicalFlows.jl and an outer damping region. Its study concerns
polar vortex regimes and sensitivity to initial energy, forcing, damping and
deformation radius. Our targeted unstable jet and mode-six diagnostics address
a narrower hexagon question, not a reproduction of that parameter survey.

## Observational tests and novelty limits

Use [Antuñano et al. (2015)](https://doi.org/10.1002/2014JE004709) for polar wind
constraints and [Sánchez-Lavega et al. (2014)](https://doi.org/10.1002/2013GL059078)
for long-term pattern motion. The [observation contract](observations.md) records
the frame, latitude, projection and filtering choices needed before comparison.
An additional modern comparison is the southern decagon reported by
[Sánchez-Lavega et al. (2026)](https://doi.org/10.1126/sciadv.aee4251), at a different
latitude from the polar southern jet in older hexagon studies. It invalidates
the blanket claim that Saturn has no southern polygon.

The useful contribution would be a reproducible regime and sensitivity study,
including the deep-flow dependence of stable meanders and quantitative
comparison with a full-sphere solver. This source list is a focused starting
point, not a claim that the proposed combination is unprecedented or that a
validated Saturn hexagon has already been obtained.

# Saturn hexagon QG

A polar gamma-plane, one-and-a-half-layer quasigeostrophic model for studying
how deformation radius, prescribed deep flow and convective PV injection affect
a Saturn-like polar jet. The research target is a finite hexagonal structure
with the observed pattern drift and upper winds.

**[Browse the results gallery](https://npizzo.github.io/saturn-hexagon-qg/)** ·
[Model](docs/model.md) · [Experiment protocol](docs/experiments-01.md) ·
[Gallery workflow](docs/gallery.md)

## Current results

The first gallery contains two completed million-second Unity integrations,
with and without convection at the same deformation radius and deep flow.
Both amplify the imposed sixfold disturbance and develop rolled-up lobes and
vortices. They do **not** demonstrate a stable hexagonal jet, and they do not
isolate the effect of the deep layer. The gallery shows the eleven actual saved
snapshots per run; animation repeats these frames without inventing intermediate
flow states. Missing jet-tracking measurements remain missing.

The next prepared campaign starts with six unforced comparisons: a barotropic
control, deformation radius 2,000 km with deep speeds 0 and 50 m/s, and radius
4,000 km with deep speeds 0, 50 and 100 m/s. Numerical checks, finite-amplitude
persistence, mode competition and matched convection stages follow the
[documented protocol](docs/experiments-01.md). These prepared cases are not
completed results. Small deformation radii are included as mathematical stress
tests where inferred thickness variations violate the standard QG assumption.

## Run the model

The pinned deployment uses Julia **1.10.12**, GeophysicalFlows.jl **0.16.5**,
Float64, and one CPU process or one CUDA GPU per realization.

```sh
git clone https://github.com/npizzo/saturn-hexagon-qg.git
cd saturn-hexagon-qg
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. scripts/run.jl configs/smoke.toml --output runs/smoke --device cpu
```

Every output directory must be new. The 400-second smoke case is a numerical
integration check, not a Saturn experiment. See [Unity](docs/unity.md) or
[DeltaAI](docs/deltaai.md) for GPU deployment; independent parameter choices
can run as an array of single-GPU jobs.

The new [Unity stabilization campaign](docs/unity-stabilization.md) adds a
prescribed deep profile that can match the entire axisymmetric upper flow,
including its central cyclone. It provides 39 staged unforced tests and a
new CPU/GPU preflight. Keep it in a separate checkout from previous runs:
source hashes have changed, so old checkpoints cannot resume with this code.
The pinned dependencies and checkpoint schema remain unchanged.
Saved outputs record the configuration, source
and dependency hashes, actual field times, diagnostics, convective events and
checkpoint state. A completed job is not a scientific validation certificate.

## Publish a new gallery

The gallery builder reads completed result directories and writes a **new**
static campaign directory. It never advances a simulation or uploads raw model
fields. Use [the gallery instructions](docs/gallery.md) to render, inspect and
publish a campaign. GitHub Pages serves only `gallery/`; raw JLD2 files,
checkpoints, runtime caches and cluster credentials are excluded.

## Scientific documentation

- [Equations, frame, boundary treatment and QG limits](docs/model.md)
- [Convective forcing](docs/convection.md)
- [Lower-layer experiments and phase-speed requirements](docs/experiment-design.md)
- [Whole-profile stabilization strategy](docs/stabilization-strategy.md)
- [Published-state reconstruction requirements](docs/published-state-readiness.md)
- [Observation comparison and goflow interface](docs/observations.md)
- [Literature context](docs/literature.md)
- [Historical initializer audit](docs/legacy-ic-audit.md)

The exact 2024 vortex-wave family has not been imported into the production
initializer. Numerical convergence, physical-regime checks and observational
agreement are separate requirements for every proposed Saturn interpretation.

# Rostami–Zeitlin–Spiga (2017) QG benchmark on Unity

This campaign starts from the paper's circular jet and central-cyclone velocity
profiles, regularizes their singular polar-core fit, and initializes the
corresponding balanced 1.5-layer QG flow. It tests formation, persistence and
drift of a hexagon. It does not import a mature shallow-water hexagon or assert
that the shallow-water and QG dynamics are equivalent. See
[the physics specification](rostami-benchmark.md).

All cases have convection off. The baseline uses a resting deep layer and
deformation radius 3,200 km. The planet radius is 55,000 km and the coordinate
rotation rate is 0.00016 s⁻¹, following the paper's dimensional parameters;
these differ from our earlier generic campaign. The upper velocity parameters
stay fixed when the deep circulation or deformation radius changes. Deep
flow is prescribed; it is not a dynamically evolving second layer.

## Download a separate immutable checkout

Keep the working Unity runtime and old outputs. Use the new directory
`/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-rostami-01` and reuse Julia
1.10.12 and its depot under
`/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity`.
Do not install packages again or resume old checkpoints into this source.

From the Mac:

```bash
ssh -i "$HOME/Downloads/privkey.key" nicholas_pizzo_uri_edu@login.unityhpc.org
bash
```

On Unity, enter the **full 40-character release commit supplied with this
campaign**. The command deliberately refuses an existing destination.

```bash
read -r -p 'Rostami campaign commit (40 hexadecimal characters): ' SATURN_COMMIT
export SATURN_COMMIT
bash <<'BASH'
set -euo pipefail
[[ ${SATURN_COMMIT:-} =~ ^[0-9a-f]{40}$ ]] || { echo 'A full release commit is required.' >&2; exit 1; }
base=/work/pi_nicholas_pizzo_uri_edu
project="$base/saturn-hexagon-qg-rostami-01"
[[ ! -e $project ]] || { echo "Refusing to overwrite $project" >&2; exit 1; }
[[ -x $base/saturn-hexagon-qg/.unity/runtime/julia-1.10.12/bin/julia ]] || { echo 'Existing Julia runtime is missing.' >&2; exit 1; }
[[ -d $base/saturn-hexagon-qg/.unity/julia-depot ]] || { echo 'Existing Julia depot is missing.' >&2; exit 1; }
download=$(mktemp -d "$base/saturn-rostami-download.XXXXXXXX")
trap 'rm -rf -- "$download"' EXIT
curl --fail --location --proto '=https' --tlsv1.2 \
    "https://github.com/npizzo/saturn-hexagon-qg/archive/${SATURN_COMMIT}.tar.gz" \
    --output "$download/source.tar.gz"
mkdir -- "$project"
tar -xzf "$download/source.tar.gz" --strip-components=1 -C "$project"
printf '%s\n' "$SATURN_COMMIT" > "$project/DEPLOYED_COMMIT"
cd -- "$project"
bash slurm/unity/rostami.sh validate configs/rostami_01 validation
BASH
```

This downloads source and checks hashes; it submits no jobs. Stop if any step
reports an error. An interrupted extraction can leave an incomplete directory;
the block intentionally refuses to overwrite it on a repeated attempt.

## Run only the validation stage first

```bash
cd /work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-rostami-01
unset SATURN_RESUME
export SATURN_ROSTAMI_RUNTIME_ROOT=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity
export SATURN_ROSTAMI_RUN_ROOT=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-rostami-01/runs/rostami-01
bash slurm/unity/rostami.sh submit configs/rostami_01 validation
```

This submits three 400-second CPU/GPU comparisons: jet only with resting deep
flow, the reference cyclone with half-matched deep jet, and the reference flow
with both deep components matched. These use 512² points and a 500-km core
regularization spanning only 3.2 cells. They are software checks; the science
cases use at least 1024² points. A passed CPU/GPU comparison is not a claim of
spatial convergence or physical stability.

The wrapper retains the tested stabilization validator and immutable launch
controller. Validation jobs therefore have the scheduler name
`saturn-stabilization-validation`; the outputs and receipts belong to this new
Rostami campaign. The default is two concurrent tasks, each requesting one A100,
eight CPUs, 32 GB and one hour on `gpu-highvram,uri-gpu`. The account uses Unity's
existing default association unless `SATURN_UNITY_ACCOUNT` is explicitly set.

```bash
bash slurm/unity/rostami.sh status configs/rostami_01 validation
bash slurm/unity/rostami.sh collect configs/rostami_01 validation
```

`status` uses the selected stage's frozen receipt, not a shell `$JOB` variable.
`collect` refuses incomplete, duplicate or configuration-mismatched outputs and
requires both validation success markers and the comparison report. It prints
the exact archive path and SHA256. Inspect `cpu_gpu_validation.toml`,
`CPU_GPU_VALIDATION_PASSED` and `COMPLETED` in all three outputs before proceeding.

`validate SUITE STAGE` is a read-only hash/provenance check. It does not execute
a model or run GPU validation. `prepare SUITE STAGE` optionally freezes inputs
without submitting; `submit` can reuse that preparation if nothing changed.

## Review the controls, then choose the next stage

The suite has 32 configurations. There is no command that submits all stages
or chains them automatically.

| Stage | Cases | Question |
|---|---:|---|
| `validation` | 3 | Does the new initializer agree on CPU and GPU? |
| `controls` | 3 | Do the unperturbed jet-only, reference-core and stronger-core equilibria remain circular numerically? |
| `formation` | 6 | What does each core do to a seeded sixfold disturbance versus competing modes 2:8? |
| `competition` | 3 | Does the reference case survive modes 1:32 with three independent phase seeds? |
| `qg_regime` | 3 | What changes when the upper cyclone speed is reduced to 10%, 25% or 50% of the paper reference? |
| `deep` | 4 | Can the deep jet, deep cyclone, or both change persistence and drift at fixed upper winds? |
| `deformation` | 2 | What changes at deformation radii 2,400 and 4,000 km with unchanged dimensional upper winds? |
| `core` | 3 | How much do inner-core regularization and spatial resolution affect the result? |
| `numerics` | 4 | Does the result depend on time step, physical viscosity or initial disturbance amplitude? |
| `boundary` | 1 | Does moving the exterior taper outward affect the result? |

After reviewing validation:

```bash
bash slurm/unity/rostami.sh submit configs/rostami_01 controls
bash slurm/unity/rostami.sh status configs/rostami_01 controls
# Once complete:
bash slurm/unity/rostami.sh collect configs/rostami_01 controls
```

Inspect the axisymmetric controls before submitting `formation`. An initially
circular flow is a continuum QG equilibrium; large changes here require
diagnosing numerical drift before interpreting perturbation growth.

The reference cyclone can have order-one vorticity/f near the pole even after
regularization. Inspect the profile and actual-grid audits before treating its
evolution as a controlled QG approximation. The `qg_regime` stage tests weaker
upper cores at fixed radius and regularization, alongside the jet-only and
reference cases. These are explicit departures from the paper; reduced
amplitude alone does not certify QG validity. Prioritize these comparisons
before making a physical claim about the reference core's stabilizing effect.

The formation and sensitivity cases last **3,125,000 seconds = 1000/f₀**. This is
the nondimensional duration illustrated in the paper, not 1000 inertial periods
and not a proof of long-term stability. Default `dt=25 s` gives 125,000 steps per
baseline case. Use validation/control timings to choose a wall time; four hours
is a starting request below, not a measured runtime guarantee.

The 50,000-second snapshot interval produces about 64 saved frames per full
run. Five Float64 fields per frame require roughly **2.5 GiB at 1024²** or
**10 GiB at 2048²**, before checkpoints and other output. Check available storage
before choosing later stages; these are estimates, not measured archive sizes.

```bash
# Only after reviewing the controls:
bash slurm/unity/rostami.sh submit configs/rostami_01 formation --time 04:00:00
bash slurm/unity/rostami.sh status configs/rostami_01 formation
# Once complete:
bash slurm/unity/rostami.sh collect configs/rostami_01 formation
```

Use the same pattern to submit each later stage only when it addresses a
remaining question. `--concurrency 1` reduces the number of simultaneous jobs.
The `numerics/nx2048_fixed_nu` case uses a 0.78125-second timestep to respect the
explicit hyperviscosity bound: **4,000,000 steps**. Budget it from measured
timings, separately from a casual one-hour request; `--time 12:00:00` is only a
possible initial request and is not a guarantee of completion. A time-limited
case is incomplete, not evidence about stability.

Every preparation freezes configurations, source/dependency/launch hashes,
campaign metadata and the scale audit. Changes after preparation are rejected.
Duplicate submissions and uncertain scheduler responses retain their reserved
attempt; inspect `submission.log` before considering `--new-attempt`. A new
attempt preserves prior outputs and does not cancel an existing job. Old
checkpoints are not resumed automatically. Bookkeeping remains in the original
runtime's `.unity/experiments-01/`, scoped by this suite's absolute path and stage;
simulation outputs are separate under `runs/rostami-01`.

## Interpret the comparisons carefully

- The smooth 500-km inner core is the baseline model. Compare 250-km and 500-km
  regularizations **at the same 2048² resolution**, then compare 500-km at 1024²
  and 2048². The 1000-km case is another core model, not a convergence target.
- `core/epsilon500km_nx2048` holds cutoff damping time fixed, reducing physical
  hyperviscosity by 256. `numerics/nx2048_fixed_nu` holds physical hyperviscosity
  fixed and uses the smaller explicit time step. Both are needed to separate
  grid and dissipation effects.
- The 80-Mm square uses planetary-PV and streamfunction tapers from 28 to 36 Mm.
  The paper used an exterior sponge; our boundary treatment is different. The
  100-Mm box retains grid spacing and physical hyperviscosity while moving the
  taper to 35–45 Mm. Inspect the zero-mean interface gauge as well as the jet.
- The 10-km disturbance setting is the coefficient-norm displacement parameter.
  Broadband modes share it with normalization 1/√N; at the jet the nominal RMS
  of δψ divided by jet speed is 10/√2 km. This is not a measured contour
  displacement. The exact perturbation differs from the paper's linear
  eigenmode. A mode-six seed does not establish spontaneous sixfold selection.
- Track mean upper winds, connected-jet coverage, competing azimuthal modes,
  phase drift, interface displacement and vorticity/f together. A persistent
  pattern with the wrong wind profile or rotation rate does not establish an
  observational fit. The core can remain marginally QG even when a numerical
  integration succeeds.

## Return results for analysis and gallery publication

On the Mac in Bash, copy the exact archive path printed by `collect`:

```bash
read -r -p 'Exact Unity collection archive path: ' SATURN_COLLECTION
download_dir=$(mktemp -d "$HOME/Downloads/saturn-rostami-01.XXXXXXXX")
scp -i "$HOME/Downloads/privkey.key" \
    "nicholas_pizzo_uri_edu@login.unityhpc.org:${SATURN_COLLECTION}" "$download_dir/"
shasum -a 256 "$download_dir"/*.tar.gz
```

Compare with Unity's checksum before extraction into a new directory. Keep raw
scientific fields and checkpoints outside the public repository. The gallery
does **not** pull Unity outputs automatically: collect, download, analyze/render,
review, and publish selected results through [the gallery workflow](gallery.md).

The radial linear driver and the profile audit run on a CPU workstation or an
allocated CPU job, not as heavy work on Unity's login node. Use the analysis
dependencies in `analysis-requirements.txt`. For example:

```bash
python3 scripts/stabilization_linear.py configs/rostami_01/formation \
    --output analysis/rostami-linear-01 --radial-nodes 128,256 --m-range 1:32
```

The circular inviscid eigenproblem and periodic nonlinear simulation have
different outer boundaries. Refine the radial grid and enclosing radius,
especially for weak growth or critical layers; do not label a branch stable
from a single quiet spectrum. After extracting nonlinear outputs, use
`scripts/summarize_experiments.py` and `scripts/field_mode_check.py` to cross-check
jet tracking against field modes before interpreting drift. The dedicated
[physics specification](rostami-benchmark.md) documents the audit command and
the balance and frame conventions.

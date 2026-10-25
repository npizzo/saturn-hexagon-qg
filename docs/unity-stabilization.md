# Staged stabilization experiments on Unity

This campaign tests whether a prescribed circular lower-layer flow changes the
stability and phase speed of the existing upper-layer jet. The new matched
profile includes the central cyclone. All initial stages are **unforced**.
The configurations are experimental proposals, not validated stable hexagons.

Use a separate checkout at
`/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-stabilization-01` and reuse
Julia 1.10.12 and the depot already installed under
`/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity`. Do not bootstrap or
install packages again. The original checkout, outputs, and campaign remain in
place. Do not resume an old checkpoint into the new source/configurations.

## Download an immutable checkout

From the Mac, connect using the existing key:

```bash
ssh -i "$HOME/Downloads/privkey.key" nicholas_pizzo_uri_edu@login.unityhpc.org
bash
```

On Unity, set `SATURN_COMMIT` to the **full 40-character commit ID supplied with
this campaign release**. A moving branch URL is not a reproducible deployment.
The following block downloads directly from the public GitHub repository; it
does not need `gh`, a GitHub token, or credentials on a compute node.

```bash
read -r -p 'Campaign commit (40 hexadecimal characters): ' SATURN_COMMIT
export SATURN_COMMIT
bash <<'BASH'
set -euo pipefail
[[ ${SATURN_COMMIT:-} =~ ^[0-9a-f]{40}$ ]] || { echo 'A full release commit is required.' >&2; exit 1; }
base=/work/pi_nicholas_pizzo_uri_edu
project="$base/saturn-hexagon-qg-stabilization-01"
[[ ! -e $project ]] || { echo "Refusing to overwrite $project" >&2; exit 1; }
[[ -x $base/saturn-hexagon-qg/.unity/runtime/julia-1.10.12/bin/julia ]] || { echo 'Existing Julia runtime is missing.' >&2; exit 1; }
[[ -d $base/saturn-hexagon-qg/.unity/julia-depot ]] || { echo 'Existing Julia depot is missing.' >&2; exit 1; }
download=$(mktemp -d "$base/saturn-stabilization-download.XXXXXXXX")
trap 'rm -rf -- "$download"' EXIT
curl --fail --location --proto '=https' --tlsv1.2 \
    "https://github.com/npizzo/saturn-hexagon-qg/archive/${SATURN_COMMIT}.tar.gz" \
    --output "$download/source.tar.gz"
mkdir -- "$project"
tar -xzf "$download/source.tar.gz" --strip-components=1 -C "$project"
printf '%s\n' "$SATURN_COMMIT" > "$project/DEPLOYED_COMMIT"
cd -- "$project"
bash slurm/unity/stabilization.sh validate configs/stabilization_01 validation
BASH
```

If the block reports an error, stop and inspect that error; it submits no jobs.
An interrupted extraction can leave the new directory incomplete, and rerunning
the block intentionally refuses to overwrite it. A successful deployment ends
with the provenance/configuration hash validation message.

## Submit only the short validation stage first

In the same Unity login shell:

```bash
cd /work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-stabilization-01
unset SATURN_RESUME
export SATURN_STABILIZATION_RUNTIME_ROOT=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity
export SATURN_STABILIZATION_RUN_ROOT=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg-stabilization-01/runs/stabilization-01
bash slurm/unity/stabilization.sh submit configs/stabilization_01 validation
```

This submits one array with three short CPU/GPU comparisons using the
`matched_upper` profile: a barotropic reference, a half-matched deep flow at
deformation radius 3.2 Mm, and a fully matched flow at deformation radius 0.5 Mm.
At most two tasks run simultaneously. Each task requests one A100 GPU, eight
CPUs, 32 GB, and one hour on the existing `gpu-highvram,uri-gpu` partitions.
The account is Unity's existing default association, unless the user explicitly
exports `SATURN_UNITY_ACCOUNT`. Inherited `SBATCH_*` options are cleared in the
submission child. Generic old Unity/Delta runtime and run-directory exports are
overridden by the dedicated stabilization settings above.

`validate SUITE STAGE` is a read-only hash check. It does **not** run the
CPU/GPU comparison; `submit SUITE validation` does. `prepare SUITE STAGE`
optionally freezes a stage without submitting it, and a subsequent `submit`
reuses that frozen preparation if nothing changed.

Check the exact recorded stage, without a reusable shell job-ID variable:

```bash
bash slurm/unity/stabilization.sh status configs/stabilization_01 validation
bash slurm/unity/stabilization.sh collect configs/stabilization_01 validation
```

`status` is read-only and reports the named receipt, each member's output, and
`squeue`/`sacct` state. `collect` refuses incomplete, duplicated, or mismatched
outputs. Once all members finish, it creates a new archive and prints its exact
path and SHA256. Inspect all three `cpu_gpu_validation.toml` reports and both
`CPU_GPU_VALIDATION_PASSED` and ordinary `COMPLETED` markers. The validator exits
nonzero on failed comparisons; failed comparisons must not be treated as
successful scientific runs. A passed short comparison does not establish
long-time stability, convergence, or observational agreement.

## Review each stage before selecting the next

The generated suite contains 39 configurations in these explicit stages:

| Stage | Cases | Purpose |
|---|---:|---|
| `validation` | 3 | Short CPU/GPU and output checks for new lower-layer profiles |
| `controls` | 5 | Four axisymmetric controls plus a perturbed barotropic reference |
| `screen` | 9 | Lower-layer matching fraction versus deformation radius |
| `transition` | 4 | Resolve the small-deformation-radius transition |
| `width` | 6 | Test sensitivity to jet width within the stated profile constraints |
| `boundary` | 2 | Enlarge the box/taper region at unchanged grid spacing and physical viscosity |
| `numerics` | 10 | Selected time-step, resolution, perturbation, and dissipation checks |

After reviewing validation, submit `controls` explicitly; after reviewing the
controls, select `screen`. The same command pattern applies separately to
`transition`, `width`, `boundary`, and `numerics` when their results will address a specific
remaining question. No command submits all stages or automatically chains them.

```bash
bash slurm/unity/stabilization.sh submit configs/stabilization_01 controls
# Later, after reviewing controls:
bash slurm/unity/stabilization.sh submit configs/stabilization_01 screen
# Later, if checking sensitivity to the outer domain:
bash slurm/unity/stabilization.sh submit configs/stabilization_01 boundary
```

Use `status` and `collect` with the same stage name each time. The default is
two concurrent tasks; `--concurrency 1` lowers that cap. `--time HH:MM:SS`
changes wall time. The fixed-viscosity 2048² cases in `numerics` require much
smaller explicit time steps and are expensive: use measured earlier timings to
budget them, with `--time 12:00:00` as a possible starting request, not a runtime
guarantee. Do not submit this stage blindly with the one-hour default.

Every submission freezes its configurations, launch/source/dependency hashes,
campaign notes, and physics audit. Frozen receipts are stored below the original
checkout's `.unity/experiments-01/` bookkeeping tree, scoped by the
new suite's absolute path and stage; outputs live in the new run directory.
Each array task verifies provenance again before starting Julia. Changes after
preparation cause a rejection. Duplicate submission is also rejected, including
an uncertain scheduler response: inspect `submission.log` before deciding on
an explicit `--new-attempt`. A new attempt preserves old attempts and does not
cancel any existing job. `status` and `collect` refer to the latest prepared
attempt for the selected suite/stage, and remain usable after source edits.

## Bring a completed stage back to the Mac

Copy the exact archive path printed by `collect`. On the Mac, in Bash:

```bash
read -r -p 'Exact Unity collection archive path: ' SATURN_COLLECTION
download_dir=$(mktemp -d "$HOME/Downloads/saturn-stabilization-01.XXXXXXXX")
scp -i "$HOME/Downloads/privkey.key" \
    "nicholas_pizzo_uri_edu@login.unityhpc.org:${SATURN_COLLECTION}" "$download_dir/"
shasum -a 256 "$download_dir"/*.tar.gz
```

Compare the printed SHA256 with Unity's collection checksum before analysis.
The archive includes frozen inputs, the physics audit, submission receipt,
available logs, and each member's output. Keep raw scientific outputs outside
the public Git repository. The separate gallery workflow renders selected
results locally for review before any explicit publication.

## Linear screening and analysis

The radial linear driver reads the same frozen configurations and verifies
their source and configuration hashes. Run it in this checkout on a CPU
workstation or a separately allocated CPU job, not as a heavy calculation on
Unity's login node. It needs NumPy, SciPy, and TOML support as listed in
`analysis-requirements.txt`. For a new workstation environment:

```bash
python3 -m venv .analysis-env
.analysis-env/bin/python -m pip install -r analysis-requirements.txt
.analysis-env/bin/python scripts/stabilization_linear.py configs/stabilization_01/screen \
    --output analysis/linear-screen-01 --radial-nodes 128,256 --m-range 1:32
```

This is a circular, inviscid eigenproblem with a circular outer boundary;
the nonlinear simulations use the periodic Cartesian domain. Compare growing
modes across radial grids and repeat with a different enclosing outer radius
using `--outer-radius` and a new output directory. For the small-LD transition,
use `configs/stabilization_01/transition` and finer radial grids. Weak growth
and critical layers need stronger convergence checks; a quiet spectrum at
one resolution is not a stable branch. Extend the mode band if growth remains
near its upper edge. The driver retains all eigenvalues and growing
eigenfunctions and does not invent a drift from a neutral spectrum.

After extracting a collected nonlinear stage into a **new** local directory,
its result directories are under `runs/`. For example, if the extraction
directory is `analysis/returned-screen-01`:

```bash
.analysis-env/bin/python scripts/summarize_experiments.py \
    analysis/returned-screen-01/runs/*/results --output analysis/screen-summary-01
```

The output records provisional growth and frame-aware phase estimates,
connected-jet coverage, and extrema of the saved thickness and Rossby-number
indicators. It does not label runs stable. Inspect snapshots and cross-check
jet-tracker modes against the field diagnostic before interpreting drift.
See [gallery.md](gallery.md) for producing a gallery from selected actual
outputs. The current public gallery remains the record of earlier runs until
new results have been collected and reviewed.

To prepare long follow-ups after identifying a candidate, generate a **new**
campaign, using an ID from `scripts/make_stabilization.jl --list-cases`:

```bash
export JULIA_BIN=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity/runtime/julia-1.10.12/bin/julia
export JULIA_DEPOT_PATH=/work/pi_nicholas_pizzo_uri_edu/saturn-hexagon-qg/.unity/julia-depot
"$JULIA_BIN" --project=. scripts/make_stabilization.jl configs/pilot.toml \
    configs/stabilization_candidates_01 --candidates ld500_a100_w1000
```

The ID above is a syntax example, **not a selected or stable candidate**.
This adds a 100-km mode-six persistence test and three 100-km broadband
competition tests, each lasting at least five million seconds, plus numerical
checks for the chosen candidate. Submit `maintenance` and `competition`
explicitly only after reviewing the candidate's growth, physical regime and
resolution. They still start from the generic perturbed circular jet, not a
reconstructed published nonlinear solution. No convection cases are generated.

## Scientific scope

This series studies circular base-flow stabilization and the nonlinear response
to specified perturbations. It does not initialize a reconstructed published
Pizzo–Salmon (2024) equilibrium or claim to reproduce Rostami et al. (2017).
The full ambient PV and rotation/frame information for the 2024 solution are
not yet available as a validated model initial state. Those routes require
reconstruction, equilibrium residual checks, and finite-deformation-radius
continuation before deployment; an RSW result also requires a checked balanced
QG counterpart. Any later candidate branch is generated and reviewed
separately. Convection calibration comes after the unforced benchmarks.

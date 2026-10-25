# Lower-layer experiments for the current Saturn QG model

This campaign asks whether finite deformation radius and prescribed circular
deep flow suppress disruptive disturbances, permit a finite hexagonal jet, and
give it an observationally acceptable drift. These are separate requirements.
The supplied MATLAB initializer is documented in [the legacy IC audit](legacy-ic-audit.md).
It helps identify the 2024 PV formula, but the experiments below use the existing,
tested annular-jet initializer. They do not claim to import the exact family.

The default prepared suite is `configs/experiments_01`: **six initial runs**,
eleven additional physics cases, and eight numerical sensitivity cases.
Only one named stage is submitted at a time. There is no automatic progression
to the next stage and no simulation was launched while preparing this suite.

## Physical controls and initial six

Hold the upper streamfunction fixed across each physics comparison: annular
jet speed coefficient 100 m/s, nominal radius 15 Mm, width parameter 1 Mm;
central cyclone speed coefficient 50 m/s and radius 3 Mm. Keep the 60-Mm box,
planetary PV continuation, zero sponge, rotation parameter and perturbation
realization fixed. Deep flow is a prescribed annular jet at the same radius
and width, with no prescribed deep central cyclone. Thus deep/upper ratio 1
matches the annular component, not the entire upper flow.

The initial `screen` contains:

| Case ID | LD (km) | Deep jet coefficient (m/s) | Purpose |
|---|---:|---:|---|
| `bt` | infinite | 0 | Barotropic control |
| `ld2000_d000` | 2000 | 0 | Finite deformation over resting deep flow |
| `ld2000_dp050` | 2000 | +50 | Current pilot's deep-flow parameters |
| `ld4000_d000` | 4000 | 0 | Weaker coupling over resting deep flow |
| `ld4000_dp050` | 4000 | +50 | Weaker coupling with coflow |
| `ld4000_dp100` | 4000 | +100 | Deep annular flow matches the upper annular component |

All six are **unforced**. They start with broadband modes 1–32 and seed 42.
The configured displacement is 10 km; the initializer divides each cosine's
coefficient by sqrt(32), so each begins near 1.77 km and the combined radial
RMS is approximately 7.07 km at the nominal ring. The actual Eulerian field and
measured centreline, not this first-order estimate, determine the initial modes.
Identical upper velocity does not imply identical total PV when LD/deep flow
changes; this is intentional and is the controlled mechanism comparison.

Use 1024², dt=50 s, duration=500,000 s, diagnostics every 2,000 s, fields every
50,000 s, and a checkpoint every 100,000 s. Record modes through 32 using
720 angles and 256 radial samples. The nominal jet circuit is 942,478 s:
this first screen lasts 5.79 Earth days, only 0.53 circuits. It can identify
rapid growth, not certify stability of a quiet case. Small-LD weak eigenvalues
are especially resolution-sensitive and may require much longer integrations.
Budget approximately 0.5 GiB per 1024² screen result before compression, or
roughly 3 GiB for the first six. A 2048² numerical-check result is approximately
2 GiB. These estimates cover the planned fields/domain/checkpoint and are not
storage quotas or runtime measurements.

The already-completed circular eigenvalue analysis motivates these brackets;
its disk boundary differs from the nonlinear periodic square. It predicts
linear circular-base modes, not a stable finite-amplitude hexagon.

## Physical-regime audit

For a standard reduced-gravity interpretation, the active-layer thickness
anomaly is h'/H=(psi-psi_d)/(f0 LD²), with f0=2 Omega. Upward displacement of
the lower interface has the opposite sign. This follows from the reduced-gravity
pressure relation in the [GFDL derivation, equation 3](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/shallow.pdf)
and the QG thickness term in [Thiry et al. (2024), equations 7 and 10](https://doi.org/10.1029/2024MS004510).
`physics-audit.tsv` reports this
quantity for the axisymmetric base, including the taper and the model's
full-square zero-mean convention for both streamfunctions. It excludes the
initial perturbation, spectral projection and subsequent evolution; it is a
base-state smallness indicator, not an online layer-thickness measurement.

At LD=500 km, max|h'/H| is 3.92–5.87 for deep speeds +50 to -50 m/s, and the
inferred base thickness becomes negative. Even matching the annular flow
leaves the unmatched central cyclone. Those runs are **mathematical mechanism
and stress tests**, not admissible thin-displacement Saturn candidates under
these profiles. They are retained in `expanded`, not the initial six. The
audit also flags problematic 1000-km cases. A stable result there would motivate
revising the physically constrained base/deep profiles or testing a model that
allows finite thickness changes; it would not establish QG stabilization of
Saturn's hexagon. Changing a streamfunction constant while holding the same
H and LD is not a physical repair of this failure.

For LD=2000 km and deep speed 0 or +50 m/s, the maximum base |h'/H| is about
0.306 or 0.245: these are already marginal for an asymptotically small
displacement. At 4000 km, the corresponding values are 0.077 and 0.061;
matching the annular flow gives about 0.046. The nominal jet-width Rossby
number U/(f0 w) is 0.303 throughout, so even the better thickness regime does
not by itself validate QG for Saturn. Evaluate actual local shear, curvature,
and the evolving interface as well as comparing against goflow. The audit
does not assign an arbitrary binary physical-validity threshold to these
finite values.

## Numerical checks before interpreting a quiet case

The `numerics` stage uses two anchors: `ld2000_dp050` and `ld4000_dp050`.
For each, compare four variants with its existing screen run:

| Variant | Changed setting | What it checks |
|---|---|---|
| `dt25` | dt=25 s | Time integration |
| `nx2048_dt25` | 2048², dt=25 s | Mesh plus the specified cutoff-damping prescription |
| `amplitude5km` | Broadband displacement 5 km | Initial linearity and diagnostic floor |
| `weaker_damping` | Cutoff damping time 1,000,000 s | Tenfold reduction in physical hyperviscosity at fixed grid |

The grid refinement keeps the damping time at the numerical cutoff fixed.
Since the model uses fourth-order hyperviscosity in the Laplacian, physical
nu decreases by approximately 256 when the resolution doubles. This is
**not fixed-nu convergence**. The separate damping test helps diagnose
sensitivity but does not fully separate the two effects. Compare normalized
mode growth over overlapping small-amplitude ranges, not chaotic final fields
or equal final displacements. If a candidate exists only under strong damping,
or its growth/phase changes materially under these checks, refine it before
making a physical claim.
An explicitly promoted candidate outside these two anchors receives its own
four numerical variants when generating the candidate suite.

The `expanded` stage contains the remaining eleven cases in LD={500,1000,2000,
4000} km × deep/upper={-0.5,0,0.5,1}, plus the single barotropic control already
in `screen`. It uses the same small broadband initialization. Run it after
reviewing the initial six; do not interpret the grid as a fitted Saturn range.

## Finite-amplitude hexagons, competition, then convection

The generator adds these stages **only for explicitly named candidate IDs**:

| Stage | Cases per candidate | Initial condition and question |
|---|---:|---|
| `maintenance` | 1 | Mode 6, 100-km displacement: does a seeded hexagon persist with the right drift? |
| `selection` | 3 | Broadband 100-km displacement, seeds 42–44: what modes survive nonlinear competition? |
| `convection` | 6 new | Same fresh mode-6 field: three unbiased and three biased forcing realizations, reusing `maintenance` as their unforced control |

The default duration is 2,000,000 s (23.15 days, 2.12 nominal circuits), with
2,000-s diagnostics and 200,000-s fields/checkpoints. Promising cases should
then be extended to at least 5,000,000 s (5.31 circuits), with longer times
set by measured growth and drift uncertainty. No finite duration proves
asymptotic stability. The generator's `--long-duration-s` prepares fresh-run
configurations; an actual checkpoint continuation must instead use the same
original config with only duration changed, as required by the existing restart
contract. Never change forcing, seed or initial amplitude while resuming.

Pure mode-six initialization restricts the excited symmetries. Its persistence
is not a nonlinear stability proof against arbitrary perturbations. Repeating
an unforced single-mode case with different seeds would produce duplicates;
the campaign deliberately does not do that. Broadband cases provide a distinct
competition test, but a strict stability test of an already formed hexagon
will need independently perturbed copies of that same state. This is a later
extension, not something the current configuration silently implements.

Convection starts from the same fresh upper field as the unforced control.
It tests development/persistence under forcing, not rescue of an equilibrated
state. The maintenance trajectory is reused as the shared unforced baseline;
it is not rerun or counted as independent evidence. Identical forcing seeds pair attempted event schedules and random
choices, but state-dependent coldness/PV limits can change actual accepted
injections and energy input. Compare event counts, caps, net energy/enstrophy
jumps and mean-flow changes. At 1024² the 250-km forcing radius spans only
4.27 cells, so any positive result needs forcing-resolution checks.

## Phase speed and criteria for advancing

Use angular pattern speed, the observed epoch and a stated longitude sign.
The production rotation parameter is 1.65e-4 s^-1 (period 10 h 34 m 39.911 s).
If it defines the model frame and +theta is eastward, the observed 2008–2014
period 10 h 39 m 23.01 s corresponds to c=-1.21761412414e-6 s^-1, or -18.2642
m/s at the model's 15-Mm radius. System III's period is 10 h 39 m 22.4 s.
These conditional values are recorded in `campaign.toml`; they are not a
replacement for registering the eventual goflow data and uncertainty.

Changing frame shifts both U/r and c, leaving their difference unchanged.
Changing `planet.rotation_rate_s` in the simulation changes the PV dynamics;
it cannot be used to cosmetically correct a phase-speed mismatch. The
observational sources and limitations are in [the experiment design](experiment-design.md#joint-requirement-a-stable-hexagon-with-the-observed-phase-speed).

Review each candidate against all of the following:

1. The base and evolved flow remain within a defensible physical regime for QG;
   review `physics-audit.tsv` before interpreting a quiet small-LD result.
   Saved fields remain finite and diagnostics/inversion/budgets are credible.
   `COMPLETED` is a numerical completion marker, not a scientific result.
2. Perturbation growth is assessed across all measured modes, including modes
   that compete with six. A quiet m=6 is insufficient if another mode grows.
3. A finite sixfold jet has a connected, unambiguous centreline, appreciable
   deformation and coherent phase. Roll-up into a vortex street or peak
   switching is not counted as maintained jet structure. Inspect snapshots
   and use the existing cubic jet analysis, including its ambiguity flags.
4. Pattern drift and upper velocity/jet width agree with the same observation
   frame within data and numerical uncertainties. A stable state with a large
   drift mismatch remains a mechanism case, not an observational candidate.
5. Apparent persistence survives the numerical checks and appropriate competing
   disturbances. Only then add convection or extend the selection experiment.

`scripts/summarize_experiments.py` reports preliminary CSV-based growth and
phase fits. Its default window is 20,000–200,000 s; it uses one contiguous
segment with at least ten samples, complete angular coverage, amplitude above
500 m and below 5% of jet width, combined measured displacement RMS below 5%
of jet width, and mean peak speed within 5% of its initial
value. It does not interpolate across gaps, assign zero to missing phase, or
label a case stable. Report actual fit windows and sensitivity to these
explicit thresholds. Runtime centreline fits remain provisional until field
and branch-continuity checks agree. A phase-step check cannot recover temporal
aliasing that has already occurred; cadence must resolve the physical mode.

For a cross-check independent of locating a velocity maximum,
`scripts/field_mode_check.py` samples saved streamfunction and tangential
velocity on the annulus Rj ± 2w and exports radial norms of their angular
Fourier coefficients. Its `2|psi_m(Rj)|/Uj` quantity is an amplitude proxy,
not a measured centreline displacement. The mode-six projection phase is
saved wrapped relative to the first snapshot; the script does not infer
speed from the sparser fields. This check recovered 99,999.983 m from the
archived 100-km seeded initial state. The weak broadband runs especially
need this check or the existing cubic jet extraction before interpreting
small centreline changes.

After downloading and extracting one completed stage into a new local folder,
pass its `runs/array-*/results` directories to the summary tool. For one result:

```sh
python scripts/summarize_experiments.py /absolute/path/to/results \
  --output analysis/NEW-summary
python scripts/field_mode_check.py /absolute/path/to/results analysis/NEW-fields
```

The field cross-check additionally needs SciPy and h5py. On the original Mac,
`/usr/local/anaconda3/bin/python` has the required analysis packages. These
tools read saved data and do not advance a model.

## Unity: submit the initial six

Clone this repository into allocated work storage, or update an existing copy
with a clean Git working tree. Install and validate the pinned environment as
shown in [the Unity instructions](unity.md). From the project root on Unity:

```sh
bash slurm/unity/experiments_01.sh submit configs/experiments_01 screen
```

This submits six array tasks with at most two running concurrently, using the
existing Unity GPU defaults and a one-hour wall limit per task. It retains
your existing default Unity Slurm association; no new account name is inferred
from the work-directory name. It freezes the
stage's configs and verifies model/dependency provenance before submission and
again in each task. A mismatched existing installation is refused. It does not
instantiate packages or run science on the login node. Source and dependencies
must match the already-validated Unity installation.

```sh
bash slurm/unity/experiments_01.sh status configs/experiments_01 screen
```

After all six finish:

```sh
bash slurm/unity/experiments_01.sh collect configs/experiments_01 screen
```

The wrapper prints the result archive path. Download that path to the Mac for
analysis. It collects only the recorded array's outputs, verifies one completed
result per member, and refuses overwriting an existing archive. A repeat
submission requires `--new-attempt`, which creates a separate receipt and
outputs. Use that only for a deliberate rerun, not as a status check.

After reviewing the screen, the next numerical checks use:

```sh
bash slurm/unity/experiments_01.sh submit configs/experiments_01 numerics
```

Submit `expanded` separately only if its additional parameter coverage is
useful. Candidate stages must first be generated after selecting IDs from the
actual results. For example, **if** `ld4000_dp050` warrants further testing:

```sh
source slurm/unity/environment.sh
"$JULIA_BIN" --project=. scripts/make_experiments.jl \
  configs/unity_followup_01/maintenance_unforced.toml configs/experiments_02 \
  --candidates ld4000_dp050
bash slurm/unity/experiments_01.sh submit configs/experiments_02 maintenance
```

This example is a conditional recipe, not a claim that this case is a good
hexagon. Submit its `selection` and, later, `convection` stages independently
after the relevant review. A longer forcing run may need `--time 02:00:00`;
measure actual elapsed time from the first stage before estimating cost.

## Reproduce preparation and tests

With Julia 1.10.12 installed, from the project directory:

```sh
julia --project=. scripts/make_experiments.jl \
  configs/unity_followup_01/maintenance_unforced.toml configs/experiments_check
julia --project=. test/experiment_tests.jl
python scripts/summarize_experiments.py --help
```

Use a new output directory. The generator uses only Julia standard libraries
and the shared configuration validator. The summary script needs NumPy and,
on Python before 3.11, `toml`; the Unity submission helper uses only Python 3
standard libraries. `cases.tsv` lists every experiment and its role;
`members.tsv` binds each config by hash; `provenance.tsv` binds the production
source, pinned environment and launch/preparation code. Protocol tests and
mocked scheduler tests accompany the update.

The separate `test/experiment_initialization_check.jl` builds each of the six
actual prepared 1024² inputs and advances one 50-s CPU step. It checks initial
upper-velocity agreement across physics controls, finite fields, angular
diagnostic coverage and the model's timestep guards. It is a short integration
check, not a run of the 500,000-s scientific protocol.

Validation during preparation: the configuration/provenance tests, analysis
tests and mocked scheduler tests passed. The expensive full-grid local check
completed the barotropic case, LD=2000 km with deep speed 0 and +50 m/s, and
the LD=500 km/-50 m/s stress case before it was deliberately interrupted during
the next initialization. The three final LD=4000-km screen inputs have not
been advanced locally; their configuration checks passed. The partial check is not a completed six-case integration test. No new Unity science
runs have been submitted, and this update leaves the production solver and
pinned dependencies unchanged.

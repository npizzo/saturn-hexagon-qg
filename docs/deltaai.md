# DeltaAI deployment and verification

These scripts prepare one simulation per GPU and distribute experiments across a Slurm job array. They have not been run on DeltaAI. The first allocation must verify the pinned software environment on the actual hardware before production use.

DeltaAI combines four GH200 superchips per node. Each combines a Grace **Linux ARM64/aarch64** CPU and an H100 GPU with 96 GB GPU memory. Its production partition is `ghx4`. The templates request one GPU, one task, eight CPU cores, 32 GB host memory, and one hour; benchmark before increasing these requests. GeophysicalFlows/FourierFlows currently uses a single GPU per simulation. Additional GPUs should run independent parameter/seed combinations, rather than being allocated to an unchanged single simulation. [NCSA architecture](https://docs.ncsa.illinois.edu/systems/deltaai/en/latest/user-guide/architecture.html), [NCSA Slurm guide](https://docs.ncsa.illinois.edu/systems/deltaai/en/latest/user-guide/running-jobs.html), [GeophysicalFlows GPU support](https://fourierflows.github.io/GeophysicalFlowsDocumentation/dev/gpu/)

## Prepare a project copy

Copy this complete source directory, including `Project.toml` and `Manifest.toml`, to your project storage. Do not copy a Mac or x86 Julia executable, compiled package cache, or system image. Source files and the pinned manifest travel; platform-specific binaries are installed on DeltaAI.

The supplied bootstrap selects the official Julia **1.10.12 Linux aarch64** archive on DeltaAI (or Linux x86_64 on Unity), verifies its SHA256, and installs only into directories you explicitly provide. It does not modify shell profiles or the global Julia depot. The runtime is intentionally pinned rather than claiming to be the newest Julia release. [Julia manual downloads](https://www.julialang.org/downloads/manual-downloads/)

Choose actual account/user paths, then run from the project copy:

```bash
bash slurm/bootstrap.sh /absolute/runtime-root /absolute/julia-depot
```

Set the three `export` lines it prints. Set a run directory on your work allocation as well:

```bash
export SATURN_PROJECT=/absolute/path/to/saturn-hexagon-qg
export JULIA_BIN=/absolute/runtime-root/julia-1.10.12/bin/julia
export JULIA_DEPOT_PATH=/absolute/julia-depot
export SATURN_RUN_ROOT=/work/nvme/ACCOUNT/USERNAME/saturn-runs
```

Replace `ACCOUNT` and `USERNAME`; no account is inferred by the code. DeltaAI provides `/work/hdd` and `/work/nvme`, rather than `/scratch`. Use project storage for source and curated results, work storage for runs, and follow your allocation's retention/quota policy. The `/projects` and `/work` areas have no file recovery snapshots; archive important results separately. [NCSA data management](https://docs.ncsa.illinois.edu/systems/deltaai/en/latest/user-guide/data-mgmt.html)

The script deliberately requires one explicit `JULIA_DEPOT_PATH` directory. Populate it before launching an array, avoiding dozens of concurrent package downloads/precompiles. Existing module settings and Slurm's GPU visibility are preserved. The scripts neither purge your modules nor guess a CUDA module name.

## Verification in a short allocation

The simplest first job is the supplied validation batch script. After bootstrap
has succeeded and all four environment variables above are exported, run this
**on DeltaAI**, from the copied project directory:

```bash
sbatch -A ACCOUNT slurm/validate_gpu.sbatch
```

Replace `ACCOUNT` with your real Slurm allocation. This requests one GPU for
30 minutes and runs the CPU test suite on the ARM host, records the CUDA
environment, compares short CPU/GPU trajectories including convection, and
executes a GPU smoke run with saved fields and a checkpoint. It stops at the first
failure. Its final `VALIDATION PASSED` line confirms these checks; it does not
certify GPU restart equivalence, production-resolution performance, or physical
hexagon stability. Complete package installation/precompilation before submission.
The script itself has only been checked locally; actual DeltaAI execution remains
to be verified. Unset `SATURN_RESUME` for this first job.

The following interactive procedure is an alternative when debugging an error.

Use your real allocation from the `accounts` command. For example:

```bash
srun -A ACCOUNT --partition=ghx4 --nodes=1 --ntasks=1 \
  --gpus-per-node=1 --cpus-per-task=8 --mem=32G \
  --time=00:30:00 --pty bash
```

Inside the allocation, first test the CPU implementation and then CUDA availability:

```bash
export JULIA_NUM_THREADS=8
"$JULIA_BIN" --startup-file=no --project="$SATURN_PROJECT" \
  -e 'using Pkg; Pkg.test()'
"$JULIA_BIN" --startup-file=no --project="$SATURN_PROJECT" \
  -e 'using CUDA; CUDA.versioninfo(); CUDA.functional(true) || error("CUDA unavailable"); CUDA.allowscalar(false); a = CUDA.ones(Float64, 16); @assert sum(a) == 16; CUDA.synchronize()'
```

Then run the small smoke configuration on CPU and GPU into separate new output directories:

```bash
"$JULIA_BIN" --startup-file=no --project="$SATURN_PROJECT" \
  "$SATURN_PROJECT/scripts/run.jl" "$SATURN_PROJECT/configs/smoke.toml" \
  --output "$SATURN_RUN_ROOT/smoke-cpu-01" --device cpu
"$JULIA_BIN" --startup-file=no --project="$SATURN_PROJECT" \
  "$SATURN_PROJECT/scripts/run.jl" "$SATURN_PROJECT/configs/smoke.toml" \
  --output "$SATURN_RUN_ROOT/smoke-gpu-01" --device gpu
```

Verify finite fields, reported device, step counts, checkpoint reading, and short-time agreement of energy and mode diagnostics within numerical tolerance. Long chaotic trajectories need statistical comparisons, not pointwise equality. Local CPU tests cannot certify GPU execution or ARM64 package compatibility.

CUDA.jl supports Linux ARM64/SBSA and Hopper, but the installed driver, pinned package version, toolkit artifacts, and Julia version must work together. Do not upgrade the manifest during a production ensemble. If toolkit discovery fails, record `CUDA.versioninfo()` and the full error; consult NCSA's available modules and CUDA.jl's installation guide. A site-local toolkit can be selected through CUDA.jl preferences when necessary, but record the resulting `LocalPreferences.toml` with the experiment. This is an explicit environment change, not an automatic fallback. [CUDA.jl installation](https://cuda.juliagpu.org/stable/installation/overview/), [CUDA.jl support matrix](https://github.com/JuliaGPU/CUDA.jl)

If package precompilation on a login node did not select GPU artifacts, finish GPU initialization once in this allocation before an array. The runner must fail when GPU execution is requested and unavailable; a CPU run should never masquerade as a GPU run.

## Single simulation and restart

From the project directory:

```bash
sbatch -A ACCOUNT slurm/single_gpu.sbatch configs/smoke.toml
```

Use the smoke configuration for deployment testing only. The batch command accepts a configuration path relative to `SATURN_PROJECT` or an absolute path. Its unique output parent is printed in the job log; a copy of the exact configuration, manifest hash, and Slurm context sit alongside the runner's `results` directory. Paths must be accessible from compute nodes. Slurm writes stdout/stderr in the submission directory unless you override `--output` and `--error`.

To resume one compatible run, submit its configuration with an absolute checkpoint path:

```bash
SATURN_RESUME=/absolute/previous/results/checkpoint.jld2 \
  sbatch -A ACCOUNT slurm/single_gpu.sbatch /absolute/compatible-config.toml
```

Use the checkpoint filename actually emitted by your run. The resumed run gets a new output directory and must preserve the scientific parameters and RNG state required by the runner's restart checks. This template does not request automatic requeue and does not promise a last-second checkpoint on timeout. Set normal checkpoint intervals well below wall time, and split long integrations into planned segments. Do not overwrite the original results.

## Prepared ensembles

Prepare scientifically reviewed configurations before queue submission. Freeze their order and hashes with:

```bash
bash slurm/prepare_array.sh /absolute/experiment-01.tsv \
  /absolute/config-a.toml /absolute/config-b.toml
```

Set the two exports printed by the command, then use the printed array length, for example:

```bash
sbatch -A ACCOUNT --array=1-2%2 slurm/ensemble.sbatch
```

The `%2` bounds simultaneous tasks. Select concurrency using your allocation rather than filling the queue by default. The job validates the manifest hash, array index, and selected configuration hash before starting. Changing a configuration after preparation causes a clear failure; create a new manifest for a revised experiment. Resume individual array members with the single-job script, not a global shared resume path.

Each launch reserves a unique directory atomically, including re-submissions. It stages the config through a temporary name before publication. The runner publishes checkpoints by renaming a `.partial` file; an interrupted `.partial` file is not a restart source. Record source revision or preserve the complete source snapshot as well as the manifest: configuration hashes do not make later source edits reproducible.

## Before increasing scale

1. Measure step time after initial compilation at representative resolution, with actual diagnostics and checkpoint cadence.
2. Check peak GPU memory and host memory, including FFT workspaces and output transfers.
3. Confirm CPU/GPU short-time agreement and continuous-versus-restarted agreement before interpreting long integrations.
4. Set wall time and ensemble concurrency from the measurements. A small problem can be dominated by compilation or data movement; no H100 speedup is assumed.

Large arrays and long integrations remain a deployment step after these checks; none are submitted by this repository.

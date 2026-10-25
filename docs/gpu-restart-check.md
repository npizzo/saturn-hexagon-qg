# Checkpoint continuation gate

**Passed on Unity in job 66043703.** The downloaded report in
`runs/from-unity/restart-check-66043703.SEWZhnZ9/results/report.toml` records
zero final field/spectral differences, with identical RNG and event schedules,
through 20,000 seconds. This was a fresh-model restoration on the same GPU in
one allocation. It is separate from spatial convergence or long-time stability.

This check compares an uninterrupted trajectory with a fresh model loaded from
a serialized checkpoint on the same device, in the same allocation. It uses the
supplied configuration unchanged, including convection and its random seed.
For the current Unity follow-up, use `configs/unity_benchmark.toml`: 1024²,
100 s nominal steps, and 20,000 s duration. The successful timestep comparison
does not itself establish checkpoint equivalence.

From the existing Unity project, after copying the updated validator and launcher
files, submit:

```bash
unset SATURN_RESUME
bash slurm/unity/submit.sh restart-check configs/unity_benchmark.toml
```

This reuses the installed Julia runtime and depot. It does not run the installer,
the complete test suite, or a separate CPU trajectory. The launcher retains the
usual Unity account selection (site default unless explicitly overridden), one
A100 GPU, eight CPUs, and 32 GB. The default wall limit is one hour; ordinary
Slurm options may follow the configuration path. As with other runs, the job ID
is recorded in `.unity/last-job-id` (under `SATURN_UNITY_ROOT` if customized).

The checkpoint is saved at the first natural integrator/event boundary at or
after half the duration. No extra shortened step is introduced at the split.
There must be at least one nonzero convective injection on each side. The check
compares:

- Restored time, step, event counter, next-event time, convection settings, and
  32 draws from a copy of the random-number generator.
- Exact spectral checkpoint serialization, followed by agreement of the restored
  retained spectral modes. Restore applies the library's usual dealias projection;
  the uninterrupted branch applies the same projection before its next step.
- Every subsequent step/event boundary, random event location, and injection flag.
- Final random stream and spectral and physical fields (`s`, `psi`, `u`, `v`,
  and total PV), with relative L2 differences at most 10⁻¹².

Results are in a unique `runs/unity/restart-check-JOBID.XXXXXXXX/results` folder,
or the corresponding custom Unity run root. They include `report.toml`, the
interior checkpoint, exact config, validator, model source, and environment
files. Success produces `RESTART_VALIDATION_PASSED` and exits zero; any failed
check exits nonzero. `FAILED.txt` records failures after the output is created.
The report identifies the device; a local CPU result is not evidence of a GPU
restart pass. This check makes no cross-node, cross-GPU, or long-time chaotic
trajectory claim. It exercises the model state restored by `restore_checkpoint!`,
not the output cadence of a resumed `run_experiment` directory.

Prepared runs can be submitted using the existing run mode with a Slurm success
dependency, without waiting interactively:

```bash
# Capture this immediately after submitting restart-check, before another job.
restart_job=$(cat "${SATURN_UNITY_ROOT:-$PWD/.unity}/last-job-id")
bash slurm/unity/submit.sh run configs/unity_followup_01/maintenance_forced.toml \
    --dependency="afterok:$restart_job"
bash slurm/unity/submit.sh run configs/unity_followup_01/maintenance_unforced.toml \
    --dependency="afterok:$restart_job"
```

The dependency holds both runs until the check succeeds. If the check fails,
Slurm will not release them; inspect/cancel or resubmit those held jobs after
addressing the failure. These examples do not submit anything automatically.

For a local CPU exercise of the validator, choose a new output directory:

```bash
bash scripts/julia-local.sh scripts/validate_gpu_restart.jl configs/smoke.toml \
    --output runs/local-restart-check --device cpu
```

The smoke configuration uses enlarged storms on a coarse grid. It checks the
validation procedure, while the Unity benchmark checks the intended spatial
resolution and forcing geometry. Keep production `src/` and `Manifest.toml`
unchanged when transferring only these validation additions: their hashes are
part of checkpoint compatibility.

# Unity deployment

Clone the project into allocated Unity work storage. Keep your SSH key and
allocation-specific paths outside the repository. The scripts use the default
Slurm association unless `SATURN_UNITY_ACCOUNT` is explicitly set.

From the project root on the Unity login node:

```sh
bash slurm/unity/submit.sh validate
```

Installation and validation run inside a scheduled GPU allocation. They install
the pinned Julia 1.10.12 Linux x86_64 runtime and packages under `.unity/`, run
CPU tests, compare short CPU/GPU trajectories and write a GPU smoke result.
The supplied defaults use one x86_64 A100 GPU and the `gpu-highvram,uri-gpu`
partitions used by the original deployment. Select partitions available to your
account when needed; these are not universal cluster access guarantees.

The production source and dependency environment have been validated on Unity.
Historical CPU/GPU relative L2 differences were approximately 10^-15 for the
short smoke comparison; this does not imply long-time chaotic trajectory
agreement or scientific convergence.

After successful validation, submit only the first six prepared experiments:

```sh
bash slurm/unity/experiments_01.sh submit configs/experiments_01 screen
bash slurm/unity/experiments_01.sh status configs/experiments_01 screen
```

The wrapper verifies source/dependency/configuration hashes and keeps a
stage-specific submission receipt. Once all six complete:

```sh
bash slurm/unity/experiments_01.sh collect configs/experiments_01 screen
```

Download the archive printed by `collect` to your analysis machine. Follow
[the gallery workflow](gallery.md) to render its completed result directories
and publish the compact figures. GitHub authentication is not required by the
simulation jobs.

See [the experiment protocol](experiments-01.md) for later stages, physical
interpretation, storage estimates and restart limitations. Consult the
[Unity documentation](https://docs.unityhpc.org/) for current site policies.

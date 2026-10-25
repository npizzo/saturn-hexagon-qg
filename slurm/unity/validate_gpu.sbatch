#!/usr/bin/env bash
# Submit via slurm/unity/submit.sh, which supplies all Slurm options.
set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Use bash slurm/unity/submit.sh validate}"
: "${SLURM_JOB_ID:?Submit through Slurm}"
source "$SATURN_PROJECT/slurm/unity/environment.sh"
[[ $# == 0 ]] || unity_die 'validation takes no configuration arguments'
[[ -z ${SATURN_RESUME:-} ]] || unity_die 'unset SATURN_RESUME for validation'
printf 'Installing/checking the Unity runtime on %s\n' "$(hostname)"
srun --ntasks=1 bash "$SATURN_PROJECT/slurm/bootstrap.sh" \
    "$SATURN_UNITY_ROOT/runtime" "$JULIA_DEPOT_PATH"
# Bash ignores the DeltaAI #SBATCH defaults in this shared workflow. The
# allocation is the Unity job already created by submit.sh, not a second job.
exec bash "$SATURN_PROJECT/slurm/validate_gpu.sbatch"

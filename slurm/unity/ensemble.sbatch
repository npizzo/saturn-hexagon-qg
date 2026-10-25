#!/usr/bin/env bash
set -euo pipefail
: "${SATURN_PROJECT:?Use bash slurm/unity/submit.sh ensemble --array=1-N}"
: "${SLURM_JOB_ID:?Submit through Slurm}"
source "$SATURN_PROJECT/slurm/unity/environment.sh"
if [[ -n ${SATURN_EXPERIMENT_ATTEMPT:-} || -n ${SATURN_EXPERIMENT_INTEGRITY_SHA256:-} ]]; then
    : "${SATURN_EXPERIMENT_ATTEMPT:?Missing experiment attempt path}"
    : "${SATURN_EXPERIMENT_INTEGRITY_SHA256:?Missing experiment integrity hash}"
    python3 "$SATURN_PROJECT/slurm/unity/experiments_01.py" verify-task
fi
exec bash "$SATURN_PROJECT/slurm/ensemble.sbatch" "$@"

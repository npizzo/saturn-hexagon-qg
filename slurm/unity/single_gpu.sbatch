#!/usr/bin/env bash
set -euo pipefail
: "${SATURN_PROJECT:?Use bash slurm/unity/submit.sh run CONFIG}"
: "${SLURM_JOB_ID:?Submit through Slurm}"
source "$SATURN_PROJECT/slurm/unity/environment.sh"
exec bash "$SATURN_PROJECT/slurm/single_gpu.sbatch" "$@"

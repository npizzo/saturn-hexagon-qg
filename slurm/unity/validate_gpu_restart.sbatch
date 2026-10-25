#!/usr/bin/env bash
# Submit via `bash slurm/unity/submit.sh restart-check CONFIG`.
set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Use bash slurm/unity/submit.sh restart-check CONFIG}"
source "$SATURN_PROJECT/slurm/unity/environment.sh"
source "$SATURN_PROJECT/slurm/common.sh"
[[ $# == 1 ]] || unity_die 'restart-check requires exactly one short convection configuration'
[[ -z ${SATURN_RESUME:-} ]] || unity_die 'unset SATURN_RESUME; restart-check creates its own checkpoint'
prepare_environment
job_dir=$(mktemp -d "$SATURN_RUN_ROOT/restart-check-${SLURM_JOB_ID}.XXXXXXXX")
stage_config "$1"
printf 'Restart validation results: %s\n' "$job_dir/results"
# Reuse the validated runtime/depot. No bootstrap, install, or CPU fallback.
srun --ntasks=1 "$JULIA_BIN" --startup-file=no --history-file=no \
    "--project=$SATURN_PROJECT" "$SATURN_PROJECT/scripts/validate_gpu_restart.jl" \
    "$job_dir/config.toml" --output "$job_dir/results" --device gpu

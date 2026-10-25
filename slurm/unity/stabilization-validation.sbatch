#!/usr/bin/env bash
# Submitted only by stabilization.sh. Scheduler options come from its wrapper.
set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Use stabilization.sh submit SUITE validation}"
: "${SLURM_JOB_ID:?Submit through Slurm}"
source "$SATURN_PROJECT/slurm/unity/environment.sh"
source "$SATURN_PROJECT/slurm/common.sh"
[[ $# == 0 ]] || die 'supply validation configs through the frozen stage manifest'
[[ -z ${SATURN_RESUME:-} ]] || die 'validation does not resume another run'
: "${SATURN_EXPERIMENT_ATTEMPT:?Missing frozen attempt path}"
: "${SATURN_EXPERIMENT_INTEGRITY_SHA256:?Missing frozen integrity hash}"
python3 "$SATURN_PROJECT/slurm/unity/experiments_01.py" verify-task
prepare_environment
: "${SATURN_CONFIG_LIST:?Missing frozen configuration manifest}"
: "${SATURN_CONFIG_LIST_SHA256:?Missing manifest hash}"
: "${SLURM_ARRAY_JOB_ID:?Submit as a Slurm array}"
: "${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"
require_absolute SATURN_CONFIG_LIST "$SATURN_CONFIG_LIST"
require_id SLURM_ARRAY_JOB_ID "$SLURM_ARRAY_JOB_ID"
require_id SLURM_ARRAY_TASK_ID "$SLURM_ARRAY_TASK_ID"
[[ $SATURN_CONFIG_LIST_SHA256 =~ ^[0-9a-f]{64}$ ]] || die 'invalid manifest SHA256'
actual=$(sha256sum "$SATURN_CONFIG_LIST")
[[ ${actual%% *} == "$SATURN_CONFIG_LIST_SHA256" ]] || die 'array manifest changed after preparation'
count=$(awk 'END { print NR }' "$SATURN_CONFIG_LIST")
(( SLURM_ARRAY_TASK_ID >= 1 && SLURM_ARRAY_TASK_ID <= count )) || die 'array task ID outside manifest bounds'
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SATURN_CONFIG_LIST")
IFS=$'\t' read -r expected config <<< "$line"
[[ $expected =~ ^[0-9a-f]{64}$ && $config == /* ]] || die 'malformed array manifest row'
job_dir=$(mktemp -d "$SATURN_RUN_ROOT/array-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.XXXXXXXX")
stage_config "$config"
actual=$(sha256sum "$job_dir/config.toml")
[[ ${actual%% *} == "$expected" ]] || die 'config changed after array preparation'
cp -- "$SATURN_CONFIG_LIST" "$job_dir/array-manifest.tsv"
printf 'Validation result directory: %s\n' "$job_dir/results"
srun --ntasks=1 "$JULIA_BIN" --startup-file=no --history-file=no \
    "--project=$SATURN_PROJECT" "$SATURN_PROJECT/scripts/validate_stabilization.jl" \
    "$job_dir/config.toml" --output "$job_dir/results"
[[ -f $job_dir/results/CPU_GPU_VALIDATION_PASSED && -f $job_dir/results/COMPLETED ]] || die 'validator did not produce both success markers'

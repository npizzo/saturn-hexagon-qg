#!/usr/bin/env bash
#SBATCH --job-name=saturn-ensemble
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=saturn-%A_%a.out
#SBATCH --error=saturn-%A_%a.err
# Array size is supplied explicitly to sbatch --array=1-N%CONCURRENCY.

set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Set SATURN_PROJECT before sbatch}"
[[ $SATURN_PROJECT == /* ]] || { printf 'SATURN_PROJECT must be absolute\n' >&2; exit 1; }
source "$SATURN_PROJECT/slurm/common.sh"
[[ $# == 0 ]] || die 'supply configs through SATURN_CONFIG_LIST, not positional arguments'
prepare_environment
: "${SATURN_CONFIG_LIST:?Prepare a manifest with slurm/prepare_array.sh}"
: "${SATURN_CONFIG_LIST_SHA256:?Use the manifest SHA256 printed by prepare_array.sh}"
: "${SLURM_ARRAY_JOB_ID:?Submit with sbatch --array=1-N}"
: "${SLURM_ARRAY_TASK_ID:?Submit with sbatch --array=1-N}"
require_absolute SATURN_CONFIG_LIST "$SATURN_CONFIG_LIST"
require_id SLURM_ARRAY_JOB_ID "$SLURM_ARRAY_JOB_ID"
require_id SLURM_ARRAY_TASK_ID "$SLURM_ARRAY_TASK_ID"
[[ -z ${SATURN_RESUME:-} ]] || die 'array resume is intentionally disabled; resume an individual config with single_gpu.sbatch'
[[ $SATURN_CONFIG_LIST_SHA256 =~ ^[0-9a-f]{64}$ ]] || die 'invalid manifest SHA256'
[[ -f $SATURN_CONFIG_LIST ]] || die 'array manifest is missing'
actual=$(sha256sum "$SATURN_CONFIG_LIST")
[[ ${actual%% *} == "$SATURN_CONFIG_LIST_SHA256" ]] || die 'array manifest has changed since preparation'
count=$(awk 'END { print NR }' "$SATURN_CONFIG_LIST")
(( SLURM_ARRAY_TASK_ID >= 1 && SLURM_ARRAY_TASK_ID <= count )) || die 'array task ID is outside manifest bounds'
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SATURN_CONFIG_LIST")
IFS=$'\t' read -r expected config <<< "$line"
[[ $expected =~ ^[0-9a-f]{64}$ && $config == /* ]] || die 'malformed array manifest row'
job_dir=$(mktemp -d "$SATURN_RUN_ROOT/array-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.XXXXXXXX")
stage_config "$config"
actual=$(sha256sum "$job_dir/config.toml")
[[ ${actual%% *} == "$expected" ]] || die 'config has changed since array preparation; refusing to launch'
cp -- "$SATURN_CONFIG_LIST" "$job_dir/array-manifest.tsv"
launch_gpu

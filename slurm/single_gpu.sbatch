#!/usr/bin/env bash
#SBATCH --job-name=saturn-qg
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=saturn-%j.out
#SBATCH --error=saturn-%j.err

set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Set SATURN_PROJECT before sbatch}"
[[ $SATURN_PROJECT == /* ]] || { printf 'SATURN_PROJECT must be absolute\n' >&2; exit 1; }
source "$SATURN_PROJECT/slurm/common.sh"
[[ $# == 1 ]] || die 'usage: sbatch -A ACCOUNT slurm/single_gpu.sbatch CONFIG.toml'
prepare_environment
job_dir=$(mktemp -d "$SATURN_RUN_ROOT/job-${SLURM_JOB_ID}.XXXXXXXX")
stage_config "$1"
launch_gpu

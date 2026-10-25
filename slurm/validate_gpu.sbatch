#!/usr/bin/env bash
# CPU checks, GPU agreement, and saved GPU output. DeltaAI scheduler defaults.
#SBATCH --job-name=saturn-validate
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=saturn-validate-%j.out
#SBATCH --error=saturn-validate-%j.err

set -euo pipefail
umask 027
: "${SATURN_PROJECT:?Set SATURN_PROJECT before sbatch}"
[[ $SATURN_PROJECT == /* ]] || { printf 'SATURN_PROJECT must be absolute\n' >&2; exit 1; }
source "$SATURN_PROJECT/slurm/common.sh"
[[ $# == 0 ]] || die 'usage: sbatch -A ACCOUNT slurm/validate_gpu.sbatch'
[[ -z ${SATURN_RESUME:-} ]] || die 'unset SATURN_RESUME for a fresh validation job'
prepare_environment
job_dir=$(mktemp -d "$SATURN_RUN_ROOT/validation-${SLURM_JOB_ID}.XXXXXXXX")
stage_config "$SATURN_PROJECT/configs/smoke.toml"
printf 'Validation directory: %s\n' "$job_dir"

printf '\nChecking the test suite on allocated CPUs\n'
srun --ntasks=1 "$JULIA_BIN" --startup-file=no --history-file=no \
  "--project=$SATURN_PROJECT" "$SATURN_PROJECT/test/runtests.jl"
printf '\nRecording CUDA environment\n'
srun --ntasks=1 "$JULIA_BIN" --startup-file=no --history-file=no \
  "--project=$SATURN_PROJECT" -e 'using CUDA; CUDA.versioninfo()'
printf '\nChecking short-time CPU/GPU agreement\n'
srun --ntasks=1 "$JULIA_BIN" --startup-file=no --history-file=no \
  "--project=$SATURN_PROJECT" "$SATURN_PROJECT/scripts/validate_gpu.jl" "$job_dir/config.toml"
printf '\nRunning the GPU smoke case with field and checkpoint output\n'
launch_gpu
printf '\nVALIDATION PASSED. Results: %s/results\n' "$job_dir"

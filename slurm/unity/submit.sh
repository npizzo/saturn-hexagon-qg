#!/usr/bin/env bash
# Submit from a Unity login node. No package installation or computation here.
set -euo pipefail
umask 027
source "$(dirname -- "${BASH_SOURCE[0]}")/environment.sh"
[[ -z ${SLURM_JOB_ID:-} ]] || unity_die 'submit from a login shell outside an existing allocation'
[[ $# -ge 1 ]] || unity_die 'usage: submit.sh validate [SBATCH_OPTIONS] | run CONFIG [SBATCH_OPTIONS] | restart-check CONFIG [SBATCH_OPTIONS] | ensemble [SBATCH_OPTIONS]'
command -v sbatch >/dev/null || unity_die 'sbatch is unavailable; run this on Unity'
mode=$1
shift
case $mode in
    validate)
        [[ -z ${SATURN_RESUME:-} ]] || unity_die 'unset SATURN_RESUME for validation'
        script=validate_gpu.sbatch
        ;;
    run|restart-check)
        [[ $# -ge 1 ]] || unity_die "$mode requires a configuration path"
        config=$1
        shift
        [[ $config == /* ]] || config=$SATURN_PROJECT/$config
        [[ -f $config ]] || unity_die "configuration is missing: $config"
        if [[ $mode == run ]]; then
            script=single_gpu.sbatch
        else
            [[ -z ${SATURN_RESUME:-} ]] || unity_die 'unset SATURN_RESUME; restart-check creates its own checkpoint'
            script=validate_gpu_restart.sbatch
        fi
        ;;
    ensemble)
        [[ -n ${SATURN_CONFIG_LIST:-} && -n ${SATURN_CONFIG_LIST_SHA256:-} ]] || unity_die 'first prepare a configuration manifest with slurm/prepare_array.sh and export its path and hash'
        [[ -z ${SATURN_RESUME:-} ]] || unity_die 'array resume is disabled; use run with one checkpoint'
        script=ensemble.sbatch
        ;;
    *) unity_die "unknown mode: $mode" ;;
esac
if [[ $mode != validate ]]; then
    [[ -x $JULIA_BIN && -d $JULIA_DEPOT_PATH ]] || unity_die 'run validate first to install and check this Unity environment'
fi
mkdir -p -- "$SATURN_UNITY_ROOT" "$SATURN_RUN_ROOT/logs"
log_pattern="$SATURN_RUN_ROOT/logs/saturn-${mode}-%j"
[[ $mode != ensemble ]] || log_pattern="$SATURN_RUN_ROOT/logs/saturn-${mode}-%A_%a"
# Explicit scheduler defaults and user options define this submission. Remove
# inherited SBATCH_* defaults only in the child (including account/QoS/array).
sbatch_args=(--parsable --export=ALL \
    --job-name="saturn-$mode" --partition=gpu-highvram,uri-gpu \
    --constraint='x86_64&a100' --nodes=1 --ntasks=1 --gpus=1 \
    --cpus-per-task=8 --mem=32G --time=01:00:00 --no-requeue \
    --chdir="$SATURN_PROJECT" --output="$log_pattern.out" --error="$log_pattern.err")
[[ -z ${SATURN_UNITY_ACCOUNT:-} ]] || sbatch_args+=(--account="$SATURN_UNITY_ACCOUNT")
sbatch_args+=("$@" "$SATURN_PROJECT/slurm/unity/$script")
[[ $mode != run && $mode != restart-check ]] || sbatch_args+=("$config")
job=$(
    for variable in "${!SBATCH_@}"; do unset "$variable"; done
    sbatch "${sbatch_args[@]}"
)
job=${job%%;*}
[[ $job =~ ^[0-9]+$ ]] || unity_die "unexpected sbatch response: $job"
printf '%s\n' "$job" > "$SATURN_UNITY_ROOT/last-job-id"
printf 'Submitted Unity %s job %s\n' "$mode" "$job"
printf 'Logs: %s\n' "$SATURN_RUN_ROOT/logs"
printf 'Monitor: squeue -j %s\n' "$job"

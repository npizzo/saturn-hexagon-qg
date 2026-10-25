#!/usr/bin/env bash
# Sourced by Slurm entry points; no module purge, device remapping, or global install.
die() { printf 'saturn job: %s\n' "$*" >&2; exit 1; }
require_absolute() { [[ $2 == /* ]] || die "$1 must be an absolute path"; }
require_id() { [[ $2 =~ ^[0-9]+$ ]] || die "$1 must contain only decimal digits"; }

prepare_environment() {
    : "${SATURN_PROJECT:?Set SATURN_PROJECT to the absolute project directory}"
    : "${SATURN_RUN_ROOT:?Set SATURN_RUN_ROOT to an absolute directory on shared work storage}"
    : "${JULIA_BIN:?Set JULIA_BIN to an absolute Julia executable}"
    : "${JULIA_DEPOT_PATH:?Set JULIA_DEPOT_PATH to the explicit project depot}"
    : "${SLURM_JOB_ID:?Submit this script through sbatch}"
    : "${SLURM_JOB_ACCOUNT:?Slurm must supply a job account; use the site default or sbatch -A ACCOUNT}"
    require_absolute SATURN_PROJECT "$SATURN_PROJECT"
    require_absolute SATURN_RUN_ROOT "$SATURN_RUN_ROOT"
    require_absolute JULIA_BIN "$JULIA_BIN"
    require_id SLURM_JOB_ID "$SLURM_JOB_ID"
    [[ -x $JULIA_BIN ]] || die "Julia is not executable: $JULIA_BIN"
    [[ -f $SATURN_PROJECT/Project.toml && -f $SATURN_PROJECT/Manifest.toml ]] || die 'project or pinned Manifest.toml is missing'
    [[ -f $SATURN_PROJECT/scripts/run.jl ]] || die 'runner is missing'
    # A single explicit depot prevents accidental writes into the default user depot.
    require_absolute JULIA_DEPOT_PATH "$JULIA_DEPOT_PATH"
    [[ $JULIA_DEPOT_PATH != *:* ]] || die 'use a single explicit depot directory, without colons'
    [[ -d $JULIA_DEPOT_PATH ]] || die 'bootstrap the Julia depot before submitting'
    export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-8}}
    export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
    export JULIA_PKG_PRECOMPILE_AUTO=${JULIA_PKG_PRECOMPILE_AUTO:-0}
    require_id JULIA_NUM_THREADS "$JULIA_NUM_THREADS"
    (( JULIA_NUM_THREADS >= 1 && JULIA_NUM_THREADS <= ${SLURM_CPUS_PER_TASK:-8} )) || die 'JULIA_NUM_THREADS must fit within the allocated CPUs'
    mkdir -p -- "$SATURN_RUN_ROOT"
    cd -- "$SATURN_PROJECT"
}

stage_config() {
    local source_config=$1
    [[ $source_config == /* ]] || source_config=$SATURN_PROJECT/$source_config
    [[ -f $source_config ]] || die "config is missing: $source_config"
    cp -- "$source_config" "$job_dir/config.toml.tmp"
    mv -- "$job_dir/config.toml.tmp" "$job_dir/config.toml"
    # Unique parent directory is reserved atomically; the runner owns its result directory.
    {
        printf 'job_id=%s\naccount=%s\nhost=%s\n' "$SLURM_JOB_ID" "$SLURM_JOB_ACCOUNT" "$(hostname)"
        printf 'project=%s\nsource_config=%s\njulia=%s\n' "$SATURN_PROJECT" "$source_config" "$JULIA_BIN"
        printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        sha256sum "$job_dir/config.toml" "$SATURN_PROJECT/Manifest.toml"
        "$JULIA_BIN" --startup-file=no --version
    } > "$job_dir/job-context.txt.tmp"
    mv -- "$job_dir/job-context.txt.tmp" "$job_dir/job-context.txt"
}

launch_gpu() {
    local -a julia_args=(--startup-file=no --history-file=no
      "--project=$SATURN_PROJECT" "$SATURN_PROJECT/scripts/run.jl"
      "$job_dir/config.toml" --output "$job_dir/results" --device gpu)
    if [[ -n ${SATURN_RESUME:-} ]]; then
        require_absolute SATURN_RESUME "$SATURN_RESUME"
        [[ -f $SATURN_RESUME ]] || die 'resume checkpoint is missing'
        julia_args+=(--resume "$SATURN_RESUME")
    fi
    printf 'Result directory: %s\n' "$job_dir/results"
    # Slurm controls CUDA_VISIBLE_DEVICES. The runner rejects unavailable CUDA.
    srun --ntasks=1 "$JULIA_BIN" "${julia_args[@]}"
}

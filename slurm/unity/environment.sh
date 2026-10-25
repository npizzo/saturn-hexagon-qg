#!/usr/bin/env bash
# Source from the submitter or a batch wrapper. Use Unity-specific paths even
# when the terminal still contains exports from a previous DeltaAI deployment.
unity_die() { printf 'saturn Unity: %s\n' "$*" >&2; exit 1; }
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || unity_die 'use Unity Linux x86_64 nodes; this is not the Mac launcher'
export SATURN_PROJECT
SATURN_PROJECT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
export SATURN_UNITY_ROOT=${SATURN_UNITY_ROOT:-$SATURN_PROJECT/.unity}
export SATURN_RUN_ROOT=${SATURN_UNITY_RUN_ROOT:-$SATURN_PROJECT/runs/unity}
[[ $SATURN_UNITY_ROOT == /* && $SATURN_RUN_ROOT == /* ]] || unity_die 'Unity runtime and run roots must be absolute paths'
export JULIA_BIN="$SATURN_UNITY_ROOT/runtime/julia-1.10.12/bin/julia"
export JULIA_DEPOT_PATH="$SATURN_UNITY_ROOT/julia-depot"
# Shared caches must also work when a later job lands on a different x86 CPU.
export JULIA_CPU_TARGET=generic
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=0

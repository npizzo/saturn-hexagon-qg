#!/usr/bin/env bash
# Reuse the installed Unity runtime, with a separate checkout and run directory.
set -euo pipefail
umask 027
stabilization_project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
stabilization_previous="$(dirname -- "$stabilization_project")/saturn-hexagon-qg"
export SATURN_UNITY_ROOT=${SATURN_STABILIZATION_RUNTIME_ROOT:-$stabilization_previous/.unity}
export SATURN_UNITY_RUN_ROOT=${SATURN_STABILIZATION_RUN_ROOT:-$stabilization_project/runs/stabilization-01}
source "$stabilization_project/slurm/unity/environment.sh"
command -v python3 >/dev/null || unity_die 'Python 3 is required for stabilization manifests'
exec python3 "$SATURN_PROJECT/slurm/unity/stabilization.py" "$@"

#!/usr/bin/env bash
# Separate 2017-benchmark checkout and outputs; reuse the working Unity runtime.
set -euo pipefail
umask 027
rostami_project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
rostami_previous="$(dirname -- "$rostami_project")/saturn-hexagon-qg"
export SATURN_UNITY_ROOT=${SATURN_ROSTAMI_RUNTIME_ROOT:-$rostami_previous/.unity}
export SATURN_UNITY_RUN_ROOT=${SATURN_ROSTAMI_RUN_ROOT:-$rostami_project/runs/rostami-01}
source "$rostami_project/slurm/unity/environment.sh"
command -v python3 >/dev/null || unity_die 'Python 3 is required for Rostami campaign manifests'
exec python3 "$SATURN_PROJECT/slurm/unity/rostami.py" "$@"

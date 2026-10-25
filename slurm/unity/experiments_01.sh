#!/usr/bin/env bash
# Explicitly prepare or submit one reviewed experiment stage; never chain stages.
set -euo pipefail
umask 027
source "$(dirname -- "${BASH_SOURCE[0]}")/environment.sh"
command -v python3 >/dev/null || unity_die 'Python 3 is required for experiment manifests'
exec python3 "$SATURN_PROJECT/slurm/unity/experiments_01.py" "$@"

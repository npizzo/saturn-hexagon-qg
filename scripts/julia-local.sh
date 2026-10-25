#!/usr/bin/env bash
# Portable wrapper; use JULIA_BIN to select an existing installation.
set -euo pipefail
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
julia_executable=${JULIA_BIN:-julia}
cd -- "$project_root"
exec "$julia_executable" --startup-file=no --project="$project_root" "$@"

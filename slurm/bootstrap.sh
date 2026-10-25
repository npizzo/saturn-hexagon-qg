#!/usr/bin/env bash
# Install the pinned Linux runtime and project dependencies in explicit paths.
set -euo pipefail
umask 027

die() { printf 'bootstrap: %s\n' "$*" >&2; exit 1; }
[[ $# == 2 ]] || die 'usage: bash slurm/bootstrap.sh /absolute/runtime-root /absolute/julia-depot'
runtime_root=$1
depot_root=$2
[[ $runtime_root == /* && $depot_root == /* ]] || die 'both directories must be absolute paths'
[[ $(uname -s) == Linux ]] || die 'this bootstrap is for Linux compute nodes, not a local Mac'
runtime_arch=$(uname -m)
case $runtime_arch in
    aarch64)
        runtime_url='https://julialang-s3.julialang.org/bin/linux/aarch64/1.10/julia-1.10.12-linux-aarch64.tar.gz'
        runtime_sha256='cc257b21967049642bfe17de1aff7dae17477eba4f423f98cb400c002080e5a2'
        ;;
    x86_64)
        runtime_url='https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.12-linux-x86_64.tar.gz'
        runtime_sha256='b0c50dfe349d76ffc8580d24e0c0be623988a6d91ed57ad8b73c0137fd27bc0d'
        ;;
    *) die "unsupported Linux architecture: $runtime_arch" ;;
esac
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
[[ -f $project_root/Project.toml && -f $project_root/Manifest.toml ]] || die 'Project.toml and the committed Manifest.toml are required'
for tool in curl sha256sum tar mktemp; do command -v "$tool" >/dev/null || die "missing command: $tool"; done
mkdir -p -- "$runtime_root" "$depot_root"
runtime_root=$(cd -- "$runtime_root" && pwd -P)
depot_root=$(cd -- "$depot_root" && pwd -P)
runtime_dir=$runtime_root/julia-1.10.12
lock_dir=$runtime_root/.saturn-julia-bootstrap.lock
mkdir -- "$lock_dir" || die "another bootstrap may be running; inspect $lock_dir before removing a stale lock"
stage_dir=''
cleanup() {
    [[ -z $stage_dir ]] || rm -rf -- "$stage_dir"
    rmdir -- "$lock_dir"
}
trap cleanup EXIT

if [[ ! -e $runtime_dir ]]; then
    stage_dir=$(mktemp -d "$runtime_root/.saturn-julia.XXXXXXXX")
    archive=$stage_dir/julia.tar.gz
    curl --fail --location --retry 3 --output "$archive" "$runtime_url"
    printf '%s  %s\n' "$runtime_sha256" "$archive" | sha256sum --check -
    tar -xzf "$archive" -C "$stage_dir"
    [[ -x $stage_dir/julia-1.10.12/bin/julia ]] || die 'archive does not contain the expected runtime'
    mv -- "$stage_dir/julia-1.10.12" "$runtime_dir"
fi
[[ -x $runtime_dir/bin/julia ]] || die "existing runtime is incomplete: $runtime_dir"

# These assignments apply to the child only. No shell profile or global depot is modified.
SATURN_BOOTSTRAP_ARCH="$runtime_arch" JULIA_DEPOT_PATH="$depot_root" JULIA_PKG_PRECOMPILE_AUTO=0 \
  "$runtime_dir/bin/julia" --startup-file=no --history-file=no --project="$project_root" -e \
  'VERSION == v"1.10.12" || error("Expected Julia 1.10.12"); string(Sys.ARCH) == ENV["SATURN_BOOTSTRAP_ARCH"] || error("Julia architecture does not match this node"); using Pkg; Pkg.instantiate(); Pkg.precompile()'
printf '\nBootstrap finished. Set these variables in the submitting shell:\n'
printf 'export SATURN_PROJECT=%q\n' "$project_root"
printf 'export JULIA_BIN=%q\n' "$runtime_dir/bin/julia"
printf 'export JULIA_DEPOT_PATH=%q\n' "$depot_root"
printf 'A GPU allocation smoke test is still required; package installation is not GPU validation.\n'

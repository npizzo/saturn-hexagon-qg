#!/usr/bin/env bash
# Freeze an ordered list of absolute config paths and their content hashes.
set -euo pipefail
umask 027
die() { printf 'prepare array: %s\n' "$*" >&2; exit 1; }
[[ $# -ge 2 ]] || die 'usage: bash slurm/prepare_array.sh /absolute/manifest.tsv CONFIG.toml ...'
manifest=$1
shift
[[ $manifest == /* ]] || die 'manifest path must be absolute'
[[ ! -e $manifest ]] || die 'manifest already exists; choose a new name'
command -v sha256sum >/dev/null || die 'sha256sum is required (run this on DeltaAI)'
mkdir -p -- "$(dirname -- "$manifest")"
temporary=$(mktemp "${manifest}.XXXXXXXX.tmp")
trap 'rm -f -- "$temporary"' EXIT
for config in "$@"; do
    [[ -f $config ]] || die "config is missing: $config"
    config=$(cd -- "$(dirname -- "$config")" && printf '%s/%s\n' "$(pwd -P)" "$(basename -- "$config")")
    [[ $config != *$'\t'* && $config != *$'\n'* && $config != *$'\r'* ]] || die 'config paths must not contain tabs or line breaks'
    digest=$(sha256sum "$config")
    digest=${digest%% *}
    printf '%s\t%s\n' "$digest" "$config" >> "$temporary"
done
# Hard-link publication is atomic and refuses to overwrite an existing manifest.
chmod 440 "$temporary"
ln -- "$temporary" "$manifest" || die 'could not publish manifest without overwriting'
digest=$(sha256sum "$manifest")
digest=${digest%% *}
printf 'export SATURN_CONFIG_LIST=%q\n' "$manifest"
printf 'export SATURN_CONFIG_LIST_SHA256=%q\n' "$digest"
printf 'Array range: 1-%s (add a concurrency cap such as %%4)\n' "$#"

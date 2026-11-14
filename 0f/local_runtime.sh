#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Shared local runtime lookup

Main options:
    workflow_local_runtime ENV EXECUTABLE FALLBACK
    Resolve a runtime from the configured Conda roots.

This is a sourced helper library, not an analysis entry point.

Examples:
    source ./0os.sh
    source ./local_runtime.sh
    workflow_local_runtime le8 python3 python3
PUBLIC_ENTRY_HELP
    exit 0
fi
# Locate an existing local runtime without activation, installation or imports.
# Callers retain explicit CLI/environment overrides before using this fallback.
workflow_local_runtime() {
    local environment=$1 executable=$2 fallback=$3 base prefix
    for base in "${OS_GU_CONDA:-}" "${DIR_SOFTWARE:-}/anaconda3" \
                "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/miniconda3"; do
        [[ -n $base ]] || continue
        prefix=$base
        [[ -z $environment ]] || prefix=$base/envs/$environment
        if [[ -x $prefix/bin/$executable ]]; then
            printf '%s\n' "$prefix/bin/$executable"
            return 0
        fi
    done
    printf '%s\n' "$fallback"
}

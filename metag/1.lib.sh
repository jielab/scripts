#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Metagenomics reference libraries

Modules:
    check   Verify existing libraries
    human   Build the human reference
    kraken  Build the Kraken library
    sylph   Build the Sylph library

Main options:
    --config FILE   Runtime configuration [config.json]
    --full          Full verification (check only)
    --plan-only     Preview verification commands (check only)

Full checks submit a metag.check job on HPC. Library operations do not run
sample analysis. Use ./1.lib.sh MODULE -h for build-specific options.

Examples:
    ./1.lib.sh check
    ./1.lib.sh check --full --plan-only
    ./1.lib.sh check --full
PUBLIC_ENTRY_HELP
    exit 0
fi
# Library operations; full checks never launch sample analysis.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module=${1:---help}
if (($#)); then shift; fi

# Only full verification needs a scheduled compute allocation.
if [[ $module == check ]]; then
    full=false; preview=false; output=""; args=(); config=false
    while (($#)); do
        case "$1" in
            -h|--help) exec bash "$HERE/f/usage.sh" check ;;
            --full) full=true; args+=("$1"); shift ;;
            --plan-only) preview=true; shift ;;
            -o|--out-dir)
                [[ $# -ge 2 && -n $2 && $2 != -* ]] || { echo "$1 requires a directory" >&2; exit 2; }
                output=$2; shift 2 ;;
            --config)
                [[ $# -ge 2 && -f $2 ]] || { echo 'Missing --config file' >&2; exit 2; }
                config=true; args+=("$1" "$(realpath -- "$2")"); shift 2 ;;
            *) args+=("$1"); shift ;;
        esac
    done
    if [[ $full == true ]]; then
        source "$HERE/../0f/0os.sh"
        [[ $config == true ]] || args+=(--config "$HERE/config.json")
        [[ $submit != local ]] || args+=(--allow-local)
        output=${output:-${ANALYSIS_ROOT}/metag/lib}
        if [[ $submit == hpc && -z ${LSB_JOBID:-} ]]; then
            printf -v command '%q ' bash "$HERE/f/analysis_entry.sh" check "${args[@]}"
            request=(bsub -J metag.check -q 33 -n 40 -R 'span[hosts=1] rusage[mem=160GB/host]' -W 120:00 -o "$output/check.%J.out" -e "$output/check.%J.err" "$command")
            if [[ $preview == true ]]; then printf '%q ' "${request[@]}"; printf '\n'; exit 0; fi
            mkdir -p -- "$output"
            exec "${request[@]}"
        fi
    fi
    if [[ $preview == true ]]; then printf '%q ' bash "$HERE/f/analysis_entry.sh" check "${args[@]}"; printf '\n'; exit 0; fi
    exec bash "$HERE/f/analysis_entry.sh" check "${args[@]}"
fi

case "$module" in
    -h|--help|--help-full) bash "$HERE/f/usage.sh" "build_lib"; exit 0 ;;
    check) exec bash "$HERE/f/analysis_entry.sh" check "$@" ;;
    human) exec bash "$HERE/f/build_human.sh" "$@" ;;
    sylph) exec bash "$HERE/f/build_sylph.sh" "$@" ;;
    kraken)
        for arg in "$@"; do
            if [[ $arg == -h || $arg == --help ]]; then bash "$HERE/f/usage.sh" "build_kraken"; exit 0; fi
        done
        if [[ ${1:-} == --help-full ]]; then shift; set -- --help "$@"; fi
        exec bash "$HERE/f/entry.sh" library "$@" ;;
    *) echo "Unknown database module: $module. Choose check, kraken, human or sylph." >&2; exit 2 ;;
esac

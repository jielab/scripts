#!/usr/bin/env bash
# Shared process-tree limits and logging for the three root entry points.

fail() { printf '%s: %s\n' "${WORKFLOW_NAME:-allele3}" "$*" >&2; exit 2; }
# Help remains available even when an explicitly configured path is missing.
_allele_help=0
for _allele_arg in "$@"; do
    case "$_allele_arg" in -h|--help|help|--help-full|--version) _allele_help=1 ;; esac
done
if ((_allele_help == 0)); then
    if [[ ${ALLELE_CONFIG+x} && ! -f $ALLELE_CONFIG ]]; then fail "Explicit ALLELE_CONFIG does not exist: $ALLELE_CONFIG"; fi
    if [[ -n ${ROOT:-} && -f ${ALLELE_CONFIG:-"$ROOT/allele3.env"} ]]; then
        set -a
        source "${ALLELE_CONFIG:-"$ROOT/allele3.env"}"
        set +a
    fi
fi
unset _allele_arg _allele_help

need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a value"; }
positive_integer() { [[ $2 =~ ^[1-9][0-9]*$ && ${#2} -le 6 ]] || fail "$1 must be a positive integer"; }

verify_limits() {
    local expected=${1:?memory bytes required}
    shift
    local cgroup
    cgroup=$(awk -F: '$1 == "0" { print $3 }' /proc/self/cgroup)
    local base=/sys/fs/cgroup$cgroup
    if [[ -z $cgroup || ! -r $base/memory.max || $(<"$base/memory.max") != "$expected" ||
          $(<"$base/memory.swap.max") != 0 || $(<"$base/memory.oom.group") != 1 ]]; then
        printf '%s: cgroup RAM/swap/OOM-group verification failed; workload was not started\n' "${WORKFLOW_NAME:-allele3}" >&2
        exit 2
    fi
    printf '[%s] verified cgroup=%s memory.max=%s swap.max=0 oom.group=1\n' "${WORKFLOW_NAME:-allele3}" "$cgroup" "$expected" >&2
    exec "$@"
}

launch_workflow() {
    local command=$1 python=$2 entry=$3 out=$4 memory=$5 screen_interval=$6
    local verbose_screen=$7 dry_run=$8 logged=$9 native_help=${10}
    shift 10
    local -a args=("$@")
    local _path_index _path_option
    for ((_path_index=0; _path_index<${#args[@]}-1; _path_index++)); do
        _path_option=${args[_path_index]}
        case "$_path_option" in
            --metadata|--coordinates|--eigenvec|--basis-dir|--reuse-basis-dir|--target-glob|--target-info|--reuse-from|--fasta|--vcf-dir|--dir-vcf|--pfile-dir|--library|--ancient-dir|--dir-ancient|--source|--store|--scan-dir|--dir-modern)
                printf '[resolved path] %s = %s\n' "$_path_option" "${args[_path_index+1]}" ;;
        esac
    done
    if ((dry_run)); then
        printf 'MemoryMax=%s GiB; MemorySwapMax=0; OOMPolicy=kill; output=%s\n' "$memory" "$out"
        if ((logged)); then
            printf 'Screen progress interval=%ss; verbose=%s; full logs=%s/logs/\n' "$screen_interval" "$verbose_screen" "$out"
        fi
        printf '%q ' "$python" "$entry" "${args[@]}"; printf '\n'
        return 0
    fi
    command -v "$python" >/dev/null 2>&1 || fail "Python not found: $python. Install environment.yml or supply --python. Help: ./allele3.sh -h"
    # Use matching native tools (bcftools/PLINK) from the selected environment.
    python=$(command -v "$python")
    PATH="$(dirname -- "$python"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export MPLCONFIGDIR=${MPLCONFIGDIR:-"$HOME/.cache/allele3-evolution/matplotlib"}
    if ((native_help)); then exec "$python" "$entry" "${args[@]}"; fi
    # Non-login WSL terminals may omit these despite a running user manager.
    export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
    export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
    command -v systemd-run >/dev/null 2>&1 || fail 'systemd-run is required for the process-tree memory cap; analysis was not started'
    [[ -r /sys/fs/cgroup/cgroup.controllers ]] || fail 'cgroup v2 is required; analysis was not started'
    timeout -k 1s 5s systemctl --user show-environment >/dev/null 2>&1 || fail 'User systemd is unavailable; analysis was not started without a memory cap'
    local bytes=$((memory * 1024 * 1024 * 1024))
    local unit
    unit="allele3-$(date +%s)-$$-$RANDOM"
    printf '[allele3] command=%s output=%s RAM=%s GiB swap=0\n' "$command" "$out" "$memory" >&2
    local -a launch=("$python" "$entry" "${args[@]}")
    if ((logged)); then
        local -a console=("$python" "${WORKFLOW_CONSOLE:-$ROOT/f/0.log.py}" --log-dir "$out/logs"
            --label "$command" --interval "$screen_interval")
        ((verbose_screen == 0)) || console+=(--verbose)
        launch=("${console[@]}" -- "${launch[@]}")
    fi
    # The kernel charges the parent AND all descendants to the same cgroup.
    # The inner guard verifies the effective limits before launching Python.
    exec systemd-run --user --scope --quiet --unit="$unit" \
        -p "MemoryMax=$bytes" -p MemorySwapMax=0 -p OOMPolicy=kill -p TasksMax=256 \
        bash -c "$(declare -f verify_limits)"$'\nset -euo pipefail\nverify_limits "$@"' \
            allele3-guard "$bytes" "${launch[@]}"
}

workflow_python() {
    local requested=$1 candidate
    shift
    if [[ -n $requested ]]; then printf '%s\n' "$requested"; return; fi
    if [[ -n ${CONDA_PREFIX:-} && ${CONDA_DEFAULT_ENV:-} != base && -x $CONDA_PREFIX/bin/python ]]; then
        printf '%s\n' "$CONDA_PREFIX/bin/python"; return
    fi
    for candidate in "$@"; do
        if [[ -x $candidate ]]; then printf '%s\n' "$candidate"; return; fi
    done
    if [[ -n ${CONDA_PREFIX:-} && -x $CONDA_PREFIX/bin/python ]]; then
        printf '%s\n' "$CONDA_PREFIX/bin/python"; return
    fi
    printf 'python3\n'
}


workflow_launch() {
    local label=$1 python=$2 entry=$3 out=$4 memory=$5 interval=$6
    local verbose=$7 dry_run=$8 logged=$9
    shift 9
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    ((memory <= 32)) || fail 'Memory cap must be at most 32 GiB on this shared host'
    export WORKFLOW_NAME=$label
    launch_workflow "$label" "$python" "$entry" "$out" "$memory" "$interval" \
        "$verbose" "$dry_run" "$logged" "$((1-logged))" "$@"
}

#!/usr/bin/env bash
# Shared launcher for independent workflows; no analysis-specific imports.

fail() { printf '%s: %s\n' "$WORKFLOW_NAME" "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a value"; }
positive_integer() { [[ $2 =~ ^[1-9][0-9]*$ && ${#2} -le 6 ]] || fail "$1 must be a positive integer"; }

workflow_python() {
    local requested=$1 candidate
    shift
    if [[ -n $requested ]]; then printf '%s\n' "$requested"; return; fi
    if [[ -n ${CONDA_PREFIX:-} && -x $CONDA_PREFIX/bin/python ]]; then
        printf '%s\n' "$CONDA_PREFIX/bin/python"; return
    fi
    for candidate in "$@"; do
        if [[ -x $candidate ]]; then printf '%s\n' "$candidate"; return; fi
    done
    printf 'python3\n'
}

workflow_launch() {
    local label=$1 python=$2 entry=$3 out=$4 memory=$5
    local interval=$6 verbose=$7 dry_run=$8 logged=$9
    shift 9
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    ((memory <= 32)) || fail 'Memory cap must be at most 32 GiB on this shared host'
    if ((dry_run)); then
        printf 'MemoryMax=%s GiB; MemorySwapMax=0; OOMPolicy=kill; output=%s\n' "$memory" "$out"
        if ((logged)); then
            printf 'Screen progress interval=%ss; verbose=%s; full logs=%s/logs/\n' "$interval" "$verbose" "$out"
        fi
        printf '%q ' "$python" "$entry" "$@"; printf '\n'
        return 0
    fi
    command -v "$python" >/dev/null 2>&1 || fail "Python not found: $python. Supply --python or activate the workflow environment."
    python=$(command -v "$python")
    PATH="$(dirname -- "$python"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export WORKFLOW_NAME=$label
    local -a launch=("$python" "$entry" "$@")
    # Native help also works on machines without systemd; it launches no analysis.
    if (( !logged )); then exec "${launch[@]}"; fi
    export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
    export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
    command -v systemd-run >/dev/null 2>&1 || fail 'systemd-run is required for the process-tree memory cap; analysis was not started'
    [[ -r /sys/fs/cgroup/cgroup.controllers ]] || fail 'cgroup v2 is required; analysis was not started'
    timeout -k 1s 5s systemctl --user show-environment >/dev/null 2>&1 || fail 'User systemd is unavailable; analysis was not started without a memory cap'
    local bytes=$((memory * 1024 * 1024 * 1024)) unit
    unit="$label-$(date +%s)-$$-$RANDOM"
    printf '[%s] output=%s RAM=%s GiB swap=0\n' "$label" "$out" "$memory" >&2
    local -a console=("$python" "$ROOT/runtime/console_log.py" --log-dir "$out/logs"
        --label "$label" --interval "$interval")
    ((verbose == 0)) || console+=(--verbose)
    exec systemd-run --user --scope --quiet --unit="$unit" \
        -p "MemoryMax=$bytes" -p MemorySwapMax=0 -p OOMPolicy=kill -p TasksMax=256 \
        bash "$ROOT/runtime/limited_exec.sh" "$bytes" "${console[@]}" -- "${launch[@]}"
}

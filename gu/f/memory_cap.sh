#!/usr/bin/env bash
# A scope accounts for the whole process tree, including detached workers.
set -euo pipefail

memory_bytes(){
  local value=${1^^} number
  [[ $value =~ ^([1-9][0-9]*)(G|GB|GIB)$ ]] || {
    echo 'ERROR: --memory-cap requires a positive whole number of GiB (for example 32G).' >&2
    return 2
  }
  number=${BASH_REMATCH[1]}
  [[ ${#number} -le 6 ]] || return 2
  printf '%s\n' "$((10#$number * 1073741824))"
}

inside_cap(){
  local bytes=$1 current expected=${GU_MEMORY_CGROUP:-}
  current=$(awk -F: '$1=="0" {print $3}' /proc/self/cgroup)
  [[ -n $expected && $expected != / && ( $current == "$expected" || $current == "$expected/"* ) ]] || return 1
  [[ $(cat "/sys/fs/cgroup$expected/memory.max") == "$bytes" &&
     $(cat "/sys/fs/cgroup$expected/memory.swap.max") == 0 &&
     $(cat "/sys/fs/cgroup$expected/memory.oom.group") == 1 ]]
}

mode=${1:-}; shift
if [[ $mode == --check ]]; then
  bytes=$(memory_bytes "${1:?}")
  inside_cap "$bytes"
  exit
fi
if [[ $mode == --inside ]]; then
  size=${1:?}; shift
  bytes=$(memory_bytes "$size")
  GU_MEMORY_CGROUP=$(awk -F: '$1=="0" {print $3}' /proc/self/cgroup)
  [[ $GU_MEMORY_CGROUP == */gu-memory-*.scope ]] || {
    echo 'ERROR: memory cap scope was not created.' >&2; exit 1;
  }
  # systemd 255 does not expose memory.oom.group as a scope property.
  # Kill the entire run on OOM instead of leaving partial workers running.
  echo 1 > "/sys/fs/cgroup$GU_MEMORY_CGROUP/memory.oom.group"
  export GU_MEMORY_CGROUP
  inside_cap "$bytes" || { echo 'ERROR: memory cap verification failed.' >&2; exit 1; }
  echo "[GU MEMORY] cap=$size (GiB units) swap=0 scope=$GU_MEMORY_CGROUP; shared by all workers; exceeding cap terminates this run" >&2
  exec "$@"
fi
size=$mode
bytes=$(memory_bytes "$size")
if inside_cap "$bytes"; then exec "$@"; fi
[[ -f /sys/fs/cgroup/cgroup.controllers ]] && command -v systemd-run >/dev/null &&
  systemctl --user show-environment >/dev/null 2>&1 || {
    echo 'ERROR: cannot enforce --memory-cap: cgroup v2 and a running systemd user manager are required. Analysis was not started.' >&2
    exit 1
  }
exec systemd-run --user --scope --quiet --unit="gu-memory-$$-$RANDOM" \
  -p "MemoryMax=$bytes" -p MemorySwapMax=0 -p KillMode=control-group -p TimeoutStopSec=10s \
  bash "${BASH_SOURCE[0]}" --inside "$size" "$@"

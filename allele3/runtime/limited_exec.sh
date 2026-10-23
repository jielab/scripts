#!/usr/bin/env bash
# Internal cgroup guard; never run a workload unless all kernel limits match.
set -euo pipefail
expected=${1:?memory bytes required}
shift
cgroup=$(awk -F: '$1 == "0" { print $3 }' /proc/self/cgroup)
base=/sys/fs/cgroup$cgroup
if [[ -z $cgroup || ! -r $base/memory.max || $(<"$base/memory.max") != "$expected" ||
      $(<"$base/memory.swap.max") != 0 || $(<"$base/memory.oom.group") != 1 ]]; then
    printf '%s: cgroup RAM/swap/OOM-group verification failed; workload was not started\n' "${WORKFLOW_NAME:-allele3}" >&2
    exit 2
fi
printf '[%s] verified cgroup=%s memory.max=%s swap.max=0 oom.group=1\n' "${WORKFLOW_NAME:-allele3}" "$cgroup" "$expected" >&2
exec "$@"

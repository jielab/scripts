#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
command=split
case "${1:-}" in
  split|classify) command="$1"; shift;;
  *)
    for arg in "$@"; do
      case "$arg" in --r1|--r1=*) command=classify; break;; esac
    done;;
esac
export METAG_PROG="./2.split.sh"
exec bash "$HERE/f/run_workflow.sh" "$command" "$@"


#!/usr/bin/env bash
# Wait for every PID explicitly; short successful jobs must not turn into wait -n=127.
refgen_children=()
refgen_stop_tree(){
  local pid=$1 child
  while read -r child; do [[ -z $child ]] || refgen_stop_tree "$child"; done < <(pgrep -P "$pid" || true)
  kill -TERM "$pid" 2>/dev/null || true
}
refgen_stop_jobs(){
  local pid
  for pid in "${refgen_children[@]}"; do refgen_stop_tree "$pid"; done
  for pid in "${refgen_children[@]}"; do wait "$pid" 2>/dev/null || true; done
  refgen_children=()
}
trap 'refgen_stop_jobs; exit 130' INT
trap 'refgen_stop_jobs; exit 143' TERM
run_parallel(){
  local fn=$1 c pid failed=0 limit=${GRID_ARG_JOBS:-${JOBS:-1}}; shift
  refgen_children=()
  for c in "$@"; do
    "$fn" "$c" & refgen_children+=("$!")
    if ((${#refgen_children[@]}>=limit)); then
      for pid in "${refgen_children[@]}"; do wait "$pid" || failed=1; done
      refgen_children=()
      ((failed==0)) || { echo "ERROR: $fn failed; no further jobs started" >&2; return 1; }
    fi
  done
  for pid in "${refgen_children[@]}"; do wait "$pid" || failed=1; done
  refgen_children=()
  ((failed==0)) || { echo "ERROR: $fn failed" >&2; return 1; }
}

#!/usr/bin/env bash
# Record exact tool invocations; verbose tool output stays in per-step logs.
grid_run_logged(){
  local logfile=$1; shift
  mkdir -p "$(dirname -- "$logfile")"
  ( flock 9; printf '%q ' "$@" >&9; printf ' > %q 2>&1\n' "$logfile" >&9; ) 9>> "$GRID_COMMAND_FILE"
  [[ $GRID_DRY_RUN == FALSE ]] || return 0
  if "$@" > "$logfile" 2>&1; then return 0; else
    local rc=$?
    echo "ERROR: command failed (exit $rc): $logfile" >&2
    tail -n 12 "$logfile" >&2
    return "$rc"
  fi
}
grid_run(){ grid_run_logged "$work/log/step.$(date +%s%N).$BASHPID.log" "$@"; }
need(){ [[ -s $1 ]] || _grid_die "Missing/empty file: $1"; }
join_comma(){ local IFS=,; echo "$*"; }
# Publish with rename on the destination filesystem; never expose partial results.
publish(){
  local src=$1 dst=$2
  need "$src"; mkdir -p "$(dirname -- "$dst")"
  cp -- "$src" "$dst.tmp.$BASHPID"
  mv -f -- "$dst.tmp.$BASHPID" "$dst"
}
# A subset run must never overwrite a genome-wide result.
suffix=''
if [[ ${CHRS[*]} != "$(seq -s ' ' 1 22)" ]]; then suffix=".chr$(join_comma "${CHRS[@]}")"; fi

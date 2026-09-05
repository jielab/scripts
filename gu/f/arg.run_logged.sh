#!/usr/bin/env bash
# Persist the whole public invocation, including startup errors and worker output.
set -euo pipefail
export PYTHONUNBUFFERED=1
action=$1; script=$2; shift 2
dir_gen= arg_dir= grch=37
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case "${args[i]}" in
    --dir-gen) dir_gen=${args[i+1]:-};;
    --arg-dir) arg_dir=${args[i+1]:-};;
    --grch) grch=${args[i+1]:-37};;
  esac
done
case "${grch,,}" in 38|b38|grch38|hg38) grch=38;; *) grch=37;; esac
dir_gen=${dir_gen:-/mnt/d/data.BIG/refGen/1kg/$grch}
dir_gen=${dir_gen%/}
case "${dir_gen##*/}" in pfile|hap|typ|vcf|imp|gen) gen_root=$(dirname -- "$dir_gen");; *) gen_root=$dir_gen;; esac
log_dir=${arg_dir:-$gen_root/arg}/log
mkdir -p "$log_dir"
log_dir=$(cd -- "$log_dir" && pwd -P)
log=$(mktemp "$log_dir/arg.$action.$(date +%Y%m%dT%H%M%S).XXXXXX.log")
finish(){
  local status=$?
  printf '\n[ARG RUN END] time=%s exit_status=%s\n' "$(date -Is)" "$status" >> "$log"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
{
  printf '[ARG RUN START] time=%s pid=%s cwd=%q\n' "$(date -Is)" "$$" "$PWD"
  printf '[ARG COMMAND] ./arg.sh %q' "$action"
  printf ' %q' "$@"
  printf '\n'
} >> "$log"
printf 'Full run log: %s\n' "$log"
set +e
bash "$script" "$@" 2>&1 | (trap '' HUP; exec tee --output-error=warn-nopipe -a "$log")
statuses=("${PIPESTATUS[@]}")
status=${statuses[0]}
if ((status==0 && statuses[1]!=0)); then status=${statuses[1]}; fi
printf 'Run exit status: %s; log: %s\n' "$status" "$log"
exit "$status"

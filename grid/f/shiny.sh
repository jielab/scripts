#!/usr/bin/env bash
# The dashboard reads only the snapshot published by final.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
final_dir=/mnt/d/analysis/grid/final
host=${GRID_SHINY_HOST:-127.0.0.1}
port=${GRID_SHINY_PORT:-3839}
check=FALSE
dry_run=FALSE
while (($#)); do
  case $1 in --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;; esac
  case $1 in
    --help|-h|help)
      cat <<'HELP'
GRID Shiny — browse PCA, CSx, Disco, PRSformer and ABM results from final.
  ./grid.sh final
  ./grid.sh shiny
  ./grid.sh shiny --port 3839 --final-dir /mnt/d/analysis/grid/final

Options: --host ADDRESS (127.0.0.1), --port NUMBER (3839), --final-dir DIR,
         --check (validate snapshot and UI without starting a server), --dry-run.
First publish a snapshot with final; --summary-only summarizes existing results.
Refresh in the page after a new final run. Shiny never trains or reads raw scores.
Stop the local server with Ctrl+C. No packages are installed automatically.
HELP
      exit 0 ;;
    --check) check=TRUE; shift ;;
    --dry-run) dry_run=TRUE; shift ;;
    --host|--port|--final-dir)
      [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "Missing value: $1" >&2; exit 2; }
      case $1 in --host) host=$2 ;; --port) port=$2 ;; --final-dir) final_dir=$2 ;; esac
      shift 2 ;;
    *) echo "Unknown Shiny option: $1" >&2; exit 2 ;;
  esac
done
[[ $port =~ ^[0-9]{1,5}$ ]] && ((10#$port>0 && 10#$port<65536)) || { echo 'Port must be 1..65535' >&2; exit 2; }
port=$((10#$port))
if [[ $dry_run == TRUE ]]; then printf 'PLAN Shiny: %s at http://%s:%s\n' "$final_dir/dashboard.rds" "$host" "$port"; exit 0; fi
source "$ROOT/f/0.common.sh"
grid_activate_environment
grid_select_r shiny,bslib,DT,plotly,data.table,ggplot2,openxlsx
[[ -s $final_dir/dashboard.rds ]] || { echo 'Missing final/dashboard.rds; run ./grid.sh final or ./grid.sh final --summary-only first.' >&2; exit 1; }
export GRID_FINAL_DIR="$final_dir" GRID_SHINY_HOST="$host" GRID_SHINY_PORT="$port"
args=()
[[ $check == FALSE ]] || args+=(--check)
echo "GRID Shiny: http://$host:$port (data: $final_dir/dashboard.rds)"
exec "${GRID_R[@]}" "$ROOT/shiny/app.R" "${args[@]}"

#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GRID: cross-ancestry unique LD boundaries.

Usage:
  ./0.ld_uniq.sh run [options]

Modules:
  run                              Build requested chromosome boundaries and validate the merge.

Options:
  --submit hpc|local|auto          Execution machine [auto]
  --chr LIST                       Requested autosomes [1-22]
  --grch 37|38                     Genome build [38]
  --dir-out DIR                    Result directory [analysis/grid/38/ld; GRCh37: grid/37/ld]
  --dir-gen DIR                    chrN.pgen/pvar[.zst]/psam inputs
  --sample-panel FILE              Sample and super_pop columns
  --replace TRUE|FALSE             Rebuild requested results [FALSE]
  --plan-only                      Write commands and plan without starting work
  --jobs-N N                       Local concurrent workers [1]; ignored on HPC
  --wait TRUE|FALSE                Wait for HPC completion [FALSE]; local always waits
  --queue NAME                     CPU queue
  --job-hours N                    Maximum worker duration [120 h]
  --threads N                      PLINK threads and maximum DP state threads [8]
  --max-d-kb N                     Largest tested common D [80000]
  --d-grid-kb CSV                  Initial increasing D grid [250,500,1000,2000,5000,10000]
  -h, --help                       Show this help.

Examples:
  ./0.ld_uniq.sh run --grch 38 --chr 1 --plan-only
  ./0.ld_uniq.sh run --grch 38 --chr 1-22 --queue 33
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; fi
# CPU ancestry × chromosome LD partitions, followed by one guarded merge.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
usage() {
bash "${BASH_SOURCE[0]}" --help
}
for arg in "$@"; do
    case $arg in --help|-h) usage; exit 0 ;; esac
done
# Runtime defaults are owned by the shared environment table.
export LD_CONFIG_LOADED=1
if [[ ${BASH_SOURCE[0]} != "$0" ]]; then
    source "$ROOT/../0f/0os.sh"
    return
fi
plan_only=FALSE
while (($#)); do
    case $1 in
        --plan-only) plan_only=TRUE; shift ;;
        --cores|--chr|--grch|--submit|--dir-out|--dir-gen|--sample-panel|--jobs-N|--wait|--queue|--job-cores|--job-memory-gb|--job-hours|--threads|--replace|--max-d-kb|--d-grid-kb)
            [[ $# -ge 2 && -n $2 && $2 != -* ]] || { echo "$1 requires a value" >&2; exit 2; }
            case $1 in
                --cores) export WORKFLOW_CORES=$2 ;;
                --chr) export CHRS=$2 ;;
                --grch) export GRCH=$2 ;;
                --submit) export SUBMIT=$2 ;;
                --dir-out) export OUT=$2 ;;
                --dir-gen) export PFILE_DIR=$2 ;;
                --sample-panel) export PANEL=$2 ;;
                --jobs-N) export LD_JOBS=$2 ;;
                --wait) export SUBMIT_WAIT=$2 ;;
                --queue) export LD_QUEUE=$2 ;;
                --job-cores) export LD_SLOTS=$2 ;;
                --job-memory-gb) export LD_MEMORY_GIB=$2 ;;
                --job-hours) export LD_HOURS=$2 ;;
                --threads) export THREADS=$2 ;;
                --replace) export LD_REPLACE=$2 ;;
                --max-d-kb) export MAX_D_KB=$2 ;;
                --d-grid-kb) export D_GRID_KB=$2 ;;
            esac
            shift 2 ;;
        *) echo "Unknown option: $1; use --help for the current interface" >&2; exit 2 ;;
    esac
done
source "$ROOT/f/env.sh"
exec "$PYTHON" "$ROOT/f/plan.py" --plan-only "$plan_only"


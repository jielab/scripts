#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
LD partitions: one independent CPU job per ancestry and chromosome.

Usage: ./ld.sh [--chr 1|1-22|1,3,5] [options]
  --submit hpc|local|auto       Execution machine [auto]
  --chr LIST                    Requested autosomes [1-22]
  --grch 37|38                  Genome build [38]
  --dir-out DIR                 Result directory
  --dir-gen DIR                 chrN.pgen/pvar[.zst]/psam inputs
  --sample-panel FILE           Sample and super_pop columns
  --replace TRUE|FALSE          Rebuild requested results [FALSE]
                                 FALSE reuses complete partitions, repairs missing ones
  --plan-only                   Write commands and plan without starting work
  --jobs-N N                    Local concurrent workers [1]; ignored on HPC
  --wait TRUE|FALSE             Wait for HPC completion [FALSE]; local always waits
  --queue NAME                  CPU queue
  --job-cores N                 Requested CPU slots per worker [THREADS+1]
  --job-memory-gb N             Host RAM per worker [64 GiB]
  --job-hours N                 Maximum worker duration [120 h]
  --threads N                   PLINK threads and maximum DP state threads [8]
  --max-d-kb N                  Largest tested common D [80000]
  --d-grid-kb CSV               Initial increasing D grid [250,500,1000,2000,5000,10000]
  --help                       Show this help

Jobs are named ld38.AFR.chr1, ld38.EAS.chr1, ld38.EUR.chr1, etc. (ld37 for GRCh37). CPU only; no GPU flag.
HPC workers have no artificial concurrency dependencies; LSF schedules resources.
All four ancestries use the same explicit D grid. HPC merge validates outputs
after all workers terminate; publication requires every successful receipt and
joint boundary stability. An unstable
grid writes NEEDS_LARGER_GRID.json; increase --max-d-kb and rerun --replace FALSE.
There is no hidden LD calculation or cross-population autoextension in merge.
Full 1-22 runs publish ALL/AFR/EAS/EUR/SAS.<build>.bed; subsets use .SUBSET.bed.
Tool defaults are in ld.sh; environment variables can override them. No sibling 0f directory is needed.

Examples:
  ./ld.sh --chr 1 --plan-only
  mkdir -p "${ANALYSIS_ROOT:-$HOME/analysis}/ld/logs"
  nohup ./ld.sh --chr 1-22 --queue 33 --job-cores 40 > "${ANALYSIS_ROOT:-$HOME/analysis}/ld/logs/submit.log" 2>&1 &
PUBLIC_ENTRY_HELP
    exit 0
fi
# CPU ancestry × chromosome LD partitions, followed by one guarded merge.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
usage() {
bash "${BASH_SOURCE[0]}" --help
}
for arg in "$@"; do
    case $arg in --help|-h) usage; exit 0 ;; esac
done
# QiMing tool paths; environment variables may override these defaults.
export LD_PYTHON=${LD_PYTHON:-${DIR_WORK:-$HOME}/.conda/envs/py311/bin/python3}
export PLINK2=${PLINK2:-${DIR_WORK:-$HOME}/software/bin/plink2}
export LD_CXX=${LD_CXX:-${DIR_WORK:-$HOME}/.conda/envs/py311/bin/x86_64-conda-linux-gnu-c++}
export LD_ZSTD_PREFIX=${LD_ZSTD_PREFIX:-${DIR_WORK:-$HOME}/.conda/envs/py311}
export LD_QUEUE=${LD_QUEUE:-33}
export LD_CONFIG_LOADED=1
# Allow the installer to load the same tool defaults without starting jobs.
if [[ ${BASH_SOURCE[0]} != "$0" ]]; then return 0; fi
plan_only=FALSE
while (($#)); do
    case $1 in
        --plan-only) plan_only=TRUE; shift ;;
        --chr|--grch|--submit|--dir-out|--dir-gen|--sample-panel|--jobs-N|--wait|--queue|--job-cores|--job-memory-gb|--job-hours|--threads|--replace|--max-d-kb|--d-grid-kb)
            [[ $# -ge 2 && -n $2 && $2 != -* ]] || { echo "$1 requires a value" >&2; exit 2; }
            case $1 in
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

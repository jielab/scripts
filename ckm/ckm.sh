#!/usr/bin/env bash
# Standalone CKM entry; scheduling is delegated to the shared project planner.
set -euo pipefail
usage() {
cat <<'HELP'
CKM: clinical phenotype, continuous within-stage prognosis, and shared/stage-dependent factors.

Usage:
  ./1.phe.sh [options]
  ./ckm.sh [all|clock|assoc|marker|enrich|report|demo] [options]

Modules:
  all       clock, assoc, marker, enrich, report (requires ckm.rds)
  clock     grouped out-of-fold clinical/prot/met/joint prognostic refinement
  assoc     stage-stratified shared effects, interactions, lag analyses and Yang
  marker    cross-sectional molecular stage markers, separately labelled
  enrich    optional offline gene-set enrichment using local annotation and GMT
  report    refresh the aggregate-only interactive ckm.html
  demo      generate SYNTHETIC data under output/test, then exercise the full pipeline

Scientific options:
  --config FILE              Default config/ckm.json
  --out-dir DIR              Default $ANALYSIS_ROOT/ckm
  --phenotype-dir DIR        UKB root containing Rdata/ and rap/
  --stage3-method PREVENT|SCORE2
  --definition aha2026|aha2023
  --domains lifestyle,prot,met
  --threads N               R/association worker budget [4]
  --rscript FILE            Explicit Rscript interpreter
  --demo-n N                Synthetic participant count [6000]

Shared execution options (scripts/0f):
  --submit local|hpc|auto --queue NAME --job-cores N --job-memory-gb N
  --job-hours N --jobs-N N --wait TRUE|FALSE --plan-only --dir-jobs DIR
  wait TRUE|FALSE           Legacy spelling of --wait; wait for job completion
  -h, --help                Show this help without loading the runtime
  --dir-dat DIR --dir-big DIR

Examples:
  ./1.phe.sh --config config/ckm.json --submit local
  ./ckm.sh all --submit hpc --queue 33 --threads 8 --job-cores 8 --job-memory-gb 64
  ./ckm.sh clock --domains prot,met --submit local
  ./ckm.sh demo --out-dir /tmp/ckm-demo --submit local

Follow-up end is read from scripts/0f/phenotype.R::date_follow_end.
Optional verified coverage tables can shorten that date, never extend it.
SCORE2 is an adapted sensitivity definition. Clock position is not years to the next stage.
GPU is not needed. Install dependencies separately with install.sh.
HELP
}

for arg in "$@"; do
    case "$arg" in -h|--help|help) usage; exit 0 ;; esac
done
CKM_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ ! -f $CKM_ROOT/../0f/project_entry.sh ]]; then
    echo 'Place ckm/ beside the current scripts/0f/ directory.' >&2
    exit 2
fi
source "$CKM_ROOT/../0f/project_entry.sh"
# Normalize the supported legacy wait spelling before shared argument parsing.
ckm_args=()
while (($#)); do
    case "$1" in
        wait|--wait)
            [[ ${2:-} == TRUE || ${2:-} == FALSE ]] || { echo 'wait requires TRUE or FALSE' >&2; exit 2; }
            ckm_args+=(--wait "$2"); shift 2 ;;
        *) ckm_args+=("$1"); shift ;;
    esac
done
project_workflow_init "${ckm_args[@]}"
set -- "${WORKFLOW_ARGS[@]}"


module=all
config=$CKM_ROOT/config/ckm.json
out_dir=${CKM_OUT_DIR:-$ANALYSIS_ROOT/ckm}
threads=4
rscript=${CKM_RSCRIPT:-Rscript}
if [[ $submit == hpc && -z ${CKM_RSCRIPT:-} ]]; then
    rscript=$DIR_WORK/.conda/envs/ckm/bin/Rscript
fi
extra=()
if [[ $# -gt 0 && $1 != -* ]]; then module=$1; shift; fi
while (($#)); do
    case "$1" in
        -h|--help|help) usage; exit 0 ;;
        --config) config=${2:?--config requires a file}; shift 2 ;;
        --out-dir|--output-dir|--analysis-root) out_dir=${2:?--out-dir requires a directory}; shift 2 ;;
        --threads) threads=${2:?--threads requires N}; shift 2 ;;
        --rscript) rscript=${2:?--rscript requires an interpreter}; shift 2 ;;
        --phenotype-dir|--stage3-method|--definition|--domains|--demo-n)
            extra+=("$1" "${2:?option requires a value}"); shift 2 ;;
        *) echo "Unknown CKM option: $1" >&2; exit 2 ;;
    esac
done
[[ $module =~ ^(phe|all|clock|assoc|marker|enrich|report|demo)$ ]] || { echo "Unknown CKM module: $module" >&2; exit 2; }
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo '--threads must be a positive integer' >&2; exit 2; }
mkdir -p -- "$out_dir/tmp"
out_dir=$(cd -- "$out_dir" && pwd -P)
export CKM_THREADS=$threads OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=$threads TMPDIR=$out_dir/tmp PYTHONDONTWRITEBYTECODE=1
if [[ -z ${CKM_PYTHON:-} && -x $(dirname -- "$rscript")/python3 ]]; then
    export CKM_PYTHON=$(dirname -- "$rscript")/python3
fi
export CKM_PYTHON=${CKM_PYTHON:-${SHARED_PYTHON:-python3}}
command -v "$rscript" >/dev/null || { echo "Rscript not found: $rscript. Run install.sh or use --rscript." >&2; exit 2; }
exec "$rscript" "$CKM_ROOT/f/run.R" --module "$module" --config "$config" --out-dir "$out_dir" --threads "$threads" "${extra[@]}"

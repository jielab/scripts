#!/usr/bin/env bash
set -e



usage() {
  cat <<'HELP'
Usage:
  ./move.sh [options]

Main options:
  -Y, --trait TRAIT          One disease outcome
      --traits CSV           Multiple outcomes, e.g. mi,stroke,depress,acd
                              Equivalent environment variable: Ys=mi,stroke,depress,acd
  -s, --start-step STEP      Start from this step. Default: all
      --transition MODE      Main transition: birth_current or home_future. Default: birth_current
      --no-prot              Skip proteomics scans
      --no-met               Skip metabolomics scans
      --max-prot N           Limit proteomic features for testing/debugging. 0 = all
      --max-met N            Limit metabolite features for testing/debugging. 0 = all
      --resume-omics         Reuse completed protein/metabolite omics layer RDS files
  -h, --help                 Show this message

Step order:
  data_prep       Build rich/poor TDI transition, individual SES, and move_core.rds
  qc              Root-level job/home/move/omics QC plots and tables
  ses_disease     Disease-specific transition x individual SES x risk analyses
  omics_signature Disease-specific proteomic/metabolomic signatures
  deep_dive       Biomarker-score, origin/destination, and attenuation analyses
  consolidate     Candidate driver evidence matrix and run index

Output structure:
  D:/analysis/move/Fig1-Fig2.* and raw/               # disease-free/root-level results
  D:/analysis/move/[Y]/Fig3-Fig7.[Y].* and [Y]/raw/   # disease-specific results

Examples:
  ./move.sh
  Ys=mi,stroke,depress,acd ./move.sh
  ./move.sh --traits mi,stroke,depress,acd --transition birth_current
  START_STEP=omics_signature ./move.sh --traits mi,stroke --max-prot 200
HELP
}

START_STEP="${START_STEP:-all}"
Y_CSV="${Ys:-${Y:-mi,stroke,depress,acd}}"
MOVE_TRANSITION="${MOVE_TRANSITION:-birth_current}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -Y|--trait) Y_CSV="$2"; shift 2 ;;
    --traits|--Ys|--ys) Y_CSV="$2"; shift 2 ;;
    -s|--start-step|--start_step) START_STEP="$2"; shift 2 ;;
    --transition) MOVE_TRANSITION="$2"; shift 2 ;;
    --no-prot) export PROT_DO=FALSE; shift ;;
    --no-met) export MET_DO=FALSE; shift ;;
    --max-prot) export MOVE_MAX_PROT_FEATURES="$2"; shift 2 ;;
    --max-met) export MOVE_MAX_MET_FEATURES="$2"; shift 2 ;;
    --resume-omics) export MOVE_RESUME_OMICS=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

export START_STEP MOVE_TRANSITION
export DIR0="${DIR0:-}"
export PROT_DO="${PROT_DO:-TRUE}"
export MET_DO="${MET_DO:-TRUE}"
export MOVE_RESUME_OMICS="${MOVE_RESUME_OMICS:-FALSE}"
export N_CORES="${N_CORES:-4}"
export SEED="${SEED:-2026}"

move_privacy_root="${MOVE_OUTDIR:-/mnt/d/analysis/move}"
privacy_cleanup_completed_run=FALSE
cleanup_ukb_individual_outputs() {
  local status=$? f
  trap - EXIT
  if (( status != 0 )) || [[ "$privacy_cleanup_completed_run" != TRUE ]]; then
    echo "MOVE did not reach full completion (exit=$status); preserving participant-level working data for resume." >&2
    exit "$status"
  fi
  case "$move_privacy_root" in ""|/|/mnt|/mnt/|/mnt/d|/mnt/d/|/mnt/d/analysis|/mnt/d/analysis/)
    echo "PRIVACY ERROR: unsafe MOVE_OUTDIR: '$move_privacy_root'" >&2; exit 70 ;;
  esac
  if [[ -d "$move_privacy_root" ]]; then
    echo "Privacy cleanup: removing MOVE UKB participant-level working data"
    while IFS= read -r -d '' f; do rm -f -- "$f"; done < <(find "$move_privacy_root" -type f \( \
      -name 'move_core*.rds' -o -name 'move_core*.tsv.gz' -o \
      -name 'move_home_long.rds' \) -print0)
  fi
  exit "$status"
}
trap cleanup_ukb_individual_outputs EXIT

RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript || true)}"
if [[ -z "${RSCRIPT_BIN}" ]]; then
  echo "Error: cannot find Rscript. Set RSCRIPT_BIN=/path/to/Rscript." >&2
  exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_FILE="${SCRIPT_DIR}/move.v6.R"
if [[ ! -s "${R_FILE}" ]]; then
  R_FILE="${SCRIPT_DIR}/move.v5.R"
fi
if [[ ! -s "${R_FILE}" ]]; then
  R_FILE="${SCRIPT_DIR}/move.R"
fi
if [[ ! -s "${R_FILE}" ]]; then
  echo "Error: cannot find move.v6.R, move.v5.R, or move.R in ${SCRIPT_DIR}" >&2
  exit 1
fi

IFS=',' read -ra TRAITS <<< "${Y_CSV}"

echo "Move v6 pipeline"
echo "Ys=${Y_CSV}"
echo "START_STEP=${START_STEP}"
echo "MOVE_TRANSITION=${MOVE_TRANSITION}"
echo "DIR0=${DIR0:-auto}"
echo "MOVE_OUTDIR=${MOVE_OUTDIR:-auto}"
echo "PROT_DO=${PROT_DO}; MET_DO=${MET_DO}"
echo "MOVE_RESUME_OMICS=${MOVE_RESUME_OMICS}"
echo "Rscript=${RSCRIPT_BIN}"
echo "R file=${R_FILE}"

for trait0 in "${TRAITS[@]}"; do
  trait=$(echo "$trait0" | xargs)
  [[ -z "$trait" ]] && continue
  export Y="$trait"
  unset MOVE_TRAIT_OUTDIR
  echo
  echo "#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "# RUN TRAIT: ${Y}"
  echo "#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  "${RSCRIPT_BIN}" "${R_FILE}" \
    --start_step "${START_STEP}" \
    --trait "${Y}" \
    --transition "${MOVE_TRANSITION}"
done

privacy_cleanup_completed_run=TRUE
echo
echo "Done. Outputs are under MOVE_OUTDIR (default D:/analysis/move or /mnt/d/analysis/move)."

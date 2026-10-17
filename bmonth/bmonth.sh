#!/usr/bin/env bash
set -euo pipefail
ORIGINAL_ARGS=("$@")
if [[ ${#ORIGINAL_ARGS[@]} -eq 0 ]]; then ORIGINAL_ARGS=(all); fi

usage(){ cat <<'HELP'
Birth-month / birth-environment analysis in UK Biobank

Examples (one line at a time; continue only after success):
  cd /mnt/d/scripts/bmonth
  ./bmonth.sh check
  ./bmonth.sh prep
  ./bmonth.sh phewas
  ./bmonth.sh assoc
  ./bmonth.sh decomp
  ./bmonth.sh med
  ./bmonth.sh final

All commands run in the foreground. The 'all' module runs the sequence above.

Optional temperature analysis (requires HadUK-Grid in --temp-dir):
  ./bmonth.sh temp
  ./bmonth.sh tempmed
  ./bmonth.sh final

The 'plot' module redraws figures from saved results.

Options:
  --out DIR                  output directory [ /mnt/d/analysis/bmonth ]
  --in DIR                   phenotype directory containing Rdata/ [ /mnt/d/data/ukb ]
  --rscript FILE             Rscript binary
  --outcomes CSV             incident outcomes [cvd_cad,ra]
  --exposures CSV            PRIMARY birth-month scalar exposures
                             [bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn]
  --context-exposures CSV    birthplace-context PheWAS exposures
                             [birth_north_100km,birth_east_100km]
  --exploratory-traits CSV   exploratory birth-order exposures; 'none' disables
                             [border12,border1n,border_rank,border_relative]
  --phewas-cores N|auto      PheWAS workers [auto: up to 8, memory-limited]
  --med-boot N               mediation/bootstrap replicates [100]
  --follow-end DATE          override reconstructed-outcome end
  --temp-dir DIR             HadUK-Grid NetCDF directory [/mnt/d/data/ukb/temperature]
  --temp-res RES             HadUK-Grid resolution [1km]
  --tempmed-phenos CSV       explicit phecodes for temperature decomposition
  --replace                  recompute scalar and joint PheWAS caches

Design:
  * PRIMARY birth-month models adjust age, sex, assessment centre, education,
    adult TDI, birth north/east and tdi.birth.
  * Birthplace-context PheWAS restores birth_north_100km and birth_east_100km
    as exposures using the OLD-compatible D2 model (age, sex, centre,
    education, adult TDI). These are context/supplementary analyses, not
    birth-month effects.
  * The decomposition stage then uses a common sample to quantify attenuation:
      D1 core
      D2 + adult SES
      D3 + birth geography
      D4 + birthplace TDI proxy
      D5 + migration distance
  * Temperature is NOT a baseline PheWAS covariate. It is handled only in
    the optional temp/tempmed modules.
  * move is post-birth and is used in D5 sensitivity / interaction analyses,
    not in the primary model.
  * Exploratory birth-order traits use age, sex, centre, education and adult TDI
    for PheWAS, Cox and omics mediation; D1-D5 assess contextual attenuation.
    Their eligibility rules apply only to birth-order traits.
HELP
}

STEP=all
OUTDIR=/mnt/d/analysis/bmonth
INDIR=/mnt/d/data/ukb
RSCRIPT_BIN=""
OUTCOMES="cvd_cad,ra"
EXPOSURES="bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn"
CONTEXT_EXPOSURES="birth_north_100km,birth_east_100km"
EXPLORATORY_TRAITS="border12,border1n,border_rank,border_relative"
PHEWAS_CORES=auto
MED_BOOT=100
FOLLOW_END=""
TEMP_DIR="${BMONTH_TEMP_DIR:-/mnt/d/data/ukb/temperature}"
TEMP_RES="1km"
TEMPMED_PHENOS=""
REPLACE=FALSE

while [[ $# -gt 0 ]]; do
  case "$1" in
    prep|assoc|phewas|decomp|med|temp|tempmed|final|plot|all|check) STEP="$1"; shift ;;
    --out) OUTDIR="${2:?}"; shift 2 ;;
    --in) INDIR="${2:?}"; shift 2 ;;
    --rscript) RSCRIPT_BIN="${2:?}"; shift 2 ;;
    --outcomes) OUTCOMES="${2:?}"; shift 2 ;;
    --exposures) EXPOSURES="${2:?}"; shift 2 ;;
    --context-exposures) CONTEXT_EXPOSURES="${2:?}"; shift 2 ;;
    --exploratory-traits) EXPLORATORY_TRAITS="${2:?}"; shift 2 ;;
    --phewas-cores) PHEWAS_CORES="${2:?}"; shift 2 ;;
    --med-boot) MED_BOOT="${2:?}"; shift 2 ;;
    --follow-end) FOLLOW_END="${2:?}"; shift 2 ;;
    --temp-dir) TEMP_DIR="${2:?}"; shift 2 ;;
    --temp-res) TEMP_RES="${2:?}"; shift 2 ;;
    --tempmed-phenos) TEMPMED_PHENOS="${2:?}"; shift 2 ;;
    --replace) REPLACE=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$PHEWAS_CORES" == auto || "$PHEWAS_CORES" =~ ^[1-9][0-9]*$ ]] || {
  echo 'ERROR: --phewas-cores must be auto or a positive integer' >&2; exit 2; }
[[ "$MED_BOOT" =~ ^(0|[1-9][0-9]*)$ ]] || {
  echo 'ERROR: --med-boot must be a nonnegative integer' >&2; exit 2; }

RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: Rscript not found" >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"

export BMONTH_OUTDIR="$OUTDIR" BMONTH_INDIR="$INDIR" BMONTH_OUTCOMES="$OUTCOMES"
export BMONTH_EXPOSURES="$EXPOSURES" BMONTH_CONTEXT_EXPOSURES="$CONTEXT_EXPOSURES"
export BMONTH_EXPLORATORY_TRAITS="$EXPLORATORY_TRAITS"
export BMONTH_PHEWAS_CORES="$PHEWAS_CORES"
export BMONTH_MED_BOOT="$MED_BOOT" BMONTH_FOLLOW_END="$FOLLOW_END"
export BMONTH_TEMP_DIR="$TEMP_DIR" BMONTH_TEMP_RES="$TEMP_RES"
export BMONTH_TEMPMED_PHENOS="$TEMPMED_PHENOS"
export BMONTH_REPLACE="$REPLACE" BMONTH_FDIR="$SCRIPT_DIR/f"
export SCRIPT_LOG_DIR="${SCRIPT_LOG_DIR:-$OUTDIR/logs}"
source /mnt/d/scripts/0f/console.sh "${ORIGINAL_ARGS[@]}" 2>/dev/null || true

[[ -f "$SCRIPT_DIR/../0f/memory_cap.sh" ]] || {
  echo 'ERROR: missing memory guard; refusing to run uncapped.' >&2; exit 1;
}
source "$SCRIPT_DIR/../0f/memory_cap.sh"
memory_cap_enter bmonth 24 0 "$SCRIPT_DIR/bmonth.sh" "${ORIGINAL_ARGS[@]}"

run_r(){
  local s="$1"
  echo "[bmonth] START ${s##*/}"
  if [[ ${s##*/} == 3.phewas.R ]]; then
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
      VECLIB_MAXIMUM_THREADS=1 "$RSCRIPT_BIN" --vanilla "$s"
  else
    "$RSCRIPT_BIN" --vanilla "$s"
  fi
  echo "[bmonth] DONE ${s##*/}"
}

case "$STEP" in
  check) run_r "$SCRIPT_DIR/f/00.check.R" ;;
  prep) run_r "$SCRIPT_DIR/f/1.prep.R" ;;
  assoc) run_r "$SCRIPT_DIR/f/2.assoc.R" ;;
  phewas) run_r "$SCRIPT_DIR/f/3.phewas.R" ;;
  decomp) run_r "$SCRIPT_DIR/f/5.decomp.R" ;;
  med) run_r "$SCRIPT_DIR/f/4.med.R" ;;
  temp) run_r "$SCRIPT_DIR/f/6.temperature.R" ;;
  tempmed) run_r "$SCRIPT_DIR/f/7.tempmed.R" ;;
  final) run_r "$SCRIPT_DIR/f/5.final.R" ;;
  plot) run_r "$SCRIPT_DIR/f/6.replot.R" ;;
  all)
    run_r "$SCRIPT_DIR/f/00.check.R"
    run_r "$SCRIPT_DIR/f/1.prep.R"
    run_r "$SCRIPT_DIR/f/3.phewas.R"
    run_r "$SCRIPT_DIR/f/2.assoc.R"
    run_r "$SCRIPT_DIR/f/5.decomp.R"
    run_r "$SCRIPT_DIR/f/4.med.R"
    run_r "$SCRIPT_DIR/f/5.final.R"
    ;;
esac

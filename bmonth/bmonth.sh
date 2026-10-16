#!/usr/bin/env bash
set -euo pipefail
ORIGINAL_ARGS=("$@")
if [[ ${#ORIGINAL_ARGS[@]} -eq 0 ]]; then ORIGINAL_ARGS=(all); fi

usage(){ cat <<'HELP'
Birth-month / birth-environment analysis in UK Biobank

cd /mnt/d/scripts/bmonth
./bmonth.sh all
./bmonth.sh prep
./bmonth.sh phewas
./bmonth.sh assoc
./bmonth.sh decomp
./bmonth.sh med
./bmonth.sh temp --temp-dir /path/to/HadUK-Grid
./bmonth.sh tempmed
./bmonth.sh final
./bmonth.sh plot
./bmonth.sh check

Options:
  --out DIR              output directory [ /mnt/d/analysis/bmonth ]
  --in DIR               phenotype directory containing Rdata/ [ /mnt/d/data/ukb/phe ]
  --rscript FILE         Rscript binary
  --outcomes CSV         incident outcomes [cvd_cad,ra]
  --exposures CSV        scalar exposures
                         [bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn]
  --phewas-cores N|auto  PheWAS workers [auto: up to 8, memory-limited]
  --med-boot N           mediation/bootstrap replicates [100]
  --follow-end DATE      override reconstructed-outcome end
  --temp-dir DIR         HadUK-Grid NetCDF directory [/mnt/d/data/ukb/temperature]
  --temp-res RES         HadUK-Grid resolution [1km]
  --tempmed-phenos CSV   explicit phecodes for temperature decomposition
  --replace              recompute scalar and joint PheWAS caches

Notes:
  * Primary birth-month models adjust birth north/east and tdi.birth in addition
    to age, sex, education, adult TDI and assessment centre.
  * move (birth-to-current distance from move.rds) is post-birth and is used in
    D4 sensitivity / interaction analyses, not in the primary model.
  * `temp` is optional and expects externally downloaded HadUK-Grid files.
HELP
}

STEP=all
OUTDIR=/mnt/d/analysis/bmonth
INDIR=/mnt/d/data/ukb/phe
RSCRIPT_BIN=""
OUTCOMES="cvd_cad,ra"
EXPOSURES="bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn"
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

[[ "$PHEWAS_CORES" == auto || "$PHEWAS_CORES" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --phewas-cores must be auto or a positive integer' >&2; exit 2; }
[[ "$MED_BOOT" =~ ^(0|[1-9][0-9]*)$ ]] || { echo 'ERROR: --med-boot must be a nonnegative integer' >&2; exit 2; }

RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: Rscript not found" >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"

export BMONTH_OUTDIR="$OUTDIR" BMONTH_INDIR="$INDIR" BMONTH_OUTCOMES="$OUTCOMES"
export BMONTH_EXPOSURES="$EXPOSURES" BMONTH_PHEWAS_CORES="$PHEWAS_CORES"
export BMONTH_MED_BOOT="$MED_BOOT" BMONTH_FOLLOW_END="$FOLLOW_END"
export BMONTH_TEMP_DIR="$TEMP_DIR" BMONTH_TEMP_RES="$TEMP_RES"
export BMONTH_TEMPMED_PHENOS="$TEMPMED_PHENOS"
export BMONTH_REPLACE="$REPLACE" BMONTH_FDIR="$SCRIPT_DIR/f"
export SCRIPT_LOG_DIR="${SCRIPT_LOG_DIR:-$OUTDIR/logs}"
source /mnt/d/scripts/0f/console.sh "${ORIGINAL_ARGS[@]}" 2>/dev/null || true

# Every stage and worker must stay inside the existing 24 GiB, zero-swap cap.
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

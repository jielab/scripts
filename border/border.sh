#!/usr/bin/env bash
set -euo pipefail
ORIGINAL_ARGS=("$@")

usage() {
cat <<'HELP'
Birth-order analysis in UK Biobank

cd /mnt/d/scripts/border
./border.sh all
./border.sh prep
./border.sh assoc
./border.sh phewas
./border.sh med
./border.sh final
./border.sh plot          # redraw saved PheWAS results only
./border.sh check

Options:
  --out DIR              output directory [ /mnt/d/analysis/border ]
  --in DIR               phenotype directory containing Rdata/ [ /mnt/d/data/ukb/phe ]
  --rscript FILE         Rscript binary
  --outcomes CSV         main outcomes [cvd_cad,ra]
  --exposures CSV        birth-order exposures [border12,border1n,border_rank,border_relative]
  --phewas-cores N       PheWAS cores [4]
  --med-boot N           mediation bootstrap replicates [100]
  --follow-end DATE      override reconstructed-outcome end [date_follow_end in 0f/phenotype.R]
  --replace              recompute matching PheWAS caches
  check                  check dependencies, inputs and outcome definitions only
HELP
}

STEP=all
OUTDIR=/mnt/d/analysis/border
INDIR=/mnt/d/data/ukb/phe
RSCRIPT_BIN=""
OUTCOMES="cvd_cad,ra"
EXPOSURES="border12,border1n,border_rank,border_relative"
PHEWAS_CORES=4
MED_BOOT=100
# Empty means read the shared phenotype.R date during prep, or the saved cohort
# setting in later stages. An explicit option is checked against that cohort.
FOLLOW_END=""
REPLACE=FALSE

while [[ $# -gt 0 ]]; do
  case "$1" in
    prep|assoc|phewas|med|final|plot|all|check) STEP="$1"; shift ;;
    --out) OUTDIR="${2:?}"; shift 2 ;;
    --in) INDIR="${2:?}"; shift 2 ;;
    --rscript) RSCRIPT_BIN="${2:?}"; shift 2 ;;
    --outcomes) OUTCOMES="${2:?}"; shift 2 ;;
    --exposures) EXPOSURES="${2:?}"; shift 2 ;;
    --phewas-cores) PHEWAS_CORES="${2:?}"; shift 2 ;;
    --med-boot) MED_BOOT="${2:?}"; shift 2 ;;
    --follow-end) FOLLOW_END="${2:?}"; shift 2 ;;
    --replace) REPLACE=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$PHEWAS_CORES" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --phewas-cores must be a positive integer' >&2; exit 2; }
[[ "$MED_BOOT" =~ ^(0|[1-9][0-9]*)$ ]] || { echo 'ERROR: --med-boot must be a nonnegative integer' >&2; exit 2; }
RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: Rscript not found" >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"

export BORDER_OUTDIR="$OUTDIR"
export BORDER_INDIR="$INDIR"
export BORDER_OUTCOMES="$OUTCOMES"
export BORDER_EXPOSURES="$EXPOSURES"
export BORDER_PHEWAS_CORES="$PHEWAS_CORES"
export BORDER_MED_BOOT="$MED_BOOT"
export BORDER_FOLLOW_END="$FOLLOW_END"
export BORDER_REPLACE="$REPLACE"
export BORDER_FDIR="$SCRIPT_DIR/f"
export SCRIPT_LOG_DIR="${SCRIPT_LOG_DIR:-$OUTDIR/logs}"
# Parse --out before starting the shared logger so isolated tests and alternate
# runs also keep their logs inside the requested output directory.
# shellcheck source=/dev/null
source /mnt/d/scripts/0f/console.sh "${ORIGINAL_ARGS[@]}" 2>/dev/null || true

STAGE_INDEX=0
STAGE_TOTAL=1
[[ "$STEP" != all ]] || STAGE_TOTAL=6
run_r() {
  local script="$1"; shift || true
  local label start=$SECONDS rc
  case "${script##*/}" in
    00.check.R) label="检查依赖和输入" ;;
    1.prep.R) label="准备分析队列" ;;
    3.phewas.R) label="PheWAS：关联分析和四张图" ;;
    2.assoc.R) label="主疾病 Cox 关联" ;;
    4.med.R) label="中介分析：筛选和 bootstrap" ;;
    5.final.R) label="汇总结果和最终绘图" ;;
    6.replot.R) label="从保存结果重画 PheWAS" ;;
    *) label="${script##*/}" ;;
  esac
  STAGE_INDEX=$((STAGE_INDEX + 1))
  echo "[border] START [$STAGE_INDEX/$STAGE_TOTAL] $label"
  if "$RSCRIPT_BIN" --vanilla "$script" "$@"; then
    echo "[border] DONE [$STAGE_INDEX/$STAGE_TOTAL] $label；耗时 $((SECONDS - start)) 秒"
  else
    rc=$?
    echo "[border] FAIL [$STAGE_INDEX/$STAGE_TOTAL] $label；退出码 $rc" >&2
    return "$rc"
  fi
}

case "$STEP" in
  check)  run_r "$SCRIPT_DIR/f/00.check.R" ;;
  prep)   run_r "$SCRIPT_DIR/f/1.prep.R" ;;
  assoc)  run_r "$SCRIPT_DIR/f/2.assoc.R" ;;
  phewas) run_r "$SCRIPT_DIR/f/3.phewas.R" ;;
  med)    run_r "$SCRIPT_DIR/f/4.med.R" ;;
  plot)   run_r "$SCRIPT_DIR/f/6.replot.R" ;;
  final)  run_r "$SCRIPT_DIR/f/5.final.R" ;;
  all)
    run_r "$SCRIPT_DIR/f/00.check.R" packages
    run_r "$SCRIPT_DIR/f/1.prep.R"
    run_r "$SCRIPT_DIR/f/3.phewas.R"
    run_r "$SCRIPT_DIR/f/2.assoc.R"
    run_r "$SCRIPT_DIR/f/4.med.R"
    run_r "$SCRIPT_DIR/f/5.final.R"
    ;;
esac

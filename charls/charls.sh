#!/usr/bin/env bash
set -euo pipefail



usage() {
  cat <<'EOF_HELP'
Usage:
  ./charls.sh [MODE]

Default:
  ./charls.sh
      Run the full CHARLS AIPFI-CMM reproduction workflow.

Modes:
  all      Run prep -> Fig1 -> Tab1 -> Tab2 -> Tab3 -> Fig2 -> Fig3 -> Fig4 -> Fig5 -> Fig6 -> supp.
  clean    Remove previous outputs under CHARLS_OUT.
  prep     Build analytic datasets from CHARLS files and write QC/variable-map files.
  fig1     Participant flow diagram and Fig1.out.xlsx.
  tab1     Baseline characteristics by baseline AIPFI quartiles.
  tab2     Cox models for baseline AIPFI and CMM incidence.
  tab3     Cox models for cumulative AIPFI and AIPFI trajectory clusters.
  fig2     Kaplan-Meier curves by baseline AIPFI quartiles and trajectory clusters.
  fig3     K-means trajectory figure and cluster summaries.
  fig4     Restricted cubic-spline / natural-spline Cox curves.
  fig5     Subgroup forest plots for baseline AIPFI and cumulative AIPFI.
  fig6     Clinical prediction model: ROC, calibration, DCA, and SHAP-like explanation.
  supp     Supplementary sensitivity analyses and supplementary figures.
  qc       Print expected output paths and key QC files.

Help:
  ./charls.sh -h
  ./charls.sh --help

Main examples:
  ./charls.sh
  ./charls.sh prep
  START_STEP=tab2 ./charls.sh
  CHARLS_CLEAN=/mnt/d/data/charls/clean CHARLS_OUT=/mnt/d/analysis/charls ./charls.sh

Common variables:
  CHARLS_CLEAN=/mnt/d/data/charls/clean
  CHARLS_OUT=/mnt/d/analysis/charls
  CHARLS_HARMONIZED=/mnt/d/data/charls/clean/Harmonized/H_CHARLS_D_Data.dta
  CHARLS_MANIFEST=/mnt/d/data/charls/clean/_extract_qc/extracted_manifest.tsv   Optional; used for manifest QC only.
  USE_MICE=1              Use mice imputation if the package is installed; otherwise median/mode fallback.
  USE_HARMONIZED=1        Prefer the Harmonized CHARLS file as the primary source.
  FAIL_ON_VAR_MISS=1      Stop if a required variable cannot be resolved.
  ALLOW_HIGH_MISSING=0    Stop if key extracted variables have >20% missingness.
  SEED=2026
  CHARLS_VAR_W1_TG=...      Optional manual variable override if auto-matching fails.
  CHARLS_VAR_W1_HDLC=...    Optional manual variable override if auto-matching fails.

Input assumptions:
  - Data are already extracted under CHARLS_CLEAN by wave folders: 2011, 2013, 2015, 2018, 2020, Harmonized.
  - The workflow uses Harmonized CHARLS as primary data and merges raw Blood/Biomarker files only through explicit ID keys.
  - Outputs are written to CHARLS_OUT, corresponding to D:\analysis\charls on Windows.

Notes:
  - Use environment variables for settings; do not pass settings as --options.
  - Unknown --options are rejected intentionally to keep the interface simple.
EOF_HELP
}

for arg in "$@"; do
  case "$arg" in
    -h|--h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "Unknown option: $arg" >&2
      usage >&2
      exit 1
      ;;
  esac
done
if [[ $# -gt 1 ]]; then
  echo "ERROR: charls.sh accepts at most one MODE, but got: $*" >&2
  usage >&2
  exit 1
fi

MODE=${1:-all}
case "$MODE" in
  all|clean|prep|fig1|tab1|tab2|tab3|fig2|fig3|fig4|fig5|fig6|supp|qc) ;;
  *) echo "ERROR invalid MODE=$MODE" >&2; usage >&2; exit 1 ;;
esac

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths and runtime settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DIRSCRIPT=${DIRSCRIPT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
CHARLS_CLEAN=${CHARLS_CLEAN:-/mnt/d/data/charls/clean}
CHARLS_OUT=${CHARLS_OUT:-/mnt/d/analysis/charls}
CHARLS_HARMONIZED=${CHARLS_HARMONIZED:-$CHARLS_CLEAN/Harmonized/H_CHARLS_D_Data.dta}
SEED=${SEED:-2026}
USE_MICE=${USE_MICE:-1}
USE_HARMONIZED=${USE_HARMONIZED:-1}
FAIL_ON_VAR_MISS=${FAIL_ON_VAR_MISS:-1}
ALLOW_HIGH_MISSING=${ALLOW_HIGH_MISSING:-0}
CHARLS_MANIFEST=${CHARLS_MANIFEST:-}
CHARLS_VERSION=${CHARLS_VERSION:-v008}
START_STEP=${START_STEP:-}

R_SCRIPT=${R_SCRIPT:-$DIRSCRIPT/charls.R}
[[ -s "$R_SCRIPT" ]] || { echo "ERROR missing R script: $R_SCRIPT" >&2; exit 1; }
[[ -d "$CHARLS_CLEAN" ]] || { echo "ERROR missing CHARLS_CLEAN directory: $CHARLS_CLEAN" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "ERROR Rscript not found in PATH" >&2; exit 1; }

mkdir -p "$CHARLS_OUT" "$CHARLS_OUT/log" "$CHARLS_OUT/dat" "$CHARLS_OUT/qc"

export CHARLS_CLEAN CHARLS_OUT CHARLS_HARMONIZED CHARLS_MANIFEST SEED USE_MICE USE_HARMONIZED FAIL_ON_VAR_MISS ALLOW_HIGH_MISSING START_STEP

if [[ "$MODE" == clean ]]; then
  echo "Removing previous outputs under $CHARLS_OUT"
  find "$CHARLS_OUT" -mindepth 1 -maxdepth 1 ! -name log -exec rm -rf {} +
  mkdir -p "$CHARLS_OUT/log" "$CHARLS_OUT/dat" "$CHARLS_OUT/qc"
  exit 0
fi

echo "CHARLS_CLEAN       = $CHARLS_CLEAN"
echo "CHARLS_HARMONIZED  = $CHARLS_HARMONIZED"
echo "CHARLS_OUT         = $CHARLS_OUT"
echo "CHARLS_MANIFEST    = ${CHARLS_MANIFEST:-<auto/optional>}"
echo "CHARLS_VERSION     = $CHARLS_VERSION"
echo "MODE               = $MODE"
echo "START_STEP         = ${START_STEP:-<none>}"
echo

Rscript "$R_SCRIPT" "$MODE" 2>&1 | tee "$CHARLS_OUT/log/${MODE}.log"

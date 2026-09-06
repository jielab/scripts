#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./maha.sh [ukb|nhanes|chns|all] [--steps <internal-step-list>]
  ./maha.sh final [ukb|nhanes|chns|all]
  ./maha.sh --help

Examples:
  ./maha.sh all
  ./maha.sh                    # same as: ./maha.sh all

  ./maha.sh ukb
  ./maha.sh nhanes
  ./maha.sh chns

  ./maha.sh ukb --steps data_qc,fig2,fig3
  ./maha.sh final   # rebuild all manuscript figures from cached cohort outputs

Modules:
  ukb, nhanes, chns       Run cohort analyses, then assemble final outputs.
  all                    Run all cohort analyses, then assemble final outputs.
  final [COHORT]         Integrate cached analyses into publication-ready figures
                         and tables for Shiny/manuscripts. Default scope: all.

Options:
  --rscript FILE          Rscript executable.
  --output-dir DIR        Publication output root. Default: /mnt/d/analysis/maha.
  --jobs N                Parallel matched-null workers. Default: 4.
  --null-max N            0 uses all direction-null scores. Default: 0.
  --helper-dir DIR        Shared R helper directory. Default: /mnt/d/scripts/0f.
  --replace-steps CSV     Completed cached steps to rerun, or all.
  --tile-proxy URL        Candidate HTTP proxy for OSM tiles.
  --chns-project-dir DIR  CHNS project root. Default: /mnt/d/data/chns.
  --chns-food-map FILE
  --chns-nutrient-map FILE
  --chns-baseline-wave YEAR
  --chns-followup-end YEAR
  --chns-min-age N
  --chns-min-diet-days N
  --chns-run-sweep 0|1        Run the focused high-event sensitivity. Default: 1.
  --chns-sweep-waves CSV      Diet waves for the cumulative-average sensitivity.
                              Default: 2004,2006,2009,2011.
  --chns-model-age N          Age threshold for Models 4-5. Default: 40.
  --chns-event-target N       Descriptive event-count flag only. Default: 500.

Output layout:
  /mnt/d/analysis/maha/
      Fig1.ukb.phewas.png
      Fig2.ukb.precision.png
      Fig3.nhanes.mortality.png
      Fig4.chns.mortality.png
      Fig5.ukb.geography.png
      FigS*.{cohort}.*.png
      matching *.out.xlsx files
      ukb/       <- UKB intermediate png/xlsx/log/tsv/rds
      nhanes/    <- NHANES intermediate files
      chns/      <- CHNS intermediate files

USAGE
}

STEP="all"
STEP_SET=0
INTERNAL_STEPS=""
FINAL_ONLY=0
RSCRIPT_BIN=""
OUTDIR=/mnt/d/analysis/maha
MAHA_JOBS=4
MAHA_NULL_MAX=0
MAHA_HELPER_DIR=/mnt/d/scripts/0f
MAHA_REPLACE_STEPS=""
MAHA_TILE_PROXY=http://127.0.0.1:7897
CHNS_PROJECT_DIR=/mnt/d/data/chns
CHNS_FOOD_MAP=""
CHNS_NUTRIENT_MAP=""
CHNS_BASELINE_WAVE=2011
CHNS_FOLLOWUP_END=2015
CHNS_MIN_AGE=40
CHNS_MIN_DIET_DAYS=2
CHNS_RUN_MORTALITY_SWEEP=1
CHNS_MORTALITY_SWEEP_WAVES=2004,2006,2009,2011
CHNS_MORTALITY_MODEL_AGE=40
CHNS_MORTALITY_EVENT_TARGET=500
while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps)
      [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: --steps requires a value." >&2; usage >&2; exit 2; }
      INTERNAL_STEPS="$2"; shift 2 ;;
    final)
      (( FINAL_ONLY == 0 )) || { echo "ERROR: final specified more than once." >&2; exit 2; }
      FINAL_ONLY=1; shift ;;
    --rscript) RSCRIPT_BIN=${2:?ERROR: --rscript requires FILE}; shift 2 ;;
    --output-dir) OUTDIR=${2:?ERROR: --output-dir requires DIR}; shift 2 ;;
    --jobs) MAHA_JOBS=${2:?ERROR: --jobs requires N}; shift 2 ;;
    --null-max) MAHA_NULL_MAX=${2:?ERROR: --null-max requires N}; shift 2 ;;
    --helper-dir) MAHA_HELPER_DIR=${2:?ERROR: --helper-dir requires DIR}; shift 2 ;;
    --replace-steps) MAHA_REPLACE_STEPS=${2:?ERROR: --replace-steps requires CSV}; shift 2 ;;
    --tile-proxy) MAHA_TILE_PROXY=${2:?ERROR: --tile-proxy requires URL}; shift 2 ;;
    --chns-project-dir) CHNS_PROJECT_DIR=${2:?ERROR: --chns-project-dir requires DIR}; shift 2 ;;
    --chns-food-map) CHNS_FOOD_MAP=${2:?ERROR: --chns-food-map requires FILE}; shift 2 ;;
    --chns-nutrient-map) CHNS_NUTRIENT_MAP=${2:?ERROR: --chns-nutrient-map requires FILE}; shift 2 ;;
    --chns-baseline-wave) CHNS_BASELINE_WAVE=${2:?ERROR: --chns-baseline-wave requires YEAR}; shift 2 ;;
    --chns-followup-end) CHNS_FOLLOWUP_END=${2:?ERROR: --chns-followup-end requires YEAR}; shift 2 ;;
    --chns-min-age) CHNS_MIN_AGE=${2:?ERROR: --chns-min-age requires N}; shift 2 ;;
    --chns-min-diet-days) CHNS_MIN_DIET_DAYS=${2:?ERROR: --chns-min-diet-days requires N}; shift 2 ;;
    --chns-run-sweep) CHNS_RUN_MORTALITY_SWEEP=${2:?ERROR: --chns-run-sweep requires 0 or 1}; shift 2 ;;
    --chns-sweep-waves) CHNS_MORTALITY_SWEEP_WAVES=${2:?ERROR: --chns-sweep-waves requires CSV}; shift 2 ;;
    --chns-model-age) CHNS_MORTALITY_MODEL_AGE=${2:?ERROR: --chns-model-age requires N}; shift 2 ;;
    --chns-event-target) CHNS_MORTALITY_EVENT_TARGET=${2:?ERROR: --chns-event-target requires N}; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    ukb|nhanes|chns|all)
      (( STEP_SET == 0 )) || { echo "ERROR: multiple modules specified: $STEP and $1" >&2; usage >&2; exit 2; }
      STEP="$1"; STEP_SET=1; shift ;;
    *)
      echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$STEP" in ukb|nhanes|chns|all) ;; *) echo "ERROR: unsupported step: $STEP" >&2; exit 2;; esac
[[ "$CHNS_RUN_MORTALITY_SWEEP" =~ ^[01]$ ]] || { echo "ERROR: --chns-run-sweep must be 0 or 1." >&2; exit 2; }
[[ "$CHNS_MORTALITY_MODEL_AGE" =~ ^[0-9]+$ ]] || { echo "ERROR: --chns-model-age must be an integer." >&2; exit 2; }
[[ "$CHNS_MORTALITY_EVENT_TARGET" =~ ^[0-9]+$ ]] || { echo "ERROR: --chns-event-target must be an integer." >&2; exit 2; }
[[ "$CHNS_MORTALITY_SWEEP_WAVES" =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "ERROR: --chns-sweep-waves must be comma-separated years." >&2; exit 2; }

RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: cannot find Rscript. Use --rscript FILE or add Rscript to PATH." >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CHNS_FOOD_MAP=${CHNS_FOOD_MAP:-$CHNS_PROJECT_DIR/MAHA_foodcode/chns_foodcode_map.tsv}
CHNS_NUTRIENT_MAP=${CHNS_NUTRIENT_MAP:-$CHNS_PROJECT_DIR/MAHA_foodcode/chns_foodcode_nutrients.tsv}
export MAHA_OUTDIR="$OUTDIR"
export RSCRIPT_BIN
export MAHA_JOBS MAHA_NULL_MAX MAHA_HELPER_DIR MAHA_REPLACE_STEPS MAHA_TILE_PROXY
export CHNS_PROJECT_DIR CHNS_FOOD_MAP CHNS_NUTRIENT_MAP CHNS_BASELINE_WAVE CHNS_FOLLOWUP_END CHNS_MIN_AGE CHNS_MIN_DIET_DAYS
export CHNS_RUN_MORTALITY_SWEEP CHNS_MORTALITY_SWEEP_WAVES
export CHNS_MORTALITY_MODEL_AGE CHNS_MORTALITY_EVENT_TARGET
mkdir -p "$OUTDIR"

configure_geospatial_env() {
  unset PROJ_LIB PROJ_DATA GDAL_DATA

  local proj_dir=""
  local gdal_dir=""
  local candidate
  for candidate in /usr/share/proj /usr/local/share/proj; do
    if [[ -f "${candidate}/proj.db" ]]; then
      proj_dir="$candidate"
      break
    fi
  done
  for candidate in /usr/share/gdal /usr/local/share/gdal; do
    if [[ -d "$candidate" ]]; then
      gdal_dir="$candidate"
      break
    fi
  done

  [[ -z "$proj_dir" ]] || export PROJ_DATA="$proj_dir"
  [[ -z "$gdal_dir" ]] || export GDAL_DATA="$gdal_dir"
}

configure_geospatial_env

configure_tile_proxy() {
  local tile_test_url="https://tile.openstreetmap.org/0/0/0.png"
  local tile_proxy="$MAHA_TILE_PROXY"

  if [[ -n "${HTTPS_PROXY:-${https_proxy:-}}" || -n "${HTTP_PROXY:-${http_proxy:-}}" ]]; then
    echo "Tile proxy             : inherited from environment"
    return
  fi

  if command -v curl >/dev/null 2>&1 && \
     curl --proxy "$tile_proxy" --silent --show-error --fail \
       --connect-timeout 5 --max-time 15 --output /dev/null "$tile_test_url"; then
    export HTTP_PROXY="$tile_proxy" HTTPS_PROXY="$tile_proxy"
    export http_proxy="$tile_proxy" https_proxy="$tile_proxy"
    echo "Tile proxy             : $tile_proxy (OSM connectivity verified)"
  else
    echo "WARNING: OSM street tiles are unreachable directly and the candidate proxy failed: $tile_proxy" >&2
    echo "         Start/enable Clash mixed-port 7897, or use --tile-proxy with a working HTTP proxy." >&2
  fi
}

configure_tile_proxy

run_one() {
  local cohort="$1"
  local script="${SCRIPT_DIR}/f/MAHA_${cohort}.R"
  [[ -f "$script" ]] || { echo "ERROR: analysis script does not exist: $script" >&2; exit 3; }

  export MAHA_AUXDIR="${OUTDIR}/${cohort}"
  mkdir -p "$MAHA_AUXDIR"

  local launcher_log="${MAHA_AUXDIR}/launcher.log"
  : > "$launcher_log"
  {
    echo "======================================================================"
    echo "MAHA pipeline: ${cohort}"
    echo "Rscript: ${RSCRIPT_BIN}"
    echo "Script : ${script}"
    echo "Publication output root: ${OUTDIR}"
    echo "Cohort auxiliary dir   : ${MAHA_AUXDIR}"
    echo "Launcher log           : ${launcher_log}"
    echo "PROJ data              : ${PROJ_DATA:-<library default>}"
    echo "GDAL data              : ${GDAL_DATA:-<library default>}"
    echo "======================================================================"

    if [[ "$cohort" == "chns" ]]; then
      export CHNS_REQUIRE_NUTRIENT_MAP=1
      export CHNS_REQUIRE_ENHANCED_MAP=1
      echo "CHNS project      : ${CHNS_PROJECT_DIR}"
      echo "Food map          : ${CHNS_FOOD_MAP}"
      echo "Nutrient map      : ${CHNS_NUTRIENT_MAP}"
      echo "Baseline wave     : ${CHNS_BASELINE_WAVE}"
      echo "Follow-up end     : ${CHNS_FOLLOWUP_END}"
      echo "Minimum age       : ${CHNS_MIN_AGE}"
      echo "Minimum diet days : ${CHNS_MIN_DIET_DAYS}"
      echo "Run high-event sensitivity : ${CHNS_RUN_MORTALITY_SWEEP}"
      echo "Repeated-exposure waves    : ${CHNS_MORTALITY_SWEEP_WAVES}"
      echo "Repeated-exposure age      : ${CHNS_MORTALITY_MODEL_AGE}"
      echo "Event-count flag  : ${CHNS_MORTALITY_EVENT_TARGET}"
    fi

    local args=()
    [[ -z "$INTERNAL_STEPS" ]] || args+=("--steps=${INTERNAL_STEPS}")
    "$RSCRIPT_BIN" "$script" "${args[@]}"
  } 2>&1 | tee -a "$launcher_log"
}

run_final() {
  local scope="$1"
  local script="${SCRIPT_DIR}/f/MAHA_publication.R"
  [[ -f "$script" ]] || { echo "ERROR: publication script does not exist: $script" >&2; exit 3; }
  echo "[MAHA FINAL] Integrating publication figures and tables for scope=${scope}"
  "$RSCRIPT_BIN" "$script" "--cohort=${scope}"
}

if [[ "$FINAL_ONLY" -eq 0 ]]; then
  case "$STEP" in
    ukb) run_one ukb ;;
    nhanes) run_one nhanes ;;
    chns) run_one chns ;;
    all) run_one ukb; run_one nhanes; run_one chns ;;
  esac
fi

if [[ "$FINAL_ONLY" -eq 1 || -z "$INTERNAL_STEPS" ]]; then
  run_final "$STEP"
else
  echo "[MAHA] Analysis step selection completed; final assembly skipped for --steps=${INTERNAL_STEPS}."
  echo "[MAHA] After all required cohort outputs exist, run: ./maha.sh final"
fi

echo "[MAHA] Done. Publication figures are in: ${OUTDIR}"

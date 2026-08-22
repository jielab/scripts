#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./maha.sh [ukb|nhanes|chns|all] [--steps <internal-step-list>] [--pub-only]
  ./maha.sh --help

Examples:
  ./maha.sh all
  ./maha.sh                    # same as: ./maha.sh all

  ./maha.sh ukb
  ./maha.sh nhanes
  ./maha.sh chns
  ./maha.sh all --pub-only

  ./maha.sh all --pub-only

  ./maha.sh ukb --steps data_qc,fig2,fig3

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

Common environment variables:
  RSCRIPT_BIN          Rscript executable (default: Rscript in PATH)
  MAHA_OUTDIR          publication output root (default: /mnt/d/analysis/maha)
  MAHA_JOBS            parallel workers for matched-null fits (default: 4 on Unix)
  MAHA_NULL_MAX        0 = all direction-null scores; >0 = deterministic subset
  MAHA_HELPER_DIR      shared R helper directory (default: /mnt/d/scripts/0f)
  MAHA_FORCE_STEPS     completed cached step(s) to rerun, e.g. fig5 or all

CHNS defaults (all may be overridden):
  CHNS_PROJECT_DIR       /mnt/d/data/chns
  CHNS_NUTR3             <project>/raw/diet/nutr3_00.sas7bdat
  CHNS_C12DIET           <project>/raw/diet/c12diet.sas7bdat
  CHNS_MASTER            <project>/CHNS2.0-TaoB/CHNS原始数据/Master_ID_201908/mast_pub_12.sas7bdat
  CHNS_ROSTER            <project>/CHNS2.0-TaoB/CHNS原始数据/Master_ID_201908/rst_12.sas7bdat
  CHNS_SURVEY            <project>/CHNS2.0-TaoB/CHNS原始数据/Master_ID_201908/surveys_pub_12.sas7bdat
  CHNS_PEXAM             <project>/raw/individual/pexampub12.sav
  CHNS_PACT              <project>/raw/individual/pact_12.sas7bdat
  CHNS_FOOD_MAP          <project>/MAHA_foodcode/chns_foodcode_map.tsv
  CHNS_NUTRIENT_MAP      <project>/MAHA_foodcode/chns_foodcode_nutrients.tsv
  CHNS_BASELINE_WAVE     2011
  CHNS_FOLLOWUP_END      2015
  CHNS_MIN_AGE           40
  CHNS_MIN_DIET_DAYS     2
USAGE
}

STEP="all"
STEP_SET=0
INTERNAL_STEPS=""
PUB_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps)
      [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: --steps requires a value." >&2; usage >&2; exit 2; }
      INTERNAL_STEPS="$2"; shift 2 ;;
    --pub-only)
      PUB_ONLY=1; shift ;;
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

RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: cannot find Rscript. Set RSCRIPT_BIN or add Rscript to PATH." >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUTDIR=${MAHA_OUTDIR:-/mnt/d/analysis/maha}
export MAHA_OUTDIR="$OUTDIR"
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
  local tile_proxy="${MAHA_TILE_PROXY:-http://127.0.0.1:7897}"

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
    echo "         Start/enable Clash mixed-port 7897, or set MAHA_TILE_PROXY to a working HTTP proxy." >&2
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
      export CHNS_PROJECT_DIR="${CHNS_PROJECT_DIR:-/mnt/d/data/chns}"
      export CHNS_FOOD_MAP="${CHNS_FOOD_MAP:-${CHNS_PROJECT_DIR}/MAHA_foodcode/chns_foodcode_map.tsv}"
      export CHNS_NUTRIENT_MAP="${CHNS_NUTRIENT_MAP:-${CHNS_PROJECT_DIR}/MAHA_foodcode/chns_foodcode_nutrients.tsv}"
      export CHNS_REQUIRE_NUTRIENT_MAP="${CHNS_REQUIRE_NUTRIENT_MAP:-1}"
      export CHNS_REQUIRE_ENHANCED_MAP="${CHNS_REQUIRE_ENHANCED_MAP:-1}"
      export CHNS_BASELINE_WAVE="${CHNS_BASELINE_WAVE:-2011}"
      export CHNS_FOLLOWUP_END="${CHNS_FOLLOWUP_END:-2015}"
      export CHNS_MIN_AGE="${CHNS_MIN_AGE:-40}"
      export CHNS_MIN_DIET_DAYS="${CHNS_MIN_DIET_DAYS:-2}"
      echo "CHNS project      : ${CHNS_PROJECT_DIR}"
      echo "Food map          : ${CHNS_FOOD_MAP}"
      echo "Nutrient map      : ${CHNS_NUTRIENT_MAP}"
      echo "Baseline wave     : ${CHNS_BASELINE_WAVE}"
      echo "Follow-up end     : ${CHNS_FOLLOWUP_END}"
      echo "Minimum age       : ${CHNS_MIN_AGE}"
      echo "Minimum diet days : ${CHNS_MIN_DIET_DAYS}"
    fi

    local args=()
    [[ -z "$INTERNAL_STEPS" ]] || args+=("--steps=${INTERNAL_STEPS}")
    "$RSCRIPT_BIN" "$script" "${args[@]}"
  } 2>&1 | tee -a "$launcher_log"
}

run_publication() {
  local scope="$1"
  local script="${SCRIPT_DIR}/f/MAHA_publication.R"
  [[ -f "$script" ]] || { echo "ERROR: publication script does not exist: $script" >&2; exit 3; }
  echo "[MAHA] Building publication figures for scope=${scope}"
  "$RSCRIPT_BIN" "$script" "--cohort=${scope}"
}

if [[ "$PUB_ONLY" -eq 0 ]]; then
  case "$STEP" in
    ukb) run_one ukb ;;
    nhanes) run_one nhanes ;;
    chns) run_one chns ;;
    all) run_one ukb; run_one nhanes; run_one chns ;;
  esac
fi

if [[ "$PUB_ONLY" -eq 1 || -z "$INTERNAL_STEPS" ]]; then
  run_publication "$STEP"
else
  echo "[MAHA] Analysis step selection completed; publication assembly skipped for --steps=${INTERNAL_STEPS}."
  echo "[MAHA] After all required cohort outputs exist, run: ./maha.sh all --pub-only"
fi

echo "[MAHA] Done. Publication figures are in: ${OUTDIR}"

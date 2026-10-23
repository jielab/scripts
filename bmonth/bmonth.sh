#!/usr/bin/env bash
# Parse canonical and legacy options locally; shared repo helpers are optional
# for console logging and mandatory for the existing actual-analysis memory cap.
set -euo pipefail
ORIGINAL_ARGS=("$@")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage(){ cat <<'HELP'
Modules:
    roster / rebuild / check / prep / phewas / assoc / decomp / followup / omics / med
    temp / tempmed / final / plot / tests
    Default: check -> prep -> phewas -> assoc -> decomp -> followup -> final -> plot
    Existing phenotype panels are reused; rebuild must be requested explicitly.
    --with-omics adds omics and med before final.
    --with-temperature adds temp -> tempmed -> final after core figures.
    Stages run in the foreground and stop on failure.

Main options:
    --dir-out DIR             Output [/mnt/d/analysis/bmonth]
    --dir-in DIR              UKB directory containing Rdata/ [/mnt/d/data/ukb]
    --phewas-input FILE       Explicit panel (otherwise discover output, then input panel)
    --outcomes CSV            Incident outcomes [cvd_cad,ra]
    --phewas-cores N|auto     Workers [4; auto: up to 8, memory-limited]
    --supplemental-covariates FILE  RDS with eid and new covariates
    --followup-extra-covars CSV      Explicit targeted sensitivity covariates
    --omics-covariates CSV    Explicit omics covariates
    --med-mode screen|bootstrap [screen]
    --med-candidates FILE     Locked CSV required for bootstrap
    --med-layers CSV          prot and/or met [prot,met]
    --med-boot N              Bootstrap replicates [100]
    --dir-temp DIR            HadUK-Grid directory [/mnt/d/data/ukb/temperature]
    --rscript FILE            Rscript executable
    --replace                Recompute compatible model checkpoints
    -h, --help               Help; all options and prerequisites: 0.README.md

Example:
    ./bmonth.sh check
    ./bmonth.sh followup --followup-extra-covars ethnic.c
    ./bmonth.sh med --med-mode screen --med-layers prot
    ./bmonth.sh med --med-mode bootstrap --med-candidates candidates.csv --med-boot 2000
    ./bmonth.sh  # Reuse an existing panel and run the core workflow
HELP
}

STEP=workflow
OUTDIR="${BMONTH_OUTDIR:-/mnt/d/analysis/bmonth}"
INDIR="${BMONTH_INDIR:-/mnt/d/data/ukb}"
RSCRIPT_BIN="${BMONTH_RSCRIPT:-}"
OUTCOMES="${BMONTH_OUTCOMES:-cvd_cad,ra}"
EXPOSURES="${BMONTH_EXPOSURES-bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn}"
CONTEXT_EXPOSURES="${BMONTH_CONTEXT_EXPOSURES-birth_north_100km,birth_east_100km}"
EXPLORATORY_TRAITS="${BMONTH_EXPLORATORY_TRAITS:-border12,border1n,border_rank,border_relative}"
PHEWAS_CORES="${BMONTH_PHEWAS_CORES:-4}"
PHEWAS_INPUT="${BMONTH_PHEWAS_INPUT:-}"
ALLOW_LEGACY="${BMONTH_ALLOW_LEGACY_PHEWAS:-FALSE}"
SUPPLEMENTAL="${BMONTH_SUPPLEMENTAL_COVARIATES:-}"
EXTRA_COVARS="${BMONTH_EXTRA_COVARS:-}"
FOLLOWUP_EXTRAS="${BMONTH_FOLLOWUP_EXTRA_COVARS-}"
FOLLOWUP_EXTRAS_SET=FALSE
if [[ ${BMONTH_FOLLOWUP_EXTRA_COVARS+x} ]]; then FOLLOWUP_EXTRAS_SET=TRUE; fi
OMICS_COVARS="${BMONTH_OMICS_COVARIATES:-}"
WITH_OMICS="${BMONTH_WITH_OMICS:-FALSE}"
WITH_TEMPERATURE="${BMONTH_WITH_TEMPERATURE:-FALSE}"
MED_MODE="${BMONTH_MED_MODE:-screen}"
MED_CANDIDATES="${BMONTH_MED_CANDIDATES:-}"
MED_LAYERS="${BMONTH_MED_LAYERS:-prot,met}"
MED_CENTER="${BMONTH_MED_CENTER:-strata}"
MED_BOOT="${BMONTH_MED_BOOT:-100}"
FOLLOW_END="${BMONTH_FOLLOW_END:-}"
TEMP_DIR="${BMONTH_TEMP_DIR:-/mnt/d/data/ukb/temperature}"
TEMP_RES="${BMONTH_TEMP_RES:-1km}"
TEMPMED_PHENOS="${BMONTH_TEMPMED_PHENOS:-}"
REPLACE="${BMONTH_REPLACE:-FALSE}"

need_value(){ if [[ $# -lt 2 ]]; then echo "ERROR: $1 requires a value" >&2; exit 2; fi; }
truthy(){ case "${1^^}" in TRUE|1|YES) return 0;; *) return 1;; esac; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    roster|rebuild|check|prep|phewas|assoc|decomp|followup|omics|med|temp|tempmed|final|plot|tests) STEP="$1"; shift ;;
    all) STEP=workflow; shift ;; # Compatibility with older commands; no module is needed.
    --dir-out|--out) need_value "$@"; OUTDIR="$2"; shift 2 ;;
    --dir-in|--in) need_value "$@"; INDIR="$2"; shift 2 ;;
    --dir-temp|--temp-dir) need_value "$@"; TEMP_DIR="$2"; shift 2 ;;
    --rscript) need_value "$@"; RSCRIPT_BIN="$2"; shift 2 ;;
    --phewas-input) need_value "$@"; PHEWAS_INPUT="$2"; shift 2 ;;
    --allow-legacy-phewas) ALLOW_LEGACY=TRUE; shift ;;
    --supplemental-covariates) need_value "$@"; SUPPLEMENTAL="$2"; shift 2 ;;
    --extra-covars) need_value "$@"; EXTRA_COVARS="$2"; shift 2 ;;
    --followup-extra-covars) need_value "$@"; FOLLOWUP_EXTRAS="$2"; FOLLOWUP_EXTRAS_SET=TRUE; shift 2 ;;
    --omics-covariates) need_value "$@"; OMICS_COVARS="$2"; shift 2 ;;
    --outcomes) need_value "$@"; OUTCOMES="$2"; shift 2 ;;
    --exposures) need_value "$@"; EXPOSURES="$2"; shift 2 ;;
    --context-exposures) need_value "$@"; CONTEXT_EXPOSURES="$2"; shift 2 ;;
    --exploratory-traits) need_value "$@"; EXPLORATORY_TRAITS="$2"; shift 2 ;;
    --phewas-cores) need_value "$@"; PHEWAS_CORES="$2"; shift 2 ;;
    --with-omics) WITH_OMICS=TRUE; shift ;;
    --with-temperature) WITH_TEMPERATURE=TRUE; shift ;;
    --med-mode) need_value "$@"; MED_MODE="$2"; shift 2 ;;
    --med-candidates) need_value "$@"; MED_CANDIDATES="$2"; shift 2 ;;
    --med-layers) need_value "$@"; MED_LAYERS="$2"; shift 2 ;;
    --med-center) need_value "$@"; MED_CENTER="$2"; shift 2 ;;
    --med-boot) need_value "$@"; MED_BOOT="$2"; shift 2 ;;
    --follow-end) need_value "$@"; FOLLOW_END="$2"; shift 2 ;;
    --temp-res) need_value "$@"; TEMP_RES="$2"; shift 2 ;;
    --tempmed-phenos) need_value "$@"; TEMPMED_PHENOS="$2"; shift 2 ;;
    --replace) REPLACE=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2 ;;
  esac
done
if ! truthy "$FOLLOWUP_EXTRAS_SET"; then FOLLOWUP_EXTRAS="$EXTRA_COVARS"; fi
PHEWAS_INPUT_EXPLICIT=FALSE
if [[ -n "$PHEWAS_INPUT" ]]; then PHEWAS_INPUT_EXPLICIT=TRUE; fi
# Rebuild writes to its own target, independent of ordinary read discovery.
REBUILD_TARGET="${PHEWAS_INPUT:-$OUTDIR/phewas/PheWAS.rds}"
if ! truthy "$PHEWAS_INPUT_EXPLICIT"; then
  if [[ -f "$OUTDIR/phewas/PheWAS.rds" ]]; then
    PHEWAS_INPUT="$OUTDIR/phewas/PheWAS.rds"
  elif [[ -f "$INDIR/Rdata/PheWAS.rds" ]]; then
    PHEWAS_INPUT="$INDIR/Rdata/PheWAS.rds"
  else
    PHEWAS_INPUT="$OUTDIR/phewas/PheWAS.rds"
  fi
fi
# Keep no-module invocation intact across console/memory-guard re-entry.
# A real option also enables the shared console logger, which bypasses empty args.
if [[ ${#ORIGINAL_ARGS[@]} -eq 0 ]]; then ORIGINAL_ARGS=(--dir-out "$OUTDIR"); fi
[[ "$PHEWAS_CORES" == auto || "$PHEWAS_CORES" =~ ^[1-9][0-9]*$ ]] || {
  echo 'ERROR: --phewas-cores must be auto or a positive integer' >&2; exit 2; }
[[ "$MED_BOOT" =~ ^(0|[1-9][0-9]*)$ ]] || {
  echo 'ERROR: --med-boot must be a nonnegative integer' >&2; exit 2; }
[[ "$MED_MODE" == screen || "$MED_MODE" == bootstrap ]] || {
  echo 'ERROR: --med-mode must be screen or bootstrap' >&2; exit 2; }
[[ "$MED_CENTER" == strata || "$MED_CENTER" == coefficient ]] || {
  echo 'ERROR: --med-center must be strata or coefficient' >&2; exit 2; }
[[ -z "$FOLLOW_END" || "$FOLLOW_END" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || {
  echo 'ERROR: --follow-end must have YYYY-MM-DD format' >&2; exit 2; }
[[ -n "$OUTDIR" && -n "$INDIR" ]] || { echo 'ERROR: input/output directories cannot be empty' >&2; exit 2; }
RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] && command -v "$RSCRIPT_BIN" >/dev/null 2>&1 || {
  echo 'ERROR: Rscript executable not found; use --rscript FILE' >&2; exit 1; }

export BMONTH_OUTDIR="$OUTDIR" BMONTH_INDIR="$INDIR" BMONTH_OUTCOMES="$OUTCOMES"
export BMONTH_EXPOSURES="$EXPOSURES" BMONTH_CONTEXT_EXPOSURES="$CONTEXT_EXPOSURES"
export BMONTH_EXPLORATORY_TRAITS="$EXPLORATORY_TRAITS" BMONTH_PHEWAS_CORES="$PHEWAS_CORES"
# Keep the explicit configuration separate across console/memory-guard re-entry.
if truthy "$PHEWAS_INPUT_EXPLICIT"; then export BMONTH_PHEWAS_INPUT="$PHEWAS_INPUT"; else unset BMONTH_PHEWAS_INPUT; fi
export BMONTH_PHEWAS_SELECTED_INPUT="$PHEWAS_INPUT"
export BMONTH_ALLOW_LEGACY_PHEWAS="$ALLOW_LEGACY" BMONTH_SUPPLEMENTAL_COVARIATES="$SUPPLEMENTAL"

export BMONTH_EXTRA_COVARS="$EXTRA_COVARS" BMONTH_FOLLOWUP_EXTRA_COVARS="$FOLLOWUP_EXTRAS"
export BMONTH_OMICS_COVARIATES="$OMICS_COVARS" BMONTH_WITH_OMICS="$WITH_OMICS" BMONTH_WITH_TEMPERATURE="$WITH_TEMPERATURE"
export BMONTH_MED_MODE="$MED_MODE" BMONTH_MED_CANDIDATES="$MED_CANDIDATES" BMONTH_MED_LAYERS="$MED_LAYERS"
export BMONTH_MED_CENTER="$MED_CENTER" BMONTH_MED_BOOT="$MED_BOOT" BMONTH_FOLLOW_END="$FOLLOW_END"
export BMONTH_TEMP_DIR="$TEMP_DIR" BMONTH_TEMP_RES="$TEMP_RES" BMONTH_TEMPMED_PHENOS="$TEMPMED_PHENOS"
export BMONTH_REPLACE="$REPLACE" BMONTH_FDIR="$SCRIPT_DIR/f" BMONTH_TEST_RSCRIPT="$RSCRIPT_BIN"

if [[ "$STEP" == tests ]]; then
  # Restore a delivered external test archive; no test sources in production.
  TEST_DIR="${BMONTH_TEST_DIR:-}"
  if [[ -z "$TEST_DIR" ]]; then
    TEST_ARCHIVE="${BMONTH_TEST_ARCHIVE:-/mnt/d/bmonth-delivery/bmonth-tests-20261020.tar.gz}"
    [[ -f "$TEST_ARCHIVE" ]] || {
      echo "ERROR: missing test archive $TEST_ARCHIVE; supply BMONTH_TEST_ARCHIVE or BMONTH_TEST_DIR (see 0.README.md)." >&2
      exit 1
    }
    TEST_DIR=$(mktemp -d /tmp/bmonth-tests-restored-XXXXXX)
    trap 'rm -rf -- "$TEST_DIR"' EXIT
    tar -xzf "$TEST_ARCHIVE" -C "$TEST_DIR"
  fi
  [[ -f "$TEST_DIR/run_tests.R" ]] || { echo "ERROR: missing $TEST_DIR/run_tests.R" >&2; exit 1; }
  "$RSCRIPT_BIN" --vanilla "$TEST_DIR/run_tests.R" "$SCRIPT_DIR/f"
  exit $?
fi

export SCRIPT_LOG_DIR="${SCRIPT_LOG_DIR:-$OUTDIR/logs}"
mkdir -p "$OUTDIR" "$SCRIPT_LOG_DIR"
if [[ -f "$SCRIPT_DIR/../0f/console.sh" ]]; then
  source "$SCRIPT_DIR/../0f/console.sh" "${ORIGINAL_ARGS[@]}" 2>/dev/null || true
fi
[[ -f "$SCRIPT_DIR/../0f/memory_cap.sh" ]] || {
  echo 'ERROR: missing ../0f/memory_cap.sh; actual analysis must run inside the scripts repository with its existing memory guard.' >&2
  exit 1
}
source "$SCRIPT_DIR/../0f/memory_cap.sh"
memory_cap_enter bmonth 24 0 "$SCRIPT_DIR/bmonth.sh" "${ORIGINAL_ARGS[@]}"

RUN_ID="$(date -u '+%Y%m%dT%H%M%SZ')-$$"
STATUS_FILE="$SCRIPT_LOG_DIR/stage_status.tsv"
if [[ ! -s "$STATUS_FILE" ]]; then printf 'run_id\ttime_utc\tstage\tstatus\texit_code\n' > "$STATUS_FILE"; fi
ACTIVE_STAGE=""
record_stage(){ printf '%s\t%s\t%s\t%s\t%s\n' "$RUN_ID" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$1" "$2" "$3" >> "$STATUS_FILE"; }
handle_exit(){
  local rc="$1"
  if [[ "$rc" -ne 0 && -n "$ACTIVE_STAGE" ]]; then record_stage "$ACTIVE_STAGE" failed "$rc"; fi
}
trap 'handle_exit $?' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
run_r(){
  local script="$1" label="${2:-${1##*/}}" rc=0
  shift 2
  ACTIVE_STAGE="$label"; record_stage "$label" started ""
  echo "[bmonth] START $label"
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 "$RSCRIPT_BIN" --vanilla "$script" "$@" || rc=$?
  if [[ "$rc" -ne 0 ]]; then
    record_stage "$label" failed "$rc"; ACTIVE_STAGE=""
    echo "[bmonth] FAILED $label (exit $rc)" >&2; return "$rc"
  fi
  record_stage "$label" complete 0; ACTIVE_STAGE=""
  echo "[bmonth] DONE $label"
}


# 🚩 Optional rebuild; ordinary analysis only discovers existing inputs.
prepare_roster(){
  local panel_dir
  panel_dir=$(dirname "$REBUILD_TARGET")
  if [[ ! -f "$panel_dir/phewas_observable_roster.rds" ]]; then
    run_r "$SCRIPT_DIR/../ukb/f/build_phewas_roster.R" roster \
      --dir-in "$INDIR" --dir-out "$panel_dir"
  fi
}
ensure_phewas_panel(){
  [[ -f "$PHEWAS_INPUT" ]] && return 0
  if truthy "$PHEWAS_INPUT_EXPLICIT"; then
    echo "ERROR: specified phenotype panel does not exist: $PHEWAS_INPUT" >&2
  else
    echo "ERROR: no existing phenotype panel. Missing: $OUTDIR/phewas/PheWAS.rds; $INDIR/Rdata/PheWAS.rds" >&2
    printf 'Optional rebuild: %q rebuild --dir-in %q --dir-out %q\n' "$SCRIPT_DIR/bmonth.sh" "$INDIR" "$OUTDIR" >&2
  fi
  return 1
}


case "$STEP" in
  roster|rebuild)
    [[ "${REBUILD_TARGET##*/}" == PheWAS.rds ]] || { echo 'ERROR: rebuild requires a PheWAS.rds basename' >&2; exit 2; }
    panel_dir=$(dirname "$REBUILD_TARGET")
    replace_args=()
    if truthy "$REPLACE"; then replace_args=(--replace); fi
    if truthy "$REPLACE"; then
      run_r "$SCRIPT_DIR/../ukb/f/build_phewas_roster.R" roster \
        --dir-in "$INDIR" --dir-out "$panel_dir" --replace
    else
      prepare_roster
    fi
    if [[ "$STEP" == rebuild ]]; then
      run_r "$SCRIPT_DIR/../ukb/f/rebuild_phewas.R" rebuild \
        --code "$INDIR/Rdata/icd10.code0.rds" --date "$INDIR/Rdata/icd10.date0.rds" \
        --roster "$panel_dir/phewas_observable_roster.rds" --dir-out "$panel_dir" "${replace_args[@]}"
    fi
    ;;
  check) ensure_phewas_panel; run_r "$SCRIPT_DIR/f/00.check.R" check ;;
  prep) run_r "$SCRIPT_DIR/f/1.prep.R" prep ;;
  phewas) ensure_phewas_panel; run_r "$SCRIPT_DIR/f/3.phewas.R" phewas ;;
  assoc) run_r "$SCRIPT_DIR/f/2.assoc.R" assoc ;;
  decomp) ensure_phewas_panel; run_r "$SCRIPT_DIR/f/5.decomp.R" decomp ;;
  followup)
    if [[ ",${BMONTH_FOLLOWUP_COMPONENTS:-targeted,cycle}," == *,targeted,* ]]; then
      ensure_phewas_panel
    fi
    run_r "$SCRIPT_DIR/f/8.followup.R" followup
    ;;
  omics) run_r "$SCRIPT_DIR/f/9.omics_followup.R" omics ;;
  med) run_r "$SCRIPT_DIR/f/4.med.R" med ;;
  temp) run_r "$SCRIPT_DIR/f/6.temperature.R" temp ;;
  tempmed) run_r "$SCRIPT_DIR/f/7.tempmed.R" tempmed ;;
  final) run_r "$SCRIPT_DIR/f/5.final.R" final ;;
  plot) run_r "$SCRIPT_DIR/f/6.replot.R" plot ;;
  workflow)
    export BMONTH_CHECK_SURVIVAL=TRUE
    ensure_phewas_panel
    run_r "$SCRIPT_DIR/f/00.check.R" check
    run_r "$SCRIPT_DIR/f/1.prep.R" prep
    run_r "$SCRIPT_DIR/f/3.phewas.R" phewas
    run_r "$SCRIPT_DIR/f/2.assoc.R" assoc
    run_r "$SCRIPT_DIR/f/5.decomp.R" decomp
    run_r "$SCRIPT_DIR/f/8.followup.R" followup
    if truthy "$WITH_OMICS"; then
      run_r "$SCRIPT_DIR/f/9.omics_followup.R" omics
      run_r "$SCRIPT_DIR/f/4.med.R" med
    fi
    run_r "$SCRIPT_DIR/f/5.final.R" final_core
    run_r "$SCRIPT_DIR/f/6.replot.R" plot_core
    if truthy "$WITH_TEMPERATURE"; then
      run_r "$SCRIPT_DIR/f/6.temperature.R" temp
      run_r "$SCRIPT_DIR/f/7.tempmed.R" tempmed
      run_r "$SCRIPT_DIR/f/5.final.R" final_temperature
    fi
    ;;
esac

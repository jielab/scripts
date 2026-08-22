#!/usr/bin/env bash
set -euo pipefail



usage() {
  cat <<'HELP'
Usage:
  ./ckm.sh [options]

Core options:
      The analysis target is fixed to incident clinical CKM stage 4:
      first CHD/HF/stroke/PAD/AF among baseline stages 0-3. Baseline stage-4
      participants are retained as Yang evidence. There is no --trait option.
      --biom X               prot, met, prot,met, or prot+met. Default: prot
      --start-step STEP      build_ckm, summarize_ckm, stage_clock, data_prep, marker_map,
                             biom_refine, yy_timeline, stage_specific, causal_triage,
                             prediction, consolidate. Default: stage_clock
      --stop-step STEP       Stop after this step. Default: consolidate.
                             Useful for: --start-step stage_clock --stop-step stage_clock
      --scan-all-biom        Scan all available biomarkers.
      --max-biom N           Debug limit after biomarker QC. 0 = no limit.
      --biomarkers CSV       Biomarkers forced into display figures.
      --ethnic X             Analysis filter, e.g. all or White. Default: all

CKM construction:
      --rebuild-ckm          Rebuild ckm.rds even when it already exists.
      --summarize_ckm        Regenerate CKM summary figures and Excel files
                             from the existing dual-stage ckm.rds.
      --stage4-source X      all or icd10. Default: all. "all" uses all dated
                             t2e-compatible sources; exact ICD-10 is sensitivity.
      --stage4-require-risk X
                             TRUE/FALSE. The 2026 staging definition requires
                             clinical CVD with CKM risk factors. Default: TRUE;
                             FALSE is a sensitivity analysis matching a broader
                             clinical-CVD-only implementation.
      --stage3-method X     SCORE2 or PREVENT. Selects which stored CKM stage
                             is analyzed. Default: SCORE2.
      --prevent-cutoff X     PREVENT-CVD 10-year cutoff for stage 3. Default: 0.20
      --prevent-model X      base or best_available. Default: base. The base
                             equation is used consistently across participants;
                             best_available is a sensitivity analysis.
      --tg-fasting-cutoff X  Fasting TG cutoff in mmol/L. Default: 1.69.
      --tg-nonfasting-cutoff X
                             Nonfasting/unknown-fasting TG cutoff. Default: 1.98.
      --statin-col COL       Statin-use column for PREVENT; auto-detects drug.statin,
                             statin, statin_use, then drug.lipid.
      --allow-unknown-fasting-glucose
                             Treat glucose with unknown fasting time as usable.
                             Default: off; HbA1c remains usable.
      --allow-unknown-fasting-tg
                             This flag enables the pragmatic nonfasting/unknown sensitivity.
      --uacr-col COL         Existing UACR column in mg/g.
      --urine-albumin-col C  Urine albumin column in mg/L.
      --urine-creatinine-col C
                             Urine creatinine column in umol/L; used to derive UACR.
      --subclinical-ascvd-col C
                             Optional baseline 0/1 subclinical ASCVD flag.
      --pre-hf-col C         Optional baseline 0/1 pre-heart-failure flag.
      --cac-col C            Optional coronary calcium score; CAC >=100 is stage 3.
      --abi-col C            Optional ankle-brachial index; ABI <0.90 is stage 3.
      --merge-ckm-to-all     Also merge compact CKM stage/date/subtype fields into all.rds.
      --audit-dir DIR        CKM construction audit directory.
                             Default: /mnt/d/analysis/ckm/audit. The directory
                             is created only when ckm.rds is built/rebuilt;
                             downstream-only runs do not create an empty audit folder.
      --icd10-list FILE      ICD-10 phenotype definition list used for audit.
                             Default: PHEDIR/common/icd10.lst
      --sbp-col COL          Baseline systolic BP column. Auto: sbp, then sbp_auto.
      --dbp-col COL          Baseline diastolic BP column. Auto: dbp, then dbp_auto.

Stage-clock options:
      --clock-reference X    healthy_le8 or population. Default: healthy_le8.
                             healthy_le8 standardizes the four behavioral LE8
                             scores (diet, activity, smoking, sleep), not the
                             physiological scores used to define CKM stages.
                             RMST, risk-age, and risk-calibrated sensitivity models use the full UKB cohort.
      --clock-adjust-mode X  crude, demog, le8_behavior, or le8_full.
                             Default: le8_behavior. le8_full is an overadjusted
                             sensitivity model because BMI/BP/HbA1c/non-HDL
                             contribute to CKM stage definitions.
      --clock-draws N        Coefficient simulations for clock uncertainty. Default: 500
      --clock-ref-n N        Reference rows for point standardization. Default: 5000
      --clock-draw-ref-n N   Reference rows used inside coefficient simulations.
                             Default: 1000; lowers RAM without changing point estimates.
      --clock-tau X          Requested RMST horizon in years. Default: 15.
                             The effective horizon is capped at supported follow-up.
      --risk-age-reference X Reference age used to display equivalent age. Default: 60.
Adjustment and prediction:
      --pred-base-mode X     demog, behavior, le8, or full. Default: behavior
      --center-mode X        none or factor. Default: none
      --pred-folds N         Outer folds. Default: 5
      --pred-screen-n N      Biomarkers retained by fold-specific screening. Default: 300
      --pred-top-n N         Biomarkers in explainable YY score. Default: 50
      --pred-min-yin-n N     Minimum eligible Yin participants. Default: 1000
      --pred-min-events N    Minimum held-out target events. Default: 50
      --lag-years CSV        Prospective lag analyses. Default: 0,2,5,10
      --lag-min-retention X  Minimum absolute lag-5/main beta ratio for a
                             lag-persistent association. Default: 0.50.
      --lag-min-hr X         Minimum HR magnitude for both main and lag-5 effects.
                             Default: 1.05 (protective equivalent <=1/1.05).
      --no-prediction        Skip cross-validated prediction.
      --benchmark-glmnet     Add an optional elastic-net benchmark using only
                             fold-specific Yin-selected biomarkers.
      --met-strict-exclude-lipid-block TRUE/FALSE
                             For biom=met, also run a strict sensitivity analysis
                             excluding the lipoprotein-lipid block. Default: TRUE.
      --seed N               Random seed. Default: 2026
      --prevent-batch-size N PREVENT rows per progress batch. Default: 5000
      --feature-chunk-size N Biomarker scan features per progress chunk.
                             Default: automatic.
      --ncore N              Parallel feature-scan workers. Default: min(2, cores-2).
                             Use 1 when RAM is limited.
  -h, --help                 Show this message.

Biomarker semantics:
  prot,met  runs proteins and metabolites independently.
  prot+met  merges both layers in the common participant cohort.

Output layout:
  Shared CKM outputs such as Fig1.stage_clock.png and, after a prot,met run,
  Fig7.cross_omic_comparison.png are written directly under CKM_OUTDIR
  (default: /mnt/d/analysis/ckm), with raw shared tables under
  CKM_OUTDIR/raw. Biomarker-specific outputs are written under
  CKM_OUTDIR/prot, CKM_OUTDIR/met, or CKM_OUTDIR/prot+met. There is no
  [biom]/[Y] subdirectory layer.

Examples:
  # Step 1 (run once): calculate and store both SCORE2 and PREVENT CKM stages.
  ./ckm.sh --rebuild-ckm

  # Regenerate the CKM summary figures and Excel workbooks.
  ./ckm.sh --summarize_ckm

  # Step 2: run protein and metabolite analyses using the existing ckm.rds and Fig1.
  # SCORE2 analysis (default):
  ./ckm.sh --biom prot,met --scan-all-biom

  # PREVENT analysis using the same ckm.rds.
  ./ckm.sh --stage3-method PREVENT --biom prot,met --scan-all-biom

  # Full analysis with commonly used optional settings. This also reuses ckm.rds.
  ./ckm.sh --biom prot,met --scan-all-biom --clock-reference healthy_le8 --clock-adjust-mode le8_behavior --clock-tau 15 --risk-age-reference 60 --pred-base-mode behavior --center-mode none --lag-years 0,2,5,10 --lag-min-retention 0.50 --lag-min-hr 1.05 --met-strict-exclude-lipid-block TRUE --benchmark-glmnet

  # Full protein-only analysis.
  ./ckm.sh --biom prot --scan-all-biom --clock-adjust-mode le8_behavior --clock-tau 15 --pred-base-mode behavior --center-mode none

  # Recreate the comparative Fig1 and CKM summary workbooks.
  ./ckm.sh --summarize_ckm

  # Resume from stage-specific analysis after existing marker/stage outputs.
  ./ckm.sh --biom prot --start-step stage_specific --scan-all-biom --pred-base-mode behavior --center-mode none

  # Resume from causal triage after a valid stage_specific.rds exists.
  ./ckm.sh --biom prot --start-step causal_triage --scan-all-biom --pred-base-mode behavior --center-mode none

  # Merge protein and metabolite layers in their common participant cohort.
  ./ckm.sh --biom prot+met --scan-all-biom --clock-adjust-mode le8_behavior --clock-tau 15 --pred-base-mode behavior --center-mode none

  # Conservative-memory run.
  ./ckm.sh --biom prot --scan-all-biom --ncore 1 --feature-chunk-size 50 --pred-base-mode behavior --center-mode none
HELP
}

if [[ -n "${START_STEP+x}" ]]; then START_STEP_EXPLICIT=TRUE; else START_STEP_EXPLICIT=FALSE; fi
if [[ -n "${END_STEP+x}" ]]; then END_STEP_EXPLICIT=TRUE; else END_STEP_EXPLICIT=FALSE; fi
START_STEP="${START_STEP:-stage_clock}"
END_STEP="${END_STEP:-consolidate}"
BIOM_REQUEST="${CKM_BIOM:-prot}"
CKM_SCAN_ALL_BIOM="${CKM_SCAN_ALL_BIOM:-FALSE}"
CKM_REBUILD="${CKM_REBUILD:-FALSE}"
CKM_STAGE4_REQUIRE_RISK="${CKM_STAGE4_REQUIRE_RISK:-TRUE}"
CKM_MERGE_ALL="${CKM_MERGE_ALL:-FALSE}"
CKM_ETHNIC="${CKM_ETHNIC:-all}"
CKM_CLOCK_REFERENCE="${CKM_CLOCK_REFERENCE:-healthy_le8}"
CKM_CLOCK_ADJUST_MODE="${CKM_CLOCK_ADJUST_MODE:-le8_behavior}"
CKM_SKIP_PREDICTION="${CKM_SKIP_PREDICTION:-FALSE}"
CKM_BENCHMARK_GLMNET="${CKM_BENCHMARK_GLMNET:-FALSE}"
CKM_PRED_BASE_MODE="${CKM_PRED_BASE_MODE:-behavior}"
SEED="${SEED:-2026}"
CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK="${CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK:-TRUE}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --biom) BIOM_REQUEST="$2"; shift 2 ;;
    --start-step|--start_step|-s) START_STEP="$2"; START_STEP_EXPLICIT=TRUE; shift 2 ;;
    --stop-step|--end-step|--stop_step|--end_step) END_STEP="$2"; END_STEP_EXPLICIT=TRUE; shift 2 ;;
    --scan-all-biom|--scan-all-biomarkers|--scan-all-proteins) CKM_SCAN_ALL_BIOM=TRUE; shift ;;
    --max-biom|--max-biomarkers|--max-proteins) export CKM_MAX_BIOM="$2"; shift 2 ;;
    --biomarkers|--proteins) export CKM_BIOMARKERS="$2"; shift 2 ;;
    --ethnic) CKM_ETHNIC="$2"; shift 2 ;;
    --rebuild-ckm) CKM_REBUILD=TRUE; shift ;;
    --summarize_ckm|--summarize-ckm) START_STEP=summarize_ckm; END_STEP=summarize_ckm; START_STEP_EXPLICIT=TRUE; END_STEP_EXPLICIT=TRUE; shift ;;
    --stage4-source) export CKM_STAGE4_SOURCE="$2"; shift 2 ;;
    --stage4-require-risk) CKM_STAGE4_REQUIRE_RISK="$2"; shift 2 ;;
    --stage3-method) export CKM_STAGE3_METHOD="$2"; shift 2 ;;
    --prevent-cutoff) export CKM_PREVENT_CUTOFF="$2"; shift 2 ;;
    --prevent-model) export CKM_PREVENT_MODEL="$2"; shift 2 ;;
    --tg-fasting-cutoff) export CKM_TG_FASTING_CUTOFF="$2"; shift 2 ;;
    --tg-nonfasting-cutoff) export CKM_TG_NONFASTING_CUTOFF="$2"; shift 2 ;;
    --statin-col) export CKM_STATIN_COL="$2"; shift 2 ;;
    --allow-unknown-fasting-glucose) export CKM_ALLOW_UNKNOWN_FASTING_GLUCOSE=TRUE; shift ;;
    --allow-unknown-fasting-tg) export CKM_ALLOW_UNKNOWN_FASTING_TG=TRUE; shift ;;
    --uacr-col) export CKM_UACR_COL="$2"; shift 2 ;;
    --urine-albumin-col) export CKM_URINE_ALBUMIN_COL="$2"; shift 2 ;;
    --urine-creatinine-col) export CKM_URINE_CREATININE_COL="$2"; shift 2 ;;
    --subclinical-ascvd-col) export CKM_SUBCLINICAL_ASCVD_COL="$2"; shift 2 ;;
    --pre-hf-col) export CKM_PREHF_COL="$2"; shift 2 ;;
    --cac-col) export CKM_CAC_COL="$2"; shift 2 ;;
    --abi-col) export CKM_ABI_COL="$2"; shift 2 ;;
    --merge-ckm-to-all) CKM_MERGE_ALL=TRUE; shift ;;
    --audit-dir) export CKM_AUDIT_DIR="$2"; shift 2 ;;
    --icd10-list) export CKM_ICD10_LIST="$2"; shift 2 ;;
    --sbp-col) export CKM_SBP_COL="$2"; shift 2 ;;
    --dbp-col) export CKM_DBP_COL="$2"; shift 2 ;;
    --clock-reference) CKM_CLOCK_REFERENCE="$2"; shift 2 ;;
    --clock-adjust-mode) CKM_CLOCK_ADJUST_MODE="$2"; shift 2 ;;
    --clock-draws) export CKM_CLOCK_DRAWS="$2"; shift 2 ;;
    --clock-ref-n) export CKM_CLOCK_REF_N="$2"; shift 2 ;;
    --clock-draw-ref-n) export CKM_CLOCK_DRAW_REF_N="$2"; shift 2 ;;
    --clock-tau) export CKM_CLOCK_TAU="$2"; shift 2 ;;
    --risk-age-reference) export CKM_RISK_AGE_REFERENCE="$2"; shift 2 ;;
    --pred-base-mode|--pred_base_mode) CKM_PRED_BASE_MODE="$2"; shift 2 ;;
    --pred-folds|--pred_folds) export CKM_PRED_FOLDS="$2"; shift 2 ;;
    --pred-screen-n|--pred_screen_n) export CKM_PRED_SCREEN_N="$2"; shift 2 ;;
    --pred-top-n|--pred_top_n) export CKM_PRED_TOP_N="$2"; shift 2 ;;
    --pred-min-yin-n|--pred_min_yin_n) export CKM_PRED_MIN_YIN_N="$2"; shift 2 ;;
    --pred-min-events|--pred_min_events) export CKM_PRED_MIN_EVENTS="$2"; shift 2 ;;
    --lag-years|--lag_years) export CKM_LAG_YEARS="$2"; shift 2 ;;
    --lag-min-retention) export CKM_LAG_MIN_RETENTION="$2"; shift 2 ;;
    --lag-min-hr) export CKM_MARKER_MIN_COX_LOGHR="$(python3 - <<PY2
import math
x=float("$2")
if x <= 1:
    raise SystemExit("--lag-min-hr must be >1")
print(math.log(x))
PY2
)"; shift 2 ;;
    --no-prediction) CKM_SKIP_PREDICTION=TRUE; shift ;;
    --benchmark-glmnet) CKM_BENCHMARK_GLMNET=TRUE; shift ;;
    --met-strict-exclude-lipid-block) CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --prevent-batch-size|--prevent_batch_size) export CKM_PREVENT_BATCH_SIZE="$2"; shift 2 ;;
    --feature-chunk-size|--feature_chunk_size) export CKM_FEATURE_CHUNK_SIZE="$2"; shift 2 ;;
    --ncore) export CKM_NCORE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# A bare rebuild is the one-time construction phase: create both CKM stage
# assignments, then stop before summaries and omic analysis.
if [[ "${CKM_REBUILD^^}" == "TRUE" && "${START_STEP_EXPLICIT}" == "FALSE" && "${END_STEP_EXPLICIT}" == "FALSE" ]]; then
  START_STEP=build_ckm
  END_STEP=build_ckm
fi

biom_key="${BIOM_REQUEST,,}"
biom_key="${biom_key//[[:space:]]/}"
case "${biom_key}" in
  prot) BIOM_RUNS=("prot") ;;
  met) BIOM_RUNS=("met") ;;
  prot,met|met,prot) BIOM_RUNS=("prot" "met") ;;
  prot+met|met+prot|prot.met|met.prot) BIOM_RUNS=("prot+met") ;;
  *) echo "Error: --biom must be prot, met, prot,met, or prot+met; got ${BIOM_REQUEST}" >&2; exit 1 ;;
esac


export START_STEP END_STEP CKM_SCAN_ALL_BIOM CKM_REBUILD CKM_STAGE4_REQUIRE_RISK CKM_MERGE_ALL
export CKM_ETHNIC CKM_CLOCK_REFERENCE CKM_CLOCK_ADJUST_MODE
export CKM_SKIP_PREDICTION CKM_BENCHMARK_GLMNET CKM_PRED_BASE_MODE SEED CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK
export DIR0="${DIR0:-/mnt/d}"
export PHEDIR="${PHEDIR:-${DIR0}/data/ukb/phe}"
export SCRIPT_DIR="${SCRIPT_DIR:-${DIR0}/scripts/0f}"
export CKM_OUTDIR="${CKM_OUTDIR:-${DIR0}/analysis/ckm}"
export CKM_AUDIT_DIR="${CKM_AUDIT_DIR:-${CKM_OUTDIR}/audit}"

privacy_cleanup_completed_run=FALSE
cleanup_ukb_individual_outputs() {
  local status=$? f
  trap - EXIT
  if (( status != 0 )) || [[ "$privacy_cleanup_completed_run" != TRUE ]]; then
    echo "CKM did not reach full completion (exit=$status); preserving participant-level outputs for resume." >&2
    exit "$status"
  fi
  case "${CKM_OUTDIR}" in ""|/|/mnt|/mnt/|/mnt/d|/mnt/d/|/mnt/d/analysis|/mnt/d/analysis/)
    echo "PRIVACY ERROR: unsafe CKM_OUTDIR: '${CKM_OUTDIR}'" >&2; exit 70 ;;
  esac
  if [[ -d "${CKM_OUTDIR}" ]]; then
    echo "Privacy cleanup: removing CKM UKB participant-level outputs"
    while IFS= read -r -d '' f; do rm -f -- "$f"; done < <(find "${CKM_OUTDIR}" -type f \( \
      -name 'ckm_build_audit.tsv.gz' -o -name 'ckm_build_audit.xlsx' -o \
      -name 'stage_clock*.rds' -o -name '*biom_refine.oof.tsv.gz' -o \
      -name '*biom_refine.rds' -o -name 'prediction.rows.tsv.gz' -o \
      -name 'prediction.rds' -o -name 'ckm.all_results.rds' \) -print0)
  fi
  exit "$status"
}
trap cleanup_ukb_individual_outputs EXIT

# Create an output directory robustly.  A common failure mode is that a regular
# file or a broken symbolic link exists at the requested directory path.  In
# that case, preserve it by moving it aside, then create the directory.
ensure_output_dir() {
  local path="$1"
  local label="$2"
  local stamp backup parent

  if [[ -z "${path}" || "${path}" == "/" ]]; then
    echo "Error: unsafe ${label} path: '${path}'" >&2
    exit 1
  fi

  parent="$(dirname -- "${path}")"
  if [[ -e "${parent}" && ! -d "${parent}" ]]; then
    echo "Error: parent of ${label} exists but is not a directory: ${parent}" >&2
    ls -ld -- "${parent}" >&2 || true
    exit 1
  fi
  mkdir -p -- "${parent}"

  if [[ -L "${path}" && ! -e "${path}" ]]; then
    stamp="$(date +%Y%m%d_%H%M%S)"
    backup="${path}.broken_symlink.${stamp}"
    echo "Warning: ${label} is a broken symlink; moving it to: ${backup}" >&2
    mv -- "${path}" "${backup}"
  elif [[ -e "${path}" && ! -d "${path}" ]]; then
    stamp="$(date +%Y%m%d_%H%M%S)"
    backup="${path}.blocking_file.${stamp}"
    echo "Warning: ${label} path exists as a non-directory; moving it to: ${backup}" >&2
    ls -ld -- "${path}" >&2 || true
    mv -- "${path}" "${backup}"
  fi

  mkdir -p -- "${path}"
  if [[ ! -d "${path}" ]]; then
    echo "Error: failed to create ${label}: ${path}" >&2
    exit 1
  fi
  if [[ ! -w "${path}" ]]; then
    echo "Error: ${label} is not writable: ${path}" >&2
    ls -ld -- "${path}" >&2 || true
    exit 1
  fi
}

ensure_output_dir "${CKM_OUTDIR}" "CKM output directory"

export CKM_ICD10_LIST="${CKM_ICD10_LIST:-${PHEDIR}/common/icd10.lst}"
export CKM_PROT_META="${CKM_PROT_META:-${DIR0}/data.BIG/gwas/ppp/ppp_3k_b38.bed}"
export CKM_MET_META="${CKM_MET_META:-${PHEDIR}/common/met.lst}"

RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript || true)}"
if [[ -z "${RSCRIPT_BIN}" ]]; then
  echo "Error: cannot find Rscript. Set RSCRIPT_BIN=/path/to/Rscript." >&2
  exit 1
fi
if ! "${RSCRIPT_BIN}" --version 2>&1 | grep -qi '^Rscript'; then
  echo "Error: RSCRIPT_BIN does not appear to be Rscript: ${RSCRIPT_BIN}" >&2
  echo "Unset RSCRIPT_BIN or set it to the full path returned by: command -v Rscript" >&2
  exit 1
fi

SCRIPT_DIR_THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CKM_SCRIPT_DIR_THIS="${SCRIPT_DIR_THIS}"
base="$(basename "$0" .sh)"
R_FILE="${SCRIPT_DIR_THIS}/${base}.R"
if [[ ! -s "${R_FILE}" ]]; then
  R_FILE="${SCRIPT_DIR_THIS}/ckm.R"
fi
if [[ ! -s "${R_FILE}" ]]; then
  echo "Error: cannot find ${base}.R or ckm.R next to this shell script" >&2
  exit 1
fi

stage3_method_display="${CKM_STAGE3_METHOD:-SCORE2}"

cat <<INFO
CKM 2026 risk-age staging and Yin-Yang pipeline
CKM_TARGET=fixed ckm4 (first CHD/HF/stroke/PAD/AF)
BIOM_REQUEST=${BIOM_REQUEST}
BIOM_RUNS=${BIOM_RUNS[*]}
START_STEP=${START_STEP}
END_STEP=${END_STEP}
CKM_SCAN_ALL_BIOM=${CKM_SCAN_ALL_BIOM}
CKM_REBUILD=${CKM_REBUILD}
CKM_STAGE4_REQUIRE_RISK=${CKM_STAGE4_REQUIRE_RISK}
CKM_STAGE3_METHOD=${stage3_method_display}
CKM_STAGE4_SOURCE=${CKM_STAGE4_SOURCE:-all}
CKM_CLOCK_REFERENCE=${CKM_CLOCK_REFERENCE}
CKM_CLOCK_ADJUST_MODE=${CKM_CLOCK_ADJUST_MODE}
CKM_PRED_BASE_MODE=${CKM_PRED_BASE_MODE}
CKM_PREVENT_MODEL=${CKM_PREVENT_MODEL:-base}
CKM_TG_FASTING_CUTOFF=${CKM_TG_FASTING_CUTOFF:-1.69}
CKM_TG_NONFASTING_CUTOFF=${CKM_TG_NONFASTING_CUTOFF:-1.98}
CKM_RISK_AGE_REFERENCE=${CKM_RISK_AGE_REFERENCE:-60}
CKM_LAG_MIN_RETENTION=${CKM_LAG_MIN_RETENTION:-0.50}
CKM_LAG_MIN_LOGHR=${CKM_MARKER_MIN_COX_LOGHR:-log(1.05)}
CKM_SKIP_PREDICTION=${CKM_SKIP_PREDICTION}
CKM_BENCHMARK_GLMNET=${CKM_BENCHMARK_GLMNET}
CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK=${CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK}
CKM_PREVENT_BATCH_SIZE=${CKM_PREVENT_BATCH_SIZE:-5000}
CKM_FEATURE_CHUNK_SIZE=${CKM_FEATURE_CHUNK_SIZE:-auto}
CKM_NCORE=${CKM_NCORE:-auto_conservative}
CKM_CLOCK_DRAW_REF_N=${CKM_CLOCK_DRAW_REF_N:-1000}
DIR0=${DIR0}
PHEDIR=${PHEDIR}
CKM_OUTDIR=${CKM_OUTDIR}
CKM_AUDIT_DIR=${CKM_AUDIT_DIR} (created only when ckm.rds is built/rebuilt)
CKM_ICD10_LIST=${CKM_ICD10_LIST}
CKM_STAGE4_COLUMNS=fod_icd10_cvd_cad,fod_icd10_mi,[fod_icd10_ihd],fod_icd10_cvd_hfail,fod_icd10_stroke,fod_icd10_cvd_pad,fod_icd10_cvd_afib
Rscript=${RSCRIPT_BIN}
R file=${R_FILE}
INFO

rebuild_requested="${CKM_REBUILD}"
run_index=0
for biom_run in "${BIOM_RUNS[@]}"; do
  export CKM_BIOM="${biom_run}"
  r_args=(--biom "${CKM_BIOM}" --start-step "${START_STEP}")
  if [[ ${run_index} -eq 0 ]]; then
    export CKM_REBUILD="${rebuild_requested}"
    export CKM_RUN_SHARED_OUTPUTS=TRUE
    if [[ "${rebuild_requested^^}" == "TRUE" ]]; then
      r_args+=(--rebuild-ckm)
    fi
  else
    # ckm.rds is omic-independent; do not rebuild it again for the second layer.
    export CKM_REBUILD=FALSE
    export CKM_RUN_SHARED_OUTPUTS=FALSE
  fi
  echo
  echo "========== RUN CKM BIOMARKER MODE: ${CKM_BIOM} =========="
  echo "CKM_RUN_SHARED_OUTPUTS=${CKM_RUN_SHARED_OUTPUTS}"
  "${RSCRIPT_BIN}" "${R_FILE}" "${r_args[@]}"
  echo "========== DONE CKM BIOMARKER MODE: ${CKM_BIOM} =========="
  run_index=$((run_index + 1))
done

privacy_cleanup_completed_run=TRUE
echo "Done. Outputs are under ${CKM_OUTDIR}."

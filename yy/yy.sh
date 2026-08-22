#!/usr/bin/env bash
set -euo pipefail



usage() {
  cat <<'HELP'
Usage:
  ./yy.sh [options]

Main options:
  -Y, --trait TRAIT        One outcome, e.g. cvd_cad
      --traits CSV         Multiple outcomes, e.g. cvd_cad,cvd_hfail,cvd_htn,stroke,t2dm,ckd,acd,depress,ra,amr,death
      --biom X             Biomarker mode:
                             prot      protein-only analysis
                             met       metabolite-only analysis
                             prot,met  run prot and met independently
                             prot+met  merge protein and metabolite data in one integrated analysis
                           Default: prot
  -s, --start-step STEP    Start from a step. Default: all
      --date-source MODE   icd10 or all. Default: icd10
      --require-aod        Stop if real AOD/aodtxt columns cannot be found
      --aod-qc-traits CSV  Root-level aodtxt QC traits. Default includes mi,cvd_cad,cvd_hfail,cvd_afib,t2dm,ckd
      --fig-standardize X  control or case. Default: control
      --timeline-adjust X  raw, age_sex_PC, age_sex_PC_tdi, le8, drugs, genetic, full. Default: age_sex_PC_tdi
      --context-proteins X all or selected for root Fig2 category composition. Default: all
      --biomarkers CSV     Biomarkers to display in QC/timeline figures; names may come from either omic layer
      --extra-biomarkers CSV Additional biomarkers for root QC figures
      --prot-meta FILE     Protein name/group BED. Default: /mnt/d/data.BIG/gwas/ppp/ppp_3k_b38.bed
      --met-meta FILE      Metabolite name/group list. Default: /mnt/d/data/ukb/phe/common/met.lst
      --scan-all-proteins  Scan all available biomarkers in disease Fig4/Fig6-Fig7; YY timeline still shows prioritized biomarkers
      --max-proteins N     Limit biomarker scan for debugging. 0 = no limit
      --yy-plot-top-n N    Number of pWAS-prioritized biomarkers in YY timeline plots. Default: 6
      --yy-plot-always CSV Always include these biomarkers when present. Default: GDF15,PCSK9
      --no-prediction      Skip YY-informed prediction
      --pred-folds N       Outer CV folds for prediction. Default: 10
      --pred-screen-n N    Incident pWAS screen size for prediction. Default: 700
      --pred-top-n N       YY-ranked biomarker set size for prediction. Default: 300
      --pred-max-proteins N Limit biomarkers entering prediction for debugging. 0 = no limit
      --pred-lambda X      glmnet lambda: 1se or min. Default: 1se
      --pred-min-yin-n N   Minimum Yin participants required for prediction. Default: 1000
      --pred-min-events N  Minimum incident events required for prediction. Default: 50
                           LE8/base variables are taken exactly from vars.basic and vars.le8 in 0phe.f.R
      --pred-base-mode X   demog, le8, or full. Default: le8
      --center-mode X      none or factor. Default: none
      --pred-adaptive-eta CSV Adaptive-penalty strengths. Default: 0,0.25,0.5,1
      --pred-layer-balance X none, proportional, sqrt, or equal. Default: sqrt for prot+met
      --pred-layer-ablation In prot+met mode, also fit protein-only and metabolite-only Yin models
      --no-pred-layer-ablation Disable layer-ablation models
      --pred-legacy-complex Also run older YY score-block/adaptive-glmnet diagnostic models. Default: off
      --pred-residualize   Residualize biomarkers on training-fold base covariates. Default: on
      --no-pred-residualize Disable training-fold biomarker residualization
      --pred-block-size N  Biomarker block size for residualization. Default: 250
      --lag-years CSV      Lag years for disease Fig5. Default: 0,2,5,10
      --specificity-outcomes CSV Outcomes for disease Fig4 specificity map
      --min-bin-n N        Minimum N per timeline bin. Default: 20
  -h, --help               Show this message

Output structure:
  --biom prot      -> YY_OUTDIR/prot/
  --biom met       -> YY_OUTDIR/met/
  --biom prot,met  -> both YY_OUTDIR/prot/ and YY_OUTDIR/met/ as independent runs
  --biom prot+met  -> YY_OUTDIR/prot+met/ using the intersection cohort and merged feature matrix

Examples:
  ./yy.sh
  ./yy.sh --biom prot,met --traits cvd_cad,stroke,mi,cvd_htn,copd,asthma,ra,depress,acd,ad,amr --require-aod --scan-all-proteins --pred-folds 10 --pred-base-mode le8 --center-mode none
  ./yy.sh --biom prot+met --trait cvd_cad --pred-layer-ablation
  YY_DATE_SOURCE=all ./yy.sh --trait cvd_cad
HELP
}

START_STEP="${START_STEP:-all}"
Y_CSV="${Ys:-${Y:-cvd_cad}}"
YY_BIOM_REQUEST="${YY_BIOM:-prot}"
YY_DATE_SOURCE="${YY_DATE_SOURCE:-icd10}"
YY_FIG_STANDARDIZE="${YY_FIG_STANDARDIZE:-${YY_FIG1_STANDARDIZE:-control}}"
YY_REQUIRE_AOD="${YY_REQUIRE_AOD:-FALSE}"
YY_TIMELINE_ADJUST="${YY_TIMELINE_ADJUST:-age_sex_PC_tdi}"
YY_CONTEXT_PROTEINS="${YY_CONTEXT_PROTEINS:-all}"
YY_SKIP_PREDICTION="${YY_SKIP_PREDICTION:-FALSE}"
YYP_RESIDUALIZE="${YYP_RESIDUALIZE:-TRUE}"
YYP_BASE_MODE="${YYP_BASE_MODE:-le8}"
YYP_ADAPTIVE_ETA="${YYP_ADAPTIVE_ETA:-0,0.25,0.5,1}"
YYP_LAYER_BALANCE="${YYP_LAYER_BALANCE:-auto}"
YYP_LAYER_ABLATION="${YYP_LAYER_ABLATION:-TRUE}"
YYP_LEGACY_COMPLEX="${YYP_LEGACY_COMPLEX:-FALSE}"

yy_privacy_root="${YY_OUTDIR:-/mnt/d/analysis/yy}"
privacy_cleanup_completed_run=FALSE
cleanup_ukb_individual_outputs() {
  local status=$? f
  trap - EXIT
  if (( status != 0 )) || [[ "$privacy_cleanup_completed_run" != TRUE ]]; then
    echo "YY did not reach full completion (exit=$status); preserving participant-level prediction data for resume." >&2
    exit "$status"
  fi
  case "$yy_privacy_root" in ""|/|/mnt|/mnt/|/mnt/d|/mnt/d/|/mnt/d/analysis|/mnt/d/analysis/)
    echo "PRIVACY ERROR: unsafe YY_OUTDIR: '$yy_privacy_root'" >&2; exit 70 ;;
  esac
  if [[ -d "$yy_privacy_root" ]]; then
    echo "Privacy cleanup: removing YY UKB participant-level prediction data"
    while IFS= read -r -d '' f; do rm -f -- "$f"; done < <(find "$yy_privacy_root" -type f \( \
      -name '*prediction_rows.tsv.gz' -o -name '*risk_deciles.*.tsv.gz' -o \
      -name '*prediction.rds' -o -name 'summary.*.rds' -o \
      -name 'yy_all_results.*.rds' -o -name '*multi_date_examples*.tsv.gz' \) -print0)
  fi
  exit "$status"
}
trap cleanup_ukb_individual_outputs EXIT

while [[ $# -gt 0 ]]; do
  case "$1" in
    -Y|--trait) Y_CSV="$2"; shift 2 ;;
    --traits|--Ys|--ys) Y_CSV="$2"; shift 2 ;;
    --biom) YY_BIOM_REQUEST="$2"; shift 2 ;;
    -s|--start-step|--start_step) START_STEP="$2"; shift 2 ;;
    --date-source) YY_DATE_SOURCE="$2"; shift 2 ;;
    --fig-standardize|--fig_standardize) YY_FIG_STANDARDIZE="$2"; shift 2 ;;
    --aod-qc-traits|--aod_qc_traits) export YY_AOD_QC_TRAITS="$2"; shift 2 ;;
    --timeline-adjust|--timeline_adjust) YY_TIMELINE_ADJUST="$2"; shift 2 ;;
    --context-proteins|--context_proteins) YY_CONTEXT_PROTEINS="$2"; shift 2 ;;
    --biomarkers|--proteins) export YY_BIOMARKERS="$2"; shift 2 ;;
    --extra-biomarkers|--extra-proteins) export YY_EXTRA_BIOMARKERS="$2"; shift 2 ;;
    --prot-meta|--prot_meta) export YY_PROT_META="$2"; shift 2 ;;
    --met-meta|--met_meta) export YY_MET_META="$2"; shift 2 ;;
    --root-balance-top-n|--root_balance_top_n) export YY_ROOT_BALANCE_TOP_N="$2"; shift 2 ;;
    --scan-all-proteins) export YY_SCAN_ALL_PROTEINS=TRUE; shift ;;
    --no-prediction) YY_SKIP_PREDICTION=TRUE; shift ;;
    --require-aod) YY_REQUIRE_AOD=TRUE; shift ;;
    --max-proteins) export YY_MAX_PROTEINS="$2"; shift 2 ;;
    --yy-plot-top-n|--yy_plot_top_n) export YY_PLOT_TOP_N="$2"; shift 2 ;;
    --yy-plot-always|--yy_plot_always) export YY_PLOT_ALWAYS_PROTEINS="$2"; shift 2 ;;
    --pred-folds|--pred_folds) export YYP_FOLDS="$2"; shift 2 ;;
    --pred-screen-n|--pred_screen_n) export YYP_SCREEN_N="$2"; shift 2 ;;
    --pred-top-n|--pred_top_n) export YYP_TOP_N="$2"; shift 2 ;;
    --pred-max-proteins|--pred_max_proteins) export YYP_MAX_PROTEINS="$2"; shift 2 ;;
    --pred-lambda|--pred_lambda) export YYP_LAMBDA="$2"; shift 2 ;;
    --pred-min-yin-n|--pred_min_yin_n) export YYP_MIN_YIN_N="$2"; shift 2 ;;
    --pred-min-events|--pred_min_events) export YYP_MIN_EVENTS="$2"; shift 2 ;;
    --pred-base-mode|--pred_base_mode) YYP_BASE_MODE="$2"; shift 2 ;;
    --pred-adaptive-eta|--pred_adaptive_eta) YYP_ADAPTIVE_ETA="$2"; shift 2 ;;
    --pred-layer-balance|--pred_layer_balance) YYP_LAYER_BALANCE="$2"; shift 2 ;;
    --pred-layer-ablation|--pred_layer_ablation) YYP_LAYER_ABLATION=TRUE; shift ;;
    --no-pred-layer-ablation|--no_pred_layer_ablation) YYP_LAYER_ABLATION=FALSE; shift ;;
    --pred-legacy-complex|--pred_legacy_complex) YYP_LEGACY_COMPLEX=TRUE; shift ;;
    --no-pred-legacy-complex|--no_pred_legacy_complex) YYP_LEGACY_COMPLEX=FALSE; shift ;;
    --pred-residualize) YYP_RESIDUALIZE=TRUE; shift ;;
    --no-pred-residualize) YYP_RESIDUALIZE=FALSE; shift ;;
    --pred-block-size|--pred_block_size) export YYP_BLOCK_SIZE="$2"; shift 2 ;;
    --lag-years|--lag_years) export YY_LAG_YEARS="$2"; shift 2 ;;
    --specificity-outcomes|--specificity_outcomes) export YY_SPECIFICITY_OUTCOMES="$2"; shift 2 ;;
    --min-bin-n|--min_bin_n) export YY_MIN_BIN_N="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# Keep comma and plus semantics distinct:
#   prot,met = two independent runs; prot+met = one merged run.
biom_key="${YY_BIOM_REQUEST,,}"
biom_key="${biom_key//[[:space:]]/}"
case "${biom_key}" in
  prot) BIOM_RUNS=("prot") ;;
  met) BIOM_RUNS=("met") ;;
  prot,met|met,prot) BIOM_RUNS=("prot" "met") ;;
  prot+met|met+prot|prot.met|met.prot) BIOM_RUNS=("prot+met") ;;
  *)
    echo "Error: --biom must be prot, met, prot,met, or prot+met; got: ${YY_BIOM_REQUEST}" >&2
    exit 1
    ;;
esac

export START_STEP
export Ys="$Y_CSV"
export YY_DATE_SOURCE YY_FIG_STANDARDIZE YY_REQUIRE_AOD YY_TIMELINE_ADJUST YY_CONTEXT_PROTEINS YY_SKIP_PREDICTION YYP_RESIDUALIZE
export YYP_BASE_MODE YYP_ADAPTIVE_ETA YYP_LAYER_BALANCE YYP_LAYER_ABLATION YYP_LEGACY_COMPLEX
export DIR0="${DIR0:-/mnt/d}"
export PHEDIR="${PHEDIR:-${DIR0}/data/ukb/phe}"
export SCRIPT_DIR="${SCRIPT_DIR:-${DIR0}/scripts/0f}"
export YY_OUTDIR="${YY_OUTDIR:-${DIR0}/analysis/yy}"
export SEED="${SEED:-2026}"

RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript || true)}"
if [[ -z "${RSCRIPT_BIN}" ]]; then
  echo "Error: cannot find Rscript. Set RSCRIPT_BIN=/path/to/Rscript." >&2
  exit 1
fi

SCRIPT_DIR_THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base="$(basename "$0" .sh)"
R_FILE="${SCRIPT_DIR_THIS}/${base}.R"
if [[ ! -s "${R_FILE}" ]]; then
  R_FILE="${SCRIPT_DIR_THIS}/yy.R"
fi
if [[ ! -s "${R_FILE}" ]]; then
  echo "Error: cannot find ${base}.R or yy.R next to this shell script." >&2
  exit 1
fi

echo "YY v18 biomarker + prediction pipeline"
echo "Ys=${Ys}"
echo "YY_BIOM_REQUEST=${YY_BIOM_REQUEST}"
echo "BIOM_RUNS=${BIOM_RUNS[*]}"
echo "START_STEP=${START_STEP}"
echo "YY_DATE_SOURCE=${YY_DATE_SOURCE}"
echo "YY_AOD=precomputed only; YY_REQUIRE_AOD=${YY_REQUIRE_AOD}"
echo "YY_AOD_QC_TRAITS=${YY_AOD_QC_TRAITS:-auto}"
echo "YY_FIG_STANDARDIZE=${YY_FIG_STANDARDIZE}"
echo "YY_TIMELINE_ADJUST=${YY_TIMELINE_ADJUST}"
echo "YY_CONTEXT_PROTEINS=${YY_CONTEXT_PROTEINS}"
echo "YY_SCAN_ALL_PROTEINS=${YY_SCAN_ALL_PROTEINS:-FALSE}"
echo "YY_SKIP_PREDICTION=${YY_SKIP_PREDICTION}"
echo "YY_PLOT_TOP_N=${YY_PLOT_TOP_N:-6}; YY_PLOT_ALWAYS_PROTEINS=${YY_PLOT_ALWAYS_PROTEINS:-GDF15,PCSK9}"
echo "YYP_FOLDS=${YYP_FOLDS:-10}; YYP_SCREEN_N=${YYP_SCREEN_N:-700}; YYP_TOP_N=${YYP_TOP_N:-300}; YYP_MAX_PROTEINS=${YYP_MAX_PROTEINS:-0}; YYP_MIN_YIN_N=${YYP_MIN_YIN_N:-1000}; YYP_MIN_EVENTS=${YYP_MIN_EVENTS:-50}"
echo "YYP_RESIDUALIZE=${YYP_RESIDUALIZE}; YYP_BLOCK_SIZE=${YYP_BLOCK_SIZE:-250}"
echo "YYP_BASE_MODE=${YYP_BASE_MODE}; YYP_ADAPTIVE_ETA=${YYP_ADAPTIVE_ETA}; YYP_LAYER_BALANCE=${YYP_LAYER_BALANCE}; YYP_LAYER_ABLATION=${YYP_LAYER_ABLATION}; YYP_LEGACY_COMPLEX=${YYP_LEGACY_COMPLEX}"
echo "YY_PROT_META=${YY_PROT_META:-${DIR0}/data.BIG/gwas/ppp/ppp_3k_b38.bed}"
echo "YY_MET_META=${YY_MET_META:-${DIR0}/data/ukb/phe/common/met.lst}"
echo "DIR0=${DIR0}"
echo "PHEDIR=${PHEDIR}"
echo "SCRIPT_DIR=${SCRIPT_DIR}"
echo "YY_OUTDIR=${YY_OUTDIR}"
echo "Rscript=${RSCRIPT_BIN}"
echo "R file=${R_FILE}"

for biom_run in "${BIOM_RUNS[@]}"; do
  export YY_BIOM="${biom_run}"
  echo
  echo "========== RUN BIOMARKER MODE: ${YY_BIOM} =========="
  echo "Output: ${YY_OUTDIR}/${YY_BIOM}"
  "${RSCRIPT_BIN}" "${R_FILE}" --traits "${Ys}" --start_step "${START_STEP}"
  echo "========== DONE BIOMARKER MODE: ${YY_BIOM} =========="
done

privacy_cleanup_completed_run=TRUE
echo "Done. Outputs are under ${YY_OUTDIR}."

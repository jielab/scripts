#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'HELP'
Usage: ./pfi.sh COMMAND [options]

Commands:
  all              Build source data and stop for manual expert labels
  after08          Continue automatically from completed steps 06/07
  download-pubmed  Download PubMed baseline/update archives
  parse-pubmed     Parse PubMed archives
  download-pmc     Download PMC Open Access XML archives
  parse-pmc        Parse PMC XML archives
  build-master     Build articles_master
  make-blinded     Build the blinded-text dataset
  features         Extract basic features
  normalize-rank   Normalize a real or explicitly allowed proxy journal rank
  review-inputs    Create blinded expert-review inputs
  auto-score       Create weak-supervision labels
  prepare-training Build the model-training dataset
  train            Train PubMedBERT/BiomedBERT
  predict          Generate article predictions
  compute          Compute the PFI index
  analysis         Run the R analysis and figures
  sample-llm       Sample an optional LLM review batch
  audit            Audit strong outliers

Common options:
  --analysis-root DIR               Analysis root [/mnt/d/analysis/pfi]
  --raw-root DIR                    Raw archive root [/mnt/g/pfi]
  --python FILE                     Python executable [python]
  --rscript FILE                    Rscript executable [Rscript]
  --test-mode TRUE|FALSE            Small engineering run [FALSE]
  --journal-rank FILE               Real journal-rank CSV
  --allow-proxy-journal-rank BOOL   Permit a non-scientific proxy rank [FALSE]
  --review-rows INT                 Review-set size [1000]
  --review-max-scan-rows INT        Maximum rows scanned for review sampling
  --max-train-rows INT              Maximum training rows [1000]
  --pred-sample-rows INT            Prediction sample size; 0 means all [1000000]
  --pred-limit INT                  Hard prediction limit; 0 means none [0]
  --pred-fresh-sample TRUE|FALSE    Rebuild the prediction sample [FALSE]
  --rows INT                        Rows for sample-llm [1000]
  --help                            Show this message

Examples:
  ./pfi.sh all
  ./pfi.sh after08 --journal-rank /mnt/d/scripts/pfi/journal_rank_real.csv
  ./pfi.sh after08 --review-rows 100 --max-train-rows 100 --pred-sample-rows 10000 --pred-limit 10000 --allow-proxy-journal-rank TRUE
HELP
}

die() { echo "ERROR: $*" >&2; exit 2; }

bool01() {
  case "${1,,}" in
    true|1|yes|y) printf '1\n' ;;
    false|0|no|n) printf '0\n' ;;
    *) return 2 ;;
  esac
}

[[ $# -gt 0 ]] || { usage; exit 2; }
case "$1" in -h|--help|help) usage; exit 0;; esac
command_name=$1
shift

analysis_root=/mnt/d/analysis/pfi
raw_root=/mnt/g/pfi
python_bin=python
rscript_bin=Rscript
test_mode=0
journal_rank=""
allow_proxy=0
review_rows=1000
review_max_scan_rows=""
max_train_rows=1000
pred_sample_rows=1000000
pred_limit=0
pred_fresh_sample=0
sample_rows=1000

while (( $# )); do
  case "$1" in
    --analysis-root) analysis_root=${2:?}; shift 2 ;;
    --raw-root) raw_root=${2:?}; shift 2 ;;
    --python) python_bin=${2:?}; shift 2 ;;
    --rscript) rscript_bin=${2:?}; shift 2 ;;
    --test-mode) test_mode=$(bool01 "${2:?}") || die "--test-mode requires TRUE or FALSE"; shift 2 ;;
    --journal-rank) journal_rank=${2:?}; shift 2 ;;
    --allow-proxy-journal-rank) allow_proxy=$(bool01 "${2:?}") || die "--allow-proxy-journal-rank requires TRUE or FALSE"; shift 2 ;;
    --review-rows) review_rows=${2:?}; shift 2 ;;
    --review-max-scan-rows) review_max_scan_rows=${2:?}; shift 2 ;;
    --max-train-rows) max_train_rows=${2:?}; shift 2 ;;
    --pred-sample-rows) pred_sample_rows=${2:?}; shift 2 ;;
    --pred-limit) pred_limit=${2:?}; shift 2 ;;
    --pred-fresh-sample) pred_fresh_sample=$(bool01 "${2:?}") || die "--pred-fresh-sample requires TRUE or FALSE"; shift 2 ;;
    --rows) sample_rows=${2:?}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done

for item in "$review_rows" "$max_train_rows" "$pred_sample_rows" "$pred_limit" "$sample_rows"; do
  [[ $item =~ ^[0-9]+$ ]] || die "row and limit options must be non-negative integers"
done
[[ -z $review_max_scan_rows || $review_max_scan_rows =~ ^[1-9][0-9]*$ ]] || die "--review-max-scan-rows must be a positive integer"

export ANALYSIS_ROOT=$analysis_root RAW_ROOT=$raw_root
export PFI_PYTHON_BIN=$python_bin PFI_RSCRIPT_BIN=$rscript_bin PFI_TEST_MODE=$test_mode
export PFI_ALLOW_PROXY_JOURNAL_RANK=$allow_proxy PFI_REVIEW_ROWS=$review_rows
export PFI_MAX_TRAIN_ROWS=$max_train_rows PFI_PRED_SAMPLE_ROWS=$pred_sample_rows
export PFI_PRED_LIMIT=$pred_limit PFI_PRED_FRESH_SAMPLE=$pred_fresh_sample
[[ -z $review_max_scan_rows ]] || export PFI_REVIEW_MAX_SCAN_ROWS=$review_max_scan_rows

cd "$SCRIPT_DIR"
case "$command_name" in
  all) bash ./00_run_pfi_all.sh "$journal_rank" ;;
  after08) bash ./00_run_pfi_after08_auto.sh "$journal_rank" ;;
  download-pubmed) bash ./02_download_pubmed.sh --analysis-root "$analysis_root" --raw-root "$raw_root" --python "$python_bin" --test-mode "$([[ $test_mode == 1 ]] && echo TRUE || echo FALSE)" ;;
  parse-pubmed) bash ./03_parse_pubmed.sh ;;
  download-pmc) bash ./04_download_pmc_oa_bulk_xml.sh ;;
  parse-pmc) bash ./05_parse_pmc_oa_bulk_xml.sh ;;
  build-master) bash ./06_build_master_dataset.sh ;;
  make-blinded) bash ./07_make_blinded_text.sh ;;
  features) bash ./08_extract_basic_features.sh ;;
  normalize-rank) bash ./08b_normalize_journal_rank.sh "$journal_rank" ;;
  review-inputs) bash ./09a_make_pfi_review_inputs.sh ;;
  auto-score) bash ./09a_auto_score_pfi_review_dataset.sh ;;
  prepare-training) bash ./09a_prepare_pfi_training_dataset.sh ;;
  train) bash ./09b_train_pubmedbert_pfi_model.sh ;;
  predict) bash ./09c_predict_pubmedbert_pfi.sh ;;
  compute) bash ./09d_compute_pfi_index.sh ;;
  analysis) bash ./10_run_pfi_analysis_R.sh ;;
  sample-llm) bash ./11_sample_llm_review_batch.sh "$sample_rows" ;;
  audit) bash ./12_audit_strong_outliers.sh ;;
  *) die "unknown command: $command_name" ;;
esac

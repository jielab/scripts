#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
PFI: publication data preparation, scoring and index analysis.

Usage:
  ./pfi.sh run|MODULE [options]

Modules:
  run                              Prepare data and review inputs; training requires completed review scores.
  after08                          Continue from prepared master data and blinded text.
  download-pubmed | parse-pubmed   Download or parse PubMed records.
  download-pmc | parse-pmc         Download or parse PMC records.
  build-master | make-blinded | features Prepare the article dataset and features.
  normalize-rank | review-inputs | auto-score Prepare ranks and review scores.
  prepare-training | train | predict Prepare training data, train or predict.
  compute | analysis | sample-llm | audit Compute the index, analyze, sample or audit.

Options:
  --plan-only                      Write a local/LSF .cmd without running
  --job-hours N --queue NAME       
  --job-gpus N                     Explicit GPU allocation for train/predict only [0]
  --dir-analysis-root DIR          Results [${ANALYSIS_ROOT}/pfi]
  --dir-raw-root DIR               Raw data [${DIR_BIG}/pfi]
  --journal-rank FILE              Real journal rank data
  --pred-sample-rows N Prediction sample size [1000000]; 0 for all 
  --python FILE --rscript FILE     Interpreters
  -h, --help                       Show this help.

Examples:
  ./pfi.sh run --journal-rank journal_rank_real.csv --plan-only
  ./pfi.sh prepare-training
  ./pfi.sh train --job-gpus 1 --queue 2v100-32-sc
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-auto}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve --dir-dat/--dir-big/--dir-tmp through the shared path table.
source "$SCRIPTS_ROOT/0f/project_entry.sh"
project_workflow_init "$@" || exit $?
set -- "${WORKFLOW_ARGS[@]}"

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
source "${SCRIPTS_ROOT}/0f/console.sh"
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)


# Command-line help
usage() {
	bash "${BASH_SOURCE[0]}" --help
}

die() {
	echo "ERROR: $*" >&2
	exit 2
}

bool01() {
	case "${1,,}" in
		true | 1 | yes | y) printf '1\n' ;;
		false | 0 | no | n) printf '0\n' ;;
		*) return 2 ;;
	esac
}

[[ $# -gt 0 ]] || {
	usage
	exit 2
}
case "$1" in -h | --help | help)
	usage
	exit 0
	;;
esac
command_name=$1
shift

analysis_root=${PFI_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/pfi}
raw_root=${RAW_ROOT:-${DIR_BIG}/pfi}
python_bin=${PFI_PYTHON_BIN:-python}
rscript_bin=${PFI_RSCRIPT_BIN:-Rscript}
test_mode=${PFI_TEST_MODE:-0}
journal_rank=""
allow_proxy=0
review_rows=1000
review_max_scan_rows=""
max_train_rows=1000
pred_sample_rows=1000000
pred_limit=0
pred_fresh_sample=0
sample_rows=1000


# Command-line arguments
while (($#)); do
	case "$1" in
		--analysis-root)
			analysis_root=${2:?}
			shift 2
			;;
		--raw-root)
			raw_root=${2:?}
			shift 2
			;;
		--python)
			python_bin=${2:?}
			shift 2
			;;
		--rscript)
			rscript_bin=${2:?}
			shift 2
			;;
		--test-mode)
			test_mode=$(bool01 "${2:?}") || die "--test-mode requires TRUE or FALSE"
			shift 2
			;;
		--journal-rank)
			journal_rank=${2:?}
			shift 2
			;;
		--allow-proxy-journal-rank)
			allow_proxy=$(bool01 "${2:?}") || die "--allow-proxy-journal-rank requires TRUE or FALSE"
			shift 2
			;;
		--review-rows)
			review_rows=${2:?}
			shift 2
			;;
		--review-max-scan-rows)
			review_max_scan_rows=${2:?}
			shift 2
			;;
		--max-train-rows)
			max_train_rows=${2:?}
			shift 2
			;;
		--pred-sample-rows)
			pred_sample_rows=${2:?}
			shift 2
			;;
		--pred-limit)
			pred_limit=${2:?}
			shift 2
			;;
		--pred-fresh-sample)
			pred_fresh_sample=$(bool01 "${2:?}") || die "--pred-fresh-sample requires TRUE or FALSE"
			shift 2
			;;
		--rows)
			sample_rows=${2:?}
			shift 2
			;;
		-h | --help)
			usage
			exit 0
			;;
		*) die "unknown option: $1" ;;
	esac
done

for item in "$review_rows" "$max_train_rows" "$pred_sample_rows" "$pred_limit" "$sample_rows"; do
	[[ $item =~ ^[0-9]+$ ]] || die "row and limit options must be non-negative integers"
done
[[ -z $review_max_scan_rows || $review_max_scan_rows =~ ^[1-9][0-9]*$ ]] || die "--review-max-scan-rows must be a positive integer"

export PFI_ANALYSIS_ROOT=$analysis_root ANALYSIS_ROOT=$analysis_root RAW_ROOT=$raw_root
export PFI_PYTHON_BIN=$python_bin PFI_RSCRIPT_BIN=$rscript_bin PFI_TEST_MODE=$test_mode
export PFI_ALLOW_PROXY_JOURNAL_RANK=$allow_proxy PFI_REVIEW_ROWS=$review_rows
export PFI_MAX_TRAIN_ROWS=$max_train_rows PFI_PRED_SAMPLE_ROWS=$pred_sample_rows
export PFI_PRED_LIMIT=$pred_limit PFI_PRED_FRESH_SAMPLE=$pred_fresh_sample
[[ -z $review_max_scan_rows ]] || export PFI_REVIEW_MAX_SCAN_ROWS=$review_max_scan_rows

cd "$SCRIPT_DIR"
case "$command_name" in
	run) bash ./00_run_pfi_all.sh run "$journal_rank" ;;
	after08) bash ./00_run_pfi_after08_auto.sh run "$journal_rank" ;;
	download-pubmed) bash ./02_download_pubmed.sh run --analysis-root "$analysis_root" --raw-root "$raw_root" --python "$python_bin" --test-mode "$([[ $test_mode == 1 ]] && echo TRUE || echo FALSE)" ;;
	parse-pubmed) bash ./03_parse_pubmed.sh run ;;
	download-pmc) bash ./04_download_pmc_oa_bulk_xml.sh run ;;
	parse-pmc) bash ./05_parse_pmc_oa_bulk_xml.sh run ;;
	build-master) bash ./06_build_master_dataset.sh run ;;
	make-blinded) bash ./07_make_blinded_text.sh run ;;
	features) bash ./08_extract_basic_features.sh run ;;
	normalize-rank) bash ./08b_normalize_journal_rank.sh run "$journal_rank" ;;
	review-inputs) bash ./09a_make_pfi_review_inputs.sh run ;;
	auto-score) bash ./09a_auto_score_pfi_review_dataset.sh run ;;
	prepare-training) bash ./09a_prepare_pfi_training_dataset.sh run ;;
	train) bash ./09b_train_pubmedbert_pfi_model.sh run ;;
	predict) bash ./09c_predict_pubmedbert_pfi.sh run ;;
	compute) bash ./09d_compute_pfi_index.sh run ;;
	analysis) bash ./10_run_pfi_analysis_R.sh run ;;
	sample-llm) bash ./11_sample_llm_review_batch.sh run "$sample_rows" ;;
	audit) bash ./12_audit_strong_outliers.sh run ;;
	*) die "unknown command: $command_name" ;;
esac

#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# This entry initializes LOCAL paths and is not a supported --submit hpc pipeline.
# The bsub example above applies only after its worker paths/runtime are ported.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-local}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve --dir-dat/--dir-big/--dir-tmp through the shared path table.
source "$SCRIPTS_ROOT/0f/local_paths.sh"
local_workflow_init "$@" || exit $?
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


# 🚩 Command-line help
usage() {
	cat <<'HELP'
Modules:
    all  Prepare data and review inputs
    prepare-training / train / predict  Model workflow; requires completed review scores
    compute / analysis / audit  Index, analysis and validation

Main options:
    --dir-analysis-root DIR  Results [${ANALYSIS_ROOT}/pfi]
    --dir-raw-root DIR       Raw data [${DIR_BIG}/pfi]
    --journal-rank FILE  Real journal rank data
    --pred-sample-rows N Prediction sample size [1000000]; 0 for all
    --python FILE --rscript FILE  Interpreters
    -h, --help           Show options and examples below

Example:
    cd ${SCRIPTS_ROOT}/0LOCAL/pfi

    ./pfi.sh all --journal-rank ${SCRIPTS_ROOT}/0LOCAL/pfi/journal_rank_real.csv

    # Fill the review XLSX; save ${ANALYSIS_ROOT}/pfi/dat/paper_score_review_1000_scored.csv before continuing.
    ./pfi.sh prepare-training
    ./pfi.sh train
    ./pfi.sh predict --pred-sample-rows 0
    ./pfi.sh compute
    ./pfi.sh analysis
    ./pfi.sh audit

    ./pfi.sh all --journal-rank journal_rank_real.csv
    ./pfi.sh predict --pred-sample-rows 0
HELP
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

analysis_root=${ANALYSIS_ROOT}/pfi
raw_root=${DIR_BIG}/pfi
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


# 🚩 Command-line arguments
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

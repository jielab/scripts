#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/output_dirs.sh"
script_output_directory_args "${BASH_SOURCE[0]}" "$@" || { return 2 2>/dev/null || exit 2; }
set -- "${SCRIPT_OUTPUT_DIRECTORY_ARGS[@]}"

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
PFI: Run the post-feature workflow.

Usage:
  ./00_run_pfi_after08_auto.sh run [options]

Modules:
  run                              Run the post-feature workflow.

Permanent output directories (legacy options remain accepted):
  --dir-out-journal_rank DIR
                                  Durable journal-rank copy [${SCRIPTS_ROOT}/pfi/journal_rank_real.csv].

Options:
  --submit auto|local|hpc          Execution mode [auto].
  --plan-only                      Write a local/LSF .cmd without running
  --job-hours N --queue NAME       
  -h, --help                       Show this help.

Examples:
  ./00_run_pfi_after08_auto.sh run --plan-only
  ./00_run_pfi_after08_auto.sh run
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
if [[ ${1:-} == run ]]; then shift; fi

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			bash "${BASH_SOURCE[0]}" --help
			exit 0
			;;
	esac
fi


# 00_run_pfi_after08_auto
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR" "$OUTPUT_DIR/analysis" "$PFI_ANALYSIS_DATA_DIR"

LOG_FILE="$LOG_DIR/00_run_pfi_after08_auto.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "===== 00_run_pfi_after08_auto.sh started at $(date -Is) ====="
echo "This starts after 06/07 have produced articles_master and blinded_text."
echo "If no expert scored CSV exists, weak-supervision pseudo labels will be generated automatically."
echo "If no real journal-rank file exists, a temporary proxy rank will be generated from articles_master."

if [[ ! -e "$PFI_ARTICLES_MASTER" ]]; then
	echo "ERROR: articles_master is missing: $PFI_ARTICLES_MASTER" >&2
	echo "Run 06_build_master_dataset.sh first." >&2
	exit 1
fi
if [[ ! -e "$PFI_BLINDED_TEXT" ]]; then
	echo "ERROR: blinded_text is missing: $PFI_BLINDED_TEXT" >&2
	echo "Run 07_make_blinded_text.sh first." >&2
	exit 1
fi

# Features are needed for review sampling and pseudo labels.
if [[ ! -e "$PFI_ARTICLE_FEATURES" ]]; then
	echo "INFO: article features are missing; running 08_extract_basic_features.sh"
	bash "$PFI_ROOT/08_extract_basic_features.sh" run
fi

bash "$PFI_ROOT/08b_normalize_journal_rank.sh" run "${1:-}"
bash "$PFI_ROOT/09a_make_pfi_review_inputs.sh" run
bash "$PFI_ROOT/09a_auto_score_pfi_review_dataset.sh" run
bash "$PFI_ROOT/09a_prepare_pfi_training_dataset.sh" run "$PFI_JOURNAL_RANK_NORMALIZED" "$PFI_PAPER_SCORE_REVIEW_PSEUDO"
bash "$PFI_ROOT/09b_train_pubmedbert_pfi_model.sh" run
bash "$PFI_ROOT/09c_predict_pubmedbert_pfi.sh" run
bash "$PFI_ROOT/09d_compute_pfi_index.sh" run "$PFI_JOURNAL_RANK_NORMALIZED"
bash "$PFI_ROOT/10_run_pfi_analysis_R.sh" run
bash "$PFI_ROOT/12_audit_strong_outliers.sh" run

echo "===== 00_run_pfi_after08_auto.sh finished at $(date -Is) ====="
echo "Main outputs:"
echo "  $PFI_TRAINING_DATASET"
echo "  $PFI_PUBMEDBERT_PREDICTIONS"
echo "  $PFI_FAIRNESS_RESULTS_CSV"
echo "  $OUTPUT_DIR/analysis"

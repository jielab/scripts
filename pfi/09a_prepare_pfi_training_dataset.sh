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

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Prepare model training data

Main options:
    RANK_FILE SCORED_FILE  Optional positional input files
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09a_prepare_pfi_training_dataset.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09a_prepare_pfi_training_dataset
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"

rank_file="${1:-$PFI_JOURNAL_RANK_NORMALIZED}"
paper_score_file="${2:-$PFI_PAPER_SCORE_REVIEW_SCORED}"

# Ensure journal rank exists. If no real rank was provided, 08b can create a temporary proxy rank.
if [[ ! -s "$rank_file" ]]; then
	echo "INFO: normalized journal-rank file is missing; running 08b_normalize_journal_rank.sh" >&2
	bash "$PFI_ROOT/08b_normalize_journal_rank.sh"
	rank_file="$PFI_JOURNAL_RANK_NORMALIZED"
fi

# Ensure review input exists.
if [[ ! -s "$PFI_PAPER_SCORE_REVIEW_CSV" || ! -s "$PFI_PAPER_SCORE_REVIEW_KEY" ]]; then
	echo "INFO: blinded review inputs are missing; running 09a_make_pfi_review_inputs.sh" >&2
	bash "$PFI_ROOT/09a_make_pfi_review_inputs.sh"
fi

# Prefer human/expert labels if present. If missing, automatically create weak-supervision labels.
if [[ ! -s "$paper_score_file" ]]; then
	pseudo_file="$PFI_PAPER_SCORE_REVIEW_PSEUDO"
	if [[ ! -s "$pseudo_file" ]]; then
		echo "INFO: scored review CSV is missing; creating weak-supervision pseudo labels automatically." >&2
		bash "$PFI_ROOT/09a_auto_score_pfi_review_dataset.sh" "$PFI_PAPER_SCORE_REVIEW_CSV" "$pseudo_file"
	fi
	paper_score_file="$pseudo_file"
	echo "INFO: using weak-supervision pseudo labels for training: $paper_score_file" >&2
fi

"$PFI_PYTHON_BIN" -u f/make_pfi_training_dataset.py \
	--review_scores "$paper_score_file" \
	--review_key "$PFI_PAPER_SCORE_REVIEW_KEY" \
	--blinded "$PFI_BLINDED_TEXT" \
	--features "$PFI_ARTICLE_FEATURES" \
	--journal_rank "$rank_file" \
	--master "$PFI_ARTICLES_MASTER" \
	--out "$PFI_TRAINING_DATASET" \
	--summary_out "$DAT_DIR/training_dataset_summary.csv" \
	--fit_summary_out "$DAT_DIR/training_label_fit_summary.csv" \
	--audit_csv_out "$DAT_DIR/training_records_preview.csv" \
	--max_rows "$PFI_MAX_TRAIN_ROWS" \
	--min_text_chars "$PFI_MIN_TEXT_CHARS" \
	--seed "$PFI_REVIEW_SEED" \
	2>&1 | tee "$LOG_DIR/09a_prepare_pfi_training_dataset.log"

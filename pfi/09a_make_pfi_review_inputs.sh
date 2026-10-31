#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
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
    Default  Prepare blinded review files

Main options:
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09a_make_pfi_review_inputs.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09a_make_pfi_review_inputs
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"
: "${PFI_HIGH_ANCHOR_CSV:=}"
: "${PFI_LOW_ANCHOR_CSV:=}"
extra_args=()
[[ -n "$PFI_HIGH_ANCHOR_CSV" && -s "$PFI_HIGH_ANCHOR_CSV" ]] && extra_args+=(--high_anchor_csv "$PFI_HIGH_ANCHOR_CSV")
[[ -n "$PFI_LOW_ANCHOR_CSV" && -s "$PFI_LOW_ANCHOR_CSV" ]] && extra_args+=(--low_anchor_csv "$PFI_LOW_ANCHOR_CSV")
if [[ -z "${PFI_REVIEW_MAX_SCAN_ROWS:-}" ]]; then
	PFI_REVIEW_MAX_SCAN_ROWS=$((PFI_REVIEW_ROWS * 2000))
	if ((PFI_REVIEW_MAX_SCAN_ROWS < 200000)); then
		PFI_REVIEW_MAX_SCAN_ROWS=200000
	fi
fi
"$PFI_PYTHON_BIN" -u f/make_pfi_review_inputs.py \
	--master "$PFI_ARTICLES_MASTER" \
	--blinded "$PFI_BLINDED_TEXT" \
	--features "$PFI_ARTICLE_FEATURES" \
	--journal_rank "$PFI_JOURNAL_RANK_NORMALIZED" \
	--out_csv "$PFI_PAPER_SCORE_REVIEW_CSV" \
	--out_xlsx "$PFI_PAPER_SCORE_REVIEW_XLSX" \
	--key_out "$PFI_PAPER_SCORE_REVIEW_KEY" \
	--summary_out "$PFI_PAPER_SCORE_REVIEW_SUMMARY" \
	--n "$PFI_REVIEW_ROWS" \
	--min_text_chars "$PFI_MIN_TEXT_CHARS" \
	--max_scan_rows "$PFI_REVIEW_MAX_SCAN_ROWS" \
	--threads "${PFI_DUCKDB_THREADS:-4}" \
	--seed "$PFI_REVIEW_SEED" \
	"${extra_args[@]}" \
	2>&1 | tee "$LOG_DIR/09a_make_pfi_review_inputs.log"

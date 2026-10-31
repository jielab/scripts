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
    Default  Sample articles for review

Main options:
    N  Number of records [1000]
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./11_sample_llm_review_batch.sh 500
    ./11_sample_llm_review_batch.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 11_sample_llm_review_batch
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"
n="${1:-1000}"
"$PFI_PYTHON_BIN" -u f/sample_llm_review_batch.py \
	--blinded "$PFI_BLINDED_TEXT" \
	--features "$PFI_ARTICLE_FEATURES" \
	--out "$DAT_DIR/llm_review_sample.jsonl" \
	--n "$n" 2>&1 | tee "$LOG_DIR/11_sample_llm_review_batch.log"

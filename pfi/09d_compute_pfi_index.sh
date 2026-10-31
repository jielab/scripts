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
    Default  Compute the PFI index

Main options:
    RANK_FILE  Optional positional normalized rank file
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09d_compute_pfi_index.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09d_compute_pfi_index
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$PFI_ANALYSIS_DATA_DIR" "$LOG_DIR"
rank_file="${1:-$PFI_JOURNAL_RANK_NORMALIZED}"
institution_rank_file="${2:-}"
args=()
[[ "$PFI_RESULTS_FULL_TEXT_ONLY" == "1" ]] && args+=(--full_text_only)
[[ -n "$institution_rank_file" && -s "$institution_rank_file" ]] && args+=(--institution_rank "$institution_rank_file")
"$PFI_PYTHON_BIN" -u f/compute_pfi_index.py \
	--master "$PFI_ARTICLES_MASTER" \
	--predictions "$PFI_PUBMEDBERT_PREDICTIONS" \
	--journal_rank "$rank_file" \
	--out_csv "$PFI_FAIRNESS_RESULTS_CSV" \
	--out_parquet "$PFI_FAIRNESS_RESULTS_PARQUET" \
	--summary_dir "$PFI_ANALYSIS_DATA_DIR" \
	"${args[@]}" \
	2>&1 | tee "$LOG_DIR/09d_compute_pfi_index.log"

#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Default  Run the complete PFI workflow

Main options:
    --plan-only  Write a local/LSF .cmd without running
    --job-cores N --job-memory-gb N --job-hours N --job-gpus N --queue NAME
    Environment is detected automatically; --submit auto|local|hpc is optional.
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Examples:
    ./00_run_pfi_all.sh
PUBLIC_ENTRY_HELP
    exit 0
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

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			bash "${BASH_SOURCE[0]}" --help
			exit 0
			;;
	esac
fi


# 00_run_pfi_all
set -euo pipefail

bash 02_download_pubmed.sh
bash 03_parse_pubmed.sh
bash 04_download_pmc_oa_bulk_xml.sh
bash 05_parse_pmc_oa_bulk_xml.sh
bash 06_build_master_dataset.sh
bash 07_make_blinded_text.sh
bash 08_extract_basic_features.sh
rank_args=()
if [[ $# -gt 0 && -n "${1:-}" ]]; then
	rank_args=("$1")
fi
bash 08b_normalize_journal_rank.sh "${rank_args[@]}"
bash 09a_make_pfi_review_inputs.sh
cat <<'MSG'
STOP: Fill the blinded review XLSX, save it as paper_score_review_1000_scored.csv, then run:
  ./pfi.sh prepare-training
  ./pfi.sh train
  ./pfi.sh predict
  ./pfi.sh compute
  ./pfi.sh analysis
MSG

#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Default  Prepare blinded article text

Main options:
    --plan-only  Write a local/LSF .cmd without running
    --job-cores N --job-memory-gb N --job-hours N --job-gpus N --queue NAME
    Environment is detected automatically; --submit auto|local|hpc is optional.
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Examples:
    ./07_make_blinded_text.sh
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
set -euo pipefail


# Shell options, paths, and shared inputs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$PFI_WORK_PARQUET_DIR" "$PFI_WORK_DUCKDB_TEMP_DIR" "$LOG_DIR"

full_text_arg=()
: "${PFI_BLINDED_FULL_TEXT_ONLY:=1}"
if [[ "$PFI_BLINDED_FULL_TEXT_ONLY" == "1" ]]; then
	full_text_arg=(--full_text_only)
fi
memory_limit_arg=()
if [[ -n "${PFI_DUCKDB_MEMORY_LIMIT:-}" ]]; then
	memory_limit_arg=(--memory_limit "$PFI_DUCKDB_MEMORY_LIMIT")
fi
resume_arg=()
: "${PFI_BLINDED_RESUME:=1}"
if [[ "$PFI_BLINDED_RESUME" == "1" ]]; then
	resume_arg=(--resume)
fi

{
	echo "[$(date -Is)] 07_make_blinded_text.sh started"
	"$PFI_PYTHON_BIN" -u f/make_blinded_text.py \
		--master "$PFI_ARTICLES_MASTER" \
		--out "$PFI_BLINDED_TEXT" \
		--temp_dir "$PFI_WORK_DUCKDB_TEMP_DIR" \
		--threads "${PFI_DUCKDB_THREADS:-4}" \
		"${memory_limit_arg[@]}" \
		--engine "${PFI_BLINDED_ENGINE:-pyarrow}" \
		--batch_size "${PFI_BLINDED_BATCH_SIZE:-1000}" \
		--rows_per_file "${PFI_BLINDED_ROWS_PER_FILE:-100000}" \
		--progress_every_rows "${PFI_BLINDED_PROGRESS_EVERY_ROWS:-250000}" \
		"${resume_arg[@]}" \
		"${full_text_arg[@]}" \
		--partition_by_year
	echo "[$(date -Is)] 07_make_blinded_text.sh finished"
} 2>&1 | tee "$LOG_DIR/07_make_blinded_text.log"

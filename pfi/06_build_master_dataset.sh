#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Default  Build the article master dataset

Main options:
    --plan-only  Write a local/LSF .cmd without running
    --job-cores N --job-memory-gb N --job-hours N --job-gpus N --queue NAME
    Environment is detected automatically; --submit auto|local|hpc is optional.
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Examples:
    ./06_build_master_dataset.sh
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
mkdir -p "$PFI_WORK_DB_DIR" "$PFI_WORK_DUCKDB_TEMP_DIR" "$PFI_WORK_PARQUET_DIR" "$LOG_DIR"

build_jobs="${PFI_BUILD_JOBS:-2}"
build_threads="${PFI_BUILD_DUCKDB_THREADS:-${PFI_DUCKDB_THREADS:-2}}"
duckdb_memory_limit="${PFI_DUCKDB_MEMORY_LIMIT:-8GB}"

echo "06 build settings:"
echo "  build jobs: $build_jobs"
echo "  DuckDB threads per build job: $build_threads"
echo "  DuckDB memory limit per connection: $duckdb_memory_limit"
echo "  reuse PMC partitions: ${PFI_REUSE_PMC:-1}"

args=(
	--pubmed_glob "$PARQUET_DIR/pubmed/**/*.parquet"
	--pmc_glob "$PARQUET_DIR/pmc_oa/articles/*.parquet"
	--work_dir "$PFI_WORK_PARQUET_DIR/06_build"
	--temp_dir "$PFI_WORK_DUCKDB_TEMP_DIR"
	--out "$PFI_ARTICLES_MASTER"
	--jobs "$build_jobs"
	--threads "$build_threads"
	--buckets "${PFI_BUILD_BUCKETS:-128}"
	--memory_limit "$duckdb_memory_limit"
	--progress_every_buckets "${PFI_BUILD_PROGRESS_EVERY_BUCKETS:-8}"
	--fresh_work
	--reuse_pubmed
)

if [[ "${PFI_REUSE_PMC:-1}" != "0" ]]; then
	args+=(--reuse_pmc)
fi

"$PFI_PYTHON_BIN" -u f/build_master_dataset.py "${args[@]}" 2>&1 | tee "$LOG_DIR/06_build_master_dataset.log"

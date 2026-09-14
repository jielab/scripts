#!/usr/bin/env bash
set -euo pipefail


# 🚩 Shell options, paths, and shared inputs
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
  --pubmed_glob "$PARQUET_DIR/pubmed/**/*.parquet" \
  --pmc_glob "$PARQUET_DIR/pmc_oa/articles/*.parquet" \
  --work_dir "$PFI_WORK_PARQUET_DIR/06_build" \
  --temp_dir "$PFI_WORK_DUCKDB_TEMP_DIR" \
  --out "$PFI_ARTICLES_MASTER" \
  --jobs "$build_jobs" \
  --threads "$build_threads" \
  --buckets "${PFI_BUILD_BUCKETS:-128}" \
  --memory_limit "$duckdb_memory_limit" \
  --progress_every_buckets "${PFI_BUILD_PROGRESS_EVERY_BUCKETS:-8}" \
  --fresh_work \
  --reuse_pubmed
)

if [[ "${PFI_REUSE_PMC:-1}" != "0" ]]; then
  args+=(--reuse_pmc)
fi

"$PFI_PYTHON_BIN" -u f/build_master_dataset.py "${args[@]}" 2>&1 | tee "$LOG_DIR/06_build_master_dataset.log"

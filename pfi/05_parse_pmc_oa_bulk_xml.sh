#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Default  Parse downloaded PMC XML

Main options:
    --plan-only  Write a local/LSF .cmd without running
    --job-cores N --job-memory-gb N --job-hours N --job-gpus N --queue NAME
    Environment is detected automatically; --submit auto|local|hpc is optional.
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Examples:
    ./05_parse_pmc_oa_bulk_xml.sh
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
mkdir -p "$PARQUET_DIR/pmc_oa/articles" "$LOG_DIR"

LOG_FILE="$LOG_DIR/05_parse_pmc_oa_bulk_xml.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== 05_parse_pmc_oa_bulk_xml.sh started at $(date -Is) ====="
echo "Logging to: $LOG_FILE"

max_files=0
if [[ "$PFI_TEST_MODE" == "1" ]]; then
	max_files=2
	echo "Test mode: parsing only the first 2 PMC OA tarballs."
fi

workers="${PMC_PARSE_WORKERS:-$N_JOBS}"
batch_size="${PMC_PARSE_BATCH_SIZE:-500}"
echo "PMC parse settings:"
echo "  workers: $workers"
echo "  batch size: $batch_size"
echo "  max chars per article: $FULL_TEXT_MAX_CHARS"

"$PFI_PYTHON_BIN" - "$PARQUET_DIR/pmc_oa/articles" <<'PY'
import re
import sys
from pathlib import Path

out_dir = Path(sys.argv[1])
removed = []

for path in out_dir.glob("*.tmp"):
    path.unlink()
    removed.append(path.name)

for path in out_dir.glob("*.parquet"):
    if path.stat().st_size == 0:
        path.unlink()
        removed.append(path.name)

part_re = re.compile(r"^(?P<prefix>.+)\.part\d+\.parquet$")
prefixes = {}
for path in out_dir.glob("*.part*.parquet"):
    m = part_re.match(path.name)
    if m:
        prefixes.setdefault(m.group("prefix"), []).append(path)

for prefix, paths in prefixes.items():
    done_file = out_dir / f"{prefix}.done"
    if done_file.exists():
        continue
    for path in paths:
        path.unlink()
        removed.append(path.name)

if removed:
    print(f"Cleaned {len(removed)} incomplete parse output files.")
    for name in removed[:20]:
        print(f"  removed: {name}")
else:
    print("No incomplete parse output files to clean.")
PY

"$PFI_PYTHON_BIN" f/parse_pmc_oa_bulk.py \
	--in_dir "$PMC_OA_RAW_DIR" \
	--out_dir "$PARQUET_DIR/pmc_oa/articles" \
	--min_year "$MIN_YEAR" \
	--max_chars "$FULL_TEXT_MAX_CHARS" \
	--workers "$workers" \
	--batch_size "$batch_size" \
	--max_files "$max_files"

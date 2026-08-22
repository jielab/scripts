#!/usr/bin/env bash
set -euo pipefail




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  echo "PFI_TEST_MODE=1: parsing only first 2 PMC OA tarballs."
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

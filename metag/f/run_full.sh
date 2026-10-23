#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# No disk-wide loops: one explicit dataset, audited database, frozen virus map.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT=/mnt/f/metag/ocean; DB=""; VIRUS_MAP=""; OUT=""; THREADS=16; EXECUTE=0; CHECK_MD5=0
usage() { cat <<'EOF'
Modules:
    Default  Audit -> check manifests -> classify -> summarize

Main options:
    --dir-db DIR --virus-map FILE  Prepared reference inputs
    --dir-root DIR --dir-out DIR  Single-dataset input/output
    --threads N --check-md5  Resources and optional MD5 recheck
    --execute  Execute; otherwise print plan only
    -h, --help  Help; prerequisites: 0.README.md

Example:
    # Single-directory form; replace reference placeholders with prepared files.
    ./metag.sh full --db /path/to/database --virus-map /path/to/virus-map.tsv --root /mnt/f/metag/ocean --out /mnt/f/metag/ocean/results_v2 --threads 16
    # Add --execute after checking the plan; --check-md5 requests a full MD5 recheck.

    ./metag.sh full emp
    ./metag.sh full emp,Weurope --execute
EOF
}
while (($#)); do
  case "$1" in
    --root) ROOT="$2"; shift 2;; --db) DB="$2"; shift 2;; --virus-map) VIRUS_MAP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;; --threads) THREADS="$2"; shift 2;; --execute) EXECUTE=1; shift;;
    --check-md5) CHECK_MD5=1; shift;;
    -h|--help) usage; exit 0;; *) usage >&2; exit 2;;
  esac
done
[[ -n "$DB" && -n "$VIRUS_MAP" && "$THREADS" =~ ^[1-9][0-9]*$ ]] || { usage >&2; exit 2; }
ROOT="$(realpath -m -- "$ROOT")"; DB="$(realpath -m -- "$DB")"; VIRUS_MAP="$(realpath -m -- "$VIRUS_MAP")"
OUT="${OUT:-$ROOT/results_v2}"
step() {
  printf '%q ' "$@"; printf '\n'
  if ((EXECUTE)); then "$@"; fi
}
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" db-audit --db "$DB" --virus-map "$VIRUS_MAP" --threads "$THREADS"
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" run --root "$ROOT" --db "$DB" --virus-map "$VIRUS_MAP" --out "$OUT" --threads "$THREADS" --dry-run
if ((CHECK_MD5)); then
  step "${METAG_PYTHON:-python3}" "$HERE/metag.py" verify --root "$ROOT" --md5
else
  echo '跳过原始 FASTQ 的重复 MD5 检查（由用户在流程外完成）；需要重查时加 --check-md5。'
fi
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" run --root "$ROOT" --db "$DB" --virus-map "$VIRUS_MAP" --out "$OUT" --threads "$THREADS"
step "${METAG_PYTHON:-python3}" "$HERE/metag.py" summary --out "$OUT"

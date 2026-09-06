#!/usr/bin/env bash
# One-time migration of existing gwas_post cis outputs to coordinate-sorted
# BGZF plus tabix indexes. New cis outputs are indexed directly by gwas_post.sh.

set -euo pipefail
export LC_ALL=C

DATA_ROOT="${CIS_INDEX_DATA_ROOT:-/mnt/d/data.BIG/gwas}"
PROJECTS="${CIS_INDEX_PROJECTS:-main,prot,met}"
JOBS="${CIS_INDEX_JOBS:-6}"
THREADS="${CIS_INDEX_THREADS:-2}"
SORT_MEMORY="${CIS_INDEX_SORT_MEMORY:-512M}"
TMP_DIR="${CIS_INDEX_TMP_DIR:-/mnt/d/tmp}"
INDEX_F="${CIS_INDEX_FUNCTIONS:-/mnt/d/scripts/0f/gwas_index.f.sh}"
REPLACE=FALSE
DRY_RUN=FALSE

usage() {
  cat <<'EOF'
Usage: ./index_cis.sh [options]

Convert existing canonical cis files
  <data-root>/<project>/common/<trait>/gwas/<trait>.cis.gz
to coordinate-sorted BGZF and create <trait>.cis.gz.tbi (or .csi when needed).

Options:
  --data-root PATH       GWAS parent directory [/mnt/d/data.BIG/gwas]
  --projects CSV         Projects to scan [main,prot,met]
  --jobs N               Files processed concurrently [6]
  --threads N            bgzip threads per file [2]
  --sort-memory SIZE     GNU sort memory per file if sorting is needed [512M]
  --tmp-dir PATH         GNU sort temporary directory [/mnt/d/tmp]
  --index-f FILE         shared BGZF/tabix helper [/mnt/d/scripts/0f/gwas_index.f.sh]
  --replace TRUE|FALSE   rebuild even a valid existing index [FALSE]
  --dry-run              list exact targets without changing them
  -h, --help             show this help

The original .cis.gz is replaced only after a temporary BGZF file and its
index pass validation. Existing valid, fresh .tbi/.csi files are skipped.
EOF
}

bool_word() {
  case "${1,,}" in
    true|1|yes|y) printf 'TRUE\n' ;;
    false|0|no|n) printf 'FALSE\n' ;;
    *) return 2 ;;
  esac
}

need_value() {
  [[ -n "${2-}" && "${2-}" != --* ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root) need_value "$1" "${2-}"; DATA_ROOT="$2"; shift 2 ;;
    --projects) need_value "$1" "${2-}"; PROJECTS="$2"; shift 2 ;;
    --jobs) need_value "$1" "${2-}"; JOBS="$2"; shift 2 ;;
    --threads) need_value "$1" "${2-}"; THREADS="$2"; shift 2 ;;
    --sort-memory) need_value "$1" "${2-}"; SORT_MEMORY="$2"; shift 2 ;;
    --tmp-dir) need_value "$1" "${2-}"; TMP_DIR="$2"; shift 2 ;;
    --index-f) need_value "$1" "${2-}"; INDEX_F="$2"; shift 2 ;;
    --replace) need_value "$1" "${2-}"; REPLACE=$(bool_word "$2") || { echo "ERROR: --replace expects TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
    --dry-run) DRY_RUN=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs must be a positive integer" >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be a positive integer" >&2; exit 2; }
[[ -d "$DATA_ROOT" ]] || { echo "ERROR: missing data root: $DATA_ROOT" >&2; exit 1; }
[[ -d "$TMP_DIR" ]] || { echo "ERROR: missing sort temp directory: $TMP_DIR" >&2; exit 1; }
[[ -s "$INDEX_F" ]] || { echo "ERROR: missing shared index helper: $INDEX_F" >&2; exit 1; }
command -v gzip >/dev/null || { echo "ERROR: gzip not found" >&2; exit 1; }
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found (install htslib)" >&2; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found (install htslib)" >&2; exit 1; }
command -v sort >/dev/null || { echo "ERROR: GNU sort not found" >&2; exit 1; }

# shellcheck source=/mnt/d/scripts/0f/gwas_index.f.sh
source "$INDEX_F"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

index_cis_one() {
  local file="$1"
  [[ -s "$file" ]] || { log "ERROR missing/empty: $file"; return 1; }
  if [[ "$REPLACE" != TRUE ]] && gwas_index_valid "$file"; then
    log "SKIP valid index: $file"
    return 0
  fi
  [[ "$REPLACE" != TRUE ]] || rm -f -- "${file}.tbi" "${file}.csi"
  log "INDEX cis: $file"
  gwas_index_file "$file"
  bgzip -t "$file" >/dev/null 2>&1 || { log "ERROR BGZF validation failed: $file"; return 1; }
  gwas_index_valid "$file" || { log "ERROR tabix validation failed: $file"; return 1; }
  log "DONE BGZF+index: $file"
}

target_list=$(mktemp "${TMPDIR:-/tmp}/index_cis.targets.XXXXXX")
trap 'rm -f -- "$target_list"' EXIT

IFS=',' read -r -a project_array <<< "$PROJECTS"
for project in "${project_array[@]}"; do
  project=${project//[[:space:]]/}
  [[ -n "$project" && "$project" != "." && "$project" != ".." && "$project" != */* ]] || {
    echo "ERROR: invalid project name: $project" >&2; exit 2;
  }
  root="$DATA_ROOT/$project/common"
  [[ -d "$root" ]] || { log "WARNING missing project root, skip: $root"; continue; }
  while IFS= read -r file; do
    trait=$(basename "$(dirname "$(dirname "$file")")")
    [[ "$(basename "$file")" == "${trait}.cis.gz" ]] && printf '%s\n' "$file"
  done < <(find "$root" -mindepth 3 -maxdepth 3 -type f -path '*/gwas/*.cis.gz' -size +0c | sort -V)
done | sort -u -V > "$target_list"

target_count=$(wc -l < "$target_list" | tr -d ' ')
log "Targets=$target_count projects=$PROJECTS jobs=$JOBS bgzip_threads=$THREADS"
if (( target_count == 0 )); then
  log "No existing canonical cis files found; nothing to do."
  exit 0
fi

if [[ "$DRY_RUN" == TRUE ]]; then
  cat "$target_list"
  exit 0
fi

export -f log gwas_index_log gwas_index_header gwas_index_col gwas_index_valid
export -f gwas_index_make_sidecar gwas_index_file index_cis_one
export REPLACE GWAS_CLEAN_COMP_THREADS="$THREADS"
export GWAS_INDEX_SORT_MEMORY="$SORT_MEMORY" GWAS_INDEX_TMP_DIR="$TMP_DIR"
export SHELL=/bin/bash
if command -v parallel >/dev/null 2>&1; then
  parallel --line-buffer --halt soon,fail=1 -j "$JOBS" index_cis_one :::: "$target_list"
else
  while IFS= read -r file; do index_cis_one "$file"; done < "$target_list"
fi

log "All existing canonical cis files are BGZF-indexed."

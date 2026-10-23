#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# Broad REFSEQ database, not a claim to cover every organism.
# Downloads can be very large. WITHOUT --execute this only prints the plan.
set -euo pipefail
DB=""; THREADS=16; DOWNLOAD_THREADS=8; MAX_GB=24; EXECUTE=0; VIRUS_MAP=""
MASK_JOBS=2; MASK_THREADS=1; PREP_MEMORY_GB=12
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
K2=("${METAG_PYTHON:-python3}" "$HERE/k2_compat.py")
usage() {
	cat <<'EOF'
Modules:
    Default  Build a new reference index or audit an existing complete index

Main options:
    --dir-db DIR --virus-map FILE  Index and frozen virus mapping
    --threads N  Workers [16]
    --mask-jobs N --masker-threads N  Parallel genomes / threads per masker [2 / 1]
    --prepare-memory-gb N  Preparation RAM ceiling, separate from index build [12]
    --max-db-gb N  Hash-table cap in GiB [24]; 0 disables cap
    --execute  Execute; otherwise print plan only
    -h, --help  Help; reference details: 0.README.md

Example:
    # Build into a NEW directory; omit --execute to inspect the plan first.
    bash metag.sh build-db --db /path/to/NEW/database --threads 16 --max-db-gb 48 --execute

    ./metag.sh build-db --dir-db /path/to/new_index
    ./metag.sh build-db --dir-db /path/to/new_index --max-db-gb 24 --execute
EOF
}
while (($#)); do
  case "$1" in
    --db) DB="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --download-threads) DOWNLOAD_THREADS="$2"; shift 2;;
    --mask-jobs) MASK_JOBS="$2"; shift 2;;
    --masker-threads) MASK_THREADS="$2"; shift 2;;
    --prepare-memory-gb) PREP_MEMORY_GB="$2"; shift 2;;
    --max-db-gb) MAX_GB="$2"; shift 2;;
    --virus-map) VIRUS_MAP="$2"; shift 2;;
    --execute) EXECUTE=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done
[[ -n "$DB" ]] || { usage; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$DOWNLOAD_THREADS" =~ ^[1-9][0-9]*$ && "$MAX_GB" =~ ^[0-9]+$ ]] || exit 2
[[ "$MASK_JOBS" =~ ^[1-9][0-9]*$ && "$MASK_THREADS" =~ ^[1-9][0-9]*$ && "$PREP_MEMORY_GB" =~ ^[1-9][0-9]*$ ]] || exit 2
((MASK_JOBS * MASK_THREADS <= THREADS)) || { echo 'Masking workers exceed --threads budget' >&2; exit 2; }
DB="$(realpath -m -- "$DB")"
AUDIT=("${METAG_PYTHON:-python3}" "$HERE/metag.py" db-audit --db "$DB" --threads "$THREADS")
if [[ -n "$VIRUS_MAP" ]]; then AUDIT+=(--virus-map "$(realpath -m -- "$VIRUS_MAP")"); else AUDIT+=(--base-only); fi
# Cleanup only follows full coverage auditing of a completed managed database.
cleanup_reference() {
  if [[ -n "$VIRUS_MAP" ]]; then
    "${METAG_PYTHON:-python3}" "$HERE/reference_cleanup.py" --db "$DB" --virus-map "$VIRUS_MAP" --execute
  fi
}
# Never mutate a pre-existing index, even one built outside this wrapper.
if [[ -f "$DB/hash.k2d" && -f "$DB/opts.k2d" && -f "$DB/taxo.k2d" ]]; then
  if ((EXECUTE)); then
    "${AUDIT[@]}"
    cleanup_reference
  else
    printf '%q ' "${AUDIT[@]}"; printf '\n'
    if [[ -n "$VIRUS_MAP" ]]; then
      printf '%q ' "${METAG_PYTHON:-python3}" "$HERE/reference_cleanup.py" --db "$DB" --virus-map "$VIRUS_MAP" --execute
      printf '\n'
    fi
  fi
  exit 0
fi
if [[ -f "$DB/hash.k2d" || -f "$DB/opts.k2d" || -f "$DB/taxo.k2d" ]]; then
  echo "Partial index exists; inspect it and use a NEW database path." >&2; exit 1
fi
LIBRARIES=(archaea bacteria viral human fungi plant protozoa vertebrate_mammalian vertebrate_other invertebrate UniVec_Core)
BUILD=("${METAG_PYTHON:-python3}" -B "$HERE/index_inputs.py" build -- --build --db "$DB" --threads "$THREADS")
if ((MAX_GB > 0)); then BUILD+=(--max-db-size "$((MAX_GB * 1024 * 1024 * 1024))"); fi
print_cmd() { printf '%q ' "$@"; printf '\n'; }
if ((EXECUTE == 0)); then
  echo '# PLAN ONLY; no downloads or database writes have occurred.'
  print_cmd "${K2[@]}" download-taxonomy --db "$DB" --skip-maps
  for lib in "${LIBRARIES[@]}"; do
    print_cmd "${METAG_PYTHON:-python3}" "$HERE/reference_limits.py" --memory-gb "$PREP_MEMORY_GB" -- \
      env METAG_BATCHED_REFERENCE=1 "METAG_MASK_JOBS=$MASK_JOBS" "${K2[@]}" download-library \
      --db "$DB" --library "$lib" --threads "$DOWNLOAD_THREADS" --masker-threads "$MASK_THREADS" \
      --assembly-source refseq --assembly-levels chromosome complete_genome --resume
  done
  print_cmd "${BUILD[@]}"
  print_cmd "${AUDIT[@]}"
  if [[ -n "$VIRUS_MAP" ]]; then
    print_cmd "${METAG_PYTHON:-python3}" "$HERE/reference_cleanup.py" --db "$DB" --virus-map "$VIRUS_MAP" --execute
  fi
  echo '# Execute with --execute only after checking available disk and RAM.'
  exit 0
fi
for tool in k2 kraken2-build; do
  command -v "$tool" >/dev/null || { echo 'Install the supplied environment with Kraken2 2.17.1.' >&2; exit 1; }
done
mkdir -p "$DB/metag_build_logs" "$DB/metag_build_state"
exec 9>"$DB/.metag_build.lock"
BUILD_PID_FILE="$DB/metag_build_state/owner.pid"
if ! flock -n 9; then
  echo 'Database build lock is held by an existing task (possibly suspended with Ctrl+Z).' >&2
  if [[ -r "$BUILD_PID_FILE" ]]; then
    owner_pid=$(cat "$BUILD_PID_FILE")
    if [[ "$owner_pid" =~ ^[1-9][0-9]*$ ]]; then
      ps -p "$owner_pid" -o pid,stat,etime,args >&2 || true
    fi
  fi
  echo 'In the original terminal, check jobs -l; use fg to resume a suspended job. Do not delete the lock file.' >&2
  exit 1
fi
printf '%s\n' "$$" > "$BUILD_PID_FILE"
trap 'rm -f -- "$BUILD_PID_FILE"' EXIT
if [[ -f "$DB/metag_build_state/all.complete" ]]; then
  echo "A completed frozen database already exists at $DB. Use a NEW path to change it." >&2
  exit 1
fi
# Record the configuration and refuse to resume a partial build with different settings.
CONFIG=$(printf 'libraries=%s\nmax_db_gb=%s\n' "${LIBRARIES[*]}" "$MAX_GB")
SELECTION_CONFIG=$("${METAG_PYTHON:-python3}" "$HERE/reference_selection.py" signature --db "$DB")
if [[ -n "$SELECTION_CONFIG" ]]; then CONFIG+=$'\n'"$SELECTION_CONFIG"; fi
if [[ -f "$DB/metag_build_state/config.txt" ]] && [[ "$(cat "$DB/metag_build_state/config.txt")" != "$CONFIG" ]]; then
  echo 'Build configuration changed; use a NEW database directory.' >&2; exit 1
fi
printf '%s\n' "$CONFIG" > "$DB/metag_build_state/config.txt"
kraken2 --version > "$DB/metag_build_logs/kraken2.version.txt" 2>&1
valid_step() {
  case "$1" in
    taxonomy) [[ -s "$DB/taxonomy/nodes.dmp" && -s "$DB/taxonomy/names.dmp" ]];;
    build) [[ -s "$DB/hash.k2d" && -s "$DB/opts.k2d" && -s "$DB/taxo.k2d" ]];;
    audit) [[ -s "$DB/metag.audit.json" ]];;
    *) "${METAG_PYTHON:-python3}" -B "$HERE/index_inputs.py" check-library --dir-db "$DB" --library "$1";;
  esac
}
step() {
  local name="$1"; shift
  if [[ -f "$DB/metag_build_state/$name.complete" ]]; then
    if valid_step "$name"; then
      echo "SKIP completed step: $name"; return
    fi
    echo "Retrying step with missing/empty products: $name"
    rm -- "$DB/metag_build_state/$name.complete"
  fi
  print_cmd "$@" | tee -a "$DB/metag_build_logs/commands.txt"
  "$@" 2>&1 | tee "$DB/metag_build_logs/$name.log"
  valid_step "$name" || { echo "Step $name produced missing/empty files; not marked complete." >&2; return 1; }
  touch "$DB/metag_build_state/$name.complete"
}
step taxonomy "${K2[@]}" download-taxonomy --db "$DB" --skip-maps
RESUME=()
if "${K2[@]}" download-library --help 2>&1 | grep -q -- '--resume'; then RESUME=(--resume); fi
for lib in "${LIBRARIES[@]}"; do
  step "$lib" "${METAG_PYTHON:-python3}" "$HERE/reference_limits.py" --memory-gb "$PREP_MEMORY_GB" -- \
    env METAG_BATCHED_REFERENCE=1 "METAG_MASK_JOBS=$MASK_JOBS" "${K2[@]}" download-library \
    --db "$DB" --library "$lib" --threads "$DOWNLOAD_THREADS" --masker-threads "$MASK_THREADS" \
    --assembly-source refseq --assembly-levels chromosome complete_genome "${RESUME[@]}"
done
# The private build wrapper does not patch the shared Conda installation.
step build "${BUILD[@]}"
step audit "${AUDIT[@]}"
"${METAG_PYTHON:-python3}" - "$DB" <<'PY'
import json, sys
from pathlib import Path
p = Path(sys.argv[1])
a = json.loads((p / 'metag.audit.json').read_text())
if a['missing_groups']:
    raise SystemExit('Missing groups after database build: ' + ','.join(a['missing_groups']))
PY
touch "$DB/metag_build_state/all.complete"
if [[ -z "$VIRUS_MAP" ]]; then echo 'Virus types NOT AUDITED: build a frozen mapping, then db-audit --virus-map before run.'; fi
echo 'Database frozen. Preserve indexes, taxonomy/*.dmp, manifests, mapping and metag.* metadata.'
echo 'Do NOT run k2 clean: it can remove the taxonomy required by this pipeline.'
# Release the build lock before the guarded cleanup acquires its own locks.
rm -f -- "$BUILD_PID_FILE"
trap - EXIT
flock -u 9
exec 9>&-
cleanup_reference

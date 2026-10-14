#!/usr/bin/env bash
# Broad REFSEQ database, not a claim to cover every organism.
# Downloads can be very large. WITHOUT --execute this only prints the plan.
set -euo pipefail
DB=""; THREADS=16; DOWNLOAD_THREADS=4; MAX_GB=48; EXECUTE=0; VIRUS_MAP=""
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
K2=("${METAG_PYTHON:-python3}" "$HERE/k2_compat.py")
usage() {
  cat <<'EOF'
bash metag.sh build-db --db /path/to/NEW/database [--threads 16] [--max-db-gb 48] [--execute]
  --max-db-gb 0 : no hash-table size cap (RAM requirements may be very large)
  --virus-map  : frozen virus_genome_types.sqlite (for existing/new compatible taxonomy)
  Existing complete index: audit only, never download into or rebuild it.
  --execute    : actually download/build; otherwise print the plan only
Uses k2 >=2.17.1 for vertebrate/invertebrate libraries. Do not substitute the older
kraken2-build --download-library interface for these additional libraries.
Index construction uses kraken2-build --build; k2 2.17.1 mishandles a size cap
above the estimated requirement and can return success after a build failure.
Reference genomes/downloads/temp files are NOT limited by --max-db-gb.
EOF
}
while (($#)); do
  case "$1" in
    --db) DB="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --download-threads) DOWNLOAD_THREADS="$2"; shift 2;;
    --max-db-gb) MAX_GB="$2"; shift 2;;
    --virus-map) VIRUS_MAP="$2"; shift 2;;
    --execute) EXECUTE=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done
[[ -n "$DB" ]] || { usage; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$DOWNLOAD_THREADS" =~ ^[1-9][0-9]*$ && "$MAX_GB" =~ ^[0-9]+$ ]] || exit 2
DB="$(realpath -m -- "$DB")"
AUDIT=("${METAG_PYTHON:-python3}" "$HERE/metag.py" db-audit --db "$DB" --threads "$THREADS")
if [[ -n "$VIRUS_MAP" ]]; then AUDIT+=(--virus-map "$(realpath -m -- "$VIRUS_MAP")"); else AUDIT+=(--base-only); fi
# Never mutate a pre-existing index, even one built outside this wrapper.
if [[ -f "$DB/hash.k2d" && -f "$DB/opts.k2d" && -f "$DB/taxo.k2d" ]]; then
  if ((EXECUTE)); then "${AUDIT[@]}"; else printf '%q ' "${AUDIT[@]}"; printf '\n'; fi
  exit 0
fi
if [[ -f "$DB/hash.k2d" || -f "$DB/opts.k2d" || -f "$DB/taxo.k2d" ]]; then
  echo "Partial index exists; inspect it and use a NEW database path." >&2; exit 1
fi
LIBRARIES=(archaea bacteria viral human fungi plant protozoa vertebrate_mammalian vertebrate_other invertebrate UniVec_Core)
BUILD=(kraken2-build --build --db "$DB" --threads "$THREADS")
if ((MAX_GB > 0)); then BUILD+=(--max-db-size "$((MAX_GB * 1024 * 1024 * 1024))"); fi
print_cmd() { printf '%q ' "$@"; printf '\n'; }
if ((EXECUTE == 0)); then
  echo '# PLAN ONLY; no downloads or database writes have occurred.'
  print_cmd "${K2[@]}" download-taxonomy --db "$DB" --skip-maps
  for lib in "${LIBRARIES[@]}"; do
    print_cmd "${K2[@]}" download-library --db "$DB" --library "$lib" --threads "$DOWNLOAD_THREADS" --assembly-source refseq --assembly-levels chromosome complete_genome
  done
  print_cmd "${BUILD[@]}"
  print_cmd "${AUDIT[@]}"
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
    *) [[ -s "$DB/library/$1/library.fna" && -s "$DB/library/$1/prelim_map.txt" ]];;
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
  step "$lib" "${K2[@]}" download-library --db "$DB" --library "$lib" --threads "$DOWNLOAD_THREADS" \
    --assembly-source refseq --assembly-levels chromosome complete_genome "${RESUME[@]}"
done
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
echo 'Database frozen. Preserve taxonomy/*.dmp, source libraries/manifests and metag.* metadata.'
echo 'Do NOT run k2 clean: it can remove the taxonomy required by this pipeline.'

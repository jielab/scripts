#!/usr/bin/env bash
set -euo pipefail

# Build the compact Sylph species-representative databases from references that
# are ALREADY present in the metag project. No raw reference is copied and the
# existing ~4.6-TB Kraken library/index is NOT rebuilt or modified.
#
# Default HPC paths resolved from metag/config.json + 0f/0os.sh:
#   representatives: /work/sph-huangj/analysis/metag/qc/reference-audit/species-representatives.tsv
#   Sylph DBs:       /work/sph-huangj/data/annot/metag/sylph/
#   sylph binary:    /work/sph-huangj/software/metag_parallel/sylph/1.0.0/sylph
#
# Usage:
#   cd /work/sph-huangj/scripts/metag
#   ./1.build-lib.sh sylph
#
# Optional:
#   SYLPH_BUCKETS=prok ./1.build-lib.sh sylph       # highest priority
#   SYLPH_BUCKETS=prok,virus,euk ./1.build-lib.sh sylph
#   SYLPH_QUEUE=38 SYLPH_THREADS=32 SYLPH_MEM_GIB=256 ./1.build-lib.sh sylph

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
for arg in "$@"; do
    if [[ $arg == -h || $arg == --help ]]; then bash "$HERE/f/usage.sh" "build_sylph"; exit 0; fi
done
CONFIG="${METAG_CONFIG:-$HERE/config.json}"
DO_SUBMIT=TRUE
ROOT_OVERRIDE=""
while (($#)); do
    case $1 in
        --config) CONFIG=$2; shift 2 ;;
        --root) ROOT_OVERRIDE=$2; shift 2 ;;
        --submit)
            [[ ${2:-} == auto || ${2:-} == local || ${2:-} == hpc ]] || { echo "--submit must be hpc, local or auto" >&2; exit 2; }
            export submit=$2 SUBMIT=$2; shift 2 ;;
        --plan-only) DO_SUBMIT=FALSE; shift ;;
        --buckets) export SYLPH_BUCKETS=$2; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
source "$HERE/f/env.sh"
if [[ $DO_SUBMIT == FALSE ]]; then
    "$METAG_PYTHON" "$HERE/f/plan_small_build.py" --module sylph --config "$CONFIG" --root "$ROOT_OVERRIDE" --submit "$submit"
    exit 0
fi

RESOLVED="$("$METAG_PYTHON" - "$HERE" "$CONFIG" "$ROOT_OVERRIDE" <<'PY'
import shlex, sys
sys.path.insert(0, sys.argv[1] + '/f')
from common import load_config
cfg = load_config(sys.argv[2], sys.argv[3] or None)
for key in ('library_dir','work_dir','tools_dir','reference_raw','reads_raw'):
    print(f'{key.upper()}={shlex.quote(str(cfg[key]))}')
PY
)"
eval "$RESOLVED"

ANNOT_METAG="$(dirname -- "$LIBRARY_DIR")"
AUDIT_DIR="${AUDIT_DIR:-$WORK_DIR/qc/reference-audit}"
REPRESENTATIVES="$AUDIT_DIR/species-representatives.tsv"
SYLPH_DIR="${SYLPH_DIR:-$ANNOT_METAG/sylph}"
SYLPH_TOOL_DIR="${SYLPH_TOOL_DIR:-$TOOLS_DIR/sylph/1.0.0}"
SYLPH="$SYLPH_TOOL_DIR/sylph"
LOG_DIR="$WORK_DIR/logs/small-db/sylph"
JOB_DIR="$WORK_DIR/jobs/small-db/sylph"
mkdir -p "$AUDIT_DIR" "$SYLPH_DIR" "$SYLPH_TOOL_DIR" "$LOG_DIR" "$JOB_DIR"
exec 9>"$SYLPH_DIR/.submission.lock"
flock -n 9 || { echo 'Another Sylph DB preparation is active' >&2; exit 2; }

SYLPH_VERSION=1.0.0
# Upstream publishes its static executable under the `latest` asset, not under
# the v1.0.0 release (the latter URL is HTTP 404). Validate the actual version.
SYLPH_URL="${SYLPH_URL:-https://github.com/bluenote-1577/sylph/releases/download/latest/sylph}"

fetch() {
    local url=$1 out=$2
    if [[ -s "$out" ]]; then return; fi
    local tmp="${out}.partial"
    if command -v aria2c >/dev/null 2>&1; then
        aria2c -c -x 8 -s 8 --file-allocation=none -d "$(dirname "$out")" -o "$(basename "$tmp")" "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl -fLsS --connect-timeout 30 --retry 6 --retry-all-errors --retry-delay 5 -C - -o "$tmp" "$url"
    elif command -v wget >/dev/null 2>&1; then
        wget -c -O "$tmp" "$url"
    else
        echo "ERROR: need aria2c, curl, or wget" >&2; exit 2
    fi
    [[ -s "$tmp" ]] || { echo "ERROR: empty download: $tmp" >&2; exit 2; }
    mv -f "$tmp" "$out"
}

# 1) Install pinned sylph v1.0.0 static binary if needed.
if [[ ! -x "$SYLPH" ]]; then
    echo "== Install sylph v${SYLPH_VERSION} =="
    BUNDLED="${DIR_SOFTWARE}/metag/sylph-1.0.0-linux-x86_64"
    if [[ -s "$BUNDLED" ]]; then
        cp "$BUNDLED" "$SYLPH"
    else
        fetch "$SYLPH_URL" "$SYLPH"
    fi
    chmod +x "$SYLPH"
fi
actual_version=$("$SYLPH" --version)
[[ $actual_version == "sylph $SYLPH_VERSION" ]] || {
    echo "Expected sylph $SYLPH_VERSION, got: $actual_version. Set SYLPH_TOOL_DIR to an existing 1.0.0 installation or SYLPH_URL to its binary." >&2
    exit 2
}
"$SYLPH" -h >/dev/null
sha256sum "$SYLPH" > "$SYLPH_TOOL_DIR/sylph.sha256"

# 2) Metadata-only audit of the EXISTING library. This scans manifests/receipts
# and NCBI summary metadata; it does not reread/rebuild the TB-scale hash files.
echo "== Select one representative assembly per NCBI species =="
ROOT_ARGS=()
[[ -z "$ROOT_OVERRIDE" ]] || ROOT_ARGS=(--root "$ROOT_OVERRIDE")
bash "$HERE/metag.sh" audit-reference --config "$CONFIG" "${ROOT_ARGS[@]}" --out "$AUDIT_DIR"
[[ -s "$REPRESENTATIVES" ]] || { echo "ERROR: missing $REPRESENTATIVES" >&2; exit 2; }

# 3) Reuse the v3 planner to make prok/virus/euk genome lists from the existing raw FASTAs.
echo "== Create Sylph build lists =="
"$METAG_PYTHON" "$HERE/f/sylph_profile.py" plan-db --config "$CONFIG" "${ROOT_ARGS[@]}" \
    --representatives "$REPRESENTATIVES" \
    --out "$SYLPH_DIR" --threads "${SYLPH_THREADS:-32}"

# Report how many representative genomes and bases will be read. No genome is copied.
"$METAG_PYTHON" - "$REPRESENTATIVES" "$REFERENCE_RAW" <<'PY'
import csv, os, sys
p, root = sys.argv[1:]
by = {}
with open(p, newline='', encoding='utf-8') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        g = r['group']
        b = 'prok' if g in {'bacteria','archaea'} else ('virus' if g == 'virus' else 'euk')
        x = by.setdefault(b, [0,0])
        x[0] += 1; x[1] += int(r.get('bases') or 0)
for b,(n,bases) in sorted(by.items()):
    print(f'{b}\trepresentatives={n:,}\treference_bases={bases:,}\t~{bases/1e9:.2f} Gb')
PY

# 4) Submit independent bucket jobs so prok/virus/euk can build in parallel.
# v1.0.0 two-stage .syl2db is used for prokaryotes; current v3 intentionally
# keeps virus/euk as .syldb.
QUEUE="${SYLPH_QUEUE:-${HPC_MEDIUM_QUEUE:-38}}"
THREADS="${SYLPH_THREADS:-32}"
MEM="${SYLPH_MEM_GIB:-256}"
HOURS="${SYLPH_HOURS:-120}"
BUCKETS="${SYLPH_BUCKETS:-prok}"
IFS=, read -ra requested_buckets <<< "$BUCKETS"
for bucket in "${requested_buckets[@]}"; do
    case $bucket in prok|virus|euk) ;; *) echo "Invalid bucket: $bucket" >&2; exit 2 ;; esac
done

want_bucket() { [[ ",${BUCKETS}," == *",$1,"* ]]; }

make_job() {
    local bucket=$1 cmd=$2 expected=$3
    local list="$SYLPH_DIR/${bucket}.genomes.txt"
    if [[ ! -s "$list" ]]; then
        echo "SKIP $bucket: empty/missing $list"
        return
    fi
    if [[ -s "$expected" && -s "$SYLPH_DIR/${bucket}.BUILD.COMPLETE" && -s "$expected.sha256" ]]; then
        echo "REUSE $bucket DB: $expected"
        return
    fi
    if [[ -s "$SYLPH_DIR/${bucket}.jobid" && ${submit:-hpc} == hpc ]]; then
        local old_id old_state
        old_id=$(cat "$SYLPH_DIR/${bucket}.jobid")
        old_state=$(bjobs -noheader -o stat "$old_id" 2>/dev/null || true)
        case $old_state in *RUN*|*PEND*|*SUSP*) echo "Existing $bucket job $old_id: $old_state"; return ;; esac
    fi
    local job run_dir run_logs
    run_dir=$("$METAG_PYTHON" "$HERE/../0f/workflow_layout.py" "$JOB_DIR/$bucket")
    run_logs="$LOG_DIR/$bucket/$(basename -- "$run_dir")"
    mkdir -p "$run_logs"
    job="$run_dir/build.cmd"
    local build="$SYLPH_DIR/.build-${bucket}"
    cat > "$job" <<JOB
#!/usr/bin/env bash
#BSUB -J metag_sylph_${bucket}
#BSUB -q $QUEUE
#BSUB -n $THREADS
#BSUB -W $HOURS:00
#BSUB -R "span[hosts=1] rusage[mem=${MEM}GB/host]"
#BSUB -o $run_logs/build.out
#BSUB -e $run_logs/build.err
set -euo pipefail
exec 8>"$SYLPH_DIR/.${bucket}.build.lock"
flock -n 8 || { echo 'Bucket build already active' >&2; exit 2; }
mkdir -p "$build"
cd -- "$build"
$cmd
test -s "$build/$(basename "$expected")"
mv "$build/$(basename "$expected")" "$expected"
test -s "$expected"
sha256sum "$expected" > "$expected.sha256"
printf 'complete\n' > "$SYLPH_DIR/${bucket}.BUILD.COMPLETE"
JOB
    chmod +x "$job"
    if [[ $DO_SUBMIT == FALSE ]]; then
        echo "Plan only; manual submission: bsub < $job"
    elif [[ ${submit:-hpc} == hpc ]]; then
        local submission job_id
        submission=$(bsub < "$job")
        echo "$submission"
        job_id=$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "$submission")
        [[ -n "$job_id" ]] || { echo 'Cannot parse submitted job ID; inspect bjobs before retrying' >&2; exit 2; }
        printf '%s\n' "$job_id" > "$SYLPH_DIR/${bucket}.jobid"
    else
        bash "$job"
    fi
}

if want_bucket prok; then
    make_job prok \
      "\"$SYLPH\" sketch -l \"$SYLPH_DIR/prok.genomes.txt\" -t \"$THREADS\" -c 200 -o prok && \"$SYLPH\" convert-db-two-screen prok.syldb -o prok_v2" \
      "$SYLPH_DIR/prok_v2.syl2db"
fi
if want_bucket virus; then
    make_job virus \
      "\"$SYLPH\" sketch -l \"$SYLPH_DIR/virus.genomes.txt\" -t \"$THREADS\" -c 100 -o virus" \
      "$SYLPH_DIR/virus.syldb"
fi
if want_bucket euk; then
    make_job euk \
      "\"$SYLPH\" sketch -l \"$SYLPH_DIR/euk.genomes.txt\" -t \"$THREADS\" -c 200 -o euk" \
      "$SYLPH_DIR/euk.syldb"
fi

cat > "$SYLPH_DIR/sylph-db.env" <<ENV
SYLPH_EXECUTABLE=$SYLPH
SYLPH_REPRESENTATIVES=$REPRESENTATIVES
SYLPH_PROK_DB=$SYLPH_DIR/prok_v2.syl2db
SYLPH_VIRUS_DB=$SYLPH_DIR/virus.syldb
SYLPH_EUK_DB=$SYLPH_DIR/euk.syldb
ENV

cat <<EOF2

Sylph mini-DB paths
-------------------
Representatives: $REPRESENTATIVES
DB directory:    $SYLPH_DIR
Sylph:           $SYLPH
Environment:     $SYLPH_DIR/sylph-db.env

Priority for ERR14374139: prok_v2.syl2db first. The sample is strongly bacteria-dominated
among the reads that already classify, so prok gives the highest immediate value.
Virus/euk are independent jobs and can finish later.

After prok.BUILD.COMPLETE exists:
  bash $HERE/metag.sh refine --profile-only --config $CONFIG \\
    --representatives $REPRESENTATIVES \\
    --db $SYLPH_DIR/prok_v2.syl2db \\
    --r1 $READS_RAW/england/ERR14374139_1.fastq.gz \\
    --r2 $READS_RAW/england/ERR14374139_2.fastq.gz \\
    --threads 32 --min-ani 95 \\
    --sylph $SYLPH --out $WORK_DIR/qc/ERR14374139-sylph --submit hpc

This script does NOT rebuild $LIBRARY_DIR.
EOF2

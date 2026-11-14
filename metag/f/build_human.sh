#!/usr/bin/env bash
set -euo pipefail

# Build the independent human alignment database used by metag validate.
# This does NOT touch/rebuild annot/metag/lib (the existing ~4.6-TB Kraken library).
#
# Default HPC paths are resolved from metag/config.json + 0f/0os.sh:
#   reference: /work/sph-huangj/data/gen/1hgp/38/
#   bwa-mem2:  /work/sph-huangj/software/metag_parallel/bwa-mem2/2.2.1/
#   logs/jobs: /work/sph-huangj/analysis/metag/{logs,jobs}/small-db/
#
# Usage:
#   cd /work/sph-huangj/scripts/metag
#   ./1.build-lib.sh human
#
# Useful overrides:
#   METAG_CONFIG=/path/config.json HUMAN_QUEUE=38 HUMAN_MEM_GIB=128 ./1.build-lib.sh human
#   GRCH38_URL=https://... ./1.build-lib.sh human

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
for arg in "$@"; do
    if [[ $arg == -h || $arg == --help ]]; then bash "$HERE/f/usage.sh" "build_human"; exit 0; fi
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
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
source "$HERE/f/env.sh"
if [[ $DO_SUBMIT == FALSE ]]; then
    "$METAG_PYTHON" "$HERE/f/plan_small_build.py" --module human --config "$CONFIG" --root "$ROOT_OVERRIDE" --submit "$submit"
    exit 0
fi

# Resolve the exact paths with the same Python/config loader used by metag.
RESOLVED="$("$METAG_PYTHON" - "$HERE" "$CONFIG" "$ROOT_OVERRIDE" <<'PY'
import shlex, sys
sys.path.insert(0, sys.argv[1] + '/f')
from common import load_config
cfg = load_config(sys.argv[2], sys.argv[3] or None)
for key in ('library_dir','work_dir','tools_dir','reads_raw','output_dir','human_reference_dir'):
    print(f'{key.upper()}={shlex.quote(str(cfg[key]))}')
PY
)"
eval "$RESOLVED"

HUMAN_DIR="${HUMAN_DIR:-$HUMAN_REFERENCE_DIR}"
BWA_DIR="${BWA_DIR:-$TOOLS_DIR/bwa-mem2/2.2.1}"
LOG_DIR="$WORK_DIR/logs/small-db"
JOB_DIR="$WORK_DIR/jobs/small-db"
mkdir -p "$HUMAN_DIR" "$BWA_DIR" "$LOG_DIR" "$JOB_DIR"
exec 9>"$HUMAN_DIR/.submission.lock"
flock -n 9 || { echo 'Another human DB preparation is active' >&2; exit 2; }


# NCBI GRCh38 no-alt analysis set: appropriate for read alignment because ALT/patch
# scaffolds that complicate non-alt-aware aligners are omitted.
REF_URL_DEFAULT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REF_URL="${GRCH38_URL:-$REF_URL_DEFAULT}"
REF_GZ="$HUMAN_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REF_FA="$HUMAN_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

"$METAG_PYTHON" "$HERE/f/data_readme.py" --directory "$HUMAN_DIR" --url "$REF_URL" \
    --reads-root "$READS_RAW" --analysis-root "$WORK_DIR" --tools-root "$TOOLS_DIR" --code-root "$HERE"
if [[ -s "$HUMAN_DIR/index.jobid" && ${submit:-hpc} == hpc ]]; then
    old_id=$(cat "$HUMAN_DIR/index.jobid")
    old_state=$(bjobs -noheader -o stat "$old_id" 2>/dev/null || true)
    case $old_state in
        *RUN*|*PEND*|*SUSP*) echo "Existing human index job $old_id: $old_state"; exit 0 ;;
    esac
fi

BWA_VERSION=2.2.1
BWA_URL="${BWA_URL:-https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_VERSION}/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2}"
BWA="$BWA_DIR/bwa-mem2"

fetch() {
    local url=$1 out=$2
    if [[ -s "$out" ]]; then
        echo "REUSE $out"
        return
    fi
    local tmp="${out}.partial"
    if command -v aria2c >/dev/null 2>&1; then
        aria2c -c -x 8 -s 8 -k 1M --file-allocation=none -d "$(dirname "$out")" -o "$(basename "$tmp")" "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl -fLsS --connect-timeout 30 --retry 6 --retry-all-errors --retry-delay 5 -C - -o "$tmp" "$url"
    elif command -v wget >/dev/null 2>&1; then
        wget -c -O "$tmp" "$url"
    else
        echo "ERROR: need aria2c, curl, or wget" >&2
        exit 2
    fi
    [[ -s "$tmp" ]] || { echo "ERROR: empty download: $tmp" >&2; exit 2; }
    mv -f "$tmp" "$out"
}

# 1) Install a pinned BWA-MEM2 copy under the existing metag software root.
if [[ ! -x "$BWA" ]]; then
    echo "== Install bwa-mem2 v${BWA_VERSION} =="
    ARCHIVE="$BWA_DIR/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2"
    BUNDLED="${DIR_SOFTWARE}/metag/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2"
    if [[ -s "$BUNDLED" ]]; then
        cp "$BUNDLED" "$ARCHIVE"
    fi
    fetch "$BWA_URL" "$ARCHIVE"
    TMP="$BWA_DIR/.extract.$$"
    rm -rf "$TMP"; mkdir -p "$TMP"
    tar -xjf "$ARCHIVE" -C "$TMP"
    SRC="$TMP/bwa-mem2-${BWA_VERSION}_x64-linux"
    [[ -x "$SRC/bwa-mem2" ]] || { echo "ERROR: bwa-mem2 binary missing after extraction" >&2; exit 2; }
    # Keep all SIMD-dispatch companion binaries beside bwa-mem2.
    cp -a "$SRC"/. "$BWA_DIR"/
    rm -rf "$TMP"
fi
cp "$HERE/f/bwa_dispatch.sh" "$BWA_DIR/bwa-mem2-launch"
chmod +x "$BWA_DIR/bwa-mem2-launch"
BWA="$BWA_DIR/bwa-mem2-launch"
"$BWA" version

# 2) Download and validate the GRCh38 no-alt FASTA.
echo "== Download GRCh38 no-alt =="
if [[ -s "$REF_GZ" ]]; then
    echo "REUSE $REF_GZ"
else
    "$METAG_PYTHON" "$HERE/f/download_reference.py" --url "$REF_URL" --output "$REF_GZ" \
        --connections "${DOWNLOAD_CONNECTIONS:-8}"
fi
# Validation, decompression and checksums run on the compute node with indexing.

# 3) BWA-MEM2 index. Official BWA-MEM2 documentation gives a high peak-memory
# requirement for index construction, so submit to a 128-GiB+ compute node.
SUFFIXES=(.0123 .amb .ann .bwt.2bit.64 .pac)
complete=0
[[ -s "$HUMAN_DIR/BUILD.COMPLETE" && -s "$REF_FA" ]] && complete=1
for s in "${SUFFIXES[@]}"; do [[ -s "${REF_FA}${s}" ]] || complete=0; done

if [[ "$complete" == 1 ]]; then
    echo "REUSE complete BWA-MEM2 index: $REF_FA"
else
    HUMAN_QUEUE="${HUMAN_QUEUE:-${HPC_MEDIUM_QUEUE:-38}}"
    HUMAN_CORES="${HUMAN_CORES:-8}"
    HUMAN_MEM_GIB="${HUMAN_MEM_GIB:-128}"
    HUMAN_HOURS="${HUMAN_HOURS:-24}"
    RUN_DIR=$("$METAG_PYTHON" "$HERE/../0f/workflow_layout.py" "$JOB_DIR/human")
    RUN_LOGS="$LOG_DIR/human/$(basename -- "$RUN_DIR")"
    mkdir -p "$RUN_LOGS"
    CMD="$RUN_DIR/index.cmd"
    cat > "$CMD" <<JOB
#!/usr/bin/env bash
#BSUB -J metag_grch38_index
#BSUB -q $HUMAN_QUEUE
#BSUB -n $HUMAN_CORES
#BSUB -W $HUMAN_HOURS:00
#BSUB -R "span[hosts=1] rusage[mem=${HUMAN_MEM_GIB}GB/host]"
#BSUB -o $RUN_LOGS/index.out
#BSUB -e $RUN_LOGS/index.err
set -euo pipefail
exec 8>"$HUMAN_DIR/.build.lock"
flock -n 8 || { echo 'Human index build already active' >&2; exit 2; }
gzip -t "$REF_GZ"
if [[ ! -s "$REF_FA" ]]; then
  gzip -dc "$REF_GZ" > "${REF_FA}.partial"
  grep -qm1 '^>' "${REF_FA}.partial"
  mv "${REF_FA}.partial" "$REF_FA"
fi
sha256sum "$REF_GZ" "$REF_FA" > "$HUMAN_DIR/reference.sha256"
printf '%s\n' "$REF_URL" > "$HUMAN_DIR/reference.url.txt"
"$BWA" index -p "${REF_FA}.building" "$REF_FA"
for s in .0123 .amb .ann .bwt.2bit.64 .pac; do
  test -s "${REF_FA}.building\${s}"
done
for s in .0123 .amb .ann .bwt.2bit.64 .pac; do
  mv "${REF_FA}.building\${s}" "${REF_FA}\${s}"
done
sha256sum "$REF_FA" "${REF_FA}".0123 "${REF_FA}".amb "${REF_FA}".ann "${REF_FA}".bwt.2bit.64 "${REF_FA}".pac > "$HUMAN_DIR/bwa-mem2-index.sha256"
printf 'complete\n' > "$HUMAN_DIR/BUILD.COMPLETE"
JOB
    chmod +x "$CMD"
    if [[ $DO_SUBMIT == FALSE ]]; then
        echo "Plan only; manual submission: bsub < $CMD"
    elif [[ ${submit:-hpc} == hpc ]]; then
        echo "== Submit BWA-MEM2 index job =="
        submission=$(bsub < "$CMD")
        echo "$submission"
        job_id=$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "$submission")
        [[ -n "$job_id" ]] || { echo 'Cannot parse submitted job ID; inspect bjobs before retrying' >&2; exit 2; }
        printf '%s\n' "$job_id" > "$HUMAN_DIR/index.jobid"
    else
        echo "== Run BWA-MEM2 index locally =="
        bash "$CMD"
    fi
fi

printf 'HUMAN_HOST_REFERENCE=%q\nHUMAN_HOST_EXECUTABLE=%q\n' "$REF_FA" "$BWA" > "$HUMAN_DIR/human-db.env"

cat <<EOF2

Human mini-DB paths
-------------------
FASTA:       $REF_FA
BWA-MEM2:    $BWA
Index dir:   $HUMAN_DIR
Environment: $HUMAN_DIR/human-db.env

After BUILD.COMPLETE exists, validate the existing human candidates:
  bash $HERE/metag.sh validate --config $CONFIG --sample england/ERR14374050 --submit hpc

Then screen the full original sample for missed human-compatible reads:
  bash $HERE/metag.sh validate --config $CONFIG --sample england/ERR14374050 --scope full --submit hpc

Use --reference and --bwa-mem2 to override the default paths when needed.

This script does NOT rebuild $LIBRARY_DIR.
EOF2

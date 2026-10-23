#!/usr/bin/env bash
# Independent individual-genome HERV census. See 0.README.md for dependencies.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=herv
# shellcheck source=runtime/standalone.sh
source "$ROOT/runtime/standalone.sh"

usage() {
    cat <<'HELP'
Usage: ./1.herv.sh [OPTIONS]

Reconstruct chr1-22 from phased VCF/PGEN and annotate reference/hap1/hap2 with
chunked RepeatMasker. Resume verified checkpoints and report diploid ERV fractions.
No gnomAD scan, selection analysis, PCA or allele3.sh is required.

  --gen NAME / --sample ID   Dataset / individual [1kg / HG00096]
  --dir-analysis-root DIR   Analysis root [/mnt/d/analysis/allele3]
  --dir-out DIR, --out DIR   Sample output [ROOT/herv/GEN/SAMPLE]
  --input-format auto|vcf|pgen  Reuse checkpoint; otherwise prefer available PGEN [auto]
  --dir-vcf DIR             VCF source [/mnt/f/gen/1kg/38/vcf]
  --dir-pfile DIR           PGEN/PVAR/PSAM source [/mnt/f/gen/1kg/38/pfile]
  --fasta FILE              Indexed reference [/mnt/f/gen/1hgp/38/GRCH38.fasta.gz]
  --library FILE            Repeat library [/mnt/f/annot/dfam_4.0_human/human.curated.fa]
  --repeatmasker FILE       Reuse checkpoint executable, then autodetect
  --pgen-python FILE        Python with numpy, pgenlib and zstandard [autodetect]
  --jobs N / --cores N      Concurrent chunks / CPU affinity limit [8 / 16]
  --chunk-bp N / --halo N    Core length / flanking context [5000000 / 20000]
  --dir-scratch DIR         Temporary sequences [/tmp/allele3-herv-assembly]
  --memory-limit-gb N       Process-tree RAM cap [16]; swap disabled; requires systemd
  --python FILE            Python interpreter [active conda / installed HERV env / python3]
  --screen-log-interval N / --verbose-screen   Console controls; full logs in OUT/logs
  --dry-run                Print the resolved command; do not create files or run analysis
  --help-full              Native Python options
  -h, --help               This help, without loading Python or configuration

Aliases: --vcf-dir, --pfile-dir and --scratch also work.
Config: HERV_CONFIG or ./herv.env; see 0.README.md (environment and configuration). CLI overrides config.
Existing output defaults are retained for checkpoint reuse; use --out elsewhere.
Portable file bundle and input conventions: 0.README.md.
HELP
}

main() {
    local arg
    for arg in "$@"; do
        case "$arg" in -h|--help|help) usage; return 0 ;; esac
    done
    local config=${HERV_CONFIG:-"$ROOT/herv.env"}
    if [[ -f $config ]]; then
        set -a
        # shellcheck disable=SC1090
        source "$config"
        set +a
    fi
    # Legacy environment names remain aliases; allele.env is never loaded here.
    local analysis_root=${HERV_ANALYSIS_ROOT:-${ALLELE_ANALYSIS_ROOT:-/mnt/d/analysis/allele3}}
    local out=${HERV_OUT:-${ALLELE_HERV_OUT:-}} gen=${HERV_GEN:-1kg} sample=${HERV_SAMPLE:-HG00096}
    local vcf=${HERV_VCF_DIR:-${ALLELE_HERV_VCF_DIR:-/mnt/f/gen/1kg/38/vcf}}
    local fasta=${HERV_FASTA:-${ALLELE_FASTA:-/mnt/f/gen/1hgp/38/GRCH38.fasta.gz}}
    local jobs=${HERV_JOBS:-${ALLELE_HERV_JOBS:-8}} cores=${HERV_CORES:-${ALLELE_HERV_CORES:-16}}
    local python=${HERV_PYTHON:-${ALLELE_PYTHON:-}}
    local memory=${HERV_MEMORY_LIMIT_GB:-${ALLELE_MEMORY_LIMIT_GB:-16}}
    local interval=${HERV_SCREEN_LOG_INTERVAL:-${ALLELE_SCREEN_LOG_INTERVAL:-300}}
    local verbose=0 dry_run=0 full_help=0
    local -a extra=() args=()
    [[ -z ${HERV_PFILE_DIR:-} ]] || extra+=(--pfile-dir "$HERV_PFILE_DIR")
    [[ -z ${HERV_LIBRARY:-} ]] || extra+=(--library "$HERV_LIBRARY")
    [[ -z ${HERV_SCRATCH:-} ]] || extra+=(--scratch "$HERV_SCRATCH")
    while (($#)); do
        if [[ $1 == --*=* ]]; then
            arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"
        fi
        case "$1" in
            --dry-run) dry_run=1; shift ;;
            --help-full) full_help=1; shift ;;
            --verbose-screen) verbose=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; shift 2 ;;
            --dir-out|--out) need_value "$@"; out=$2; shift 2 ;;
            --gen) need_value "$@"; gen=$2; shift 2 ;;
            --sample) need_value "$@"; sample=$2; shift 2 ;;
            --dir-vcf|--vcf-dir) need_value "$@"; vcf=$2; shift 2 ;;
            --fasta|--fasta38) need_value "$@"; fasta=$2; shift 2 ;;
            --jobs) need_value "$@"; jobs=$2; shift 2 ;;
            --cores) need_value "$@"; cores=$2; shift 2 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --input-format|--pfile-dir|--dir-pfile|--pgen-python|--library|--repeatmasker|--chunk-bp|--halo|--scratch|--dir-scratch)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    [[ $gen =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail '--gen and --sample must be single names (letters, digits, dots, underscores or hyphens)'
    positive_integer --jobs "$jobs"
    positive_integer --cores "$cores"
    ((jobs <= 8)) || fail 'Use at most 8 workers on this shared host'
    local conda_root=${HERV_CONDA_ROOT:-${ALLELE_CONDA_ROOT:-${METAG_CONDA_ROOT:-$HOME/anaconda3}}}
    python=$(workflow_python "$python" "$conda_root/envs/herv-pgen/bin/python" \
        "$conda_root/envs/allele3-herv-pgen/bin/python" "$conda_root/envs/allele3-evolution/bin/python")
    out=${out:-$analysis_root/herv/$gen/$sample}
    args=(--out "$out" --gen "$gen" --sample "$sample" --fasta "$fasta"
          --vcf-dir "$vcf" --jobs "$jobs" --cores "$cores" "${extra[@]}")
    if ((full_help)); then args=(--help); fi
    workflow_launch herv "$python" "$ROOT/f/herv/herv_assembly_job.py" "$out" \
        "$memory" "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
}

main "$@"

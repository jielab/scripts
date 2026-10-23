#!/usr/bin/env bash
# Standalone fixed-basis PCA and PC1/PC2/age visualization.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=pca3d
# shellcheck source=runtime/standalone.sh
source "$ROOT/runtime/standalone.sh"

usage() {
    cat <<'HELP'
Usage: ./0.pca3D.sh [run|plan|render] [OPTIONS]

  run (default)   Prepare metadata if needed, fit/project PCA, then render HTML/PNG/PDF
  plan            Write the PLINK2 workflow without executing it
  render          Render existing shared PCA coordinates; no PGEN needed without overlays

The 3D axes are PC1, PC2 and age (BP), not PC1/PC2/PC3. PCA planning currently
uses human autosomal GRCh37 PGEN data. No gnomAD/selection/Atlas stage is required.

Inputs and output:
  --pgen-prefix PREFIX      Input .pgen/.pvar/.psam prefix
  --metadata TSV           sample_id, cohort, age_bp, include, stratum; see 0.README.md
  --hq-variants TSV         Background variants: #CHROM,POS,ID,REF,ALT (GRCh37)
  --eigenvec FILE           Shared .sscore or .eigenvec (required for render)
  --out DIR, --dir-out DIR  Output root [/mnt/d/analysis/pca3d]; PCA files under OUT/pca
  --dir-analysis-root DIR   Parent of default pca3d output [/mnt/d/analysis]
  --html FILE               HTML output [OUT/pca/pca_time.html]
  --dir-ancient DIR         Original dataset layout [/mnt/f/gen/ancient]
  --supplement-xlsx FILE    Original dataset metadata adapter; custom panels use --metadata
  --qpadm TSV / --ancestry-model 3way|4way   Optional original metadata adapter inputs
  --dir-pgen-cache DIR / --no-pgen-cache    Optional verified execution copy

PCA settings:
  --pca-basis joint_balanced|modern_projection [joint_balanced]
  --pca-maf F / --pca-geno F [0.05 / 0.1]
  --pca-components N / --pca-max-per-group N [10 / 200]
  --pca-time-bin-years N / --seed N [2000 / 20261021]
  --pca-exclude-bed FILE    Additional GRCh37 exclusions
  --target-variants TSV     Optional background exclusions with chrom,pos[,build]
  --pca-target-flank-bp N   Target exclusion flank [1000000]
  --pca-threads N / --pca-memory-mb N [4 / 8192]; --plink2 FILE [plink2]

Optional plot settings:
  --selected-variants TSV  Verified PVAR-indexed targets; at most 12 are displayed
                           Also excluded from PCA unless --target-variants is supplied
  --pca-max-variants N / --count-kind hardcall|dosage [12 / hardcall]
  --pca-basis-label TEXT / --title TEXT / --simulated

Execution:
  --python FILE / --memory-limit-gb N [16; entire process tree, systemd required]
  --screen-log-interval N / --verbose-screen [300 seconds / off]
  --dry-run                Print command only; no files, input reads or computation
  --help-full              Native Python options
  -h, --help               Help without Python or configuration

Config: PCA3D_CONFIG or ./pca3d.env; see 0.README.md (environment and configuration).
Existing allele3 pipeline and pca-plan/pca-time commands use the same calculation library.
HELP
}

main() {
    local arg
    for arg in "$@"; do
        case "$arg" in -h|--help|help) usage; return 0 ;; esac
    done
    local config=${PCA3D_CONFIG:-"$ROOT/pca3d.env"}
    if [[ -f $config ]]; then
        set -a
        # shellcheck disable=SC1090
        source "$config"
        set +a
    fi
    local command=run
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    case "$command" in
        pca-plan) command=plan ;; pca-time) command=render ;;
        run|plan|render) ;; *) fail "Unknown command: $command" ;;
    esac
    local analysis_root=${PCA3D_ANALYSIS_ROOT:-/mnt/d/analysis} out=${PCA3D_OUT:-}
    local ancient=${PCA3D_ANCIENT_DIR:-/mnt/f/gen/ancient} cache=${PCA3D_PGEN_CACHE_DIR:-}
    local python=${PCA3D_PYTHON:-} memory=${PCA3D_MEMORY_LIMIT_GB:-16}
    local threads=${PCA3D_THREADS:-4} pca_memory=${PCA3D_MEMORY_MB:-8192}
    local interval=${PCA3D_SCREEN_LOG_INTERVAL:-300} verbose=0 dry_run=0 full_help=0
    local -a extra=() args=()
    [[ -z ${PCA3D_PGEN_PREFIX:-} ]] || extra+=(--pgen-prefix "$PCA3D_PGEN_PREFIX")
    [[ -z ${PCA3D_METADATA:-} ]] || extra+=(--metadata "$PCA3D_METADATA")
    [[ -z ${PCA3D_HQ_VARIANTS:-} ]] || extra+=(--hq-variants "$PCA3D_HQ_VARIANTS")
    [[ -z ${PCA3D_PLINK2:-} ]] || extra+=(--plink2 "$PCA3D_PLINK2")
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
            --dir-ancient|--ancient-dir) need_value "$@"; ancient=$2; shift 2 ;;
            --dir-pgen-cache) need_value "$@"; cache=$2; shift 2 ;;
            --no-pgen-cache) cache=''; shift ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --pca-threads|--threads|--cores) need_value "$@"; threads=$2; shift 2 ;;
            --pca-memory-mb) need_value "$@"; pca_memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --simulated) extra+=("$1"); shift ;;
            --pgen-prefix|--metadata|--hq-variants|--eigenvec|--html|--supplement-xlsx|--qpadm|--ancestry-model|--selected-variants|--target-variants|--pca-basis|--pca-basis-label|--pca-maf|--pca-geno|--pca-max-per-group|--pca-components|--pca-target-flank-bp|--pca-time-bin-years|--pca-exclude-bed|--pca-max-variants|--count-kind|--plink2|--seed|--title)
                need_value "$@"; extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help)" ;;
        esac
    done
    positive_integer --memory-limit-gb "$memory"
    positive_integer --pca-threads "$threads"
    positive_integer --pca-memory-mb "$pca_memory"
    ((threads <= 8)) || fail 'Use at most 8 threads on this shared host'
    if ((pca_memory > memory * 512)); then pca_memory=$((memory * 512)); fi
    ((pca_memory >= 640)) || fail 'PLINK needs at least 640 MiB; increase --memory-limit-gb / --pca-memory-mb'
    local conda_root=${PCA3D_CONDA_ROOT:-$HOME/anaconda3}
    python=$(workflow_python "$python" "$conda_root/envs/pca3d/bin/python" \
        "$conda_root/envs/allele3-evolution/bin/python")
    out=${out:-$analysis_root/pca3d}
    args=("$command" --out "$out" --ancient-dir "$ancient" --pca-threads "$threads"
        --pca-memory-mb "$pca_memory" "${extra[@]}")
    [[ -z $cache ]] || args+=(--dir-pgen-cache "$cache")
    if ((full_help)); then args=(--help); fi
    workflow_launch pca3d "$python" "$ROOT/f/pca3d/cli.py" "$out" "$memory" \
        "$interval" "$verbose" "$dry_run" "$((1-full_help))" "${args[@]}"
}

main "$@"

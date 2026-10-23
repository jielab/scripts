#!/usr/bin/env bash
# allele3: modern gnomAD scans and ancient-DNA / PCA / sparse Atlas analysis.
# All paths and common settings are visible here; allele.env may override them.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<'HELP'
Usage: ./allele3.sh [COMMAND] [OPTIONS]

Commands:
  pipeline (default)  Resume modern scan/summary -> ancient analysis -> PCA -> report
                      Completed results are reused; AlphaGenome is off by default.
  run                 Scan modern VCFs and generate local summaries only
  scan / summarize    Run either modern stage separately; resume completed shards
  check               Validate modern inputs without scanning variants
  herv                Compatibility alias for ./1.herv.sh (independent HERV workflow)
  evolve              Ancient metadata -> selection join -> candidates -> PGEN
                      trajectories -> Atlas/PCA plans -> evolution report
  prepare-ancient / join-selection / select-candidates / panel-index
  trajectories / audit-ancient / audit-counts
  pca-plan / pca-time / report-evolution   Standalone PCA: ./0.pca3D.sh
  atlas-plan / atlas-query / atlas-summarize / compare-functions
  validate-ld / prepare-erv / summarize-architecture / ukb-mfi
  build-universes / select-u3 / final-report
  annotate / make-callable   Local score annotation / callable BED utilities

Directories (input data remain on F:, all default analysis output is on D:):
  --dir-analysis-root DIR  Analysis root [/mnt/d/analysis/allele3]
  --dir-out DIR            Output [pipeline/evolution: ROOT/evolution; run: ROOT/gnomad_local]
  --dir-vcf DIR            Modern chromosome VCFs [/mnt/f/annot/gnomad/sites]
  --dir-cache DIR          Modern shard cache [OUT/shards]
  --dir-ancient DIR        Ancient panel and metadata [/mnt/f/gen/ancient]
  --dir-modern DIR         Modern scan/cache for pipeline or evolve [ROOT/gnomad_local]
  --dir-atlas-cache DIR    Existing Atlas results/cache, if supplied
  --dir-pgen-cache DIR     Verified PGEN execution copy [local config: $HOME/.cache/allele3-evolution/pgen]
  --no-pgen-cache          Read the original PGEN directly for PCA

Modern scan settings:
  --vcf FILE --region CHR:START-END   Scan one indexed file / region
  --chr CSV               Directory chromosomes, e.g. 1,2,22,X [autosomes 1-22]
  --fasta FILE            GRCh38 reference [/mnt/f/gen/1hgp/38/GRCH38.fasta.gz]
  --jobs N                Concurrent chromosome scans [2]
  --backend NAME          auto, bcftools, python [auto]
  --min-an N              Minimum per-ALT AN [100]
  --min-af F              Minimum ALT frequency [0; retains singletons]
  --thresholds CSV        ALT count thresholds [1,2,5,10]
  --vep                   Read existing VEP consequences [off]
  --with-architecture     LD / joined ERV sidecar counters (original QC preserved)
  --ld-block-dir DIR / --architecture-build 37|38 / --erv-joined FILE
  --architecture-prefix NAME / --architecture-min-sites N [10000]
  --architecture-segdup-bed FILE / --architecture-lcr-bed FILE
  --architecture-low-mappability-bed FILE   Excluded low-mappability intervals
  --candidate-storage MODE  deferred or eager [deferred]
  --no-plots              Skip optional static plots

Ancient / optional analysis settings:
  --fasta37 FILE          GRCh37 reference [/mnt/f/gen/1hgp/37/GRCH37.fasta.gz]
  --chain38to37 FILE / --chain37to38 FILE   Verified reciprocal chain files
  --join-strategy NAME    selection-first or modern-first [selection-first]
  --max-sites N / --max-variants N         Candidate / Atlas budgets [200 / 600]
  --with-universes        U0 full / U1 stratified / U2 all strong / independent U3
  --u1-sites N / --u3-sites N / --controls-per-case N [20000 / 3000 / 4]
  --with-ukb / --ukb-mfi-dir DIR   Full U0 exact GRCh37 MFI scan before PGEN
  --metag-census FILE      Reuse the completed hg38 census read-only
  --with-pca / --no-pca   Execute PCA [on for pipeline; off for explicit evolve]
  --pca-threads N / --pca-memory-mb N      PLINK budgets [4 / 8192]
  --with-atlas / --no-atlas  Execute bounded AlphaGenome queries [off]
  --api-key-file FILE     External secret [local config: /mnt/d/files/alphaG.api]
  --workers N / --requests-per-second N   Atlas concurrency / rate [2 / 2]

Execution and safety:
  --screen-log-interval N  Repeated screen progress interval [300 seconds]
  --verbose-screen       Show full output for pipeline/run/scan/evolve
                         Full output is always saved in OUT/logs/ for these commands.
  --memory-limit-gb N     Hard RAM cap for the ENTIRE process tree [16 GiB]
                         Swap is disabled in this cgroup; OOM kills this job.
                         Requires user systemd + cgroup v2; fails closed otherwise.
  --python FILE          Python interpreter [allele3-evolution, then python3]
  --preflight            Check environment, resource isolation and input headers
  --dry-run              Print the resolved command and limits; execute nothing
  -h, --help             This Bash help; Python/environment is not required
  --help-full            Full native options for the selected command
  --version              Installed code version

Legacy --out, --vcf-dir, --cache-dir, --ancient-dir and --modern-out also work.
All native Python options remain available. Config: ALLELE_CONFIG or ./allele.env.

Examples:
  ./allele3.sh -h
  ./allele3.sh herv --sample HG00096 --dry-run
  ./allele3.sh --preflight
  ./allele3.sh            # Resume modern results, ancient analysis and PCA; no API
  ./allele3.sh --with-atlas  # Also query AlphaGenome and compare functional effects
  ./allele3.sh run --vcf /mnt/f/annot/gnomad/sites/chr22.vcf.bgz \
      --region chr22:20000001-22000000 --vep \
      --dir-out /mnt/d/analysis/allele3/pilot_chr22
  ./allele3.sh evolve --dir-modern /mnt/d/analysis/allele3/gnomad_local \
      --dir-out /mnt/d/analysis/allele3/evolution --with-pca
  ./allele3.sh atlas-query --dir-out /mnt/d/analysis/allele3/evolution \
      --api-key-file /mnt/d/files/alphaG.api
HELP
}

fail() { printf 'allele3: %s\n' "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || fail "$1 requires a value"; }
positive_integer() { [[ $2 =~ ^[1-9][0-9]*$ && ${#2} -le 6 ]] || fail "$1 must be a positive integer"; }

main() {
    # Help must work even when the configured Python is missing or config is broken.
    local arg
    for arg in "$@"; do
        case "$arg" in -h|--help|help)
            if [[ ${1:-} == herv ]]; then exec bash "$ROOT/1.herv.sh" --help; else usage; fi
            return 0 ;; esac
    done
    local config=${ALLELE_CONFIG:-"$ROOT/allele.env"}
    if [[ -f $config ]]; then
        set -a
        # shellcheck disable=SC1090
        source "$config"
        set +a
    fi

    # The independent HERV entry owns all its options and execution.
    if [[ ${1:-} == herv ]]; then shift; exec bash "$ROOT/1.herv.sh" "$@"; fi

    # Public defaults. Explicit CLI options override configuration and environment.
    local command=pipeline family=modern
    local analysis_root=${ALLELE_ANALYSIS_ROOT:-/mnt/d/analysis/allele3}
    local out=${ALLELE_LOCAL_OUT:-${ALLELE_OUT:-}} modern_out=${ALLELE_MODERN_OUT:-}
    local vcf_dir=${ALLELE_VCF_DIR:-/mnt/f/annot/gnomad/sites} cache=${ALLELE_CACHE_DIR:-}
    local ancient_dir=${ALLELE_ANCIENT_DIR:-/mnt/f/gen/ancient}
    local pgen_cache=${ALLELE_PGEN_CACHE_DIR:-}
    local fasta=${ALLELE_FASTA:-/mnt/f/gen/1hgp/38/GRCH38.fasta.gz}
    local fasta37=${ALLELE_FASTA37:-/mnt/f/gen/1hgp/37/GRCH37.fasta.gz}
    local chain38=${ALLELE_CHAIN38TO37:-/mnt/f/gen/chains/hg38ToHg19.over.chain.gz}
    local chain37=${ALLELE_CHAIN37TO38:-/mnt/f/gen/chains/hg19ToHg38.over.chain.gz}
    local python=${ALLELE_PYTHON:-} jobs=${ALLELE_JOBS:-2}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} pca_memory=8192 pca_threads=4
    local screen_interval=${ALLELE_SCREEN_LOG_INTERVAL:-300} verbose_screen=0
    local min_an=100 min_af=0 thresholds=1,2,5,10 backend=auto storage=deferred
    local max_sites=200 max_variants=600 strategy=selection-first
    local pca=auto with_atlas=0 api_key_file=${ALLELE_API_KEY_FILE:-}
    local vcf='' chromosomes='' region='' dry_run=0 preflight=0 full_help=0 version=0
    local explicit_vcf_dir=0
    local -a final_args=() extra=() modern_extra=() evolution_extra=() args=() modern_args=() evolution_args=()
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    case "$command" in
        validate-ld|prepare-erv|summarize-architecture|ukb-mfi|build-universes|select-u3|final-report) family=final; final_args=("$command" "$@"); set -- ;;
        pipeline) family=pipeline; out=${ALLELE_EVOLUTION_OUT:-} ;;
        run|scan|summarize|check|annotate|make-callable) ;;
        evolve|prepare-ancient|join-selection|select-candidates|panel-index|trajectories|atlas-plan|atlas-query|atlas-summarize|compare-functions|pca-plan|pca-time|report-evolution|audit-ancient|audit-counts)
            family=evolution; out=${ALLELE_EVOLUTION_OUT:-} ;;
        *) fail "Unknown command: $command (see --help)" ;;
    esac

    # FINAL local profile is opt-in through config or explicit switches.
    if [[ ${ALLELE_WITH_ARCHITECTURE:-0} == 1 && ( $family == modern || $family == pipeline ) ]]; then
        modern_extra+=(--with-architecture --ld-block-dir "${ALLELE_LD_BLOCK_DIR:-/mnt/f/LD/block}"
            --architecture-build 38 --erv-joined "${ALLELE_ERV_JOINED:-/mnt/f/annot/repeats/hg38/rmskJoinedCurrent.txt.gz}"
            --architecture-annotation-cache "${ALLELE_ANNOTATION_CACHE:-/mnt/d/analysis/allele3/final/annotations}")
        if [[ $family == modern ]]; then extra+=("${modern_extra[@]}"); fi
    fi
    if [[ ${ALLELE_WITH_UNIVERSES:-0} == 1 && ( $family == evolution || $family == pipeline ) ]]; then
        evolution_extra+=(--with-universes --with-ukb --ld-block-dir "${ALLELE_LD_BLOCK_DIR:-/mnt/f/LD/block}"
            --u1-sites "${ALLELE_U1_SITES:-20000}" --u3-sites "${ALLELE_U3_SITES:-3000}"
            --ukb-mfi-dir "${ALLELE_UKB_MFI_DIR:-/mnt/f/gen/ukb/37/imp_mfi}")
        [[ -z ${ALLELE_METAG_CENSUS:-} ]] || evolution_extra+=(--metag-census "$ALLELE_METAG_CENSUS")
        if [[ $family == evolution ]]; then extra+=("${evolution_extra[@]}"); fi
        max_variants=${ALLELE_MAX_VARIANTS:-9000}
    fi
    while (($#)); do
        # Normalize --option=value without losing spaces or shell metacharacters.
        if [[ $1 == --*=* ]]; then
            arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"
        fi
        case "$1" in
            --dry-run) dry_run=1; shift ;;
            --preflight) preflight=1; shift ;;
            --help-full) full_help=1; shift ;;
            --version) version=1; shift ;;
            --screen-log-interval) need_value "$@"; screen_interval=$2; shift 2 ;;
            --verbose-screen) verbose_screen=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; shift 2 ;;
            --dir-out|--out) need_value "$@"; out=$2; shift 2 ;;
            --dir-vcf|--vcf-dir) need_value "$@"; vcf_dir=$2; explicit_vcf_dir=1; shift 2 ;;
            --dir-cache|--cache-dir) need_value "$@"; cache=$2; shift 2 ;;
            --dir-ancient|--ancient-dir) need_value "$@"; ancient_dir=$2; shift 2 ;;
            --dir-modern|--modern-out) need_value "$@"; modern_out=$2; shift 2 ;;
            --dir-pgen-cache) need_value "$@"; pgen_cache=$2; shift 2 ;;
            --no-pgen-cache) pgen_cache=""; shift ;;
            --dir-atlas-cache|--atlas-cache) need_value "$@"; extra+=(--atlas-cache "$2"); shift 2 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --jobs) need_value "$@"; jobs=$2; shift 2 ;;
            --cores) need_value "$@"; jobs=$2; shift 2 ;;
            --pca-threads) need_value "$@"; pca_threads=$2; shift 2 ;;
            --pca-memory-mb) need_value "$@"; pca_memory=$2; shift 2 ;;
            --vcf) need_value "$@"; vcf=$2; shift 2 ;;
            --chr) need_value "$@"; chromosomes=$2; shift 2 ;;
            --region) need_value "$@"; region=$2; shift 2 ;;
            --fasta|--fasta38) need_value "$@"; fasta=$2; shift 2 ;;
            --fasta37) need_value "$@"; fasta37=$2; shift 2 ;;
            --chain38to37) need_value "$@"; chain38=$2; shift 2 ;;
            --chain37to38) need_value "$@"; chain37=$2; shift 2 ;;
            --min-an) need_value "$@"; min_an=$2; shift 2 ;;
            --min-af) need_value "$@"; min_af=$2; shift 2 ;;
            --thresholds) need_value "$@"; thresholds=$2; shift 2 ;;
            --backend) need_value "$@"; backend=$2; shift 2 ;;
            --candidate-storage) need_value "$@"; storage=$2; shift 2 ;;
            --max-sites) need_value "$@"; max_sites=$2; shift 2 ;;
            --max-variants) need_value "$@"; max_variants=$2; shift 2 ;;
            --join-strategy) need_value "$@"; strategy=$2; shift 2 ;;
            --with-pca) pca=1; shift ;;
            --no-pca) pca=0; shift ;;
            --with-atlas) with_atlas=1; shift ;;
            --no-atlas) with_atlas=0; shift ;;
            --api-key-file) need_value "$@"; api_key_file=$2; shift 2 ;;
            --api-key-env)
                need_value "$@"; api_key_file=''; extra+=("$1" "$2"); evolution_extra+=("$1" "$2"); shift 2 ;;
            --with-architecture)
                extra+=("$1"); modern_extra+=("$1"); shift ;;
            --with-universes|--with-ukb)
                extra+=("$1"); evolution_extra+=("$1"); shift ;;
            --ld-block-dir)
                need_value "$@"; extra+=("$1" "$2"); modern_extra+=("$1" "$2"); evolution_extra+=("$1" "$2"); shift 2 ;;
            --architecture-build|--erv-joined|--architecture-prefix|--architecture-min-sites|--architecture-top-fraction|--architecture-annotation-cache|--architecture-segdup-bed|--architecture-lcr-bed|--architecture-low-mappability-bed)
                need_value "$@"; extra+=("$1" "$2"); modern_extra+=("$1" "$2"); shift 2 ;;
            --u1-sites|--u3-sites|--controls-per-case|--ukb-mfi-dir|--metag-census)
                need_value "$@"; extra+=("$1" "$2"); evolution_extra+=("$1" "$2"); shift 2 ;;
            --include-unfiltered|--vep|--write-all-sites)
                extra+=("$1"); modern_extra+=("$1"); shift ;;
            --no-plots)
                extra+=("$1"); modern_extra+=("$1"); evolution_extra+=("$1"); shift ;;
            --include-mhc|--allow-selection-warnings|--skip-trajectories|--simulated)
                extra+=("$1"); evolution_extra+=("$1"); shift ;;
            --window-size|--populations|--min-pop-an|--min-pop-af|--exclude-flags|--callable-bed|--build|--source-label|--min-training-sites|--scores|--score-column|--score-format|--score-label|--score-version|--score-build|--score-backend|--all-sites-an|--callable-min-an|--output-bed)
                need_value "$@"; extra+=("$1" "$2"); modern_extra+=("$1" "$2"); shift 2 ;;
            --pgen-prefix|--selection|--counts-tsv|--hq-variants|--supplement-xlsx|--qpadm|--ancestry-model|--metadata|--intersection|--candidates|--selected-variants|--strong-posterior|--weak-posterior|--min-extra-ac|--block-bp|--seed|--bins|--min-samples|--weights|--count-kind|--atlas-plan|--scorers|--ontology-terms|--workers|--requests-per-second|--retries|--timeout|--geometry|--eigenvec|--pca-basis|--pca-basis-label|--pca-maf|--pca-geno|--pca-max-per-group|--pca-components|--pca-target-flank-bp|--pca-exclude-bed|--pca-max-variants|--plink2|--html)
                need_value "$@"; extra+=("$1" "$2"); evolution_extra+=("$1" "$2"); shift 2 ;;
            *) fail "Unknown option: $1 (see --help or --help-full)" ;;
        esac
    done

    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$screen_interval"
    positive_integer --jobs "$jobs"
    positive_integer --pca-threads "$pca_threads"
    positive_integer --pca-memory-mb "$pca_memory"
    ((memory <= 32)) || fail 'Memory cap must be at most 32 GiB on this shared WSL host'
    ((jobs <= 8 && pca_threads <= 8)) || fail 'Use at most 8 workers/threads on this shared host'
    # Leave half of the cgroup budget for Python, PGEN metadata and other children.
    if ((pca_memory > memory * 512)); then pca_memory=$((memory * 512)); fi
    [[ -z $vcf || $explicit_vcf_dir == 0 ]] || fail '--vcf and --dir-vcf are mutually exclusive'
    [[ -z $vcf || -z $chromosomes ]] || fail 'Use --region with --vcf; --chr selects directory inputs'
    [[ -z $region || -n $vcf ]] || fail '--region requires --vcf'

    if [[ -z $python ]]; then
        if [[ -x $HOME/anaconda3/envs/allele3-evolution/bin/python3 ]]; then
            python=$HOME/anaconda3/envs/allele3-evolution/bin/python3
        else python=python3; fi
    fi
    modern_out=${modern_out:-${ALLELE_LOCAL_OUT:-${ALLELE_OUT:-$analysis_root/gnomad_local}}}
    if ((preflight)); then
        [[ $family != evolution ]] || fail '--preflight uses modern input check; use audit-ancient for the ancient panel'
        if [[ $family == pipeline ]]; then out=$modern_out; fi
        command=check; family=modern
    fi
    local entry=$ROOT/f/cli.py modern_command=$command evolution_command=$command
    if [[ $family == pipeline ]]; then modern_command=run; evolution_command=evolve; fi
    if [[ $family == modern || $family == pipeline ]]; then
        if [[ $family == modern ]]; then modern_out=${out:-$analysis_root/gnomad_local}; fi
        modern_args=("$modern_command" --dir-out "$modern_out" --fasta "$fasta" --jobs "$jobs"
              --min-an "$min_an" --min-af "$min_af" --thresholds "$thresholds"
              --backend "$backend" --candidate-storage "$storage")
        if [[ -n $vcf ]]; then modern_args+=(--vcf "$vcf"); else modern_args+=(--dir-vcf "$vcf_dir"); fi
        [[ -z $chromosomes ]] || modern_args+=(--chr "$chromosomes")
        [[ -z $region ]] || modern_args+=(--region "$region")
        [[ -z $cache ]] || modern_args+=(--dir-cache "$cache")
    fi
    if [[ $family == evolution || $family == pipeline ]]; then
        entry=$ROOT/f/evolve_cli.py
        out=${out:-$analysis_root/evolution}
        evolution_args=("$evolution_command" --out "$out" --ancient-dir "$ancient_dir" --modern-out "$modern_out"
              --fasta38 "$fasta" --fasta37 "$fasta37" --chain38to37 "$chain38" --chain37to38 "$chain37"
              --join-strategy "$strategy" --max-sites "$max_sites" --max-variants "$max_variants"
              --pca-threads "$pca_threads" --pca-memory-mb "$pca_memory")
        [[ -z $pgen_cache ]] || evolution_args+=(--dir-pgen-cache "$pgen_cache")
        if [[ $pca == 1 || ( $pca == auto && $family == pipeline ) ]]; then evolution_args+=(--with-pca); fi
        if ((with_atlas)); then evolution_args+=(--with-atlas); fi
        [[ -z $api_key_file ]] || evolution_args+=(--api-key-file "$api_key_file")
    fi
    case "$family" in
        pipeline)
            entry=$ROOT/runtime/pipeline.py
            args=(pipeline "${modern_args[@]}" "${modern_extra[@]}" -- "${evolution_args[@]}" "${evolution_extra[@]}") ;;
        modern)
            out=$modern_out
            ((preflight)) || [[ $pca == auto && $with_atlas == 0 ]] || fail 'PCA/Atlas options require pipeline or evolve'
            args=("${modern_args[@]}" "${extra[@]}") ;;
        evolution) args=("${evolution_args[@]}" "${extra[@]}") ;;
        final) entry=$ROOT/f/final_cli.py; args=("${final_args[@]}") ;;
    esac
    if ((full_help)); then
        args=("$command" --help)
    fi
    if ((version)); then args=(--version); fi
    local logged=0
    if (( !full_help && !version )); then
        case "$command" in pipeline|run|scan|evolve) logged=1 ;; esac
    fi
    if ((dry_run)); then
        printf 'MemoryMax=%s GiB; MemorySwapMax=0; OOMPolicy=kill; output=%s\n' "$memory" "$out"
        if ((logged)); then
            printf 'Screen progress interval=%ss; verbose=%s; full logs=%s/logs/\n' "$screen_interval" "$verbose_screen" "$out"
        fi
        printf '%q ' "$python" "$entry" "${args[@]}"; printf '\n'
        return 0
    fi
    command -v "$python" >/dev/null 2>&1 || fail "Python not found: $python. Install environment.yml or supply --python. Help: ./allele3.sh -h"
    # Use matching native tools (bcftools/PLINK) from the selected environment.
    python=$(command -v "$python")
    PATH="$(dirname -- "$python"):$PATH"
    export PATH
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export MPLCONFIGDIR=${MPLCONFIGDIR:-"$HOME/.cache/allele3-evolution/matplotlib"}
    # Non-login WSL terminals may omit these despite a running user manager.
    export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
    export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
    command -v systemd-run >/dev/null 2>&1 || fail 'systemd-run is required for the process-tree memory cap; analysis was not started'
    [[ -r /sys/fs/cgroup/cgroup.controllers ]] || fail 'cgroup v2 is required; analysis was not started'
    timeout -k 1s 5s systemctl --user show-environment >/dev/null 2>&1 || fail 'User systemd is unavailable; analysis was not started without a memory cap'
    local bytes=$((memory * 1024 * 1024 * 1024))
    local unit
    unit="allele3-$(date +%s)-$$-$RANDOM"
    printf '[allele3] command=%s output=%s RAM=%s GiB swap=0\n' "$command" "$out" "$memory" >&2
    local -a launch=("$python" "$entry" "${args[@]}")
    if ((logged)); then
        local -a console=("$python" "$ROOT/runtime/console_log.py" --log-dir "$out/logs"
            --label "$command" --interval "$screen_interval")
        ((verbose_screen == 0)) || console+=(--verbose)
        launch=("${console[@]}" -- "${launch[@]}")
    fi
    # The kernel charges the parent AND all descendants to the same cgroup.
    # The inner guard verifies the effective limits before launching Python.
    exec systemd-run --user --scope --quiet --unit="$unit" \
        -p "MemoryMax=$bytes" -p MemorySwapMax=0 -p OOMPolicy=kill -p TasksMax=256 \
        bash "$ROOT/runtime/limited_exec.sh" "$bytes" "${launch[@]}"
}

main "$@"

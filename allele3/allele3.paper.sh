#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
allele3: region enrichment and within-site allele tests.

Usage:
  ./allele3.paper.sh run|enrich_test|allele_test [options]

Modules:
  run                             Launch enrich_test and allele_test in separate directories.
  enrich_test                     Region enrichment against the selected reference.
  allele_test                     Within-site allele tests using existing Atlas VCFs.

Options:
  --grch 38                       Target build [38]; uses the shared build-specific root.
  --out DIR                       Output root/directory for the selected workflow.
  --replace TRUE|FALSE            Recompute / reuse completed results [FALSE].
  --python FILE                   Scientific interpreter [ALLELE_PYTHON].
  --job-memory-gb N               Allocation RAM [run: HPC 160 / local 32 GiB].
  --memory-limit-gb N             Scientific module RAM (individual modules only).
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    LSF queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --job-hours N                   Job duration limit.
  --job-dir DIR                   Generated command/job directory.
  --write-cmd                     Write a reviewable command file.
  --wait TRUE|FALSE               Wait for submitted work.
  --plan-only                     Write plan without executing/submitting.
  --help-full                     Native full help for the selected scientific module.
  --target gnomad|allofus|custom  Target collection [gnomad].
  --ref alphagenome|herv|archaic  Reference for enrichment [alphagenome].
  --scan-dir DIR                  Existing target scan/summary directory.
  --atlas-dir DIR                 Existing AlphaGenome Atlas resources.
  --avi-dir DIR                   Existing AVI resources.
  --herv-scan-dir DIR             Existing HERV scan reference.
  --herv-bed FILE --stable-herv-bed FILEHERV interval overrides.
  --erv-joined FILE               Joined ERV annotations.
  --archaic-dir DIR               Existing IBDmix segment reference.
  --archaic-population NAME       Archaic population [ALL].
  --archaic-reference NAME        Archaic reference [all].
  --archaic-min-carriers N        Minimum archaic carriers [2].
  --archaic-min-slod N            Minimum segment LOD [0].
  --archaic-min-bp N              Minimum segment length [0].
  --chr LIST                      Chromosomes [autosomes].
  --population NAME               Target population [ALL].
  --fasta FILE --no-context       Sequence context source / disable context matching.
  --exclude-bed FILE              Excluded intervals; repeatable.
  --outcomes CSV                  Tested outcomes (see --help-full for current defaults).
  --threshold-kind auto|AC|AF     Eligibility threshold [auto].
  --min-ac N --min-an N --min-af NEligibility thresholds [5, 100, 0.0001].
  --window-bp N --block-bp N      Analysis window/block sizes [10000, 1000000].
  --min-background-sites N        Minimum background sites [100].
  --region-windows CSV            Enrichment windows [10000,100000,1000000 bp].
  --highlight-fold N              Export O/E candidate threshold [3].
  --bootstrap N --seed N          Bootstrap/seed [1000, 20261115].
  --profile-sites N               Profile site budget [10000].
  --controls-per-stratum N        Controls per stratum [100].
  --controls-per-case N           Controls per case [3].
  --radii CSV                     Context radii [1000,10000,100000,1000000].
  --native-binary FILE            Native worker override.
  --jobs N                        Enrichment chromosome workers [0: allocated cores].
  -h, --help                      Show this help.

Examples:
  ./allele3.paper.sh run --target gnomad --grch 38 --replace FALSE --plan-only
  ./allele3.paper.sh enrich_test --target gnomad --ref alphagenome --grch 38 --plan-only
  ./allele3.paper.sh allele_test --target gnomad --grch 38 --submit hpc --queue 33 --cores 40
  ./allele3.paper.sh enrich_test --help-full
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
case "$1" in
    run|enrich_test|allele_test) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ $1 == run ]]; then
    shift
    for module in enrich_test allele_test; do
        args=("$@")
        for ((i=0; i<${#args[@]}; i++)); do
            case ${args[i]} in
                --out|--job-dir) args[i+1]=${args[i+1]:?missing directory}/$module; ((i+=1)) ;;
                --out=*|--job-dir=*) args[i]=${args[i]}/$module ;;
            esac
        done
        bash "$ROOT/allele3.paper.sh" "$module" "${args[@]}"
    done
    exit 0
fi
WORKFLOW_NAME=allele3
source "$ROOT/f/0.run.sh"
research_main() {
    local module=$1; shift
    local target=gnomad reference=alphagenome out= python=${ALLELE_PYTHON:-}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local -a args=()
    while (($#)); do
        case "$1" in
            --target) need_value "$@";target=$2;args+=("$1" "$2");shift 2;;
            --ref) need_value "$@";reference=$2;args+=("$1" "$2");shift 2;;
            --out) need_value "$@";out=$2;shift 2;;
            --python) need_value "$@";python=$2;shift 2;;
            --memory-limit-gb) need_value "$@";memory=$2;shift 2;;
            --screen-log-interval) need_value "$@";interval=$2;shift 2;;
            --help-full|-h|--help) exec "${python:-python3}" "$ROOT/f/research_tests.py" --help;;
            *) args+=("$1");shift;;
        esac
    done
    if [[ -z $out ]]; then
        out="$ALLELE_ANALYSIS_ROOT/$module/$target"
        [[ $module != enrich_test ]] || out="$ALLELE_ANALYSIS_ROOT/$module/${reference}_${target}"
    fi
    python=$(workflow_python "$python" "$HOME/anaconda3/envs/allele3-evolution/bin/python")
    export WORKFLOW_LOG_DIR="$out/logs"
    workflow_launch "$module" "$python" "$ROOT/f/research_tests.py" "$out" "$memory" "$interval" 0 0 1 \
        "$module" --out "$out" --replace "${ALLELE_REPLACE:-FALSE}" "${args[@]}"
}

research_main "$@"

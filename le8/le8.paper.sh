#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
LE8: final reports, portable exports and Shiny preparation.

Usage:
  ./le8.paper.sh run|MODULE[,MODULE...] [options]

Modules:
  run                             final -> share -> shiny for selected scopes.
  final                           Aggregate completed analysis results.
  share                           Export a portable viewer bundle.
  shiny                           Start local viewer; prepare viewer inputs on HPC.

Options:
  --Y CSV                         Outcomes [Y or cvd_cad,ra].
  --biom CSV                      Molecular layers: prot,met [BIOM or prot,met].
  --analysis-root DIR             Results root [LE8_ANALYSIS_ROOT or ANALYSIS_ROOT/le8].
  --replace TRUE|FALSE            Recompute / reuse complete results [TRUE].
  --config FILE                   C1 method configuration [config/c1_methods.json].
  --memory-limit-gb N             Per-job host RAM budget [LE8_MEMORY_LIMIT_GB or 32 GiB].
  --hours N                       Job time limit [LE8_HOURS or 72 hours].
  --python FILE                   Scientific Python [shared LE8 runtime].
  --r-bin FILE                    Rscript; HPC workers use shared 0f/0os.r.sh.
  --plan-dir DIR                  Plan root [analysis-root/logs/workflows].
  --seed N                        Scientific random seed [2026].
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --jobs-N N                      Shared workflow concurrency setting.
  --wait TRUE|FALSE               Wait for submitted workflow [FALSE].
  --plan-only                     Write commands/plan without executing/submitting.
  --fit-joint                     Enable optional final-stage joint fitting.
  --fit-genetic                   Enable optional final-stage genetic fitting.
  --details                       Run detailed native final-stage processing.
  --index-only                    Request final-stage index-only processing.
  --host ADDRESS                  Local viewer bind address [127.0.0.1].
  --port N                        Local viewer port [3839].
  -h, --help                      Show this help.

Examples:
  ./le8.paper.sh run --Y cvd_cad,ra --biom prot,met --plan-only
  ./le8.paper.sh final --Y cvd_cad --biom prot --submit hpc --queue 33 --cores 40
  ./le8.paper.sh final,share --Y cvd_cad --biom prot --analysis-root /work/sph-huangj/analysis/le8 --plan-only
  ./le8.paper.sh shiny --Y cvd_cad --biom prot --submit local --host 127.0.0.1 --port 3839
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run ]]; then shift; set -- final,share,shiny "$@"; fi
if [[ $1 == *, ]]; then usage >&2; exit 2; fi
IFS=, read -ra selected_modules <<< "$1"
for module in "${selected_modules[@]}"; do
    case "$module" in
        final|share|shiny) ;;
        *) printf 'Unknown module: %s\n' "$module" >&2; usage >&2; exit 2 ;;
    esac
done
LE8_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$LE8_ROOT/f/entry.sh"
set -- "${WORKFLOW_ARGS[@]}"
source "$LE8_ROOT/../0f/console.sh"
{
    le8_main "$@"
    exit $?
}

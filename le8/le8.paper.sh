#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd -P)/paper_name.sh"
workflow_paper_name "${BASH_SOURCE[0]}" "$@"
set -euo pipefail

usage() {
    cat <<'HELP'
LE8: paper reports, portable exports and Shiny preparation.

Usage:
  ./le8.paper.sh run|MODULE[,MODULE...] [options]

Modules:
  run | report                    Summarize saved results, led by first-year ICD profiles and individual predictivity.
  paper                           Detailed native paper processing for selected scopes.
  share                           Export a portable viewer bundle.
  shiny                           Start local viewer; prepare viewer inputs on HPC.

Output directory: analysis/le8/paper/{main/results,main/reports,logs,jobs,qc,tmp}/

Results contract: only this code release is accepted (RELEASE.json).
Start a new release in an empty analysis directory; no migration, repair or old-result fallback.

Options:
  --Y CSV                         Outcomes [run/report: discover saved traits; paper: Y or cvd_cad,ra].
  --dx-date CSV                   Dates for current results [icd,any1]. Optional sensitivity CSV.
  --biom CSV                      Molecular layers: prot,met [BIOM or prot,met].
  --analysis-root DIR             Results root [LE8_ANALYSIS_ROOT or ANALYSIS_ROOT/le8].
  --replace TRUE|FALSE            Detailed modules recompute / reuse [TRUE]; run refreshes the summary.
  --config FILE                   C1 method configuration [config/c1_methods.json].
  --python FILE                   Scientific Python [shared LE8 runtime].
  --r-bin FILE                    Rscript; HPC workers use shared 0f/0os.r.sh.
  --plan-dir DIR                  Plan root [analysis-root/paper/jobs].
  --seed N                        Scientific random seed [2026].
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --jobs-N N                      Shared workflow concurrency setting.
  --wait TRUE|FALSE               Wait for submitted workflow [FALSE].
  --plan-only                     Write commands/plan without executing/submitting.
  --fit-joint                     Enable optional paper-stage joint fitting.
  --fit-genetic                   Enable optional paper-stage genetic fitting.
  --details                       Run detailed native paper-stage processing.
  --index-only                    Request paper-stage index-only processing.
  --host ADDRESS                  Local viewer bind address [127.0.0.1].
  --port N                        Local viewer port [3839].
  -h, --help                      Show this help.

Examples:
  ./le8.paper.sh run --Y cvd_cad,cvd_stroke_i,cvd_pad,cvd_afib,cvd_hfail,mi --biom prot,met
  ./le8.paper.sh run --Y cvd_cad,ra --biom prot,met --plan-only
  ./le8.paper.sh paper --Y cvd_cad --biom prot --queue 33 --cores 40
  ./le8.paper.sh paper,share --Y cvd_cad --biom prot --analysis-root /work/sph-huangj/analysis/le8 --plan-only
  ./le8.paper.sh shiny --Y cvd_cad --biom prot --host 127.0.0.1 --port 3839
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run || $1 == report ]]; then
    shift
    # Keep the detailed native fit switches supported by the former run command.
    action=paper_summary
    for arg in "$@"; do
        case "$arg" in --details|--fit-joint|--fit-genetic|--fit-reference|--index-only) action=paper ;; esac
    done
    set -- "$action" "$@"
fi
if [[ $1 == *, ]]; then usage >&2; exit 2; fi
IFS=, read -ra selected_modules <<< "$1"
for module in "${selected_modules[@]}"; do
    case "$module" in
        paper_summary|paper|share|shiny) ;;
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

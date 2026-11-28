#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd -P)/paper_name.sh"
workflow_paper_name "${BASH_SOURCE[0]}" "$@"
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/output_dirs.sh"
script_output_directory_args "${BASH_SOURCE[0]}" "$@" || { return 2 2>/dev/null || exit 2; }
set -- "${SCRIPT_OUTPUT_DIRECTORY_ARGS[@]}"
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
case "${1:-}" in
    ""|-h|--help|help)
        cat <<'HELP'
GU reports from completed chromosome exports.

Usage:
  ./gu.paper.sh ibdmix --batches ygb001-049
  ./gu.paper.sh ibdmix --target NAME --dir-ibdmix RAW [options]
  ./gu.paper.sh export --report DIR/asv --matrix --out ASV.h5
  ./gu.paper.sh export --report DIR/asv --locus ASV_ID --out genotypes.tsv.gz

IBDmix batch report:
  --batches LIST                  YGB batches, e.g. ygb001 or ygb001-049 [ygb001].
                                  Includes 1kg and chromosomes 1-22 and X.
  --dir-ibdmix DIR                Root containing 1kg/ and ygbNNN/ [ANALYSIS_ROOT/gu/ibdmix].
  --window-bp N --min-carriers N  ASV window and carrier threshold [50000, 5].
  --callable-bed FILE             Optional shared GRCh37 callable regions.
  --roster-dir DIR                Optional tested rosters: ygbNNN.tsv (sample, x_tested).
  Main report                     ANALYSIS_ROOT/gu/ibdmix/paper/asv/asv.html

IBDmix single-cohort report:
  --target NAME --dir-ibdmix DIR  Report one cohort from its raw/ directory.
  --grch 37|38 --window-bp N      Coordinate build and summary window [37, 100000].
  --sample-panel FILE --x-males   Optional tested roster and male-only X denominator.
  Output                          ANALYSIS_ROOT/gu/ibdmix/paper/NAME.segments.html

Shared options:
  --dir-out DIR                   Report directory [ANALYSIS_ROOT/gu/ibdmix/paper].
                                  Also accepts --dir-out-ibdmix_report.
  --source TEXT                   Record input provenance.
  --submit auto|local|hpc         Select execution environment.
  --plan-only                     Write the job plan without running the analysis.
  --cores N --queue NAME          Shared job CPU budget and queue.
  -h, --help                      Show help.

Examples:
  ./gu.paper.sh ibdmix --batches ygb001-049
  ./gu.paper.sh ibdmix --target 1kg --dir-ibdmix /work/sph-huangj/analysis/gu/ibdmix/1kg/raw
  ./gu.paper.sh export --report /work/sph-huangj/analysis/gu/ibdmix/paper/asv --matrix --out /work/sph-huangj/analysis/gu/ibdmix/paper/ASV.h5
HELP
        exit 0 ;;
    run|ibdmix|asv|export) ;;
    *) echo "Unknown report: $1. See ./gu.paper.sh -h." >&2; exit 2 ;;
esac
module=$1
# The same IBDmix module accepts a batch atlas or a named single-cohort report.
if [[ $module == ibdmix ]]; then
    module=asv
    for arg in "$@"; do
        case "$arg" in
            --target|--target=*) module=run; break ;;
        esac
    done
fi

source "$root/../0f/project_entry.sh"
# The shared planner otherwise places logs under gu/gu.paper. Supply the GU
# default only on submission; an explicit output or job directory takes priority.
if [[ ${WORKFLOW_WORKER:-0} != 1 && ${PROJECT_JOB_ACTIVE:-0} != 1 ]]; then
    default_logs=1
    for arg in "$@"; do
        case "${arg%%=*}" in
            --dir-out|--dir-jobs) default_logs=0; break ;;
        esac
    done
    if (( default_logs )); then
        source "$root/../0f/0os.sh"
        workflow_args "$@"
        set -- "${WORKFLOW_ARGS[@]}" --dir-jobs "$ANALYSIS_ROOT/gu/ibdmix/paper/logs"
    fi
fi
project_workflow_init "$@"
set -- "${WORKFLOW_ARGS[@]}"
if [[ ${1:-} == --help ]]; then
    set -- "$module" --help
fi
shift
source "$root/f/runtime.env"
python=${GU_ASV_PYTHON:-$GU_ENV_PREFIX/bin/python3}
[[ -x $python ]] || python=${SHARED_PYTHON:-python3}
if [[ $module == export ]]; then
    exec "$python" "$root/f/asv_export.py" "$@"
elif [[ $module == asv ]]; then
    exec "$python" "$root/f/asv_scan.py" "$@"
fi
exec "$python" "$root/f/ibdmix_paper.py" "$@"

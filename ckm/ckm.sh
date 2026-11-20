#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
CKM: within-stage prognosis and shared/stage-dependent factors.

Usage:
  ./ckm.sh run|clock|assoc|marker|enrich [options]

Modules:
  run                             clock -> assoc -> marker -> enrich; writes caches/tables.
  clock                           Grouped out-of-fold clinical/prot/met/joint prognosis.
  assoc                           Stage-stratified effects, interactions, lag analyses and Yang.
  marker                          Cross-sectional molecular stage markers, separately labelled.
  enrich                          Optional offline gene-set enrichment with local annotation/GMT.
  Publication                     Use ./ckm.paper.sh to refresh combined ckm.html.

Options:
  --config FILE                   Model/staging configuration [config/ckm.json].
  --out-dir DIR                   Results and ckm.rds [CKM_OUT_DIR or ANALYSIS_ROOT/ckm].
                                  Aliases: --output-dir, --analysis-root.
  --phenotype-dir DIR             UKB root containing Rdata/ and rap/.
  --stage3-method PREVENT|SCORE2  Clinical stage-3 definition.
  --definition aha2026|aha2023    CKM staging definition.
  --domains CSV                   Selected domains, e.g. lifestyle,prot,met.
  --threads N                     R/association workers [CKM_THREADS or WORKFLOW_CORES];
                                  must be positive and no greater than --cores.
  --rscript FILE                  Explicit Rscript interpreter [CKM_RSCRIPT or Rscript].
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --jobs-N N                      Shared workflow concurrency setting.
  --wait TRUE|FALSE               Wait for submitted workflow [FALSE].
  --plan-only                     Write commands/plan without executing/submitting.
  --job-hours N                   Shared job duration limit.
  --dir-jobs DIR                  Shared command/plan directory.
  --dir-dat DIR --dir-big DIR     Shared data-root overrides.
  wait TRUE|FALSE                 Legacy spelling of --wait.
  -h, --help                      Show this help.

Examples:
  ./1.phe.sh --config config/ckm.json --submit local
  # Wait for phenotype completion before analysis.
  ./ckm.sh run --submit hpc --queue 33 --cores 40 --threads 8 --plan-only
  ./ckm.sh run --submit hpc --queue 33 --cores 40 --threads 8
  ./ckm.sh clock --domains prot,met --submit local
  # Wait for analysis completion before the combined report.
  ./ckm.paper.sh report --out-dir /work/sph-huangj/analysis/ckm --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
case "$1" in
    run|clock|assoc|marker|enrich) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
CKM_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$CKM_ROOT/../0f/project_entry.sh"
# Normalize the supported legacy wait spelling before shared argument parsing.
ckm_args=()
while (($#)); do
    case "$1" in
        wait|--wait)
            [[ ${2:-} == TRUE || ${2:-} == FALSE ]] || { echo 'wait requires TRUE or FALSE' >&2; exit 2; }
            ckm_args+=(--wait "$2"); shift 2 ;;
        *) ckm_args+=("$1"); shift ;;
    esac
done
project_workflow_init "${ckm_args[@]}"
set -- "${WORKFLOW_ARGS[@]}"
source "$CKM_ROOT/f/entry.sh"

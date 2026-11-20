#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
Diet: cohort phenotypes, associations and downstream analysis.

Usage:
  ./diet.sh run|MODULE [options]

Modules:
  run                             Cohort phenotypes -> UKB DOSE foodwas/ml -> assoc ->
                                  UKB DOSE mechanism/signature -> bayes -> extensions.
  ukb | nhanes | chns | charls    Build selected cohort phenotypes; not diet scoring.
  assoc                           Cohort association analysis.
  bayes                           Bayesian analysis, then configured extensions.
  foodwas | ml                    UKB DOSE derivation / machine-learning stages.
  mechanism | signature           UKB DOSE molecular follow-up stages.
  extensions                      Run named configured extensions.
  status                          Display recorded stage status.
  Publication                     Use ./diet.paper.sh for final and web.

Options:
  --cohort NAME                   ukb|nhanes|chns|charls|all [all].
  --study dose|maha|all           Study selection [all]; this all remains supported.
  --config FILE                   Pipeline configuration [config/diet.json].
  --dir-output DIR                Results root [configured output directory].
  --dir-files DIR                 Dictionary/resource directory.
  --dir-ukb DIR                   UKB input directory.
  --dir-nhanes DIR                NHANES input directory.
  --dir-nhanes-dose DIR           DOSE NHANES root with cdc/ and FPED/MPED inputs.
  --dir-chns DIR                  CHNS input directory.
  --dir-charls DIR                CHARLS input directory.
  --charls-input FILE             CHARLS input override.
  --charls-phe FILE               CHARLS phenotype override.
  --charls-config FILE            CHARLS configuration override.
  --foods FILE                    Food input override.
  --replace TRUE|FALSE            Recompute / reuse verified outputs [TRUE].
  --rscript FILE                  R interpreter; alias --r-bin.
  --python FILE                   Python interpreter; alias --python-bin.
  --dry-run                       Print scientific commands without running them.
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --jobs-N N                      Shared workflow concurrency setting.
  --wait TRUE|FALSE               Wait for submitted workflow [FALSE].
  --plan-only                     Write commands/plan without executing/submitting.
  --job-hours N                   Shared job duration limit.
  --dir-jobs DIR                  Shared command/plan directory.
  --core                          Omit the four UKB DOSE extensions from run.
  --extension NAME                Select configured extension; repeatable.
  --steps CSV                     Original DOSE steps in registered execution order.
  --list-steps                    List registered steps for the selected module.
  --jobs N                        Compatibility alias for scientific --cores.
  -h, --help                      Show this help.

Examples:
  ./diet.sh run --cohort ukb --study dose --core --plan-only
  # Generate scores, then wait for successful completion.
  ./diet.sh run --cohort ukb --study dose --core --replace FALSE
  ./diet.sh run --cohort ukb --study dose --submit hpc --queue 33 --cores 40 --plan-only
  ./diet.sh foodwas --cohort ukb --study dose --list-steps
  ./diet.sh assoc --cohort nhanes --study dose --dry-run
  # After the analyses finish:
  ./diet.paper.sh run --cohort all --study all --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
case "$1" in
    run|ukb|nhanes|chns|charls|assoc|bayes|foodwas|ml|mechanism|signature|extensions|status) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
diet_entry="$ROOT/diet.sh"
source "$ROOT/f/entry.sh"

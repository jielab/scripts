#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
Diet: final figures, collected results and offline report.

Usage:
  ./diet.paper.sh run|final|web [options]

Modules:
  run                             final -> web (also accepted as final,web).
  final                           Render/collect publication figures from existing results.
  web                             Collect existing results and refresh offline viewer.

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
  --collect-only                  Collect existing final figures without rerendering composites.
  -h, --help                      Show this help.

Examples:
  ./diet.paper.sh run --cohort all --study all --plan-only
  ./diet.paper.sh final --cohort ukb --study dose --collect-only --dry-run
  ./diet.paper.sh web --cohort ukb --study dose --dir-output /work/sph-huangj/analysis/diet --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run ]]; then shift; set -- final,web "$@"; fi
case "$1" in
    final|web|final,web) ;;
    *) printf 'Unknown module: %s\n' "$1" >&2; usage >&2; exit 2 ;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
diet_entry="$ROOT/diet.paper.sh"
source "$ROOT/f/entry.sh"

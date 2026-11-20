#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
FLUX: HealthFlux model training and evaluation.

Usage:
  ./flux.sh run|MODULE [options]

Modules:
  run                              Run the configured analysis stages in order.
  prepare | train | readout        Prepare inputs, train or extract model readouts.
  evaluate | baseline | report     Evaluate models, baselines or reports.
  demo | rollout | ablate          Run the selected experiment.

Options:
  -h, --help                       Show this help.

Examples:
  ./flux.sh run --plan-only
  ./flux.sh demo --machine local --device cpu
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- all "$@"; fi
# HealthFlux independent reproduction. All analysis code lives in f/.
set -euo pipefail
usage() {
bash "${BASH_SOURCE[0]}" --help
}
for arg in "$@"; do
    case "$arg" in -h|--help) usage; exit 0 ;; esac
done
FLUX_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export FLUX_ROOT PYTHONDONTWRITEBYTECODE=1
source "$FLUX_ROOT/f/runtime.sh"

flux_runtime "$@"
"$FLUX_PYTHON" -c 'import yaml' 2>/dev/null || {
    echo "ERROR: PyYAML missing in $FLUX_PYTHON. Run ./install.sh --machine hpc --device gpu first." >&2
    exit 2
}
exec "$FLUX_PYTHON" "$FLUX_ROOT/f/flux.py" "$@"

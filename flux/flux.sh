#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Usage: ./flux.sh [prepare|train|readout|evaluate|baseline|report|all|demo|rollout|ablate|test] [OPTIONS]
  --config FILE                 YAML overrides [config/flux.yml]
  --machine auto|local|hpc       Execution machine [hpc]
  --device auto|cpu|gpu          Processor [gpu]
  --dir-ukb DIR --dir-out DIR    Input and output roots
  --cores N --gpus N             Worker resources
  --resume --replace            Resume training / rebuild selected outputs
  --plan-only                   Write commands without analysis/submission

Examples:
    ./flux.sh --plan-only
    ./flux.sh
    # Optional local CPU demonstration.
    ./flux.sh demo --machine local --device cpu
PUBLIC_ENTRY_HELP
    exit 0
fi
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

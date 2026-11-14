#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
module=$1
shift
if [[ ${1:-} == --help-full ]]; then
    shift
    source "$HERE/f/env.sh"
    exec "$METAG_PYTHON" "$HERE/f/${module}.py" "$@" --help
fi
for arg in "$@"; do
    if [[ $arg == -h || $arg == --help ]]; then
        bash "$HERE/f/usage.sh" "$module"
        exit 0
    fi
done
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/qc_submit.py" "$module" "$@"

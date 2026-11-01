#!/usr/bin/env bash
set -euo pipefail
METAG_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
source "$METAG_CODE_ROOT/f/scratch.sh"
source "$METAG_CODE_ROOT/f/env.sh"
"$METAG_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else "metag HPC needs Python >= 3.10; load a site Python module or set METAG_PYTHON in env.sh")'
exec "$METAG_PYTHON" "$METAG_CODE_ROOT/f/hpc_cli.py" "$@"

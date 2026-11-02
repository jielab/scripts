#!/usr/bin/env bash
set -euo pipefail
METAG_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
source "$METAG_CODE_ROOT/f/env.sh"
"$METAG_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else "metag HPC needs Python >= 3.10; load a site Python module or set METAG_PYTHON in env.sh")'
"$METAG_PYTHON" -c 'import yaml' || { echo 'PyYAML is required in METAG_PYTHON (python -m pip install PyYAML).' >&2; exit 1; }
exec "$METAG_PYTHON" "$METAG_CODE_ROOT/f/parallel_cli.py" "$@"

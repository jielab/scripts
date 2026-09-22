#!/usr/bin/env bash
set -euo pipefail
panome_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PANOME_PYTHON:-}" ]]; then
  panome_python="$PANOME_PYTHON"
elif [[ -x "$panome_dir/.venv/bin/python" ]]; then
  panome_python="$panome_dir/.venv/bin/python"
else
  panome_python="python3"
fi
exec "$panome_python" "$panome_dir/f/panome.py" "$@"

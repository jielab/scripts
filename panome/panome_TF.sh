#!/usr/bin/env bash
set -euo pipefail
panome_tf_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PANOME_PYTHON:-}" ]]; then
  panome_tf_python="$PANOME_PYTHON"
elif [[ -x "$HOME/venvs/panome-v3/bin/python" ]]; then
  panome_tf_python="$HOME/venvs/panome-v3/bin/python"
else
  panome_tf_python=python3
fi
exec "$panome_tf_python" "$panome_tf_dir/f/panome_tf.py" --device cuda --cores 16 "$@"

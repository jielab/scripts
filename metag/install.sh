#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in
    doctor|setup|status|audit) exec bash "$HERE/f/entry.sh" "$@" ;;
esac
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/install_tools.py" "$@"

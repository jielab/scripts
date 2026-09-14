#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
set -euo pipefail


# 🚩 Command-line help
usage(){ cat <<'HELP'
cd /mnt/d/scripts/grid

./disco.sh --trait height
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/environment.sh"
exec bash "$ROOT/f/disco.sh" "$@"

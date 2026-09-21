#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
set -euo pipefail


# 🚩 Command-line help
usage(){ cat <<'HELP'
cd /mnt/d/scripts/grid

# GRID needs native Needle trees and four-population features (reuses prepared UKB inputs).
./grid.sh grid --trait height --chrs 1-22 --jobs 4 --threads 8
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/environment.sh"
module=${1:-help}
shift || true
case "$module" in
  -h|--help|help) usage; exit 0 ;;
  grid) exec bash "$ROOT/f/grid.sh" "$@" ;;
  *) echo "ERROR: unknown GRID module '$module'" >&2; usage >&2; exit 2 ;;
esac

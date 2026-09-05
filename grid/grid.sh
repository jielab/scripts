#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/environment.sh"
module=${1:-help}
shift || true
case "$module" in
  -h|--help|help) exec bash "$ROOT/f/help.sh" ;;
  pca) exec bash "$ROOT/pca.sh" "$@" ;;
  grid|eval) exec bash "$ROOT/f/$module.sh" "$@" ;;
  csx|disco) echo "ERROR: use ./$module.sh $*" >&2; exit 2 ;;
  ancestry) echo 'ERROR: ancestry is included in ./pca.sh (or ./grid.sh pca).' >&2; exit 2 ;;
  *) echo "ERROR: unknown GRID module '$module'" >&2; bash "$ROOT/f/help.sh" >&2; exit 2 ;;
esac

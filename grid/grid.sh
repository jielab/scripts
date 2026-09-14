#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
set -euo pipefail


# 🚩 Command-line help
usage(){ cat <<'HELP'
cd /mnt/d/scripts/grid

./pca.sh --jobs 4 --threads 8
./csx.sh --trait height --chrs 1-22 --jobs 4 --threads 8
./disco.sh --trait height

# GRID needs native Needle trees and four-population features (reuses prepared UKB inputs).
/mnt/d/scripts/gu/arg.sh build --dir-gen /mnt/i/ukbGen/37 --dir-pfile /mnt/i/ukbGen/37/hap \
  --method needle --format native --features TRUE --chr 1-22 --threads 8 --jobs 1
./grid.sh grid --trait height --chrs 1-22 --jobs 4 --threads 8
./grid.sh eval --trait height
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/environment.sh"
module=${1:-help}
shift || true
case "$module" in
  -h|--help|help) usage; exit 0 ;;
  pca) exec bash "$ROOT/pca.sh" "$@" ;;
  grid|eval) exec bash "$ROOT/f/$module.sh" "$@" ;;
  csx|disco) echo "ERROR: use ./$module.sh $*" >&2; exit 2 ;;
  ancestry) echo 'ERROR: ancestry is included in ./pca.sh (or ./grid.sh pca).' >&2; exit 2 ;;
  *) echo "ERROR: unknown GRID module '$module'" >&2; usage >&2; exit 2 ;;
esac

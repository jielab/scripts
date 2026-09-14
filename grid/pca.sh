#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
set -euo pipefail


# 🚩 Command-line help
usage(){ cat <<'HELP'
cd /mnt/d/scripts/grid

./pca.sh --dir-imp /mnt/i/ukbGen/37/imp --jobs 4 --threads 8
HELP
}
case "${1:-}" in -h|--help|help) usage; exit 0;; esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/environment.sh"
exec bash "$ROOT/f/pca.sh" "$@"

#!/usr/bin/env bash
set -euo pipefail
SRC=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
TARGET=${1:-/mnt/d/scripts/grid}
YES=${2:-}
[[ -s "$SRC/grid.sh" && -s "$SRC/csx/PRScsx.py" && -s "$SRC/DiscoDivas/files/g1k_hm3_maf5_woamb_wolr.pca.weight" ]] || { echo 'ERROR: release is missing required bundled assets.' >&2; exit 2; }
[[ $YES == --yes ]] || { echo "This will replace $TARGET after making a complete backup."; echo "Run: $0 '$TARGET' --yes"; exit 2; }
ts=$(date +%Y%m%d-%H%M%S); backup="${TARGET}.backup.${ts}"
if [[ -e $TARGET ]]; then mv "$TARGET" "$backup"; echo "Backup: $backup"; fi
mkdir -p "$(dirname "$TARGET")"
if command -v rsync >/dev/null 2>&1; then rsync -a "$SRC/" "$TARGET/"; else cp -a "$SRC" "$TARGET"; fi
bash -n "$TARGET/grid.sh"; python3 -m py_compile "$TARGET"/f/*.py
cat > "$TARGET/INSTALL_RECORD.txt" <<META
installed=$(date -Is)
source=$SRC
target=$TARGET
backup=$backup
rollback=rm -rf '$TARGET' && mv '$backup' '$TARGET'
META
echo "Installed GRID at $TARGET"
echo "Rollback: rm -rf '$TARGET' && mv '$backup' '$TARGET'"

#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$HERE/config.json"
ROOT_OVERRIDE=""
while (($#)); do
    case $1 in
        -h|--help)
            cat <<'HELP'
Usage: ./install.sh --qc-tools [--config FILE] [--root DIR]
Install minimap2 2.28 (Linux x86_64, SSE2 binary) for metag.sh refine.
Example: ./install.sh --qc-tools --config config.json
Example: submit=local ./install.sh --qc-tools --config config.json#local --root /data/metag-project
Output: tools_dir/minimap2/2.28/minimap2; metag.sh refine uses it by default.
Override downloads with MINIMAP2_URL; override the installer destination with MINIMAP2_DIR.
HELP
            exit 0 ;;
        --config) CONFIG=$2; shift 2 ;;
        --root) ROOT_OVERRIDE=$2; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
source "$HERE/f/env.sh"
TOOL_ROOT=$("$METAG_PYTHON" - "$HERE" "$CONFIG" "$ROOT_OVERRIDE" <<'PY'
import sys
sys.path.insert(0, sys.argv[1] + '/f')
from common import load_config
print(load_config(sys.argv[2], sys.argv[3] or None)['tools_dir'])
PY
)
DEST="${MINIMAP2_DIR:-$TOOL_ROOT/minimap2/2.28}"
mkdir -p "$DEST"
exec 9>"$DEST/.install.lock"
flock -n 9 || { echo 'minimap2 installation is already active' >&2; exit 2; }
if [[ ! -x "$DEST/minimap2" ]]; then
    [[ $(uname -m) == x86_64 && $(uname -s) == Linux ]] || {
        echo 'Install minimap2 2.28 for your platform and pass --minimap2 PATH' >&2; exit 2;
    }
    archive="$DEST/minimap2-2.28_x64-linux.tar.bz2"
    url="${MINIMAP2_URL:-https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2}"
    bundled="${DIR_SOFTWARE}/metag/minimap2-2.28_x64-linux.tar.bz2"
    if [[ -s "$bundled" ]]; then
        cp "$bundled" "$archive.partial"
    else
        curl -fLsS --connect-timeout 30 --retry 4 --retry-delay 5 -C - -o "$archive.partial" "$url"
    fi
    tar -tjf "$archive.partial" >/dev/null
    mv "$archive.partial" "$archive"
    tar -xjf "$archive" -C "$DEST" --strip-components=1
    # The SSE2 build works across ordinary x86_64 compute nodes.
    if [[ -x "$DEST/minimap2.sse2" ]]; then cp "$DEST/minimap2.sse2" "$DEST/minimap2"; fi
fi
"$DEST/minimap2" --version | grep -q '^2\.28'
sha256sum "$DEST/minimap2" > "$DEST/minimap2.sha256"
echo "minimap2 ready: $DEST/minimap2"

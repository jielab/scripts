#!/usr/bin/env bash
set -euo pipefail
# Official SMR 1.4.3 is an AppImage. WSL and batch nodes may not allow FUSE.
export APPIMAGE_EXTRACT_AND_RUN=1
exec "$(dirname -- "${BASH_SOURCE[0]}")/smr.AppImage" "$@"

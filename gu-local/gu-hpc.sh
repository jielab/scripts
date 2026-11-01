#!/usr/bin/env bash
# Compatibility entry point; use ./gu.sh hpc for new commands.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$root/gu.sh" hpc "$@"

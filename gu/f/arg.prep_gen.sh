#!/usr/bin/env bash
# Public input-preparation entry point; the shared parser keeps the two steps compatible.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export REFGEN_GEN4ARG_ONLY=1
exec bash "$HERE/arg.cli.sh" "$@"

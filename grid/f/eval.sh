#!/usr/bin/env bash
# Legacy internal entry: evaluation now lives in Yeval.sh.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
exec bash "$ROOT/Yeval.sh" "$@"

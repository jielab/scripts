#!/usr/bin/env bash
set -euo pipefail
shared=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../0f" && pwd -P)
exec bash "$shared/0os.r.sh" "$@"

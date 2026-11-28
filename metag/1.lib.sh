#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
METAG: reference-library construction and validation.

Usage:
  ./1.lib.sh run|MODULE [options]

Modules:
  run | build                      Build selected group index from a completed Kraken library.
  setup | kraken                   Prepare tools/references, or build full category indexes.
  check | status                   Validate a library or show index availability.
  human | sylph                    Build an optional independent library.

Options:
  --config FILE --root DIR         Library configuration and root.
  --index-use NAME                 optimized for build; also full for check.
  --cores N --threads N            CPU budget and build threads.
  --mem-gib N --hours N            Job memory and duration.
  --max-gib N                      Optimized hash memory cap.
  --estimate auto|reuse|native     Optimized-build size estimation.
  --full                           Hash full union files during check.
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME                     Build queue.
  --plan-only                      Write commands without submitting.
  -h, --help                       Show this help.

Examples:
  ./1.lib.sh run --plan-only
  ./1.lib.sh check
  ./1.lib.sh status
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- build "$@"; fi
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/library_cli.py" "$@"

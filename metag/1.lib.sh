#!/usr/bin/env bash
# One reference library: setup -> kraken -> named group builds -> check/status.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
case "${1:-}" in
    -h|--help|--help-full)
        cat <<'HELP'
Usage: ./1.lib.sh [setup|kraken|build|check|status|human|sylph] [OPTIONS]
Default: build (optimized group only; requires a completed Kraken library).

Stages, in order:
  setup   Install pinned tools and plan references from frozen metadata.
  kraken  Prepare references and build/publish the full category indexes.
  build   Build the selected group index; reuse matching completed output.
  check   Validate the selected library; --full also hashes full union files.
  status  Show legacy, optimized and full index availability.

Current HPC library:
  ./1.lib.sh build
  # Wait for successful completion before checking.
  ./1.lib.sh check
  ./1.lib.sh status

Fresh machine (wait for each HPC build to finish successfully):
  ./1.lib.sh setup --offline
  ./1.lib.sh kraken
  ./1.lib.sh build --index-use legacy
  ./1.lib.sh build --index-use optimized
  ./1.lib.sh check --index-use legacy
  ./1.lib.sh check --index-use optimized
  ./1.lib.sh status
Copy frozen metadata and pinned software bundles as well as raw inputs.
See REBUILD.txt for paths, requirements and verification limits.

Common: --config FILE --root DIR --submit hpc|local|auto --plan-only
Build:  --index-use optimized|legacy [optimized] --queue NAME
        --cores N --threads N --mem-gib N --hours N
        --max-gib N --estimate auto|reuse|native (optimized only)
        --cores defaults to 40 on HPC, 16 locally; optimized build threads
        follow --cores unless --threads is explicit. Legacy requires 24 threads.
Check:  --index-use optimized|legacy|full [optimized] --full
Existing legacy is preserved; optimized may be replaced when settings change.
HPC submission is asynchronous; local execution waits for completion.

Optional independent libraries (not required for the main lib/):
  ./1.lib.sh human
  ./1.lib.sh sylph --buckets prok,virus,euk
HELP
        exit 0 ;;
esac
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/library_cli.py" "$@"

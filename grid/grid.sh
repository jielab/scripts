#!/usr/bin/env bash

# Public entry for GRID twin analysis, benchmark dispatch, and the viewer.
# Copy-ready workflow and common paths are kept in this entry script.
set -euo pipefail
umask 077
usage() {
	cat <<'HELP'
Modules:
    prsformer   PRSformer benchmark
    prsformer2  PRSformer with separate prediction calibration
    abm         Twin-based ABM analysis (Clinical, CSx-linear, MoE, BlockNet)
    workflow    prsformer -> prsformer2 -> abm -> final (all is equivalent)
    final       Final comparison and result index
    cleanup     Preview unused large caches; add --apply to release
    shiny       Result viewer

Usage:
    ./grid.sh [abm] --traits height,ldl,t2dm
    ./grid.sh module1,module2 --traits height,ldl
    Default module: abm (seven independent ABM phases)
    Comma-separated modules run in order with shared options; stop on failure.

Options:
    --traits CSV          [height,ldl,t2dm]
    --module-only         Accepted compatibility option; abm is already standalone
    --predictor NAME      [blocknet] or moe
    --representation NAME [csx_ldblocks_fixed] or csx_scores
    --device NAME         [auto], cpu or cuda
    --config FILE_OR_DIR  [config/abm.json]
    --analysis-root DIR   [/mnt/d/analysis/grid]
    --threads N           [4]
    --dry-run             Preview without running
    <module> --help       Module options
    -h, --help            Help

Example:
    # Run both benchmarks sequentially; stop if the first fails:
    ./grid.sh prsformer,prsformer2 --traits height,ldl
    # Run benchmarks, ABM and the final summary:
    ./grid.sh workflow --traits height,ldl,t2dm
    ./grid.sh --traits height,ldl,t2dm
    ./grid.sh shiny
HELP
}

case "${1:-abm}" in
	-h | --help | help) usage; exit 0 ;;
esac

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# One foreground logger covers the dispatcher and all children, as in le8.sh.
# Screen: stage/progress summaries. Full output: SCRIPT_LOG_DIR.
# SCRIPT_VERBOSE=1 restores raw terminal output for troubleshooting.
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-/mnt/d/analysis/grid/logs}
source "$ROOT/../0f/console.sh"
# Runtime paths (override with GRID_CONDA_ENV or GRID_PYTHON).
# Analysis paths: /mnt/d/analysis/grid; config: $ROOT/config/abm.json.
# Genotypes: /mnt/f/gen/ukb/37/hap; LD blocks: /mnt/f/LD/block.
# Module CLI options override their own defaults; see 0.README.md.
runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then
	default_python="$runtime/bin/python3"
else
	default_python=python3
fi
PYTHON=${GRID_PYTHON:-$default_python}
PYTHON_PATH=$(command -v -- "$PYTHON")
export GRID_PYTHON="$PYTHON_PATH"
export PATH="$(dirname -- "$PYTHON_PATH"):$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
exec "$PYTHON_PATH" "$ROOT/f/twin.py" "$@"

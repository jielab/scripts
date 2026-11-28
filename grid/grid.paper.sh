#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
case "${1:-}" in
    ""|--help|-h|help)
        cat <<'HELP'
./grid.paper.sh run --grch 38 --traits height,ldl,t2dm
  Evaluate available predictions and refresh the dashboard and manuscript.
  Existing run options include --analysis-root, --ldl-threshold, --bootstrap,
  --cores, --submit, --wait and --plan-only. Use run --help for all options.

./grid.paper.sh render --results /path/to/grid/paper/main/results
  Redraw Table 1, Figures 1–5 and the supplement from published aggregate
  TSVs, without a scheduler, genotypes, individual outcomes or model fitting.
  Default output: manuscript/ beside results/. Use render --help for options.
HELP
        exit 0 ;;
    run) exec bash "$root/f/run_entry.sh" paper "$@" ;;
    render)
        shift
        if [[ -n ${GRID_PYTHON:-} ]]; then python=$GRID_PYTHON
        elif [[ -x ${DIR_WORK:-/work/${USER:-}}/.conda/envs/grid/bin/python ]]; then python=${DIR_WORK:-/work/${USER:-}}/.conda/envs/grid/bin/python
        elif [[ -x ${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}/bin/python ]]; then python=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}/bin/python
        else python=python3; fi
        exec "$python" "$root/f/paper_figures.py" "$@" ;;
    *) echo 'Usage: ./grid.paper.sh {run|render} [options]' >&2; exit 2 ;;
esac

#!/usr/bin/env bash
set -euo pipefail
HOME_ARG=${1:-${ARG_NEEDLE_HOME:-}}
[[ -n $HOME_ARG ]] || { echo 'Usage: test_argneedle_install.sh /path/to/arg-needle-scripts' >&2; exit 2; }
GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid-argneedle}
if [[ -x "$GRID_CONDA_ENV/bin/python3" && ${CONDA_PREFIX:-} != "$GRID_CONDA_ENV" ]]; then
  export PATH="$GRID_CONDA_ENV/bin:$PATH"
fi
export PYTHONNOUSERSITE=1
export PYTHONPATH="$HOME_ARG${PYTHONPATH:+:$PYTHONPATH}"
python3 - <<'PY'
import arg_needle,arg_needle_lib,tskit
print('ARG-Needle and tskit imports: PASS')
PY
s=$(find "$HOME_ARG" -type f -path '*/arg_needle/*infer_args_advanced.py' -print -quit)
[[ -s $s ]] || { echo 'infer_args_advanced.py not found' >&2; exit 1; }
help=$(python3 "$s" --help)
for o in --hap_gz --map --out --mode --step --chromosome --num_snp_samples; do grep -Fq -- "$o" <<<"$help" || { echo "missing upstream option $o" >&2; exit 1; }; done
echo 'ARG-Needle advanced CLI: PASS'

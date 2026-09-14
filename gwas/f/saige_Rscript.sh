#!/usr/bin/env bash
# Keep the conda R ABI separate from the user's global R library/profile.
set -euo pipefail
saige_env=${SAIGE_ENV_DIR:-$HOME/anaconda3/envs/saige}
if [[ ! -x "$saige_env/bin/Rscript" ]]; then
  saige_env=$HOME/miniforge3/envs/saige
fi
[[ -x "$saige_env/bin/Rscript" ]] || { echo 'SAIGE R environment not found' >&2; exit 1; }
export R_LIBS="$saige_env/lib/R/library"
export R_LIBS_USER="$saige_env/lib/R/library"
export R_LIBS_SITE="$saige_env/lib/R/library"
exec "$saige_env/bin/Rscript" --vanilla "$@"

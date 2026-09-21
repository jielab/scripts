#!/usr/bin/env bash
set -euo pipefail
args=(--trait "$GRID_TRAIT" --mode "$csx_model" --gwas-dir "$GRID_GWAS_DIR"
 --target-dir "$GRID_TARGET_DIR" --snpinfo "$GRID_CSX_SNPINFO" --ref-dir "$GRID_CSX_REF_DIR"
 --bim "$GRID_CSX_BIM_PREFIX" --work "$work" --score-home "$GRID_SCORE_DIR/$GRID_TRAIT${suffix:+/${suffix#.}}"
 --chrs "${CHRS[*]}" --jobs "$GRID_JOBS" --threads "$GRID_THREADS" --iterations "$GRID_MCMC_ITER"
 --burnin "$GRID_MCMC_BURNIN" --thin "$GRID_MCMC_THIN" --seed "$GRID_SEED" --phi "$GRID_PHI"
 --n-gwas "$GRID_N_GWAS" --remove "$GRID_REMOVE" --keep "$GRID_KEEP" --stage "$GRID_STAGE" --replace "$GRID_REPLACE")
[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || args+=(--check)
printf 'python3 %q ' "$ROOT/f/csx_combined.py" >> "$GRID_COMMAND_FILE"
printf '%q ' "${args[@]}" >> "$GRID_COMMAND_FILE"; printf '\n' >> "$GRID_COMMAND_FILE"
python3 "$ROOT/f/csx_combined.py" "${args[@]}"

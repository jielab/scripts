#!/usr/bin/env bash
set -uo pipefail
[[ -n ${LSB_JOBID:-} ]] || exit 2
cd /work/sph-huangj/scripts/grid
source ../0f/workflow_env.sh
workflow_init grid
AUDIT=/work/sph-huangj/scripts/grid/.preflight/20261105
export WORKFLOW_WORKER=1 WORKFLOW_RESOURCE_MANAGER=lsf WORKFLOW_MEMORY_GIB=16 GRID_MEMORY_CAP_GB=16 GRID_THREADS=2 GRID_JOBS=1 GRID_CACHE_ROOT="$AUDIT/cache"
export PATH="$GRID_CONDA_ENV/bin:$PATH" PYTHONPATH="$PWD/f:$PWD/../0f" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
hostname
for script in 1.pca.sh 2.csx.sh 3.disco.sh; do
 args=(--check --output-root "$AUDIT/native")
 if [[ $script == 2.csx.sh ]]; then args+=(--bim-prefix "$AUDIT/ukb_array"); fi
 bash "$script" "${args[@]}" >"$AUDIT/$script.check.log" 2>&1
 printf '%s\t%s\n' "$script" "$?" >> "$AUDIT/native-exits.tsv"
done
python3 f/prsformer.py prepare --check --cache-dir "$AUDIT/prsformer-cache" --dir-gen "$GRID_TARGET_DIR" --pheno-file "$GRID_PHE_FILE" --ancestry-file "$(dirname "$GRID_PCA_FILE")/ukb.ancestry.auto.tsv.gz" --remove none > "$AUDIT/shared-cohort.check.log" 2>&1
printf 'shared-cohort\t%s\n' "$?" >> "$AUDIT/native-exits.tsv"
cat "$AUDIT/native-exits.tsv"

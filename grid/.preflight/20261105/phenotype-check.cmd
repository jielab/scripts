#!/usr/bin/env bash
set -uo pipefail
[[ -n ${LSB_JOBID:-} ]] || exit 2
cd /work/sph-huangj/scripts/grid
source ../0f/workflow_env.sh
workflow_init grid
AUDIT=/work/sph-huangj/scripts/grid/.preflight/20261105
export PATH="$GRID_CONDA_ENV/bin:$PATH" PYTHONPATH="$PWD/f:$PWD/../0f" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla - <<'R' > "$AUDIT/phenotype-columns.log" 2>&1
for (file in c("phe.rds", "all.rds")) {
 d <- readRDS(file.path("/work/sph-huangj/phe/ukb/Rdata",file))
 cat(file, "rows", nrow(d), "columns", ncol(d), "\n")
 wanted <- c("eid","age","sex","height","ldl","ethnicity", "PC1","PC2","t2dm.Yr2e","t2dm.Yt2e","cnt_icd10_t2dm","date_death","fod_icd10_t2dm")
 for (column in wanted) cat(column, column %in% names(d), "\n")
 cat("matching columns:",paste(grep("t2dm|ldl|cholesterol|ethnic",names(d),value=TRUE),collapse=", "),"\n")
 rm(d); gc()
}
R
python3 f/prsformer.py prepare --check --cache-dir "$AUDIT/prsformer-cache" --dir-gen "$GRID_TARGET_DIR" --pheno-file "$DIR_UKB/Rdata/all.rds" --ancestry-file "$(dirname "$GRID_PCA_FILE")/ukb.ancestry.auto.tsv.gz" --remove none > "$AUDIT/shared-cohort-all.check.log" 2>&1
printf 'shared-cohort-all\t%s\n' "$?" > "$AUDIT/phenotype-exits.tsv"
cat "$AUDIT/phenotype-exits.tsv"

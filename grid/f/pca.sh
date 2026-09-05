#!/usr/bin/env bash
# Cohort preparation: reference-PC projection -> distances/QC -> ancestry.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
for n in "$GRID_COV_PCS" "$GRID_DISTANCE_PCS"; do
  [[ $n =~ ^[1-9][0-9]*$ ]] || _grid_die 'PC counts must be positive integers'
done
((GRID_DISTANCE_PCS >= 5 && GRID_DISTANCE_PCS <= GRID_COV_PCS)) || _grid_die 'Require 5 <= distance PCs <= covariate PCs'
pca=$GRID_PCA_FILE
qc=$(dirname -- "$pca")
raw="$GRID_OUTPUT_ROOT/pca"
mkdir -p "$qc" "$raw/log"
need(){ [[ -s $1 ]] || _grid_die "Missing/empty file: $1"; }
log(){ echo "[PCA] $*" >&2; }
# The projection cache is independent of distances and ancestry outputs.
projected=FALSE
if [[ $GRID_REPLACE == FALSE ]] && python3 "$ROOT/f/pca_cache.py" "$pca" "$GRID_COV_PCS"; then
  log "SKIP projection: existing PC1-PC$GRID_COV_PCS in $pca"
else
  need "$GRID_PCA_WEIGHT"
  available_mb=$(awk '/^MemAvailable:/ {print int($2/1024)}' /proc/meminfo)
  pca_jobs=$((available_mb * 80 / 100 / 4096))
  ((pca_jobs >= 1)) || pca_jobs=1
  ((pca_jobs <= 4)) || pca_jobs=4
  ((pca_jobs <= GRID_JOBS)) || pca_jobs=$GRID_JOBS
  project_chr(){
    local c=$1 prefix="$GRID_IMP_DIR/chr$1" dest="$raw/chr$1" marker="$raw/chr$1.pca${GRID_COV_PCS}.done"
    if [[ $GRID_REPLACE == FALSE && -s $dest.sscore && -s $dest.sscore.vars && -e $marker ]]; then
      log "SKIP projection chr$c: cached chromosome scores"
      return
    fi
    need "$prefix.pgen"; need "$prefix.psam"
    local -a input=(--pfile "$prefix")
    if [[ -s $prefix.pvar.zst ]]; then input+=(vzs); else need "$prefix.pvar"; fi
    grid_run_logged "$raw/log/chr$c.log" plink2 "${input[@]}" \
      --score "$GRID_PCA_WEIGHT" 2 6 header-read no-mean-imputation list-variants cols=nallele,scoresums \
      --score-col-nums "7-$((6 + GRID_COV_PCS))" --rm-dup exclude-mismatch \
      --memory 4096 --threads 1 --out "$dest"
    if [[ $GRID_DRY_RUN == FALSE ]]; then need "$dest.sscore"; need "$dest.sscore.vars"; touch "$marker"; fi
  }
  running=0; status=0
  for c in {1..22}; do
    project_chr "$c" & ((++running))
    if ((running >= pca_jobs)); then wait -n || status=1; running=$((running-1)); fi
  done
  while ((running)); do wait -n || status=1; running=$((running-1)); done
  ((status == 0)) || _grid_die 'PCA projection failed'
  # Persist the projection before distance/QC work, so that either can resume alone.
  grid_run_logged "$raw/log/combine.log" Rscript "$ROOT/f/combine_disco_pca.R" \
    --dir "$raw" --out "$pca" --outdir "$qc" --n-pc "$GRID_COV_PCS" \
    --distance-pcs "$GRID_DISTANCE_PCS" --projection-only TRUE
  projected=TRUE
fi

distance_updated=FALSE
if [[ $GRID_REPLACE == TRUE || $projected == TRUE || ! -s $qc/ukb_reference_distances.tsv.gz || ! -s $qc/pca_qc.xlsx || ! -s $qc/Fig1.PCA_QC.png || $qc/ukb_reference_distances.tsv.gz -ot $pca || $qc/pca_qc.xlsx -ot $pca || $qc/Fig1.PCA_QC.png -ot $pca ]]; then
  need "$GRID_MED_FILE"
  grid_run_logged "$raw/log/distance.log" Rscript "$ROOT/f/combine_disco_pca.R" \
    --pca "$pca" --med "$GRID_MED_FILE" --outdir "$qc" \
    --n-pc "$GRID_COV_PCS" --distance-pcs "$GRID_DISTANCE_PCS"
  distance_updated=TRUE
else
  log 'SKIP reference distances/QC: outputs exist'
fi

if [[ $GRID_REPLACE == TRUE || $distance_updated == TRUE || ! -s $GRID_ANCESTRY_FILE || ! -s $qc/ancestry_auto_qc.xlsx || ! -s $qc/Fig2.Ancestry_QC.png || $GRID_ANCESTRY_FILE -ot $qc/ukb_reference_distances.tsv.gz || $qc/ancestry_auto_qc.xlsx -ot $pca || $qc/Fig2.Ancestry_QC.png -ot $pca ]]; then
  need "$GRID_MED_FILE"; need "$GRID_PHE_FILE"
  grid_run_logged "$raw/log/ancestry.log" Rscript "$ROOT/f/prepare_ancestry_auto.R" \
    --phe "$GRID_PHE_FILE" --pca "$pca" --med "$GRID_MED_FILE" \
    --out "$GRID_ANCESTRY_FILE" --outdir "$qc" --ethnicity-col "$GRID_ETHNICITY_COL" \
    --n-pc "$GRID_COV_PCS" --distance-pcs "$GRID_DISTANCE_PCS" \
    --prob-threshold "$GRID_ANCESTRY_PROB_MIN" --anchor-max-per-group "$GRID_ANCHOR_MAX_PER_GROUP"
else
  log 'SKIP ancestry: outputs exist'
fi
log "Completed: $pca; $qc/ukb_reference_distances.tsv.gz; $GRID_ANCESTRY_FILE"

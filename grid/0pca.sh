#!/usr/bin/env bash
# UKB preparation: reference PCA projection -> distance QC -> ancestry assignment.
# Projection and ancestry: f/0.pca.R; cache and discovery centres: f/0.pca.py.
set -euo pipefail
usage() {
	cat <<'HELP'
PCA projection — project UKB into the reference space used by DiscoDivas

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./0pca.sh --check
  ./0pca.sh --dir-imp /mnt/f/gen/ukb/37/imp --jobs 4
  ./0pca.sh --replace TRUE --jobs 4

Modules (run in order; existing completed outputs are reused):
  projection  Score chromosomes 1-22 with reference PCA loadings; sum PC scores.
  distance    Reference distances, PCA QC workbook and figure.
  ancestry    Phenotype-assisted ancestry assignment and QC (separate from Disco).

Important parameters:
  --dir-imp DIR                 /mnt/f/gen/ukb/37/imp (chr*.pgen/pvar/psam)
  --pca-file FILE               /mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz
  --pca-weight FILE             /mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight
  --med-file FILE               /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --cov-pcs N / --distance-pcs N 20 / 10; require 5 <= distance <= covariate PCs.
  --phe-file FILE               /mnt/d/data/ukb/phe/Rdata/phe.rds
  --ancestry-file FILE          Default: ukb.ancestry.auto.tsv.gz beside PCA file.
  --ancestry-prob-min VALUE     0.90; --anchor-max-per-group 10000
  --jobs N                     At most 4 concurrent chromosomes, memory-limited.
  --threads N                  Accepted shared option; projection uses 1 thread/job.
  --output-root DIR             /mnt/d/analysis/grid (temporary).
  --replace TRUE|FALSE          Default FALSE. TRUE rebuilds projection and QC.
  --check                       Validate existing projection or genotype inputs.
  --dry-run TRUE                Same preflight/plan behavior; no analysis.

Permanent PCA, distance, ancestry and QC files: /mnt/d/data/ukb/pca_proj/
Temporary chromosome scores/commands/logs: /mnt/d/analysis/grid/pca/
Existing PCA is reused independently of temporary files. If changing reference
loadings/genotypes, explicitly use --replace TRUE or a different --pca-file.
HELP
}
case "${1:-}" in -h | --help | help)
	usage
	exit 0
	;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# 🚩 pca_run
pca_run() {
	# Cohort preparation: reference-PC projection -> distances/QC -> ancestry.
	set -euo pipefail
	for n in "$GRID_COV_PCS" "$GRID_DISTANCE_PCS"; do
		[[ $n =~ ^[1-9][0-9]*$ ]] || _grid_die 'PC counts must be positive integers'
	done
	((GRID_DISTANCE_PCS >= 5 && GRID_DISTANCE_PCS <= GRID_COV_PCS && GRID_COV_PCS <= 20)) || _grid_die 'Require 5 <= distance PCs <= covariate PCs'
	pca=$GRID_PCA_FILE
	qc=$(dirname -- "$pca")
	raw="$GRID_OUTPUT_ROOT/pca"
	echo "input genotype: $GRID_IMP_DIR/chr{1..22}"
	echo "input PCA weights: $GRID_PCA_WEIGHT"
	echo "output PCA files: $pca"
	echo "output distance/ancestry/QC: $qc"
	# Fingerprint projection inputs; do not reuse a cache solely because a file exists.
	projection_inputs=("$GRID_PCA_WEIGHT" "$ROOT/f/0.pca.R" "$ROOT/0pca.sh")
	for c in {1..22}; do
		projection_inputs+=("$GRID_IMP_DIR/chr$c.pgen" "$GRID_IMP_DIR/chr$c.psam")
		if [[ -s $GRID_IMP_DIR/chr$c.pvar ]]; then projection_inputs+=("$GRID_IMP_DIR/chr$c.pvar"); else projection_inputs+=("$GRID_IMP_DIR/chr$c.pvar.zst"); fi
	done
	projection_sig=$(python3 "$ROOT/f/0.common.py" signature "$GRID_COV_PCS" --files "${projection_inputs[@]}")
	projection_valid() {
		[[ -f $pca.signature && $(cat "$pca.signature") == "$projection_sig" ]] && python3 "$ROOT/f/0.pca.py" cache "$pca" "$GRID_COV_PCS"
	}
	raw="$raw/$projection_sig"

	grid_select_r data.table,ggplot2,openxlsx,patchwork,MASS
	if [[ $GRID_CHECK == TRUE || $GRID_DRY_RUN == TRUE ]]; then
		if projection_valid; then
			echo "CHECK PASS: existing PCA projection contains PC1-PC$GRID_COV_PCS"
		else
			need "$GRID_PCA_WEIGHT"
			for c in {1..22}; do
				need "$GRID_IMP_DIR/chr$c.pgen"
				need "$GRID_IMP_DIR/chr$c.psam"
				[[ -s $GRID_IMP_DIR/chr$c.pvar || -s $GRID_IMP_DIR/chr$c.pvar.zst ]] || _grid_die "Missing chr$c PVAR"
			done
			echo 'CHECK PASS: genotype inputs exist; projection is needed'
		fi
		return 0
	fi
	mkdir -p "$qc" "$raw/log"
	exec {pca_lock}>"$GRID_OUTPUT_ROOT/pca/run.lock"
	flock -n "$pca_lock" || _grid_die 'Another PCA run is active'
	need() { [[ -s $1 ]] || _grid_die "Missing/empty file: $1"; }
	log() { echo "[PCA] $*" >&2; }
	# The projection cache is independent of distances and ancestry outputs.
	projected=FALSE
	if [[ $GRID_REPLACE == FALSE ]] && projection_valid; then
		log "SKIP projection: existing PC1-PC$GRID_COV_PCS in $pca"
	else
		need "$GRID_PCA_WEIGHT"
		available_mb=$(awk '/^MemAvailable:/ {print int($2/1024)}' /proc/meminfo)
		pca_jobs=$((available_mb * 80 / 100 / 4096))
		((pca_jobs >= 1)) || pca_jobs=1
		((pca_jobs <= 4)) || pca_jobs=4
		((pca_jobs <= GRID_JOBS)) || pca_jobs=$GRID_JOBS
		project_chr() {
			local c=$1 prefix="$GRID_IMP_DIR/chr$1" dest="$raw/chr$1" marker="$raw/chr$1.pca${GRID_COV_PCS}.done"
			if [[ $GRID_REPLACE == FALSE && -s $dest.sscore && -s $dest.sscore.vars && -e $marker ]]; then
				log "SKIP projection chr$c: cached chromosome scores"
				return
			fi
			rm -f -- "$marker"
			need "$prefix.pgen"
			need "$prefix.psam"
			local -a input=(--pfile "$prefix")
			if [[ ! -s $prefix.pvar && -s $prefix.pvar.zst ]]; then input+=(vzs); else need "$prefix.pvar"; fi
			grid_run_logged "$raw/log/chr$c.log" plink2 "${input[@]}" \
				--score "$GRID_PCA_WEIGHT" 2 6 header-read no-mean-imputation list-variants cols=nallele,scoresums \
				--score-col-nums "7-$((6 + GRID_COV_PCS))" --rm-dup exclude-mismatch \
				--memory 4096 --threads 1 --out "$dest"
			if [[ $GRID_DRY_RUN == FALSE ]]; then
				need "$dest.sscore"
				need "$dest.sscore.vars"
				touch "$marker"
			fi
		}
		pids=()
		status=0
		for c in {1..22}; do
			project_chr "$c" &
			pids+=("$!")
			if ((${#pids[@]} >= pca_jobs)); then
				for pid in "${pids[@]}"; do wait "$pid" || status=1; done
				pids=()
				((status == 0)) || _grid_die 'PCA chromosome projection failed'
			fi
		done
		for pid in "${pids[@]}"; do wait "$pid" || status=1; done
		((status == 0)) || _grid_die 'PCA projection failed'
		# Persist the projection before distance/QC work, so that either can resume alone.
		grid_run_logged "$raw/log/combine.log" "${GRID_R[@]}" "$ROOT/f/0.pca.R" projection \
			--dir "$raw" --out "$pca" --outdir "$qc" --n-pc "$GRID_COV_PCS" \
			--distance-pcs "$GRID_DISTANCE_PCS" --projection-only TRUE
		printf '%s\n' "$projection_sig" >"$pca.signature"
		projected=TRUE
	fi

	distance_updated=FALSE
	need "$GRID_MED_FILE"
	distance_sig=$(python3 "$ROOT/f/0.common.py" signature "$GRID_COV_PCS" "$GRID_DISTANCE_PCS" --files "$pca" "$GRID_MED_FILE" "$ROOT/f/0.pca.R")
	if [[ $GRID_REPLACE == TRUE || $projected == TRUE || ! -f $qc/distance.signature || $(cat "$qc/distance.signature" 2>/dev/null) != "$distance_sig" || ! -s $qc/ukb_reference_distances.tsv.gz || ! -s $qc/pca_qc.xlsx || ! -s $qc/Fig1.PCA_QC.png || $qc/ukb_reference_distances.tsv.gz -ot $pca || $qc/pca_qc.xlsx -ot $pca || $qc/Fig1.PCA_QC.png -ot $pca ]]; then
		need "$GRID_MED_FILE"
		grid_run_logged "$raw/log/distance.log" "${GRID_R[@]}" "$ROOT/f/0.pca.R" projection \
			--pca "$pca" --med "$GRID_MED_FILE" --outdir "$qc" \
			--n-pc "$GRID_COV_PCS" --distance-pcs "$GRID_DISTANCE_PCS"
		printf '%s\n' "$distance_sig" >"$qc/distance.signature"
		distance_updated=TRUE
	else
		log 'SKIP reference distances/QC: outputs exist'
	fi

	need "$GRID_PHE_FILE"
	ancestry_sig=$(python3 "$ROOT/f/0.common.py" signature "$GRID_COV_PCS" "$GRID_DISTANCE_PCS" "$GRID_ETHNICITY_COL" "$GRID_ANCESTRY_PROB_MIN" "$GRID_ANCHOR_MAX_PER_GROUP" --files "$pca" "$GRID_MED_FILE" "$GRID_PHE_FILE" "$ROOT/f/0.pca.R")
	if [[ $GRID_REPLACE == TRUE || $distance_updated == TRUE || ! -f $qc/ancestry.signature || $(cat "$qc/ancestry.signature" 2>/dev/null) != "$ancestry_sig" || ! -s $GRID_ANCESTRY_FILE || ! -s $qc/ancestry_auto_qc.xlsx || ! -s $qc/Fig2.Ancestry_QC.png || $GRID_ANCESTRY_FILE -ot $qc/ukb_reference_distances.tsv.gz || $qc/ancestry_auto_qc.xlsx -ot $pca || $qc/Fig2.Ancestry_QC.png -ot $pca ]]; then
		need "$GRID_MED_FILE"
		need "$GRID_PHE_FILE"
		grid_run_logged "$raw/log/ancestry.log" "${GRID_R[@]}" "$ROOT/f/0.pca.R" ancestry \
			--phe "$GRID_PHE_FILE" --pca "$pca" --med "$GRID_MED_FILE" \
			--out "$GRID_ANCESTRY_FILE" --outdir "$qc" --ethnicity-col "$GRID_ETHNICITY_COL" \
			--n-pc "$GRID_COV_PCS" --distance-pcs "$GRID_DISTANCE_PCS" \
			--prob-threshold "$GRID_ANCESTRY_PROB_MIN" --anchor-max-per-group "$GRID_ANCHOR_MAX_PER_GROUP"
		printf '%s\n' "$ancestry_sig" >"$qc/ancestry.signature"
	else
		log 'SKIP ancestry: outputs exist'
	fi
	log "Completed: $pca; $qc/ukb_reference_distances.tsv.gz; $GRID_ANCESTRY_FILE"
}

source "$ROOT/f/0.common.sh"
grid_pipeline pca "$@"

#!/usr/bin/env bash
# Official DiscoDivas: align PCA/PRS -> genetic-distance interpolation -> publish.
# Workflow orchestration is here; pinned upstream program: f/2disco.R.
# Reference: https://github.com/YunfengRuan/DiscoDivas
set -euo pipefail
usage() {
	cat <<'HELP'
DiscoDivas — combine the four ancestry-specific CSx scores for each UKB individual

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./0pca.sh --check
  ./1csx.sh --traits height,ldl,t2dm --jobs 4 --threads 1
  ./2disco.sh --traits height,ldl,t2dm --check
  ./2disco.sh --traits height,ldl,t2dm
  ./2disco.sh --trait height --a-list 1,1,1,1 --regress-pca TRUE

Prerequisites:
  Completed reference PCA projection and the permanent merged 1csx.scores.rds file.
  Existing projection in pca_proj can be reused; no need to rerun it per trait.
  Scored samples absent from PCA are filtered out of all four PRS inputs.
  Negative IDs and /mnt/d/files/ukb.exclude.id are excluded.
  Logs report excluded, filtered and retained sample counts.
  DiscoDivas combines individual scores; it does not directly analyze GWAS or
  produce a single universal SNP weight file (weights vary by individual).

Important parameters:
  --trait NAME / --traits LIST   Default: height,ldl,t2dm; traits run sequentially.
  --pca-file FILE               /mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz
  --med-file FILE               /mnt/d/files/DiscoDivas/med.g1000.4pop.tsv
  --distance-pcs N              10; reference and target use the same PCs (>=5).
  --a-list A1,A2,A3,A4          Default 1,1,1,1, in AFR,EAS,EUR,SAS order.
  --regress-pca TRUE|FALSE      Default TRUE (official recommendation).
  --score-dir DIR               /mnt/d/data/ukb/pgs (input and permanent output).
  --output-root DIR             /mnt/d/analysis/grid (result namespace).
  --replace TRUE|FALSE          Default FALSE; reuse matching outputs.
  --check                       Check prerequisite files/packages; no computation.
  --dry-run TRUE                Same preflight/plan behavior.
  --chrs LIST                   Match the chromosome subset used in CSx, if any.

Input scores:  <score-dir>/<trait>/1csx.scores.rds
Output scores: <score-dir>/<trait>/2disco.scores.rds (IID,disco)
Coefficients:  <score-dir>/<trait>/disco.coef.tsv.gz
Commands/logs/intermediates: /tmp/grid-cache/

Uses official DiscoDivas distance-matrix interpolation, not inverse-square mixing.
Defaults use reference centers and equal quality factors, without phenotype tuning.
HELP
}
case "${1:-}" in -h | --help | help)
	usage
	exit 0
	;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# 🚩 disco_run
disco_run() {
	# Official DiscoDivas: population PRS + reference-projected PCs -> individual PRS.
	set -euo pipefail
	trait=$GRID_TRAIT
	score_home="$GRID_SCORE_DIR/$trait${suffix:+/${suffix#.}}"
	io="$ROOT/f/0.common.py"
	need "$GRID_PCA_FILE"
	need "$GRID_MED_FILE"
	echo "input PCA: $GRID_PCA_FILE"
	echo "input reference centers: $GRID_MED_FILE"
	inputs=("$score_home/1csx.scores.rds")
	need "${inputs[0]}"
	score_lock=$(python3 "$ROOT/f/0.common.py" cache-path "$score_home/1csx.scores.rds")
	mkdir -p "$score_lock"
	exec {score_read_lock}>>"$score_lock/write.lock"
	flock -s "$score_read_lock"
	[[ -z $GRID_REMOVE || -f $GRID_REMOVE ]] || _grid_die "Missing withdrawal file: $GRID_REMOVE"
	echo "input score file: ${inputs[0]}"
	echo "output score files: $score_home/2disco.scores.rds; $score_home/2disco.coefficients.rds"
	[[ -z $GRID_REMOVE ]] || inputs+=("$GRID_REMOVE")
	python3 - "$GRID_DISCO_A" "$GRID_DISTANCE_PCS" <<'PY'
import math,sys
x=list(map(float,sys.argv[1].split(',')))
assert len(x)==4 and all(math.isfinite(v) and v>=0 for v in x) and any(x), '--a-list requires four nonnegative values, at least one positive'
assert 5<=int(sys.argv[2])<=20, '--distance-pcs must be 5..20'
PY

	grid_select_r data.table,dplyr,stringr,rio,optparse
	disco_r=("${GRID_R[@]}")
	echo "Disco R runtime: ${disco_r[*]}"
	[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || {
		echo 'CHECK/PLAN complete; no Disco calculation executed'
		return 0
	}
	sig=$(python3 "$io" signature "$GRID_DISCO_A" "$GRID_REGRESS_PCA" "$GRID_DISTANCE_PCS" --files "${inputs[@]}" "$GRID_PCA_FILE" "$GRID_MED_FILE" "$io" "$ROOT/2disco.sh" "$ROOT/f/2disco.R" "$ROOT/f/0.common.py")
	run="$work/$trait/$sig"
	mkdir -p "$run" "$logdir/$trait"
	cache_sig="$work/$trait/disco.signature"
	exec {lock}>"$work/$trait/run.lock"
	flock -n "$lock" || _grid_die "Another Disco run is active for $trait"
	if [[ -s $score_home/2disco.scores.rds && -s $score_home/2disco.coefficients.rds && -s $cache_sig && $(cat "$cache_sig") == "$sig" && $GRID_REPLACE == FALSE ]]; then
		echo "SKIP $trait: matching permanent Disco results"
		return 0
	fi
	grid_run_logged "$logdir/$trait/prepare.log" python3 "$io" disco-inputs "$GRID_PCA_FILE" "$GRID_MED_FILE" "$score_home" "$run" "$GRID_DISTANCE_PCS" "$GRID_REMOVE"
	grep '^Disco ' "$logdir/$trait/prepare.log"
	prs=()
	for p in "${POPS[@]}"; do prs+=("$run/$p.tsv"); done
	grid_run_logged "$logdir/$trait/disco.log" "${disco_r[@]}" "$ROOT/f/2disco.R" -m "$run/centers.tsv" -p "$run/pca.tsv" --prs.list "$(join_comma "${prs[@]}")" -s IID,PRS -A "$GRID_DISCO_A" --regress.PCA "$GRID_REGRESS_PCA" --print.coef TRUE -o "$run/disco"
	grid_run_logged "$logdir/$trait/validate.log" python3 "$io" disco-output "$run/disco.tsv.gz" "$run"
	grid_run python3 "$ROOT/f/0.common.py" publish disco "$run/disco.tsv.gz" "$score_home/2disco.scores.rds" --remove "$GRID_REMOVE"
	grid_run "${disco_r[@]}" "$ROOT/../0f/results.R" import-table "$run/disco.coef.tsv.gz" "$score_home/2disco.coefficients.rds"
	printf '%s\n' "$sig" >"$run/signature"
	publish "$run/signature" "$cache_sig"
	# Permanent provenance stays usable even after deleting the scratch directory.
	{
		printf 'key\tvalue\n'
		printf 'PCA\t%s\ncenters\t%s\nPCs\t%s\nA\t%s\nregress_PCA\t%s\n' "$GRID_PCA_FILE" "$GRID_MED_FILE" "$GRID_DISTANCE_PCS" "$GRID_DISCO_A" "$GRID_REGRESS_PCA"
		for p in "${POPS[@]}"; do printf '%s\t%s\n' "$p" "$score_home/1csx.scores.rds"; done
	} >"$run/manifest.tsv"
	echo "DONE $trait: $score_home/2disco.scores.rds"
}

source "$ROOT/f/0.common.sh"
grid_pipeline disco "$@"

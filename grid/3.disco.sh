#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
GRID DiscoDivas interpolation and held-out evaluation, CPU only.
Modules: all [default], score, evaluate, prepare.
Inputs: --traits CSV --pca-file FILE --dir-score DIR --analysis-root DIR
  --distance-pcs N --a-list A1,A2,A3,A4 --regress-pca TRUE|FALSE
No GPU processor option applies to Disco.
Machine: --submit hpc|local|auto [auto]
Results: --replace TRUE|FALSE [TRUE]; FALSE reuses complete outputs and rebuilds missing outputs.
Planning: --plan-only --jobs-N N --wait TRUE|FALSE
CPU resources: --queue NAME --job-cores N --job-memory-gb N --job-hours N --threads N

Examples:
  ./3.disco.sh all --traits height,ldl,t2dm
PUBLIC_ENTRY_HELP
    exit 0
fi

# Raw ancestry-score/reference-center interpolation; source identity:
# csx-disco-legacy/raw-interpolation. Pinned score interpolation followed by a same-roster outcome adaptation.
# Workflow orchestration is here; pinned upstream program: f/3.disco.R.
# Reference: https://github.com/YunfengRuan/DiscoDivas
set -euo pipefail
usage() {
	bash "${BASH_SOURCE[0]}" --help
}
case "${1:-}" in -h | --help | help)
	usage
	exit 0
	;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry disco "$@"
set -- "${GRID_ENTRY_ARGS[@]}"
GRID_DISCO_EVALUATE=TRUE
case "${1:-}" in
    prepare|evaluate)
        module=$1; shift
        source "$ROOT/f/0.common.sh"
        grid_activate_environment
        export PYTHONPATH="$ROOT/f${PYTHONPATH:+:$PYTHONPATH}"
        if [[ $module == prepare ]]; then
            exec "${GRID_PYTHON:-python3}" "$ROOT/f/shared_cohort.py" prepare "$@"
        fi
        exec "${GRID_PYTHON:-python3}" "$ROOT/f/disco_evaluate.py" "$@"
        ;;
    score) GRID_DISCO_EVALUATE=FALSE; shift ;;
    all) shift ;;
esac


# disco_run
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
	exec {score_read_lock}<>"$score_lock/write.lock"
	"${GRID_PYTHON:-python3}" "$ROOT/f/file_lock.py" -s "$score_read_lock"
	[[ -z $GRID_REMOVE || -f $GRID_REMOVE ]] || _grid_die "Missing withdrawal file: $GRID_REMOVE"
	echo "input score file: ${inputs[0]}"
	echo "output score files: $score_home/2disco.scores.rds; $score_home/2disco.coefficients.xlsx"
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
	sig=$(python3 "$io" signature "$GRID_DISCO_A" "$GRID_REGRESS_PCA" "$GRID_DISTANCE_PCS" --files "${inputs[@]}" "$GRID_PCA_FILE" "$GRID_MED_FILE" "$io" "$ROOT/3.disco.sh" "$ROOT/f/3.disco.R" "$ROOT/f/0.common.py")
	run="$work/$trait${suffix:+/${suffix#.}}"
	trait_log="$logdir/$trait${suffix:+/${suffix#.}}"
	mkdir -p "$run" "$trait_log"
	cache_sig="$run/disco.signature"
	exec {lock}<>"$work/$trait/run.lock"
	"${GRID_PYTHON:-python3}" "$ROOT/f/file_lock.py" -n "$lock" || _grid_die "Another Disco run is active for $trait"
	if [[ -s $score_home/2disco.scores.rds && -s $score_home/2disco.coefficients.xlsx && -s $cache_sig && $(cat "$cache_sig") == "$sig" && $GRID_REPLACE == FALSE ]]; then
		echo "SKIP $trait: matching permanent Disco results"
		disco_evaluate
		return 0
	fi
	# Invalidate completion before replacing files in this stable workspace.
	# A failed rebuild must not make a later invocation reuse partial publication.
	rm -f -- "$cache_sig" "$run/signature"
	grid_run_logged "$trait_log/prepare.log" python3 "$io" disco-inputs "$GRID_PCA_FILE" "$GRID_MED_FILE" "$score_home" "$run" "$GRID_DISTANCE_PCS" "$GRID_REMOVE"
	grep '^Disco ' "$trait_log/prepare.log"
	prs=()
	for p in "${POPS[@]}"; do prs+=("$run/$p.tsv"); done
	grid_run_logged "$trait_log/disco.log" "${disco_r[@]}" "$ROOT/f/3.disco.R" -m "$run/centers.tsv" -p "$run/pca.tsv" --prs.list "$(join_comma "${prs[@]}")" -s IID,PRS -A "$GRID_DISCO_A" --regress.PCA "$GRID_REGRESS_PCA" --print.coef TRUE -o "$run/disco"
	grid_run_logged "$trait_log/validate.log" python3 "$io" disco-output "$run/disco.tsv.gz" "$run"
	grid_run python3 "$ROOT/f/0.common.py" publish disco "$run/disco.tsv.gz" "$score_home/2disco.scores.rds" --remove "$GRID_REMOVE"
	grid_run "${disco_r[@]}" "$ROOT/../0f/results.R" import-table "$run/disco.coef.tsv.gz" "$score_home/2disco.coefficients.xlsx"
	# Permanent provenance stays usable even after deleting the scratch directory.
	{
		printf 'key\tvalue\n'
		printf 'PCA\t%s\ncenters\t%s\nPCs\t%s\nA\t%s\nregress_PCA\t%s\n' "$GRID_PCA_FILE" "$GRID_MED_FILE" "$GRID_DISTANCE_PCS" "$GRID_DISCO_A" "$GRID_REGRESS_PCA"
		for p in "${POPS[@]}"; do printf '%s\t%s\n' "$p" "$score_home/1csx.scores.rds"; done
	} >"$run/manifest.tsv"
	printf '%s\n' "$sig" >"$run/signature"
	publish "$run/signature" "$cache_sig"
	disco_evaluate
	echo "DONE $trait: $score_home/2disco.scores.rds"
}

disco_evaluate() {
    [[ $GRID_DISCO_EVALUATE == TRUE ]] || return 0
    python3 "$ROOT/f/disco_evaluate.py" --analysis-root "$GRID_OUTPUT_ROOT" --trait "$trait" \
        --scores "$score_home/2disco.scores.rds" "${GRID_BENCHMARK_ARGS[@]}"
}

source "$ROOT/f/0.common.sh"
grid_pipeline disco "$@"

#!/usr/bin/env bash
# Official DiscoDivas: population PRS + reference-projected PCs -> individual PRS.
if [[ -z ${GRID_COMMAND_FILE:-} ]]; then
  exec bash "$(dirname -- "${BASH_SOURCE[0]}")/pipeline.sh" disco "$@"
fi
set -euo pipefail
trait=$GRID_TRAIT
score_home="$GRID_SCORE_DIR/$trait${suffix:+/${suffix#.}}"
io="$ROOT/f/pipeline_io.py"
need "$GRID_PCA_FILE"; need "$GRID_MED_FILE"
echo "input PCA: $GRID_PCA_FILE"
echo "input reference centers: $GRID_MED_FILE"
inputs=("$score_home/csx.pgs.gz")
need "${inputs[0]}"
exec {score_read_lock}>>"$score_home/csx.pgs.gz.lock"
flock -s "$score_read_lock"
[[ -z $GRID_REMOVE || -f $GRID_REMOVE ]] || _grid_die "Missing withdrawal file: $GRID_REMOVE"
echo "input score file: ${inputs[0]}"
echo "output score files: $score_home/disco.pgs.gz; $score_home/disco.coef.tsv.gz"
[[ -z $GRID_REMOVE ]] || inputs+=("$GRID_REMOVE")
python3 - "$GRID_DISCO_A" "$GRID_DISTANCE_PCS" <<'PY'
import math,sys
x=list(map(float,sys.argv[1].split(',')))
assert len(x)==4 and all(math.isfinite(v) and v>=0 for v in x) and any(x), '--a-list requires four nonnegative values, at least one positive'
assert 5<=int(sys.argv[2])<=20, '--distance-pcs must be 5..20'
PY
source "$ROOT/f/r_runtime.sh"
grid_select_r data.table,dplyr,stringr,rio,optparse
disco_r=("${GRID_R[@]}")
echo "Disco R runtime: ${disco_r[*]}"
[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || { echo 'CHECK/PLAN complete; no Disco calculation executed'; return 0; }
sig=$(python3 "$io" signature "$GRID_DISCO_A" "$GRID_REGRESS_PCA" "$GRID_DISTANCE_PCS" --files "${inputs[@]}" "$GRID_PCA_FILE" "$GRID_MED_FILE" "$io" "$ROOT/f/disco.sh" "$ROOT/f/disco/DiscoDivas.R" "$ROOT/f/score_output.py")
run="$work/$trait/$sig"; mkdir -p "$run" "$work/log/$trait"
cache_sig="$work/$trait/disco.signature"
exec {lock}>"$work/$trait/run.lock"; flock -n "$lock" || _grid_die "Another Disco run is active for $trait"
if [[ -s $score_home/disco.pgs.gz && -s $score_home/disco.coef.tsv.gz && -s $cache_sig && $(cat "$cache_sig") == "$sig" && $GRID_REPLACE == FALSE ]]; then
  echo "SKIP $trait: matching permanent Disco results"; return 0
fi
grid_run_logged "$work/log/$trait/prepare.log" python3 "$io" disco-inputs "$GRID_PCA_FILE" "$GRID_MED_FILE" "$score_home" "$run" "$GRID_DISTANCE_PCS" "$GRID_REMOVE"
grep '^Disco ' "$work/log/$trait/prepare.log"
prs=(); for p in "${POPS[@]}"; do prs+=("$run/$p.tsv"); done
grid_run_logged "$work/log/$trait/disco.log" "${disco_r[@]}" "$ROOT/f/disco/DiscoDivas.R" -m "$run/centers.tsv" -p "$run/pca.tsv" --prs.list "$(join_comma "${prs[@]}")" -s IID,PRS -A "$GRID_DISCO_A" --regress.PCA "$GRID_REGRESS_PCA" --print.coef TRUE -o "$run/disco"
grid_run_logged "$work/log/$trait/validate.log" python3 "$io" disco-output "$run/disco.tsv.gz" "$run"
grid_run python3 "$ROOT/f/score_output.py" disco "$run/disco.tsv.gz" "$score_home/disco.pgs.gz" --remove "$GRID_REMOVE"
publish "$run/disco.coef.tsv.gz" "$score_home/disco.coef.tsv.gz"
printf '%s\n' "$sig" > "$run/signature"; publish "$run/signature" "$cache_sig"
# Permanent provenance stays usable even after deleting the scratch directory.
{ printf 'key\tvalue\n'; printf 'PCA\t%s\ncenters\t%s\nPCs\t%s\nA\t%s\nregress_PCA\t%s\n' "$GRID_PCA_FILE" "$GRID_MED_FILE" "$GRID_DISTANCE_PCS" "$GRID_DISCO_A" "$GRID_REGRESS_PCA"; for p in "${POPS[@]}"; do printf '%s\t%s\n' "$p" "$score_home/csx.pgs.gz"; done; } > "$run/manifest.tsv"
publish "$run/manifest.tsv" "$score_home/disco.manifest.tsv"
echo "DONE $trait: $score_home/disco.pgs.gz"

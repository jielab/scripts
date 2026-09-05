#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/grid_common.sh" "$@"
score_one(){
  local label=$1 c=$2 mode prefix w o
  mode=$(grid_target_mode "$c" || true); [[ -n $mode ]] || _grid_die "Missing target genotype chr$c"
  prefix="$GRID_TARGET_DIR/chr$c"; w="$gdir/weights/$label.chr$c.tsv"; o="$out/scores/tmp/$label.chr$c"; mkdir -p "$out/scores/tmp"
  [[ -s $w ]] || _grid_die "Missing $w"
  [[ -s $o.sscore && $GRID_REPLACE == FALSE ]] && return
  cmd=(plink2 "--$mode" "$prefix" --score "$w" 1 2 3 header-read no-mean-imputation cols=+scoresums --threads "$GRID_THREADS" --out "$o")
  [[ -z $GRID_KEEP ]] || cmd+=(--keep "$GRID_KEEP")
  [[ -z $GRID_REMOVE ]] || cmd+=(--remove "$GRID_REMOVE")
  grid_run_logged "$out/log/grid.score.$label.chr$c.log" "${cmd[@]}"
}
make_scores(){
  labels=(GRID_shared); for p in "${POPS[@]}"; do labels+=("GRID_${p^^}"); done
  files=()
  for label in "${labels[@]}"; do
    running=0;status=0
    for c in "${CHRS[@]}"; do score_one "$label" "$c" & ((++running)); if ((running>=GRID_JOBS)); then wait -n || status=1; running=$((running-1)); fi; done
    while ((running)); do wait -n || status=1; running=$((running-1)); done; ((status==0)) || _grid_die "Scoring failed for $label"
    inp=(); for c in "${CHRS[@]}"; do inp+=("$out/scores/tmp/$label.chr$c.sscore"); done
    grid_run python3 "$ROOT/f/combine_scores.py" --inputs "${inp[@]}" --name "$label" --output "$out/scores/$label.tsv.gz"; files+=("$out/scores/$label.tsv.gz")
  done
  grid_run python3 "$ROOT/f/merge_scores.py" --inputs "${files[@]}" --output "$out/scores/grid_population.tsv.gz"
  grid_run python3 "$ROOT/f/mix_population_scores.py" --scores "$out/scores/grid_population.tsv.gz" --ancestry "$GRID_ANCESTRY_FILE" --prefix GRID --out "$out/scores/grid_mixed.tsv.gz"
  grid_run python3 "$ROOT/f/merge_scores.py" --inputs "$out/scores/grid_mixed.tsv.gz" "$out/scores/GRID_shared.tsv.gz" --output "$out/scores/grid.tsv.gz"
  # Matched/posterior PRS-CSx baselines use the same genotype-only ancestry probabilities.
  if [[ -s $out/scores/csx.tsv.gz ]]; then grid_run python3 "$ROOT/f/mix_population_scores.py" --scores "$out/scores/csx.tsv.gz" --ancestry "$GRID_ANCESTRY_FILE" --prefix CSX --out "$out/scores/csx_mixed.tsv.gz"; fi
}

make_scores
echo "GRID score completed: $gdir"

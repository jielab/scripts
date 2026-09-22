#!/usr/bin/env bash
# Score every UKB individual with each population's permanent SNP weights.
for f in "${finals[@]}"; do
  need "$f"; need "$f.signature"
  [[ $(cat "$f.signature") == "$sig" ]] || _grid_die "Weights/settings mismatch: $f; rerun --stage weights with these settings"
done
score_inputs=("${finals[@]}" "$ROOT/f/csx_score.sh" "$ROOT/f/combine_scores.py" "$ROOT/f/score_output.py" "$ROOT/f/merge_scores.py")
for c in "${CHRS[@]}"; do
  mode=$(grid_target_mode "$c"); prefix="$GRID_TARGET_DIR/chr$c"
  if [[ $mode == pfile ]]; then
    score_inputs+=("$prefix.pgen" "$prefix.psam")
    if [[ -s $prefix.pvar ]]; then score_inputs+=("$prefix.pvar"); else score_inputs+=("$prefix.pvar.zst"); fi
  else score_inputs+=("$prefix.bed" "$prefix.bim" "$prefix.fam"); fi
done
[[ -z $GRID_KEEP ]] || score_inputs+=("$GRID_KEEP")
[[ -z $GRID_REMOVE ]] || score_inputs+=("$GRID_REMOVE")
score_sig=$(python3 "$io" signature "$sig" "$GRID_KEEP" "$GRID_REMOVE" --files "${score_inputs[@]}")
score_run="$run/scores/$score_sig"; mkdir -p "$score_run" "$score_home"
score_cache="$work/scores/$trait${suffix:+/${suffix#.}}"; mkdir -p "$score_cache"
for i in "${!POPS[@]}"; do
  p=${POPS[$i]}; dest="$score_cache/csx.$p.tsv.gz"
  if [[ -s $dest && -s $dest.signature && $(cat "$dest.signature") == "$score_sig" && $GRID_REPLACE == FALSE ]]; then echo "SKIP $trait $p: matching cached scores"; continue; fi
  score_chr(){
    local c=$1 mode prefix o
    mode=$(grid_target_mode "$c"); prefix="$GRID_TARGET_DIR/chr$c"; o="$score_run/$p.chr$c"
    if [[ -s $o.sscore && -f $o.done && $GRID_REPLACE == FALSE ]]; then return; fi
    rm -f -- "$o.done"
    local cmd=(plink2 "--$mode" "$prefix")
    if [[ $mode == pfile && ! -s $prefix.pvar ]]; then cmd+=(vzs); fi
    cmd+=(--score "${finals[$i]}" 1 2 3 header-read no-mean-imputation list-variants cols=+scoresums --threads "$GRID_THREADS" --memory 4096 --out "$o")
    [[ -z $GRID_KEEP ]] || cmd+=(--keep "$GRID_KEEP")
    [[ ! -s $GRID_REMOVE ]] || cmd+=(--remove "$GRID_REMOVE")
    grid_run_logged "$work/log/$trait/score.$p.chr$c.log" "${cmd[@]}" || return $?
    need "$o.sscore"; touch "$o.done"
  }
  echo "RUN $trait $p: UKB scoring"
  pids=()
  for c in "${CHRS[@]}"; do
    score_chr "$c" & pids+=("$!")
    if ((${#pids[@]} >= GRID_JOBS)); then status=0; for pid in "${pids[@]}"; do wait "$pid" || status=1; done; ((status==0)) || _grid_die 'CSx scoring failed'; pids=(); fi
  done
  status=0; for pid in "${pids[@]}"; do wait "$pid" || status=1; done; ((status==0)) || _grid_die 'CSx scoring failed'
  inputs=(); for c in "${CHRS[@]}"; do inputs+=("$score_run/$p.chr$c.sscore"); done
  grid_run_logged "$work/log/$trait/combine.$p.log" python3 "$ROOT/f/combine_scores.py" --inputs "${inputs[@]}" --name "CSX_$p" --output "$score_run/$p.tsv.gz"
  publish "$score_run/$p.tsv.gz" "$dest"
  printf '%s\n' "$score_sig" > "$score_run/signature"; publish "$score_run/signature" "$dest.signature"
done
merge=(); for p in "${POPS[@]}"; do merge+=("$score_cache/csx.$p.tsv.gz"); done
grid_run python3 "$ROOT/f/merge_scores.py" --inputs "${merge[@]}" --output "$score_run/csx.tsv.gz"
grid_run python3 "$ROOT/f/score_output.py" csx "$score_run/csx.tsv.gz" "$score_home/csx.pgs.gz" --remove "$GRID_REMOVE"
{
  printf 'trait\tpopulation\tinput_gwas\tweights\tukb_score\n'
  for i in "${!POPS[@]}"; do printf '%s\t%s\t%s\t%s\t%s\n' "$trait" "${POPS[$i]}" "${gwas[$i]}" "${finals[$i]}" "$score_home/csx.pgs.gz"; done
} > "$score_run/manifest.tsv"
publish "$score_run/manifest.tsv" "$score_home/csx.manifest.tsv"
echo "DONE $trait: $score_home/csx.pgs.gz"

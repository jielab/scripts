#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; POPS=($(grid_csv_words "$GRID_POPS")); CHRS=($(grid_expand_chrs "$GRID_CHRS"))
out="$GRID_OUTPUT_ROOT/$trait"; gdir="$out/grid"; mkdir -p "$gdir" "$out/log" "$out/scores" "$gdir/weights" "$gdir/model"
case "$GRID_GRID_ACTION" in ld|transport|weights|score|all) ;; *) _grid_die "bad --grid-action=$GRID_GRID_ACTION";; esac

ld_hdf5_for(){
  local p=${1,,} c=$2
  find "$GRID_CSX_REF_DIR" -type f -iname "*chr${c}*.hdf5" -ipath "*${p}*" -print | sort | head -n1
}
make_ld(){
  local p c f dest
  for p in "${POPS[@]}"; do p=${p^^}; mkdir -p "$GRID_LDSCORE_DIR/$p"
    for c in "${CHRS[@]}"; do
      dest="$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz"; [[ -s $dest && $GRID_REPLACE == FALSE ]] && continue
      f=$(ld_hdf5_for "$p" "$c"); [[ -s $f ]] || _grid_die "No PRS-CSx HDF5 LD file for $p chr$c below $GRID_CSX_REF_DIR"
      grid_run_logged "$out/log/grid.ld.$p.chr$c.log" python3 "$ROOT/f/extract_ld_scores.py" --hdf5 "$f" --pop "$p" --chr "$c" --out "$dest"
    done
  done
}
make_transport(){
  for c in "${CHRS[@]}"; do
    [[ -s $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz ]] || _grid_die "Missing ARG feature file $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
    if grid_is_true "$GRID_REQUIRE_LD"; then for p in "${POPS[@]}"; do [[ -s $GRID_LDSCORE_DIR/${p^^}/chr$c.ldscore.tsv.gz ]] || _grid_die "Missing LD score for ${p^^} chr$c"; done; fi
  done
  t="$gdir/transport.tsv.gz"
  cmd=(python3 "$ROOT/f/build_transport_table.py" --trait "$trait" --pops "$GRID_POPS" --chrs "${CHRS[*]}" --sumstats-dir "$out/sumstats/bychr" --arg-dir "$GRID_ARG_TREES_DIR" --ldscore-dir "$GRID_LDSCORE_DIR" --centers "$GRID_MED_FILE" --max-snps-per-chr "$GRID_MAX_SNPS_PER_CHR" --out "$t")
  [[ -z $GRID_EXTERNAL_AGE ]] || cmd+=(--external-age "$GRID_EXTERNAL_AGE")
  [[ -s $t && $GRID_REPLACE == FALSE ]] || grid_run_logged "$out/log/grid.transport.log" "${cmd[@]}"
  grid_run_logged "$out/log/grid.transport.fit.log" python3 "$ROOT/f/fit_transport_model.py" --input "$t" --out-dir "$gdir/model" --ridge-alpha "$GRID_RIDGE_ALPHA" --conservation-min "$GRID_CONSERVATION_MIN" --conservation-max "$GRID_CONSERVATION_MAX" --seed "$GRID_SEED"
}
make_weights(){
  [[ -s $gdir/model/variant_conservation.tsv.gz ]] || _grid_die 'Run --grid-action transport first'
  [[ -s $out/csx/manifest.tsv ]] || _grid_die 'Run grid.sh csx first'
  grid_run_logged "$out/log/grid.weights.log" python3 "$ROOT/f/make_grid_weights.py" --weights-dir "$out/csx/weights" --conservation "$gdir/model/variant_conservation.tsv.gz" --manifest "$out/csx/manifest.tsv" --out-dir "$gdir/weights" --pops "$GRID_POPS" --chrs "${CHRS[*]}"
}
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
case "$GRID_GRID_ACTION" in
 ld) make_ld;;
 transport) make_transport;;
 weights) make_weights;;
 score) make_scores;;
 all) make_ld; make_transport; make_weights; make_scores;;
esac
cat > "$gdir/GRID_RUN.txt" <<META
created=$(date -Is)
trait=$trait
chromosomes=${CHRS[*]}
method=genealogy-informed shrinkage of PRS-CSx population effects toward a shared effect
validation=blocked out-of-fold transportability model
local_genealogy=$GRID_ARG_TREES_DIR/chrCHR.variants.tsv.gz
ld_source=PRS-CSx HDF5 reference panels
participant_phenotype_used_for_weights=FALSE
META
echo "GRID module completed: $out/scores/grid.tsv.gz"

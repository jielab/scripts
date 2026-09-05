#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
SCRIPTS_ROOT=$(cd -- "$HERE/../.." && pwd -P)
GRID_ROOT=${REFGEN_GRID_ROOT:-$SCRIPTS_ROOT/grid}
GRID_F=$GRID_ROOT/f
[[ -s $GRID_F/common.sh ]] || { echo "ERROR: GRID helpers are missing: $GRID_F" >&2; exit 2; }

# Use the validated ARG-Needle environment without requiring an interactive
# `conda activate` in data-preparation jobs.
GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid-argneedle}
if [[ ! -x $GRID_CONDA_ENV/bin/python3 && -x $HOME/miniforge3/envs/grid/bin/python3 ]]; then
  GRID_CONDA_ENV=$HOME/miniforge3/envs/grid
fi
if [[ -x $GRID_CONDA_ENV/bin/python3 ]]; then export PATH="$GRID_CONDA_ENV/bin:$PATH"; fi
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null
[[ ! -d $GRID_CONDA_ENV/lib/R/library ]] || export R_LIBS_USER=$GRID_CONDA_ENV/lib/R/library

export GRID_ARG_ACTION=${REFGEN_ARG_ACTION:?}
export GRID_ARG_FORMAT=${REFGEN_ARG_FORMAT:-native}
export GRID_ARG_HAP_DIR=${REFGEN_PFILE_DIR:?}
export GRID_ARG_MAP_DIR=${REFGEN_MAP_DIR:-}
export GRID_ARG_MAP_PATTERN=${REFGEN_MAP_PATTERN:-}
export GRID_ARG_OUT=${REFGEN_ARG_DIR:?}
export GRID_ARG_TREES_DIR=$GRID_ARG_OUT/trees
export GRID_ARG_SCRATCH=${REFGEN_ARG_SCRATCH:-}
export GRID_ARG_FULL=${REFGEN_ARG_FULL:-FALSE}
export GRID_ARG_MAX_INDIVIDUALS=${REFGEN_MAX_INDIVIDUALS:-20000}
export GRID_ARG_SEED_HAPLOTYPES=${REFGEN_SEED_HAPLOTYPES:-4000}
export GRID_ARG_THREADS=${REFGEN_THREADS:-8}
export GRID_ARG_JOBS=${REFGEN_JOBS:-1}
export GRID_ARG_ANCHORS_PER_POP=${REFGEN_ANCHORS_PER_POP:-1000}
export GRID_ARG_KEEP=${REFGEN_ARG_KEEP:-}
export GRID_ARG_AFFINITY=${REFGEN_ARG_AFFINITY:-FALSE}
export GRID_ARG_NEEDLE_HOME=${REFGEN_ARG_NEEDLE_HOME:-}
export GRID_ANCESTRY_FILE=${REFGEN_ANCESTRY_FILE:?}
export GRID_CHRS=${REFGEN_ARG_CHRS:-1-22}
export GRID_REPLACE=${REFGEN_REPLACE:-FALSE}
export GRID_DRY_RUN=FALSE
export GRID_OUTPUT_ROOT=$GRID_ARG_OUT
export GRID_SEED=${REFGEN_SEED:-20260904}
export GRID_ANCHOR_PROB_MIN=${REFGEN_ANCHOR_PROB_MIN:-0.999}
export GRID_ARG_MAF_MIN=${REFGEN_ARG_MAF_MIN:-0.001}
export GRID_ARG_GENO_MAX=${REFGEN_ARG_GENO_MAX:-0.05}
export GRID_ARG_NORMALIZE=${REFGEN_ARG_NORMALIZE:-TRUE}
export GRID_ARG_WINDOW_BP=${REFGEN_ARG_WINDOW_BP:-1000000}

# shellcheck source=/dev/null
source "$GRID_F/common.sh"
grid_run(){ printf '[REFGEN ARG]'; printf ' %q' "$@"; printf '\n' >&2; "$@"; }
grid_run_logged(){
  local log=$1; shift
  mkdir -p "$(dirname -- "$log")"
  printf '[REFGEN ARG]'; printf ' %q' "$@"; printf '\n' | tee -a "$log" >&2
  "$@" >>"$log" 2>&1
}
if [[ -z $GRID_ARG_NEEDLE_HOME ]]; then
  GRID_ARG_NEEDLE_HOME=$(python3 - <<'PY' 2>/dev/null || true
import importlib.util
from pathlib import Path
s=importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
if s is not None and s.origin is not None: print(Path(s.origin).resolve().parents[2])
PY
)
fi
if [[ -z $GRID_ARG_NEEDLE_HOME && -s $HOME/software/arg-needle/arg-needle-scripts/arg_needle/infer_args_advanced.py ]]; then
  GRID_ARG_NEEDLE_HOME=$HOME/software/arg-needle/arg-needle-scripts
fi
read -r -a CHRS <<< "$(grid_expand_chrs "$GRID_CHRS")"
(( ${#CHRS[@]} )) || _grid_die 'No chromosomes selected'

WORK=${GRID_ARG_SCRATCH:+$GRID_ARG_SCRATCH/refgen_argneedle}
PREP=${REFGEN_GEN4ARG_DIR:?}
INF=${WORK:+$WORK/infer}
KEEP=$PREP/panel/arg.keep.txt
PANEL=$PREP/panel/arg.panel.tsv
ORDER=$PREP/panel/arg.order.txt
TRACE_DIR=$GRID_ARG_OUT/trace/needle

check_chr(){
  local c=$1 f
  for f in "$GRID_ARG_OUT/argn/chr$c.argn" "$GRID_ARG_TREES_DIR/chr$c.trees" \
           "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"; do
    [[ -s $f ]] || { echo "ERROR: missing ARG-Needle output for chr$c: $f" >&2; return 1; }
  done
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE ]]; then gzip -t "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"; fi
  python3 - "$GRID_ARG_TREES_DIR/chr$c.trees" "$c" <<'PY'
import sys,tskit
ts=tskit.load(sys.argv[1])
if min(ts.num_trees,ts.num_samples,ts.num_sites,ts.num_mutations) < 1:
    raise SystemExit(f"ERROR: empty ARG-Needle tree for chr{sys.argv[2]}: trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
print(f"ARG-Needle CHECK PASS chr{sys.argv[2]} trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
PY
  if [[ $c == X ]]; then
    python3 "$HERE/arg.needle_x.py" check --tree "$GRID_ARG_TREES_DIR/chrX.trees" --identities "$PREP/chrX.haplotypes.tsv" --output "$GRID_ARG_TREES_DIR/chrX.sample_map.tsv" --psam "$(grid_sample_for X)" --keep "$PREP/panel/X/arg.keep.txt"
  fi
}

trace_chr(){
  local c=$1
  local -a cmd=(python3 "$HERE/arg.trace.py" --method needle
    --tree "$GRID_ARG_TREES_DIR/chr$c.trees"
    --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
    --out-dir "$TRACE_DIR" --chr "$c")
  [[ $GRID_ARG_ACTION != check ]] || cmd+=(--check)
  [[ $GRID_REPLACE == FALSE ]] || cmd+=(--replace)
  grid_run "${cmd[@]}"
}

if [[ $GRID_ARG_ACTION == check ]]; then
  for c in "${CHRS[@]}"; do check_chr "$c"; done
  if [[ $GRID_ARG_FORMAT == trace ]]; then
    for c in "${CHRS[@]}"; do trace_chr "$c"; done
  fi
  exit 0
fi

[[ -n $GRID_ARG_SCRATCH ]] || _grid_die '--scratch is required for ARG-Needle construction'
case "$GRID_ARG_SCRATCH" in /mnt/[a-zA-Z]|/mnt/[a-zA-Z]/*) _grid_die '--scratch must use the native Linux filesystem, not /mnt/*';; esac
mkdir -p "$PREP/panel" "$PREP/log" "$INF"
exec 7>"$WORK/.lock"
flock -n 7 || _grid_die "another task uses scratch: $WORK"
owner=$(printf '%s\n' "$(realpath "$PREP")" "$(realpath -m "$GRID_ARG_OUT")")
if [[ -e $WORK/owner.txt ]]; then
  [[ $(cat "$WORK/owner.txt") == "$owner" ]] || _grid_die "scratch belongs to a different dataset/output; choose another --scratch: $WORK"
elif find "$INF" -type f -print -quit | grep -q . || [[ -d $WORK/prepare ]]; then
  _grid_die "unverified legacy scratch; choose a NEW --scratch directory: $WORK"
fi
printf '%s\n' "$owner" > "$WORK/owner.txt"
if [[ $GRID_ARG_ACTION != prepare ]]; then mkdir -p "$GRID_ARG_TREES_DIR" "$GRID_ARG_OUT/log"; fi

select_panel(){
  local dir=$PREP/panel
  [[ $1 != X ]] || dir=$dir/X
  KEEP=$dir/arg.keep.txt; PANEL=$dir/arg.panel.tsv; ORDER=$dir/arg.order.txt
}

make_panel_one(){
  local c=$1
  select_panel "$c"
  mkdir -p "$(dirname -- "$KEEP")"
  local sample; sample=$(grid_sample_for "$c" || true)
  [[ -s $sample ]] || _grid_die "Missing phased sample/PSAM for chr${CHRS[0]} under $GRID_ARG_HAP_DIR"
  local -a cmd=(python3 "$GRID_F/make_arg_keep.py" --sample "$sample" --ancestry "$GRID_ANCESTRY_FILE" --anchors-per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN" --seed "$GRID_SEED")
  if grid_is_true "$GRID_ARG_FULL"; then cmd+=(--full --max-individuals 0); else cmd+=(--max-individuals "$GRID_ARG_MAX_INDIVIDUALS"); fi
  [[ -z $GRID_ARG_KEEP ]] || cmd+=(--custom-keep "$GRID_ARG_KEEP")
  # Regenerate deterministically, but preserve timestamps when contents match.
  cmd+=(--out-keep "$KEEP.next" --out-panel "$PANEL.next")
  grid_run "${cmd[@]}"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.needle_x.py" panel --keep "$KEEP.next" --psam "$sample" --policy "${REFGEN_X_ODD_MALE:-drop-last}" --qc "$PREP/panel/X/selection.qc.json"
    python3 - "$KEEP.next" "$PANEL.next" <<'PY'
import sys,pandas as pd
from pathlib import Path
wanted={r.split()[-1] for r in Path(sys.argv[1]).read_text().splitlines()[1:] if r.strip()}
p=pd.read_csv(sys.argv[2],sep='\t',dtype={'eid':str}); p[p.eid.isin(wanted)].to_csv(sys.argv[2],sep='\t',index=False)
PY
  fi
  for f in "$KEEP" "$PANEL"; do
    if cmp -s "$f.next" "$f"; then rm -f "$f.next"; else mv -f "$f.next" "$f"; fi
  done
  tail -n +2 "$KEEP" > "$ORDER.next"
  if cmp -s "$ORDER.next" "$ORDER"; then rm -f "$ORDER.next"; else mv -f "$ORDER.next" "$ORDER"; fi
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE || $GRID_ARG_AFFINITY == TRUE ]]; then
    python3 - "$PANEL" "$GRID_ARG_ANCHORS_PER_POP" "$GRID_ANCHOR_PROB_MIN" <<'PY'
import sys,pandas as pd
p=pd.read_csv(sys.argv[1],sep='\t'); prob=pd.to_numeric(p['prob'],errors='coerce')
for pop in ('AFR','EAS','EUR','SAS'):
    n=min(int(sys.argv[2]),int(((p['pop']==pop)&(prob.isna()|(prob>=float(sys.argv[3])))).sum()))
    if n<10: raise SystemExit(f'GRID features need >=10 anchors for {pop} (have {n}); use --features FALSE for ARG/TRACE only')
PY
  fi
}

make_panel(){
  local c made=0
  for c in "${CHRS[@]}"; do
    if [[ $c == X ]]; then make_panel_one X
    elif ((made==0)); then make_panel_one "$c"; made=1; fi
  done
}

prepare_record(){
  local c=$1 pf s m
  select_panel "$c"
  pf=$(grid_pfile_for "$c" || true); s=$(grid_sample_for "$c" || true); m=$(grid_map_for "$c" || true)
  local -a args=(--value "build=${REFGEN_GRCH};maf=$GRID_ARG_MAF_MIN;geno=$GRID_ARG_GENO_MAX;mode=${REFGEN_NEEDLE_MODE};seed=$GRID_SEED"
    --text-file "$KEEP" --text-file "$PANEL" --file "$m" --file "$s"
    --text-file "$HERE/arg.needle.sh" --text-file "$HERE/arg.unique_haps.py" --text-file "$GRID_F/build_argneedle_map.py" --text-file "$GRID_F/validate_haps.py")
  if [[ -n $pf ]]; then
    args+=(--file "$pf.pgen")
    if [[ -s $pf.pvar.zst ]]; then args+=(--file "$pf.pvar.zst"); else args+=(--file "$pf.pvar"); fi
  else args+=(--file "$(grid_bgen_for "$c")"); fi
  [[ -z ${REFGEN_EXTRACT:-} ]] || args+=(--text-file "$REFGEN_EXTRACT")
  if [[ $c == X ]]; then
    args+=(--text-file "$HERE/arg.needle_x.py" --text-file "$HERE/arg_x.py" --text-file "$PREP/panel/X/selection.qc.json" --value "x_odd_male=${REFGEN_X_ODD_MALE:-drop-last}")
  fi
  args+=(--file "$(command -v plink2)")
  python3 "$HERE/arg.cache.py" "${args[@]}"
}
prepared_outputs(){
  local p=$PREP/chr$1
  local -a args=(--file "$p.haps.gz" --text-file "$p.sample" --text-file "$p.map" --text-file "$p.qc.json" --text-file "$p.positions.qc.json")
  if [[ $1 == X ]]; then args+=(--text-file "$p.haplotypes.tsv" --text-file "$p.packed.keep"); fi
  python3 "$HERE/arg.cache.py" "${args[@]}"
}
require_prepared(){
  local c=$1 p=$PREP/chr$1
  prepare_record "$c" > "$p.request.next"
  [[ -s $p.request.json ]] && cmp -s "$p.request.next" "$p.request.json" || _grid_die "prepared chr$c does not match this request; run prep_gen with the SAME input/panel/filter options"
  prepared_outputs "$c" > "$p.outputs.next"
  cmp -s "$p.outputs.next" "$p.outputs.json" || _grid_die "prepared chr$c was modified or incomplete; run prep_gen"
  rm -f "$p.request.next" "$p.outputs.next"
}

export_chr(){
  local c=$1 b pf s srcmap prefix subset haps sampleout gz
  select_panel "$c"
  b=$(grid_bgen_for "$c" || true); pf=$(grid_pfile_for "$c" || true)
  s=$(grid_sample_for "$c" || true); srcmap=$(grid_map_for "$c" || true)
  [[ -s $b || -n $pf ]] || _grid_die "Missing phased BGEN or PGEN for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $s ]] || _grid_die "Missing Oxford sample or PSAM for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $srcmap ]] || _grid_die "Missing GRCh${REFGEN_GRCH:-37} genetic map for chr$c"
  prefix=$PREP/chr$c; subset=$prefix.subset
  prepare_record "$c" > "$prefix.request.next"
  if [[ $GRID_REPLACE == FALSE && -s $prefix.outputs.json ]] && cmp -s "$prefix.request.next" "$prefix.request.json"; then
    if prepared_outputs "$c" > "$prefix.outputs.next" && cmp -s "$prefix.outputs.next" "$prefix.outputs.json"; then
      echo "SKIP verified genetic input: $prefix"; rm -f "$prefix.request.next" "$prefix.outputs.next"; return
    fi
  fi
  rm -f "$prefix.request.json" "$prefix.outputs.json"
  rm -f "$prefix.haps" "$prefix.haps.gz" "$prefix.sample" "$prefix.map" "$prefix.log"
  local -a cmd
  if [[ -n $pf ]]; then
    cmd=(plink2 --pfile "$pf"); [[ ! -s $pf.pvar.zst ]] || cmd+=(vzs)
  else cmd=(plink2 --bgen "$b" ref-first --sample "$s"); fi
  cmd+=(--threads "$GRID_ARG_THREADS")
  if [[ $c == X ]]; then
    [[ -n $pf ]] || _grid_die 'X preparation requires a PGEN/PSAM input with known sex'
    local lo=2699521 hi=154931043
    if [[ $REFGEN_GRCH == 38 ]]; then lo=2781480; hi=155701382; fi
    cmd+=(--chr X --from-bp "$lo" --to-bp "$hi")
  fi
  [[ -z ${REFGEN_EXTRACT:-} ]] || cmd+=(--extract "$REFGEN_EXTRACT")
  cmd+=(--keep "$KEEP" --indiv-sort f "$ORDER" --snps-only just-acgt --max-alleles 2 --maf "$GRID_ARG_MAF_MIN" --geno "$GRID_ARG_GENO_MAX" --sort-vars --make-pgen vzs --out "$subset")
  grid_run_logged "$PREP/log/arg.prepare.chr$c.log" "${cmd[@]}"
  grid_run_logged "$PREP/log/arg.prepare.chr$c.log" plink2 --pfile "$subset" vzs --threads "$GRID_ARG_THREADS" --export haps --out "$prefix"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.needle_x.py" pack --prefix "$prefix" --psam "$subset.psam" --keep "$KEEP" --build "$REFGEN_GRCH"
  fi
  grid_run python3 "$HERE/arg.unique_haps.py" --haps "$prefix.haps" --chr "$c" --out "$prefix.positions.qc.json"
  haps=$(find "$PREP" -maxdepth 1 -type f \( -name "chr$c.haps" -o -name "chr$c.haps.gz" \) -print -quit)
  sampleout=$(find "$PREP" -maxdepth 1 -type f -name "chr$c.sample" -print -quit)
  [[ -s $haps && -s $sampleout ]] || _grid_die "PLINK2 did not create chr$c Oxford HAPS/sample"
  if [[ $haps != *.gz ]]; then
    if command -v pigz >/dev/null 2>&1; then grid_run pigz -f -p "$GRID_ARG_THREADS" "$haps"; else grid_run gzip -f "$haps"; fi
    gz=$haps.gz
  else gz=$haps; fi
  [[ $gz == "$prefix.haps.gz" ]] || mv -f "$gz" "$prefix.haps.gz"
  rm -f "$subset.pgen" "$subset.pvar" "$subset.pvar.zst" "$subset.psam" "$subset.log"
  grid_run python3 "$GRID_F/build_argneedle_map.py" --source "$srcmap" --haps "$prefix.haps.gz" --chr "$c" --out "$prefix.map"
  local validation_keep=$KEEP
  [[ $c != X ]] || validation_keep=$prefix.packed.keep
  grid_run python3 "$GRID_F/validate_haps.py" --haps "$prefix.haps.gz" --sample "$prefix.sample" --keep "$validation_keep" --map "$prefix.map" --out "$prefix.qc.json"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.needle_x.py" qc --qc "$prefix.qc.json" --identities "$prefix.haplotypes.tsv"
  fi
  prepared_outputs "$c" > "$prefix.outputs.next"
  mv -f "$prefix.outputs.next" "$prefix.outputs.json"
  mv -f "$prefix.request.next" "$prefix.request.json"
}

infer_chr(){
  local c=$1 prefix nind seed final
  select_panel "$c"
  prefix=$PREP/chr$c
  require_prepared "$c"
  mkdir -p "$WORK/input"
  # Stage verified inputs on native Linux storage without changing identical files.
  for ext in haps.gz sample map; do
    if [[ ! -s $WORK/input/chr$c.$ext ]] || ! cmp -s "$prefix.$ext" "$WORK/input/chr$c.$ext"; then
      cp "$prefix.$ext" "$WORK/input/chr$c.$ext.next"; mv -f "$WORK/input/chr$c.$ext.next" "$WORK/input/chr$c.$ext"
    fi
  done
  prefix=$WORK/input/chr$c
  [[ -s $prefix.haps.gz && -s $prefix.sample && -s $prefix.map ]] || _grid_die "Run arg.sh prep_gen first for chr$c"
  nind=$(( $(wc -l < "$prefix.sample") - 2 )); seed=$GRID_ARG_SEED_HAPLOTYPES
  ((seed<=2*nind)) || seed=$((2*nind))
  ((2*nind>=300)) || _grid_die "ARG-Needle ASMC decoding requires at least 300 total haplotypes (selected $((2*nind)))"
  if [[ ${REFGEN_NEEDLE_MODE} == sequence && $seed -ge $((2*nind)) ]]; then seed=$((2*nind-1)); fi
  ((seed>=2)) || _grid_die "Too few scaffold haplotypes for chr$c"
  local chromosome=$c
  [[ $c != X ]] || chromosome=23
  if [[ $c == X ]]; then echo "Needle X: $(wc -l < "$PREP/chrX.haplotypes.tsv") identity rows (including header); normalization=$(chr_normalize "$c")"; fi
  local -a cmd=(python3 "$GRID_F/run_argneedle_advanced.py" --home "$GRID_ARG_NEEDLE_HOME" --haps "$prefix.haps.gz" --map "$prefix.map" --out "$INF/chr$c" --chr "$chromosome" --seed-haplotypes "$seed" --threads "$GRID_ARG_THREADS" --normalize "$(chr_normalize "$c")" --random-seed "$GRID_SEED" --log "$GRID_ARG_OUT/log/arg.infer.chr$c.log")
  [[ $GRID_REPLACE == FALSE ]] || cmd+=(--replace)
  cmd+=(--mode "$REFGEN_NEEDLE_MODE")
  [[ -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ]] || cmd+=(--normalize-demography "$REFGEN_NORMALIZE_DEMOGRAPHY")
  grid_run "${cmd[@]}"
  final=$INF/chr$c.argn; [[ -s $final ]] || _grid_die "Missing final ARG-Needle file $final"
  mkdir -p "$GRID_ARG_OUT/argn"
  if ! cmp -s "$final" "$GRID_ARG_OUT/argn/chr$c.argn"; then
    cp "$final" "$GRID_ARG_OUT/argn/chr$c.argn.next"; mv -f "$GRID_ARG_OUT/argn/chr$c.argn.next" "$GRID_ARG_OUT/argn/chr$c.argn"
  fi
  inferred_record "$c" > "$GRID_ARG_OUT/argn/chr$c.input.json.next"
  mv -f "$GRID_ARG_OUT/argn/chr$c.input.json.next" "$GRID_ARG_OUT/argn/chr$c.input.json"
}

chr_normalize(){
  # The bundled CEU normalization is not an X-specific demography.
  if [[ $GRID_ARG_NORMALIZE != TRUE || ( $1 == X && -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ) ]]; then echo 0; else echo 1; fi
}

inferred_record(){
  local c=$1
  local -a args=(--file "$GRID_ARG_OUT/argn/chr$c.argn" --text-file "$PREP/chr$c.request.json" --text-file "$PREP/chr$c.outputs.json"
    --value "mode=$REFGEN_NEEDLE_MODE;normalize=$(chr_normalize "$c");scaffold=$GRID_ARG_SEED_HAPLOTYPES;seed=$GRID_SEED")
  [[ -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ]] || args+=(--file "$REFGEN_NORMALIZE_DEMOGRAPHY")
  python3 "$HERE/arg.cache.py" "${args[@]}"
}
require_inferred(){
  local c=$1
  require_prepared "$c"
  inferred_record "$c" > "$WORK/chr$c.inferred.next"
  cmp -s "$WORK/chr$c.inferred.next" "$GRID_ARG_OUT/argn/chr$c.input.json" || _grid_die "ARG chr$c does not match prepared input/options; run infer first"
  rm -f "$WORK/chr$c.inferred.next"
}

convert_chr(){
  local c=$1 argn=$GRID_ARG_OUT/argn/chr$1.argn trees=$GRID_ARG_TREES_DIR/chr$1.trees sampleout=$PREP/chr$1.sample
  select_panel "$c"
  [[ -s $argn && -s $sampleout ]] || _grid_die "Missing chr$c ARG/HAPS preparation output"
  require_inferred "$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.convert.chr$c.log" python3 "$GRID_F/argn_to_trees.py" --argn "$argn" --haps "$PREP/chr$c.haps.gz" --out "$trees" --home "$GRID_ARG_NEEDLE_HOME"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.needle_x.py" map --tree "$trees" --identities "$PREP/chrX.haplotypes.tsv" --output "$GRID_ARG_TREES_DIR/chrX.sample_map.tsv"
  else
    grid_run python3 "$GRID_F/make_sample_map.py" --trees "$trees" --sample "$sampleout" --out "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
  fi
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE || $GRID_ARG_AFFINITY == TRUE ]]; then
    grid_run python3 "$GRID_F/make_anchors.py" --panel "$PANEL" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --out "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN"
  fi
}

features_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.variants.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $PREP/chr$c.haps.gz && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing tree/HAPS/anchors for chr$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.features.chr$c.log" python3 "$GRID_F/arg_features.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --haps "$PREP/chr$c.haps.gz" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --chr "$c" --window-bp "$GRID_ARG_WINDOW_BP" --out "$out"
}

affinity_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.affinity.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $GRID_ARG_TREES_DIR/chr$c.sample_map.tsv && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing ARG conversion files for chr$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.affinity.chr$c.log" python3 "$GRID_F/arg_affinity.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --out "$out"
}

source "$HERE/arg.jobs.sh"

case "$GRID_ARG_ACTION" in
  prepare) make_panel; run_parallel export_chr "${CHRS[@]}";;
  infer) make_panel; run_parallel infer_chr "${CHRS[@]}";;
  convert) make_panel; run_parallel convert_chr "${CHRS[@]}";;
  features) make_panel; run_parallel features_chr "${CHRS[@]}";;
  affinity) run_parallel affinity_chr "${CHRS[@]}";;
  all)
    make_panel
    run_parallel export_chr "${CHRS[@]}"
    run_parallel infer_chr "${CHRS[@]}"
    run_parallel convert_chr "${CHRS[@]}"
    if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE ]]; then run_parallel features_chr "${CHRS[@]}"; fi
    if grid_is_true "$GRID_ARG_AFFINITY"; then run_parallel affinity_chr "${CHRS[@]}"; fi
    ;;
esac

if [[ $GRID_ARG_FORMAT == trace && ( $GRID_ARG_ACTION == all || $GRID_ARG_ACTION == convert ) ]]; then
  run_parallel trace_chr "${CHRS[@]}"
  printf 'method\tneedle\nformat\ttrace\nbuild\tGRCh%s\nchromosomes\t%s\n' \
    "${REFGEN_GRCH:-37}" "${CHRS[*]}" > "$TRACE_DIR/ARG_TRACE_BUILD.tsv"
fi

if [[ $GRID_ARG_ACTION == prepare ]]; then echo "Genetic input ready: $PREP"; exit 0; fi
cat > "$GRID_ARG_OUT/ARG_NEEDLE_BUILD.txt" <<META
created=$(date -Is)
producer=/mnt/d/scripts/gu/arg.sh build
method=needle
build=GRCh${REFGEN_GRCH:-37}
genotype_dir=$GRID_ARG_HAP_DIR
chromosomes=${CHRS[*]}
full=$GRID_ARG_FULL
max_individuals=$GRID_ARG_MAX_INDIVIDUALS
seed_haplotypes=$GRID_ARG_SEED_HAPLOTYPES
scratch=$GRID_ARG_SCRATCH
META

echo "ARG-Needle data preparation complete: $GRID_ARG_OUT"

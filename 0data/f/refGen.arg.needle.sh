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
read -r -a CHRS <<< "$(grid_expand_chrs "$GRID_CHRS")"
(( ${#CHRS[@]} )) || _grid_die 'No chromosomes selected'

WORK=${GRID_ARG_SCRATCH:+$GRID_ARG_SCRATCH/refgen_argneedle}
PREP=${WORK:+$WORK/prepare}
INF=${WORK:+$WORK/infer}
KEEP=$GRID_ARG_OUT/panel/arg.keep.txt
PANEL=$GRID_ARG_OUT/panel/arg.panel.tsv
ORDER=$GRID_ARG_OUT/panel/arg.order.txt
TRACE_DIR=$GRID_ARG_OUT/trace/needle

check_chr(){
  local c=$1 f
  for f in "$GRID_ARG_OUT/argn/chr$c.argn" "$GRID_ARG_TREES_DIR/chr$c.trees" \
           "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" \
           "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"; do
    [[ -s $f ]] || { echo "ERROR: missing ARG-Needle output for chr$c: $f" >&2; return 1; }
  done
  gzip -t "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
  python3 - "$GRID_ARG_TREES_DIR/chr$c.trees" "$c" <<'PY'
import sys,tskit
ts=tskit.load(sys.argv[1])
if min(ts.num_trees,ts.num_samples,ts.num_sites,ts.num_mutations) < 1:
    raise SystemExit(f"ERROR: empty ARG-Needle tree for chr{sys.argv[2]}: trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
print(f"ARG-Needle CHECK PASS chr{sys.argv[2]} trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
PY
}

trace_chr(){
  local c=$1
  local -a cmd=(python3 "$HERE/refGen.arg.trace.py" --method needle
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
mkdir -p "$PREP" "$INF" "$GRID_ARG_OUT" "$GRID_ARG_TREES_DIR" "$GRID_ARG_OUT/panel" "$GRID_ARG_OUT/log"

REQUEST_FILE=$GRID_ARG_OUT/ARG_NEEDLE_REQUEST.tsv
REQUEST_EXPECTED=$WORK/ARG_NEEDLE_REQUEST.expected.tsv
write_request(){
  local c p s m
  printf 'schema\t1\nmethod\tneedle\nbuild\tGRCh%s\nfull\t%s\nmax_individuals\t%s\nseed_haplotypes\t%s\nanchors_per_pop\t%s\nancestry_prob_min\t%s\nmaf_min\t%s\ngeno_max\t%s\nnormalize\t%s\n' \
    "${REFGEN_GRCH:-37}" "$GRID_ARG_FULL" "$GRID_ARG_MAX_INDIVIDUALS" "$GRID_ARG_SEED_HAPLOTYPES" \
    "$GRID_ARG_ANCHORS_PER_POP" "$GRID_ANCHOR_PROB_MIN" "$GRID_ARG_MAF_MIN" "$GRID_ARG_GENO_MAX" "$GRID_ARG_NORMALIZE"
  stat -c 'ancestry\t%n:%s:%Y' "$GRID_ANCESTRY_FILE"
  [[ -z $GRID_ARG_KEEP ]] || stat -c 'custom_keep\t%n:%s:%Y' "$GRID_ARG_KEEP"
  stat -c 'helper\t%n:%s:%Y' "$GRID_F/make_arg_keep.py" "$GRID_F/run_argneedle_advanced.py" "$GRID_F/argn_to_trees.py" "$GRID_F/arg_features.py"
  for c in "${CHRS[@]}"; do
    p=$(grid_pfile_for "$c" || true); s=$(grid_sample_for "$c" || true); m=$(grid_map_for "$c" || true)
    if [[ -n $p ]]; then stat -c "genotype_chr$c\t%n:%s:%Y" "$p.pgen" "$p.psam"; else p=$(grid_bgen_for "$c" || true); stat -c "genotype_chr$c\t%n:%s:%Y" "$p" "$s"; fi
    stat -c "map_chr$c\t%n:%s:%Y" "$m"
  done
}
if [[ $GRID_ARG_ACTION == prepare || $GRID_ARG_ACTION == all ]]; then
  write_request > "$REQUEST_EXPECTED"
  if [[ -s $PANEL ]] && { [[ ! -s $REQUEST_FILE ]] || ! cmp -s "$REQUEST_EXPECTED" "$REQUEST_FILE"; }; then
    echo "WARNING: ARG-Needle request differs from the existing panel; rebuilding requested outputs." >&2
    GRID_REPLACE=TRUE
  fi
fi

make_panel(){
  if [[ -s $KEEP && -s $PANEL && $GRID_REPLACE == FALSE ]]; then tail -n +2 "$KEEP" > "$ORDER"; return; fi
  local sample; sample=$(grid_sample_for "${CHRS[0]}" || true)
  [[ -s $sample ]] || _grid_die "Missing phased sample/PSAM for chr${CHRS[0]} under $GRID_ARG_HAP_DIR"
  local -a cmd=(python3 "$GRID_F/make_arg_keep.py" --sample "$sample" --ancestry "$GRID_ANCESTRY_FILE" --out-keep "$KEEP" --out-panel "$PANEL" --anchors-per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN" --seed "$GRID_SEED")
  if grid_is_true "$GRID_ARG_FULL"; then cmd+=(--full --max-individuals 0); else cmd+=(--max-individuals "$GRID_ARG_MAX_INDIVIDUALS"); fi
  [[ -z $GRID_ARG_KEEP ]] || cmd+=(--custom-keep "$GRID_ARG_KEEP")
  grid_run "${cmd[@]}"
  tail -n +2 "$KEEP" > "$ORDER"
}

export_chr(){
  local c=$1 b pf s srcmap prefix subset haps sampleout gz
  b=$(grid_bgen_for "$c" || true); pf=$(grid_pfile_for "$c" || true)
  s=$(grid_sample_for "$c" || true); srcmap=$(grid_map_for "$c" || true)
  [[ -s $b || -n $pf ]] || _grid_die "Missing phased BGEN or PGEN for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $s ]] || _grid_die "Missing Oxford sample or PSAM for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $srcmap ]] || _grid_die "Missing GRCh${REFGEN_GRCH:-37} genetic map for chr$c"
  prefix=$PREP/chr$c; subset=$prefix.subset
  if [[ -s $prefix.haps.gz && -s $prefix.sample && -s $prefix.map && $GRID_REPLACE == FALSE ]]; then return; fi
  rm -f "$prefix.haps" "$prefix.haps.gz" "$prefix.sample" "$prefix.map" "$prefix.log"
  if [[ -s $b ]]; then cmd=(plink2 --bgen "$b" ref-first --sample "$s"); else cmd=(plink2 --pfile "$pf" vzs); fi
  cmd+=(--keep "$KEEP" --indiv-sort f "$ORDER" --snps-only just-acgt --max-alleles 2 --maf "$GRID_ARG_MAF_MIN" --geno "$GRID_ARG_GENO_MAX" --make-pgen vzs --out "$subset")
  grid_run_logged "$GRID_ARG_OUT/log/arg.prepare.chr$c.log" "${cmd[@]}"
  grid_run_logged "$GRID_ARG_OUT/log/arg.prepare.chr$c.log" plink2 --pfile "$subset" vzs --export haps --out "$prefix"
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
  grid_run python3 "$GRID_F/validate_haps.py" --haps "$prefix.haps.gz" --sample "$prefix.sample" --out "$GRID_ARG_OUT/log/arg.prepare.chr$c.qc.json"
}

infer_chr(){
  local c=$1 prefix nind seed final
  prefix=$PREP/chr$c
  [[ -s $prefix.haps.gz && -s $prefix.sample && -s $prefix.map ]] || _grid_die "Run make-arg --action prepare first for chr$c"
  nind=$(( $(wc -l < "$prefix.sample") - 2 )); seed=$GRID_ARG_SEED_HAPLOTYPES
  ((seed<=2*nind)) || seed=$((2*nind))
  ((2*nind>=300)) || _grid_die "ARG-Needle normalization needs at least 300 total haplotypes (selected $((2*nind)))"
  ((seed>=2)) || _grid_die "Too few scaffold haplotypes for chr$c"
  local -a cmd=(python3 "$GRID_F/run_argneedle_advanced.py" --home "$GRID_ARG_NEEDLE_HOME" --haps "$prefix.haps.gz" --map "$prefix.map" --out "$INF/chr$c" --chr "$c" --seed-haplotypes "$seed" --threads "$GRID_ARG_THREADS" --normalize "$([[ $GRID_ARG_NORMALIZE == TRUE ]] && echo 1 || echo 0)" --random-seed "$GRID_SEED" --log "$GRID_ARG_OUT/log/arg.infer.chr$c.log")
  [[ $GRID_REPLACE == FALSE ]] || cmd+=(--replace)
  grid_run "${cmd[@]}"
  final=$INF/chr$c.argn; [[ -s $final ]] || _grid_die "Missing final ARG-Needle file $final"
  mkdir -p "$GRID_ARG_OUT/argn"; cp -f "$final" "$GRID_ARG_OUT/argn/chr$c.argn"
}

convert_chr(){
  local c=$1 argn=$GRID_ARG_OUT/argn/chr$1.argn trees=$GRID_ARG_TREES_DIR/chr$1.trees sampleout=$PREP/chr$1.sample
  [[ -s $argn && -s $sampleout ]] || _grid_die "Missing chr$c ARG/HAPS preparation output"
  [[ ! -s $trees || $GRID_REPLACE == TRUE ]] || return
  grid_run_logged "$GRID_ARG_OUT/log/arg.convert.chr$c.log" python3 "$GRID_F/argn_to_trees.py" --argn "$argn" --haps "$PREP/chr$c.haps.gz" --out "$trees" --home "$GRID_ARG_NEEDLE_HOME"
  grid_run python3 "$GRID_F/make_sample_map.py" --trees "$trees" --sample "$sampleout" --out "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
  grid_run python3 "$GRID_F/make_anchors.py" --panel "$PANEL" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --out "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN"
}

features_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.variants.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $PREP/chr$c.haps.gz && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing tree/HAPS/anchors for chr$c"
  [[ ! -s $out || $GRID_REPLACE == TRUE ]] || return
  grid_run_logged "$GRID_ARG_OUT/log/arg.features.chr$c.log" python3 "$GRID_F/arg_features.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --haps "$PREP/chr$c.haps.gz" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --chr "$c" --window-bp "$GRID_ARG_WINDOW_BP" --out "$out"
}

affinity_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.affinity.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $GRID_ARG_TREES_DIR/chr$c.sample_map.tsv && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing ARG conversion files for chr$c"
  [[ ! -s $out || $GRID_REPLACE == TRUE ]] || return
  grid_run_logged "$GRID_ARG_OUT/log/arg.affinity.chr$c.log" python3 "$GRID_F/arg_affinity.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --out "$out"
}

run_parallel(){
  local fn=$1 c running=0 status=0; shift
  for c in "$@"; do "$fn" "$c" & ((++running)); if ((running>=GRID_ARG_JOBS)); then wait -n || status=1; running=$((running-1)); fi; done
  while ((running)); do wait -n || status=1; running=$((running-1)); done
  ((status==0)) || _grid_die "$fn failed"
}

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
    run_parallel features_chr "${CHRS[@]}"
    grid_is_true "$GRID_ARG_AFFINITY" && run_parallel affinity_chr "${CHRS[@]}"
    ;;
esac

if [[ $GRID_ARG_FORMAT == trace ]]; then
  run_parallel trace_chr "${CHRS[@]}"
  printf 'method\tneedle\nformat\ttrace\nbuild\tGRCh%s\nchromosomes\t%s\n' \
    "${REFGEN_GRCH:-37}" "${CHRS[*]}" > "$TRACE_DIR/ARG_TRACE_BUILD.tsv"
fi

cat > "$GRID_ARG_OUT/ARG_NEEDLE_BUILD.txt" <<META
created=$(date -Is)
producer=/mnt/d/scripts/0data/refGen.sh make-arg
method=needle
build=GRCh${REFGEN_GRCH:-37}
genotype_dir=$GRID_ARG_HAP_DIR
chromosomes=${CHRS[*]}
full=$GRID_ARG_FULL
max_individuals=$GRID_ARG_MAX_INDIVIDUALS
seed_haplotypes=$GRID_ARG_SEED_HAPLOTYPES
scratch=$GRID_ARG_SCRATCH
META
if [[ $GRID_ARG_ACTION == prepare || $GRID_ARG_ACTION == all ]]; then mv -f "$REQUEST_EXPECTED" "$REQUEST_FILE"; fi
echo "ARG-Needle data preparation complete: $GRID_ARG_OUT"

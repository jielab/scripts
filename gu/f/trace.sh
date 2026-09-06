#!/usr/bin/env bash
# TRACE computation driver. Shiny owns reporting and aggregation.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd); F=$ROOT/f
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }; }
ACTION=${1:-run}; case "$ACTION" in run|extract|infer|segments|check) ;; *) echo "ERROR: invalid TRACE action: $ACTION" >&2; exit 2;; esac
ARG_DIR=${GU_ARG_DIR:-${GU_TARGET_ROOT:-.}/arg}; TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-}
OUT=${TRACE_OUT:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu}/trace/${GU_SCOPE_ID:-genome}}
BUILD=GRCh${GU_BUILD:-37}; CHRS=${TRACE_CHRS:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}}; LOCI_FILE=${TRACE_LOCI_FILE:-${GU_LOCI_FILE:-}}
LOCI_MODE=${TRACE_LOCI_MODE:-posthoc}
case "$LOCI_MODE" in posthoc|extract) ;; *) echo "ERROR: internal trace-loci-mode must be posthoc or extract" >&2; exit 2;; esac
GROUP_SIZE=${TRACE_GROUP_SIZE:-50}; JEX=${TRACE_JOB_EXTRACT:-1}; JIN=${TRACE_JOB_INFER:-8}; JSUM=${TRACE_JOB_SUMMARIZE:-8}
TARC=${TRACE_T_ARCHAIC:-15000}; POST=${TRACE_POSTERIOR_THRESHOLD:-0.90}; MINBP=${TRACE_PHYSICAL_LENGTH_THRESHOLD:-50000}; MINCM=${TRACE_GENETIC_DISTANCE_THRESHOLD:-0.05}; WIN=${TRACE_WINDOW_SIZE:-1000}
SAMPLE_MAP=${TRACE_SAMPLE_MAP:-$OUT/samples/trace_sample_map.tsv}; TARGET_SAMPLES=${TRACE_TARGET_SAMPLES:-}
TRACE_ARG_METHOD=$(awk -F'\t' '$1=="method"{print $2;exit}' "$ARG_DIR/ARG_TRACE_BUILD.tsv" 2>/dev/null || true)
REPLACE=${TRACE_REPLACE:-0}; [[ $REPLACE == 0 || $REPLACE == 1 ]] || { echo "ERROR: internal replace-trace state must be 0 or 1" >&2; exit 2; }
TRACE_FINAL=$OUT/final/trace_haplotype_segments.tsv.gz
if [[ $ACTION == run && $REPLACE == 0 && -s $TRACE_FINAL ]]; then
  echo "[GU TRACE] SKIP reason=output_exists final=$TRACE_FINAL"
  exit 0
fi
for x in python3 bcftools awk sort find gzip stat cmp grep diff comm; do need "$x"; done
if [[ $ACTION != check ]]; then for x in trace-extract trace-infer trace-summarize; do need "$x"; done; fi
mkdir -p "$OUT"/{manifest,samples/groups,regions,extract,datafiles,infer,calls,final,log}
if [[ $ACTION != check && $REPLACE == 1 ]]; then
  rm -rf -- "$OUT/extract" "$OUT/datafiles" "$OUT/infer" "$OUT/calls" "$OUT/final"
  mkdir -p "$OUT"/{extract,datafiles,infer,calls,final}
fi

# GU TRACE consumes mutation-bearing, TRACE-formatted trees stored directly in
# ARG_DIR. Native method artifacts stay in their method directories and are
# exported explicitly with `arg.sh build --format trace`.
find_trees(){ local c=$1; find -L "$ARG_DIR" -maxdepth 1 -type f \( -name '*.trees' -o -name '*.tsz' \) -print | awk -v c="$c" 'BEGIN{IGNORECASE=1}{b=$0;gsub(/.*\//,"",b);if(b~("(^|[^A-Za-z0-9])chr"c"([^A-Za-z0-9]|$)"))print;else if(b~("(^|[^A-Za-z0-9])"c"([^A-Za-z0-9]|$)"))print}' | sort; }
validate_tree_scope(){
  local tree=$1 chr=$2
  python3 - "$tree" "$chr" "${TRACE_MIN_FULL_ARG_SITES:-10000}" "${TRACE_MIN_FULL_ARG_SPAN_FRACTION:-0.5}" <<'PY'
import sys
from pathlib import Path
tree, chrom, min_sites, min_span = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
if tree.suffix == ".tsz":
    import tszip
    ts = tszip.decompress(tree)
else:
    import tskit
    ts = tskit.load(tree)
positions = [site.position for site in ts.sites()]
if not positions:
    raise SystemExit(f"ERROR: TRACE ARG chr{chrom} has no sites: {tree}")
span_fraction = (max(positions) - min(positions)) / max(1.0, float(ts.sequence_length))
print(f"TRACE ARG scope chr{chrom}: sites={ts.num_sites} span_fraction={span_fraction:.6f} tree={tree}")
if ts.num_sites < min_sites or span_fraction < min_span:
    raise SystemExit(
        f"ERROR: TRACE requires a full-chromosome ARG for chr{chrom}; this tree looks loci-only "
        f"(sites={ts.num_sites}, span_fraction={span_fraction:.6f}, required sites>={min_sites}, span>={min_span}): {tree}"
    )
PY
}
write_manifest(){
  local mf=$OUT/manifest/tree_files.tsv c f idx; printf 'chr\tposterior\ttree_file\n' > "$mf"
  for c in $CHRS; do mapfile -t fs < <(find_trees "$c"); (( ${#fs[@]} )) || { echo "ERROR: no GU TRACE-format ARG tree for chr$c under $ARG_DIR. Build it with arg.sh build (add --format trace for needle), then run gu.sh trace." >&2; exit 1; }
    if printf '%s\n' "${fs[@]}" | grep -Eqi '(_part[0-9]+|chunk[._-]?[0-9]+)'; then echo "ERROR: coordinate-split ARG detected for chr$c; chunks cannot be posterior replicates" >&2; exit 1; fi
    idx=0; for f in "${fs[@]}"; do validate_tree_scope "$f" "$c"; idx=$((idx+1)); printf '%s\t%s\t%s\n' "$c" "$idx" "$f" >> "$mf"; done
  done
}
prepare_regions(){
  local c out; [[ -n $LOCI_FILE && $LOCI_MODE == extract ]] || return 0; [[ -s $LOCI_FILE ]] || { echo "ERROR: TRACE loci file missing: $LOCI_FILE" >&2; exit 1; }
  # trace-extract expects exactly three BED columns. The locus ID remains in
  # GU's staged/map files and is restored when final calls are clipped.
  for c in $CHRS; do out=$OUT/regions/chr$c.bed; awk -F'\t' -v c="$c" 'BEGIN{OFS="\t"}$1==c{print "chr"c,$2,$3}' "$LOCI_FILE" > "$out"; [[ -s $out ]] || { echo "ERROR: no loci remain for TRACE chr$c" >&2; exit 1; }; done
}
prepare_sample_map(){
  local tree chr vcf arg_map chr_map work first=1 validate_all=${TRACE_VALIDATE_ALL_TREES:-0}
  work=$OUT/samples/map_check; rm -rf "$work"; mkdir -p "$work" "$(dirname "$SAMPLE_MAP")"
  for chr in $CHRS; do
    vcf=
    tree=$(awk -F'\t' -v c="$chr" 'NR>1&&$1==c{print $3;exit}' "$OUT/manifest/tree_files.tsv")
    [[ -n $tree ]] || { echo "ERROR: manifest has no tree for chr$chr" >&2; return 1; }
    arg_map=$ARG_DIR/chr${chr}.sample_map.tsv; chr_map=$work/chr${chr}.sample_map.tsv
    [[ -n $TARGET_VCF_DIR && -s $TARGET_VCF_DIR/chr${chr}.vcf.gz ]] && vcf=$TARGET_VCF_DIR/chr${chr}.vcf.gz
    if [[ -s $arg_map ]]; then
      cp "$arg_map" "$chr_map"
    else
      python3 "$F/make_sample_map.py" --tree "$tree" --out "$chr_map"
    fi
    if [[ -n $vcf ]]; then
      bcftools query -l "$vcf" | sort -u > "$work/chr${chr}.target.samples"
      awk -F'\t' 'NR>1&&NF>=2{print $2}' "$chr_map" | sort -u > "$work/chr${chr}.arg.samples"
      if [[ $TRACE_ARG_METHOD == needle || $TRACE_ARG_METHOD == threads ]]; then
        comm -13 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" > "$work/chr${chr}.arg_not_target.samples"
        if [[ -s $work/chr${chr}.arg_not_target.samples ]]; then
          echo "ERROR: TRACE $TRACE_ARG_METHOD ARG contains chr$chr samples absent from the target cohort" >&2
          head -10 "$work/chr${chr}.arg_not_target.samples" >&2 || true
          return 1
        fi
        echo "TRACE $TRACE_ARG_METHOD panel sample check passed chr$chr: ARG=$(wc -l < "$work/chr${chr}.arg.samples") target=$(wc -l < "$work/chr${chr}.target.samples")"
      elif ! cmp -s "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples"; then
        echo "ERROR: TRACE target samples differ from the chr$chr ARG sample map; a target cohort cannot be projected onto a 1KG-only ARG" >&2
        echo "Only in target VCF:" >&2; comm -23 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" | head -10 >&2 || true
        echo "Only in ARG:" >&2; comm -13 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" | head -10 >&2 || true
        return 1
      fi
    fi
    if (( first )); then
      if [[ ! -s $SAMPLE_MAP ]]; then cp "$chr_map" "$SAMPLE_MAP"; fi
      first=0
    fi
    awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$SAMPLE_MAP" | sort > "$work/canonical.tsv"
    awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$chr_map" | sort > "$work/current.tsv"
    if ! cmp -s "$work/canonical.tsv" "$work/current.tsv"; then
      echo "ERROR: TRACE node/sample/haplotype map differs from the cached map for chr$chr; use a fresh chromosome/locus output or rebuild its cache" >&2
      diff -u "$work/canonical.tsv" "$work/current.tsv" | head -40 >&2 || true
      return 1
    fi
    if [[ $validate_all == 1 ]]; then
      while IFS=$'\t' read -r _ _ tree; do
        [[ $tree == tree_file ]] && continue
        python3 "$F/make_sample_map.py" --tree "$tree" --out "$work/posterior.sample_map.tsv"
        awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$work/posterior.sample_map.tsv" | sort > "$work/posterior.tsv"
        cmp -s "$work/canonical.tsv" "$work/posterior.tsv" || { echo "ERROR: TRACE posterior tree node map differs: $tree" >&2; return 1; }
      done < <(awk -F'\t' -v c="$chr" 'NR==1||$1==c' "$OUT/manifest/tree_files.tsv")
    fi
  done
  echo "TRACE sample-map check passed for chromosomes: $CHRS"
}
split_groups(){
  if [[ -n $TARGET_SAMPLES ]]; then [[ -s $TARGET_SAMPLES ]] || { echo "ERROR: TRACE_TARGET_SAMPLES missing: $TARGET_SAMPLES" >&2; exit 1; }; awk 'BEGIN{FS="\t"}NR==FNR{if(NF)want[$1]=1;next}FNR>1&&($2 in want){print $1}' "$TARGET_SAMPLES" "$SAMPLE_MAP" > "$OUT/samples/tree_nodes.txt"
  else awk 'BEGIN{FS="\t"}NR>1{print $1}' "$SAMPLE_MAP" > "$OUT/samples/tree_nodes.txt"; fi
  [[ -s $OUT/samples/tree_nodes.txt ]] || { echo "ERROR: no TRACE target nodes" >&2; exit 1; }; rm -rf "$OUT/samples/groups"; mkdir -p "$OUT/samples/groups"
  awk -v n="$GROUP_SIZE" -v d="$OUT/samples/groups" 'NF{g=int((NR-1)/n)+1;f=sprintf("%s/group%05d.nodes",d,g);print $1>>f}' "$OUT/samples/tree_nodes.txt"
  : > "$OUT/samples/node_group.tsv"; for g in "$OUT"/samples/groups/*.nodes; do b=$(basename "$g" .nodes); awk -v b="$b" '{print $1"\t"b}' "$g" >> "$OUT/samples/node_group.tsv"; done
}
run_record(){
  local c post tree cmd
  printf 'build\t%s\nchromosomes\t%s\nloci_mode\t%s\nt_archaic\t%s\nposterior_threshold\t%s\nphysical_length_threshold\t%s\ngenetic_distance_threshold\t%s\nwindow_size\t%s\ngroup_size\t%s\n' \
    "$BUILD" "$CHRS" "$LOCI_MODE" "$TARC" "$POST" "$MINBP" "$MINCM" "$WIN" "$GROUP_SIZE"
  [[ -z $LOCI_FILE ]] || stat -c 'loci\t%n:%s:%Y' "$LOCI_FILE"
  stat -c 'sample_map\t%n:%s:%Y' "$SAMPLE_MAP"
  [[ -z $TARGET_SAMPLES ]] || stat -c 'target_samples\t%n:%s:%Y' "$TARGET_SAMPLES"
  while IFS=$'\t' read -r c post tree; do [[ $c == chr ]] && continue; stat -c "tree_chr${c}_p${post}\\t%n:%s:%Y" "$tree"; done < "$OUT/manifest/tree_files.tsv"
  for cmd in trace-extract trace-infer trace-summarize; do
    command -v "$cmd" >/dev/null 2>&1 && stat -c 'software\t%n:%s:%Y' "$(command -v "$cmd")"
  done
  stat -c 'software\t%n:%s:%Y' "$F/trace.sh" "$F/trace_combine.py" "$F/make_sample_map.py"
}
check_run_provenance(){
  local current=$OUT/run.meta.tsv
  [[ -s $current ]] || run_record > "$current"
}
prepare_run(){ write_manifest; prepare_regions; prepare_sample_map; split_groups; check_run_provenance; }
extract_one(){
  local c=$1 gf=$2 group=$3 post=$4 tree=$5 prefix ids npost; local -a opts=()
  prefix=$OUT/extract/${group}.chr${c}.p${post}; [[ -s $prefix.npz ]] && return 0; ids=$(paste -sd, "$gf"); npost=$(awk -v c="$c" 'NR>1&&$1==c{n++}END{print n+0}' "$OUT/manifest/tree_files.tsv")
  (( npost > 1 )) && opts+=(--window-size "$WIN")
  [[ -n $LOCI_FILE && $LOCI_MODE == extract ]] && opts+=(--include-regions "$OUT/regions/chr$c.bed" --chrom "chr$c")
  trace-extract --tree-file "$tree" -t "$TARC" --individuals "$ids" "${opts[@]}" -o "$prefix" > "$OUT/log/extract.${group}.chr${c}.p${post}.log" 2>&1
}
run_extract(){
  local running=0 status=0 gf group c post tree
  for gf in "$OUT"/samples/groups/*.nodes; do group=$(basename "$gf" .nodes); while IFS=$'\t' read -r c post tree; do [[ $c == chr ]] && continue; extract_one "$c" "$gf" "$group" "$post" "$tree" & running=$((running+1)); if (( running >= JEX )); then wait -n || status=1; running=$((running-1)); fi; done < "$OUT/manifest/tree_files.tsv"; done
  while (( running )); do wait -n || status=1; running=$((running-1)); done; (( status == 0 )) || exit 1
  for gf in "$OUT"/samples/groups/*.nodes; do group=$(basename "$gf" .nodes); for c in $CHRS; do find "$OUT/extract" -maxdepth 1 -name "${group}.chr${c}.p*.npz" -print | sort > "$OUT/datafiles/${group}.chr${c}.txt"; done; done
}
gmap_for_chr(){ local c=$1 d=${TRACE_GENETIC_MAP_DIR:-}; [[ -n $d && -d $d ]] || return 1; find "$d" -maxdepth 1 -type f \( -iname "*chr${c}*" -o -iname "*_${c}.*" \) -print | sort | head -1; }
infer_one(){
  local node=$1 group=$2 c npost gmap C D G out any_multi=0 complete=1; local -a chroms=() datafiles=() npzs=() gmaps=() opt=()
  out=$OUT/infer/hap$node
  for c in $CHRS; do chroms+=("chr$c"); npost=$(wc -l < "$OUT/datafiles/${group}.chr${c}.txt"); (( npost > 1 )) && any_multi=1; datafiles+=("$OUT/datafiles/${group}.chr${c}.txt"); npzs+=("$(head -1 "$OUT/datafiles/${group}.chr${c}.txt")"); if gmap=$(gmap_for_chr "$c"); then gmaps+=("$gmap"); fi; done
  # A prior run may have produced only the final chromosome before failing.
  # Reuse is valid only when every requested chromosome output exists.
  for c in $CHRS; do [[ -s $out.chr${c}.xss.npz ]] || complete=0; done
  (( complete )) && return 0
  C=$(IFS=,;echo "${chroms[*]}")
  if (( any_multi )); then D=$(IFS=,;echo "${datafiles[*]}"); opt=(--data-files "$D"); else D=$(IFS=,;echo "${npzs[*]}"); opt=(--npz-files "$D"); fi
  if (( ${#gmaps[@]} == ${#chroms[@]} )); then G=$(IFS=,;echo "${gmaps[*]}"); opt+=(--genetic-maps "$G"); fi
  trace-infer -i "$node" "${opt[@]}" --chroms "$C" -o "$out" > "$OUT/log/infer.hap${node}.log" 2>&1
}
run_infer(){ [[ -s $OUT/samples/node_group.tsv ]] || run_extract; local running=0 status=0 node group; while IFS=$'\t' read -r node group; do infer_one "$node" "$group" & running=$((running+1)); if (( running >= JIN )); then wait -n || status=1; running=$((running-1)); fi; done < "$OUT/samples/node_group.tsv"; while (( running )); do wait -n || status=1; running=$((running-1)); done; (( status == 0 )) || exit 1; }
summarize_one(){ local node=$1 c f out FSTR CSTR; local -a files=() chroms=(); out=$OUT/calls/hap$node; [[ -s $out.summary.txt ]] && return 0; for c in $CHRS; do f=$OUT/infer/hap${node}.chr${c}.xss.npz; [[ -s $f ]] || return 1; files+=("$f"); chroms+=("chr$c"); done; FSTR=$(IFS=,;echo "${files[*]}"); CSTR=$(IFS=,;echo "${chroms[*]}"); trace-summarize -f "$FSTR" -c "$CSTR" --posterior-threshold "$POST" --physical-length-threshold "$MINBP" --genetic-distance-threshold "$MINCM" -o "$out" > "$OUT/log/summarize.hap${node}.log" 2>&1; }
run_segments(){
  [[ -s $OUT/samples/tree_nodes.txt ]] || run_extract; local running=0 status=0 node; while read -r node; do summarize_one "$node" & running=$((running+1)); if (( running >= JSUM )); then wait -n || status=1; running=$((running-1)); fi; done < "$OUT/samples/tree_nodes.txt"; while (( running )); do wait -n || status=1; running=$((running-1)); done; (( status == 0 )) || exit 1
  local -a args=(--root "$OUT" --sample-map "$SAMPLE_MAP" --build "$BUILD"); [[ -n $LOCI_FILE ]] && args+=(--loci "$LOCI_FILE"); python3 "$F/trace_combine.py" "${args[@]}"
}
prepare_run
case "$ACTION" in
  run) run_extract; run_infer; run_segments ;;
  extract) run_extract ;;
  infer) run_infer; run_segments ;;
  segments) run_segments ;;
  check) echo "TRACE CHECK PASSED: manifest=$OUT/manifest/tree_files.tsv sample_map=$SAMPLE_MAP" ;;
esac

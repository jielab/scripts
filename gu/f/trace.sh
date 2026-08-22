#!/usr/bin/env bash
# Standalone TRACE driver. It intentionally does not share the IBDmix engine.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd); F=$ROOT/f
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }; }
for x in python3 trace-extract trace-infer trace-summarize awk sort find gzip; do need "$x"; done
ACTION=${1:-run}; shift || true
ARG_DIR=${GU_ARG_DIR:?set GU_ARG_DIR}; TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-}; OUT=${TRACE_OUT:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu}/trace}
BUILD=GRCh${GU_BUILD:-37}; CHRS=${TRACE_CHRS:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}}
GROUP_SIZE=${TRACE_GROUP_SIZE:-50}; JEX=${TRACE_JOB_EXTRACT:-1}; JIN=${TRACE_JOB_INFER:-8}; JSUM=${TRACE_JOB_SUMMARIZE:-8}
TARC=${TRACE_T_ARCHAIC:-15000}; POST=${TRACE_POSTERIOR_THRESHOLD:-0.90}; MINBP=${TRACE_PHYSICAL_LENGTH_THRESHOLD:-50000}; MINCM=${TRACE_GENETIC_DISTANCE_THRESHOLD:-0.05}; WIN=${TRACE_WINDOW_SIZE:-1000}
SAMPLE_MAP=${TRACE_SAMPLE_MAP:-$OUT/samples/trace_sample_map.tsv}; TARGET_SAMPLES=${TRACE_TARGET_SAMPLES:-}
mkdir -p "$OUT"/{manifest,samples/groups,extract,datafiles,infer,summary,final,report,log}
usage(){ cat <<USAGE
Usage: ./gu.sh trace [run|extract|infer|segments|report]
Environment: GU_ARG_DIR, GU_TARGET_VCF_DIR, TRACE_CHRS, TRACE_SAMPLE_MAP.
Optional TRACE_TARGET_SAMPLES is a one-column list of target people to call. This is
recommended for UKB anchor+target ARGs: anchors shape the genealogy but are not called.
The default is chromosomes 1..22 and X. Explicit TRACE_CHRS or GU_CHRS values
are used unchanged.
USAGE
}
[[ $ACTION == help || $ACTION == -h || $ACTION == --help ]] && { usage; exit 0; }
find_trees(){ local c=$1; find "$ARG_DIR" -type f \( -name '*.trees' -o -name '*.tsz' \) -print | awk -v c="$c" 'BEGIN{IGNORECASE=1}{b=$0;gsub(/.*\//,"",b);if(b~("(^|[^A-Za-z0-9])chr"c"([^A-Za-z0-9]|$)"))print $0;else if(b~("(^|[^A-Za-z0-9])"c"([^A-Za-z0-9]|$)"))print $0}' | sort; }
write_manifest(){ local mf=$OUT/manifest/tree_files.tsv c f idx n; printf 'chr\tposterior\ttree_file\n' > "$mf"; for c in $CHRS; do mapfile -t fs < <(find_trees "$c"); n=${#fs[@]}; ((n>0)) || { echo "ERROR: no ARG tree for chr$c under $ARG_DIR" >&2; exit 1; }; if printf '%s\n' "${fs[@]}"|grep -Eqi '(_part[0-9]+|chunk[._-]?[0-9]+)'; then echo "ERROR coordinate-split ARG detected for chr$c; coordinate chunks are not posterior replicates." >&2; exit 1; fi; idx=0; for f in "${fs[@]}"; do idx=$((idx+1)); printf '%s\t%s\t%s\n' "$c" "$idx" "$f" >> "$mf"; done; done; }
prepare_sample_map(){ [[ -s $SAMPLE_MAP ]] && return 0; local first_tree first_chr vcf; first_tree=$(awk 'NR==2{print $3}' "$OUT/manifest/tree_files.tsv"); first_chr=$(awk 'NR==2{print $1}' "$OUT/manifest/tree_files.tsv"); [[ -n $TARGET_VCF_DIR ]] || { echo "ERROR TRACE_SAMPLE_MAP absent and GU_TARGET_VCF_DIR unset" >&2; exit 1; }; vcf=$TARGET_VCF_DIR/chr${first_chr}.vcf.gz; [[ -s $vcf ]] || { echo "ERROR missing $vcf" >&2; exit 1; }; python3 "$F/make_sample_map.py" --tree "$first_tree" --vcf "$vcf" --out "$SAMPLE_MAP"; }
split_groups(){
  if [[ -n $TARGET_SAMPLES ]]; then
    [[ -s $TARGET_SAMPLES ]] || { echo "ERROR TRACE_TARGET_SAMPLES missing: $TARGET_SAMPLES" >&2; exit 1; }
    awk 'BEGIN{FS="\t"} NR==FNR{if(NF)want[$1]=1;next} FNR>1 && ($2 in want){print $1}' "$TARGET_SAMPLES" "$SAMPLE_MAP" > "$OUT/samples/tree_nodes.txt"
    npeople=$(awk 'NF{n++}END{print n+0}' "$TARGET_SAMPLES"); nnodes=$(wc -l < "$OUT/samples/tree_nodes.txt"); echo "TRACE target filter: people=$npeople nodes=$nnodes"
    ((nnodes>0)) || { echo "ERROR none of TRACE_TARGET_SAMPLES mapped to ARG" >&2; exit 1; }
  else awk 'BEGIN{FS="\t"} NR>1{print $1}' "$SAMPLE_MAP" > "$OUT/samples/tree_nodes.txt"; fi
  rm -rf "$OUT/samples/groups"; mkdir -p "$OUT/samples/groups"; awk -v n="$GROUP_SIZE" -v d="$OUT/samples/groups" 'NF{g=int((NR-1)/n)+1;f=sprintf("%s/group%05d.nodes",d,g);print $1>>f}' "$OUT/samples/tree_nodes.txt"
  : > "$OUT/samples/node_group.tsv"; for g in "$OUT"/samples/groups/*.nodes; do [[ -e $g ]] || continue; bn=$(basename "$g" .nodes); awk -v b="$bn" '{print $1"\t"b}' "$g" >> "$OUT/samples/node_group.tsv"; done
}
extract_one(){ local c=$1 gf=$2 gname=$3 post=$4 tree=$5 prefix npz opts=(); prefix=$OUT/extract/${gname}.chr${c}.p${post}; npz=$prefix.npz; [[ -s $npz ]]&&return 0; ids=$(paste -sd, "$gf"); npost=$(awk -v c="$c" 'NR>1&&$1==c{n++}END{print n+0}' "$OUT/manifest/tree_files.tsv"); ((npost>1))&&opts+=(--window-size "$WIN"); trace-extract --tree-file "$tree" -t "$TARC" --individuals "$ids" "${opts[@]}" -o "$prefix" > "$OUT/log/extract.${gname}.chr${c}.p${post}.log" 2>&1; }
run_extract(){ write_manifest; prepare_sample_map; split_groups; local running=0 status=0 gf gname c post tree; for gf in "$OUT"/samples/groups/*.nodes; do gname=$(basename "$gf" .nodes); while IFS=$'\t' read -r c post tree; do [[ $c == chr ]]&&continue; extract_one "$c" "$gf" "$gname" "$post" "$tree" & running=$((running+1)); if ((running>=JEX)); then wait -n||status=1; running=$((running-1)); fi; done < "$OUT/manifest/tree_files.tsv"; done; while ((running>0)); do wait -n||status=1; running=$((running-1)); done; ((status==0))||exit 1; for gf in "$OUT"/samples/groups/*.nodes; do gname=$(basename "$gf" .nodes); for c in $CHRS; do ls "$OUT/extract/${gname}.chr${c}.p"*.npz 2>/dev/null|sort > "$OUT/datafiles/${gname}.chr${c}.txt"; done; done; }
gmap_for_chr(){ local c=$1 d=${TRACE_GENETIC_MAP_DIR:-}; [[ -n $d && -d $d ]]||return 1; find "$d" -maxdepth 1 -type f \( -iname "*chr${c}*" -o -iname "*_${c}.*" \) -print|sort|head -1; }
infer_one(){ local node=$1 group=$2 chroms=() datafiles=() npzs=() gmaps=() c npost gmap outprefix=$OUT/infer/hap$node any_multi=0; for c in $CHRS; do chroms+=("chr$c"); npost=$(wc -l < "$OUT/datafiles/${group}.chr${c}.txt"); ((npost>1))&&any_multi=1; datafiles+=("$OUT/datafiles/${group}.chr${c}.txt"); npzs+=("$(head -1 "$OUT/datafiles/${group}.chr${c}.txt")"); if gmap=$(gmap_for_chr "$c"); then gmaps+=("$gmap"); fi; done; lastc=${chroms[$((${#chroms[@]}-1))]}; [[ -s "$outprefix.${lastc}.xss.npz" ]]&&return 0; C=$(IFS=,;echo "${chroms[*]}"); if ((any_multi)); then D=$(IFS=,;echo "${datafiles[*]}"); opt=(--data-files "$D"); else D=$(IFS=,;echo "${npzs[*]}"); opt=(--npz-files "$D"); fi; if ((${#gmaps[@]}==${#chroms[@]})); then G=$(IFS=,;echo "${gmaps[*]}"); opt+=(--genetic-maps "$G"); fi; trace-infer -i "$node" "${opt[@]}" --chroms "$C" -o "$outprefix" > "$OUT/log/infer.hap${node}.log" 2>&1; }
run_infer(){ [[ -s $OUT/samples/node_group.tsv ]]||run_extract; local running=0 status=0 node group; while IFS=$'\t' read -r node group; do infer_one "$node" "$group" & running=$((running+1)); if ((running>=JIN)); then wait -n||status=1; running=$((running-1)); fi; done < "$OUT/samples/node_group.tsv"; while ((running>0)); do wait -n||status=1; running=$((running-1)); done; ((status==0))||exit 1; }
summarize_one(){ local node=$1 c files=() chroms=() f out=$OUT/summary/hap$node; [[ -s $out.summary.txt ]]&&return 0; for c in $CHRS; do f=$OUT/infer/hap${node}.chr${c}.xss.npz; [[ -s $f ]]||return 1; files+=("$f"); chroms+=("chr$c"); done; FSTR=$(IFS=,;echo "${files[*]}"); CSTR=$(IFS=,;echo "${chroms[*]}"); trace-summarize -f "$FSTR" -c "$CSTR" --posterior-threshold "$POST" --physical-length-threshold "$MINBP" --genetic-distance-threshold "$MINCM" -o "$out" > "$OUT/log/summarize.hap${node}.log" 2>&1; }
run_segments(){ [[ -s $OUT/samples/tree_nodes.txt ]]||run_extract; local running=0 status=0 node; while read -r node; do summarize_one "$node" & running=$((running+1)); if ((running>=JSUM)); then wait -n||status=1; running=$((running-1)); fi; done < "$OUT/samples/tree_nodes.txt"; while ((running>0)); do wait -n||status=1; running=$((running-1)); done; ((status==0))||exit 1; python3 "$F/trace_combine.py" --root "$OUT" --sample-map "$SAMPLE_MAP" --build "$BUILD"; }
case "$ACTION" in run)run_extract;run_infer;run_segments;; extract)run_extract;; infer)run_infer;run_segments;; segments)run_segments;; report)python3 "$F/trace_combine.py" --root "$OUT" --sample-map "$SAMPLE_MAP" --build "$BUILD";; *)usage;exit 2;; esac

#!/usr/bin/env bash

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Internal TRACE driver; public entry: gu.sh trace

Main options:
    Use ./gu.sh trace with --chr, --dir-arg and target options.
    -h, --help  Show help

Example:
    ./gu.sh -h
HELP
			exit 0
			;;
	esac
fi


# 🚩 trace
# TRACE computation driver. Shiny owns reporting and aggregation.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
F=$ROOT/f
source "$F/0.common.sh"
gu_load_config
need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERROR: missing command: $1" >&2
	exit 1
}; }
ACTION=${1:-run}
case "$ACTION" in run | extract | infer | segments | check) ;; *)
	echo "ERROR: invalid TRACE action: $ACTION" >&2
	exit 2
	;;
esac
ARG_DIR=${GU_ARG_DIR:-${GU_TARGET_ROOT:-.}/arg.tsinfer}
TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-}
TRACE_RESULT_MODULE=trace
[[ ${TRACE_BACKEND:-cpu} != cuda ]] || TRACE_RESULT_MODULE=trace-GPU
OUT=${TRACE_OUT:-${GU_ANALYSIS_ROOT:-$GU_ANALYSIS_BASE/1kg}/${TRACE_RESULT_MODULE:-trace}/${GU_SCOPE_ID:-genome}}
BUILD=GRCh${GU_BUILD:-37}
CHRS=${TRACE_CHRS:-${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}}
LOCI_FILE=${TRACE_LOCI_FILE:-${GU_LOCI_FILE:-}}
LOCI_MODE=${TRACE_LOCI_MODE:-posthoc}
case "$LOCI_MODE" in posthoc | extract) ;; *)
	echo "ERROR: internal trace-loci-mode must be posthoc or extract" >&2
	exit 2
	;;
esac
GROUP_SIZE=${TRACE_GROUP_SIZE:-50}
JEX=${TRACE_JOB_EXTRACT:-1}
JIN=${TRACE_JOB_INFER:-8}
JSUM=${TRACE_JOB_SUMMARIZE:-8}
TARC=${TRACE_T_ARCHAIC:-15000}
POST=${TRACE_POSTERIOR_THRESHOLD:-0.90}
MINBP=${TRACE_PHYSICAL_LENGTH_THRESHOLD:-50000}
MINCM=${TRACE_GENETIC_DISTANCE_THRESHOLD:-0.05}
WIN=${TRACE_WINDOW_SIZE:-1000}
SAMPLE_MAP=${TRACE_SAMPLE_MAP:-$OUT/samples/trace_sample_map.tsv}
TARGET_SAMPLES=${TRACE_TARGET_SAMPLES:-}
FINAL_NAME=${TRACE_FINAL_NAME:-trace_haplotype_segments.tsv.gz}
export TRACE_BACKEND=${TRACE_BACKEND:-cpu}
TRACE_SEED=${TRACE_SEED:-42}
INFER_CMD=(python3 "$F/trace.py" infer)
case "$TRACE_BACKEND" in
	cpu) ;;
	cuda)
		python3 "$F/trace.py" infer --gpu-check
		# One device-owning inference process; extraction and summary stay CPU.
		JIN=1
		;;
	*) echo "ERROR: unknown TRACE backend: $TRACE_BACKEND" >&2; exit 2 ;;
esac
TRACE_ARG_METHOD=$(awk -F'\t' '$1=="method"{print $2;exit}' "$ARG_DIR/ARG_TRACE_BUILD.tsv" 2>/dev/null || true)
for value in "$GROUP_SIZE" "$JEX" "$JIN" "$JSUM"; do
	[[ $value =~ ^[1-9][0-9]*$ ]] || {
		echo "ERROR: TRACE group size and worker counts must be positive integers" >&2
		exit 2
	}
done
# One native-library thread per worker; parallelism comes from the job counts.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
REPLACE=${TRACE_REPLACE:-0}
[[ $REPLACE == 0 || $REPLACE == 1 ]] || {
	echo "ERROR: internal replace-trace state must be 0 or 1" >&2
	exit 2
}
TRACE_FINAL=$OUT/final/$FINAL_NAME
python3 - "$TARC" "$POST" "$MINBP" "$MINCM" "$WIN" "$TRACE_SEED" <<'PYVALID'
import math,sys
t,p,b,c,w,seed=map(float,sys.argv[1:])
if not all(map(math.isfinite,(t,p,b,c,w,seed))) or not (t>0 and 0<p<=1 and b>=0 and c>=0 and w>=1 and seed>=0 and b.is_integer() and w.is_integer() and seed.is_integer()):
    raise SystemExit('ERROR: invalid TRACE thresholds/window/seed')
PYVALID
for x in python3 awk sort find gzip stat cmp grep diff comm flock pgrep; do need "$x"; done
if [[ $ACTION != check ]]; then for x in trace-extract trace-infer trace-summarize; do need "$x"; done; fi
TRACE_LOG_DIR=${DIR_RUNTIME_TMP}/gu-logs/trace/$(printf '%s' "$OUT" | sha256sum | cut -c1-16)
mkdir -p "$OUT"/{manifest,samples/groups,regions,maps,extract,datafiles,infer,calls,final} "$TRACE_LOG_DIR"
# A persistent lock in the shared result directory also protects jobs on
# different HPC nodes. Never unlink a live flock inode.
exec 9>"$OUT/.trace.lock"
flock -n 9 || {
	echo "ERROR: another TRACE run is using $OUT" >&2
	exit 1
}
trace_children=()
trace_stop_tree() {
	local pid=$1 child
	while read -r child; do [[ -z $child ]] || trace_stop_tree "$child"; done < <(pgrep -P "$pid" 2>/dev/null || true)
	kill -TERM "$pid" 2>/dev/null || true
}
trace_exit() {
	local status=$? pid
	trap - EXIT
	if ((status != 0)); then
		for pid in "${trace_children[@]}"; do trace_stop_tree "$pid"; done
		for pid in "${trace_children[@]}"; do wait "$pid" 2>/dev/null || true; done
	fi
	exit "$status"
}
trap trace_exit EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

# Autosome and male non-PAR X sample nodes are different biological copies.
# Keep joint autosomal HMM training, and fit X separately with its own map.
# Internal final filenames are deliberately distinct to avoid double import
# by the recursive Shiny result collector.
autosomes=()
has_x=0
for c in $CHRS; do
	if [[ $c == X ]]; then has_x=1; else autosomes+=("$c"); fi
done
if ((has_x && ${#autosomes[@]})); then
	[[ $ACTION == check ]] || rm -f -- "$TRACE_FINAL"
	for class in autosomes chrX; do
		if [[ $class == autosomes ]]; then class_chrs="${autosomes[*]}"; else class_chrs=X; fi
		env TRACE_CHRS="$class_chrs" TRACE_OUT="$OUT/chromosome_classes/$class" TRACE_SAMPLE_MAP= \
			TRACE_FINAL_NAME=trace_haplotype_segments.internal.tsv.gz bash "$F/trace.sh" "$ACTION" &
		trace_children=("$!")
		wait "${trace_children[0]}"
		trace_children=()
	done
	if [[ $ACTION == run || $ACTION == infer || $ACTION == segments ]]; then
		python3 - "$OUT" "$FINAL_NAME" <<'PYSPLIT'
from pathlib import Path
import pandas as pd, sys
root=Path(sys.argv[1])
frames=[]
for kind in ('autosomes','chrX'):
    f=root/'chromosome_classes'/kind/'final/trace_haplotype_segments.internal.tsv.gz'
    d=pd.read_csv(f,sep='\t',dtype={'sample':str,'tree_node_id':str})
    d['trace_chromosome_class']=kind
    frames.append(d)
stage=root/'final'/(sys.argv[2]+'.part')
pd.concat(frames,ignore_index=True).to_csv(stage,sep='\t',index=False,compression='gzip')
stage.replace(root/'final'/sys.argv[2])
PYSPLIT
	fi
	{
		printf 'build\t%s\nchromosomes\t%s\nchromosome_classes\tautosomes,chrX\nbackend\t%s\n' "$BUILD" "$CHRS" "$TRACE_BACKEND"
		printf 'runtime\tgu-trace-split-nonPAR-X-v2\n'
		for class in autosomes chrX; do sha256sum "$OUT/chromosome_classes/$class/run.meta.tsv"; done
	} >"$OUT/run.meta.tsv.next"
	mv -f "$OUT/run.meta.tsv.next" "$OUT/run.meta.tsv"
	echo "TRACE $ACTION completed with joint autosomes and a separate male non-PAR X model: $OUT"
	exit 0
fi
if [[ $ACTION != check && $REPLACE == 1 ]]; then
	rm -rf -- "$OUT/extract" "$OUT/datafiles" "$OUT/infer" "$OUT/calls" "$OUT/final"
	mkdir -p "$OUT"/{extract,datafiles,infer,calls,final}
	if [[ $SAMPLE_MAP == "$OUT/samples/trace_sample_map.tsv" ]]; then rm -f -- "$SAMPLE_MAP"; fi
fi

# GU TRACE consumes mutation-bearing, TRACE-formatted trees stored directly in
# ARG_DIR. Native method artifacts stay in their method directories and are
# exported explicitly with `arg.sh build --format trace`.
find_trees() {
	local c=$1
	find -L "$ARG_DIR" -maxdepth 1 -type f \( -name '*.trees' -o -name '*.tsz' \) -print | awk -v c="$c" 'BEGIN{IGNORECASE=1}{b=$0;gsub(/.*\//,"",b);if(b~("(^|[^A-Za-z0-9])chr"c"([^A-Za-z0-9]|$)"))print;else if(b~("(^|[^A-Za-z0-9])"c"([^A-Za-z0-9]|$)"))print}' | sort
}
validate_tree_scope() {
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
if ts.time_units != 'generations':
    raise SystemExit(f'ERROR: TRACE requires ARG node times in generations; got {ts.time_units!r}: {tree}')
if isinstance(ts.metadata, dict) and ts.metadata.get('method') == 'singer' and ts.metadata.get('regional'):
    raise SystemExit('ERROR: regional SINGER test ARG cannot be used as full-chromosome TRACE evidence')
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
write_manifest() {
	local mf=$OUT/manifest/tree_files.tsv c f idx
	printf 'chr\tposterior\ttree_file\n' >"$mf"
	for c in $CHRS; do
		mapfile -t fs < <(find_trees "$c")
		((${#fs[@]})) || {
			echo "ERROR: no GU TRACE-format ARG tree for chr$c under $ARG_DIR. Build it with arg.sh build (add --format trace for needle), then run gu.sh trace." >&2
			exit 1
		}
		if printf '%s\n' "${fs[@]}" | grep -Eqi '(_part[0-9]+|chunk[._-]?[0-9]+)'; then
			echo "ERROR: coordinate-split ARG detected for chr$c; chunks cannot be posterior replicates" >&2
			exit 1
		fi
		idx=0
		for f in "${fs[@]}"; do
			validate_tree_scope "$f" "$c"
			idx=$((idx + 1))
			printf '%s\t%s\t%s\n' "$c" "$idx" "$f" >>"$mf"
		done
	done
}
prepare_regions() {
	local c tree
	for c in $CHRS; do
		tree=$(awk -F'\t' -v c="$c" 'NR>1&&$1==c{print $3;exit}' "$OUT/manifest/tree_files.tsv")
		python3 - "$TRACE_ARG_METHOD" "$ARG_DIR" "$LOCI_FILE" "$LOCI_MODE" "$c" "${GU_BUILD:-37}" "$tree" "$OUT/regions/chr$c.bed" "${TRACE_INCLUDE_REGIONS:-}" <<'PYMASK'
import sys
from pathlib import Path
method,root,loci,mode,c,build,tree,out,include=sys.argv[1:]
def read(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        r=line.split()
        if r[0].removeprefix('chr')==c:
            a,b=int(r[1]),int(r[2])
            if a<0 or b<=a:raise SystemExit(f'ERROR: invalid BED interval in {path}: {line}')
            rows.append((a,b))
    return rows
def intersect(a,b):
    return [(max(x,l),min(y,r)) for x,y in a for l,r in b if max(x,l)<min(y,r)]
regions=None
if method=='singer':regions=read(Path(root)/f'chr{c}.inferred.bed')
if c=='X':
    start,end={'37':(2699520,154931043),'38':(2781479,155701382)}[build]
    x=[(start,end)]
    regions=x if regions is None else intersect(regions,x)
if include:
    mask=read(include)
    regions=mask if regions is None else intersect(regions,mask)
if mode=='extract' and loci:
    mask=read(loci)
    regions=mask if regions is None else intersect(regions,mask)
if regions is None:
    Path(out).unlink(missing_ok=True)
else:
    if not regions:raise SystemExit(f'ERROR: no accessible TRACE intervals for chr{c}')
    merged=[]
    for a,b in sorted(regions):
        if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
        else:merged.append((a,b))
    Path(out).write_text(''.join(f'chr{c}\t{a}\t{b}\n' for a,b in merged))
PYMASK
	done
}

prepare_sample_map() {
	local tree chr vcf arg_map chr_map work first=1 validate_all=${TRACE_VALIDATE_ALL_TREES:-1}
	work=$OUT/samples/map_check
	rm -rf "$work"
	mkdir -p "$work" "$(dirname "$SAMPLE_MAP")"
	for chr in $CHRS; do
		vcf=
		tree=$(awk -F'\t' -v c="$chr" 'NR>1&&$1==c{print $3;exit}' "$OUT/manifest/tree_files.tsv")
		[[ -n $tree ]] || {
			echo "ERROR: manifest has no tree for chr$chr" >&2
			return 1
		}
		arg_map=$ARG_DIR/chr${chr}.sample_map.tsv
		chr_map=$work/chr${chr}.sample_map.tsv
		if [[ $chr == X && -n $TARGET_VCF_DIR && -s $TARGET_VCF_DIR/chrX.male.vcf.gz ]]; then
			vcf=$TARGET_VCF_DIR/chrX.male.vcf.gz
		elif [[ -n $TARGET_VCF_DIR && -s $TARGET_VCF_DIR/chr${chr}.vcf.gz ]]; then
			vcf=$TARGET_VCF_DIR/chr${chr}.vcf.gz
		fi
		if [[ -s $arg_map ]]; then
			cp "$arg_map" "$chr_map"
		else
			python3 "$F/0.common.py" sample-map --tree "$tree" --out "$chr_map"
		fi
		if [[ -n $vcf ]]; then
			need bcftools
			bcftools query -l "$vcf" | sort -u >"$work/chr${chr}.target.samples"
			awk -F'\t' 'NR>1&&NF>=2{print $2}' "$chr_map" | sort -u >"$work/chr${chr}.arg.samples"
			if [[ $TRACE_ARG_METHOD == needle || $TRACE_ARG_METHOD == threads || $TRACE_ARG_METHOD == singer ]]; then
				comm -13 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" >"$work/chr${chr}.arg_not_target.samples"
				if [[ -s $work/chr${chr}.arg_not_target.samples ]]; then
					echo "ERROR: TRACE $TRACE_ARG_METHOD ARG contains chr$chr samples absent from the target cohort" >&2
					head -10 "$work/chr${chr}.arg_not_target.samples" >&2 || true
					return 1
				fi
				echo "TRACE $TRACE_ARG_METHOD panel sample check passed chr$chr: ARG=$(wc -l <"$work/chr${chr}.arg.samples") target=$(wc -l <"$work/chr${chr}.target.samples")"
			elif ! cmp -s "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples"; then
				echo "ERROR: TRACE target samples differ from the chr$chr ARG sample map; a target cohort cannot be projected onto a 1KG-only ARG" >&2
				echo "Only in target VCF:" >&2
				comm -23 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" | head -10 >&2 || true
				echo "Only in ARG:" >&2
				comm -13 "$work/chr${chr}.target.samples" "$work/chr${chr}.arg.samples" | head -10 >&2 || true
				return 1
			fi
		fi
		if ((first)); then
			if [[ ! -s $SAMPLE_MAP ]]; then cp "$chr_map" "$SAMPLE_MAP"; fi
			first=0
		fi
		awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$SAMPLE_MAP" | sort >"$work/canonical.tsv"
		awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$chr_map" | sort >"$work/current.tsv"
		if ! cmp -s "$work/canonical.tsv" "$work/current.tsv"; then
			echo "ERROR: TRACE node/sample/haplotype map differs from the cached map for chr$chr; use a fresh chromosome/locus output or rebuild its cache" >&2
			diff -u "$work/canonical.tsv" "$work/current.tsv" | head -40 >&2 || true
			return 1
		fi
		if [[ $validate_all == 1 ]]; then
			while IFS=$'\t' read -r _ _ tree; do
				[[ $tree == tree_file ]] && continue
				python3 "$F/trace.py" validate-map --tree "$tree" --sample-map "$SAMPLE_MAP"

			done < <(awk -F'\t' -v c="$chr" 'NR==1||$1==c' "$OUT/manifest/tree_files.tsv")
		fi
	done
	echo "TRACE sample-map check passed for chromosomes: $CHRS"
}
split_groups() {
	if [[ -n $TARGET_SAMPLES ]]; then
		[[ -s $TARGET_SAMPLES ]] || {
			echo "ERROR: TRACE_TARGET_SAMPLES missing: $TARGET_SAMPLES" >&2
			exit 1
		}
		awk 'BEGIN{FS="\t"}NR==FNR{if(NF)want[$1]=1;next}FNR>1&&($2 in want){print $1}' "$TARGET_SAMPLES" "$SAMPLE_MAP" >"$OUT/samples/tree_nodes.txt"
	else awk 'BEGIN{FS="\t"}NR>1{print $1}' "$SAMPLE_MAP" >"$OUT/samples/tree_nodes.txt"; fi
	[[ -s $OUT/samples/tree_nodes.txt ]] || {
		echo "ERROR: no TRACE target nodes" >&2
		exit 1
	}
	rm -rf "$OUT/samples/groups"
	mkdir -p "$OUT/samples/groups"
	awk -v n="$GROUP_SIZE" -v d="$OUT/samples/groups" 'NF{g=int((NR-1)/n)+1;f=sprintf("%s/group%05d.nodes",d,g);print $1>>f}' "$OUT/samples/tree_nodes.txt"
	: >"$OUT/samples/node_group.tsv"
	for g in "$OUT"/samples/groups/*.nodes; do
		b=$(basename "$g" .nodes)
		awk -v b="$b" '{print $1"\t"b}' "$g" >>"$OUT/samples/node_group.tsv"
	done
}
gmap_for_chr() {
	local c=$1 d=${TRACE_GENETIC_MAP_DIR:-}
	[[ -n $d ]] || return 1
	[[ -d $d ]] || { echo "ERROR: missing TRACE genetic-map directory: $d" >&2; return 2; }
	python3 - "$d" "$c" <<'PYGMAP'
from pathlib import Path
import re,sys
root,chrom=Path(sys.argv[1]),sys.argv[2]
pat=re.compile(r'(?:^|[^A-Za-z0-9])(?:chr)?'+re.escape(chrom)+r'(?:[^A-Za-z0-9]|$)',re.I)
files=sorted(p for p in root.iterdir() if p.is_file() and pat.search(p.name) and not p.name.endswith(('.tbi','.csi','.md5','.sha256')))
if len(files)!=1:
    raise SystemExit(f'ERROR: expected exactly one genetic map for chr{chrom} in {root}; found {len(files)}: '+', '.join(p.name for p in files))
print(files[0])
PYGMAP
}
prepare_maps() {
	local c source
	if [[ -z ${TRACE_GENETIC_MAP_DIR:-} ]]; then
		echo "TRACE genetic map: constant 1 cM/Mb (set TRACE_GENETIC_MAP_DIR for build-matched maps)"
		return
	fi
	for c in $CHRS; do
		source=$(gmap_for_chr "$c") || return 2
		python3 "$F/trace.py" normalize-map --input "$source" --output "$OUT/maps/chr$c.txt" --chr "$c"
	done
}
run_record() {
	local c post tree cmd source
	printf 'build\t%s\nchromosomes\t%s\nloci_mode\t%s\nt_archaic\t%s\nposterior_threshold\t%s\nphysical_length_threshold\t%s\ngenetic_distance_threshold\t%s\nwindow_size\t%s\ngroup_size\t%s\n' \
		"$BUILD" "$CHRS" "$LOCI_MODE" "$TARC" "$POST" "$MINBP" "$MINCM" "$WIN" "$GROUP_SIZE"
	printf 'backend\t%s\nseed\t%s\nruntime\tgu-trace-input-provenance-v2\n' "$TRACE_BACKEND" "$TRACE_SEED"
	printf 'genetic_map_mode\t%s\n' "${TRACE_GENETIC_MAP_DIR:-constant_1_cM_per_Mb}"
	sha256sum "$F/trace.py" "$F/trace.sh" "$F/0.common.py"
	if [[ $TRACE_BACKEND == cuda ]]; then
		sha256sum "$F/trace_gpu/backend.py" "${TRACE_GPU_HOME:-${GU_SOFT:-$HOME/.local/share/gu}/trace-gpu}/libtrace_gpu.so"
	fi
	[[ -z $LOCI_FILE ]] || sha256sum "$LOCI_FILE"
	sha256sum "$SAMPLE_MAP"
	[[ -z $TARGET_SAMPLES ]] || sha256sum "$TARGET_SAMPLES"
	[[ -z ${TRACE_INCLUDE_REGIONS:-} ]] || sha256sum "$TRACE_INCLUDE_REGIONS"
	for c in $CHRS; do
		[[ ! -s $ARG_DIR/chr$c.sample_map.tsv ]] || sha256sum "$ARG_DIR/chr$c.sample_map.tsv"
		[[ ! -s $OUT/regions/chr$c.bed ]] || sha256sum "$OUT/regions/chr$c.bed"
		if [[ -n ${TRACE_GENETIC_MAP_DIR:-} ]]; then
			source=$(gmap_for_chr "$c") || return 2
			sha256sum "$source" "$OUT/maps/chr$c.txt"
		fi
	done
	while IFS=$'\t' read -r c post tree; do
		[[ $c == chr ]] && continue
		stat -Lc "tree_chr${c}_p${post}\\t%n:%s:%y" "$tree"
	done <"$OUT/manifest/tree_files.tsv"
	for cmd in trace-extract trace-infer trace-summarize; do
		command -v "$cmd" >/dev/null 2>&1 && sha256sum "$(command -v "$cmd")"
	done
	python3 - <<'PYSOURCE'
import hashlib,importlib.util
for name in ('tracehmm','tracehmm.trace','tracehmm.infer_cli','tracehmm.extract_cli','tracehmm.summarize_cli','trace_utils'):
    try: spec=importlib.util.find_spec(name)
    except ModuleNotFoundError: spec=None
    if spec is not None and spec.origin:
        from pathlib import Path
        p=Path(spec.origin)
        print('software_sha256\t'+name+'\t'+hashlib.sha256(p.read_bytes()).hexdigest())
PYSOURCE
}
check_run_provenance() {
	local current=$OUT/run.meta.tsv candidate=$OUT/run.meta.next.tsv
	run_record >"$candidate"
	if [[ ! -s $current && $REPLACE == 0 ]] && find "$OUT/extract" "$OUT/infer" "$OUT/calls" -type f \( -name '*.npz' -o -name '*.summary.txt' \) -print -quit | grep -q .; then
		rm -f "$candidate"
		echo "ERROR: cached TRACE outputs have no scientific provenance; use a fresh output or --replace-trace TRUE" >&2
		return 1
	fi
	if [[ -s $current && $REPLACE == 0 ]] && ! cmp -s "$current" "$candidate"; then
		rm -f "$candidate"
		echo "ERROR: TRACE scientific settings, ARG inputs, sample map/selection, genetic map, mask or runtime changed; use a fresh output directory or --replace-trace TRUE" >&2
		return 1
	fi
	mv "$candidate" "$current"
}
prepare_run() {
	write_manifest
	prepare_regions
	prepare_sample_map
	split_groups
	prepare_maps
	check_run_provenance
}

extract_one() {
	local c=$1 gf=$2 group=$3 post=$4 tree=$5 prefix ids npost
	local -a opts=()
	prefix=$OUT/extract/${group}.chr${c}.p${post}
	ids=$(paste -sd, "$gf")
	npost=$(awk -v c="$c" 'NR>1&&$1==c{n++}END{print n+0}' "$OUT/manifest/tree_files.tsv")
	((npost > 1)) && opts+=(--window-size "$WIN")
	[[ ! -s $OUT/regions/chr$c.bed ]] || opts+=(--include-regions "$OUT/regions/chr$c.bed" --chrom "chr$c")
	python3 "$F/trace.py" output --prefix "$prefix" --suffix .npz -- trace-extract --tree-file "$tree" -t "$TARC" --individuals "$ids" "${opts[@]}" 2>&1 | tee "$TRACE_LOG_DIR/extract.${group}.chr${c}.p${post}.log"
}
run_extract() {
	local status=0 gf group c post tree pid
	local -a pids=()
	for gf in "$OUT"/samples/groups/*.nodes; do
		group=$(basename "$gf" .nodes)
		while IFS=$'\t' read -r c post tree; do
			[[ $c == chr ]] && continue
			extract_one "$c" "$gf" "$group" "$post" "$tree" &
			pids+=("$!")
			trace_children=("${pids[@]}")
			if ((${#pids[@]} >= JEX)); then
				wait "${pids[0]}" || return 1
				pids=("${pids[@]:1}")
				trace_children=("${pids[@]}")
			fi
		done <"$OUT/manifest/tree_files.tsv"
	done
	for pid in "${pids[@]}"; do wait "$pid" || status=1; done
	trace_children=()
	((status == 0)) || exit 1
	for gf in "$OUT"/samples/groups/*.nodes; do
		group=$(basename "$gf" .nodes)
		for c in $CHRS; do find "$OUT/extract" -maxdepth 1 -name "${group}.chr${c}.p*.npz" -print | sort >"$OUT/datafiles/${group}.chr${c}.txt"; done
	done
}

infer_one() {
	local node=$1 group=$2 c npost gmap C D G out any_multi=0
	local -a chroms=() datafiles=() npzs=() gmaps=() opt=()
	out=$OUT/infer/hap$node
	for c in $CHRS; do
		chroms+=("chr$c")
		npost=$(wc -l <"$OUT/datafiles/${group}.chr${c}.txt")
		((npost > 1)) && any_multi=1
		datafiles+=("$OUT/datafiles/${group}.chr${c}.txt")
		npzs+=("$(head -1 "$OUT/datafiles/${group}.chr${c}.txt")")
		if [[ -n ${TRACE_GENETIC_MAP_DIR:-} ]]; then gmaps+=("$OUT/maps/chr$c.txt"); fi
	done
	# A prior run may have produced only the final chromosome before failing.
	# Reuse is valid only when every requested chromosome output exists.
	local -a suffix_args=()
	for c in $CHRS; do suffix_args+=(--suffix ".chr${c}.xss.npz"); done
	C=$(
		IFS=,
		echo "${chroms[*]}"
	)
	if ((any_multi)); then
		D=$(
			IFS=,
			echo "${datafiles[*]}"
		)
		opt=(--data-files "$D")
	else
		D=$(
			IFS=,
			echo "${npzs[*]}"
		)
		opt=(--npz-files "$D")
	fi
	if ((${#gmaps[@]} == ${#chroms[@]})); then
		G=$(
			IFS=,
			echo "${gmaps[*]}"
		)
		opt+=(--genetic-maps "$G")
	fi
	python3 "$F/trace.py" output --prefix "$out" "${suffix_args[@]}" -- "${INFER_CMD[@]}" -i "$node" "${opt[@]}" --chroms "$C" --seed "$TRACE_SEED" 2>&1 | tee "$TRACE_LOG_DIR/infer.hap${node}.log"
}
run_infer() {
	[[ -s $OUT/samples/node_group.tsv ]] || run_extract
	local status=0 node group pid
	local -a pids=()
	while IFS=$'\t' read -r node group; do
		infer_one "$node" "$group" &
		pids+=("$!")
		trace_children=("${pids[@]}")
		if ((${#pids[@]} >= JIN)); then
			wait "${pids[0]}" || return 1
			pids=("${pids[@]:1}")
			trace_children=("${pids[@]}")
		fi
	done <"$OUT/samples/node_group.tsv"
	for pid in "${pids[@]}"; do wait "$pid" || status=1; done
	trace_children=()
	((status == 0)) || exit 1
}
summarize_one() {
	local node=$1 c f out FSTR CSTR
	local -a files=() chroms=()
	out=$OUT/calls/hap$node
	for c in $CHRS; do
		f=$OUT/infer/hap${node}.chr${c}.xss.npz
		[[ -s $f ]] || return 1
		files+=("$f")
		chroms+=("chr$c")
	done
	FSTR=$(
		IFS=,
		echo "${files[*]}"
	)
	CSTR=$(
		IFS=,
		echo "${chroms[*]}"
	)
	python3 "$F/trace.py" output --prefix "$out" --suffix .summary.txt -- trace-summarize -f "$FSTR" -c "$CSTR" --posterior-threshold "$POST" --physical-length-threshold "$MINBP" --genetic-distance-threshold "$MINCM" 2>&1 | tee "$TRACE_LOG_DIR/summarize.hap${node}.log"
}
run_segments() {
	[[ -s $OUT/samples/tree_nodes.txt ]] || run_extract
	local status=0 node pid
	local -a pids=()
	while read -r node; do
		summarize_one "$node" &
		pids+=("$!")
		trace_children=("${pids[@]}")
		if ((${#pids[@]} >= JSUM)); then
			wait "${pids[0]}" || return 1
			pids=("${pids[@]:1}")
			trace_children=("${pids[@]}")
		fi
	done <"$OUT/samples/tree_nodes.txt"
	for pid in "${pids[@]}"; do wait "$pid" || status=1; done
	trace_children=()
	((status == 0)) || exit 1
	local -a args=(--root "$OUT" --sample-map "$SAMPLE_MAP" --build "$BUILD" --filename "$FINAL_NAME")
	[[ -n $LOCI_FILE ]] && args+=(--loci "$LOCI_FILE")
	python3 "$F/trace.py" combine "${args[@]}"
}
prepare_run
case "$ACTION" in
	run)
		run_extract
		run_infer
		run_segments
		;;
	extract) run_extract ;;
	infer)
		run_infer
		run_segments
		;;
	segments) run_segments ;;
	check) echo "TRACE CHECK PASSED: manifest=$OUT/manifest/tree_files.tsv sample_map=$SAMPLE_MAP" ;;
esac

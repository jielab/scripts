#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
GU_F=$HERE
for f in arg.infer_tsinfer.py arg.check_tsinfer.py; do
  [[ -s $GU_F/$f ]] || { echo "ERROR: ARG helper is missing: $GU_F/$f" >&2; exit 2; }
done

GU_ENV=${GU_ENV_PREFIX:-$HOME/anaconda3/envs/gu}
[[ ! -x $GU_ENV/bin/python ]] || export PATH="$GU_ENV/bin:$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ACTION=${REFGEN_ARG_ACTION:-build}
IN_DIR=${REFGEN_VCF_DIR:?}
OUT_DIR=${REFGEN_ARG_DIR:?}
ARG_VCF_DIR=${REFGEN_GEN4ARG_DIR:?}
THREADS=${REFGEN_THREADS:-8}
export JOBS=${REFGEN_JOBS:-1}
SAMPLE_MATCH_MIN_WORK=${REFGEN_SAMPLE_MATCH_MIN_WORK:-50000000}
export MISS=${REFGEN_MISSING_MAX:-0.05}
export MAC=${REFGEN_MAC_MIN:-2}
REPLACE=${REFGEN_REPLACE:-FALSE}
FORMAT=${REFGEN_ARG_FORMAT:-native}
TRACE_DIR=$OUT_DIR/trace/tsinfer
TMP_PARENT=${REFGEN_TSINFER_TMP_ROOT:-/tmp/refgen/arg/tsinfer/$(basename -- "${REFGEN_GEN_ROOT:?}").b${REFGEN_GRCH:-37}}
WORK=
if [[ $ACTION == check ]]; then
  [[ -d $OUT_DIR ]] || { echo "ERROR: ARG directory is missing: $OUT_DIR" >&2; exit 2; }
else
  mkdir -p "$ARG_VCF_DIR" "$TMP_PARENT"
  [[ $ACTION == prepare ]] || mkdir -p "$OUT_DIR"
  TMP_PARENT=$(cd -- "$TMP_PARENT" && pwd -P)
  WORK=$(mktemp -d "$TMP_PARENT/run.XXXXXX")
  cleanup(){ case "$WORK" in "$TMP_PARENT"/run.*) rm -rf -- "$WORK";; *) echo "WARNING: refusing to remove unexpected temporary path: $WORK" >&2;; esac; }
  trap cleanup EXIT
fi

CHRS=$(python3 - "${REFGEN_ARG_CHRS:-1-22,X}" <<'PY'
import re,sys
out=[]
for x in re.split(r"[,;\s]+",sys.argv[1].strip()):
    x=re.sub(r"^chr","",x,flags=re.I)
    m=re.fullmatch(r"(\d+)-(\d+)",x)
    if m: out += [str(i) for i in range(int(m[1]),int(m[2])+1)]
    elif x: out.append(x.upper())
print(" ".join(dict.fromkeys(out)))
PY
)

prepared_vcf(){ printf '%s/chr%s.vcf.gz\n' "$ARG_VCF_DIR" "$1"; }
tmp_for(){ printf '%s/chr%s\n' "$WORK" "$1"; }
source "$HERE/arg.prep_tsinfer.sh"
build_record(){
  local c=$1
  prepared_record "$c"
  python3 "$HERE/arg.cache.py" --text-file "$HERE/arg.tsinfer.sh" --text-file "$GU_F/arg.infer_tsinfer.py" --text-file "$GU_F/arg.check_tsinfer.py" --text-file "$HERE/arg_x.py" --file "$(prepared_vcf "$c")" --file "$(prepared_vcf "$c").input.meta.tsv"
  printf 'sample_match_min_work\t%s\nrecombination_rate\t%s\nmutation_rate\t%s\none_bit\t%s\n' \
    "$SAMPLE_MATCH_MIN_WORK" "${REFGEN_RECOMB_RATE:-1e-8}" "${REFGEN_MUT_RATE:-1.25e-8}" "${REFGEN_ARG_ONE_BIT:-0}"
}
build_one(){
  local c=$1 in out meta expected
  local -a extra=()
  in=$(prepared_vcf "$c"); out=$OUT_DIR/chr$c.trees; meta=$out.input.meta.tsv; expected=$WORK/chr$c.arg.meta.tsv
  build_record "$c" > "$expected"
  [[ ${REFGEN_ARG_ONE_BIT:-0} != 1 ]] || extra+=(--one-bit)
  if [[ $REPLACE == FALSE && -s $out && -s $meta && -s $OUT_DIR/chr$c.sample_map.tsv && -s $OUT_DIR/chr$c.arg_qc.json ]] && cmp -s "$expected" "$meta"; then
    python3 "$HERE/arg.cache.py" --file "$out" --file "$OUT_DIR/chr$c.sample_map.tsv" --file "$OUT_DIR/chr$c.arg_qc.json" > "$expected.outputs"
    if cmp -s "$expected.outputs" "$out.outputs.json"; then check_one "$c"; echo "SKIP verified tsinfer ARG: $out"; return; fi
  fi
  rm -f "$meta" "$out.outputs.json"
  python3 "$GU_F/arg.infer_tsinfer.py" --chr "$c" --vcf "$in" --out "$out" --sample-map "$OUT_DIR/chr$c.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr$c.arg_qc.json" --threads "$THREADS" --tmp-root "$(tmp_for "$c")" \
    --sample-match-min-work "$SAMPLE_MATCH_MIN_WORK" --recombination-rate "${REFGEN_RECOMB_RATE:-1e-8}" \
    --mutation-rate "${REFGEN_MUT_RATE:-1.25e-8}" "${extra[@]}"
  python3 "$GU_F/arg.check_tsinfer.py" --trees "$out" --sample-map "$OUT_DIR/chr$c.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr$c.arg_qc.json" --input-vcf "$in" --chr "$c"
  python3 "$HERE/arg.cache.py" --file "$out" --file "$OUT_DIR/chr$c.sample_map.tsv" --file "$OUT_DIR/chr$c.arg_qc.json" > "$out.outputs.json.next"
  mv -f "$out.outputs.json.next" "$out.outputs.json"
  mv -f "$expected" "$meta"
}
check_one(){
  local c=$1 tree=$OUT_DIR/chr$1.trees map=$OUT_DIR/chr$1.sample_map.tsv qc=$OUT_DIR/chr$1.arg_qc.json
  for f in "$tree" "$map" "$qc"; do [[ -s $f ]] || { echo "ERROR: missing GU TRACE ARG output for chr$c: $f" >&2; return 1; }; done
  if [[ $c == X ]]; then
    python3 "$GU_F/arg.check_tsinfer.py" --trees "$tree" --sample-map "$map" --qc-json "$qc" --input-vcf "$(prepared_vcf "$c")" --chr "$c"
    return
  fi
  python3 - "$tree" "$map" "$qc" "$c" <<'PY'
import csv,json,sys,tskit
tree,sample_map,qc,chrom=sys.argv[1:]
ts=tskit.load(tree)
with open(sample_map,newline="") as f: n=sum(1 for _ in csv.DictReader(f,delimiter="\t"))
with open(qc) as f: json.load(f)
if min(ts.num_trees,ts.num_sites,ts.num_samples,n)<1: raise SystemExit(f"ERROR: empty GU TRACE ARG for chr{chrom}")
if n != ts.num_samples: raise SystemExit(f"ERROR: sample map rows {n} != tree samples {ts.num_samples} for chr{chrom}")
print(f"tsinfer CHECK PASS chr{chrom} trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites}")
PY
}
trace_one(){
  local c=$1
  local -a cmd=(python3 "$HERE/arg.trace.py" --method tsinfer
    --tree "$OUT_DIR/chr$c.trees" --sample-map "$OUT_DIR/chr$c.sample_map.tsv"
    --out-dir "$TRACE_DIR" --chr "$c")
  [[ $ACTION != check ]] || cmd+=(--check)
  [[ $REPLACE == FALSE ]] || cmd+=(--replace)
  "${cmd[@]}"
}

source "$HERE/arg.jobs.sh"
build_chr(){ prep_one "$1"; build_one "$1"; }
read -r -a CHR_ARRAY <<< "$CHRS"

case "$ACTION" in
  prepare) run_parallel prep_one "${CHR_ARRAY[@]}";;
  build)
    run_parallel build_chr "${CHR_ARRAY[@]}"
    cat > "$OUT_DIR/ARG_TSINFER_BUILD.txt" <<META
created=$(date -Is)
producer=/mnt/d/scripts/gu/arg.sh build
method=tsinfer
build=GRCh${REFGEN_GRCH:-37}
vcf_dir=$IN_DIR
chromosomes=$CHRS
META
    ;;
  check) run_parallel check_one "${CHR_ARRAY[@]}";;
esac
if [[ $FORMAT == trace && $ACTION != prepare ]]; then
  run_parallel trace_one "${CHR_ARRAY[@]}"
  if [[ $ACTION != check ]]; then
    printf 'method\ttsinfer\nformat\ttrace\nbuild\tGRCh%s\nchromosomes\t%s\n' \
      "${REFGEN_GRCH:-37}" "$CHRS" > "$TRACE_DIR/ARG_TRACE_BUILD.tsv"
  fi
fi
echo "tsinfer ARG data preparation $ACTION complete: $OUT_DIR"

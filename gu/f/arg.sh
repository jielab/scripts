#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd); F=$ROOT/f
ACTION=${1:-build}; shift || true
BACKEND=${GU_ARG_BACKEND:-tsinfer}
IN_DIR=${GU_TARGET_VCF_DIR:?use ./gu.sh arg build --target-vcf-dir DIR}
GU_TARGET_TMP_DIR=${GU_TARGET_TMP_DIR:?}
OUT_DIR=${GU_ARG_DIR:-${GU_TARGET_ROOT:-.}/arg}
ARG_VCF_DIR=${GU_ARG_VCF_DIR:-${GU_TARGET_ROOT:-.}/vcf.4arg}
THREADS=${GU_ARG_THREADS:-8}
SAMPLE_MATCH_MIN_WORK=${GU_ARG_SAMPLE_MATCH_MIN_WORK:-50000000}
CHRS=${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}
MISS=${GU_ARG_MISSING_MAX:-0.05}; MAC=${GU_ARG_MAC_MIN:-2}
REPLACE=${GU_ARG_REPLACE:-0}
[[ $REPLACE == 0 || $REPLACE == 1 ]] || { echo "ERROR: internal replace-arg state must be 0 or 1" >&2; exit 2; }
mkdir -p "$OUT_DIR" "$ARG_VCF_DIR"
ARG_TMP_PARENT=${GU_ARG_TMP_ROOT:-/tmp/gu/arg/${GU_TARGET_NAMESPACE:-target}/${GU_BUILD:-${GRCH:-unknown}}}
mkdir -p "$ARG_TMP_PARENT"
ARG_TMP_PARENT=$(cd -- "$ARG_TMP_PARENT" && pwd -P)
ARG_WORK_ROOT=$(mktemp -d "$ARG_TMP_PARENT/run.XXXXXX")
cleanup_arg_work() {
  case "$ARG_WORK_ROOT" in
    "$ARG_TMP_PARENT"/run.*) rm -rf -- "$ARG_WORK_ROOT" ;;
    *) echo "WARNING: refusing to remove unexpected ARG temporary path: $ARG_WORK_ROOT" >&2 ;;
  esac
}
trap cleanup_arg_work EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

usage(){ cat <<EOF2
Usage: ./gu.sh arg build
GU_ARG_BACKEND=tsinfer|external

Recommended:
  1KG exploratory/engineering: GU_ARG_BACKEND=tsinfer
  TRACE sensitivity validation: build posterior ARGs with SINGER, convert to .trees/.tsz,
                                then use GU_ARG_BACKEND=external and GU_ARG_DIR=...

The ARG workflow does not automatically split chr1/chrX into coordinate chunks for TRACE.
Coordinate chunks are different genomic regions, not posterior ARG replicates.
The default chromosome set is 1..22 and X. An explicit GU_CHRS value is used unchanged.
EOF2
}
[[ $ACTION == help || $ACTION == -h || $ACTION == --help ]] && { usage; exit 0; }

prep_one(){
  local c=$1 in out tmp meta expected
  in=$IN_DIR/chr${c}.vcf.gz; out=$(prepared_vcf "$c"); tmp=$(arg_tmp "$c")
  meta=$out.input.meta.tsv; expected=$ARG_WORK_ROOT/chr${c}.prepared.meta.tsv
  prepared_record "$c" > "$expected"
  if [[ $REPLACE == 0 && -s $out && -s $out.tbi && -s $meta ]] && cmp -s "$expected" "$meta"; then
    echo "SKIP: prepared ARG VCF already exists: $out"
    return 0
  fi
  [[ -s $in ]] || { echo "ERROR missing $in" >&2; return 1; }
  local replace_word=FALSE; [[ $REPLACE == 0 && ! -e $out ]] || replace_word=TRUE
  bash "$F/arg_vcf_prep.sh" --in "$in" --out "$out" --threads "$THREADS" --missing-max "$MISS" --mac-min "$MAC" --tmp-root "$tmp" --replace "$replace_word"
  mv -f -- "$expected" "$meta"
}

prepared_vcf(){ printf '%s/chr%s.vcf.gz\n' "$ARG_VCF_DIR" "$1"; }
arg_tmp(){ printf '%s/chr%s\n' "$ARG_WORK_ROOT" "$1"; }

prepared_record(){
  local c=$1 in=$IN_DIR/chr${c}.vcf.gz index
  [[ -s $in ]] || { echo "ERROR missing $in" >&2; return 1; }
  if [[ -s $in.tbi ]]; then index=$in.tbi; elif [[ -s $in.csi ]]; then index=$in.csi; else echo "ERROR missing VCF index for $in" >&2; return 1; fi
  stat -c 'input\t%n:%s:%Y' "$in" "$index" "$F/arg_vcf_prep.sh"
  printf 'missing_max\t%s\nmac_min\t%s\nthreads\t%s\n' "$MISS" "$MAC" "$THREADS"
}

arg_build_record(){
  local c=$1 source_meta=$GU_TARGET_TMP_DIR/chr${c}/source.tsv line
  [[ -s $source_meta ]] || { echo "ERROR: target source provenance is missing: $source_meta" >&2; return 1; }
  while IFS= read -r line; do printf 'target_source\tchr%s\t%s\n' "$c" "$line"; done < "$source_meta"
  stat -c 'script\t%n:%s:%Y' "$F/arg_vcf_prep.sh" "$F/arg_tsinfer.py"
  stat -c 'prepared_vcf\t%n:%s:%Y' "$(prepared_vcf "$c")" "$(prepared_vcf "$c").input.meta.tsv"
  printf 'prep_missing_max\t%s\nprep_mac_min\t%s\nthreads\t%s\nsample_match_min_work\t%s\nrecombination_rate\t%s\nmutation_rate\t%s\none_bit\t%s\n' \
    "$MISS" "$MAC" "$THREADS" "$SAMPLE_MATCH_MIN_WORK" \
    "${GU_RECOMB_RATE:-1e-8}" "${GU_MUT_RATE:-1.25e-8}" "${GU_ARG_ONE_BIT:-0}"
}

build_one(){
  local c=$1 in out meta sample_map qc expected
  in=$(prepared_vcf "$c"); out=$OUT_DIR/chr${c}.trees
  meta=$out.input.meta.tsv; sample_map=$OUT_DIR/chr${c}.sample_map.tsv; qc=$OUT_DIR/chr${c}.arg_qc.json
  expected=$ARG_WORK_ROOT/chr${c}.arg.meta.tsv
  arg_build_record "$c" > "$expected"
  local extra=()
  [[ ${GU_ARG_ONE_BIT:-0} == 1 ]] && extra+=(--one-bit)
  if [[ $REPLACE != 1 && -s $out && -s $meta && -s $sample_map && -s $qc ]] && cmp -s "$expected" "$meta"; then echo "SKIP: ARG outputs exist: $out"; return 0; fi
  python3 "$F/arg_tsinfer.py" --vcf "$in" --out "$out" --sample-map "$OUT_DIR/chr${c}.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr${c}.arg_qc.json" --threads "$THREADS" --tmp-root "$(arg_tmp "$c")" \
    --sample-match-min-work "$SAMPLE_MATCH_MIN_WORK" \
    --recombination-rate "${GU_RECOMB_RATE:-1e-8}" --mutation-rate "${GU_MUT_RATE:-1.25e-8}" \
    "${extra[@]}"
  mv -f -- "$expected" "$meta"
}

check_one(){
  local c=$1 in out
  in=$(prepared_vcf "$c"); out=$OUT_DIR/chr${c}.trees
  local meta=$out.input.meta.tsv sample_map=$OUT_DIR/chr${c}.sample_map.tsv qc=$OUT_DIR/chr${c}.arg_qc.json
  for required in "$in" "$out" "$meta" "$sample_map" "$qc"; do
    [[ -s $required ]] || { echo "ERROR: ARG post-build check missing or empty: $required" >&2; return 1; }
  done
  python3 "$F/arg_check.py" --trees "$out" --sample-map "$sample_map" --qc-json "$qc" --input-vcf "$in" --chr "$c"
}

check_external(){
  local c n bad=0
  for c in $CHRS; do
    n=$(find "$OUT_DIR" -type f \( -name '*.trees' -o -name '*.tsz' \) -print | awk -v c="$c" 'BEGIN{IGNORECASE=1}{b=$0;gsub(/.*\//,"",b);if(b ~ ("(^|[^A-Za-z0-9])chr" c "([^A-Za-z0-9]|$)")) n++}END{print n+0}')
    ((n>0)) || { echo "ERROR: external ARG missing chr$c under $OUT_DIR" >&2; bad=1; }
  done
  ((bad==0))
}

case "$BACKEND:$ACTION" in
  tsinfer:build)
    for c in $CHRS; do prep_one "$c"; build_one "$c"; check_one "$c"; done
    echo "ARG BUILD COMPLETE: chromosomes=$CHRS post_build_check=PASS output=$OUT_DIR"
    ;;
  external:build) check_external ;;
  *) echo "ERROR unsupported GU_ARG_BACKEND/ACTION: $BACKEND/$ACTION" >&2; usage; exit 2 ;;
esac

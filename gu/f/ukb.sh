#!/usr/bin/env bash
# UK Biobank helpers for GU.
# Primary supported UKB input in this installation: Field 22438 WTCHG phased haplotype BGENs (GRCh37, chr1-22).
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
F=$ROOT/f
ACTION=${1:-help}; shift || true

HAP=${UKB_HAP_ROOT:-/mnt/h/ukbGen/hap}
TYPED=${UKB_TYPED_ROOT:-/mnt/h/ukbGen/typ}
OUT=${UKB_WORK:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu}/ukb}
CHRS=${GU_CHRS:-22}
THREADS=${UKB_THREADS:-8}
REF_FASTA=${UKB_REF_FASTA:-}
mkdir -p "$OUT"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }; }
die(){ echo "ERROR: $*" >&2; exit 1; }
bgen_for(){ find "$HAP" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2.bgen" -print -quit; }
sample_for(){ find "$HAP" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2_s*.sample" -print -quit; }

sample_hash(){
  # Hash IDs in BGEN order. Oxford .sample files have two header/type lines.
  awk 'NR>2 && NF{print $2}' "$1" | sha256sum | awk '{print $1}'
}

variant_count(){
  local b=$1
  if command -v sqlite3 >/dev/null 2>&1 && [[ -s $b.bgi ]]; then
    local n
    n=$(sqlite3 "$b.bgi" 'SELECT count(*) FROM Variant;' 2>/dev/null || true)
    [[ $n =~ ^[0-9]+$ ]] && { echo "$n"; return; }
  fi
  # Fallback is slower but robust.
  bgenix -g "$b" -list 2>/dev/null | awk 'BEGIN{n=0} /^[#]/ {next} NF{n++} END{print n}'
}

phase_check_vcf(){
  local vcf=$1
  python3 - "$vcf" <<'PY'
import sys
import pysam
path=sys.argv[1]
v=pysam.VariantFile(path)
samples=list(v.header.samples)[:20]
checked=het=bad=0
examples=[]
for rec in v:
    checked += 1
    for sid in samples:
        call=rec.samples[sid]
        gt=call.get("GT")
        if gt is None or len(gt)!=2 or any(x is None for x in gt):
            continue
        if gt[0] != gt[1]:
            het += 1
            if not call.phased:
                bad += 1
                if len(examples)<5:
                    examples.append(f"{rec.chrom}:{rec.pos}:{sid}:{gt}")
    if checked >= 2000:
        break
if bad:
    for x in examples:
        print("UNPHASED", x, file=sys.stderr)
    raise SystemExit(f"ERROR: {bad} sampled heterozygous GTs were unphased in {path}")
if het == 0:
    print(f"WARNING: phase QC saw no heterozygotes in {checked} records x {len(samples)} samples", file=sys.stderr)
print(f"phase_check records={checked} samples={len(samples)} heterozygotes={het} unphased={bad}")
PY
}

usage(){ cat <<'HELP'
Usage: ./gu.sh ukb ACTION

Actions:
  inspect-hap    Inspect Field 22438 phased BGENs, index them if necessary, and verify sample order
  make-panel     Create GU sample panel from an Oxford .sample file (pop=UKB, super_pop=ALL)
  batches        Make deterministic ancestry-stratified target batches and optional fixed UKB anchors
  hap-vcf        Export one selected phased-BGEN sample set to phased, reference-oriented VCF
  hap-arg-vcf    Run hap-vcf and then transfer INFO/AA from build-matched 1KG by exact allele match
  inspect-typed  Summarize the older /typ BED/BIM/FAM array dataset

Important variables:
  UKB_HAP_ROOT       default /mnt/h/ukbGen/hap
  UKB_REF_FASTA      GRCh37 FASTA used to set REF during BGEN -> VCF conversion (recommended/required)
  UKB_KEEP           one-column UKB sample list to export; usually a *.joint.txt batch file
  UKB_1KG_VCF_DIR    GRCh37 1KG VCF directory containing INFO/AA
  UKB_VCF_OUT        raw phased VCF output directory
  UKB_ARG_VCF_OUT    phased VCF + transferred INFO/AA output directory
  GU_CHRS            chromosomes; start with 22

Minimal chr22 preparation:
  export GU_BUILD=37
  export GU_CHRS=22
  export UKB_HAP_ROOT=/mnt/h/ukbGen/hap
  export UKB_REF_FASTA=/path/to/human_g1k_v37.fasta
  ./gu.sh ukb inspect-hap
  export GU_SAMPLE_PANEL=$GU_ANALYSIS_ROOT/ukb/ukb.sample_panel.tsv
  ./gu.sh ukb make-panel
  GU_ANCHORS_PER_GROUP=1000 GU_BATCH_SIZE=1000 ./gu.sh ukb batches
  export UKB_KEEP=$GU_ANALYSIS_ROOT/ukb/batches/ALL.b0001.joint.txt
  export UKB_1KG_VCF_DIR=/mnt/d/data.BIG/refGen/1kg/GRCH37/vcf
  export UKB_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/ALL.b0001/vcf
  export UKB_ARG_VCF_OUT=$GU_ANALYSIS_ROOT/ukb/work/ALL.b0001/vcf.aa
  ./gu.sh ukb hap-arg-vcf
HELP
}

case "$ACTION" in
  inspect-hap)
    need bgenix; need sha256sum
    report="$OUT/hap.inspect.tsv"
    printf 'chr\tn_samples\tn_variants\tsample_order_sha256\tbgen\tsample\n' > "$report"
    first_hash=""; first_chr=""
    for c in $CHRS; do
      b=$(bgen_for "$c"); s=$(sample_for "$c")
      [[ -s $b ]] || die "missing chr$c BGEN under $HAP"
      [[ -s $s ]] || die "missing chr$c Oxford .sample under $HAP"
      if [[ ! -s $b.bgi ]]; then
        echo "Indexing chr$c BGEN: $b" >&2
        bgenix -g "$b" -index >/dev/null
      fi
      ns=$(( $(wc -l < "$s") - 2 ))
      (( ns > 0 )) || die "bad sample count in $s"
      h=$(sample_hash "$s")
      nv=$(variant_count "$b")
      [[ $nv =~ ^[0-9]+$ && $nv -gt 0 ]] || die "could not count variants in $b"
      if [[ -z $first_hash ]]; then first_hash=$h; first_chr=$c
      elif [[ $h != "$first_hash" ]]; then
        die "sample IDs/order differ between chr$first_chr and chr$c; do not use one panel across chromosomes until resolved"
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$c" "$ns" "$nv" "$h" "$b" "$s" | tee -a "$report"
    done
    echo "PASS: selected chromosomes have identical sample ID order. Report: $report"
    ;;

  make-panel)
    c=${UKB_PANEL_CHR:-22}; s=$(sample_for "$c")
    [[ -s $s ]] || die "missing sample file for chr$c under $HAP"
    o=${GU_SAMPLE_PANEL:-$OUT/ukb.sample_panel.tsv}
    mkdir -p "$(dirname "$o")"
    python3 - "$s" "$o" <<'PY'
import sys
from pathlib import Path
src=Path(sys.argv[1]); out=Path(sys.argv[2])
lines=[x.rstrip('\n') for x in src.open() if x.strip()]
if len(lines)<3: raise SystemExit(f"Bad Oxford sample file: {src}")
header=lines[0].split(); rows=[x.split() for x in lines[2:]]
idx={x.lower():i for i,x in enumerate(header)}
id_col=idx.get('id_2', idx.get('id', idx.get('id_1', 1)))
sex_col=idx.get('sex')
seen=set()
with out.open('w') as h:
    h.write('sample\tpop\tsuper_pop\tsex\n')
    for r in rows:
        if id_col >= len(r): continue
        sid=r[id_col]
        if sid in {'0','NA',''}: continue
        if sid in seen: raise SystemExit(f"Duplicate sample ID in {src}: {sid}")
        seen.add(sid)
        sex='unknown'
        if sex_col is not None and sex_col < len(r):
            sex={'1':'male','2':'female'}.get(r[sex_col], 'unknown')
        h.write(f'{sid}\tUKB\tALL\t{sex}\n')
print(f'Wrote {len(seen)} samples to {out}')
PY
    ;;

  batches)
    : "${GU_SAMPLE_PANEL:?set GU_SAMPLE_PANEL; run ./gu.sh ukb make-panel or supply an ancestry-labelled panel}"
    [[ -s $GU_SAMPLE_PANEL ]] || die "missing GU_SAMPLE_PANEL=$GU_SAMPLE_PANEL"
    args=(--panel "$GU_SAMPLE_PANEL" --outdir "$OUT/batches" --batch-size "${GU_BATCH_SIZE:-1000}" --sample-col "${GU_SAMPLE_COL:-sample}" --group-col "${GU_GROUP_COL:-super_pop}" --anchors-per-group "${GU_ANCHORS_PER_GROUP:-0}")
    [[ -n ${GU_ANCHOR_LIST:-} ]] && args+=(--anchor-list "$GU_ANCHOR_LIST")
    python3 "$F/make_batches.py" "${args[@]}"
    ;;

  hap-vcf)
    need plink2; need tabix; need bcftools; need python3
    : "${UKB_KEEP:?set UKB_KEEP to a one-column sample list; exporting all ~487k samples to VCF is intentionally blocked}"
    [[ -s $UKB_KEEP ]] || die "missing UKB_KEEP=$UKB_KEEP"
    if [[ -z $REF_FASTA || ! -s $REF_FASTA ]]; then
      if [[ ${UKB_ALLOW_PROVISIONAL_REF:-0} != 1 ]]; then
        die "set UKB_REF_FASTA to the build-matched GRCh37 FASTA. Set UKB_ALLOW_PROVISIONAL_REF=1 only for engineering tests."
      fi
      echo "WARNING: no UKB_REF_FASTA; REF/ALT remain provisional and AA transfer may retain few sites." >&2
    fi
    ODIR=${UKB_VCF_OUT:-$OUT/hap_vcf}; mkdir -p "$ODIR"
    keep_tmp=$(mktemp "${TMPDIR:-/tmp}/gu-ukb-keep.XXXXXX")
    trap 'rm -f "$keep_tmp"' EXIT
    { echo '#IID'; awk 'NF && $1!~/^#/{print $1}' "$UKB_KEEP"; } > "$keep_tmp"
    expected=$(awk 'NR>1{n++}END{print n+0}' "$keep_tmp")
    ((expected>0)) || die "UKB_KEEP contains no sample IDs"
    for c in $CHRS; do
      b=$(bgen_for "$c"); s=$(sample_for "$c")
      [[ -s $b && -s $s ]] || die "missing chr$c BGEN/sample"
      p="$ODIR/chr$c"
      cmd=(plink2 --bgen "$b" ref-unknown --sample "$s" --oxford-single-chr "$c" --keep "$keep_tmp" --snps-only just-acgt --max-alleles 2)
      if [[ -n $REF_FASTA && -s $REF_FASTA ]]; then cmd+=(--fa "$REF_FASTA" --ref-from-fa); fi
      cmd+=(--export vcf bgz id-paste=iid --out "$p")
      printf 'Running:' >&2; printf ' %q' "${cmd[@]}" >&2; echo >&2
      "${cmd[@]}"
      tabix -f -p vcf "$p.vcf.gz"
      observed=$(bcftools query -l "$p.vcf.gz" | wc -l)
      [[ $observed -eq $expected ]] || die "chr$c VCF sample count $observed != requested $expected; inspect PLINK keep/sample IDs"
      phase_check_vcf "$p.vcf.gz"
      nv=$(bcftools index -n "$p.vcf.gz")
      printf 'metric\tvalue\nchromosome\t%s\nsamples\t%s\nvariants\t%s\nreference_fasta\t%s\n' "$c" "$observed" "$nv" "${REF_FASTA:-NONE}" > "$p.qc.tsv"
      echo "chr$c phased VCF: $p.vcf.gz samples=$observed variants=$nv"
    done
    ;;

  hap-arg-vcf)
    : "${UKB_KEEP:?set UKB_KEEP to the joint target+anchor sample list}"
    : "${UKB_1KG_VCF_DIR:?set UKB_1KG_VCF_DIR to build-matched 1KG VCF directory containing INFO/AA}"
    RAW=${UKB_VCF_OUT:-$OUT/hap_vcf}
    AAOUT=${UKB_ARG_VCF_OUT:-$OUT/hap_vcf.aa}
    mkdir -p "$AAOUT"
    UKB_VCF_OUT="$RAW" bash "$0" hap-vcf
    for c in $CHRS; do
      ref="$UKB_1KG_VCF_DIR/chr$c.vcf.gz"
      [[ -s $ref ]] || die "missing 1KG AA source: $ref"
      bash "$F/add_aa_from_1kg.sh" --in "$RAW/chr$c.vcf.gz" --ref "$ref" --out "$AAOUT/chr$c.vcf.gz" --threads "$THREADS"
    done
    echo "AA-tagged phased VCFs are ready under: $AAOUT"
    echo "Next set GU_TARGET_VCF_DIR=$AAOUT and run: ./gu.sh arg build"
    ;;

  inspect-typed)
    printf 'chr\tn_variants\tn_samples\n'
    for bim in "$TYPED"/chr*.bim; do
      [[ -e $bim ]] || continue
      c=$(basename "$bim" .bim | sed 's/^chr//'); fam="$TYPED/chr$c.fam"
      printf '%s\t%s\t%s\n' "$c" "$(wc -l < "$bim")" "$([[ -s $fam ]] && wc -l < "$fam" || echo NA)"
    done
    ;;

  help|-h|--help) usage ;;
  *) echo "ERROR: unknown UKB action: $ACTION" >&2; usage >&2; exit 2 ;;
esac

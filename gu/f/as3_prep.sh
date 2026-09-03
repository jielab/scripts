#!/usr/bin/env bash
# Validate direct ArchaicSeeker3 Ref1028 inputs and write a lightweight manifest.
set -euo pipefail

need() { command -v "$1" >/dev/null 2>&1 || { echo "FATAL ERROR: missing command: $1" >&2; exit 1; }; }
for cmd in bcftools awk sort comm cmp mktemp readlink stat; do need "$cmd"; done
AS3_BCFTOOLS_BIN=$(type -P bcftools)
bcftools() {
    "$AS3_BCFTOOLS_BIN" "$@" 2> >(sed '/^\[W::bcf_hdr_check_sanity\] AC should be declared as Number=A$/d' >&2)
}

dir0=${dir0:-/mnt/d}
dir_ref=${dir_ref:-$dir0/data.BIG/refGen}
reference_dir=${AS3_REFERENCE_PANEL_DIR:-$dir_ref/archaic/38/vcf}
reference_map=${AS3_REFERENCE_MAP:-}
out=${AS3_DATA_OUT:?AS3_DATA_OUT is required; run this internal helper through gu.sh}
target_dir=${AS3_TARGET_VCF_DIR:-${AS3_MODERN_HUMAN_VCF_DIR:-$dir_ref/1kg/${GRCH:-38}/vcf}}
target_template=${AS3_TARGET_VCF_TEMPLATE:-}
reference_template=${AS3_REFERENCE_VCF_TEMPLATE:-}
[[ -n $target_template ]] || target_template='chr{CHR}.vcf.gz'
[[ -n $reference_template ]] || reference_template='Ref_Panel.chr{CHR}.vcf.gz'
panel_id=${AS3_REFERENCE_PANEL_ID:-Ref1028}
genome_build=${AS3_GENOME_BUILD:-b${GRCH:-38}}
requested=${AS3_CHRS:-"$(seq 1 22)"}
check_only=${AS3_PREP_CHECK_ONLY:-0}

while (( $# )); do
    case "$1" in
        --reference-panel-dir) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --reference-panel-dir requires DIR" >&2; exit 2; }; reference_dir=$2; shift 2 ;;
        --reference-map) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --reference-map requires FILE" >&2; exit 2; }; reference_map=$2; shift 2 ;;
        --out) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --out requires DIR" >&2; exit 2; }; out=$2; shift 2 ;;
        --chr) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --chr requires LIST" >&2; exit 2; }; requested=${2//,/ }; shift 2 ;;
        --target-vcf-dir) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --target-vcf-dir requires DIR" >&2; exit 2; }; target_dir=$2; shift 2 ;;
        --target-vcf-template) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --target-vcf-template requires TEMPLATE" >&2; exit 2; }; target_template=$2; shift 2 ;;
        --reference-vcf-template) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --reference-vcf-template requires TEMPLATE" >&2; exit 2; }; reference_template=$2; shift 2 ;;
        --genome-build) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --genome-build requires b38" >&2; exit 2; }; genome_build=$2; shift 2 ;;
        --check-only) check_only=1; shift ;;
        --threads|--max-procs|--tmp-root) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a value" >&2; exit 2; }; shift 2 ;;
        --allow-chrx) echo "FATAL ERROR: official Ref1028 supports GRCh38 chr1-22 only" >&2; exit 2 ;;
        *) echo "FATAL ERROR: unknown AS3 input option: $1" >&2; exit 2 ;;
    esac
done

case "${genome_build,,}" in
    b38|38|grch38|hg38) genome_build=b38 ;;
    *) echo "FATAL ERROR: official Ref1028 is GRCh38 and supports chr1-22 only; got $genome_build" >&2; exit 2 ;;
esac
[[ $check_only == 0 || $check_only == 1 ]] || { echo "FATAL ERROR: AS3_PREP_CHECK_ONLY must be 0 or 1" >&2; exit 2; }
[[ -d $reference_dir ]] || { echo "FATAL ERROR: Ref1028 directory missing: $reference_dir" >&2; exit 1; }
[[ -d $target_dir ]] || { echo "FATAL ERROR: target VCF directory missing: $target_dir" >&2; exit 1; }
reference_map=${reference_map:-$reference_dir/Ref_Panel.map.txt}
[[ -s $reference_map ]] || { echo "FATAL ERROR: Ref1028 map missing or empty: $reference_map" >&2; exit 1; }

render() {
    local dir=$1 template=$2 chr=$3 path
    path=$(awk -v value="$template" -v chr="$chr" 'BEGIN{
      token="{CHR}";p=index(value,token)
      if(!p){token="{chr}";p=index(value,token)}
      if(p)value=substr(value,1,p-1) chr substr(value,p+length(token))
      print value
    }')
    [[ $path == /* ]] && printf '%s\n' "$path" || printf '%s/%s\n' "$dir" "$path"
}

valid_vcf() {
    local file=$1 count
    [[ -s $file && ( -s $file.tbi || -s $file.csi ) ]] || return 1
    bcftools view -h "$file" >/dev/null 2>&1 || return 1
    count=$(bcftools index -n "$file" 2>/dev/null || true)
    [[ $count =~ ^[1-9][0-9]*$ ]]
}

validate_single_record_contig() {
    local file=$1 chr=$2 role=$3
    bcftools index -s "$file" | awk -v expected="$chr" -v role="$role" -v file="$file" '
      $3+0>0 {
        c=$1; sub(/^chr/,"",c); n++
        if(c!=expected){print "FATAL ERROR: " role " VCF has records on unexpected contig " $1 ": " file > "/dev/stderr"; bad=1}
      }
      END{
        if(n!=1){print "FATAL ERROR: " role " VCF must contain records on exactly one chromosome (found " n "): " file > "/dev/stderr"; bad=1}
        exit bad
      }'
}

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/gu-as3-ref1028.XXXXXX")
trap 'rm -rf -- "$tmp_root"' EXIT
awk 'BEGIN{FS="[[:space:]]+"}
  NF {
    if(NF!=2){print "FATAL ERROR: Ref1028 map must have exactly two columns at line " NR > "/dev/stderr"; bad=1}
    if($2!="AFR"&&$2!="DEN"&&$2!="NEAN"){print "FATAL ERROR: invalid Ref1028 ancestry label " $2 " at line " NR > "/dev/stderr"; bad=1}
    print $1
    n[$2]++
  }
  END{
    if(n["AFR"]!=146||n["DEN"]!=1||n["NEAN"]!=3){print "FATAL ERROR: Ref1028 map requires exactly AFR=146, DEN=1, and NEAN=3" > "/dev/stderr"; bad=1}
    exit bad
  }' "$reference_map" > "$tmp_root/map.samples.unsorted"
sort "$tmp_root/map.samples.unsorted" > "$tmp_root/map.samples"
duplicate=$(uniq -d "$tmp_root/map.samples" | sed -n '1p')
[[ -z $duplicate ]] || { echo "FATAL ERROR: duplicate Ref1028 map sample: $duplicate" >&2; exit 1; }
read -r n_afr n_den n_nean < <(awk 'BEGIN{FS="[[:space:]]+"}NF{n[$2]++}END{print n["AFR"]+0,n["DEN"]+0,n["NEAN"]+0}' "$reference_map")

mapfile -t chromosomes < <(tr ', ' '\n\n' <<< "$requested" | awk 'NF{c=$1;sub(/^chr/,"",c);sub(/^CHR/,"",c);if(c==23)c="X";if(!seen[c]++)print c}')
(( ${#chromosomes[@]} > 0 )) || { echo "FATAL ERROR: AS3 chromosome list is empty" >&2; exit 2; }
for chr in "${chromosomes[@]}"; do
    [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo "FATAL ERROR: official Ref1028 supports GRCh38 chr1-22 only; got chr$chr" >&2; exit 2; }
done

mkdir -p "$out"
manifest_tmp=$out/.manifest.$$.tsv
release_tmp=$out/.release.$$.tsv
printf 'panel_id\tunit\ttarget_vcf\treference_vcf\treference_map\treference_target_overlap\n' > "$manifest_tmp"
printf 'key\tvalue\npanel_id\t%s\ngenome_build\tGRCh38\nreference_dir\t%s\nreference_map\t%s\nmap_AFR\t%s\nmap_DEN\t%s\nmap_NEAN\t%s\n' \
    "$panel_id" "$reference_dir" "$reference_map" "$n_afr" "$n_den" "$n_nean" > "$release_tmp"

for chr in "${chromosomes[@]}"; do
    target=$(render "$target_dir" "$target_template" "$chr")
    reference=$(render "$reference_dir" "$reference_template" "$chr")
    valid_vcf "$target" || { echo "FATAL ERROR: target chr$chr VCF/index is missing, empty, or unreadable: $target" >&2; exit 1; }
    valid_vcf "$reference" || { echo "FATAL ERROR: Ref1028 chr$chr VCF/index is missing, empty, or unreadable: $reference" >&2; exit 1; }
    # A run directory may expose a native VCF through a short-lived symlink.
    # Persist only the canonical source in the manifest so the AS3 worker never
    # depends on that directory after input preparation.
    target=$(readlink -f -- "$target")
    reference=$(readlink -f -- "$reference")
    validate_single_record_contig "$target" "$chr" target
    validate_single_record_contig "$reference" "$chr" Ref1028

    bcftools query -l "$reference" | sort > "$tmp_root/reference.chr${chr}.samples"
    cmp -s "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" || {
        echo "FATAL ERROR: Ref1028 chr$chr VCF samples do not exactly match $reference_map" >&2
        echo "Only in VCF:" >&2; comm -23 "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" | sed -n '1,10p' >&2
        echo "Only in map:" >&2; comm -13 "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" | sed -n '1,10p' >&2
        exit 1
    }
    target_samples_tmp=$tmp_root/target.chr${chr}.samples
    bcftools query -l "$target" | sort > "$target_samples_tmp"
    [[ -s $target_samples_tmp ]] || { echo "FATAL ERROR: target chr$chr VCF has no samples: $target" >&2; exit 1; }
    target_duplicate=$(uniq -d "$target_samples_tmp" | sed -n '1p')
    [[ -z $target_duplicate ]] || { echo "FATAL ERROR: duplicate target chr$chr sample: $target_duplicate" >&2; exit 1; }
    overlap=$(comm -12 "$target_samples_tmp" "$tmp_root/map.samples" | awk 'END{print NR+0}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$panel_id" "$chr" "$target" "$reference" "$reference_map" "$overlap" >> "$manifest_tmp"
    target_n=$(awk 'END{print NR+0}' "$target_samples_tmp")
    ref_n=$(awk 'END{print NR+0}' "$tmp_root/map.samples")
    target_variants=$(bcftools index -n "$target")
    ref_variants=$(bcftools index -n "$reference")
    printf '[AS3 INPUT chr%s] target_samples=%s reference_samples=%s target_variants=%s reference_variants=%s target_reference_sample_overlap=%s\n' \
        "$chr" "$target_n" "$ref_n" "$target_variants" "$ref_variants" "$overlap"
done

if [[ ! -s $out/manifest.tsv ]] || ! cmp -s "$manifest_tmp" "$out/manifest.tsv"; then mv -f "$manifest_tmp" "$out/manifest.tsv"; else rm -f "$manifest_tmp"; fi
if [[ ! -s $out/release.tsv ]] || ! cmp -s "$release_tmp" "$out/release.tsv"; then mv -f "$release_tmp" "$out/release.tsv"; else rm -f "$release_tmp"; fi
printf 'AS3 REF1028 INPUT CHECK PASSED: build=GRCh38 panel=%s chromosomes=%s map(AFR=%s,DEN=%s,NEAN=%s) manifest=%s mode=%s\n' \
    "$panel_id" "${chromosomes[*]}" "$n_afr" "$n_den" "$n_nean" "$out/manifest.tsv" "$([[ $check_only == 1 ]] && printf check || printf prepare)"

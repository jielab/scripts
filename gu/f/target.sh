#!/usr/bin/env bash
# Target genotype discovery/validation shared by GU entry points and tests.

# Some published 1KG VCFs declare INFO/AC as Number=1. GU analyses use GT and
# do not consume INFO/AC, so hide only this htslib header warning while keeping
# every other bcftools diagnostic visible.
if GU_BCFTOOLS_BIN=$(type -P bcftools 2>/dev/null); then
  bcftools() {
    "$GU_BCFTOOLS_BIN" "$@" 2> >(sed '/^\[W::bcf_hdr_check_sanity\] AC should be declared as Number=A$/d' >&2)
  }
fi

gu_target_genotype_base() {
  local prefix=${1:?gu_target_genotype_base requires PREFIX}
  local chr=${2:?gu_target_genotype_base requires CHR}
  local x_mode=${3:-${GU_CHRX_TARGET_MODE:-male}}
  chr=${chr#chr}; chr=${chr#CHR}; chr=${chr^^}; [[ $chr == 23 ]] && chr=X
  if [[ $chr == X ]]; then
    case "$x_mode" in
      par|xy) printf '%sXY\n' "$prefix" ;;
      male|pure|nonpar|'') printf '%sX.male\n' "$prefix" ;;
      *) echo "ERROR: invalid chrX target mode: $x_mode" >&2; return 2 ;;
    esac
  else
    printf '%s%s\n' "$prefix" "$chr"
  fi
}

gu_chrx_par_bounds_bed() {
  local build=${1:?gu_chrx_par_bounds_bed requires BUILD}
  case "${build,,}" in
    37|b37|grch37|hg19) printf '60000\t2699520\t154931043\t155260560\n' ;;
    38|b38|grch38|hg38) printf '10000\t2781479\t155701382\t156030895\n' ;;
    *) echo "ERROR: chrX/PAR routing needs GRCh37 or GRCh38, got: $build" >&2; return 2 ;;
  esac
}

gu_classify_chrx_loci() {
  local bed=${1:?gu_classify_chrx_loci requires BED}
  local build=${2:?gu_classify_chrx_loci requires BUILD}
  local p1s p1e p2s p2e
  read -r p1s p1e p2s p2e < <(gu_chrx_par_bounds_bed "$build") || return 2
  awk -F '\t' -v p1s="$p1s" -v p1e="$p1e" -v p2s="$p2s" -v p2e="$p2e" '
    BEGIN{OFS="\t"}
    {
      c=toupper($1);sub(/^CHR/,"",c);if(c=="23")c="X";if(c!="X")next
      n++;s=$2+0;e=$3+0
      if(e<=s){print "ERROR: invalid chrX BED interval: "$0 > "/dev/stderr";bad=1;next}
      in1=(s>=p1s&&e<=p1e);in2=(s>=p2s&&e<=p2e)
      overlaps=(s<p1e&&e>p1s)||(s<p2e&&e>p2s)
      if(in1||in2)par++
      else if(overlaps){print "ERROR: chrX locus crosses a PAR/non-PAR boundary: "$0 > "/dev/stderr";boundary=1}
      else pure++
    }
    END{
      if(bad||boundary)exit 2
      if(!n){print "none";exit}
      if(par&&pure){print "ERROR: one --loci request cannot mix pure-X and PAR intervals; run them separately" > "/dev/stderr";exit 2}
      print par?"par":"male"
    }' "$bed"
}

gu_make_chrxy_extract_bed() {
  local input=${1:?gu_make_chrxy_extract_bed requires INPUT}
  local output=${2:?gu_make_chrxy_extract_bed requires OUTPUT}
  local build=${3:?gu_make_chrxy_extract_bed requires BUILD}
  local p1s p1e p2s p2e
  read -r p1s p1e p2s p2e < <(gu_chrx_par_bounds_bed "$build") || return 2
  awk -F '\t' -v p1s="$p1s" -v p1e="$p1e" -v p2s="$p2s" -v p2e="$p2e" '
    BEGIN{OFS="\t"}
    {
      c=toupper($1);sub(/^CHR/,"",c);if(c=="23")c="X"
      if(c!="X"){print;next}
      if($2>=p1s&&$3<=p1e)$1="PAR1"
      else if($2>=p2s&&$3<=p2e)$1="PAR2"
      else {print "ERROR: non-PAR interval reached chrXY extraction: "$0 > "/dev/stderr";bad=1;next}
      print
    }
    END{exit bad}' "$input" > "$output"
}

gu_pvar_cat() {
  local path=${1:?gu_pvar_cat requires PVAR}
  case "$path" in
    *.zst)
      if command -v zstdcat >/dev/null 2>&1; then zstdcat -- "$path"
      elif command -v zstd >/dev/null 2>&1; then zstd -dc -- "$path"
      else echo "ERROR: zstdcat or zstd is required to read $path" >&2; return 2
      fi
      ;;
    *.gz) gzip -dc -- "$path" ;;
    *) cat -- "$path" ;;
  esac
}

gu_link_if_needed() {
  local source=${1:?gu_link_if_needed requires SOURCE}
  local dest=${2:?gu_link_if_needed requires DEST}
  local source_resolved dest_resolved
  source_resolved=$(readlink -f -- "$source") || return 2
  if [[ -L $dest ]]; then
    dest_resolved=$(readlink -f -- "$dest" 2>/dev/null || true)
    [[ $dest_resolved != "$source_resolved" ]] || return 0
  fi
  rm -f -- "$dest"
  ln -s -- "$source" "$dest"
}

gu_validate_chrx_male_pfile() {
  local base=${1:?gu_validate_chrx_male_pfile requires PFILE_PREFIX}
  local build=${2:?gu_validate_chrx_male_pfile requires BUILD}
  local pvar="" p1s p1e p2s p2e variants bad_chr par_variants samples nonmale
  local hethap=NA stats=${base}.basic-stats.txt
  case "${build,,}" in
    37|b37|grch37|hg19) build=37; p1s=60001; p1e=2699520; p2s=154931044; p2e=155260560 ;;
    38|b38|grch38|hg38) build=38; p1s=10001; p1e=2781479; p2s=155701383; p2e=156030895 ;;
    *) echo "ERROR: chrX.male validation needs GRCh37 or GRCh38, got: $build" >&2; return 2 ;;
  esac
  [[ -s ${base}.pgen && -s ${base}.psam ]] || {
    echo "ERROR: incomplete chrX.male pfile: expected ${base}.{pgen,psam,pvar[.zst]}" >&2
    return 2
  }
  for pvar in "${base}.pvar" "${base}.pvar.zst"; do [[ -s $pvar ]] && break; done
  [[ -s $pvar ]] || { echo "ERROR: chrX.male pvar is missing for $base" >&2; return 2; }

  read -r variants bad_chr par_variants < <(
    gu_pvar_cat "$pvar" | awk -v p1s="$p1s" -v p1e="$p1e" -v p2s="$p2s" -v p2e="$p2e" '
      BEGIN{FS="[[:space:]]+"}
      /^#/ {next}
      NF>=2 {
        c=toupper($1); sub(/^CHR/,"",c); if(c=="23")c="X"
        n++; if(c!="X")bad++
        if(($2>=p1s&&$2<=p1e)||($2>=p2s&&$2<=p2e))par++
      }
      END{print n+0,bad+0,par+0}'
  )
  (( variants > 0 && bad_chr == 0 && par_variants == 0 )) || {
    echo "ERROR: $base is not pure non-PAR chrX: variants=$variants nonX=$bad_chr PAR=$par_variants" >&2
    return 2
  }

  read -r samples nonmale < <(
    awk 'BEGIN{FS="[[:space:]]+"}
      /^#/ {
        for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);if(x=="IID"||x=="SAMPLE")id=i;if(x=="SEX")sex=i}
        next
      }
      NF {
        n++; value=tolower($sex)
        if(!id||!sex||!(value=="1"||value=="m"||value=="male"))bad++
      }
      END{print n+0,bad+0}' "${base}.psam"
  )
  (( samples > 0 && nonmale == 0 )) || {
    echo "ERROR: $base is not male-only or its PSAM lacks IID/SEX: samples=$samples nonmale_or_unknown=$nonmale" >&2
    return 2
  }

  if [[ -s ${base}.smiss ]]; then
    hethap=$(awk 'BEGIN{FS="[[:space:]]+"}
      NR==1{for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);if(x=="HETHAP_CT")col=i};next}
      col{sum+=$col} END{if(col)print sum+0;else print "NA"}' "${base}.smiss")
    [[ $hethap == NA || $hethap == 0 ]] || {
      echo "ERROR: $base has $hethap heterozygous-haploid calls" >&2
      return 2
    }
  fi

  if [[ -s $stats ]]; then
    awk -F= -v build="GRCh$build" '
      {sub(/\r$/,"");v[$1]=$2}
      END{
        if(("build" in v)&&v["build"]!=build)exit 2
        if(("chromosome" in v)&&v["chromosome"]!="X")exit 2
        if(("par_regions_present" in v)&&toupper(v["par_regions_present"])!="FALSE")exit 2
        if(("nonmale_samples" in v)&&(v["nonmale_samples"]+0)!=0)exit 2
        if(("plink2_encoding" in v)&&v["plink2_encoding"]!="male_chrX_haploid")exit 2
      }' "$stats" || {
        echo "ERROR: chrX.male basic statistics disagree with the requested input/build: $stats" >&2
        return 2
      }
  fi
  printf 'chrX.male validation=PASS build=GRCh%s variants=%s male_samples=%s PAR_variants=0 heterozygous_haploid=%s source=%s\n' \
    "$build" "$variants" "$samples" "$hethap" "$base"
}

gu_validate_chrxy_pfile() {
  local base=${1:?gu_validate_chrxy_pfile requires PFILE_PREFIX}
  local build=${2:?gu_validate_chrxy_pfile requires BUILD}
  local pvar="" p1s p1e p2s p2e variants bad_chr outside_par samples
  read -r p1s p1e p2s p2e < <(gu_chrx_par_bounds_bed "$build") || return 2
  [[ -s ${base}.pgen && -s ${base}.psam ]] || {
    echo "ERROR: incomplete chrXY pfile: expected ${base}.{pgen,psam,pvar[.zst]}" >&2
    return 2
  }
  for pvar in "${base}.pvar" "${base}.pvar.zst"; do [[ -s $pvar ]] && break; done
  [[ -s $pvar ]] || { echo "ERROR: chrXY pvar is missing for $base" >&2; return 2; }
  read -r variants bad_chr outside_par < <(
    gu_pvar_cat "$pvar" | awk -v p1s="$p1s" -v p1e="$p1e" -v p2s="$p2s" -v p2e="$p2e" '
      BEGIN{FS="[[:space:]]+"}
      /^#/ {next}
      NF>=2 {
        c=toupper($1);sub(/^CHR/,"",c);if(c=="23")c="X"
        n++
        if(c!="PAR1"&&c!="PAR2"&&c!="X"&&c!="XY")bad++
        in1=($2>p1s&&$2<=p1e);in2=($2>p2s&&$2<=p2e)
        if(!(in1||in2))outside++
        if(c=="PAR1"&&!in1)outside++
        if(c=="PAR2"&&!in2)outside++
      }
      END{print n+0,bad+0,outside+0}'
  )
  (( variants > 0 && bad_chr == 0 && outside_par == 0 )) || {
    echo "ERROR: $base is not a pure PAR pfile: variants=$variants invalid_chrom=$bad_chr outside_PAR=$outside_par" >&2
    return 2
  }
  samples=$(awk 'BEGIN{FS="[[:space:]]+"}
    /^#/ {for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);if(x=="IID"||x=="SAMPLE")id=i};next}
    NF&&id{n++} END{print n+0}' "${base}.psam")
  (( samples > 0 )) || { echo "ERROR: chrXY PSAM lacks samples/IID: ${base}.psam" >&2; return 2; }
  printf 'chrXY validation=PASS build=GRCh%s variants=%s samples=%s outside_PAR=0 source=%s\n' \
    "${build#GRCh}" "$variants" "$samples" "$base"
}

gu_validate_chrx_male_vcf() {
  local vcf=${1:?gu_validate_chrx_male_vcf requires VCF}
  local build=${2:?gu_validate_chrx_male_vcf requires BUILD}
  local panel=${3:?gu_validate_chrx_male_vcf requires SAMPLE_PANEL}
  local p1s p1e p2s p2e contig region first work bad variants bad_record
  command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools is required to validate $vcf" >&2; return 2; }
  [[ -s $vcf && -s $panel ]] || { echo "ERROR: chrX.male VCF validation needs VCF and sample/sex metadata" >&2; return 2; }
  case "${build,,}" in
    37|b37|grch37|hg19) build=37; p1s=60001; p1e=2699520; p2s=154931044; p2e=155260560 ;;
    38|b38|grch38|hg38) build=38; p1s=10001; p1e=2781479; p2s=155701383; p2e=156030895 ;;
    *) echo "ERROR: chrX.male validation needs GRCh37 or GRCh38, got: $build" >&2; return 2 ;;
  esac
  contig=$(bcftools view -h "$vcf" | awk -F'[=,>]' '
    /^##contig=<ID=/{id=$3;if(id=="X"||id=="chrX"||id=="23"){print id;exit}}')
  [[ -n $contig ]] || { echo "ERROR: no X/chrX/23 contig in $vcf" >&2; return 2; }
  variants=$(bcftools index -n "$vcf") || { echo "ERROR: could not count indexed variants in $vcf" >&2; return 2; }
  [[ $variants =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: chrX.male VCF has no variants: $vcf" >&2; return 2; }
  bad_record=$(bcftools index -s "$vcf" | awk -v expected="$contig" '$3+0>0&&$1!=expected{print "non-X contig "$1;exit}') || {
    echo "ERROR: could not inspect indexed chromosome records in $vcf" >&2; return 2;
  }
  [[ -z $bad_record ]] || { echo "ERROR: chrX.male VCF violates the pure haploid-X contract: $bad_record" >&2; return 2; }
  # bcftools stats computes ploidy counts without emitting every GT. PSC columns
  # 4-6 count diploid SNP calls and PSI columns 8-11 count diploid indel calls;
  # a pure male X export must have zero in all of them and nonzero haploid SNPs.
  bad_record=$(bcftools stats -s - "$vcf" | awk '
    $1=="PSC"{n++;dip=$4+$5+$6;hap+=$12+$13;if(dip>0&&bad=="")bad="diploid SNP calls for sample "$3}
    $1=="PSI"{dip=$8+$9+$10+$11;if(dip>0&&bad=="")bad="diploid indel calls for sample "$3}
    END{if(!n&&bad=="")bad="no per-sample genotype statistics";if(!hap&&bad=="")bad="no nonmissing haploid SNP calls";if(bad!="")print bad}') || {
    echo "ERROR: could not inspect genotype ploidy statistics in $vcf" >&2; return 2;
  }
  [[ -z $bad_record ]] || { echo "ERROR: chrX.male VCF violates the pure haploid-X contract: $bad_record" >&2; return 2; }
  for region in "$contig:$p1s-$p1e" "$contig:$p2s-$p2e"; do
    # bcftools -r uses record-overlap semantics and can return a non-PAR SV
    # whose POS is outside PAR but INFO/END crosses the boundary. PLINK2 pfile
    # chromosome partitioning and --extract bed0 are POS-based, so test POS.
    first=$(bcftools query -r "$region" -f '%POS\n' "$vcf" | awk -F '[:-]' -v r="$region" '
      BEGIN{split(r,a,/[:-]/);lo=a[2]+0;hi=a[3]+0}$1>=lo&&$1<=hi{print;exit}')
    [[ -z $first ]] || { echo "ERROR: chrX.male VCF contains a PAR variant in $region" >&2; return 2; }
  done

  work=$(mktemp -d "${TMPDIR:-/tmp}/gu-x-vcf-check.XXXXXX")
  bcftools query -l "$vcf" | sort -u > "$work/vcf.samples"
  awk 'BEGIN{FS="[[:space:]]+"}
    NR==1{for(i=1;i<=NF;i++){x=tolower($i);sub(/^#/,"",x);if(x=="sample"||x=="iid")id=i;if(x=="sex")sex=i};next}
    NF{v=tolower($sex);if(v=="1"||v=="m"||v=="male")print $id}' "$panel" | sort -u > "$work/male.samples"
  bad=$(comm -3 "$work/vcf.samples" "$work/male.samples" | sed -n '1p')
  if [[ -n $bad || ! -s $work/vcf.samples ]]; then
    echo "ERROR: chrX.male VCF samples are not exactly the male sample set from $panel" >&2
    rm -rf -- "$work"
    return 2
  fi
  local samples
  samples=$(awk 'END{print NR+0}' "$work/vcf.samples")
  rm -rf -- "$work"
  printf 'chrX.male VCF validation=PASS build=GRCh%s variants=%s male_samples=%s PAR_variants=0 source=%s\n' \
    "$build" "$variants" "$samples" "$vcf"
}

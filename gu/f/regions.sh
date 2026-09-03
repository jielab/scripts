#!/usr/bin/env bash
# Shared GU region parsing. BED coordinates are 0-based, half-open.

gu_normalize_chrs() {
  local raw=${1:-} token
  raw=${raw//,/ }
  for token in $raw; do
    token=${token#chr}; token=${token#CHR}; token=${token^^}
    [[ $token == 23 ]] && token=X
    case "$token" in
      [1-9]|1[0-9]|2[0-2]|X) printf '%s\n' "$token" ;;
      *) echo "ERROR: unsupported chromosome in --chr: $token" >&2; return 2 ;;
    esac
  done | awk '!seen[$0]++' | paste -sd' ' -
}

gu_normalize_loci() {
  local input=${1:?gu_normalize_loci requires INPUT} output=${2:?gu_normalize_loci requires OUTPUT}
  [[ -s $input ]] || { echo "ERROR: --loci file is missing or empty: $input" >&2; return 2; }
  mkdir -p "$(dirname -- "$output")"
  awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"; bad=0; n=0}
    /^[[:space:]]*#/ || /^[[:space:]]*$/ {next}
    {
      chr=toupper($1); sub(/^CHR/,"",chr); if(chr=="23") chr="X"
      if(!(chr=="X" || (chr ~ /^[0-9]+$/ && (chr+0)>=1 && (chr+0)<=22))){print "ERROR: invalid chromosome at loci line " NR ": " $1 > "/dev/stderr"; bad=1; next}
      if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || ($3+0)<=($2+0)){
        print "ERROR: invalid BED interval at loci line " NR ": " $0 > "/dev/stderr"; bad=1; next
      }
      name=(NF>=4 && $4!="" && $4!=".") ? $4 : sprintf("locus_%s_%s_%s",chr,$2,$3)
      key=chr SUBSEP $2 SUBSEP $3 SUBSEP name
      if(!seen[key]++){print chr,$2+0,$3+0,name; n++}
    }
    END{
      if(bad) exit 2
      if(n==0){print "ERROR: no valid loci found" > "/dev/stderr"; exit 2}
    }' "$input" > "$output"
}

gu_loci_chrs() {
  awk 'BEGIN{FS="\t"} $1!~/^#/ && NF>=3 && !seen[$1]++{print $1}' "$1" | paste -sd' ' -
}

gu_loci_regions_for_chr() {
  local file=${1:?} chr=${2:?} prefix=${3:-}
  awk -F'\t' -v chr="$chr" -v prefix="$prefix" '$1==chr{printf "%s%s:%d-%d%s",sep,prefix chr,$2+1,$3;sep=","} END{print ""}' "$file"
}

gu_scope_id() {
  local loci=${1:-} chrs=${2:-} name_source=${3:-} base sum
  if [[ -n $loci ]]; then
    [[ -n $name_source ]] || name_source=$loci
    base=$(basename -- "$name_source"); base=${base%.*}; base=$(printf '%s' "$base" | tr -c 'A-Za-z0-9._-' '_')
    sum=$(cksum < "$loci" | awk '{print $1}')
    printf 'loci-%s-%s\n' "$base" "$sum"
  elif [[ -n $chrs ]]; then
    printf 'chr-%s\n' "$(printf '%s' "$chrs" | tr ' ' '_')"
  else
    printf 'genome\n'
  fi
}

# Human-readable label for public result directories and log names. The
# content checksum from gu_scope_id remains in hidden cache paths/provenance.
gu_scope_label() {
  local loci=${1:-} chrs=${2:-} name_source=${3:-} base
  if [[ -n $loci ]]; then
    [[ -n $name_source ]] || name_source=$loci
    base=$(basename -- "$name_source"); base=${base%.*}
    base=$(printf '%s' "$base" | tr -c 'A-Za-z0-9._-' '_')
    printf '%s\n' "${base:-loci}"
  elif [[ -n $chrs ]]; then
    printf 'chr%s\n' "$(printf '%s' "$chrs" | tr ' ' '_')"
  else
    printf 'genome\n'
  fi
}

# Build one public result directory per analysis unit.
gu_chr_result_dir() {
  local analysis_root=${1:?} method=${2:?} unit_label=${3:?}
  local target_namespace=${4:-} override=${5:-} request_unit_count=${6:-1} out
  if [[ -n $override ]]; then
    out=$override
    (( request_unit_count <= 1 )) || out=$out/$unit_label
  else
    out=$analysis_root/$method
    [[ -z $target_namespace ]] || out=$out/$target_namespace
    out=$out/$unit_label
  fi
  printf '%s\n' "$out"
}

gu_locus_unit_label() {
  local map_file=${1:?} chr=${2:?} core_start=${3:?} core_end=${4:?} locus=${5:?}
  local safe duplicates
  safe=$(printf '%s' "$locus" | tr -c 'A-Za-z0-9._-' '_'); safe=${safe##_}; safe=${safe%%_}
  [[ -n $safe ]] || safe=locus_${core_start}_${core_end}
  duplicates=$(awk -F'\t' -v c="$chr" -v l="$locus" 'NR>1&&$1==c&&$6==l{n++}END{print n+0}' "$map_file")
  if (( duplicates > 1 )); then printf 'chr%s_%s_%s_%s\n' "$chr" "$safe" "$core_start" "$core_end"
  else printf 'chr%s_%s\n' "$chr" "$safe"
  fi
}

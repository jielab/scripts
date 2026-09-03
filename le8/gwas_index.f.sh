#!/usr/bin/env bash
# Shared BGZF/tabix helpers for standardized GWAS files.

gwas_index_log() {
  if declare -F gwas_clean_log >/dev/null 2>&1; then
    gwas_clean_log "$*"
  else
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
  fi
}

gwas_index_header() {
  local file="$1"
  if declare -F gwas_clean_zcat >/dev/null 2>&1; then
    (set +o pipefail; gwas_clean_zcat "$file" | head -n 1 | tr -d '\r' || true)
  else
    (set +o pipefail; gzip -cd -- "$file" | head -n 1 | tr -d '\r' || true)
  fi
}

gwas_index_col() {
  local header="$1" wanted="$2"
  awk -F '\t' -v wanted="$wanted" 'NR==1{for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);gsub(/[^A-Z0-9]/,"",x);if(x==wanted){print i;exit}}}' <<< "$header"
}

gwas_index_valid() {
  local file="$1" index=""
  [[ -s "${file}.tbi" ]] && index="${file}.tbi"
  [[ -z "$index" && -s "${file}.csi" ]] && index="${file}.csi"
  [[ -n "$index" && ! "$file" -nt "$index" ]] || return 1
  tabix -l "$file" >/dev/null 2>&1
}

gwas_index_make_sidecar() {
  local file="$1" chr_col="$2" pos_col="$3"
  rm -f -- "${file}.tbi" "${file}.csi"
  if tabix -f -s "$chr_col" -b "$pos_col" -e "$pos_col" -S 1 "$file" >/dev/null 2>&1; then
    return 0
  fi
  rm -f -- "${file}.tbi" "${file}.csi"
  tabix -f -C -s "$chr_col" -b "$pos_col" -e "$pos_col" -S 1 "$file" >/dev/null 2>&1
}

gwas_index_file() {
  local target="$1" header chr_col pos_col tmp index_ext sort_tmp
  [[ -s "$target" ]] || { echo "ERROR: cannot index missing/empty GWAS: $target" >&2; return 1; }
  command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip is required to index GWAS files" >&2; return 127; }
  command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix is required to index GWAS files" >&2; return 127; }

  header=$(gwas_index_header "$target")
  [[ -n "$header" ]] || { echo "ERROR: unreadable GWAS header: $target" >&2; return 1; }
  chr_col=$(gwas_index_col "$header" CHR)
  [[ -n "$chr_col" ]] || chr_col=$(gwas_index_col "$header" CHROM)
  pos_col=$(gwas_index_col "$header" POS)
  [[ -n "$pos_col" ]] || pos_col=$(gwas_index_col "$header" BP)
  [[ -n "$chr_col" && -n "$pos_col" ]] || {
    echo "ERROR: GWAS needs CHR and POS columns for tabix: $target" >&2
    return 1
  }

  # Normal gwas_post output is already sorted BGZF, so this only writes the
  # small sidecar. Plain gzip or an unsorted legacy file takes the fallback.
  if gwas_index_make_sidecar "$target" "$chr_col" "$pos_col" && gwas_index_valid "$target"; then
    gwas_index_log "tabix index: $target"
    return 0
  fi
  rm -f -- "${target}.tbi" "${target}.csi"

  tmp="${target%.gz}.index.tmp.${BASHPID:-$$}.gz"
  sort_tmp="${GWAS_POST_TMP:-${GWAS_INDEX_TMP_DIR:-${TMPDIR:-/tmp}}}"
  mkdir -p "$sort_tmp"
  rm -f -- "$tmp" "${tmp}.tbi" "${tmp}.csi"
  gwas_index_log "normalize coordinate order/BGZF before tabix: $target"
  if declare -F gwas_clean_zcat >/dev/null 2>&1; then
    gwas_clean_zcat "$target"
  else
    gzip -cd -- "$target"
  fi | {
    IFS= read -r first_header
    printf '%s\n' "$first_header"
    sort -T "$sort_tmp" -S "${GWAS_INDEX_SORT_MEMORY:-512M}" -t $'\t' \
      -k"${chr_col}","${chr_col}"V -k"${pos_col}","${pos_col}"n
  } | bgzip -@ "${GWAS_CLEAN_COMP_THREADS:-1}" -c > "$tmp"

  bgzip -t "$tmp"
  gwas_index_make_sidecar "$tmp" "$chr_col" "$pos_col"
  gwas_index_valid "$tmp" || {
    echo "ERROR: tabix validation failed: $target" >&2
    rm -f -- "$tmp" "${tmp}.tbi" "${tmp}.csi"
    return 1
  }
  chmod --reference="$target" "$tmp" 2>/dev/null || true
  if [[ -s "${tmp}.tbi" ]]; then index_ext=tbi; else index_ext=csi; fi
  mv -f -- "$tmp" "$target"
  rm -f -- "${target}.tbi" "${target}.csi"
  mv -f -- "${tmp}.${index_ext}" "${target}.${index_ext}"
  gwas_index_valid "$target"
  gwas_index_log "BGZF+$index_ext ready: $target"
}

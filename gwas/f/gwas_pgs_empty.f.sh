#!/usr/bin/env bash
# A successful COJO run may have no score variants. Never infer this from a
# missing score file alone: validate the trait, settings and explicit status.
if declare -F gwas_post_pgs >/dev/null 2>&1 &&
   ! declare -F gwas_post_pgs_with_scores >/dev/null 2>&1; then
  eval "$(declare -f gwas_post_pgs | sed '1s/gwas_post_pgs/gwas_post_pgs_with_scores/')"
fi

gwas_post_pgs() {
  [[ ",$PGS_STEPS," == *,pgs,* || "$PGS_STEPS" == all ]] || return 0
  local reason=""
  if [[ ! -s "$PGS_SCORE_FILE" ]] &&
     gwas_lead_marker_matches "$COJO_DONE" "$GWAS" cojo "$P_LEAD" "$LEAD_WINDOW" "$CHRS"; then
    reason=$(awk -F '\t' 'NR==2{print $6}' "$COJO_DONE")
    case "$reason" in
      no_significant_variants)
        if [[ ! -s "$AWK_SNP" ]] || ! awk 'NR>1{found=1}END{exit found?1:0}' "$AWK_SNP"; then
          reason=""
        fi;;
      no_reference_matched_variants|no_snps_selected) ;;
      *) reason="";;
    esac
  fi
  if [[ -z "$reason" ]]; then
    gwas_post_pgs_with_scores
    return
  fi
  mkdir -p "$PGS_DIR"
  rm -f -- "$PGS_DONE"
  # Match the established no-matched-variants schema; no invented zero scores.
  printf '#IID\t%s.ALLELE_CT\t%s.SCORE_SUM\n' "$GWAS" "$GWAS" | gzip -c > "${PGS_OUTPUT}.tmp.$$"
  mv -f -- "${PGS_OUTPUT}.tmp.$$" "$PGS_OUTPUT"
  {
    printf 'key\tvalue\ngwas\t%s\nsource\t%s\nmatched_variants\tnone\nreason\t%s\n' "$GWAS" "$PGS_SCORE_FILE" "$reason"
    printf 'pfile_dir\t%s\nchromosomes\tnone\noutput\t%s\nchromosome_intermediates\tnone\n' "$PGS_PFILE_DIR" "$PGS_OUTPUT"
  } > "${PGS_META}.tmp.$$"
  mv -f -- "${PGS_META}.tmp.$$" "$PGS_META"
  printf 'GWAS\tSTATUS\tTIME\n%s\tcomplete_no_score_variants\t%s\n' "$GWAS" "$(date '+%F %T')" > "${PGS_DONE}.tmp.$$"
  mv -f -- "${PGS_DONE}.tmp.$$" "$PGS_DONE"
  gwas_post_log "PGS skipped: $reason (header-only output; no individual scores)"
}

#!/usr/bin/env bash
# Shared lead reference storage and indexed lookup. Loaded once by --worker;
# generated commands must not redefine these functions after sourcing it.
# GWAS_LEAD_* values are consumed by the generated chromosome loop.
# shellcheck disable=SC2034

gwas_post_reference_file() {
	local ref="$1" kind="$2" ext
	case "$kind:$ref" in
		pvar:*.pvar|pvar:*.pvar.gz|pvar:*.pvar.bgz|pvar:*.pvar.zst|bim:*.bim|bim:*.bim.gz|bim:*.bim.bgz)
			gwas_post_need_file "$ref" || return
			printf '%s\n' "$ref"; return ;;
	esac
	local -a suffixes
	case "$kind" in
		pvar) suffixes=(.pvar .pvar.gz .pvar.bgz .pvar.zst) ;;
		bim) suffixes=(.bim .bim.gz .bim.bgz) ;;
		*) echo "ERROR: invalid lead reference kind: $kind" >&2; return 2 ;;
	esac
	for ext in "${suffixes[@]}"; do
		if [[ -s "${ref}${ext}" ]]; then printf '%s\n' "${ref}${ext}"; return; fi
	done
	echo "ERROR: no $kind reference found: $ref" >&2
	return 1
}

# Stage immutable groups lazily, once a chromosome has candidates. The Python
# helper locks each version and publishes it atomically; all traits/populations
# can share the same full pfile. Cache the population keep file separately.
gwas_post_stage_lead_references() {
	local ref="$1" cref="$2" assoc="$3" ma="$4" result pvar
	local -a files staged
	GWAS_LEAD_CLUMP_REF="$ref"
	GWAS_LEAD_COJO_REF="$cref"
	[[ "${GWAS_LOCAL_IO:-FALSE}" == TRUE ]] || return 0
	if gwas_post_has_data_rows "$assoc"; then
		pvar=$(gwas_post_reference_file "$ref" pvar) || return
		files=("${ref}.pgen" "$pvar" "${ref}.psam")
		[[ ! -f "${ref}.pgen.pgi" ]] || files+=("${ref}.pgen.pgi")
		result=$(python3 "$GWAS_STORAGE_HELPER" stage --root "$GWAS_LOCAL_CACHE/references/lead-genotypes" "${files[@]}") || return
		mapfile -t staged <<<"$result"
		GWAS_LEAD_CLUMP_REF="${staged[0]%.pgen}"
		if [[ -n "${REFGEN_KEEP:-}" && "${GWAS_LEAD_KEEP_STAGED:-FALSE}" != TRUE ]]; then
			REFGEN_KEEP=$(python3 "$GWAS_STORAGE_HELPER" stage --root "$GWAS_LOCAL_CACHE/references/lead-keep" "$REFGEN_KEEP") || return
			GWAS_LEAD_KEEP_STAGED=TRUE
		fi
	fi
	if gwas_post_has_data_rows "$ma"; then
		result=$(python3 "$GWAS_STORAGE_HELPER" stage --root "$GWAS_LOCAL_CACHE/references/lead-genotypes" \
			"${cref}.bed" "${cref}.bim" "${cref}.fam") || return
		mapfile -t staged <<<"$result"
		GWAS_LEAD_COJO_REF="${staged[0]%.bed}"
	fi
}

gwas_post_reference_cache() {
	python3 "$GWAS_STORAGE_HELPER" lead-index --root "$LEAD_REF_CACHE" \
		--reference "$1" --kind "$2" --sort-memory "${GWAS_POST_SORT_MEMORY:-512M}"
}

gwas_post_subset_match_reference() {
	local ref="$1" query="$2" out="$3" kind="$4" source cached regions tmp
	source=$(gwas_post_reference_file "$ref" "$kind") || return
	cached=$(gwas_post_reference_cache "$source" "$kind") || return
	regions=$(mktemp "$GWAS_POST_TMP/reference.$kind.XXXXXXXX.tsv") || return
	tmp="${out}.tmp.$$"
	# A file of unique 1-based intervals avoids ARG_MAX for large hit lists and
	# avoids returning the same reference row repeatedly for duplicate queries.
	if ! awk -v FS='[ \t]+' -v OFS='\t' '
		function normchr(x){x=toupper(x);sub(/^CHR/,"",x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x==""?"0":x}
		NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);if(h=="CHR")c=i;else if(h=="POS")p=i}next}
		c&&p&&$p~/^[0-9]+$/&&$p+0>0{pos=sprintf("%.0f",$p+0);print normchr($c),pos,pos}
		END{if(!c||!p){print "ERROR: lead query requires CHR and POS columns" > "/dev/stderr";exit 2}}
	' "$query" | LC_ALL=C sort -u -T "$GWAS_POST_TMP" -S "${GWAS_POST_SORT_MEMORY:-512M}" >"$regions"; then
		rm -f "$regions"; return 1
	fi
	if [[ "$kind" == pvar ]]; then printf '#CHROM\tPOS\tID\tREF\tALT\n' >"$tmp"; else : >"$tmp"; fi
	if [[ -s "$regions" ]] && ! tabix -R "$regions" "$cached" >>"$tmp"; then
		rm -f "$tmp" "$regions"; return 1
	fi
	# match_SNP expects a nonempty BIM even when every query is unmatched.
	if [[ "$kind" == bim && ! -s "$tmp" ]]; then printf '0\t.\t0\t0\tN\tN\n' >"$tmp"; fi
	mv -f -- "$tmp" "$out"
	rm -f -- "$regions"
}

gwas_post_filter_clump_multiallelic() {
	local ref="$1" tag="$2" assoc="$3" pvar cached bad filtered n_before n_after n_removed
	pvar=$(gwas_post_reference_file "$ref" pvar) || return
	cached=$(gwas_post_reference_cache "$pvar" pvar) || return
	bad="${cached%/*}/ambiguous.snp"
	n_before=$(awk 'NR>1{n++} END{print n+0}' "$assoc")
	filtered="${assoc}.filtered.$$"
	awk -v FS='\t' -v OFS='\t' '
		FILENAME==ARGV[1]{bad[$1]=1;next}
		FNR==1{print;next}
		!($1 in bad){print}
	' "$bad" "$assoc" >"$filtered" || return
	mv -f -- "$filtered" "$assoc"
	n_after=$(awk 'NR>1{n++} END{print n+0}' "$assoc")
	n_removed=$((n_before-n_after))
	printf 'GWAS\tCHR\tN_BEFORE\tN_REMOVED_AMBIGUOUS_REF_ID\tN_AFTER\n%s\t%s\t%s\t%s\t%s\n' \
		"$GWAS" "$tag" "$n_before" "$n_removed" "$n_after" >"${QC_PREFIX}.${tag}.clump_multiallelic.tsv"
	gwas_post_log "clump ambiguous-reference-ID filter $tag: before=$n_before removed=$n_removed after=$n_after cache=$bad"
}

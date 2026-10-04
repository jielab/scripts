#!/usr/bin/env bash


# 🚩 GWAS coordinate indexing
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
		(
			set +o pipefail
			gwas_clean_zcat "$file" | head -n 1 | tr -d '\r' || true
		)
	else
		(
			set +o pipefail
			gzip -cd -- "$file" | head -n 1 | tr -d '\r' || true
		)
	fi
}

gwas_index_col() {
	local header="$1" wanted="$2"
	awk -F '\t' -v wanted="$wanted" 'NR==1{for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);gsub(/[^A-Z0-9]/,"",x);if(x==wanted){print i;exit}}}' <<<"$header"
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
	[[ -s "$target" ]] || {
		echo "ERROR: cannot index missing/empty GWAS: $target" >&2
		return 1
	}
	command -v bgzip >/dev/null 2>&1 || {
		echo "ERROR: bgzip is required to index GWAS files" >&2
		return 127
	}
	command -v tabix >/dev/null 2>&1 || {
		echo "ERROR: tabix is required to index GWAS files" >&2
		return 127
	}

	header=$(gwas_index_header "$target")
	[[ -n "$header" ]] || {
		echo "ERROR: unreadable GWAS header: $target" >&2
		return 1
	}
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
	} | bgzip -@ "${GWAS_CLEAN_COMP_THREADS:-1}" -c >"$tmp"

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

# LE8 loads only the indexing functions; format workers load the full module.
if [[ ${1:-} == --index ]]; then return 0; fi


# 🚩 Shared completion checks
# Shared thin completion checks for the coordinator and workers.
# Fingerprints use paths, sizes and nanosecond mtimes; no full GWAS scan.

# Plot-only runs consume an existing artifact, independently of the producer's
# path-sensitive cache marker. Never regenerate data as a plotting side effect.
gwas_thin_plot_only() {
	[[ ",$1," == *,mplot,* && "$1" != all && ! ",$1," =~ ,(format|thin|liftover), ]]
}

gwas_thin_plot_ready() {
	local input="$1" output="$2"
	[[ -s "$input" && -s "$output" && -s "$output.tbi" && ! "$input" -nt "$output" && ! "$output" -nt "$output.tbi" ]] || return 1
	tabix -l "$output" >/dev/null 2>&1
}

gwas_thin_signature() {
	local input="$1" output="$2" grch="$3" cap="$4" hm3="$5" pos="$6" script="$7" phe="$8" f
	local facts
	facts=$(
		printf 'thin_state_v1\t%s\t%s\n' "$grch" "$cap"
		for f in "$input" "$output" "$output.tbi" "$hm3" "$script" "${script%/*}/format.thin.awk" "$phe"; do
			[[ -s "$f" ]] || exit 1
			stat -Lc '%n|%s|%y' -- "$f" || exit 1
		done
		if [[ -n "$pos" ]]; then
			[[ -s "$pos" ]] || exit 1
			stat -Lc '%n|%s|%y' -- "$pos" || exit 1
		fi
	) || return 1
	printf '%s\n' "$facts" | sha256sum | cut -d ' ' -f 1
}

gwas_thin_complete() {
	local input="$1" output="$2" grch="$3" cap="$4" hm3="$5" pos="$6" script="$7" phe="$8"
	local signature saved f
	[[ -s "$input" && -s "$output" && -s "$output.tbi" && "$output" -nt "$input" && ! "$output" -nt "$output.tbi" ]] || return 1
	[[ -s "$output.done" ]] || return 1
	signature=$(gwas_thin_signature "$input" "$output" "$grch" "$cap" "$hm3" "$pos" "$script" "$phe") || return 1
	saved=$(cat -- "$output.done") || return 1
	[[ "$saved" == "$signature" ]] || return 1
	tabix -l "$output" >/dev/null 2>&1
}


# 🚩 GWAS formatting and tabix helpers


# 🚩 GWAS logging and file helpers
gwas_clean_log() { echo "[$(date '+%F %T')] [$GWAS] $*"; }
gwas_clean_gzip_ok() { [[ -s "$1" ]] && gzip -t "$1" >/dev/null 2>&1; }
gwas_clean_need_file() { [[ -s "$1" ]] || {
	echo "ERROR: missing or empty file: $1" >&2
	exit 1
}; }

gwas_clean_zcat() {
	local f="$1"
	case "$f" in
		*.gz | *.bgz)
			if command -v pigz >/dev/null 2>&1; then
				pigz -dc -p "${GWAS_CLEAN_DECOMP_THREADS:-1}" "$f"
			else
				gzip -dc "$f"
			fi
			;;
		*.zst)
			if command -v zstdcat >/dev/null 2>&1; then
				zstdcat -- "$f"
			elif command -v zstd >/dev/null 2>&1; then
				zstd -q -dc -- "$f"
			else
				echo "ERROR: zstd or zstdcat is required to read: $f" >&2
				return 127
			fi
			;;
		*)
			cat "$f"
			;;
	esac
}

gwas_clean_compress() {
	if command -v bgzip >/dev/null 2>&1; then
		bgzip -@ "${GWAS_CLEAN_COMP_THREADS:-1}" -c
	elif command -v pigz >/dev/null 2>&1; then
		pigz -c -p "${GWAS_CLEAN_COMP_THREADS:-1}"
	else
		gzip -c
	fi
}

gwas_clean_detect_fs() {
	local f="$1" line
	line=$(
		set +o pipefail
		gwas_clean_zcat "$f" | head -n 1 | tr -d '\r' || true
	)
	case "$line" in
		*$'\t'*) printf '\t' ;;
		*) printf '[ \t]+' ;;
	esac
}

gwas_clean_header_line() {
	local f="$1"
	(
		set +o pipefail
		gwas_clean_zcat "$f" | head -n 1 | tr -d '\r' || true
	)
}

gwas_clean_col() {
	[[ "${1:-}" =~ ^[0-9]+$ ]] && printf '%s\n' "$1" || printf '0\n'
}

gwas_clean_load_phef() {
	gwas_clean_need_file "$PHEF"
	if ! command -v dos2unix >/dev/null 2>&1; then
		dos2unix() { sed 's/\r$//'; }
	fi
	source "$PHEF"
}

gwas_clean_header_names() {
	local src="$1" out="$2" pipefail_on=FALSE errexit_on=FALSE rc=0
	gwas_clean_load_phef
	if set -o | awk '$1=="errexit" && $2=="on"{found=1} END{exit !found}'; then
		errexit_on=TRUE
	fi
	if set -o | awk '$1=="pipefail" && $2=="on"{found=1} END{exit !found}'; then
		pipefail_on=TRUE
	fi
	set +e
	set +o pipefail
	if phe_header_names "$src" >"$out"; then
		rc=0
	else
		rc=$?
	fi
	[[ "$pipefail_on" == "TRUE" ]] && set -o pipefail
	[[ "$errexit_on" == "TRUE" ]] && set -e
	return "$rc"
}


# 🚩 Small GWAS formatting
std_small() {
	local src="$1" out="$2" header_out="$3" tmp bad_snp cand_raw cand_uniq input_fs t0 n_raw n_uniq
	gwas_clean_need_file "$src"
	gwas_clean_need_file "$HM3"
	t0=$(date +%s)
	gwas_clean_log "Small GWAS: $src -> $out"

	gwas_clean_header_names "$src" "$header_out"
	SNP_col=$(gwas_clean_col "${SNP_col:-}")
	CHR_col=$(gwas_clean_col "${CHR_col:-}")
	POS_col=$(gwas_clean_col "${POS_col:-}")
	EA_col=$(gwas_clean_col "${EA_col:-}")
	NEA_col=$(gwas_clean_col "${NEA_col:-}")
	EAF_col=$(gwas_clean_col "${EAF_col:-}")
	N_col=$(gwas_clean_col "${N_col:-}")
	BETA_col=$(gwas_clean_col "${BETA_col:-}")
	SE_col=$(gwas_clean_col "${SE_col:-}")
	P_col=$(gwas_clean_col "${P_col:-}")
	LOG10P_col=$(gwas_clean_col "${LOG10P_col:-}")

	[[ "$CHR_col" -gt 0 && "$POS_col" -gt 0 ]] || {
		echo "ERROR: CHR/POS not detected: $src" >&2
		exit 1
	}
	[[ "$P_col" -gt 0 || "$LOG10P_col" -gt 0 ]] || {
		echo "ERROR: no P or LOG10P detected: $src" >&2
		exit 1
	}

	tmp="${out}.tmp.$$"
	bad_snp="${out}.bad.snp.$$"
	cand_raw="${out}.candidate.raw.$$"
	cand_uniq="${out}.candidate.uniq.$$"
	: >"$bad_snp"
	if [[ -n "${REFGEN_CLUMP:-}" ]]; then
		while read -r _ _ _ bad; do
			[[ -s "$bad" ]] && cat "$bad" >>"$bad_snp"
		done < <(ref_clump_pfiles "$REFGEN_CLUMP")
	fi

	input_fs=$(gwas_clean_detect_fs "$src")
	gwas_clean_header_line "$src" >"$cand_raw"

	gwas_clean_log "prefilter HM3/P/cis candidates with one-pass awk: $src"
	gwas_clean_zcat "$src" | awk -v FS="$input_fs" \
		-v hm3_file="$HM3" -v pthr="$P_SMALL" \
		-v cis_bed="${CIS_BED:-}" -v cis_name="$GWAS" -v flank="${CIS_FLANK:-0}" \
		-v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
		-v p_col="$P_col" -v logp_col="$LOG10P_col" '
    function get(c, x){x=(c>0 ? $c : "");gsub(/^[[:space:]]+|[[:space:]]+$/,"",x);return x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x); if(x=="X") x="23"; if(x=="Y") x="24"; if(x=="MT"||x=="M") x="25"; return x}
    function keep_p(  p,lp){
      if(pthr=="" || pthr+0<=0) return 0
      if(p_col>0 && isnum(get(p_col))) return get(p_col)+0 <= pthr+0
      if(logp_col>0 && isnum(get(logp_col))) return get(logp_col)+0 >= logpthr
      return 0
    }
    function in_cis(chr,pos,  i){for(i=1;i<=n_cis;i++){if(chr==rchr[i] && pos>=rstart[i] && pos<=rend[i]) return 1} return 0}
    BEGIN{
      logpthr=(pthr+0>0 ? -log(pthr+0)/log(10) : 0)
      while((getline line < hm3_file)>0){
        gsub(/\r/,"",line)
        split(line,a,/[ \t]+/)
        if(a[1]!="" && a[1]!="SNP") hm3[a[1]]=1
      }
      close(hm3_file)
      if(cis_bed!=""){
        while((getline line < cis_bed)>0){
          gsub(/\r/,"",line); if(line=="" || line ~ /^#/) continue
          split(line,a,/[ \t]+/)
          if(a[4]==cis_name){n_cis++; rchr[n_cis]=normchr(a[1]); rstart[n_cis]=a[2]-flank; if(rstart[n_cis]<0) rstart[n_cis]=0; rend[n_cis]=a[3]+flank}
        }
        close(cis_bed)
      }
    }
    NR==1{next}
    {
      snp=get(snp_col)
      if((snp!="" && (snp in hm3)) || keep_p()){
        print
        next
      }
      if(n_cis>0){
        chr=normchr(get(chr_col)); pos=get(pos_col)+0
        if(chr!="" && pos>0 && in_cis(chr,pos)) print
      }
    }' >>"$cand_raw"

	awk 'NR==1{print; next} !seen[$0]++' "$cand_raw" >"$cand_uniq"
	n_raw=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$cand_raw")
	n_uniq=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$cand_uniq")
	gwas_clean_log "prefilter candidates: raw=$n_raw unique=$n_uniq"

	awk -v FS="$input_fs" -v OFS='\t' \
		-v hm3_file="$HM3" -v pthr="$P_SMALL" \
		-v bad_snp_file="$bad_snp" \
		-v cis_bed="${CIS_BED:-}" -v cis_name="$GWAS" -v flank="${CIS_FLANK:-0}" \
		-v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
		-v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
		-v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" '
    function get(c, x){x=(c>0 ? $c : "");gsub(/^[[:space:]]+|[[:space:]]+$/,"",x);return x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x); if(x=="X") x="23"; if(x=="Y") x="24"; if(x=="MT"||x=="M") x="25"; return x}
    function snpid(chr,pos,  s,ea,nea){s=get(snp_col); if(s==""||s=="NA"||s=="."){ea=get(ea_col); nea=get(nea_col); s=chr":"pos; if(ea!="") s=s":"ea; if(nea!="") s=s":"nea} return s}
    function pvalue(){if(p_col>0 && isnum(get(p_col))) return get(p_col)+0; if(logp_col>0 && isnum(get(logp_col))) return 10^(-(get(logp_col)+0)); return ""}
    function logpvalue(p){if(logp_col>0 && isnum(get(logp_col))) return get(logp_col); if(p!="" && p>0) return -log(p)/log(10); return ""}
    function in_cis(chr,pos,  i){for(i=1;i<=n_cis;i++){if(chr==rchr[i] && pos>=rstart[i] && pos<=rend[i]) return 1} return 0}
    BEGIN{
      while((getline line < hm3_file)>0){gsub(/\r/,"",line); split(line,a,/[ \t]+/); if(a[1]!="" && a[1]!="SNP") hm3[a[1]]=1}
      close(hm3_file)
      while((getline line < bad_snp_file)>0){gsub(/\r/,"",line); split(line,a,/[ \t]+/); if(a[1]!="") bad_snp[a[1]]=1}
      close(bad_snp_file)
      if(cis_bed!=""){
        while((getline line < cis_bed)>0){
          gsub(/\r/,"",line); if(line=="" || line ~ /^#/) continue
          split(line,a,/[ \t]+/)
          if(a[4]==cis_name){n_cis++; rchr[n_cis]=normchr(a[1]); rstart[n_cis]=a[2]-flank; if(rstart[n_cis]<0) rstart[n_cis]=0; rend[n_cis]=a[3]+flank}
        }
        close(cis_bed)
      }
    }
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P"; next}
    {
      chr=normchr(get(chr_col)); pos=get(pos_col)+0
      if(chr=="" || pos<=0) next
      snp=snpid(chr,pos)
      if(snp in bad_snp) next
      p=pvalue(); lp=logpvalue(p)
      if((snp in hm3) || (p!="" && p<=pthr) || (n_cis>0 && in_cis(chr,pos))){
        print snp,chr,pos,get(ea_col),get(nea_col),get(eaf_col),get(n_col),get(beta_col),get(se_col),p,lp
      }
    }' "$cand_uniq" |
		{
			IFS= read -r format_header
			printf '%s\n' "$format_header"
			sort -t $'\t' -k2,2V -k3,3n
		} | gwas_clean_compress >"$tmp"

	gzip -t "$tmp"
	rm -f "$bad_snp" "$cand_raw" "$cand_uniq"
	mv -f "$tmp" "$out"
	gwas_clean_zcat "$out" | awk -v g="$GWAS" 'NR==1{next} END{print g"\t"NR-1}' >"${QC_PREFIX}.small.nrow.tsv"
	gwas_clean_log "Small GWAS done in $(($(date +%s) - t0)) sec"
}

gwas_clean_copy_small_to_final() {
	if [[ "$REPLACE" != "TRUE" ]] && gwas_clean_gzip_ok "$FINAL"; then
		gwas_clean_log "Final clean GWAS exists: $FINAL"
		return 0
	fi
	gwas_clean_need_file "$SMALL"
	gwas_clean_log "No liftOver: $SMALL -> $FINAL"
	tmp="${FINAL}.tmp.$$"
	cp -f "$SMALL" "$tmp"
	gzip -t "$tmp"
	mv -f "$tmp" "$FINAL"
}


# 🚩 Genome-build conversion
gwas_clean_liftover_small_to_final() {
	local work n_small n_lift n_unmap
	if [[ "$REPLACE" != "TRUE" ]] && gwas_clean_gzip_ok "$FINAL"; then
		gwas_clean_log "Lifted final clean GWAS exists: $FINAL"
		return 0
	fi
	gwas_clean_need_file "$SMALL"
	gwas_clean_need_file "$CHAIN"
	command -v "$LIFTOVER_BIN" >/dev/null 2>&1 || {
		echo "ERROR: liftOver not found: $LIFTOVER_BIN" >&2
		exit 1
	}
	gwas_clean_log "liftOver Small GWAS: $SMALL -> $FINAL"

	work="$(dirname "$FINAL")/.${GWAS}.liftover.$$"
	rm -rf "$work"
	mkdir -p "$work"
	trap 'rm -rf "$work"' EXIT

	gwas_clean_zcat "$SMALL" | awk -v FS='\t' -v OFS='\t' -v bed="$work/to_lift.bed" '
    NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "row_id","SNP","CHR_OLD","POS_OLD","EA","NEA","EAF","N","BETA","SE","P","LOG10P"; next}
    NR>1 && ("CHR" in c) && ("POS" in c){
      row=NR-1; chr=$(c["CHR"]); pos=$(c["POS"])+0
      if(chr=="" || pos<=0) next
      print "chr"chr,pos-1,pos,row > bed
      print row,$(c["SNP"]),chr,pos,$(c["EA"]),$(c["NEA"]),$(c["EAF"]),$(c["N"]),$(c["BETA"]),$(c["SE"]),$(c["P"]),$(c["LOG10P"])
    }' >"$work/std.old.tsv"

	"$LIFTOVER_BIN" "$work/to_lift.bed" "$CHAIN" "$work/lifted.bed" "$work/unmapped.bed" >"${QC_PREFIX}.liftover.log" 2>&1

	awk -v OFS='\t' '{chr=$1; sub(/^chr/,"",chr); if(chr=="X") chr="23"; if(chr=="Y") chr="24"; if(chr=="MT"||chr=="M") chr="25"; print $4,chr,$3}' "$work/lifted.bed" | sort -k1,1n >"$work/lifted.pos.tsv"
	tail -n +2 "$work/std.old.tsv" | sort -k1,1n >"$work/std.old.sorted.tsv"

	printf 'SNP\tCHR\tPOS\tEA\tNEA\tEAF\tN\tBETA\tSE\tP\tLOG10P\n' >"$work/final.tsv"
	join -t $'\t' -1 1 -2 1 "$work/lifted.pos.tsv" "$work/std.old.sorted.tsv" |
		awk 'BEGIN{FS=OFS="\t"}{print $4,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14}' >>"$work/final.tsv"

	{
		IFS= read -r final_header
		printf '%s\n' "$final_header"
		sort -t $'\t' -k2,2V -k3,3n
	} <"$work/final.tsv" | gwas_clean_compress >"${FINAL}.tmp.$$"
	gzip -t "${FINAL}.tmp.$$"
	mv -f "${FINAL}.tmp.$$" "$FINAL"

	{
		echo -e "GWAS\tN_SMALL\tN_LIFTED\tN_UNMAPPED"
		n_small=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$work/std.old.tsv")
		n_lift=$(awk 'END{print NR+0}' "$work/lifted.bed")
		n_unmap=$(grep -vc '^#' "$work/unmapped.bed" 2>/dev/null || echo 0)
		echo -e "$GWAS\t$n_small\t$n_lift\t$n_unmap"
	} >"${QC_PREFIX}.liftover.n.tsv"
}

gwas_clean_make_cis() {
	if [[ -z "${CIS_BED:-}" ]]; then
		echo "ERROR: --cis-bed is required for cis output" >&2
		exit 2
	fi
	if [[ "$REPLACE" != "TRUE" ]] && gwas_clean_gzip_ok "$CIS_OUT"; then
		gwas_clean_log "cis file exists: $CIS_OUT"
		return 0
	fi
	gwas_clean_need_file "$FINAL"
	gwas_clean_need_file "$CIS_BED"
	gwas_clean_log "cis subset: $FINAL -> $CIS_OUT"

	tmp="${CIS_OUT}.tmp.$$"
	gwas_clean_zcat "$FINAL" | awk -v FS='\t' -v OFS='\t' \
		-v cis_bed="$CIS_BED" -v cis_name="$GWAS" -v flank="$CIS_FLANK" '
    function normchr(x){gsub(/^chr/,"",x); if(x=="X") x="23"; if(x=="Y") x="24"; if(x=="MT"||x=="M") x="25"; return x}
    BEGIN{
      while((getline line < cis_bed)>0){
        gsub(/\r/,"",line); if(line=="" || line ~ /^#/) continue
        split(line,a,/[ \t]+/)
        if(a[4]==cis_name){n++; rchr[n]=normchr(a[1]); rstart[n]=a[2]-flank; if(rstart[n]<0) rstart[n]=0; rend[n]=a[3]+flank}
      }
      close(cis_bed)
    }
    NR==1{for(i=1;i<=NF;i++) c[$i]=i; print; next}
    n>0 && NR>1{
      chr=normchr($(c["CHR"])); pos=$(c["POS"])+0; if(chr=="" || pos<=0) next
      ok=0; for(i=1;i<=n;i++){if(chr==rchr[i] && pos>=rstart[i] && pos<=rend[i]){ok=1; break}}
      if(ok) print
    }' | gwas_clean_compress >"$tmp"

	gzip -t "$tmp"
	mv -f "$tmp" "$CIS_OUT"
	gwas_clean_zcat "$CIS_OUT" | awk -v g="$GWAS" 'NR==1{next} END{print g"\t"NR-1}' >"${QC_PREFIX}.cis.nrow.tsv"
}

run_lead() {
	[[ "$DO_STEP" == "all" || "$DO_STEP" == "lead" ]]
}

awk_lead() {
	cat <<'AWK'
function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
NR==1{
  for(i=1;i<=NF;i++) c[$i]=i
  header=$0
  print header
  next
}
NR>1 && ("SNP" in c) && ("CHR" in c) && ("POS" in c) && ("P" in c) {
  p=$(c["P"])
  if(!isnum(p) || p+0 > pthr) next
  chr=$(c["CHR"])
  pos=$(c["POS"])+0
  if(chr=="" || pos<=0) next
  bin=int((pos-1)/win)
  key=chr SUBSEP bin
  if(!(key in best_p) || p+0 < best_p[key]){
    best_p[key]=p+0
    best_line[key]=$0
    seen[key]=1
  }
}
END{
  for(key in seen) print best_line[key]
}
AWK
}

ref_chr() {
	local x="$1" chr
	x=$(basename "$x")
	if [[ "$x" =~ [Cc][Hh][Rr]([0-9]+|X|Y|MT|M)([^0-9A-Za-z]|$) ]]; then
		chr="${BASH_REMATCH[1]}"
		case "$chr" in
			X) echo 23 ;;
			Y) echo 24 ;;
			MT | M) echo 25 ;;
			*) echo "$chr" ;;
		esac
		return 0
	fi
	return 1
}

chr_label() {
	local chr="$1"
	case "$chr" in
		23) echo chrX ;;
		24) echo chrY ;;
		25) echo chrMT ;;
		*) echo "chr$chr" ;;
	esac
}

want_chr() {
	local chr="$1" spec=",${2:-all},"
	spec=${spec,,}
	[[ "$spec" == *,all,* ]] && return 0
	[[ "$spec" == *",$chr,"* ]] && return 0
	[[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ && "$spec" == *",autosome,"* ]] && return 0
	[[ "$chr" == 23 && "$spec" == *",x,"* ]] && return 0
	[[ "$chr" == 24 && "$spec" == *",y,"* ]] && return 0
	return 1
}

ref_clump_pfiles() {
	local ref="$1" f prefix label chr found=FALSE
	if [[ -f "${ref}.pgen" && -f "${ref}.psam" ]] &&
		[[ -f "${ref}.pvar" || -f "${ref}.pvar.zst" ]]; then
		label=$(basename "$ref")
		chr=$(ref_chr "$label") || {
			echo "ERROR: cannot infer chromosome from refGen_clump prefix: $ref" >&2
			return 1
		}
		echo "$ref $label $chr ${ref}.bad.snp"
		return 0
	fi
	if [[ -d "$ref" ]]; then
		while read -r f; do
			prefix=${f%.pgen}
			[[ -f "${prefix}.psam" && (-f "${prefix}.pvar" || -f "${prefix}.pvar.zst") ]] || continue
			label=$(basename "$prefix")
			# Directory discovery uses one full reference per chromosome. Sex-only
			# chrX.male/chrX.female subsets would repeat and overwrite chrX outputs;
			# chrXY is a separate PAR reference, not a second chrX reference.
			[[ "$label" =~ ^[Cc][Hh][Rr]([0-9]+|X|Y|MT|M)$ ]] || continue
			chr=$(ref_chr "$label") || {
				echo "ERROR: cannot infer chromosome from refGen_clump prefix: $prefix" >&2
				return 1
			}
			echo "$prefix $label $chr ${prefix}.bad.snp"
			found=TRUE
		done < <(find "$ref" -maxdepth 1 -type f -name 'chr*.pgen' | sort -V)
		[[ "$found" == "TRUE" ]] && return 0
	fi
	while read -r f; do
		prefix=${f%.pgen}
		[[ -f "${prefix}.psam" && (-f "${prefix}.pvar" || -f "${prefix}.pvar.zst") ]] || continue
		label=$(basename "$prefix")
		[[ "$label" =~ ^[Cc][Hh][Rr]([0-9]+|X|Y|MT|M)$ ]] || continue
		chr=$(ref_chr "$label") || {
			echo "ERROR: cannot infer chromosome from refGen_clump prefix: $prefix" >&2
			return 1
		}
		echo "$prefix $label $chr ${prefix}.bad.snp"
		found=TRUE
	done < <(compgen -G "${ref}chr*.pgen" | sort -V)
	[[ "$found" == "TRUE" ]] && return 0
	echo "ERROR: no pfile reference found from --refGen_clump: $ref" >&2
	return 1
}

cojo_bfile() {
	local ref="$1" chr="$2" tag
	tag=$(chr_label "$chr")
	if [[ -d "$ref" ]]; then
		echo "${ref%/}/$tag"
	else
		echo "${ref}${tag}"
	fi
}

run_tool() {
	local prefix="$1" tmp log rc
	shift
	tmp="${prefix}.tmp.out.$$"
	log="${prefix}.log"
	rm -f "${prefix}.err"
	gwas_clean_log "CMD: $*"
	if "$@" >"$tmp" 2>&1; then
		if [[ -s "$log" ]]; then
			rm -f "$tmp"
		else
			mv -f "$tmp" "$log"
		fi
		return 0
	else
		rc=$?
	fi
	if [[ -s "$log" ]]; then
		rm -f "$tmp"
	else
		mv -f "$tmp" "$log"
	fi
	{
		echo "ERROR: command failed with exit=$rc"
		echo "ERROR: log=$log"
		grep -Ei 'error|failed|invalid' "$log" || true
	} >"${prefix}.err"
	return "$rc"
}

cat_tables() {
	local out="$1" first=TRUE file
	shift
	[[ "$#" -gt 0 ]] || return 1
	: >"$out"
	for file in "$@"; do
		[[ -s "$file" ]] || continue
		if [[ "$first" == "TRUE" ]]; then
			cat "$file" >>"$out"
			first=FALSE
		else
			tail -n +2 "$file" >>"$out"
		fi
	done
	[[ "$first" == "FALSE" ]]
}

concat_chr_outputs() {
	local chr_dir="$1" out_prefix="$2" labels="$3" suffix label file out files=()
	shift 3
	[[ -s "$labels" ]] || return 0

	for suffix in "$@"; do
		files=()
		while read -r label; do
			[[ -n "$label" ]] || continue
			if [[ -d "$chr_dir" ]]; then
				file="${chr_dir}/${label}.${suffix}"
			else
				file="${chr_dir}.${label}.${suffix}"
			fi
			[[ -s "$file" ]] && files+=("$file")
		done <"$labels"

		[[ "${#files[@]}" -gt 0 ]] || continue
		out="${out_prefix}.${suffix}"
		cat_tables "$out" "${files[@]}" || true
	done
}

gwas_clean_run_core() {
	local run_small=FALSE run_liftover=FALSE run_cis=FALSE
	mkdir -p "$(dirname "$SMALL")" "$(dirname "$FINAL")" "$(dirname "$QC_PREFIX")"

	case "$DO_STEP" in
		small) run_small=TRUE ;;
		liftover) run_liftover=TRUE ;;
		cis) run_cis=TRUE ;;
		lead) ;;
		all)
			run_small=TRUE
			[[ "$DO_LIFTOVER" == "TRUE" ]] && run_liftover=TRUE
			[[ -n "${CIS_BED:-}" ]] && run_cis=TRUE
			;;
	esac

	if [[ "$run_small" == "TRUE" ]]; then
		gwas_clean_need_file "$RAW"
		if [[ "$REPLACE" != "TRUE" ]] && gwas_clean_gzip_ok "$SMALL"; then
			gwas_clean_log "Small GWAS exists: $SMALL"
		else
			std_small "$RAW" "$SMALL" "${QC_PREFIX}.header.small.txt"
		fi

		if [[ "$DO_STEP" == "small" || "$DO_STEP" == "all" ]]; then
			if [[ "$DO_LIFTOVER" == "TRUE" ]]; then
				gwas_clean_log "Final clean GWAS will be generated by liftOver."
			fi
		fi

		if [[ "$DELETE_RAW" == "TRUE" && -n "$RAW" && -f "$RAW" ]]; then
			gzip -t "$SMALL"
			gwas_clean_log "delete raw: $RAW"
			rm -f "$RAW"
		fi
	fi

	if [[ "$run_liftover" == "TRUE" ]]; then
		gwas_clean_liftover_small_to_final
	fi

	if [[ ! -s "$FINAL" && -s "$SMALL" && "$DO_LIFTOVER" != "TRUE" ]]; then
		gwas_clean_copy_small_to_final
	fi

	if [[ "$run_cis" == "TRUE" && "${DEFER_CIS:-FALSE}" != "TRUE" ]]; then
		gwas_clean_make_cis
	fi

	gwas_clean_log "DONE"
}

# The generated worker loads definitions before its run-specific functions.
if [[ ${1:-} == --base ]]; then return 0; fi

if [[ ${1:-} == --worker ]]; then


	# 🚩 Worker scoring completion
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
					fi
					;;
				no_reference_matched_variants | no_snps_selected) ;;
				*) reason="" ;;
			esac
		fi
		if [[ -z "$reason" ]]; then
			gwas_post_pgs_with_scores
			return
		fi
		mkdir -p "$PGS_DIR"
		rm -f -- "$PGS_DONE"
		# Match the established no-matched-variants schema; no invented zero scores.
		printf '#IID\t%s.ALLELE_CT\t%s.SCORE_SUM\n' "$GWAS" "$GWAS" | gzip -c >"${PGS_OUTPUT}.tmp.$$"
		mv -f -- "${PGS_OUTPUT}.tmp.$$" "$PGS_OUTPUT"
		{
			printf 'key\tvalue\ngwas\t%s\nsource\t%s\nmatched_variants\tnone\nreason\t%s\n' "$GWAS" "$PGS_SCORE_FILE" "$reason"
			printf 'pfile_dir\t%s\nchromosomes\tnone\noutput\t%s\nchromosome_intermediates\tnone\n' "$PGS_PFILE_DIR" "$PGS_OUTPUT"
		} >"${PGS_META}.tmp.$$"
		mv -f -- "${PGS_META}.tmp.$$" "$PGS_META"
		printf 'GWAS\tSTATUS\tTIME\n%s\tcomplete_no_score_variants\t%s\n' "$GWAS" "$(date '+%F %T')" >"${PGS_DONE}.tmp.$$"
		mv -f -- "${PGS_DONE}.tmp.$$" "$PGS_DONE"
		gwas_post_log "PGS skipped: $reason (header-only output; no individual scores)"
	}


	# 🚩 Worker formatting and post-processing
	# Performance overrides for format_gwas.sh at multi-thousand-GWAS scale.
	# This file is sourced by each generated per-GWAS command after format.f.sh.

	# Preserve the original fill function so integrated format-time N filling can
	# skip only its redundant rewrite while retaining match_EAF behavior.
	if declare -F gwas_post_fill_missing_fields >/dev/null 2>&1 &&
		! declare -F gwas_post_fill_missing_fields_legacy >/dev/null 2>&1; then
		eval "$(declare -f gwas_post_fill_missing_fields | sed '1s/gwas_post_fill_missing_fields/gwas_post_fill_missing_fields_legacy/')"
	fi
	if declare -F gwas_clean_liftover_small_to_final >/dev/null 2>&1 &&
		! declare -F gwas_clean_liftover_small_to_final_legacy >/dev/null 2>&1; then
		eval "$(declare -f gwas_clean_liftover_small_to_final | sed '1s/gwas_clean_liftover_small_to_final/gwas_clean_liftover_small_to_final_legacy/')"
	fi

	: "${GWAS_POST_SORT_MEMORY:=512M}"
	: "${GWAS_CLEAN_COMP_THREADS:=1}"
	GWAS_POST_N_FILLED_DURING_FORMAT=FALSE
	GWAS_POST_VIEWS_READY=FALSE
	GWAS_POST_MAGMA_ROWS=""
	GWAS_POST_MPLOT_INPUT=""
	GWAS_POST_LEAD_VIEW=""

	# Check the indexed input, not just chromosomes that happened to have an LD
	# panel. Missing panels must never produce a whole-GWAS completion marker.
	gwas_post_check_chromosome_coverage() {
		local phase="$1" refs="$2" audit="${QC_PREFIX}.${1}.chromosomes.log"
		local contigs="$GWAS_POST_TMP/${phase}.contigs" label chr row cref ext missing=0 reason
		tabix -l "$FINAL" >"$contigs" || return 1
		printf 'CHR\tSTATUS\tREASON\n' >"$audit"
		while IFS= read -r label; do
			chr=${label#chr}
			chr=${chr^^}
			case "$chr" in X) chr=23 ;; Y) chr=24 ;; M | MT) chr=25 ;; esac
			if [[ "$phase" == lead ]] && ! want_chr "$chr" "$CHRS"; then
				printf '%s\tEXCLUDED\texplicit --chr %s\n' "$label" "$CHRS" >>"$audit"
				continue
			fi
			reason=''
			if [[ "$phase" == lead ]]; then
				row=$(awk -F '\t' -v chr="$chr" '$3==chr{print;exit}' "$refs")
				if [[ -z "$row" ]]; then
					reason="missing clump reference for input chromosome $label"
				else
					cref=$(cut -f6 <<<"$row")
					for ext in bed bim fam; do
						[[ -s "${cref}.${ext}" ]] || reason="missing COJO reference ${cref}.${ext}"
					done
				fi
			else
				awk -v chr="$chr" '$1==chr{ok=1}END{exit !ok}' "$refs" ||
					reason="missing MAGMA LD reference for input chromosome $label"
				# The current MAGMA annotation path supports chromosomes 1-23. Do not
				# silently drop Y/MT even when a custom LD reference contains them.
				if [[ "$chr" == 24 || "$chr" == 25 ]]; then
					reason="chromosome $label is not supported by the current MAGMA annotation pipeline"
				fi
			fi
			if [[ -n "$reason" ]]; then
				printf '%s\tUNAVAILABLE\t%s\n' "$label" "$reason" >>"$audit"
				gwas_post_log "ERROR: $phase: $reason; see $audit"
				missing=1
			else
				printf '%s\tREFERENCE_AVAILABLE\tinput chromosome included\n' "$label" >>"$audit"
			fi
		done <"$contigs"
		((missing == 0))
	}

	gwas_clean_compress() {
		if command -v bgzip >/dev/null 2>&1; then
			bgzip -@ "$GWAS_CLEAN_COMP_THREADS" -c
		elif command -v pigz >/dev/null 2>&1; then
			pigz -c -p "$GWAS_CLEAN_COMP_THREADS"
		else
			gzip -c
		fi
	}

	gwas_clean_gzip_ok() {
		local file="$1"
		[[ -s "$file" ]] || return 1
		if { [[ -s "${file}.tbi" ]] || [[ -s "${file}.csi" ]]; } &&
			command -v tabix >/dev/null 2>&1; then
			tabix -l "$file" >/dev/null 2>&1
		else
			gzip -t "$file" >/dev/null 2>&1
		fi
	}

	gwas_post_sort_bgzf() {
		local output="$1" chr_col="$2" pos_col="$3" tmp
		tmp="${output}.tmp.$$"
		{
			IFS= read -r header
			printf '%s\n' "$header"
			sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' \
				-k"${chr_col}","${chr_col}"n -k"${pos_col}","${pos_col}"n -k1,1
		} | gwas_clean_compress >"$tmp"
		[[ -s "$tmp" ]] || {
			echo "ERROR: empty BGZF output: $output" >&2
			return 1
		}
		mv -f -- "$tmp" "$output"
	}

	# Full-format mode: one raw scan does schema conversion, coordinate filtering,
	# optional N filling, and QC counting. No post-format counting scan is needed.
	std_format() {
		local src="$1" out="$2" header_out="$3" input_fs tmp audit miss_col=0
		gwas_post_need_file "$src"
		gwas_post_log "Format full GWAS (single scan): $src -> $out"
		gwas_clean_header_names "$src" "$header_out"
		SNP_col=$(gwas_clean_col "${SNP_col:-}")
		CHR_col=$(gwas_clean_col "${CHR_col:-}")
		POS_col=$(gwas_clean_col "${POS_col:-}")
		EA_col=$(gwas_clean_col "${EA_col:-}")
		NEA_col=$(gwas_clean_col "${NEA_col:-}")
		EAF_col=$(gwas_clean_col "${EAF_col:-}")
		N_col=$(gwas_clean_col "${N_col:-}")
		BETA_col=$(gwas_clean_col "${BETA_col:-}")
		SE_col=$(gwas_clean_col "${SE_col:-}")
		P_col=$(gwas_clean_col "${P_col:-}")
		LOG10P_col=$(gwas_clean_col "${LOG10P_col:-}")
		[[ "$SNP_col" -gt 0 && "$CHR_col" -gt 0 && "$POS_col" -gt 0 && "$EA_col" -gt 0 &&
			"$NEA_col" -gt 0 && "$BETA_col" -gt 0 && "$SE_col" -gt 0 && "$P_col" -gt 0 ]] || {
			echo "ERROR: required GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: $src" >&2
			return 1
		}
		input_fs=$(gwas_clean_detect_fs "$src")
		if [[ -n "${N_TOTAL:-}" ]]; then
			miss_col=$(gwas_clean_header_line "$src" | awk -v FS="$input_fs" '{for(i=1;i<=NF;i++)if(toupper($i)=="F_MISS"){print i;exit}}')
			[[ "$miss_col" =~ ^[1-9][0-9]*$ ]] || {
				echo 'ERROR: --n-total requires an F_MISS column' >&2
				return 1
			}
		fi
		tmp="${out}.tmp.$$"
		audit="${QC_PREFIX}.format.audit.tsv"
		gwas_post_zcat "$src" | awk -v FS="$input_fs" -v OFS='\t' \
			-v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
			-v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
			-v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" \
			-v fill_n="$FILL_N" -v n_total="${N_TOTAL:-}" -v miss_col="$miss_col" -v audit="$audit" -v gwas="$GWAS" '
    function get(c, x){x=(c>0 ? $c : "");gsub(/^[[:space:]]+|[[:space:]]+$/,"",x);return x}
    function val(c, x){x=get(c);return x=="" ? "NA" : x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function validn(x){return x ~ /^[0-9]+([.][0-9]+)?$/ && x+0>0}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P";next}
    {
      chr=normchr(get(chr_col));pos=get(pos_col)
      if(chr!~/^[0-9]+$/ || chr+0<1 || chr+0>25 || pos!~/^[0-9]+$/ || pos+0<1){badcoord++;next}
      snp=get(snp_col);ea=get(ea_col);nea=get(nea_col)
      if(snp==""||snp=="NA"||snp=="."){snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
      p=(p_col>0 ? get(p_col) : "");lp=(logp_col>0 ? get(logp_col) : "")
      if(p=="" && isnum(lp))p=10^(-lp);if(lp=="" && isnum(p) && p>0)lp=-log(p)/log(10)
      if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
      n=val(n_col)
      if(n_total!="" && !validn(n)){
        miss=get(miss_col)
        if(!isnum(miss)||miss+0<0||miss+0>=1){print "ERROR: invalid F_MISS at row " NR > "/dev/stderr";exit 2}
        n=int(n_total*(1-miss)+0.5);nfilled++
      }else if(fill_n!="" && !validn(n)){n=fill_n;nfilled++}else if(validn(n))nkept++
      print snp,chr+0,pos+0,ea,nea,val(eaf_col),n,val(beta_col),val(se_col),p,lp
      kept++
    }
    END{
      print "GWAS\tN_OUTPUT\tN_DROPPED_BAD_COORD\tN_EXISTING\tN_FILLED" > audit
      print gwas "\t" kept+0 "\t" badcoord+0 "\t" nkept+0 "\t" nfilled+0 >> audit
    }' | {
			IFS= read -r format_header
			printf '%s\n' "$format_header"
			sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k2,2n -k3,3n -k1,1
		} | gwas_clean_compress >"$tmp"
		[[ -s "$tmp" && -s "$audit" ]] || {
			echo "ERROR: format failed: $src" >&2
			rm -f -- "$tmp"
			return 1
		}
		mv -f -- "$tmp" "$out"
		rm -f -- "${out}.tbi" "${out}.csi"
		awk -F '\t' 'NR==2{print $1 "\t" $2}' "$audit" >"${QC_PREFIX}.format.nrow.tsv"
		[[ -z "$FILL_N" ]] || GWAS_POST_N_FILLED_DURING_FORMAT=TRUE
	}

	# Skip auxiliary panels such as chrXY, chrX.female, and chrX.male while
	# enumerating a reference directory.  lead/MAGMA use only the canonical
	# chr1..chr22/chrX/chrY panels; otherwise several prefixes can collapse to the
	# same normalized chromosome and overwrite the same per-chromosome files.
	ref_clump_pfiles() {
		local ref="$1" f prefix label chr found=FALSE
		if [[ -f "${ref}.pgen" && -f "${ref}.psam" && (-f "${ref}.pvar" || -f "${ref}.pvar.zst") ]]; then
			label=$(basename "$ref")
			chr=$(ref_chr "$label") || return 1
			printf '%s %s %s\n' "$ref" "$label" "$chr"
			return 0
		fi
		if [[ -d "$ref" ]]; then
			while IFS= read -r f; do
				prefix=${f%.pgen}
				[[ -f "${prefix}.psam" && (-f "${prefix}.pvar" || -f "${prefix}.pvar.zst") ]] || continue
				label=$(basename "$prefix")
				chr=$(ref_chr "$label") || continue
				[[ "$label" == "$(chr_label "$chr")" ]] || continue
				printf '%s %s %s\n' "$prefix" "$label" "$chr"
				found=TRUE
			done < <(find "$ref" -maxdepth 1 -type f -name 'chr*.pgen' | sort -V)
		else
			while IFS= read -r f; do
				prefix=${f%.pgen}
				[[ -f "${prefix}.psam" && (-f "${prefix}.pvar" || -f "${prefix}.pvar.zst") ]] || continue
				label=$(basename "$prefix")
				chr=$(ref_chr "$label") || continue
				[[ "$label" == "$(chr_label "$chr")" ]] || continue
				printf '%s %s %s\n' "$prefix" "$label" "$chr"
				found=TRUE
			done < <(compgen -G "${ref}chr*.pgen" | sort -V)
		fi
		[[ "$found" == TRUE ]] || {
			echo "ERROR: no chromosome pfile reference found: $ref" >&2
			return 1
		}
	}

	# --hm3 TRUE: filter and standardize in one raw scan. Keep exactly the union
	# of HapMap3 IDs/coordinates and variants meeting --p-hm3, then deduplicate
	# exact rows.  Coordinates are the fallback for files whose rsID column is NA.
	std_hm3() {
		local src="$1" out="$2" header_out="$3" input_fs tmp_tsv tmp audit n_out
		gwas_clean_need_file "$src"
		gwas_clean_need_file "$HM3"
		[[ -z "${HM3_POS:-}" ]] || gwas_clean_need_file "$HM3_POS"
		gwas_post_log "HM3 GWAS (single scan): $src -> $out"
		gwas_clean_header_names "$src" "$header_out"
		SNP_col=$(gwas_clean_col "${SNP_col:-}")
		CHR_col=$(gwas_clean_col "${CHR_col:-}")
		POS_col=$(gwas_clean_col "${POS_col:-}")
		EA_col=$(gwas_clean_col "${EA_col:-}")
		NEA_col=$(gwas_clean_col "${NEA_col:-}")
		EAF_col=$(gwas_clean_col "${EAF_col:-}")
		N_col=$(gwas_clean_col "${N_col:-}")
		BETA_col=$(gwas_clean_col "${BETA_col:-}")
		SE_col=$(gwas_clean_col "${SE_col:-}")
		P_col=$(gwas_clean_col "${P_col:-}")
		LOG10P_col=$(gwas_clean_col "${LOG10P_col:-}")
		[[ "$CHR_col" -gt 0 && "$POS_col" -gt 0 && ("$P_col" -gt 0 || "$LOG10P_col" -gt 0) ]] || {
			echo "ERROR: CHR/POS and P/LOG10P are required: $src" >&2
			return 1
		}
		input_fs=$(gwas_clean_detect_fs "$src")
		tmp_tsv="$GWAS_POST_TMP/${GWAS}.hm3.sorted.tsv"
		tmp="${out}.tmp.$$"
		audit="${QC_PREFIX}.hm3.audit.tsv"
		gwas_clean_zcat "$src" | awk -v FS="$input_fs" -v OFS='\t' \
			-v hm3_file="$HM3" -v hm3_pos_file="${HM3_POS:-}" -v pthr="$P_HM3" \
			-v snp_col="$SNP_col" -v chr_col="$CHR_col" -v pos_col="$POS_col" \
			-v ea_col="$EA_col" -v nea_col="$NEA_col" -v eaf_col="$EAF_col" -v n_col="$N_col" \
			-v beta_col="$BETA_col" -v se_col="$SE_col" -v p_col="$P_col" -v logp_col="$LOG10P_col" \
			-v fill_n="$FILL_N" -v audit="$audit" -v gwas="$GWAS" '
    function get(c,x){x=(c>0?$c:"");gsub(/^[[:space:]]+|[[:space:]]+$/,"",x);return x}
    function val(c,x){x=get(c);return x==""?"NA":x}
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function validn(x){return x~/^[0-9]+([.][0-9]+)?$/&&x+0>0}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    BEGIN{
      while((getline x<hm3_file)>0){gsub(/\r/,"",x);split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)
      if(hm3_pos_file!="")while((getline x<hm3_pos_file)>0){
        gsub(/\r/,"",x);split(x,a,/[ \t]+/);ch=normchr(a[1]);bp=a[4]
        key=ch SUBSEP (bp+0)
        if(ch~/^[0-9]+$/&&bp~/^[0-9]+$/&&bp+0>0&&!(key in hm3pos))hm3pos[key]=a[2]
      }
      if(hm3_pos_file!="")close(hm3_pos_file)
    }
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P";next}
    {
      chr=normchr(get(chr_col));pos=get(pos_col);if(chr!~/^[0-9]+$/||chr+0<1||chr+0>25||pos!~/^[0-9]+$/||pos+0<1){badcoord++;next}
      key=chr SUBSEP (pos+0);in_hm3_pos=(key in hm3pos)
      snp=get(snp_col);ea=get(ea_col);nea=get(nea_col);if(snp==""||snp=="NA"||snp=="."){
        if(in_hm3_pos&&hm3pos[key]!="")snp=hm3pos[key]
        else{snp=chr ":" pos;if(ea!="")snp=snp ":" ea;if(nea!="")snp=snp ":" nea}
      }
      p=(p_col>0?get(p_col):"");lp=(logp_col>0?get(logp_col):"");if(p==""&&isnum(lp))p=10^(-lp);if(lp==""&&isnum(p)&&p>0)lp=-log(p)/log(10)
      keep=(snp in hm3)||in_hm3_pos||(isnum(p)&&p+0<pthr+0);if(!keep)next
      if(ea=="")ea="NA";if(nea=="")nea="NA";if(p=="")p="NA";if(lp=="")lp="NA"
      n=val(n_col);if(fill_n!=""&&!validn(n)){n=fill_n;nfilled++}else if(validn(n))nkept++
      print snp,chr+0,pos+0,ea,nea,val(eaf_col),n,val(beta_col),val(se_col),p,lp;kept++
    }
    END{print "GWAS\tN_OUTPUT\tN_DROPPED_BAD_COORD\tN_EXISTING\tN_FILLED">audit;print gwas "\t" kept+0 "\t" badcoord+0 "\t" nkept+0 "\t" nfilled+0>>audit}
  ' | {
			IFS= read -r header
			printf '%s\n' "$header"
			sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' \
				-k2,2n -k3,3n -k1,1 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 -k11,11 -u
		} >"$tmp_tsv"
		n_out=$(($(wc -l <"$tmp_tsv") - 1))
		((n_out < 0)) && n_out=0
		gwas_clean_compress <"$tmp_tsv" >"$tmp"
		[[ -s "$tmp" ]] || {
			echo "ERROR: HM3 format failed: $src" >&2
			return 1
		}
		mv -f -- "$tmp" "$out"
		rm -f -- "${out}.tbi" "${out}.csi"
		printf '%s\t%s\n' "$GWAS" "$n_out" >"${QC_PREFIX}.hm3.nrow.tsv"
		[[ -z "$FILL_N" ]] || GWAS_POST_N_FILLED_DURING_FORMAT=TRUE
	}

	# Adapter for the unchanged format.f.sh formatter API.
	std_small() {
		if [[ "${HM3_MODE:-FALSE}" == TRUE ]]; then
			std_hm3 "$1" "$2" "${QC_PREFIX}.header.hm3.txt"
		else
			std_format "$1" "$2" "${QC_PREFIX}.header.format.txt"
		fi
	}

	gwas_post_fill_missing_fields() {
		local target="$1" saved_fill_n="$FILL_N"
		if [[ "$GWAS_POST_N_FILLED_DURING_FORMAT" == TRUE ]]; then FILL_N=""; fi
		gwas_post_fill_missing_fields_legacy "$target"
		FILL_N="$saved_fill_n"
		if [[ -n "$saved_fill_n" || "$FILL_EAF" == TRUE ]]; then
			rm -f -- "${target}.tbi" "${target}.csi"
		fi
	}

	gwas_clean_liftover_small_to_final() {
		python3 "$LIFTOVER_HELPER" liftover --input "$SMALL" --output "$FINAL" \
			--chain "$CHAIN" --liftOver "$LIFTOVER_BIN" --qc-prefix "$QC_PREFIX" \
			--source-build 37 --target-build 38 --replace "$REPLACE"
	}

	gwas_post_ensure_index() {
		local file="$1"
		[[ -s "$file" ]] || return 0
		if ! gwas_index_valid "$file"; then gwas_index_file "$file"; fi
	}

	gwas_post_cis_regions() {
		local output="$1"
		awk -v name="$GWAS" -v flank="$CIS_FLANK" 'BEGIN{OFS="\t"}
    function chr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    /^[[:space:]]*($|#)/{next}
    $4==name{s=$2-flank;if(s<1)s=1;e=$3+flank;print chr($1),s,e}
  ' "$CIS_BED" | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN{OFS="\t"}
    NR==1{c=$1;s=$2;e=$3;next}
    $1==c&&$2<=e+1{if($3>e)e=$3;next}
    {print c,s,e;c=$1;s=$2;e=$3}
    END{if(NR)print c,s,e}' >"$output"
	}

	# Indexed cis extraction is proportional to the requested loci, not to the
	# whole GWAS. A full-scan fallback remains for legacy/non-indexable inputs.
	gwas_clean_make_cis() {
		local regions="$GWAS_POST_TMP/${GWAS}.cis.regions.tsv" body="$GWAS_POST_TMP/${GWAS}.cis.body.tsv" tmp header
		[[ -n "${CIS_BED:-}" ]] || {
			echo "ERROR: --cis-bed is required for cis output" >&2
			return 2
		}
		if [[ "$REPLACE" != TRUE && -s "$CIS_OUT" ]]; then
			if gwas_index_valid "$CIS_OUT"; then
				gwas_clean_log "indexed cis file exists: $CIS_OUT"
				return 0
			fi
			gwas_clean_log "repair missing/stale cis BGZF index: $CIS_OUT"
			if gwas_index_file "$CIS_OUT" && gwas_index_valid "$CIS_OUT"; then
				gwas_clean_log "indexed cis file ready: $CIS_OUT"
				return 0
			fi
			echo "ERROR: failed to repair cis BGZF/index: $CIS_OUT" >&2
			return 1
		fi
		gwas_clean_need_file "$FINAL"
		gwas_clean_need_file "$CIS_BED"
		gwas_post_cis_regions "$regions"
		tmp="${CIS_OUT}.tmp.$$"
		header=$(gwas_index_header "$FINAL")
		[[ -n "$header" ]] || {
			echo "ERROR: unreadable GWAS header: $FINAL" >&2
			return 1
		}
		: >"$body"
		if [[ -s "$regions" ]]; then
			gwas_post_ensure_index "$FINAL"
			while IFS=$'\t' read -r chr start end; do
				tabix "$FINAL" "${chr}:${start}-${end}"
			done <"$regions" >"$body"
		fi
		{
			printf '%s\n' "$header"
			sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k2,2n -k3,3n -k1,1 -u "$body"
		} |
			gwas_clean_compress >"$tmp"
		mv -f -- "$tmp" "$CIS_OUT"
		rm -f -- "${CIS_OUT}.tbi" "${CIS_OUT}.csi"
		gwas_index_file "$CIS_OUT"
		printf '%s\t%s\n' "$GWAS" "$(wc -l <"$body" | tr -d ' ')" >"${QC_PREFIX}.cis.nrow.tsv"
		gwas_clean_log "cis subset via tabix: $FINAL -> $CIS_OUT"
	}

	gwas_post_thin() {
		[[ "${THIN_MODE:-FALSE}" == TRUE ]] || return 0
		[[ "$DO_STEP" == all || ",$DO_STEP," =~ ,(format|thin|mplot|liftover), ]] || return 0
		local signature marker="${THIN_OUT}.done" before after
		if gwas_thin_plot_only "$DO_STEP"; then
			if ! gwas_thin_plot_ready "$FINAL" "$THIN_OUT"; then
				echo "ERROR: mplot requires an existing thin GWAS with a valid, current index: $THIN_OUT" >&2
				echo "ERROR: run thin first (or thin,mplot); mplot does not regenerate thin data." >&2
				return 1
			fi
			gwas_post_log "Reuse existing thin GWAS for plotting: $THIN_OUT"
			return 0
		fi
		if [[ "$REPLACE" != TRUE ]] && gwas_thin_complete "$FINAL" "$THIN_OUT" "$GRCH" "$THIN_CHR_MAX" "$HM3" "${HM3_POS:-}" "$THIN_R" "$PHE_R"; then
			gwas_post_log "SKIP completed thin GWAS: $THIN_OUT"
			return 0
		fi
		# Invalidate legacy plot completion before regenerating the thin data.
		rm -f -- "$marker" "${MH_META}"
		before=$(stat -Lc '%s|%y' -- "$FINAL") || return 1
		Rscript --vanilla "$THIN_R" thin --input "$FINAL" --output "$THIN_OUT" \
			--grch "$GRCH" --chr-max "$THIN_CHR_MAX" --hm3-file "$HM3" --hm3-pos "${HM3_POS:-}" \
			--phe-r "$PHE_R" --replace "$REPLACE" || return $?
		after=$(stat -Lc '%s|%y' -- "$FINAL") || return 1
		[[ "$before" == "$after" ]] || {
			echo "ERROR: source changed while thinning: $FINAL" >&2
			return 1
		}
		signature=$(gwas_thin_signature "$FINAL" "$THIN_OUT" "$GRCH" "$THIN_CHR_MAX" "$HM3" "${HM3_POS:-}" "$THIN_R" "$PHE_R") || return 1
		printf '%s\n' "$signature" >"$marker.tmp.$$"
		mv -f -- "$marker.tmp.$$" "$marker"
	}

	gwas_post_prepare_views() {
		local want_magma=FALSE want_mplot=FALSE want_lead=FALSE raw_plot plot_tsv
		[[ "$GWAS_POST_VIEWS_READY" == TRUE ]] && return 0
		[[ ",${DO_STEP}," == *,magma,* ]] && want_magma=TRUE
		[[ "$DO_STEP" == all || ",${DO_STEP}," == *,mplot,* ]] && want_mplot=TRUE
		[[ "$DO_STEP" == all || ",${DO_STEP}," == *,lead,* ]] && want_lead=TRUE
		if [[ "$want_mplot" == TRUE && "${THIN_MODE:-FALSE}" == TRUE ]]; then
			gwas_post_need_file "$THIN_OUT"
			GWAS_POST_MPLOT_INPUT="$THIN_OUT"
			want_mplot=FALSE
		fi
		[[ "$want_magma" == TRUE || "$want_mplot" == TRUE || "$want_lead" == TRUE ]] || {
			GWAS_POST_VIEWS_READY=TRUE
			return 0
		}
		gwas_post_need_file "$FINAL"
		if [[ "$want_mplot" == TRUE ]]; then
			gwas_post_need_file "$HM3"
			[[ -z "${HM3_POS:-}" ]] || gwas_post_need_file "$HM3_POS"
		fi
		GWAS_POST_MAGMA_ROWS="$GWAS_POST_TMP/${GWAS}.magma.rows.tsv"
		GWAS_POST_LEAD_VIEW="$GWAS_POST_TMP/${GWAS}.lead.tsv"
		raw_plot="$GWAS_POST_TMP/${GWAS}.mplot.tsv"
		plot_tsv="${QC_PREFIX}.mplot.hm3.gz"
		: >"$GWAS_POST_MAGMA_ROWS"
		: >"$GWAS_POST_LEAD_VIEW"
		: >"$raw_plot"
		gwas_post_log "single full-GWAS scan for requested module views: magma=$want_magma mplot=$want_mplot lead=$want_lead"
		gwas_post_zcat "$FINAL" | awk -v FS='\t' -v OFS='\t' \
			-v want_magma="$want_magma" -v want_mplot="$want_mplot" -v want_lead="$want_lead" \
			-v magma_out="$GWAS_POST_MAGMA_ROWS" -v plot_out="$raw_plot" -v lead_out="$GWAS_POST_LEAD_VIEW" \
			-v hm3_file="$HM3" -v hm3_pos_file="${HM3_POS:-}" -v plot_p='0.001' -v lead_p="$P_LEAD" '
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x}
    BEGIN{
      if(want_mplot=="TRUE"){
        while((getline x<hm3_file)>0){split(x,a,/[ \t]+/);if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1}close(hm3_file)
        if(hm3_pos_file!="")while((getline x<hm3_pos_file)>0){
          gsub(/\r/,"",x);split(x,a,/[ \t]+/);ch=normchr(a[1]);bp=a[4]
          if(ch~/^[0-9]+$/&&bp~/^[0-9]+$/&&bp+0>0)hm3pos[ch SUBSEP (bp+0)]=1
        }
        if(hm3_pos_file!="")close(hm3_pos_file)
      }
    }
    NR==1{
      for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}
      if(want_magma=="TRUE"&&(!("SNP" in c)||!("P" in c)||!("CHR" in c)||!("POS" in c)||!("EA" in c)||!("NEA" in c))){print "ERROR: MAGMA view requires SNP/P/CHR/POS/EA/NEA" > "/dev/stderr";exit 2}
      if(want_mplot=="TRUE")print > plot_out
      if(want_lead=="TRUE")print > lead_out
      next
    }
    {
      snp=("SNP" in c?$(c["SNP"]):"");p=("P" in c?$(c["P"]):"");n=("N" in c?$(c["N"]):"NA")
      chr=("CHR" in c?normchr($(c["CHR"])):"");pos=("POS" in c?$(c["POS"]):"")
      if(want_magma=="TRUE"&&snp!=""&&snp!="NA"&&snp!="."&&isnum(p)&&p+0>0&&p+0<=1)print snp,p,n,chr,pos,$(c["EA"]),$(c["NEA"]) > magma_out
      if(want_mplot=="TRUE"&&((snp in hm3)||((chr SUBSEP (pos+0)) in hm3pos)||(isnum(p)&&p+0<plot_p)))print > plot_out
      if(want_lead=="TRUE"&&isnum(p)&&p+0<=lead_p)print > lead_out
    }'
		if [[ "$want_mplot" == TRUE ]]; then
			gwas_clean_compress <"$raw_plot" >"$plot_tsv"
			GWAS_POST_MPLOT_INPUT="$plot_tsv"
		fi
		GWAS_POST_VIEWS_READY=TRUE
	}

	gwas_post_magma_annotation() {
		local ref_tag window_tag cache_key cache_dir annot meta lock_file lock_fd tmp_prefix
		local snploc_raw snploc pvar prefix nloc hash meta_tmp
		ref_tag=$(printf '%s' "$REFGEN_CLUMP" | sha256sum | awk '{print substr($1,1,16)}')
		window_tag=$(printf '%s' "$MAGMA_WINDOW" | tr -c 'A-Za-z0-9._-' '_')
		cache_key="v2.GRCh${GRCH}.window_${window_tag}.ref_${ref_tag}"
		cache_dir="$MAGMA_ANNOT_CACHE/v2/GRCh$GRCH/window_$window_tag/ref_$ref_tag"
		annot="$cache_dir/genes.annot"
		meta="$cache_dir/annotation.meta.tsv"
		mkdir -p "$(dirname "$cache_dir")"
		lock_file="${cache_dir}.lock"
		exec {lock_fd}>"$lock_file"
		flock "$lock_fd"
		if [[ ! -s "$annot" || ! -s "$meta" ]]; then
			mkdir -p "$cache_dir"
			snploc_raw="$GWAS_POST_TMP/reference.snp.loc.raw"
			snploc="$GWAS_POST_TMP/reference.snp.loc"
			: >"$snploc_raw"
			while read -r prefix _rest; do
				if [[ ! -s "${prefix}.pvar" && -s "${prefix}.pvar.zst" ]]; then pvar="${prefix}.pvar.zst"; else pvar="${prefix}.pvar"; fi
				gwas_clean_zcat "$pvar" | awk -v FS='[ \t]+' -v OFS='\t' '
        /^##/{next}
        /^#CHROM/{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
        {id=$(c["ID"]);ch=$(c["CHROM"]);pos=$(c["POS"]);gsub(/^chr/,"",ch);if(id!=""&&id!="."&&ch~/^([1-9]|1[0-9]|2[0-5])$/&&pos~/^[0-9]+$/)print id,ch,pos}' >>"$snploc_raw"
			done < <(ref_clump_pfiles "$REFGEN_CLUMP")
			sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2n -k3,3n "$snploc_raw" | awk '$1!=last{print;last=$1}' >"$snploc"
			nloc=$(wc -l <"$snploc" | tr -d ' ')
			((nloc > 1000)) || {
				echo "ERROR: too few shared MAGMA reference SNPs: $nloc" >&2
				return 1
			}
			hash=$(sha256sum "$snploc" | awk '{print $1}')
			tmp_prefix="$GWAS_POST_TMP/annotation.$$"
			gwas_post_log "Build shared MAGMA annotation from GRCh$GRCH reference coordinates: $annot"
			if [[ "$MAGMA_WINDOW" == "0,0" || "$MAGMA_WINDOW" == 0 ]]; then
				magma --annotate --snp-loc "$snploc" --gene-loc "$GENE_LOC" --out "$tmp_prefix"
			else
				magma --annotate window="$MAGMA_WINDOW" --snp-loc "$snploc" --gene-loc "$GENE_LOC" --out "$tmp_prefix"
			fi
			gwas_post_need_file "$tmp_prefix.genes.annot"
			mv -f -- "$tmp_prefix.genes.annot" "${annot}.tmp.$$"
			mv -f -- "${annot}.tmp.$$" "$annot"
			meta_tmp="${meta}.tmp.$$"
			printf 'key\tvalue\ngrch\t%s\nwindow_kb\t%s\ngene_loc\t%s\nrefgen_clump\t%s\nsnp_loc_n\t%s\nsnp_loc_sha256\t%s\n' \
				"$GRCH" "$MAGMA_WINDOW" "$GENE_LOC" "$REFGEN_CLUMP" "$nloc" "$hash" >"$meta_tmp"
			mv -f -- "$meta_tmp" "$meta"
		else
			gwas_post_log "Reuse shared MAGMA reference annotation: $annot"
		fi
		MAGMA_ANNOT="$annot"
		MAGMA_ANNOT_KEY="$cache_key"
		MAGMA_SNPLOC_N=$(awk -F '\t' '$1=="snp_loc_n"{print $2;exit}' "$meta")
		MAGMA_SNPLOC_HASH=$(awk -F '\t' '$1=="snp_loc_sha256"{print $2;exit}' "$meta")
		flock -u "$lock_fd"
		exec {lock_fd}>&-
	}

	# Scan a shared reference once, using large reads across Windows mounts.
	gwas_post_magma_reference_chromosomes() {
		local bim="$1" out="$2" key cache lock_fd
		key=$({
			printf '%s\n' "$bim"
			stat -Lc '%s %y' -- "$bim"
		} | sha256sum | cut -d ' ' -f 1)
		cache="${GWAS_POST_TMP_BASE:-${TMPDIR:-/tmp}}/shared/magma/${key}.chromosomes"
		mkdir -p "${cache%/*}"
		exec {lock_fd}>"${cache}.lock"
		flock "$lock_fd"
		if [[ ! -s "$cache" ]]; then
			if ! cat -- "$bim" | awk '{ch=$1;sub(/^chr/,"",ch);if(ch=="X")ch=23;if(ch=="Y")ch=24;if(ch=="MT"||ch=="M")ch=25;if(!seen[ch]++)print ch}' >"${cache}.tmp.$$"; then
				rm -f -- "${cache}.tmp.$$"
				flock -u "$lock_fd"
				exec {lock_fd}>&-
				return 1
			fi
			mv -f -- "${cache}.tmp.$$" "$cache"
		fi
		cp -- "$cache" "$out"
		flock -u "$lock_fd"
		exec {lock_fd}>&-
	}

	gwas_post_magma() {
		local pval header narg has_usable_n nloc npval meta_tmp
		[[ ",${DO_STEP}," == *,magma,* ]] || return 0
		local reference_chromosomes="$GWAS_POST_TMP/magma.reference.chromosomes"
		gwas_post_magma_reference_chromosomes "${MAGMA_REF}.bim" "$reference_chromosomes" || return 1
		if ! gwas_post_check_chromosome_coverage magma "$reference_chromosomes"; then
			rm -f "$MAGMA_DIR/magma.done"
			return 1
		fi
		if [[ "$REPLACE" != TRUE && -s "$MAGMA_DIR/magma.done" && -s "$MAGMA_PREFIX.genes.out" && -s "$MAGMA_PREFIX.genes.raw" ]]; then
			gwas_post_prune_magma_dir
			gwas_post_log "MAGMA exists: $MAGMA_PREFIX.genes.out"
			return 0
		fi
		command -v magma >/dev/null 2>&1 || {
			echo "ERROR: magma not found in PATH" >&2
			return 1
		}
		gwas_post_prepare_views
		mkdir -p "$MAGMA_DIR"
		rm -f "$MAGMA_DIR/magma.done" "$MAGMA_DIR/magma.meta.tsv"
		local mapped_rows="$GWAS_POST_TMP/$GWAS.magma.rsids.tsv"
		local mapped_loc="$GWAS_POST_TMP/$GWAS.magma.rsids.loc"
		gwas_post_log "Prepare compatible MAGMA rsIDs and GRCh$GRCH annotation coordinates"
		python3 "${PERF_F%/*}/format.py" magma-ids \
			--rows "$GWAS_POST_MAGMA_ROWS" --output "$mapped_rows" --snploc "$mapped_loc.raw" \
			--audit "${QC_PREFIX}.magma.ids.tsv" --bim "${MAGMA_REF}.bim" \
			--dbsnp "/mnt/f/annot/dbsnp/rsids-v154-hg${GRCH/37/19}.tsv.gz" \
			--cache "$MAGMA_ANNOT_CACHE/rsid-map-v1/GRCh$GRCH" || return 1
		sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2n -k3,3n -u "$mapped_loc.raw" >"$mapped_loc"
		GWAS_POST_MAGMA_ROWS="$mapped_rows"
		[[ -z "$MAGMA_N" || "$MAGMA_N" =~ ^[1-9][0-9]*$ ]] || {
			echo "ERROR: invalid MAGMA N: $MAGMA_N" >&2
			return 1
		}
		if [[ -z "$MAGMA_N" && -s "$GWAS_DIR/$GWAS.magma.N" ]]; then MAGMA_N=$(awk 'NF{print $1;exit}' "$GWAS_DIR/$GWAS.magma.N"); fi
		pval="$GWAS_POST_TMP/$GWAS.pval"
		if [[ -n "$MAGMA_N" ]]; then
			{
				printf 'SNP\tP\n'
				sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last{print $1 "\t" $2;last=$1}'
			} >"$pval"
			narg="N=$MAGMA_N"
		else
			has_usable_n=$(awk -F '\t' '$3~/^[0-9]+([.][0-9]+)?$/&&$3+0>=50{n++}END{print n+0}' "$GWAS_POST_MAGMA_ROWS")
			if ((has_usable_n > 1000)); then
				{
					printf 'SNP\tP\tN\n'
					sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last&&$3~/^[0-9]+([.][0-9]+)?$/&&$3+0>=50{print;last=$1}'
				} >"$pval"
				narg='ncol=N'
			else
				MAGMA_N="$GWAS_N"
				{
					printf 'SNP\tP\n'
					sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -k1,1 -k2,2g "$GWAS_POST_MAGMA_ROWS" | awk -F '\t' '$1!=last{print $1 "\t" $2;last=$1}'
				} >"$pval"
				narg="N=$MAGMA_N"
			fi
		fi
		npval=$(($(wc -l <"$pval") - 1))
		((npval > 1000)) || {
			echo "ERROR: too few MAGMA SNPs: pval=$npval" >&2
			return 1
		}
		gwas_post_magma_annotation "$mapped_loc"
		nloc="$MAGMA_SNPLOC_N"
		magma --bfile "$MAGMA_REF" synonyms="$SYNONYMS" --pval "$pval" "$narg" --gene-annot "$MAGMA_ANNOT" --out "$MAGMA_PREFIX"
		[[ -s "$MAGMA_PREFIX.genes.out" && -s "$MAGMA_PREFIX.genes.raw" ]] || {
			echo "ERROR: MAGMA output missing: $MAGMA_PREFIX" >&2
			return 1
		}
		meta_tmp="$MAGMA_DIR/magma.meta.tsv.tmp.$$"
		printf 'key\tvalue\ngwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\nannotation_cache\t%s\nannotation_key\t%s\nsnp_loc_sha256\t%s\n' \
			"$GWAS" "$GRCH" "$GENE_LOC" "$MAGMA_REF" "$SYNONYMS" "$MAGMA_WINDOW" "$nloc" "$npval" "$MAGMA_ANNOT" "$MAGMA_ANNOT_KEY" "$MAGMA_SNPLOC_HASH" >"$meta_tmp"
		mv -f -- "$meta_tmp" "$MAGMA_DIR/magma.meta.tsv"
		gwas_post_prune_magma_dir
		date '+%F %T' >"$MAGMA_DIR/magma.done"
		gwas_post_log "MAGMA done: $MAGMA_PREFIX.genes.out"
	}

	gwas_post_mplot_current() {
		local genes="${MAGMA_PREFIX}.genes.out" signal_input=none
		[[ -z "$ADD_SIGNAL" ]] || signal_input="$ADD_SIGNAL"
		[[ -s "$MH_PNG" && -s "$MH_META" && -s "$MH_FLAG" && "$MH_PNG" -nt "$FINAL" ]] || return 1
		if [[ "${THIN_MODE:-FALSE}" == TRUE ]]; then
			[[ -s "$THIN_OUT" && "$MH_PNG" -nt "$THIN_OUT" ]] || return 1
		fi
		awk -F '\t' -v thin="${THIN_MODE:-FALSE}" -v cap="${THIN_CHR_MAX:-10000}" -v method="$PLOT_METHOD" -v panel="$ADD_PANEL" -v grch="$GRCH" \
			-v signal="$signal_input" -v match_col="$SIGNAL_MATCH_COL" -v match_value="$SIGNAL_MATCH_VALUE" \
			-v locus_pos="$SIGNAL_LOCUS_POS" -v display_col="$SIGNAL_DISPLAY_COL" -v write_sig="$WRITE_SIG" \
			-v plot_width="$PLOT_WIDTH" -v plot_height="$PLOT_HEIGHT" -v plot_res="$PLOT_RES" '
    $1=="plot_method"&&$2==method{m=1}
    $1=="add_panel"&&$2==panel{p=1}
    $1=="grch"&&$2==grch{g=1}
    $1=="magma_threshold"&&$2+0==2.5e-6{t=1}
    $1=="mplot_style"&&$2==9{s=1}
    $1=="thin_mode"&&$2==thin{tn=1}
    $1=="thin_chr_max"&&$2==cap{tc=1}
    $1=="plot_width"&&$2+0==plot_width+0{x=1}
    $1=="plot_height"&&$2==plot_height{y=1}
    $1=="plot_res"&&$2+0==plot_res+0{r=1}
    $1=="write_sig"&&$2==write_sig{w=1}
    $1=="signal_input"&&$2==signal{a=1}
    $1=="signal_match_col"&&$2==match_col{b=1}
    $1=="signal_match_value"&&$2==match_value{c=1}
    $1=="signal_locus_pos"&&$2==locus_pos{d=1}
    $1=="signal_display_col"&&$2==display_col{e=1}
    END{exit !(m&&p&&g&&t&&s&&x&&y&&r&&w&&a&&b&&c&&d&&e&&tn&&tc)}' "$MH_META" || return 1
		if [[ "$ADD_PANEL" == magma ]]; then
			[[ -s "$genes" && "$MH_PNG" -nt "$genes" ]] || return 1
		fi
		if [[ -n "$ADD_SIGNAL" ]]; then
			[[ -s "$ADD_SIGNAL" && "$MH_PNG" -nt "$ADD_SIGNAL" ]] || return 1
		fi
		if [[ "$WRITE_SIG" == TRUE ]]; then
			[[ -s "$MH_SIG" && "$MH_SIG" -nt "$FINAL" ]] || return 1
			[[ ! -s "$COJO_FILE" || "$MH_SIG" -nt "$COJO_FILE" ]] || return 1
			[[ "$ADD_PANEL" != magma || "$MH_SIG" -nt "$genes" ]] || return 1
		fi
	}

	gwas_post_append_mplot_flag() {
		local lock_file lock_fd fragment_header aggregate_header
		[[ -s "$MH_FLAG" ]] || {
			echo "ERROR: mplot flag fragment is missing: $MH_FLAG" >&2
			return 1
		}
		mkdir -p "$(dirname "$MH_FLAG")" "$(dirname "$MPLOT_FLAG_FILE")"
		command -v flock >/dev/null 2>&1 || {
			echo "ERROR: flock is required to append $MPLOT_FLAG_FILE" >&2
			return 1
		}
		lock_file="$(dirname "$MPLOT_FLAG_FILE")/.0flag.lock"
		exec {lock_fd}>"$lock_file"
		flock "$lock_fd"

		fragment_header=$(awk 'NR==1 {sub(/\r$/, ""); print; exit}' "$MH_FLAG")
		if [[ "$fragment_header" != $'GWAS\tFLAG' ]]; then
			echo "ERROR: invalid mplot flag fragment header: $MH_FLAG" >&2
			flock -u "$lock_fd"
			exec {lock_fd}>&-
			return 1
		fi

		if [[ -s "$MPLOT_FLAG_FILE" ]]; then
			aggregate_header=$(awk 'NR==1 {sub(/\r$/, ""); print; exit}' "$MPLOT_FLAG_FILE")
			if [[ "$aggregate_header" != $'GWAS\tFLAG' ]]; then
				echo "ERROR: invalid mplot flag header; refusing to overwrite or append: $MPLOT_FLAG_FILE" >&2
				flock -u "$lock_fd"
				exec {lock_fd}>&-
				return 1
			fi
		else
			printf 'GWAS\tFLAG\n' >>"$MPLOT_FLAG_FILE"
		fi

		# Append only rows whose GWAS is not already present.  The aggregate is never
		# rebuilt, truncated, renamed over, or rewritten, so completed reruns preserve it.
		if ! awk -F '\t' '
    FNR==NR {if (FNR>1 && $1!="") seen[$1]=1; next}
    FNR>1 && $1!="" && $2!="" && !seen[$1]++ {print $1 "\t" $2}
  ' "$MPLOT_FLAG_FILE" "$MH_FLAG" >>"$MPLOT_FLAG_FILE"; then
			echo "ERROR: failed to append mplot flag: $MPLOT_FLAG_FILE" >&2
			flock -u "$lock_fd"
			exec {lock_fd}>&-
			return 1
		fi
		flock -u "$lock_fd"
		exec {lock_fd}>&-
	}

	gwas_post_mplot() {
		local genes="${MAGMA_PREFIX}.genes.out" meta_tmp panel_input signal_input sig_input cojo_input
		[[ "$DO_STEP" == all || ",${DO_STEP}," == *,mplot,* ]] || return 0
		if [[ "$REPLACE" != TRUE ]] && gwas_post_mplot_current; then
			gwas_post_append_mplot_flag
			gwas_post_log "Manhattan plot exists and is current: $MH_PNG"
			return 0
		fi
		if [[ "$ADD_PANEL" == magma && ! -s "$genes" ]]; then
			echo "ERROR: --add-panel magma requires MAGMA output: $genes" >&2
			echo "ERROR: run the magma,mplot modules together, or create MAGMA output first" >&2
			return 1
		fi
		gwas_post_prepare_views
		[[ -s "$GWAS_POST_MPLOT_INPUT" ]] || {
			echo "ERROR: Manhattan subset missing" >&2
			return 1
		}
		mkdir -p "$(dirname "$MH_PNG")" "$(dirname "$MH_META")"
		rm -f -- "$MH_META" "${MH_PNG}.meta.tsv"
		[[ "$WRITE_SIG" != TRUE ]] || rm -f -- "$MH_SIG"
		gwas_post_log "Manhattan plot method=$PLOT_METHOD panel=$ADD_PANEL write-sig=$WRITE_SIG: $GWAS_POST_MPLOT_INPUT -> $MH_PNG"
		GWAS_FULL="$FINAL" GWAS_THIN="${THIN_MODE:-FALSE}" Rscript "$MPLOT_R" mplot "$GWAS_POST_MPLOT_INPUT" "$MH_PNG" "$GWAS" "$PLOT_METHOD" "$ADD_PANEL" \
			"$genes" "$PLOT_F" "$CIS_BED" "$MH_PLOT_BED" "$GRCH" "$MH_FLAG" \
			"$ADD_SIGNAL" "$SIGNAL_MATCH_COL" "$SIGNAL_MATCH_VALUE" "$SIGNAL_LOCUS_POS" "$SIGNAL_DISPLAY_COL" \
			"$WRITE_SIG" "$MH_SIG" "$COJO_FILE" "$PLOT_WIDTH" "$PLOT_HEIGHT" "$PLOT_RES"
		[[ -s "$MH_PNG" && -s "$MH_FLAG" ]] || {
			echo "ERROR: Manhattan plot or flag fragment was not created: $MH_PNG $MH_FLAG" >&2
			return 1
		}
		[[ "$WRITE_SIG" != TRUE || -s "$MH_SIG" ]] || {
			echo "ERROR: significant-hit summary was not created: $MH_SIG" >&2
			return 1
		}
		gwas_post_append_mplot_flag
		panel_input=none
		[[ "$ADD_PANEL" != magma ]] || panel_input="$genes"
		meta_tmp="${MH_META}.tmp.$$"
		signal_input=none
		[[ -z "$ADD_SIGNAL" ]] || signal_input="$ADD_SIGNAL"
		sig_input=none
		cojo_input=none
		if [[ "$WRITE_SIG" == TRUE ]]; then
			sig_input="$MH_SIG"
			[[ ! -s "$COJO_FILE" ]] || cojo_input="$COJO_FILE"
		fi
		printf 'key\tvalue\nplot_method\t%s\nadd_panel\t%s\ngrch\t%s\nmagma_threshold\t2.5e-6\nmplot_style\t9\nplot_width\t%s\nplot_height\t%s\nplot_res\t%s\nwrite_sig\t%s\nsig_output\t%s\ncojo_input\t%s\nsignal_input\t%s\nsignal_match_col\t%s\nsignal_match_value\t%s\nsignal_locus_pos\t%s\nsignal_display_col\t%s\ngwas\t%s\ngwas_input\t%s\nmagma_input\t%s\nflag_fragment\t%s\ncreated\t%s\n' \
			"$PLOT_METHOD" "$ADD_PANEL" "$GRCH" "$PLOT_WIDTH" "$PLOT_HEIGHT" "$PLOT_RES" "$WRITE_SIG" "$sig_input" "$cojo_input" \
			"$signal_input" "$SIGNAL_MATCH_COL" "$SIGNAL_MATCH_VALUE" "$SIGNAL_LOCUS_POS" \
			"$SIGNAL_DISPLAY_COL" "$GWAS" "$FINAL" "$panel_input" "$MH_FLAG" "$(date '+%F %T')" >"$meta_tmp"
		printf 'thin_mode\t%s\nthin_chr_max\t%s\ndisplay_input\t%s\n' "${THIN_MODE:-FALSE}" "${THIN_CHR_MAX:-10000}" "$GWAS_POST_MPLOT_INPUT" >>"$meta_tmp"
		mv -f -- "$meta_tmp" "$MH_META"
	}

	# Lead extraction now reads the tiny P-threshold view made by the shared scan.
	gwas_post_prep_lead_inputs() {
		local refs="$1" suffix=".tmp.$$" assoc ma source="${GWAS_POST_LEAD_VIEW:-$FINAL}"
		rm -f -- "${QC_PREFIX}".*.cojo_skip.log
		while IFS=$'\t' read -r _ _ _ assoc ma _; do rm -f "${assoc}${suffix}" "${ma}${suffix}"; done <"$refs"
		gwas_post_zcat "$source" | awk -v FS='\t' -v OFS='\t' -v refs="$refs" -v suffix="$suffix" -v default_n="$GWAS_N" '
    function isnum(x){return x~/^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;return x}
    function get(x){return x in c?$(c[x]):"NA"}
    BEGIN{while((getline line<refs)>0){split(line,a,"\t");chr=a[3];if(chr=="")continue;assoc[chr]=a[4] suffix;ma[chr]=a[5] suffix;print "SNP","CHR","POS","EA","NEA","P">assoc[chr];print "SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","N">ma[chr]}close(refs)}
    NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}next}
    {s=get("SNP");ch=normchr(get("CHR"));pos=get("POS");ea=get("EA");nea=get("NEA");p=get("P");if(!(ch in assoc)||pos!~/^[0-9]+$/||s=="NA"||!isnum(p))next;print s,ch,pos,ea,nea,p>assoc[ch];b=get("BETA");se=get("SE");if(!isnum(b)||!isnum(se))next;n=get("N");if(!isnum(n))n=default_n;print s,ch,pos,ea,nea,get("EAF"),b,se,p,n>ma[ch]}
  '
		while IFS=$'\t' read -r _ _ _ assoc ma _; do
			mv -f "${assoc}${suffix}" "$assoc"
			mv -f "${ma}${suffix}" "$ma"
		done <"$refs"
	}

	gwas_post_reference_cache() {
		local source="$1" kind="$2" stat_key key dir cached lock fd tmp
		stat_key=$(stat -c '%n|%s|%Y' "$source")
		key=$(printf '%s|%s' "$kind" "$stat_key" | sha256sum | awk '{print substr($1,1,20)}')
		dir="$LEAD_REF_CACHE/$kind"
		cached="$dir/$key.bgz"
		lock="$cached.lock"
		mkdir -p "$dir"
		exec {fd}>"$lock"
		flock "$fd"
		if ! gwas_index_valid "$cached"; then
			tmp="$dir/$key.tmp.$$.bgz"
			rm -f -- "$tmp" "$tmp.tbi" "$tmp.csi"
			if [[ "$kind" == pvar ]]; then
				gwas_clean_zcat "$source" | awk -v FS='[ \t]+' -v OFS='\t' '
        /^##/{next} /^#CHROM/{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);c[h]=i}print "#CHROM","POS","ID","REF","ALT";next}
        {ch=$(c["CHROM"]);gsub(/^chr/,"",ch);pos=$(c["POS"]);if(ch!=""&&pos~/^[0-9]+$/)print ch,pos,$(c["ID"]),$(c["REF"]),$(c["ALT"])}' |
					{
						IFS= read -r h
						printf '%s\n' "$h"
						sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k1,1V -k2,2n
					} | bgzip -c >"$tmp"
				tabix -f -s 1 -b 2 -e 2 -S 1 "$tmp"
			else
				gwas_clean_zcat "$source" | awk -v FS='[ \t]+' -v OFS='\t' 'NF>=6{ch=$1;gsub(/^chr/,"",ch);if($4~/^[0-9]+$/)print ch,$2,$3,$4,$5,$6}' |
					sort -T "$GWAS_POST_TMP" -S "$GWAS_POST_SORT_MEMORY" -t $'\t' -k1,1V -k4,4n | bgzip -c >"$tmp"
				tabix -f -s 1 -b 4 -e 4 "$tmp"
			fi
			mv -f -- "$tmp" "$cached"
			[[ ! -s "$tmp.tbi" ]] || mv -f -- "$tmp.tbi" "$cached.tbi"
			[[ ! -s "$tmp.csi" ]] || mv -f -- "$tmp.csi" "$cached.csi"
		fi
		flock -u "$fd"
		exec {fd}>&-
		printf '%s\n' "$cached"
	}

	gwas_post_subset_match_reference() {
		local ref="$1" query="$2" out="$3" kind="$4" source="" ext cached regions
		if [[ "$kind" == pvar ]]; then
			for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst; do [[ -s "${ref}${ext}" ]] && {
				source="${ref}${ext}"
				break
			}; done
		else
			for ext in .bim .bim.gz .bim.bgz; do [[ -s "${ref}${ext}" ]] && {
				source="${ref}${ext}"
				break
			}; done
		fi
		gwas_post_need_file "$source"
		cached=$(gwas_post_reference_cache "$source" "$kind")
		regions="$GWAS_POST_TMP/reference.$kind.regions.txt"
		awk -v FS='\t' 'NR==1{for(i=1;i<=NF;i++){h=toupper($i);sub(/^#/,"",h);if(h=="CHR")c=i;else if(h=="POS")p=i}next} c&&p&&$p~/^[0-9]+$/{ch=$c;gsub(/^chr/,"",ch);print ch ":" $p "-" $p}' "$query" | sort -u >"$regions"
		if [[ "$kind" == pvar ]]; then printf '#CHROM\tPOS\tID\tREF\tALT\n' >"$out"; else : >"$out"; fi
		if [[ -s "$regions" ]]; then
			mapfile -t region_args <"$regions"
			tabix "$cached" "${region_args[@]}" >>"$out"
		fi
		if [[ "$kind" == bim && ! -s "$out" ]]; then printf '0\t.\t0\t0\tN\tN\n' >"$out"; fi
	}

else


	# 🚩 Coordinator validation and discovery
	# Sourced by format.sh: validation, completion checks, discovery and execution.
	# Configuration comes from the entry point; sourcing only defines functions.

	need_arg_value() {
		local opt="$1" val="${2-}"
		if [[ -z "$val" || "$val" == --* ]]; then
			echo "ERROR: missing value for $opt" >&2
			usage >&2
			exit 2
		fi
	}

	upper() { echo "$1" | tr '[:lower:]' '[:upper:]'; }
	has_step() { [[ ",$step," == *",$1,"* || "$step" == "all" ]]; }
	wants_magma() { [[ ",$step," == *",magma,"* ]]; }

	gwas_format_validate_options() {
		replace=$(upper "$replace")
		fill_eaf=$(upper "$fill_eaf")
		run_cmd=$(upper "$run_cmd")
		is_bsub=$(upper "$is_bsub")
		refGen_pop=$(upper "$refGen_pop")
		if [[ -z "$refGen_pop" || "$refGen_pop" == *[!A-Z0-9_-]* ]]; then
			echo "ERROR: --refgen-pop must be a simple population label: $refGen_pop" >&2
			exit 2
		fi
		foreground=$(upper "$foreground")
		liftOver=$(upper "$liftOver")
		delete_raw=$(upper "$delete_raw")
		hm3_mode=$(upper "$hm3_mode")
		thin=$(upper "$thin")
		write_sig=$(upper "$write_sig")
		step=$(echo "$step" | tr '[:upper:]' '[:lower:]')
		category=$(echo "$category" | tr '[:upper:]' '[:lower:]')
		add_panel=$(echo "$add_panel" | tr '[:upper:]' '[:lower:]')
		plot_height=$(echo "$plot_height" | tr '[:upper:]' '[:lower:]')
		[[ "$add_panel" == none || "$add_panel" == magma ]] || {
			echo "ERROR: --add-panel must be none or magma: $add_panel" >&2
			exit 2
		}

		[[ "$hm3_mode" == "TRUE" || "$hm3_mode" == "FALSE" ]] || {
			echo "ERROR: --hm3 must be TRUE or FALSE; use --hm3-file FILE for the reference list" >&2
			exit 2
		}
		[[ "$thin" == TRUE || "$thin" == FALSE ]] || {
			echo "ERROR: --thin must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$thin_chr_max" =~ ^[1-9][0-9]*$ ]] || {
			echo "ERROR: --thin-chr-max must be a positive integer" >&2
			exit 2
		}
		[[ "$replace" == "TRUE" || "$replace" == "FALSE" ]] || {
			echo "ERROR: --replace must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$fill_eaf" == "TRUE" || "$fill_eaf" == "FALSE" ]] || {
			echo "ERROR: --fill-eaf must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$run_cmd" == "TRUE" || "$run_cmd" == "FALSE" ]] || {
			echo "ERROR: --run-cmd must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$is_bsub" == "TRUE" || "$is_bsub" == "FALSE" ]] || {
			echo "ERROR: --submit-bsub must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$foreground" == "TRUE" || "$foreground" == "FALSE" ]] || {
			echo "ERROR: --foreground must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$liftOver" == "TRUE" || "$liftOver" == "FALSE" ]] || {
			echo "ERROR: --liftover must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$delete_raw" == "TRUE" || "$delete_raw" == "FALSE" ]] || {
			echo "ERROR: --delete-raw must be TRUE or FALSE" >&2
			exit 2
		}
		[[ "$write_sig" == "TRUE" || "$write_sig" == "FALSE" ]] || {
			echo "ERROR: --write-sig must be TRUE or FALSE" >&2
			exit 2
		}
		[[ -n "$category" && "$category" != "." && "$category" != ".." && "$category" != *[!a-z0-9._-]* ]] || {
			echo "ERROR: --category must be a simple folder name such as common or rare: $category" >&2
			exit 2
		}
		if [[ -n "$fill_n" ]]; then
			awk -v n="$fill_n" 'BEGIN{exit !(n ~ /^[0-9]+([.][0-9]+)?$/ && n+0>0)}' || {
				echo "ERROR: --fill-n must be a positive number: $fill_n" >&2
				exit 2
			}
		fi
		IFS=',' read -r -a requested_steps <<<"$step"
		for requested_step in "${requested_steps[@]}"; do
			case "$requested_step" in
				format | thin | magma | liftover | cis | lead | mplot | h2 | pgs | all) ;;
				*)
					echo "ERROR: step/module must contain format|thin|magma|liftover|cis|lead|mplot|h2|pgs, or all" >&2
					exit 2
					;;
			esac
		done

		[[ -z "$raw_file_arg" || (-f "$raw_file_arg" && -s "$raw_file_arg" && -n "$gwas_arg" && "$gwas_arg" != *,*) ]] || {
			echo 'ERROR: --raw-file requires an existing file and one --gwas NAME' >&2
			exit 2
		}
		[[ -z "$n_total" || "$n_total" =~ ^[1-9][0-9]*$ ]] || {
			echo 'ERROR: --n-total must be a positive integer' >&2
			exit 2
		}
		[[ -z "$n_total" || "$hm3_mode" == FALSE ]] || {
			echo 'ERROR: --n-total currently requires --hm3 FALSE' >&2
			exit 2
		}
		case "$h2_sex" in unknown | male | female | mixed) ;; *)
			echo 'ERROR: invalid --h2-sex' >&2
			exit 2
			;;
		esac
		if has_step liftover && [[ "$liftOver" == TRUE ]]; then
			[[ "$step" == liftover || "$step" == format,liftover ]] || {
				echo 'ERROR: run format,liftover first, then downstream modules in a separate invocation with --grch auto' >&2
				exit 2
			}
		fi

		if ! [[ "$jobs" =~ ^[1-9][0-9]*$ ]]; then
			echo "ERROR: --jobs must be a positive integer: $jobs" >&2
			exit 2
		fi
		if ! [[ "$pgs_threads" =~ ^[1-9][0-9]*$ ]]; then
			echo "ERROR: --pgs-threads must be a positive integer: $pgs_threads" >&2
			exit 2
		fi
		awk -v x="$plot_width" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+$/ && x+0>0)}' || {
			echo "ERROR: --plot-width must be a positive number of inches: $plot_width" >&2
			exit 2
		}
		if [[ "$plot_height" != auto ]]; then
			awk -v x="$plot_height" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+$/ && x+0>0)}' || {
				echo "ERROR: --plot-height must be auto or a positive number of inches: $plot_height" >&2
				exit 2
			}
		fi
		if ! [[ "$plot_res" =~ ^[1-9][0-9]*$ ]]; then
			echo "ERROR: --plot-res must be a positive integer: $plot_res" >&2
			exit 2
		fi
		awk -v p="$p_lead" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' || {
			echo "ERROR: --p-lead must be a number in (0,1]: $p_lead" >&2
			exit 2
		}
		awk -v p="$p_hm3" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' || {
			echo "ERROR: --p-hm3 must be a number in (0,1]: $p_hm3" >&2
			exit 2
		}
		if has_step format && [[ "$hm3_mode" == TRUE && "$thin" == TRUE ]]; then
			awk -v p="$p_hm3" 'BEGIN{exit !(p+0>=0.001)}' || {
				echo 'ERROR: --hm3 TRUE --thin TRUE requires --p-hm3 >= 1e-3 so formatting retains all thin candidates.' >&2
				exit 2
			}
		fi
		if ! [[ "$lead_window" =~ ^[1-9][0-9]*$ ]]; then
			echo "ERROR: --lead-window must be a positive integer: $lead_window" >&2
			exit 2
		fi

	}

	label_from_dir_raw() {
		local p="$1" b
		p="${p%/}"
		b=$(basename "$p")
		[[ "$b" != "raw" ]] || b=$(basename "$(dirname "$p")")
		echo "$b"
	}

	label_from_dir_clean() {
		local p="$1"
		p="${p%/}"
		basename "$(dirname "$(dirname "$(dirname "$p")")")"
	}

	gwas_format_maybe_background() {
		if [[ "$run_cmd" == "TRUE" && "$is_bsub" != "TRUE" && "$foreground" != "TRUE" ]]; then
			background_log="$dir_cmd/gwas_post.background.log"
			[[ -s "$phef" ]] || {
				echo "ERROR: missing or empty file: $phef" >&2
				exit 1
			}
			# shellcheck source=/mnt/d/scripts/0f/phenotype.sh
			source "$phef"
			if declare -F phe_check_existing_tasks >/dev/null 2>&1; then
				# Only reject another run of the same module set.  Different modules use
				# separate coordinator files, generated commands, and per-module logs.
				phe_check_existing_tasks "format_gwas.sh $step"
			fi
			phe_run_background --title "background gwas_post $label/$category/$step" "$background_log" bash "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" --foreground TRUE
			exit 0
		fi

	}

	log() { echo "[$(date '+%F %T')] $*" >&2; }
	need_file() { [[ -s "$1" ]] || {
		echo "ERROR: missing or empty file: $1" >&2
		exit 1
	}; }
	need_dir() { [[ -d "$1" ]] || {
		echo "ERROR: missing directory: $1" >&2
		exit 1
	}; }
	need_refgen_clump() {
		local d="$1"
		if [[ -f "${d}.pgen" && -f "${d}.psam" ]] &&
			[[ -f "${d}.pvar" || -f "${d}.pvar.zst" ]]; then
			return 0
		fi
		if [[ -d "$d" ]]; then
			compgen -G "$d/chr*.pgen" >/dev/null || {
				echo "ERROR: no chr*.pgen files found in refGen_clump: $d" >&2
				exit 1
			}
			{ compgen -G "$d/chr*.pvar" >/dev/null || compgen -G "$d/chr*.pvar.zst" >/dev/null; } || {
				echo "ERROR: no chr*.pvar[.zst] files found in refGen_clump: $d" >&2
				exit 1
			}
			compgen -G "$d/chr*.psam" >/dev/null || {
				echo "ERROR: no chr*.psam files found in refGen_clump: $d" >&2
				exit 1
			}
			return 0
		fi
		compgen -G "${d}chr*.pgen" >/dev/null || {
			echo "ERROR: no chr*.pgen files found in refGen_clump prefix: $d" >&2
			exit 1
		}
		{ compgen -G "${d}chr*.pvar" >/dev/null || compgen -G "${d}chr*.pvar.zst" >/dev/null; } || {
			echo "ERROR: no chr*.pvar[.zst] files found in refGen_clump prefix: $d" >&2
			exit 1
		}
		compgen -G "${d}chr*.psam" >/dev/null || {
			echo "ERROR: no chr*.psam files found in refGen_clump prefix: $d" >&2
			exit 1
		}
	}

	need_refgen_cojo() {
		local d="$1"
		if [[ -f "${d}.bed" && -f "${d}.bim" && -f "${d}.fam" ]]; then
			return 0
		fi
		if [[ -d "$d" ]]; then
			compgen -G "$d/chr*.bed" >/dev/null || {
				echo "ERROR: no chr*.bed files found in refGen_cojo: $d" >&2
				exit 1
			}
			compgen -G "$d/chr*.bim" >/dev/null || {
				echo "ERROR: no chr*.bim files found in refGen_cojo: $d" >&2
				exit 1
			}
			compgen -G "$d/chr*.fam" >/dev/null || {
				echo "ERROR: no chr*.fam files found in refGen_cojo: $d" >&2
				exit 1
			}
			return 0
		fi
		compgen -G "${d}chr*.bed" >/dev/null || {
			echo "ERROR: no chr*.bed files found in refGen_cojo prefix: $d" >&2
			exit 1
		}
		compgen -G "${d}chr*.bim" >/dev/null || {
			echo "ERROR: no chr*.bim files found in refGen_cojo prefix: $d" >&2
			exit 1
		}
		compgen -G "${d}chr*.fam" >/dev/null || {
			echo "ERROR: no chr*.fam files found in refGen_cojo prefix: $d" >&2
			exit 1
		}
	}

	need_pgs_pfiles() {
		local d="${1%/}"
		[[ -d "$d" ]] || {
			echo "ERROR: missing PGS pfile directory: $d" >&2
			exit 1
		}
		compgen -G "$d/chr*.pgen" >/dev/null || {
			echo "ERROR: no chr*.pgen files found in PGS pfile directory: $d" >&2
			exit 1
		}
		{ compgen -G "$d/chr*.pvar" >/dev/null || compgen -G "$d/chr*.pvar.zst" >/dev/null; } || {
			echo "ERROR: no chr*.pvar[.zst] files found in PGS pfile directory: $d" >&2
			exit 1
		}
		compgen -G "$d/chr*.psam" >/dev/null || {
			echo "ERROR: no chr*.psam files found in PGS pfile directory: $d" >&2
			exit 1
		}
	}

	ensure_magma_resources() {
		for ext in bed bim fam; do need_file "${magma_ref}.${ext}"; done
		if [[ -n "$gene_loc" ]]; then
			need_file "$gene_loc"
		elif [[ "$grch" == auto ]]; then
			need_file "$dir0/files/NCBI.37.gene.loc"
			need_file "$dir0/files/NCBI.38.gene.loc"
		else
			echo "ERROR: no MAGMA gene-location file resolved for GRCh$grch" >&2
			exit 1
		fi
		need_file "$synonyms"
	}
	q() { printf '%q' "$1"; }

	# A result is reusable only for the requested trait, phase and lead settings.
	# In particular, an autosome-only marker must not satisfy all chromosomes.
	gwas_lead_marker_matches() {
		local marker="$1" trait="$2" phase="$3" p="$4" window="$5" chromosomes="$6"
		[[ -s "$marker" ]] || return 1
		awk -F '\t' -v trait="$trait" -v phase="$phase" -v p="$p" \
			-v window="$window" -v chromosomes="$chromosomes" '
    NR==1 {header=($0=="GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS")}
    NR==2 {ok=(header && NF==6 && $1==trait && $2==phase &&
      $3==p && $4==window && $5==chromosomes &&
      ($6=="complete" || $6=="no_significant_variants" || $6=="no_reference_matched_variants" ||
       (phase=="cojo" && $6=="no_snps_selected")))}
    END {exit !(NR==2 && ok)}
  ' "$marker"
	}

	# Embedded in each worker. GCTA reports an empty reference-QC intersection as
	# exit 1; retain that evidence without treating an empty chromosome as a crash.
	gwas_post_run_cojo() {
		local prefix="$1" audit="$2" rc ext empty_reason=""
		shift 2
		GCTA_COJO_HAS_MATCH=FALSE
		GCTA_COJO_NO_SNPS_SELECTED=FALSE
		mkdir -p "$(dirname "$audit")"
		rm -f -- "${audit}.tsv" "${audit}.log" "${audit}.freq.badsnps" "${audit}.badsnps"
		rm -f -- "${prefix}.freq.badsnps" "${prefix}.badsnps"
		if run_tool "$prefix" "$@"; then
			GCTA_COJO_HAS_MATCH=TRUE
			# GCTA can finish selection successfully without producing a .jma.cojo.
			# Require both the explicit empty result and normal completion; missing
			# output alone must never turn a truncated or broken run into success.
			if grep -Fxq 'No SNPs have been selected.' "${prefix}.log" &&
				grep -Eq '^Analysis finished at ' "${prefix}.log"; then
				cp -- "${prefix}.log" "${audit}.log"
				printf 'STATUS\tREASON\tSELECTED_SNPS\nEMPTY\tno_snps_selected\t0\n' >"${audit}.tsv"
				rm -f -- "${prefix}.jma.cojo" "${prefix}.cma.cojo" "${prefix}.ldr.cojo"
				GCTA_COJO_NO_SNPS_SELECTED=TRUE
				gwas_post_log "COJO completed with no SNPs selected: $prefix; evidence: ${audit}.tsv"
			elif [[ ! -s "${prefix}.jma.cojo" && ! -s "${prefix}.ldr.cojo" ]]; then
				echo "ERROR: COJO returned success without results or an explicit completed empty selection: ${prefix}.log" >&2
				return 1
			fi
			return 0
		else
			rc=$?
		fi

		# Recognize only explicit empty results from known QC stages. In a small
		# chromosome subset every candidate can be fixed in the reference, causing
		# the MAF filter to stop before GCTA even reads the summary statistics.
		if [[ "$rc" == 1 ]] &&
			grep -Fxq 'Error: none of the SNPs in the GWAS summary data can be found in the genotype data.' "${prefix}.log" &&
			grep -Fxq 'Matching the GWAS meta-analysis results to the genotype data ...' "${prefix}.log" &&
			grep -Eq '^GWAS summary statistics of [1-9][0-9]* SNPs read from ' "${prefix}.log"; then
			empty_reason=no_reference_matched_variants_after_gcta_qc
		elif [[ "$rc" == 1 ]] && awk '
    /^Genotype data for [1-9][0-9]* individuals and [1-9][0-9]* SNPs to be included from / {loaded=1}
    /^Calculating allele frequencies \.\.\.$/ {frequencies=loaded}
    $0=="Error: no SNP is retained for analysis." && frequencies &&
      previous ~ /^Filtering SNPs with MAF > [0-9.eE+-]+ \.\.\.$/ {empty=1}
    {previous=$0}
    END {exit !empty}
  ' "${prefix}.log"; then
			empty_reason=no_reference_variants_after_gcta_maf_filter
		fi
		if [[ -n "$empty_reason" ]]; then
			cp -- "${prefix}.log" "${audit}.log"
			for ext in freq.badsnps badsnps; do
				if [[ -s "${prefix}.${ext}" ]]; then cp -- "${prefix}.${ext}" "${audit}.${ext}"; fi
			done
			printf 'STATUS\tREASON\tUSABLE_SNPS\nEMPTY\t%s\t0\n' "$empty_reason" >"${audit}.tsv"
			rm -f -- "${prefix}.err" "${prefix}.jma.cojo" "${prefix}.cma.cojo" "${prefix}.ldr.cojo"
			gwas_post_log "COJO has no usable reference SNPs after GCTA QC: $prefix; evidence: ${audit}.tsv"
			return 0
		fi
		return "$rc"
	}

	write_pgs_step2_cmd() {
		local root="$dir_out/pgs" script="$dir_out/pgs/pgs.step2.cmd"
		mkdir -p "$root"
		rm -f -- "$root/pgs.step.cmd"
		{
			cat <<'PGS_STEP_HEADER'
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
PGS_STEP_HEADER
			printf 'DEFAULT_PROJECT=%q\nDEFAULT_CATEGORY=%q\nDEFAULT_LABEL=%q\n' "$dir_out" "$category" "$label"
			cat <<'PGS_STEP_BODY'
PROJECT=${1:-$DEFAULT_PROJECT}
CATEGORY=${2:-$DEFAULT_CATEGORY}
LABEL=${3:-$DEFAULT_LABEL}
BATCH_SIZE=${4:-32}
[[ "$BATCH_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: batch size must be a positive integer: $BATCH_SIZE" >&2; exit 2; }
ROOT="$PROJECT/pgs"
OUTPUT="$ROOT/$LABEL.pgs.gz"
ALLELE_OUTPUT="$ROOT/$LABEL.ALLELE_CT.tsv"
DONE="$ROOT/$LABEL.pgs.done"
MANIFEST="$ROOT/$LABEL.pgs.files.tsv"
TMP_ROOT="$ROOT/.tmp"

mkdir -p "$ROOT" "$TMP_ROOT"
TMP_DIR=$(mktemp -d "$TMP_ROOT/merge.XXXXXX")
case "$TMP_DIR" in "$TMP_ROOT"/merge.*) ;; *) echo "ERROR: unsafe temporary directory: $TMP_DIR" >&2; exit 1;; esac
declare -a ACTIVE_STREAM_DIRS=()
cleanup(){
  local d
  for d in "${ACTIVE_STREAM_DIRS[@]}"; do
    case "$d" in /tmp/gwas-post-pgs-streams.*)
      [[ ! -d "$d" ]] || { rm -f -- "$d"/*; rmdir -- "$d" 2>/dev/null || true; }
      ;;
    esac
  done
  rm -rf -- "$TMP_DIR"
  rmdir -- "$TMP_ROOT" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

for tool in awk cp find gzip head mkfifo mktemp mv paste rm rmdir sort; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required command not found: $tool" >&2; exit 1; }
done
compress=(gzip -c)
if command -v pigz >/dev/null 2>&1; then compress=(pigz -c -p 4); fi

manifest_tmp="$TMP_DIR/files.tsv"
allele_tmp="$TMP_DIR/$LABEL.ALLELE_CT.tsv"
printf 'GWAS\tSTATUS\tPGS\n' > "$manifest_tmp"
printf 'trait\tALLELE_CT\n' > "$allele_tmp"
declare -a inputs=() preview=()
missing=0
empty=0
expected=0

while IFS= read -r -d '' score_file; do
  trait=${score_file##*/}
  trait=${trait%.jma.cojo}
  trait_dir="$PROJECT/$CATEGORY/$trait"
  pgs_file="$trait_dir/pgs/$trait.pgs.gz"
  done_file="$trait_dir/pgs/$trait.pgs.done"
  err_file="$trait_dir/$trait.pgs.err"
  ((expected+=1))

  if [[ ! -s "$pgs_file" || ! -s "$done_file" || -s "$err_file" ]]; then
    printf '%s\tincomplete\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    echo "ERROR: incomplete PGS for $trait: output=$pgs_file done=$done_file err=$err_file" >&2
    ((missing+=1))
    continue
  fi
  gzip -t -- "$pgs_file"
  preview=()
  mapfile -t preview < <(set +o pipefail; gzip -cd -- "$pgs_file" | head -n 2)
  expected_header=$(printf '#IID\t%s.ALLELE_CT\t%s.SCORE_SUM' "$trait" "$trait")
  [[ "${preview[0]:-}" == "$expected_header" ]] || {
    echo "ERROR: unexpected PGS header for $trait: ${preview[0]:-<empty>}" >&2
    exit 1
  }
  if [[ -z "${preview[1]:-}" ]]; then
    printf '%s\tempty_no_matched_variants\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    ((empty+=1))
  else
    IFS=$'\t' read -r first_iid allele_ct first_score extra <<< "${preview[1]}"
    [[ -n "$first_iid" && "$allele_ct" =~ ^[0-9]+$ && -n "$first_score" && -z "${extra:-}" ]] || {
      echo "ERROR: invalid first PGS data row for $trait: ${preview[1]}" >&2
      exit 1
    }
    printf '%s\t%s\n' "$trait" "$allele_ct" >> "$allele_tmp"
    printf '%s\tincluded\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    inputs+=("$pgs_file")
  fi
done < <(find "$PROJECT/$CATEGORY" -mindepth 3 -maxdepth 3 -type f -path '*/gwas/*.jma.cojo' -print0 | sort -z)

(( expected > 0 )) || { echo "ERROR: no .jma.cojo inputs found under $PROJECT/$CATEGORY" >&2; exit 1; }
if (( missing > 0 )); then
  echo "ERROR: $missing of $expected expected PGS files are incomplete; merge not written." >&2
  exit 1
fi
(( ${#inputs[@]} > 0 )) || { echo "ERROR: every completed PGS is empty; merge not written." >&2; exit 1; }

merge_trait_batch(){
  local out="$1"; shift
  local n=$# f fifo pid rc producer_rc=0 i=0
  local stream_dir
  local -a fifos=() producer_pids=()
  stream_dir=$(mktemp -d /tmp/gwas-post-pgs-streams.XXXXXX)
  case "$stream_dir" in /tmp/gwas-post-pgs-streams.*) ;; *) echo "ERROR: unsafe stream directory: $stream_dir" >&2; return 1;; esac
  ACTIVE_STREAM_DIRS+=("$stream_dir")
  for f in "$@"; do
    fifo="$stream_dir/in.$(printf '%04d' "$i")"
    mkfifo -- "$fifo"
    fifos+=("$fifo")
    gzip -cd -- "$f" > "$fifo" &
    producer_pids+=("$!")
    ((i+=1))
  done
  set +e
  paste "${fifos[@]}" | awk -v FS='\t' -v OFS='\t' -v n="$n" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    {
      expected=3*n
      if(NF!=expected){print "ERROR: batch field count mismatch at row " NR ": expected " expected ", got " NF > "/dev/stderr";exit 2}
      if(NR==1){
        printf "eid"
        for(i=0;i<n;i++){
          base=3*i
          if($(base+1)!="#IID" || $(base+3)!~ /[.]SCORE_SUM$/){print "ERROR: invalid PGS header in batch input " i+1 > "/dev/stderr";exit 2}
          printf "%s%s",OFS,$(base+3)
        }
        printf "\n"
        next
      }
      eid=$1
      printf "%s",eid
      for(i=0;i<n;i++){
        base=3*i
        if($(base+1)!=eid){print "ERROR: IID mismatch in batch at row " NR ", input " i+1 > "/dev/stderr";exit 2}
        score=$(base+3)
        if(!isnum(score)){print "ERROR: non-numeric SCORE_SUM in batch at row " NR ", input " i+1 ": " score > "/dev/stderr";exit 2}
        score+=0
        if(score>-0.0005 && score<0.0005)score=0
        printf "%s%.3f",OFS,score
      }
      printf "\n"
    }
  ' | "${compress[@]}" > "$out"
  rc=$?
  for pid in "${producer_pids[@]}"; do wait "$pid" || producer_rc=1; done
  set -e
  rm -f -- "${fifos[@]}"
  rmdir -- "$stream_dir"
  (( rc == 0 && producer_rc == 0 )) || return 1
  gzip -t -- "$out"
}

declare -a chunks=() chunk_counts=() batch=()
chunk_index=0
for f in "${inputs[@]}"; do
  batch+=("$f")
  if (( ${#batch[@]} == BATCH_SIZE )); then
    chunk="$TMP_DIR/chunk.$(printf '%04d' "$chunk_index").gz"
    merge_trait_batch "$chunk" "${batch[@]}"
    chunks+=("$chunk"); chunk_counts+=("${#batch[@]}")
    batch=(); ((chunk_index+=1))
  fi
done
if (( ${#batch[@]} > 0 )); then
  chunk="$TMP_DIR/chunk.$(printf '%04d' "$chunk_index").gz"
  merge_trait_batch "$chunk" "${batch[@]}"
  chunks+=("$chunk"); chunk_counts+=("${#batch[@]}")
fi

final_tmp="$TMP_DIR/$LABEL.pgs.gz"
if (( ${#chunks[@]} == 1 )); then
  cp -- "${chunks[0]}" "$final_tmp"
else
  stream_dir=$(mktemp -d /tmp/gwas-post-pgs-streams.XXXXXX)
  case "$stream_dir" in /tmp/gwas-post-pgs-streams.*) ;; *) echo "ERROR: unsafe stream directory: $stream_dir" >&2; exit 1;; esac
  ACTIVE_STREAM_DIRS+=("$stream_dir")
  fifos=(); producer_pids=(); i=0; producer_rc=0
  for f in "${chunks[@]}"; do
    fifo="$stream_dir/in.$(printf '%04d' "$i")"
    mkfifo -- "$fifo"
    fifos+=("$fifo")
    gzip -cd -- "$f" > "$fifo" &
    producer_pids+=("$!")
    ((i+=1))
  done
  counts=$(IFS=,; echo "${chunk_counts[*]}")
  set +e
  paste "${fifos[@]}" | awk -v FS='\t' -v OFS='\t' -v counts="$counts" '
    BEGIN{nchunk=split(counts,n,","); expected=0; for(j=1;j<=nchunk;j++)expected+=1+n[j]}
    {
      if(NF!=expected){print "ERROR: final field count mismatch at row " NR ": expected " expected ", got " NF > "/dev/stderr";exit 2}
      eid=$1; pos=1
      printf "%s",eid
      for(j=1;j<=nchunk;j++){
        if($(pos)!=eid){print "ERROR: IID mismatch between chunks at row " NR ", chunk " j > "/dev/stderr";exit 2}
        for(k=1;k<=n[j];k++)printf "%s%s",OFS,$(pos+k)
        pos+=1+n[j]
      }
      printf "\n"
    }
  ' | "${compress[@]}" > "$final_tmp"
  rc=$?
  for pid in "${producer_pids[@]}"; do wait "$pid" || producer_rc=1; done
  set -e
  rm -f -- "${fifos[@]}"
  rmdir -- "$stream_dir"
  (( rc == 0 && producer_rc == 0 )) || { echo "ERROR: final PGS merge pipeline failed" >&2; exit 1; }
fi
gzip -t -- "$final_tmp"
mv -f -- "$final_tmp" "$OUTPUT"
mv -f -- "$allele_tmp" "$ALLELE_OUTPUT"
mv -f -- "$manifest_tmp" "$MANIFEST"
done_tmp="$TMP_DIR/$LABEL.pgs.done"
printf 'LABEL\tSTATUS\tEXPECTED\tINCLUDED\tEMPTY\tTIME\n%s\tcomplete\t%s\t%s\t%s\t%s\n' \
  "$LABEL" "$expected" "${#inputs[@]}" "$empty" "$(date '+%F %T')" > "$done_tmp"
mv -f -- "$done_tmp" "$DONE"
echo "PGS merge done: $OUTPUT; allele counts: $ALLELE_OUTPUT (expected=$expected included=${#inputs[@]} empty=$empty)"
PGS_STEP_BODY
		} >"$script"
		# This merge command is invoked with bash, including on mounted drives.
		log "Wrote $label PGS merge command: $script"
	}

	# A tabix sidecar makes completion checks O(1). Legacy gzip files retain the
	# full integrity fallback until the one-time migration creates their index.
	gzip_ok() {
		local file="$1"
		[[ -s "$file" ]] || return 1
		if { [[ -s "${file}.tbi" ]] || [[ -s "${file}.csi" ]]; } && command -v tabix >/dev/null 2>&1; then
			tabix -l "$file" >/dev/null 2>&1
		else
			gzip -t "$file" >/dev/null 2>&1
		fi
	}

	cis_output_complete() {
		local file="$1" index=""
		[[ -s "$file" ]] || return 1
		[[ -s "${file}.tbi" ]] && index="${file}.tbi"
		[[ -z "$index" && -s "${file}.csi" ]] && index="${file}.csi"
		[[ -n "$index" && ! "$file" -nt "$index" ]] || return 1
		tabix -l "$file" >/dev/null 2>&1
	}

	magma_output_complete() {
		local magma_dir="$1" magma_prefix="$2"
		[[ -s "$magma_dir/magma.done" && -s "$magma_prefix.genes.out" && -s "$magma_prefix.genes.raw" &&
			-s "${magma_dir%/*}/qc/$(basename "$magma_prefix").magma.chromosomes.log" ]]
	}

	mplot_output_complete() {
		local png="$1" final="$2" genes="$3" meta="$4" expected_grch="$5" flag="$6" aggregate="$7" sig="$8" cojo="$9"
		local signal_input=none thin_output="${final%.gz}.thin.gz"
		if [[ "${thin:-FALSE}" == TRUE ]]; then
			if gwas_thin_plot_only "$step"; then
				gwas_thin_plot_ready "$final" "$thin_output" || return 1
			else
				gwas_thin_complete "$final" "$thin_output" "$expected_grch" "$thin_chr_max" "$hm3_file" "${hm3_pos//\{grch\}/$expected_grch}" "$thin_r" "$phe_r" "$meta" || return 1
			fi
			[[ "$png" -nt "$thin_output" ]] || return 1
		fi
		[[ -z "$add_signal" ]] || signal_input="$add_signal"
		[[ -s "$png" && -s "$meta" && -s "$flag" && -s "$aggregate" && "$png" -nt "$final" ]] || return 1
		awk -F '\t' -v panel="$add_panel" -v grch="$expected_grch" -v signal="$signal_input" \
			-v match_col="$signal_match_col" -v match_value="$signal_match_value" \
			-v locus_pos="$signal_locus_pos" -v display_col="$signal_display_col" -v write_sig="$write_sig" \
			-v thin="${thin:-FALSE}" -v cap="${thin_chr_max:-10000}" \
			-v plot_width="$plot_width" -v plot_height="$plot_height" -v plot_res="$plot_res" '
    $1=="plot_method"&&$2=="self"{m=1}
    $1=="add_panel"&&$2==panel{p=1}
    $1=="grch"&&$2==grch{g=1}
    $1=="magma_threshold"&&$2+0==2.5e-6{t=1}
    $1=="mplot_style"&&$2==9{s=1}
    $1=="thin_mode"&&$2==thin{tn=1}
    $1=="thin_chr_max"&&$2==cap{tc=1}
    $1=="plot_width"&&$2+0==plot_width+0{x=1}
    $1=="plot_height"&&$2==plot_height{y=1}
    $1=="plot_res"&&$2+0==plot_res+0{r=1}
    $1=="write_sig"&&$2==write_sig{w=1}
    $1=="signal_input"&&$2==signal{a=1}
    $1=="signal_match_col"&&$2==match_col{b=1}
    $1=="signal_match_value"&&$2==match_value{c=1}
    $1=="signal_locus_pos"&&$2==locus_pos{d=1}
    $1=="signal_display_col"&&$2==display_col{e=1}
    END{exit !(m&&p&&g&&t&&s&&x&&y&&r&&w&&a&&b&&c&&d&&e&&tn&&(thin!="TRUE"||tc))}' "$meta" || return 1
		if [[ "$add_panel" == magma ]]; then
			[[ -s "$genes" && "$png" -nt "$genes" ]] || return 1
		fi
		if [[ -n "$add_signal" ]]; then [[ -s "$add_signal" && "$png" -nt "$add_signal" ]] || return 1; fi
		if [[ "$write_sig" == TRUE ]]; then
			[[ -s "$sig" && "$sig" -nt "$final" ]] || return 1
			[[ ! -s "$cojo" || "$sig" -nt "$cojo" ]] || return 1
			[[ "$add_panel" != magma || "$sig" -nt "$genes" ]] || return 1
		fi
	}

	pgs_output_complete() {
		local output="$1" done_file="$2"
		[[ -s "$output" && -s "$done_file" ]]
	}

	prune_pgs_dir() {
		local d="$1" gwas="$2" f keep_logs=FALSE meta="$1/$2.pgs.meta.tsv"
		[[ -d "$d" ]] || return 0
		if [[ -s "$meta" ]] && awk -F '\t' '$1=="matched_variants"&&$2=="none"{found=1} END{exit !found}' "$meta"; then
			keep_logs=TRUE
		fi
		for f in "$d/${gwas}.chr"*; do
			[[ -f "$f" ]] || continue
			if [[ "$keep_logs" == TRUE && "$f" == *.log ]]; then continue; fi
			rm -f -- "$f"
		done
	}

	prune_magma_dir() {
		local d="$1"
		[[ -d "$d" ]] || return 0
		find "$d" -mindepth 1 -maxdepth 1 -type f \
			! -name '*.genes.out' ! -name '*.genes.raw' ! -name '*.log' \
			! -name 'magma.meta.tsv' ! -name 'magma.done' -delete
	}

	cleanup_failed_output_dirs() {
		[[ "$step" == "all" ]] || return 0
		log "Automatic failed-folder deletion is disabled for the category-first layout."
	}

	# Strip common GWAS file extensions without destroying phenotype names containing dots.
	gwas_name_from_file() {
		local b
		b=$(basename "$1")
		b=${b%.gz}
		b=${b%.bgz}
		b=${b%.tsv}
		b=${b%.txt}
		b=${b%.sumstats}
		b=${b%.assoc}
		echo "$b"
	}

	list_raw_files() {
		[[ -d "$dir_raw" ]] || return 0
		# Bash globbing is more reliable than find on /mnt/* (DrvFS) immediately
		# after a Windows-side download/rename becomes visible to WSL.
		(
			shopt -s nullglob
			local f
			find "$dir_raw" -mindepth 4 -maxdepth 4 -type f -path "*/$category/*/raw/*" \
				\( -name '*.gz' -o -name '*.bgz' -o -name '*.tsv' -o -name '*.txt' -o -name '*.sumstats' -o -name '*.assoc' \) \
				-size +0c 2>/dev/null
			for f in "$dir_raw"/*.gz "$dir_raw"/*.bgz "$dir_raw"/*.tsv "$dir_raw"/*.txt "$dir_raw"/*.sumstats "$dir_raw"/*.assoc; do
				[[ -f "$f" && -s "$f" && "$f" != *.aria2 ]] && printf '%s\n' "$f"
			done
		) | sort -u -V
	}

	list_names_from_dir() {
		local d="$1"
		[[ -d "$d" ]] || return 0
		find "$d" -type f \( -name '*.gz' -o -name '*.bgz' \) \
			! -name '*.small.gz' ! -name '*.hm3.gz' ! -name '*.thin.gz' ! -name '*.cis.gz' ! -name '*.sig.tsv.gz' ! -name '*.lead.tsv.gz' \
			-size +0c 2>/dev/null |
			while read -r f; do
				[[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
				gwas_name_from_file "$f"
			done | sort -u -V
	}

	list_liftover_names() {
		[[ -d "$dir_clean" ]] || return 0
		find "$dir_clean" -type f \( -name '*.hm3.gz' -o -name '*.small.gz' -o -name '*.source.grch37.gz' \) -size +0c 2>/dev/null |
			while read -r f; do
				[[ "$(basename "$(dirname "$f")")" == "gwas" || "$(basename "$(dirname "$f")")" == qc ]] || continue
				b=$(basename "$f")
				b=${b%.gz}
				b=${b%.hm3}
				b=${b%.small}
				b=${b%.source.grch37}
				echo "$b"
			done | sort -u -V
	}

	list_pgs_names() {
		[[ -d "$dir_clean" ]] || return 0
		find "$dir_clean" -type f -name '*.jma.cojo' -size +0c 2>/dev/null |
			while read -r f; do
				[[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
				b=$(basename "$f")
				echo "${b%.jma.cojo}"
			done | sort -u -V
	}

	raw_file_for_name() {
		local g="$1" f
		for f in \
			"$dir_out/$category/$g/raw/$g.gz" "$dir_out/$category/$g/raw/$g.bgz" \
			"$dir_out/$category/$g/raw/$g.tsv.gz" "$dir_out/$category/$g/raw/$g.txt.gz" \
			"$dir_raw/$g.gz" "$dir_raw/$g.bgz" "$dir_raw/$g.tsv.gz" "$dir_raw/$g.txt.gz" \
			"$dir_raw/$g.sumstats.gz" "$dir_raw/$g.assoc.gz" "$dir_raw/$g.tsv" "$dir_raw/$g.txt" \
			"$dir_raw/$g.sumstats" "$dir_raw/$g.assoc" "$dir_raw/$g"; do
			[[ -s "$f" ]] && {
				echo "$f"
				return 0
			}
		done
		list_raw_files | while read -r f; do [[ "$(gwas_name_from_file "$f")" == "$g" ]] && {
			echo "$f"
			break
		}; done
	}

	collect_gwas_names() {
		if [[ -n "$gwas_arg" ]]; then
			tr ',' '\n' <<<"$gwas_arg" | sed '/^[[:space:]]*$/d' | sort -u -V
		elif has_step format; then
			list_raw_files | while read -r f; do gwas_name_from_file "$f"; done | sort -u -V
		elif [[ "$step" == "liftover" ]]; then
			list_liftover_names
		elif [[ "$step" == "pgs" ]]; then
			list_pgs_names
		else
			list_names_from_dir "$dir_clean"
		fi
	}

	run_cmds() {
		local list="$1" n rc=0 base joblog
		[[ -s "$list" ]] || {
			log "No command files in $list"
			return 0
		}
		n=$(wc -l <"$list" | tr -d ' ')

		if [[ "$is_bsub" == "TRUE" ]]; then
			log "Submitting $n command files to bsub; per-GWAS logs: $dir_log/$category/<GWAS>/<GWAS>.$step_key.log"
			while read -r cmd; do
				[[ -s "$cmd" ]] || continue
				base=$(basename "$cmd" .cmd)
				bsub -J "gwas_post_${step_key}_$base" -oo "$dir_log/$category/$base/$base.$step_key.bsub.log" -eo "$dir_log/$category/$base/$base.$step_key.bsub.err" \
					"mkdir -p '$dir_log/$category/$base'; bash '$cmd' > '$dir_log/$category/$base/$base.$step_key.log' 2>&1"
			done <"$list"
			return 0
		fi

		log "Running $n command files locally (jobs=$jobs); per-GWAS logs: $dir_log/$category/<GWAS>/<GWAS>.$step_key.log"
		run_one_cmd() {
			local cmd="$1" base log_dir log_file err_file tool_err rc
			[[ -s "$cmd" ]] || return 0
			base=$(basename "$cmd" .cmd)
			log_dir="$dir_log/$category/$base"
			log_file="$log_dir/$base.$step_key.log"
			err_file="$log_dir/$base.$step_key.err"
			mkdir -p "$log_dir"
			echo "[$(date '+%F %T')] START [$step_key] $base" >&2
			rm -f "$err_file"
			if bash "$cmd" >"$log_file" 2>&1; then
				echo "[$(date '+%F %T')] DONE  [$step_key] $base" >&2
			else
				rc=$?
				{
					echo "ERROR: [$step_key] $base failed with exit=$rc"
					echo "ERROR: log=$log_file"
					while IFS= read -r tool_err; do
						[[ -s "$tool_err" ]] || continue
						echo "ERROR: detail=$tool_err"
						grep -Ei 'error|failed|invalid' "$tool_err" || true
					done < <(find "$log_dir" -mindepth 2 -type f -name '*.err' -size +0c 2>/dev/null | sort -V)
				} >"$err_file"
				echo "[$(date '+%F %T')] FAIL  [$step_key] $base exit=$rc log=$err_file" >&2
				return "$rc"
			fi
		}
		export -f run_one_cmd
		export dir_log
		export category
		export step_key
		if command -v parallel >/dev/null 2>&1; then
			joblog=$(mktemp "$list.joblog.XXXXXXXX") || return 1
			log "Parallel job log: $joblog"
			parallel --line-buffer -j "$jobs" --joblog "$joblog" run_one_cmd {} :::: "$list" || rc=$?
			if [[ "$rc" -ne 0 ]]; then
				log "ERROR: [$step_key] one or more command files failed. First failed jobs from $joblog:"
				awk -F '\t' 'NR>1 && $7 != 0 {print "  exit="$7" cmd="$9; n++; if(n>=10) exit}' "$joblog" >&2 || true
				return "$rc"
			fi
		else
			xargs -I{} -P "$jobs" bash -c 'run_one_cmd "$1"' _ {} <"$list" || rc=$?
			return "$rc"
		fi
	}

	gwas_format_check_resources() {
		need_file "$phef"
		need_file "$perf_f"
		command -v bgzip >/dev/null 2>&1 || {
			echo "ERROR: bgzip not found; install htslib" >&2
			exit 1
		}
		command -v tabix >/dev/null 2>&1 || {
			echo "ERROR: tabix not found; install htslib" >&2
			exit 1
		}
		if wants_magma; then
			ensure_magma_resources
			command -v magma >/dev/null 2>&1 || {
				echo "ERROR: magma not found in PATH ($PATH)" >&2
				exit 1
			}
		fi
		if has_step pgs; then
			command -v plink2 >/dev/null 2>&1 || {
				echo "ERROR: plink2 not found in PATH ($PATH)" >&2
				exit 1
			}
		fi
		if [[ "$step" == "all" ]]; then
			cleanup_failed_output_dirs
		fi
		if has_step mplot; then
			need_file "$mplot_r"
			need_file "$plot_f"
			[[ -z "$add_signal" ]] || need_file "$add_signal"
		fi
		if [[ "$thin" == TRUE ]] && { has_step format || has_step thin || has_step mplot || has_step liftover; }; then
			need_file "$thin_r"
			need_file "$phe_r"
			need_file "$hm3_file"
			command -v Rscript >/dev/null 2>&1 || {
				echo "ERROR: Rscript is required for --thin TRUE" >&2
				exit 1
			}
		fi
		if { has_step format && [[ "$hm3_mode" == "TRUE" ]]; } || has_step mplot; then
			need_file "$hm3_file"
		fi
		if [[ "$liftOver" == "TRUE" ]] && has_step liftover; then
			need_file "$chain"
		fi
		if has_step cis && [[ -z "$cis_bed" ]]; then
			echo "ERROR: --cis-bed is required for the cis module" >&2
			exit 2
		fi
		if [[ -n "$cis_bed" ]] && has_step cis; then
			need_file "$cis_bed"
		fi
		if has_step mplot; then
			# Create project-level plot destinations before the potentially long
			# per-GWAS command-generation pass, so background progress is visible at once.
			mkdir -p "$dir_out/mplot" "$dir_out/.project/$category/mplot/flag"
			[[ -z "$cis_bed" ]] || need_file "$cis_bed"
			if [[ -n "$mh_plot_bed" ]]; then
				need_file "$mh_plot_bed"
			else
				need_file "$dir0/files/glist.37.bed"
				need_file "$dir0/files/glist.38.bed"
			fi
		fi
		if { has_step lead || wants_magma; } && [[ "$grch" != auto ]]; then
			need_refgen_clump "$refGen_clump"
		fi
		if has_step lead && [[ "$grch" != auto ]]; then
			need_refgen_cojo "$refGen_cojo"
		fi

	}


	# 🚩 Worker command generation
	# Sourced by format.sh: build one self-contained worker command per GWAS.
	# Uses entry-point configuration and the helpers in format.f.sh.
	# Keep the worker heredoc literal/escaped as written: expansion happens at plan time.

	write_gwas_cmd() {
		local gwas="$1" raw trait_dir gwas_dir source_gwas final cmd awk_snp cis_out qc_prefix clump_dir cojo_dir merged_prefix clump_done cojo_done mh_png mh_meta mh_flag mh_sig mplot_flag_file magma_dir magma_prefix clump_kb
		local pgs_dir pgs_score_file pgs_output pgs_done pgs_meta gwas_pgs_pfile_dir
		local gwas_grch gwas_grch_cache gwas_refGen_clump gwas_refGen_cojo gwas_refGen_id_dir gwas_refGen_keep gwas_gene_loc gwas_mh_plot_bed gwas_hm3_pos detection_input cached_grch rc
		raw=""
		# Only formatting consumes the raw GWAS.  Downstream-only steps use files in
		# the trait's gwas directory, so resolving raw here would rescan the entire
		# project once per trait when no matching raw file exists.
		if has_step format; then
			if [[ -n "$raw_file_arg" ]]; then raw="$raw_file_arg"; else raw=$(raw_file_for_name "$gwas" || true); fi
		fi
		if [[ -n "$dir_clean_arg" ]]; then
			gwas_dir="$dir_clean"
			trait_dir=$(dirname "$gwas_dir")
		else
			trait_dir="$dir_out/$category/$gwas"
			gwas_dir="$trait_dir/gwas"
		fi
		final="$gwas_dir/$gwas.gz"
		if [[ "$liftOver" == "TRUE" ]]; then
			if [[ "$hm3_mode" == FALSE ]]; then source_gwas="$trait_dir/qc/$gwas.source.grch37.gz"; else source_gwas="$gwas_dir/$gwas.hm3.gz"; fi
		else
			source_gwas="$final"
		fi
		# Existing liftOver inputs may still have the old suffix.
		if [[ "$liftOver" == TRUE && "$hm3_mode" == TRUE && ! -s "$source_gwas" ]] &&
			! has_step format && [[ -s "$gwas_dir/$gwas.small.gz" ]]; then
			source_gwas="$gwas_dir/$gwas.small.gz"
		fi
		awk_snp="$gwas_dir/$gwas.awk.snp"
		cis_out="$gwas_dir/$gwas.cis.gz"
		qc_prefix="$trait_dir/qc/$gwas"
		gwas_grch_cache="${qc_prefix}.grch"
		clump_dir="$gwas_dir/clump"
		cojo_dir="$gwas_dir/cojo"
		merged_prefix="$gwas_dir/$gwas"
		clump_done="$gwas_dir/$gwas.clump.done"
		cojo_done="$gwas_dir/$gwas.cojo.done"
		mh_png="$dir_out/mplot/$gwas.png"
		mh_meta="$dir_out/.project/$category/mplot/$gwas.state"
		mh_flag="$dir_out/.project/$category/mplot/flag/$gwas.tsv"
		mh_sig="$dir_out/mplot/$gwas.sig.txt"
		mplot_flag_file="$dir_out/mplot/0flag.tsv"
		pgs_dir="$trait_dir/pgs"
		pgs_score_file="${merged_prefix}.jma.cojo"
		pgs_output="$pgs_dir/$gwas.pgs.gz"
		pgs_done="$pgs_dir/$gwas.pgs.done"
		pgs_meta="$pgs_dir/$gwas.pgs.meta.tsv"
		if [[ -n "$dir_magma" ]]; then magma_dir="$dir_magma"; else magma_dir="$trait_dir/magma"; fi
		magma_prefix="$magma_dir/$gwas"
		clump_kb=$(((lead_window + 999) / 1000))
		cmd="$dir_cmd/$gwas.cmd"
		mkdir -p "$gwas_dir" "$trait_dir/qc" "$(dirname "$mh_meta")"

		# Keep the standalone PGS fast path intentionally minimal: the two non-empty
		# result/marker files alone define completion.  Skip build detection, pfile
		# validation, gzip tests, metadata reads, and timestamp comparisons.
		if [[ "$replace" != "TRUE" && "$step" == pgs ]] &&
			pgs_output_complete "$pgs_output" "$pgs_done"; then
			prune_pgs_dir "$pgs_dir" "$gwas"
			log "SKIP completed PGS: $gwas"
			rm -f "$cmd"
			return 0
		fi

		# A completed lead result does not need GWAS-build detection or reference
		# validation.  Require independent clump/COJO markers plus their artifacts,
		# except for the explicit no-significant-row case.
		if [[ "$replace" != "TRUE" && "$step" == lead ]] && gzip_ok "$final" &&
			gwas_lead_marker_matches "$clump_done" "$gwas" clump "$p_lead" "$lead_window" "$chrs" &&
			gwas_lead_marker_matches "$cojo_done" "$gwas" cojo "$p_lead" "$lead_window" "$chrs" &&
			[[ -s "$awk_snp" && ! -d "$clump_dir" && ! -d "$cojo_dir" ]] &&
			{ { [[ -s "$clump_done" && -s "${merged_prefix}.clumps" && -s "$cojo_done" ]] &&
				{ [[ -s "${merged_prefix}.jma.cojo" || -s "${merged_prefix}.ldr.cojo" ]] ||
					awk -F '\t' 'NR==2&&$6=="no_snps_selected"{ok=1}END{exit !ok}' "$cojo_done"; }; } ||
				{ [[ -s "$clump_done" && -s "$cojo_done" ]] && awk 'NR>1{exit 1}' "$awk_snp"; }; }; then
			log "SKIP completed lead GWAS: $gwas"
			rm -f "$cmd"
			return 0
		fi

		# Standalone thin does not need source indexing or any plotting setup.
		# Resolve builds in workers so an uncached input cannot serialize planning.
		if [[ "$step" == thin ]]; then
			cat >"$cmd" <<THIN_CMD
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
FINAL=$(q "$final")
THIN_OUT=$(q "${final%.gz}.thin.gz")
GRCH=$(q "$grch")
GRCH_CACHE=$(q "$gwas_grch_cache")
HM3=$(q "$hm3_file")
HM3_POS_TEMPLATE=$(q "$hm3_pos")
THIN_R=$(q "$thin_r")
PHE_R=$(q "$phe_r")
THIN_CHR_MAX=$(q "$thin_chr_max")
THIN_MODE=$(q "$thin")
REPLACE=$(q "$replace")
MH_META=$(q "$mh_meta")
DO_STEP=thin
source $(q "$phef")
source $(q "$perf_f") --worker
gwas_post_log(){ printf '[thin] %s\n' "\$*" >&2; }
[[ "\$THIN_MODE" == TRUE ]] || exit 0
[[ -s "\$FINAL" ]] || { echo "ERROR: missing GWAS: \$FINAL" >&2; exit 1; }
if [[ "\$GRCH" == auto ]]; then
  cached=""
  if [[ -s "\$GRCH_CACHE" && "\$GRCH_CACHE" -nt "\$FINAL" ]]; then
    cached=\$(awk 'NR==1 && (\$1==37 || \$1==38){print \$1}' "\$GRCH_CACHE")
  fi
  if [[ -n "\$cached" ]]; then
    GRCH="\$cached"
  else
    check_GRCH "\$FINAL" >&2
    GRCH="\$CHECK_GRCH_RESULT"
    printf '%s\n' "\$GRCH" > "\$GRCH_CACHE.tmp.\$\$"
    mv -f "\$GRCH_CACHE.tmp.\$\$" "\$GRCH_CACHE"
  fi
fi
HM3_POS=\${HM3_POS_TEMPLATE//\{grch\}/\$GRCH}
gwas_post_thin
THIN_CMD
			# Workers run via bash; do not require chmod on mounted project drives.
			printf '%s\n' "$cmd"
			return 0
		fi

		gwas_grch="$grch"
		# With no fixed build (omitted or --grch auto), resolve every GWAS independently
		# from the 39 sentinel rsIDs.  Use RAW before format, SMALL before a standalone
		# liftOver, and the standardized FINAL for downstream-only requests.
		if [[ "$gwas_grch" == auto ]]; then
			if has_step format; then
				detection_input="$raw"
				[[ -s "$detection_input" ]] || {
					echo "ERROR: missing raw GWAS for GRCh detection: $gwas" >&2
					exit 1
				}
			elif [[ "$step" == liftover ]]; then
				detection_input="$source_gwas"
				[[ -s "$detection_input" ]] || {
					echo "ERROR: missing source GWAS for GRCh detection: $gwas" >&2
					exit 1
				}
			else
				detection_input="$final"
				[[ -s "$detection_input" ]] || {
					echo "ERROR: missing standardized GWAS for GRCh detection: $gwas" >&2
					exit 1
				}
			fi
			cached_grch=""
			if [[ "$liftOver" != TRUE && -s "$gwas_grch_cache" && "$gwas_grch_cache" -nt "$detection_input" ]]; then
				cached_grch=$(awk 'NR==1 && ($1==37 || $1==38){print $1}' "$gwas_grch_cache")
			fi
			if [[ -n "$cached_grch" ]]; then
				gwas_grch="$cached_grch"
				log "Reuse cached GRCh$gwas_grch: $gwas"
			elif check_GRCH "$detection_input" >&2; then
				gwas_grch="$CHECK_GRCH_RESULT"
				printf '%s\n' "$gwas_grch" >"${gwas_grch_cache}.tmp.$$"
				mv -f "${gwas_grch_cache}.tmp.$$" "$gwas_grch_cache"
			else
				rc=$?
				echo "ERROR: automatic GRCh detection failed for $gwas; the input must contain usable rsIDs from ${CHECK_GRCH_SNP_LIST:-$dir0/data/ukb/phe/common/snp.lst}, or specify --grch 37/38." >&2
				return "$rc"
			fi
		fi
		gwas_hm3_pos=${hm3_pos//\{grch\}/$gwas_grch}
		gwas_pgs_pfile_dir="${pgs_pfile_dir:-/mnt/f/gen/ukb/${gwas_grch}/imp}"
		if has_step pgs; then
			need_pgs_pfiles "$gwas_pgs_pfile_dir"
		fi
		gwas_gene_loc="$gene_loc"
		if [[ -z "$gwas_gene_loc" ]]; then
			[[ "$gwas_grch" == 37 ]] && gwas_gene_loc="$dir0/files/NCBI.37.gene.loc" || gwas_gene_loc="$dir0/files/NCBI.38.gene.loc"
		fi
		gwas_mh_plot_bed=""
		gwas_mh_plot_bed="$mh_plot_bed"
		[[ -n "$gwas_mh_plot_bed" ]] || gwas_mh_plot_bed="$dir0/files/glist.${gwas_grch}.bed"
		if has_step format && [[ "$liftOver" == TRUE && "$gwas_grch" != 37 ]]; then
			echo "ERROR: --liftover TRUE uses the GRCh37-to-GRCh38 chain, but $gwas was detected/configured as GRCh$gwas_grch" >&2
			exit 1
		fi
		gwas_refGen_clump="${refGen_clump:-/mnt/f/gen/1kg/${gwas_grch}/pfile/}"
		gwas_refGen_cojo="${refGen_cojo:-/mnt/f/gen/1kg/${gwas_grch}/pfile/${refGen_pop}/}"
		gwas_refGen_id_dir="${refGen_id_dir:-/mnt/f/gen/1kg/${gwas_grch}/id}"
		gwas_refGen_keep=""
		if has_step lead || wants_magma || { has_step format && [[ "$fill_eaf" == TRUE ]]; }; then
			if [[ "$refGen_pop" != ALL ]]; then
				gwas_refGen_keep="${gwas_refGen_id_dir%/}/${refGen_pop}.id.2col"
				if ! awk 'BEGIN{FS="[ \t]+"}{a=$1;b=$2;gsub(/\r/,"",a);gsub(/\r/,"",b);if(NF!=2||a!=b)bad=1;n++}END{exit bad||n==0?2:0}' "$gwas_refGen_keep"; then
					echo "ERROR: invalid PLINK two-column keep file: $gwas_refGen_keep" >&2
					exit 1
				fi
			fi
			need_refgen_clump "$gwas_refGen_clump"
		fi
		if has_step lead; then
			need_refgen_cojo "$gwas_refGen_cojo"
		fi

		if [[ "$replace" != "TRUE" ]]; then
			case "$step" in
				mplot | thin,mplot | mplot,thin)
					if mplot_output_complete "$mh_png" "$final" "$magma_prefix.genes.out" "$mh_meta" "$gwas_grch" "$mh_flag" "$mplot_flag_file" "$mh_sig" "${merged_prefix}.jma.cojo"; then
						log "SKIP completed Manhattan plot: $gwas panel=$add_panel"
						rm -f "$cmd"
						return 0
					fi
					;;
				magma)
					if magma_output_complete "$magma_dir" "$magma_prefix"; then
						prune_magma_dir "$magma_dir"
						log "SKIP completed MAGMA: $gwas"
						rm -f "$cmd"
						return 0
					fi
					;;
				cis)
					if cis_output_complete "$cis_out"; then
						log "SKIP completed indexed cis GWAS: $gwas"
						rm -f "$cmd"
						return 0
					fi
					;;

			esac
		fi

		cat >"$cmd" <<CMD_TOP
#!/bin/bash
set -euo pipefail
export LC_ALL=C
export PATH=$(q "$dir0/software/bin"):\$PATH

gwas_post_need_file(){
  [[ -s "\$1" ]] || { echo "ERROR: missing or empty file: \$1" >&2; exit 1; }
}
gwas_post_log(){
  printf '[%s] %s\n' "\$(date '+%F %T')" "\$*"
}
gwas_post_zcat(){
  case "\$1" in *.gz|*.bgz) gzip -cd -- "\$1";; *) cat -- "\$1";; esac
}
gwas_post_has_data_rows(){
  [[ -s "\$1" ]] && awk 'NR>1{found=1; exit} END{exit found ? 0 : 1}' "\$1"
}

GWAS=$(q "$gwas")
GWAS_DIR=$(q "$gwas_dir")
if [[ -n "\${GWAS_POST_TMP_BASE:-}" ]]; then
  GWAS_POST_TMP_ROOT="\${GWAS_POST_TMP_BASE%/}/\$GWAS"
else
  GWAS_POST_TMP_ROOT="\$GWAS_DIR/.tmp"
fi
mkdir -p "\$GWAS_POST_TMP_ROOT"
# Isolate temporary files per process and remove them on exit.
for stale_tmp in "\$GWAS_POST_TMP_ROOT"/run.*; do
  [[ -d "\$stale_tmp" ]] || continue
  stale_pid=\${stale_tmp##*.}
  if [[ ! "\$stale_pid" =~ ^[1-9][0-9]*\$ ]] || ! kill -0 "\$stale_pid" 2>/dev/null; then
    rm -rf -- "\$stale_tmp"
  fi
done
GWAS_POST_TMP="\$GWAS_POST_TMP_ROOT/run.\$\$"
mkdir -p "\$GWAS_POST_TMP"
export TMPDIR="\$GWAS_POST_TMP"
gwas_post_cleanup_tmp(){
  rm -rf -- "\$GWAS_POST_TMP"
  rmdir -- "\$GWAS_POST_TMP_ROOT" 2>/dev/null || true
}
trap gwas_post_cleanup_tmp EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
RAW=$(q "$raw")
# SMALL and DO_STEP=small are the unchanged format.f.sh API.
# The project uses --hm3/HM3_MODE and std_hm3 for the filtering choice.
SMALL=$(q "$source_gwas")
FINAL=$(q "$final")
AWK_SNP=$(q "$awk_snp")
CIS_OUT=$(q "$cis_out")
QC_PREFIX=$(q "$qc_prefix")
CLUMP=$(q "$clump_dir")
COJO=$(q "$cojo_dir")
MERGED=$(q "$merged_prefix")
CLUMP_DONE=$(q "$clump_done")
COJO_DONE=$(q "$cojo_done")
REF_BAD_DIR=$(q "$dir_out/.project/$category/qc/ref_bad")

PHEF=$(q "$phef")
PERF_F=$(q "$perf_f")
PLOT_F=$(q "$plot_f")
MPLOT_R=$(q "$mplot_r")
PLOT_METHOD=self
ADD_PANEL=$(q "$add_panel")
PLOT_WIDTH=$(q "$plot_width")
PLOT_HEIGHT=$(q "$plot_height")
PLOT_RES=$(q "$plot_res")
WRITE_SIG=$(q "$write_sig")
ADD_SIGNAL=$(q "$add_signal")
SIGNAL_MATCH_COL=$(q "$signal_match_col")
SIGNAL_MATCH_VALUE=$(q "$signal_match_value")
SIGNAL_LOCUS_POS=$(q "$signal_locus_pos")
SIGNAL_DISPLAY_COL=$(q "$signal_display_col")
PGS_STEPS=$(q "$step")
LEAD_REF_CACHE=$(q "$dir_out/.project/$category/lead_reference")
MH_PLOT_BED=$(q "$gwas_mh_plot_bed")
HM3=$(q "$hm3_file")
HM3_POS=$(q "$gwas_hm3_pos")
P_HM3=$(q "$p_hm3")
HM3_MODE=$(q "$hm3_mode")
THIN_MODE=$(q "$thin")
THIN_CHR_MAX=$(q "$thin_chr_max")
THIN_R=$(q "$thin_r")
PHE_R=$(q "$phe_r")
THIN_OUT=$(q "${final%.gz}.thin.gz")
MH_PNG=$(q "$mh_png")
MH_META=$(q "$mh_meta")
MH_FLAG=$(q "$mh_flag")
MH_SIG=$(q "$mh_sig")
COJO_FILE=$(q "${merged_prefix}.jma.cojo")
MPLOT_FLAG_FILE=$(q "$mplot_flag_file")
MAGMA_DIR=$(q "$magma_dir")
MAGMA_PREFIX=$(q "$magma_prefix")
MAGMA_ANNOT_CACHE=$(q "$magma_annot_cache")
DELETE_RAW_AFTER_SUCCESS=$(q "$delete_raw")
# The shared formatter can delete raw immediately after formatting. Defer it
# until all requested modules in this worker have returned successfully.
DELETE_RAW=FALSE

DO_STEP=$(q "$step")
REPLACE=$(q "$replace")
FILL_EAF=$(q "$fill_eaf")
FILL_N=$(q "$fill_n")
N_TOTAL=$(q "$n_total")
H2_HELPER=$(q "${SCRIPT_PATH%/*}/f/format.py")
LIFTOVER_HELPER=$(q "${SCRIPT_PATH%/*}/f/format.py")
H2_SEX=$(q "$h2_sex")
H2_REF_LD=$(q "$h2_ref_ld")
H2_W_LD=$(q "$h2_w_ld")
H2_MERGE_ALLELES=$(q "$h2_merge_alleles")
H2_PYTHON=$(q "$h2_python")
H2_CONDA_ENV=$(q "$h2_conda_env")
H2_REQUESTED=$(q "$step")
DO_LIFTOVER=$(q "$liftOver")
CHAIN=$(q "$chain")
LIFTOVER_BIN=$(q "$liftover_bin")

CIS_BED=$(q "$cis_bed")
CIS_FLANK=$(q "$cis_flank")

P_LEAD=$(q "$p_lead")
LEAD_WINDOW=$(q "$lead_window")
CLUMP_KB=$(q "$clump_kb")
CHRS=$(q "$chrs")
REFGEN_CLUMP=$(q "$gwas_refGen_clump")
REFGEN_POP=$(q "$refGen_pop")
REFGEN_KEEP=$(q "$gwas_refGen_keep")
REFGEN_COJO=$(q "$gwas_refGen_cojo")

GRCH=$(q "$gwas_grch")
MAGMA_REF=$(q "$magma_ref")
GENE_LOC=$(q "$gwas_gene_loc")
SYNONYMS=$(q "$synonyms")
MAGMA_WINDOW=$(q "$magma_window")
MAGMA_N=$(q "$magma_N")
GWAS_N=$(q "$gwas_N")
PGS_DIR=$(q "$pgs_dir")
PGS_SCORE_FILE=$(q "$pgs_score_file")
PGS_OUTPUT=$(q "$pgs_output")
PGS_DONE=$(q "$pgs_done")
PGS_META=$(q "$pgs_meta")
PGS_PFILE_DIR=$(q "${gwas_pgs_pfile_dir%/}")
PGS_THREADS=$(q "$pgs_threads")

# Do not delete trait-wide error files here: another module may be using the
# same trait concurrently.  Each coordinator now owns a module-specific error.

if [[ ",\$DO_STEP," == *,format,* || "\$DO_STEP" == "all" ]]; then
  if [[ ! -s "\$RAW" ]]; then
    echo "ERROR: missing/empty raw GWAS for \$GWAS: \$RAW" >&2
    echo "ERROR: removing incomplete outputs for \$GWAS before re-run." >&2
    find "\$GWAS_DIR" -mindepth 1 -depth ! -path "\$GWAS_DIR/\$GWAS.log" -delete 2>/dev/null || true
    rm -f -- "\${QC_PREFIX}"* "\$MH_PNG" "\$MH_META" "\${MH_PNG}.meta.tsv" "\$MH_SIG"
    exit 1
  fi
fi

source "\$PERF_F" --base

# Format every raw data row into the standard 11-column schema.  Unlike
# std_hm3(), this intentionally performs no variant filtering.
std_format(){
  local src="\$1" out="\$2" header_out="\$3" tmp input_fs
  gwas_post_need_file "\$src"
  gwas_post_log "Format full GWAS: \$src -> \$out"
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); EAF_col=\$(gwas_clean_col "\${EAF_col:-}")
  N_col=\$(gwas_clean_col "\${N_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  LOG10P_col=\$(gwas_clean_col "\${LOG10P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  input_fs=\$(gwas_clean_detect_fs "\$src"); tmp="\${out}.tmp.\$\$"
  gwas_post_zcat "\$src" | awk -v FS="\$input_fs" -v OFS='\t' \
    -v snp_col="\$SNP_col" -v chr_col="\$CHR_col" -v pos_col="\$POS_col" \
    -v ea_col="\$EA_col" -v nea_col="\$NEA_col" -v eaf_col="\$EAF_col" -v n_col="\$N_col" \
    -v beta_col="\$BETA_col" -v se_col="\$SE_col" -v p_col="\$P_col" -v logp_col="\$LOG10P_col" '
    function get(c, x){x=(c>0 ? \$c : "");gsub(/^[[:space:]]+|[[:space:]]+\$/,"",x);return x}
    function val(c, x){x=get(c); return x=="" ? "NA" : x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function normchr(x){gsub(/^chr/,"",x); if(x=="X")x="23"; if(x=="Y")x="24"; if(x=="MT"||x=="M")x="25"; return x}
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P"; next}
    {chr=normchr(get(chr_col)); pos=get(pos_col); snp=get(snp_col); ea=get(ea_col); nea=get(nea_col)
     if(snp==""||snp=="NA"||snp=="."){snp=chr":"pos; if(ea!="")snp=snp":"ea; if(nea!="")snp=snp":"nea}
     p=(p_col>0 ? get(p_col) : ""); lp=(logp_col>0 ? get(logp_col) : "")
     if(p=="" && isnum(lp))p=10^(-lp); if(lp=="" && isnum(p) && p>0)lp=-log(p)/log(10)
     if(snp=="")snp="NA"; if(chr=="")chr="NA"; if(pos=="")pos="NA"; if(ea=="")ea="NA"; if(nea=="")nea="NA"
     if(p=="")p="NA"; if(lp=="")lp="NA"
     print snp,chr,pos,ea,nea,val(eaf_col),val(n_col),val(beta_col),val(se_col),p,lp}' | \
    { IFS= read -r format_header; printf '%s\n' "\$format_header"; sort -t \$'\t' -k2,2n -k3,3n; } | \
    gwas_clean_compress > "\$tmp"
  gzip -t "\$tmp"; mv -f "\$tmp" "\$out"
  gwas_post_zcat "\$out" | awk -v g="\$GWAS" 'NR==1{next} END{print g"\t"NR-1}' > "\${QC_PREFIX}.format.nrow.tsv"
}

gwas_post_validate_format_columns(){
  local src="\$1" header_out="\${QC_PREFIX}.header.required.txt"
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required raw GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  gwas_post_log "format column map: SNP=\$SNP_col CHR=\$CHR_col POS=\$POS_col EA=\$EA_col NEA=\$NEA_col BETA=\$BETA_col SE=\$SE_col P=\$P_col"
}

gwas_post_fill_missing_fields(){
  local target="\$1" tmp eaf_tmp
  local -a eaf_keep_args=()
  [[ -s "\$target" ]] || { echo "ERROR: formatted GWAS is missing: \$target" >&2; return 1; }

  if [[ -n "\$FILL_N" ]]; then
    tmp="\${target}.fill_n.tmp.\$\$"
    gwas_post_log "fill missing/invalid N with \$FILL_N: \$target"
    gwas_post_zcat "\$target" | awk -v FS='\t' -v OFS='\t' -v fill_n="\$FILL_N" -v audit="\${QC_PREFIX}.fill_n.tsv" '
      function isnum(x){return x ~ /^[0-9]+([.][0-9]+)?\$/ && x+0>0}
      NR==1{for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}if(!("N" in c)){print "ERROR: standardized GWAS lacks N" > "/dev/stderr";exit 2}print;next}
      {if(isnum(\$(c["N"])))kept++;else{\$(c["N"])=fill_n;filled++}print}
      END{print "STATUS\tN" > audit;print "existing\t" kept+0 >> audit;print "filled\t" filled+0 >> audit}
    ' | gzip -c > "\$tmp"
    gzip -t "\$tmp"; mv -f "\$tmp" "\$target"
  fi

  if [[ "\$FILL_EAF" == TRUE ]]; then
    declare -F match_EAF >/dev/null 2>&1 || gwas_clean_load_phef
    eaf_tmp="\${target}.fill_eaf.tmp.\$\$.gz"
    gwas_post_log "fill missing/invalid EAF from 1KG GRCh\$GRCH: \$target"
    [[ -z "\$REFGEN_KEEP" ]] || eaf_keep_args=(--keep "\$REFGEN_KEEP")
    match_EAF --reference "\$REFGEN_CLUMP" "\${eaf_keep_args[@]}" --output "\$eaf_tmp" \
      --audit "\${QC_PREFIX}.fill_eaf.tsv" "\$target"
    gzip -t "\$eaf_tmp"; mv -f "\$eaf_tmp" "\$target"
  fi
}

gwas_post_prune_magma_dir(){
  [[ -d "\$MAGMA_DIR" ]] || return 0
  find "\$MAGMA_DIR" -mindepth 1 -maxdepth 1 -type f \
    ! -name '*.genes.out' ! -name '*.genes.raw' ! -name '*.log' \
    ! -name 'magma.meta.tsv' ! -name 'magma.done' -delete
}

gwas_post_magma_annotation(){
  local snploc="\${1:-}" cache_key window_tag resource_tag cache_dir annot meta coordinate_tag
  local lock_file lock_fd tmp_prefix annot_tmp meta_tmp nloc snploc_hash
  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock not found; required for the shared MAGMA annotation cache" >&2; exit 1; }

  window_tag=\$(printf '%s' "\$MAGMA_WINDOW" | tr -c 'A-Za-z0-9._-' '_')
  resource_tag=\$(printf '%s\n%s\n' "\$MAGMA_REF" "\$GENE_LOC" | sha256sum | awk '{print substr(\$1,1,16)}')
  # Include the actual annotation input: a project can contain different SNP
  # sets, and v2 may contain coordinate IDs incompatible with the MAGMA BIM.
  if [[ -n "\$snploc" ]]; then
    coordinate_tag=\$(sha256sum "\$snploc" | awk '{print \$1}')
  else
    coordinate_tag=\$(stat -c '%n:%s:%y' "\$FINAL" | sha256sum | awk '{print \$1}')
  fi
  cache_key="v3.GRCh\${GRCH}.window_\${window_tag}.magma_\${resource_tag}.\${coordinate_tag}"
  cache_dir="\$MAGMA_ANNOT_CACHE/v3/GRCh\$GRCH/window_\$window_tag/magma_\$resource_tag/\$coordinate_tag"
  annot="\$cache_dir/genes.annot"
  meta="\$cache_dir/annotation.meta.tsv"
  mkdir -p "\$(dirname "\$cache_dir")"

  lock_file="\${cache_dir}.lock"
  exec {lock_fd}> "\$lock_file"
  flock "\$lock_fd"
  if [[ ! -s "\$annot" || ! -s "\$meta" ]]; then
    mkdir -p "\$cache_dir"
    if [[ -z "\$snploc" ]]; then
      snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"
      gwas_post_log "Prepare MAGMA annotation coordinates from \$FINAL"
      gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
        NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
        {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
         if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
        sort -T "\$GWAS_POST_TMP" -k1,1 -k2,2n -k3,3n -u > "\$snploc"
    fi
    gwas_post_need_file "\$snploc"
    nloc=\$(wc -l < "\$snploc" | tr -d ' ')
    (( nloc > 1000 )) || { echo "ERROR: too few MAGMA annotation SNPs: \$nloc" >&2; return 1; }
    snploc_hash=\$(sha256sum "\$snploc" | awk '{print \$1}')
    tmp_prefix="\$GWAS_POST_TMP/annotation.\$\$"
    rm -f -- "\${tmp_prefix}"*
    gwas_post_log "Build shared MAGMA annotation from GWAS rsID coordinates: \$annot"
    if [[ "\$MAGMA_WINDOW" == "0,0" || "\$MAGMA_WINDOW" == 0 ]]; then
      magma --annotate --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    else
      magma --annotate window="\$MAGMA_WINDOW" --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    fi
    gwas_post_need_file "\$tmp_prefix.genes.annot"
    annot_tmp="\${annot}.tmp.\$\$"
    mv -f "\$tmp_prefix.genes.annot" "\$annot_tmp"
    mv -f "\$annot_tmp" "\$annot"
    meta_tmp="\${meta}.tmp.\$\$"
    {
      printf 'key\tvalue\n'
      printf 'grch\t%s\nwindow_kb\t%s\ngene_loc\t%s\nld_reference\t%s\nsource_gwas\t%s\nsnp_loc_n\t%s\nsnp_loc_sha256\t%s\n' \
        "\$GRCH" "\$MAGMA_WINDOW" "\$GENE_LOC" "\$MAGMA_REF" "\$GWAS" "\$nloc" "\$snploc_hash"
    } > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$meta"
  else
    gwas_post_log "Reuse shared MAGMA annotation: \$annot"
  fi
  flock -u "\$lock_fd"
  exec {lock_fd}>&-

  MAGMA_ANNOT="\$annot"
  MAGMA_ANNOT_KEY="\$cache_key"
  MAGMA_SNPLOC_N=\$(awk -F '\t' '\$1=="snp_loc_n"{print \$2;exit}' "\$meta")
  MAGMA_SNPLOC_HASH=\$(awk -F '\t' '\$1=="snp_loc_sha256"{print \$2;exit}' "\$meta")
}

gwas_post_magma(){
  [[ ",\$DO_STEP," == *,magma,* ]] || return 0
  [[ -s "\$FINAL" ]] || { echo "ERROR: missing or empty file: \$FINAL" >&2; exit 1; }

  if [[ "\$GRCH" == "38" ]]; then
    echo " [\$GWAS] GRCh build 38" >&2
  else
    echo "[\$GWAS] GRCh build 37" >&2
  fi
  echo "   gene location : \$GENE_LOC" >&2
  echo "   LD reference  : \$MAGMA_REF" >&2

  if [[ "\$REPLACE" != TRUE && -s "\$MAGMA_DIR/magma.done" && -s "\$MAGMA_PREFIX.genes.out" && -s "\$MAGMA_PREFIX.genes.raw" ]]; then
    gwas_post_prune_magma_dir
    echo "[\$(date '+%F %T')] MAGMA exists: \$MAGMA_PREFIX.genes.out (GRCh\$GRCH)" >&2
    return 0
  fi
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH" >&2; exit 1; }
  for ext in bed bim fam; do [[ -s "\$MAGMA_REF.\$ext" ]] || { echo "ERROR: missing or empty file: \$MAGMA_REF.\$ext" >&2; exit 1; }; done
  [[ -s "\$GENE_LOC" ]] || { echo "ERROR: missing or empty file: \$GENE_LOC" >&2; exit 1; }
  [[ -s "\$SYNONYMS" ]] || { echo "ERROR: missing or empty file: \$SYNONYMS" >&2; exit 1; }
  mkdir -p "\$MAGMA_DIR"
  # Never leave an old completion marker/meta file behind during a rerun.
  rm -f "\$MAGMA_DIR/magma.done" "\$MAGMA_DIR/magma.meta.tsv"

  if [[ -z "\$MAGMA_N" && -s "\$GWAS_DIR/\$GWAS.magma.N" ]]; then
    MAGMA_N=\$(awk 'NF{print \$1; exit}' "\$GWAS_DIR/\$GWAS.magma.N")
  fi
  if [[ -z "\$MAGMA_N" && ! "\$GWAS_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid default gwas_N: \$GWAS_N" >&2; exit 1
  fi
  if [[ -n "\$MAGMA_N" && ! "\$MAGMA_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid MAGMA sample size for \$GWAS: \$MAGMA_N" >&2; exit 1
  fi

  header=\$(set +o pipefail; gwas_clean_zcat "\$FINAL" | head -1)
  for col in SNP CHR POS P; do
    awk -F '\t' -v c="\$col" '{for(i=1;i<=NF;i++)if(\$i==c)exit 0;exit 1}' <<< "\$header" ||
      { echo "ERROR: \$FINAL lacks required column \$col" >&2; exit 1; }
  done
  snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"; pval="\$GWAS_POST_TMP/\$GWAS.pval"
  gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
    NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
    {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
     if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
    sort -k1,1 -k2,2n -k3,3n -u > "\$snploc"

  has_n=FALSE
  awk -F '\t' '{for(i=1;i<=NF;i++)if(\$i=="N")exit 0;exit 1}' <<< "\$header" && has_n=TRUE || true
  if [[ -n "\$MAGMA_N" ]]; then
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  elif [[ "\$has_n" == TRUE ]]; then
    { printf 'SNP\tP\tN\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);n=\$(c["N"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1&&n+0>=50)print s,p,n}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg='ncol=N'
  else
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  fi
  nloc=\$(wc -l < "\$snploc" | tr -d ' '); npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  if (( npval <= 1000 )) && [[ "\$narg" == ncol=N ]]; then
    echo "WARNING: no usable per-SNP N for \$GWAS; falling back to gwas_N=\$GWAS_N" >&2
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"; npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  fi
  (( nloc > 1000 && npval > 1000 )) || { echo "ERROR: too few MAGMA SNPs: snploc=\$nloc pval=\$npval" >&2; exit 1; }
  gwas_post_magma_annotation "\$snploc"
  magma --bfile "\$MAGMA_REF" synonyms="\$SYNONYMS" --pval "\$pval" "\$narg" --gene-annot "\$MAGMA_ANNOT" --out "\$MAGMA_PREFIX"
  [[ -s "\$MAGMA_PREFIX.genes.out" ]] || { echo "ERROR: MAGMA output missing: \$MAGMA_PREFIX.genes.out" >&2; exit 1; }
  [[ -s "\$MAGMA_PREFIX.genes.raw" ]] || { echo "ERROR: MAGMA intermediate gene result missing: \$MAGMA_PREFIX.genes.raw" >&2; exit 1; }
  meta_tmp="\$MAGMA_DIR/magma.meta.tsv.tmp.\$\$"
  { printf 'key\tvalue\n'; printf 'gwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\nannotation_cache\t%s\nannotation_key\t%s\nsnp_loc_sha256\t%s\n' \
      "\$GWAS" "\$GRCH" "\$GENE_LOC" "\$MAGMA_REF" "\$SYNONYMS" "\$MAGMA_WINDOW" "\$nloc" "\$npval" "\$MAGMA_ANNOT" "\$MAGMA_ANNOT_KEY" "\$MAGMA_SNPLOC_HASH"; } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$MAGMA_DIR/magma.meta.tsv"
  gwas_post_prune_magma_dir
  date '+%F %T' > "\$MAGMA_DIR/magma.done"
  echo "[\$(date '+%F %T')] MAGMA done: \$MAGMA_PREFIX.genes.out" >&2
}

gwas_post_mplot(){
  [[ ",\$DO_STEP," == *,mplot,* || "\$DO_STEP" == "all" ]] || return 0
  gwas_post_need_file "\$FINAL"
  if [[ "\$REPLACE" != "TRUE" && -s "\$MH_PNG" ]]; then
    gwas_post_log "Manhattan plot exists: \$MH_PNG"
    return 0
  fi
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found; required for Manhattan plotting" >&2; exit 1; }
  mkdir -p "\$(dirname "\$MH_PNG")"
  plot_input="\$FINAL"
  if [[ "\$HM3_MODE" == "FALSE" ]]; then
    plot_input="\${QC_PREFIX}.mplot.hm3.gz"
    gwas_post_log "Subset Manhattan input to HM3 or P < 0.001: \$plot_input"
    gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v hm3_file="\$HM3" -v pthr='0.001' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{while((getline x<hm3_file)>0){split(x,a,/[ \t]+/); if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1} close(hm3_file)}
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i; print; next}
      (("SNP" in c)&&\$(c["SNP"]) in hm3) || (("P" in c)&&isnum(\$(c["P"]))&&\$(c["P"])+0<pthr){print}' | gwas_clean_compress > "\$plot_input"
    gzip -t "\$plot_input"
  fi
  gwas_post_log "Manhattan plot: \$plot_input -> \$MH_PNG"
  Rscript - "\$plot_input" "\$MH_PNG" "\$GWAS" "\$PLOT_F" "\$CIS_BED" "\$MH_PLOT_BED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
output <- args[[2]]
gwas <- args[[3]]
plot_f <- args[[4]]
cis_bed <- if (nzchar(args[[5]])) args[[5]] else NULL
mh_plot_bed <- if (nzchar(args[[6]])) args[[6]] else NULL
source(plot_f)

cis_gene_arg <- NULL
if (!is.null(cis_bed) && file.exists(cis_bed)) {
  cis_rows <- tryCatch(utils::read.table(cis_bed, header=FALSE, comment.char="#", stringsAsFactors=FALSE),
                       error=function(e) NULL)
  if (!is.null(cis_rows) && ncol(cis_rows) >= 4 && any(as.character(cis_rows[[4]]) == gwas))
    cis_gene_arg <- gwas
}

png_args <- list(filename = output, width = 13.333, height = 7.5, units = "in", res = 300)
if (capabilities("cairo")) png_args[["type"]] <- "cairo"
do.call(grDevices::png, png_args)
on.exit(grDevices::dev.off(), add = TRUE)
graphics::par(mar = c(5, 5, 3, 1) + 0.1)
print(mh_plot(input, col=c("gray", "darkgray"), cis_gene=cis_gene_arg, cis_bed=cis_bed, mh_plot_bed=mh_plot_bed,
              cis_color="red", other_top_color="green",
              other_top_max_per_chr=2, locus_size=1e6, main=gwas))
RSCRIPT
  [[ "\$plot_input" == "\$FINAL" ]] || rm -f "\$plot_input"
  [[ -s "\$MH_PNG" ]] || { echo "ERROR: Manhattan plot was not created: \$MH_PNG" >&2; exit 1; }
}

gwas_post_pgs_current(){
  [[ -s "\$PGS_OUTPUT" && -s "\$PGS_DONE" ]]
}

gwas_post_prune_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    rm -f -- "\$f"
  done
}

gwas_post_prune_empty_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    case "\$f" in *.log) continue;; esac
    rm -f -- "\$f"
  done
}

gwas_post_pgs(){
  [[ ",\$PGS_STEPS," == *,pgs,* || "\$PGS_STEPS" == "all" ]] || return 0
  if [[ "\$REPLACE" != "TRUE" ]] && gwas_post_pgs_current; then
    if [[ -s "\$PGS_META" ]] && awk -F '\t' '\$1=="matched_variants"&&\$2=="none"{found=1} END{exit !found}' "\$PGS_META"; then
      gwas_post_prune_empty_pgs_dir
    else
      gwas_post_prune_pgs_dir
    fi
    gwas_post_log "PGS exists and is current: \$PGS_OUTPUT"
    return 0
  fi
  gwas_post_need_file "\$PGS_SCORE_FILE"
  gwas_clean_load_phef
  declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from \$PHEF" >&2; exit 1; }
  mkdir -p "\$PGS_DIR"
  rm -f -- "\$PGS_DONE"

  local status_file="\$GWAS_POST_TMP/\$GWAS.pgs.status.tsv"
  local pgs_match_status score_chrs intermediate_status meta_tmp done_tmp done_status=complete
  pgs_plink_calc \
    --input "\$PGS_SCORE_FILE" \
    --pfile-dir "\$PGS_PFILE_DIR" \
    --output "\$PGS_OUTPUT" \
    --label "\$GWAS" \
    --work-dir "\$PGS_DIR" \
    --status-file "\$status_file" \
    --threads "\$PGS_THREADS"
  pgs_match_status=\$(awk -F '\t' '\$1=="matched_variants"{print \$2}' "\$status_file")
  score_chrs=\$(awk -F '\t' '\$1=="chromosomes"{print \$2}' "\$status_file")
  intermediate_status=\$(awk -F '\t' '\$1=="chromosome_intermediates"{print \$2}' "\$status_file")
  [[ "\$pgs_match_status" == scored || "\$pgs_match_status" == none ]] || { echo "ERROR: invalid PGS status: \$pgs_match_status" >&2; exit 1; }
  [[ -n "\$score_chrs" && -n "\$intermediate_status" ]] || { echo "ERROR: incomplete PGS status file: \$status_file" >&2; exit 1; }
  [[ "\$pgs_match_status" != none ]] || done_status=complete_no_matched_variants

  meta_tmp="\${PGS_META}.tmp.\$\$"
  {
    printf 'key\tvalue\n'
    printf 'gwas\t%s\nsource\t%s\neffect_allele\trefA\nscore_weight\tbJ\nscore_stat\tdosage_weighted_sum\nmissing_genotype\tno_mean_imputation\nmatched_variants\t%s\npfile_dir\t%s\nchromosomes\t%s\noutput\t%s\nchromosome_intermediates\t%s\n' \
      "\$GWAS" "\$PGS_SCORE_FILE" "\$pgs_match_status" "\$PGS_PFILE_DIR" "\$score_chrs" "\$PGS_OUTPUT" "\$intermediate_status"
  } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$PGS_META"
  done_tmp="\${PGS_DONE}.tmp.\$\$"
  printf 'GWAS\tSTATUS\tTIME\n%s\t%s\t%s\n' "\$GWAS" "\$done_status" "\$(date '+%F %T')" > "\$done_tmp"
  mv -f "\$done_tmp" "\$PGS_DONE"
  if [[ "\$pgs_match_status" == none ]]; then
    gwas_post_log "PGS done with no matched variants: \$PGS_OUTPUT (header only; PLINK2 logs retained)"
  else
    gwas_post_log "PGS done: \$PGS_OUTPUT (per-chromosome working files removed from \$PGS_DIR)"
  fi
}

# Load the index and performance overrides after all legacy definitions so
# the optimized implementations replace them deliberately.
gwas_post_magma_annotation_impl=\$(declare -f gwas_post_magma_annotation)
source "\$PERF_F" --worker
# The performance helper's MAGMA override derives annotation IDs from
# REFGEN_CLUMP.  Those IDs are not guaranteed to match MAGMA_REF (for example,
# GRCh38 pvars use CHR:POS:REF:ALT while g1000_eur.bim uses rsIDs).  Keep the
# optimized p-value preparation, but restore the compatible annotation builder.
eval "\$gwas_post_magma_annotation_impl"
unset gwas_post_magma_annotation_impl

if [[ "\$DO_STEP" == "all" ]]; then
  gwas_post_validate_format_columns "\$RAW"
  DO_STEP=small; gwas_clean_run_core
  gwas_post_fill_missing_fields "\$SMALL"
  gwas_post_ensure_index "\$SMALL"
  printf '%s\n' "\$GRCH" > "\${SMALL}.grch"
  if [[ "\$DO_LIFTOVER" == TRUE ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP=all
  gwas_post_ensure_index "\$FINAL"
  gwas_post_thin
  gwas_post_prepare_views
  gwas_post_mplot
  [[ -z "\$CIS_BED" ]] || gwas_clean_make_cis
else
  REQUESTED_STEP="\$DO_STEP"
  if [[ ",\$REQUESTED_STEP," == *,format,* ]]; then
    gwas_post_validate_format_columns "\$RAW"
    DO_STEP=small; gwas_clean_run_core
    gwas_post_fill_missing_fields "\$SMALL"
    gwas_post_ensure_index "\$SMALL"
  printf '%s\n' "\$GRCH" > "\${SMALL}.grch"
  fi
  if [[ ",\$REQUESTED_STEP," == *,liftover,* ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP="\$REQUESTED_STEP"
  gwas_post_ensure_index "\$FINAL"
  gwas_post_thin
  gwas_post_prepare_views
  gwas_post_magma
  gwas_post_mplot
  if [[ ",\$REQUESTED_STEP," == *,cis,* ]]; then gwas_clean_make_cis; fi
  if [[ ",\$REQUESTED_STEP," == *,lead,* ]]; then DO_STEP=lead; else DO_STEP=none; fi
fi

$(declare -f gwas_lead_marker_matches)
$(declare -f gwas_post_run_cojo)

gwas_post_clump_complete(){
  gwas_lead_marker_matches "\$CLUMP_DONE" "\$GWAS" clump "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" || return 1
  [[ -s "\$CLUMP_DONE" ]] && {
    [[ -s "\${MERGED}.clumps" ]] || gwas_post_lead_awk_has_no_rows ||
      awk -F '\t' 'NR==2&&\$2=="clump"&&\$6=="no_reference_matched_variants"{ok=1}END{exit !ok}' "\$CLUMP_DONE"
  }
}

gwas_post_cojo_complete(){
  gwas_lead_marker_matches "\$COJO_DONE" "\$GWAS" cojo "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" || return 1
  [[ -s "\$COJO_DONE" ]] && {
    [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]] || gwas_post_lead_awk_has_no_rows ||
      awk -F '\t' 'NR==2&&\$2=="cojo"&&(\$6=="no_reference_matched_variants"||\$6=="no_snps_selected"){ok=1}END{exit !ok}' "\$COJO_DONE"
  }
}

gwas_post_lead_awk_has_no_rows(){
  [[ -s "\$AWK_SNP" ]] || return 1
  awk 'NR>1{found=1} END{exit found ? 1 : 0}' "\$AWK_SNP"
}

gwas_post_lead_done(){
  [[ "\$REPLACE" != "TRUE" ]] || return 1
  [[ -s "\$AWK_SNP" ]] || return 1
  [[ ! -d "\$CLUMP" && ! -d "\$COJO" ]] || return 1
  gwas_post_clump_complete && gwas_post_cojo_complete
}

gwas_post_mark_phase_done(){
  local phase="\$1" status="\${2:-complete}" marker tmp
  case "\$phase" in
    clump) marker="\$CLUMP_DONE" ;;
    cojo) marker="\$COJO_DONE" ;;
    *) echo "ERROR: unknown lead phase: \$phase" >&2; return 2 ;;
  esac
  tmp="\${marker}.tmp.\$\$"
  mkdir -p "\$(dirname "\$marker")"
  {
    printf 'GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$GWAS" "\$phase" "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" "\$status"
  } > "\$tmp"
  mv -f "\$tmp" "\$marker"
}

# Exclude ambiguous PLINK IDs. Reuse the build/chromosome output until deleted.
gwas_post_filter_clump_multiallelic(){
  local ref="\$1" tag="\$2" assoc="\$3" pvar pvar_real bad_dir bad ref_meta tmp meta_tmp all_ids multi_ids dup_ids filtered n_before n_after n_removed
  if [[ ! -s "\${ref}.pvar" && -s "\${ref}.pvar.zst" ]]; then
    pvar="\${ref}.pvar.zst"
  else
    pvar="\${ref}.pvar"
  fi
  gwas_post_need_file "\$pvar"
  pvar_real=\$(readlink -f -- "\$pvar")
  [[ -n "\$pvar_real" ]] || { echo "ERROR: cannot resolve reference pvar: \$pvar" >&2; exit 1; }
  bad_dir="\$REF_BAD_DIR/GRCh\${GRCH}"
  bad="\$bad_dir/\${tag}.ambiguous.snp"
  ref_meta="\${bad}.reference.tsv"
  mkdir -p "\$bad_dir"

  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock is required for shared lead reference caches" >&2; exit 1; }
  lock_file="\${bad}.lock"
  exec {bad_lock_fd}> "\$lock_file"
  flock "\$bad_lock_fd"
  if [[ ! -e "\$bad" ]]; then
    tmp="\${bad}.tmp.\$\$"
    meta_tmp="\${ref_meta}.tmp.\$\$"
    all_ids="\${tmp}.all"
    multi_ids="\${tmp}.multi"
    dup_ids="\${tmp}.dup"
    : > "\$multi_ids"
    : > "\$dup_ids"
    gwas_clean_zcat "\$pvar" | awk -v FS='\t' -v multi="\$multi_ids" '
      \$1=="#CHROM" {for(i=1;i<=NF;i++){x=\$i; sub(/^#/,"",x); if(x=="ID")id=i; if(x=="ALT")alt=i} next}
      /^##/ {next}
      id>0 && \$id!="" && \$id!="." {print \$id; if(alt>0 && index(\$alt,",")>0) print \$id > multi}
    ' > "\$all_ids"
    sort "\$all_ids" | uniq -d > "\$dup_ids"
    cat "\$multi_ids" "\$dup_ids" 2>/dev/null | sort -u > "\$tmp"
    mv -f "\$tmp" "\$bad"
    printf 'GRCH\tPVAR\n%s\t%s\n' "\$GRCH" "\$pvar_real" > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$ref_meta"
    rm -f "\$all_ids" "\$multi_ids" "\$dup_ids"
  fi
  flock -u "\$bad_lock_fd"
  exec {bad_lock_fd}>&-

  n_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  filtered="\${assoc}.filtered.\$\$"
  awk -v FS='\t' -v OFS='\t' '
    FILENAME==ARGV[1] {bad[\$1]=1; next}
    FNR==1 {print; next}
    !(\$1 in bad) {print}
  ' "\$bad" "\$assoc" > "\$filtered"
  mv -f "\$filtered" "\$assoc"
  n_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  n_removed=\$((n_before-n_after))
  printf 'GWAS\tCHR\tN_BEFORE\tN_REMOVED_AMBIGUOUS_REF_ID\tN_AFTER\n%s\t%s\t%s\t%s\t%s\n' \
    "\$GWAS" "\$tag" "\$n_before" "\$n_removed" "\$n_after" > "\${QC_PREFIX}.\${tag}.clump_multiallelic.tsv"
  gwas_post_log "clump ambiguous-reference-ID filter \$tag: before=\$n_before removed=\$n_removed after=\$n_after cache=GRCh\${GRCH}/\${tag}.ambiguous.snp"
}

gwas_post_prep_lead_inputs(){
  local refs="\$1" suffix=".tmp.\$\$" tag assoc ma
  rm -f -- "\${QC_PREFIX}".*.cojo_skip.log
  while IFS=\$'\t' read -r _ tag _ assoc ma _; do
    rm -f "\${assoc}\${suffix}" "\${ma}\${suffix}"
  done < "\$refs"
  gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v refs="\$refs" \
    -v suffix="\$suffix" -v default_n="\$GWAS_N" -v p_lead="\$P_LEAD" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function valid(x){return x!="" && x!="NA" && x!="."}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;return x}
    function get(name){return name in c ? \$(c[name]) : "NA"}
    BEGIN{
      while((getline line < refs)>0){split(line,a,"\t");chr=a[3];if(chr=="")continue;
        assoc[chr]=a[4] suffix;ma[chr]=a[5] suffix
        print "SNP","CHR","POS","EA","NEA","P" > assoc[chr]
        print "SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","N" > ma[chr]}
      close(refs)
    }
    NR==1{
      for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}
      required[1]="SNP";required[2]="CHR";required[3]="POS";required[4]="EA"
      required[5]="NEA";required[6]="BETA";required[7]="SE";required[8]="P"
      for(i=1;i<=8;i++)if(!(required[i] in c)){print "ERROR: required GWAS column is missing: " required[i] > "/dev/stderr";fatal=1}
      if(fatal)exit 2
      next
    }
    {
      snp=get("SNP");chr=normchr(get("CHR"));pos=get("POS");ea=get("EA");nea=get("NEA")
      beta=get("BETA");se=get("SE");p=get("P");eaf=get("EAF");n=get("N")
      # Filter significant variants before the in-memory ID match.
      if(!(chr in assoc)||pos!~/^[0-9]+\$/||pos+0<=0||!valid(snp)||!valid(ea)||!valid(nea)||!isnum(p)||p+0>p_lead)next
      print snp,chr,pos,ea,nea,p > assoc[chr]
      if(!isnum(beta)||!isnum(se))next
      if(!isnum(n))n=default_n
      print snp,chr,pos,ea,nea,eaf,beta,se,p,n > ma[chr]
    }
  '
  while IFS=\$'\t' read -r _ _ _ assoc ma _; do
    mv -f "\${assoc}\${suffix}" "\$assoc";mv -f "\${ma}\${suffix}" "\$ma"
  done < "\$refs"
}

gwas_post_subset_match_reference(){
  local ref="\$1" query="\$2" out="\$3" kind="\$4" ref_file="" ext
  case "\$kind" in
    pvar)
      case "\$ref" in
        *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst) ref_file="\$ref" ;;
        *)
          for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    bim)
      case "\$ref" in
        *.bim|*.bim.gz|*.bim.bgz) ref_file="\$ref" ;;
        *)
          for ext in .bim .bim.gz .bim.bgz; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    *) echo "ERROR: invalid reference subset kind: \$kind" >&2; return 2 ;;
  esac
  gwas_post_need_file "\$ref_file"

  # Restrict the reference to query coordinates before match_SNP.
  awk -v FS='[ \t]+' -v kind="\$kind" '
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x=="" ? "0" : x}
    NR==FNR{
      if(FNR==1){for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);if(h=="CHR")qchr=i;else if(h=="POS")qpos=i}next}
      if(!qchr||!qpos){fatal=1;next}
      if(\$(qpos)~/^[0-9]+\$/)wanted[normchr(\$(qchr)) SUBSEP \$(qpos)]=1
      next
    }
    kind=="pvar" && /^##/{print;next}
    kind=="pvar" && /^#CHROM/{print;header=1;next}
    kind=="pvar"{
      if((normchr(\$1) SUBSEP \$2) in wanted){print;kept++}
      next
    }
    kind=="bim"{
      if(NF>=4 && (normchr(\$1) SUBSEP \$4) in wanted){print;kept++}
      next
    }
    END{
      if(fatal){print "ERROR: lead query requires CHR and POS columns" > "/dev/stderr";exit 2}
      if(kind=="pvar"&&!header){print "ERROR: invalid .pvar header" > "/dev/stderr";exit 2}
      # Keep an empty BIM subset non-empty so match_SNP can report all query
      # rows as unmatched instead of rejecting the reference argument.
      if(kind=="bim"&&kept==0)print "0\t.\t0\t0\tN\tN"
    }
  ' "\$query" <(gwas_clean_zcat "\$ref_file") > "\$out"
}

# GCTA treats only numeric chromosomes up to --autosome-num as usable.  PLINK
# references commonly encode chrX/chrY as X/Y, which makes GCTA silently load
# just one unusable record from an otherwise complete sex-chromosome BIM.  For
# sex chromosomes, make a small per-trait BED containing only the already
# matched COJO candidates and emit numeric chromosome codes (X=23, Y=24).
gwas_post_prepare_gcta_bfile(){
  local source="\$1" chr="\$2" ma="\$3" out="\$4"
  GCTA_BFILE="\$source"
  GCTA_CHR_ARGS=()
  (( chr > 22 )) || return 0
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 is required to prepare chr\$chr for GCTA" >&2; return 1; }
  gwas_post_need_file "\${source}.bed"
  gwas_post_need_file "\${source}.bim"
  gwas_post_need_file "\${source}.fam"
  gwas_post_need_file "\$ma"
  rm -f -- "\${out}.bed" "\${out}.bim" "\${out}.fam" "\${out}.log" "\${out}.err"
  run_tool "\$out" plink2 --bfile "\$source" --extract "\$ma" --make-bed --output-chr 26 --out "\$out"
  gwas_post_need_file "\${out}.bed"
  gwas_post_need_file "\${out}.bim"
  gwas_post_need_file "\${out}.fam"
  GCTA_BFILE="\$out"
  GCTA_CHR_ARGS=(--autosome-num "\$chr")
}

gwas_post_match_lead_inputs(){
  local ref="\$1" cref="\$2" tag="\$3" assoc="\$4" ma="\$5"
  local matched ref_subset n_assoc_before n_assoc_after n_ma_before n_ma_after
  declare -F match_SNP >/dev/null 2>&1 || gwas_clean_load_phef

  n_assoc_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  if (( n_assoc_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.clump.reference.pvar"
    gwas_post_subset_match_reference "\$ref" "\$assoc" "\$ref_subset" pvar
    matched="\${assoc}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.clump_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv" "\$assoc"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$6+0)<best[\$1]{best[\$1]=\$6+0;line[\$1]=\$1 OFS \$6}
      END{print "SNP","P";for(snp in line)print line[snp]}' "\$matched" > "\${assoc}.tmp.\$\$"
    mv -f "\${assoc}.tmp.\$\$" "\$assoc"; rm -f "\$matched"
    n_assoc_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  else
    rm -f "\${QC_PREFIX}.\${tag}.clump_match.tsv" "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv"
    n_assoc_after=0
  fi

  n_ma_before=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  if (( n_ma_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.cojo.reference.bim"
    gwas_post_subset_match_reference "\$cref" "\$ma" "\$ref_subset" bim
    matched="\${ma}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.cojo_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv" "\$ma"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$9+0)<best[\$1]{best[\$1]=\$9+0;line[\$1]=\$1 OFS \$4 OFS \$5 OFS \$6 OFS \$7 OFS \$8 OFS \$9 OFS \$10}
      END{print "SNP","A1","A2","freq","b","se","p","N";for(snp in line)print line[snp]}' "\$matched" > "\${ma}.tmp.\$\$"
    mv -f "\${ma}.tmp.\$\$" "\$ma"; rm -f "\$matched"
    n_ma_after=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  else
    rm -f "\${QC_PREFIX}.\${tag}.cojo_match.tsv" "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv"
    n_ma_after=0
  fi

  gwas_post_log "reference-ID match \$tag: clump=\$n_assoc_after/\$n_assoc_before cojo=\$n_ma_after/\$n_ma_before"
}

if run_lead; then
  if gwas_post_lead_done; then
    gwas_post_log "lead SNP discovery exists: \$AWK_SNP"
  else
    lead_t0=\$(date +%s)
    clump_was_done=FALSE
    cojo_was_done=FALSE
    cojo_skipped=FALSE
    clump_has_usable_input=FALSE
    cojo_has_matched_input=FALSE
    cojo_has_empty_selection=FALSE
    if [[ "\$REPLACE" != "TRUE" ]]; then
      gwas_post_clump_complete && clump_was_done=TRUE
      gwas_post_cojo_complete && cojo_was_done=TRUE
    fi
    gwas_post_need_file "\$FINAL"
    labels="\${GWAS_POST_TMP}/labels.\$\$.txt"
    refs="\${GWAS_POST_TMP}/refs.\$\$.tsv"
    : > "\$labels"
    : > "\$refs"

  # 1) Resolve selected reference chromosomes. GWAS coordinates, alleles, and
  # frequencies are never filled from the reference panel.
  while read -r ref lab chr unused; do
    [[ -n "\$ref" ]] || continue
    want_chr "\$chr" "\$CHRS" || continue
    tag=\$(chr_label "\$chr")
    echo "\$tag" >> "\$labels"
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    assoc="\${cp}.assoc"
    ma="\${jp}.ma"
    cref=\$(cojo_bfile "\$REFGEN_COJO" "\$chr")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$ref" "\$tag" "\$chr" "\$assoc" "\$ma" "\$cref" >> "\$refs"
  done < <(ref_clump_pfiles "\$REFGEN_CLUMP")
  gwas_post_check_chromosome_coverage lead "\$refs"
  rm -f -- "\${MERGED}.lead.done"
  if [[ "\$clump_was_done" != "TRUE" ]]; then
    rm -f -- "\$CLUMP_DONE" "\${MERGED}.clumps"
  fi
  if [[ "\$cojo_was_done" != "TRUE" ]]; then
    rm -f -- "\$COJO_DONE" "\${MERGED}.jma.cojo" "\${MERGED}.cma.cojo" "\${MERGED}.ldr.cojo"
  fi
  rm -rf -- "\$CLUMP" "\$COJO"
  mkdir -p "\$CLUMP" "\$COJO"
  gwas_post_log "prepare P<=\$P_LEAD lead inputs from required SNP/CHR/POS/EA/NEA/BETA/SE/P fields: \$FINAL"
  gwas_post_prep_lead_inputs "\$refs"

  # 2) No-LD top SNP per chromosome/window and per-chromosome tool inputs.
  gwas_post_log "awk distance lead SNPs: \${GWAS_POST_LEAD_VIEW:-\$FINAL} -> \$AWK_SNP"
  gwas_post_zcat "\${GWAS_POST_LEAD_VIEW:-\$FINAL}" | awk -v FS='\t' -v OFS='\t' -v pthr="\$P_LEAD" -v win="\$LEAD_WINDOW" -f <(awk_lead) > "\${AWK_SNP}.tmp.\$\$"
  mv -f "\${AWK_SNP}.tmp.\$\$" "\$AWK_SNP"
  awk -v g="\$GWAS" 'NR==1{next} END{print g"\\t"NR-1}' "\$AWK_SNP" > "\${QC_PREFIX}.awk.nrow.tsv"

  # 3) plink2 --clump and gcta --cojo per chromosome.
  while IFS=\$'\t' read -r ref tag chr assoc ma cref; do
    [[ -n "\$ref" ]] || continue
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    gwas_post_match_lead_inputs "\$ref" "\$cref" "\$tag" "\$assoc" "\$ma"
    # Inputs were already restricted to P<=P_LEAD before reference matching,
    # which bounds memory independently of the full GWAS size.
    cojo_skip_log="\${QC_PREFIX}.\${tag}.cojo_skip.log"
    awk -v FS='\t' -v OFS='\t' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{print "STATUS","REASON","SNP","A1","A2","EAF","P"}
      NR>1 && (!isnum(\$4)||\$4+0<0||\$4+0>1){print "SKIP_GCTA","missing_or_invalid_EAF_after_reference_match",\$1,\$2,\$3,\$4,\$7}
    ' "\$ma" > "\${cojo_skip_log}.tmp.\$\$"
    if gwas_post_has_data_rows "\${cojo_skip_log}.tmp.\$\$"; then
      mv -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    else
      rm -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    fi
    if [[ "\$clump_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP plink2 chr\$chr: clump phase already complete"
    elif gwas_post_has_data_rows "\$assoc"; then
      gwas_post_filter_clump_multiallelic "\$ref" "\$tag" "\$assoc"
      if gwas_post_has_data_rows "\$assoc"; then
        clump_has_usable_input=TRUE
        pfile_modifier=()
        [[ -s "\${ref}.pvar" || ! -s "\${ref}.pvar.zst" ]] || pfile_modifier=(vzs)
        pfile_keep=()
        [[ -z "\$REFGEN_KEEP" ]] || pfile_keep=(--keep "\$REFGEN_KEEP")
        run_tool "\$cp" plink2 --pfile "\$ref" "\${pfile_modifier[@]}" "\${pfile_keep[@]}" --clump "\$assoc" --clump-snp-field SNP --clump-p-field P --clump-p1 "\$P_LEAD" --clump-kb "\$CLUMP_KB" --out "\$cp"
      else
        gwas_post_log "skip plink2 chr\$chr: all significant IDs are ambiguous in the reference"
      fi
    else
      gwas_post_log "skip plink2 chr\$chr: no SNPs in \$assoc"
    fi

    # 3) gcta --cojo: COJO selection with the chromosome-specific bed/bim/fam reference
    if [[ "\$cojo_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP gcta chr\$chr: COJO phase already complete"
    elif [[ -s "\${QC_PREFIX}.\${tag}.cojo_skip.log" ]]; then
      cojo_skipped=TRUE
      gwas_post_log "SKIP gcta chr\$chr: significant variant has missing/invalid EAF; see \${QC_PREFIX}.\${tag}.cojo_skip.log"
    elif gwas_post_has_data_rows "\$ma"; then
      gwas_post_prepare_gcta_bfile "\$cref" "\$chr" "\$ma" "\${jp}.gcta_ref"
      # A matched SNP can be monomorphic in the reference population. GCTA
      # 1.94.1 can segfault on its zero genotype variance during COJO selection.
      # This floor excludes fixed alleles while retaining rare 1000G variants.
      gwas_post_run_cojo "\$jp" "\${QC_PREFIX}.\${tag}.cojo_empty" gcta --bfile "\$GCTA_BFILE" "\${GCTA_CHR_ARGS[@]}" \
        --maf 1e-6 --cojo-file "\$ma" --cojo-slct --cojo-p "\$P_LEAD" --out "\$jp"
      if [[ "\$GCTA_COJO_HAS_MATCH" == TRUE ]]; then cojo_has_matched_input=TRUE; fi
      if [[ "\$GCTA_COJO_NO_SNPS_SELECTED" == TRUE ]]; then cojo_has_empty_selection=TRUE; fi
    else
      gwas_post_log "skip gcta chr\$chr: no SNPs in \$ma"
    fi
  done < "\$refs"

    if gwas_post_lead_awk_has_no_rows; then
      rm -f "\$labels" "\$refs"
      rm -rf -- "\$CLUMP" "\$COJO"
      gwas_post_mark_phase_done clump no_significant_variants
      gwas_post_mark_phase_done cojo no_significant_variants
    else
      if [[ "\$clump_was_done" != "TRUE" ]]; then
        concat_chr_outputs "\$CLUMP" "\$MERGED" "\$labels" clumps
        if [[ -s "\${MERGED}.clumps" ]]; then
          gwas_post_mark_phase_done clump
        elif [[ "\$clump_has_usable_input" != "TRUE" ]]; then
          gwas_post_log "clump complete with no reference-matched significant variants"
          gwas_post_mark_phase_done clump no_reference_matched_variants
        fi
      fi
      if [[ "\$cojo_was_done" != "TRUE" && "\$cojo_skipped" != "TRUE" ]]; then
        concat_chr_outputs "\$COJO" "\$MERGED" "\$labels" jma.cojo cma.cojo ldr.cojo
        if [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]]; then
          gwas_post_mark_phase_done cojo
        elif [[ "\$cojo_has_empty_selection" == TRUE ]]; then
          gwas_post_log "COJO complete with no SNPs selected"
          gwas_post_mark_phase_done cojo no_snps_selected
        elif [[ "\$cojo_has_matched_input" != "TRUE" ]]; then
          gwas_post_log "COJO complete with no reference-matched significant variants"
          gwas_post_mark_phase_done cojo no_reference_matched_variants
        fi
      fi
      rm -f "\$labels" "\$refs"
      if ! gwas_post_clump_complete; then
        echo "ERROR: clump phase did not create \${MERGED}.clumps; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      if ! gwas_post_cojo_complete; then
        echo "ERROR: COJO phase is incomplete for \$GWAS; GCTA skip details are under \${QC_PREFIX}.*.cojo_skip.log; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      gwas_post_log "delete lead intermediates: \$CLUMP \$COJO \${MERGED}.cma.cojo"
      rm -rf -- "\$CLUMP" "\$COJO"
      rm -f -- "\${MERGED}.cma.cojo"
    fi
    gwas_post_log "lead SNP discovery done in \$((\$(date +%s)-lead_t0)) sec"
  fi
fi

# PGS depends on the genome-wide .jma.cojo created by the lead/COJO phase, so
# it is deliberately invoked after lead even though it is declared after mplot.
gwas_post_pgs

if [[ ",\$H2_REQUESTED," == *,h2,* || "\$H2_REQUESTED" == all ]]; then
  h2_args=()
  [[ -z "\$H2_PYTHON" ]] || h2_args+=(--python "\$H2_PYTHON")
  python3 "\$H2_HELPER" h2 --gwas-file "\$FINAL" --trait "\$GWAS" \
    --output-dir "\$(dirname "\$GWAS_DIR")/h2" --sex "\$H2_SEX" \
    --ref-ld-chr "\$H2_REF_LD" --w-ld-chr "\$H2_W_LD" --merge-alleles "\$H2_MERGE_ALLELES" \
    --conda-env "\$H2_CONDA_ENV" --replace "\$REPLACE" "\${h2_args[@]}"
fi

# Delete raw only after every requested module has succeeded.
if [[ "\$DELETE_RAW_AFTER_SUCCESS" == TRUE && -n "\$RAW" && -f "\$RAW" ]]; then
  if [[ "\$RAW" -ef "\$FINAL" ]]; then
    echo "ERROR: raw and final GWAS are the same file; refusing to delete \$RAW" >&2
    exit 1
  fi
  gzip -t -- "\$FINAL"
  rm -- "\$RAW"
  gwas_post_log "All requested modules succeeded; deleted raw: \$RAW"
fi

CMD_TOP

		# Workers run via bash; do not require chmod on mounted project drives.
		echo "$cmd"
	}

fi

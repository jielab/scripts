#!/usr/bin/env bash
# Variables are provided by the format.sh coordinator.
# shellcheck disable=SC2154

# Configure paths cheaply. Copy/cache work starts only for the first GWAS that
# failed the existing completion check, so a finished project remains fast.
gwas_storage_configure() {
	local key
	local_io=${local_io^^}
	case "$local_io" in
		AUTO)
			local_io=FALSE
			if [[ "$is_bsub" != TRUE ]] && grep -qi microsoft /proc/sys/kernel/osrelease 2>/dev/null; then
				local_io=TRUE
			fi ;;
		TRUE|FALSE) ;;
		*) echo 'ERROR: --local-io must be AUTO, TRUE, or FALSE' >&2; return 2 ;;
	esac
	gwas_storage_ready=FALSE
	[[ "$local_io" == TRUE ]] || return 0
	local_cache_root=$(realpath -m "${local_cache_arg:-${XDG_CACHE_HOME:-$HOME/.cache}/gwas/format}")
	key=$(printf '%s\n%s\n' "$(realpath -m "$dir_out")" "$category" | sha256sum | cut -c1-20)
	local_project_cache="$local_cache_root/projects/$key"
	: "${GWAS_POST_TMP_BASE:=/tmp/gwas-post-$UID/$key/$step_key}"
	GWAS_POST_TMP_BASE=$(realpath -m "$GWAS_POST_TMP_BASE")
	legacy_rsid_cache=""
	if [[ "$rsid_cache_explicit" == FALSE ]]; then
		legacy_rsid_cache="$rsid_cache"
		rsid_cache="$local_project_cache/rsid"
	fi
	log "Local I/O: scratch=$GWAS_POST_TMP_BASE rsid-cache=$rsid_cache references=$local_cache_root/references"
}

gwas_storage_prepare() {
	[[ "$local_io" == TRUE && "$gwas_storage_ready" != TRUE ]] || return 0
	local helper="${SCRIPT_PATH%/*}/f/format.storage.py" result
	local -a staged inputs
	mkdir -p "$GWAS_POST_TMP_BASE" "$local_cache_root"
	if has_step format && [[ "$rsid_mode" == TRUE && -n "$legacy_rsid_cache" ]]; then
		python3 "$helper" seed-rsid --source "$legacy_rsid_cache" --target "$rsid_cache" || return
	fi
	if wants_magma; then
		result=$(python3 "$helper" stage --root "$local_cache_root/references" \
			"${magma_ref}.bed" "${magma_ref}.bim" "${magma_ref}.fam" "$synonyms") || return
		mapfile -t staged <<<"$result"
		magma_ref="${staged[0]%.bed}"
		synonyms="${staged[3]}"
	fi
	if { has_step format && [[ "$hm3_mode" == TRUE || "$thin" == TRUE ]]; } || has_step thin || has_step mplot; then
		inputs=("$hm3_file")
		if [[ -n "$hm3_pos" && "$grch" == auto ]]; then
			inputs+=("${hm3_pos//\{grch\}/37}" "${hm3_pos//\{grch\}/38}")
		elif [[ -n "$hm3_pos" ]]; then
			inputs+=("${hm3_pos//\{grch\}/$grch}")
		fi
		result=$(python3 "$helper" stage --root "$local_cache_root/references" "${inputs[@]}") || return
		mapfile -t staged <<<"$result"
		hm3_file="${staged[0]}"
		[[ -z "$hm3_pos" ]] || hm3_pos="${staged[1]%/*}/${hm3_pos##*/}"
	fi
	gwas_storage_ready=TRUE
}

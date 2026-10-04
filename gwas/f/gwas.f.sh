#!/usr/bin/env bash
set -euo pipefail


# 🚩 gwas_extract
gwas_extract() (
	set -euo pipefail
	here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
	if [[ ${1:-} == --intersect ]]; then
		requested=$2 shared=$3 output=$4
	else
		counts=$1 output=$2 maf=$3
	fi
	tmp=$(mktemp "${output}.tmp.XXXXXX")
	trap 'rm -f -- "$tmp" "$tmp.counts.tsv"' EXIT
	if [[ ${1:-} == --intersect ]]; then
		awk 'FILENAME == ARGV[1] { if (NF) wanted[$1]=1; next } $1 in wanted { print $1 }' "$requested" "$shared" >"$tmp"
	else
		LC_ALL=C awk -v maf="$maf" -v summary="$tmp.counts.tsv" -f "$here/gwas.extract.awk" "$counts" >"$tmp"
	fi
	if [[ ! -s "$tmp" ]]; then
		echo "ERROR: No variants pass for $output; stopping before association" >&2
		exit 2
	fi
	[[ ! -f "$tmp.counts.tsv" ]] || mv -f -- "$tmp.counts.tsv" "$output.counts.tsv"
	mv -f -- "$tmp" "$output"
)


# 🚩 gwas_saige
gwas_saige() (
	# Keep the conda R ABI separate from the user's global R library/profile.
	set -euo pipefail
	saige_env=${SAIGE_ENV_DIR:-$HOME/anaconda3/envs/saige}
	if [[ ! -x "$saige_env/bin/Rscript" ]]; then
		saige_env=$HOME/miniforge3/envs/saige
	fi
	[[ -x "$saige_env/bin/Rscript" ]] || {
		echo 'SAIGE R environment not found' >&2
		exit 1
	}
	export R_LIBS="$saige_env/lib/R/library"
	export R_LIBS_USER="$saige_env/lib/R/library"
	export R_LIBS_SITE="$saige_env/lib/R/library"
	exec "$saige_env/bin/Rscript" --vanilla "$@"
)


# 🚩 Dispatch
case "${1:-help}" in
	extract)
		shift
		gwas_extract "$@"
		;;
	saige | saige-rscript)
		shift
		gwas_saige "$@"
		;;
	*.R) gwas_saige "$@" ;;
	-h | --help | help) echo 'Usage: gwas.f.sh extract|saige-rscript [arguments]' ;;
	*)
		echo "Unknown GWAS helper: $1" >&2
		exit 2
		;;
esac

#!/usr/bin/env bash


# 🚩 c2.cause
# C2 utilities: native MR-link-2 runner and existing cis-file indexing.
set -euo pipefail

index_cis_main() (
	# One-time migration of existing gwas_post cis outputs to coordinate-sorted
	# BGZF plus tabix indexes. New cis outputs are indexed directly by the GWAS formatting module.

	set -euo pipefail
	export LC_ALL=C

	DATA_ROOT="${CIS_INDEX_DATA_ROOT:-/mnt/d/data/gwas}"
	PROJECTS="${CIS_INDEX_PROJECTS:-main,prot,met}"
	JOBS="${CIS_INDEX_JOBS:-6}"
	THREADS="${CIS_INDEX_THREADS:-2}"
	SORT_MEMORY="${CIS_INDEX_SORT_MEMORY:-512M}"
	TMP_DIR="${CIS_INDEX_TMP_DIR:-/tmp}"
	INDEX_F="${CIS_INDEX_FUNCTIONS:-/mnt/d/scripts/gwas/f/format.f.sh}"
	REPLACE=FALSE
	DRY_RUN=FALSE

	usage() {
		cat <<'EOF'
Usage: ./c2.cause.sh index-cis [options]

Convert existing canonical cis files
  <data-root>/<project>/common/<trait>/gwas/<trait>.cis.gz
to coordinate-sorted BGZF and create <trait>.cis.gz.tbi (or .csi when needed).

Options:
  --data-root PATH       GWAS parent directory [/mnt/d/data/gwas]
  --projects CSV         Projects to scan [main,prot,met]
  --jobs N               Files processed concurrently [6]
  --threads N            bgzip threads per file [2]
  --sort-memory SIZE     GNU sort memory per file if sorting is needed [512M]
  --tmp-dir PATH         GNU sort temporary directory [/tmp]
  --index-f FILE         GWAS formatting module (loads --index) [/mnt/d/scripts/gwas/f/format.f.sh]
  --replace TRUE|FALSE   rebuild even a valid existing index [FALSE]
  --dry-run              list exact targets without changing them
  -h, --help             show this help

The original .cis.gz is replaced only after a temporary BGZF file and its
index pass validation. Existing valid, fresh .tbi/.csi files are skipped.
EOF
	}

	bool_word() {
		case "${1,,}" in
			true | 1 | yes | y) printf 'TRUE\n' ;;
			false | 0 | no | n) printf 'FALSE\n' ;;
			*) return 2 ;;
		esac
	}

	need_value() {
		[[ -n "${2-}" && "${2-}" != --* ]] || {
			echo "ERROR: $1 requires a value" >&2
			exit 2
		}
	}

	while [[ $# -gt 0 ]]; do
		case "$1" in
			--data-root)
				need_value "$1" "${2-}"
				DATA_ROOT="$2"
				shift 2
				;;
			--projects)
				need_value "$1" "${2-}"
				PROJECTS="$2"
				shift 2
				;;
			--jobs)
				need_value "$1" "${2-}"
				JOBS="$2"
				shift 2
				;;
			--threads)
				need_value "$1" "${2-}"
				THREADS="$2"
				shift 2
				;;
			--sort-memory)
				need_value "$1" "${2-}"
				SORT_MEMORY="$2"
				shift 2
				;;
			--tmp-dir)
				need_value "$1" "${2-}"
				TMP_DIR="$2"
				shift 2
				;;
			--index-f)
				need_value "$1" "${2-}"
				INDEX_F="$2"
				shift 2
				;;
			--replace)
				need_value "$1" "${2-}"
				REPLACE=$(bool_word "$2") || {
					echo "ERROR: --replace expects TRUE or FALSE" >&2
					exit 2
				}
				shift 2
				;;
			--dry-run)
				DRY_RUN=TRUE
				shift
				;;
			-h | --help)
				usage
				exit 0
				;;
			*)
				echo "ERROR: unknown option: $1" >&2
				usage >&2
				exit 2
				;;
		esac
	done

	[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || {
		echo "ERROR: --jobs must be a positive integer" >&2
		exit 2
	}
	[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
		echo "ERROR: --threads must be a positive integer" >&2
		exit 2
	}
	[[ -d "$DATA_ROOT" ]] || {
		echo "ERROR: missing data root: $DATA_ROOT" >&2
		exit 1
	}
	[[ -d "$TMP_DIR" ]] || {
		echo "ERROR: missing sort temp directory: $TMP_DIR" >&2
		exit 1
	}
	[[ -s "$INDEX_F" ]] || {
		echo "ERROR: missing GWAS index module: $INDEX_F" >&2
		exit 1
	}
	command -v gzip >/dev/null || {
		echo "ERROR: gzip not found" >&2
		exit 1
	}
	command -v bgzip >/dev/null || {
		echo "ERROR: bgzip not found (install htslib)" >&2
		exit 1
	}
	command -v tabix >/dev/null || {
		echo "ERROR: tabix not found (install htslib)" >&2
		exit 1
	}
	command -v sort >/dev/null || {
		echo "ERROR: GNU sort not found" >&2
		exit 1
	}

	# shellcheck source=/mnt/d/scripts/gwas/f/format.f.sh
	source "$INDEX_F" --index

	log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

	index_cis_one() {
		local file="$1"
		[[ -s "$file" ]] || {
			log "ERROR missing/empty: $file"
			return 1
		}
		if [[ "$REPLACE" != TRUE ]] && gwas_index_valid "$file"; then
			log "SKIP valid index: $file"
			return 0
		fi
		[[ "$REPLACE" != TRUE ]] || rm -f -- "${file}.tbi" "${file}.csi"
		log "INDEX cis: $file"
		gwas_index_file "$file"
		bgzip -t "$file" >/dev/null 2>&1 || {
			log "ERROR BGZF validation failed: $file"
			return 1
		}
		gwas_index_valid "$file" || {
			log "ERROR tabix validation failed: $file"
			return 1
		}
		log "DONE BGZF+index: $file"
	}

	target_list=$(mktemp "${TMPDIR:-/tmp}/index_cis.targets.XXXXXX")
	trap 'rm -f -- "$target_list"' EXIT

	IFS=',' read -r -a project_array <<<"$PROJECTS"
	for project in "${project_array[@]}"; do
		project=${project//[[:space:]]/}
		[[ -n "$project" && "$project" != "." && "$project" != ".." && "$project" != */* ]] || {
			echo "ERROR: invalid project name: $project" >&2
			exit 2
		}
		root="$DATA_ROOT/$project/common"
		[[ -d "$root" ]] || {
			log "WARNING missing project root, skip: $root"
			continue
		}
		while IFS= read -r file; do
			trait=$(basename "$(dirname "$(dirname "$file")")")
			[[ "$(basename "$file")" == "${trait}.cis.gz" ]] && printf '%s\n' "$file"
		done < <(find "$root" -mindepth 3 -maxdepth 3 -type f -path '*/gwas/*.cis.gz' -size +0c | sort -V)
	done | sort -u -V >"$target_list"

	target_count=$(wc -l <"$target_list" | tr -d ' ')
	log "Targets=$target_count projects=$PROJECTS jobs=$JOBS bgzip_threads=$THREADS"
	if ((target_count == 0)); then
		log "No existing canonical cis files found; nothing to do."
		exit 0
	fi

	if [[ "$DRY_RUN" == TRUE ]]; then
		cat "$target_list"
		exit 0
	fi

	export -f log gwas_index_log gwas_index_header gwas_index_col gwas_index_valid
	export -f gwas_index_make_sidecar gwas_index_file index_cis_one
	export REPLACE GWAS_CLEAN_COMP_THREADS="$THREADS"
	export GWAS_INDEX_SORT_MEMORY="$SORT_MEMORY" GWAS_INDEX_TMP_DIR="$TMP_DIR"
	export SHELL=/bin/bash
	if command -v parallel >/dev/null 2>&1; then
		parallel --line-buffer --halt soon,fail=1 -j "$JOBS" index_cis_one :::: "$target_list"
	else
		while IFS= read -r file; do index_cis_one "$file"; done <"$target_list"
	fi

	log "All existing canonical cis files are BGZF-indexed."
)

mr_link2_main() (
	set -euo pipefail

	# MR-link-2 runner for the 5C CAD framework.
	# The job TSV is generated by c2.cause.R.

	if [[ -n "${CONDA_SH:-}" ]]; then
		conda_sh="$CONDA_SH"
	elif command -v conda >/dev/null 2>&1; then
		conda_sh="$(conda info --base)/etc/profile.d/conda.sh"
	else
		conda_sh="/home/huangj/anaconda3/etc/profile.d/conda.sh"
	fi
	if [[ ! -r "$conda_sh" ]]; then
		echo "Conda initialization script not found: $conda_sh" >&2
		exit 1
	fi
	# shellcheck disable=SC1090
	set +u
	source "$conda_sh"
	conda activate le8 || {
		echo "Failed to activate conda environment: le8" >&2
		exit 1
	}
	set -u
	cleanup_conda() {
		set +u
		conda deactivate >/dev/null 2>&1 || true
	}
	trap cleanup_conda EXIT

	if [[ -d /mnt/d/software/bin ]]; then
		export PATH="/mnt/d/software/bin:${PATH}"
	fi

	script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
	prep_py="${script_dir}/c2.mr_link2.py"

	jobs=""
	cad_gwas=""
	ref_bed="${MRLINK2_REF_BED:-}"
	ref_pfile_dir="${MRLINK2_REF_PFILE_DIR:-/mnt/f/gen/1kg/${LE8_GRCH:-37}/pfile}"
	ref_pfile_pop="${MRLINK2_REF_POP:-EUR}"
	ref_samples="${MRLINK2_REF_SAMPLES:-}"
	ref_id_dir="${MRLINK2_REF_ID_DIR:-}"
	mrlink2="${MRLINK2_SCRIPT:-${script_dir}/c2.mr_link2.py}"
	outdir="c2_cause/link2"
	python_bin="${PYTHON_BIN:-python3}"
	verbose="${MRLINK2_VERBOSE:-0}"
	p_threshold="${MRLINK2_P_THRESHOLD:-5e-8}"
	region_padding="${MRLINK2_REGION_PADDING:-500000}"
	maf_threshold="${MRLINK2_MAF_THRESHOLD:-0.01}"
	max_correlation="${MRLINK2_MAX_CORRELATION:-0.99}"
	continue_analysis="${MRLINK2_CONTINUE:-true}"

	usage() {
		cat <<'EOF'
Usage:
  ./f/c2.cause.sh mr-link2 \
    --jobs c2_cause/link2/c2.link2.jobs.tsv \
    --cad-gwas /mnt/d/data/gwas/main/common/cvd_cad/gwas/cvd_cad.gz \
    --reference-pfile-dir /mnt/f/gen/1kg/38/pfile \
    --reference-pop EUR \
    --reference-id-dir /mnt/f/gen/1kg/38/id \
    --reference-samples /mnt/f/gen/1kg/38/samples.txt \
    --mrlink2 f/c2.mr_link2.py \
    --outdir c2_cause/link2

The --jobs file is generated by c2.cause.R with columns:
  omics  trait  exposure  outcome  region  [audit columns ...]

When a cleaned *.cis.gz exists, exposure and region are that dense marginal
cis file and its observed bounds. --region-padding is only a fallback for jobs
without a prespecified region.
EOF
	}

	while [[ $# -gt 0 ]]; do
		case "$1" in
			--jobs)
				jobs="$2"
				shift 2
				;;
			--cad-gwas)
				cad_gwas="$2"
				shift 2
				;;
			--reference-bed)
				ref_bed="$2"
				shift 2
				;;
			--reference-pfile-dir)
				ref_pfile_dir="$2"
				shift 2
				;;
			--reference-pop)
				ref_pfile_pop="$2"
				shift 2
				;;
			--reference-id-dir)
				ref_id_dir="$2"
				shift 2
				;;
			--reference-samples)
				ref_samples="$2"
				shift 2
				;;
			--mrlink2)
				mrlink2="$2"
				shift 2
				;;
			--outdir)
				outdir="$2"
				shift 2
				;;
			--python)
				python_bin="$2"
				shift 2
				;;
			--p-threshold)
				p_threshold="$2"
				shift 2
				;;
			--region-padding)
				region_padding="$2"
				shift 2
				;;
			--maf-threshold)
				maf_threshold="$2"
				shift 2
				;;
			--help | -h)
				usage
				exit 0
				;;
			*)
				echo "Unknown argument: $1" >&2
				usage
				exit 1
				;;
		esac
	done

	if [[ "${LE8_MRLINK2_WORKER:-0}" != 1 ]]; then
		exec "$python_bin" "$script_dir/c2.parallel.py" --worker-script "$script_dir/c2.cause.sh" \
			--jobs "$jobs" --cad-gwas "$cad_gwas" --outdir "$outdir" \
			--reference-bed "$ref_bed" --reference-pfile-dir "$ref_pfile_dir" --reference-pop "$ref_pfile_pop" \
			--reference-id-dir "$ref_id_dir" --reference-samples "$ref_samples" --mrlink2 "$mrlink2" \
			--p-threshold "$p_threshold" --region-padding "$region_padding" --maf-threshold "$maf_threshold"
	fi
	ref_pfile_pop="${ref_pfile_pop^^}"
	if [[ -z "$ref_pfile_pop" || "$ref_pfile_pop" == *[!A-Z0-9_-]* ]]; then
		echo "Invalid reference population: $ref_pfile_pop" >&2
		exit 2
	fi
	if [[ -z "$ref_id_dir" && -n "$ref_pfile_dir" ]]; then
		ref_id_dir="${ref_pfile_dir%/}/../id"
	fi
	if [[ -z "$ref_samples" && -n "$ref_pfile_dir" ]]; then
		ref_samples="${ref_pfile_dir%/}/../samples.txt"
	fi

	mkdir -p "$outdir" "$outdir/prepared" "$outdir/logs" "$outdir/tmp" "$outdir/results"
	status_file="$outdir/mrlink2.status.tsv"
	all_file="$outdir/mrlink2.all.tsv"
	complete_file="$outdir/mrlink2.complete"
	phe_f="${PHE_F:-/mnt/d/scripts/0f/phenotype.sh}"
	[[ -s "$phe_f" ]] || {
		echo "Missing shared phenotype functions: $phe_f" >&2
		exit 1
	}
	# shellcheck source=/mnt/d/scripts/0f/phenotype.sh
	source "$phe_f"
	declare -F match_GRCH >/dev/null 2>&1 || {
		echo "match_GRCH is missing from $phe_f" >&2
		exit 1
	}
	replace_run=false
	case "${LE8_REPLACE:-FALSE}" in TRUE | true | 1 | yes | YES) replace_run=true ;; esac
	if [[ -s "$complete_file" && "$replace_run" != true ]]; then
		echo "MR-link-2 already completed, skip MR-link-2 step: $status_file. Delete $complete_file to resume or delete $outdir to re-run."
		exit 0
	fi
	if [[ "$replace_run" == true ]]; then
		rm -f "$complete_file"
		# A changed exposure, region or LD reference must not be combined with an
		# old MR-link-2 --continue_analysis result or old standardized summary files.
		rm -rf "$outdir/prepared" "$outdir/results" "$outdir/tmp" "$outdir/reference_bed"
		mkdir -p "$outdir/prepared" "$outdir/results" "$outdir/tmp"
		: >"$all_file"
		printf "omics\ttrait\texposure\toutcome\tstatus\tout_prefix\tmessage\n" >"$status_file"
	else
		[[ -e "$all_file" ]] || : >"$all_file"
		[[ -s "$status_file" ]] || printf "omics\ttrait\texposure\toutcome\tstatus\tout_prefix\tmessage\n" >"$status_file"
	fi

	if [[ ! -s "$prep_py" ]]; then
		echo "Missing summary-stat standardizer: $prep_py" >&2
		exit 1
	fi
	if [[ -z "$jobs" || ! -s "$jobs" ]]; then
		echo "Missing job file: $jobs" >&2
		exit 1
	fi
	if [[ -z "$cad_gwas" || ! -s "$cad_gwas" ]]; then
		echo "Missing CAD GWAS file: $cad_gwas" >&2
		exit 1
	fi
	if [[ ! -s "$mrlink2" ]]; then
		echo "MR-link-2 script not found: $mrlink2" >&2
		exit 1
	fi
	have_ref_bed=true
	[[ -n "$ref_bed" ]] || have_ref_bed=false
	for ext in bed bim fam; do
		[[ -s "${ref_bed}.${ext}" ]] || have_ref_bed=false
	done
	if [[ "$have_ref_bed" != "true" ]]; then
		if [[ -z "$ref_pfile_dir" || ! -d "$ref_pfile_dir" ]]; then
			echo "Missing reference bed prefix (${ref_bed}.bed/.bim/.fam) and pfile directory: $ref_pfile_dir" >&2
			exit 1
		fi
	fi
	if ! command -v plink >/dev/null 2>&1; then
		echo "plink 1.9 is required by MR-link-2. Install plink in WSL and ensure it is in PATH." >&2
		exit 1
	fi
	if ! "$python_bin" -c 'import bitarray, duckdb, numpy, pandas, pyarrow, scipy' >/dev/null 2>&1; then
		echo "MR-link-2 Python dependencies are unavailable via $python_bin." >&2
		echo "Synchronize the environment first: conda env update -n le8 -f environment.yml" >&2
		exit 1
	fi

	region_chrom() {
		local region="$1"
		local chrom="${region%%:*}"
		chrom="${chrom#chr}"
		echo "$chrom"
	}

	bed_prefix_complete() {
		local prefix="$1"
		[[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]]
	}

	pfile_prefix_complete() {
		local prefix="$1"
		[[ -s "${prefix}.pgen" && -s "${prefix}.psam" ]] &&
			[[ -s "${prefix}.pvar" || -s "${prefix}.pvar.zst" ]]
	}

	prepare_population_keep() {
		local bed_dir="$1"
		local persistent_keep="${ref_id_dir%/}/${ref_pfile_pop}.id.2col"
		local keep_file="${bed_dir}/${ref_pfile_pop}.keep"
		local keep_tmp="${keep_file}.tmp.$$"

		if [[ -s "$persistent_keep" ]]; then
			if ! awk 'BEGIN{FS="[ \t]+"}
        {a=$1;b=$2;gsub(/\r/,"",a);gsub(/\r/,"",b);if(NF!=2||a!=b){bad=1;exit}}
        END{exit bad?2:0}' "$persistent_keep"; then
				echo "Invalid PLINK two-column keep file: $persistent_keep" >&2
				return 1
			fi
			echo "$persistent_keep"
			return 0
		fi
		if [[ ! -s "$ref_samples" ]]; then
			echo "Population ${ref_pfile_pop} requires $persistent_keep or a readable 1KG sample table: $ref_samples" >&2
			return 1
		fi
		if [[ ! -s "$keep_file" ]]; then
			if ! awk -v target="$ref_pfile_pop" 'BEGIN{FS="[ \t]+";OFS="\t"}
        NR==1{h=toupper($NF);gsub(/\r/,"",h);if(h!="SUPER_POP")exit 2;next}
        {pop=toupper($NF);gsub(/\r/,"",pop);if(pop==target)print $1,$1}' \
				"$ref_samples" >"$keep_tmp"; then
				echo "Invalid 1KG sample table (last column must be super_pop): $ref_samples" >&2
				rm -f "$keep_tmp"
				return 1
			fi
			if [[ ! -s "$keep_tmp" ]]; then
				echo "No ${ref_pfile_pop} samples found in: $ref_samples" >&2
				rm -f "$keep_tmp"
				return 1
			fi
			mv -f "$keep_tmp" "$keep_file"
		fi
		echo "$keep_file"
	}

	prepare_ref_bed() {
		local region="$1"
		if [[ "$have_ref_bed" == "true" ]]; then
			echo "$ref_bed"
			return 0
		fi
		local chrom pop_prefix global_prefix pop_bfile_prefix
		chrom="$(region_chrom "$region")"
		if bed_prefix_complete "${ref_pfile_dir}/${ref_pfile_pop}/chr${chrom}"; then
			printf '%s\n' "${ref_pfile_dir}/${ref_pfile_pop}/chr${chrom}"
			return 0
		fi
		pop_prefix="${ref_pfile_dir}/${ref_pfile_pop}.chr${chrom}"
		global_prefix="${ref_pfile_dir}/chr${chrom}"
		pop_bfile_prefix="${ref_pfile_dir%/}/../bfile/${ref_pfile_pop}/chr${chrom}"
		# Prefer the persistent population-specific BED set when it exists.  The
		# EUR references under 1kg/<build>/bfile/EUR have already been subset and
		# validated, so no per-run --keep conversion is needed.
		if bed_prefix_complete "$pop_bfile_prefix"; then
			echo "$pop_bfile_prefix"
			return 0
		fi
		if [[ "$ref_pfile_pop" != "ALL" ]] && bed_prefix_complete "$pop_prefix"; then
			echo "$pop_prefix"
			return 0
		fi
		if [[ "$ref_pfile_pop" == "ALL" ]] && bed_prefix_complete "$global_prefix"; then
			echo "$global_prefix"
			return 0
		fi

		local source_prefix="" source_kind="" use_keep=false
		if [[ "$ref_pfile_pop" != "ALL" ]] && pfile_prefix_complete "$pop_prefix"; then
			source_prefix="$pop_prefix"
			source_kind=pfile
		elif pfile_prefix_complete "$global_prefix"; then
			source_prefix="$global_prefix"
			source_kind=pfile
			[[ "$ref_pfile_pop" == "ALL" ]] || use_keep=true
		elif bed_prefix_complete "$global_prefix"; then
			source_prefix="$global_prefix"
			source_kind=bfile
			[[ "$ref_pfile_pop" == "ALL" ]] || use_keep=true
		fi
		if [[ -z "$source_prefix" ]]; then
			echo "Missing PLINK reference for chromosome ${chrom}: expected ${pop_prefix} or ${global_prefix} with .bed/.bim/.fam or .pgen/.pvar[.zst]/.psam" >&2
			return 1
		fi

		local bed_dir="$outdir/reference_bed"
		mkdir -p "$bed_dir"
		local keep_file=""
		if [[ "$use_keep" == true ]]; then
			keep_file="$(prepare_population_keep "$bed_dir")" || return 1
		fi

		local -a cache_args=(--source "$source_prefix" --kind "$source_kind"
			--cache-root "${MRLINK2_BED_CACHE:-$ref_pfile_dir/.mr-link2-bed-cache}")
		[[ -z "$keep_file" ]] || cache_args+=(--keep "$keep_file")
		"$python_bin" "$script_dir/c2.mr_link2.py" prepare-reference "${cache_args[@]}"

	}

	standardize_sumstats() {
		local input="$1"
		local output="$2"
		local default_n="${3:-nan}"
		local region="${4:-}"
		if [[ -s "$output" && "$replace_run" != true ]]; then
			if gzip -t "$output" 2>/dev/null; then
				return 0
			fi
			echo "Invalid cached gzip; rebuilding atomically: $output" >&2
		fi
		local output_tmp="${output}.tmp.${BASHPID:-$$}"
		rm -f "$output_tmp"
		if [[ -n "$region" ]]; then
			"$python_bin" "$script_dir/c2.parallel.py" prepare "$input" "$output_tmp" "$default_n" "$region" || {
				rm -f "$output_tmp"
				return 1
			}
		else
			"$python_bin" "$script_dir/c2.parallel.py" prepare "$input" "$output_tmp" "$default_n" || {
				rm -f "$output_tmp"
				return 1
			}
		fi
		gzip -t "$output_tmp" 2>/dev/null || {
			echo "Standardized output failed gzip validation: $output_tmp" >&2
			rm -f "$output_tmp"
			return 1
		}
		mv -f "$output_tmp" "$output"
	}

	total_jobs=$(awk -F '\t' 'NR > 1 && $2 != "" {n++} END {print n+0}' "$jobs")
	job_index=0
	tail -n +2 "$jobs" | while IFS=$'\t' read -r omics trait exposure outcome region rest; do
		[[ -z "${trait:-}" ]] && continue
		job_index=$((job_index + 1))
		outcome_file="${outcome:-$cad_gwas}"
		[[ -z "$outcome_file" ]] && outcome_file="$cad_gwas"
		if [[ ! -s "$exposure" ]]; then
			printf "%s\t%s\t%s\t%s\tnot_run\t\tmissing exposure file\n" "$omics" "$trait" "$exposure" "$outcome_file" >>"$status_file"
			continue
		fi
		safe_trait=$(echo "${omics}_${trait}" | sed 's/[^A-Za-z0-9_.-]/_/g')
		echo "MR-link-2 [$job_index/$total_jobs]: $safe_trait"
		if awk -F '\t' -v om="$omics" -v tr="$trait" 'NR > 1 && $1 == om && $2 == tr && ($5 == "ok" || $5 == "no_estimate") { found=1 } END { exit !found }' "$status_file"; then
			echo "MR-link-2 [$job_index/$total_jobs]: $safe_trait already completed, skip."
			continue
		fi
		exp_std="$outdir/prepared/${safe_trait}.exposure.mrlink2.tsv.gz"
		out_std="$outdir/prepared/${safe_trait}.outcome.mrlink2.tsv.gz"
		exp_matched="$outdir/prepared/${safe_trait}.exposure.ref.tsv"
		out_matched="$outdir/prepared/${safe_trait}.outcome.ref.tsv"
		out_prefix="$outdir/results/${safe_trait}.mrlink2"
		log_file="$outdir/logs/${safe_trait}.log"
		err_file="$outdir/logs/${safe_trait}.err"
		if ! standardize_sumstats "$exposure" "$exp_std" "${MRLINK2_EXPOSURE_N:-nan}" "$region" >"$log_file.prepare" 2>"$err_file.prepare"; then
			status=not_run; grep -q not_run_missing_N "$err_file.prepare" && status=not_run_missing_N
			printf "%s\t%s\t%s\t%s\t$status\t%s\tstandardization failed; see %s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$err_file.prepare" >>"$status_file"
			continue
		fi
		if ! standardize_sumstats "$outcome_file" "$out_std" "${MRLINK2_OUTCOME_N:-nan}" "$region" >>"$log_file.prepare" 2>>"$err_file.prepare"; then
			status=not_run; grep -q not_run_missing_N "$err_file.prepare" && status=not_run_missing_N
			printf "%s\t%s\t%s\t%s\t$status\t%s\toutcome standardization failed; see %s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$err_file.prepare" >>"$status_file"
			continue
		fi
		ref_for_trait="$(prepare_ref_bed "$region" 2>>"$err_file.prepare" || true)"
		if [[ -z "$ref_for_trait" ]]; then
			printf "%s\t%s\t%s\t%s\tnot_run\t%s\treference conversion failed; see %s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$err_file.prepare" >>"$status_file"
			continue
		fi
		if ! "$python_bin" "$script_dir/c2.parallel.py" match "$ref_for_trait" "$exp_std" "$exp_matched" "$outdir/prepared/${safe_trait}.exposure.match.tsv" >>"$log_file.prepare" 2>>"$err_file.prepare" ||
            ! "$python_bin" "$script_dir/c2.parallel.py" match "$ref_for_trait" "$out_std" "$out_matched" "$outdir/prepared/${safe_trait}.outcome.match.tsv" >>"$log_file.prepare" 2>>"$err_file.prepare"; then
			printf "%s\t%s\t%s\t%s\tnot_run\t%s\treference-ID matching failed; see %s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$err_file.prepare" >>"$status_file"
			continue
		fi
		if [[ "$(wc -l <"$exp_matched")" -le 1 || "$(wc -l <"$out_matched")" -le 1 ]]; then
			printf "%s\t%s\t%s\t%s\tnot_run\t%s\tno variants matched the build-specific reference\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" >>"$status_file"
			continue
		fi
		cmd=("$python_bin" "$mrlink2"
			--reference_bed "$ref_for_trait"
			--sumstats_exposure "$exp_matched"
			--sumstats_outcome "$out_matched"
			--out "$out_prefix"
			--p_threshold "$p_threshold"
			--region_padding "$region_padding"
			--maf_threshold "$maf_threshold"
			--max_correlation "$max_correlation"
			--verbose "$verbose")
		# c2.cause.R records the biologically defined cis/local region in the jobs
		# TSV. Pass that region through so MR-link-2 does not replace it by a second
		# clumping-derived region.
		if [[ -n "${region:-}" ]]; then
			cmd+=(--prespecified_regions "$region")
		fi
		if [[ "$continue_analysis" == "true" || "$continue_analysis" == "TRUE" || "$continue_analysis" == "1" ]]; then
			cmd+=(--continue_analysis)
		fi
		if "${cmd[@]}" --tmp "$outdir/tmp/${safe_trait}" >"$log_file" 2>"$err_file"; then
			if [[ -s "$out_prefix" && "$(wc -l <"$out_prefix")" -gt 1 ]]; then
				printf "%s\t%s\t%s\t%s\tok\t%s\tcompleted\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" >>"$status_file"
				echo "MR-link-2 [$job_index/$total_jobs]: $safe_trait completed (ok)."
			else
				detail="MR-link-2 exited normally but produced no regional estimate"
				[[ -s "${out_prefix}_no_estimate" ]] && detail="${detail}; see ${out_prefix}_no_estimate"
				printf "%s\t%s\t%s\t%s\tno_estimate\t%s\t%s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$detail" >>"$status_file"
				echo "MR-link-2 [$job_index/$total_jobs]: $safe_trait completed (no estimate)."
			fi
		else
			printf "%s\t%s\t%s\t%s\tfailed\t%s\tMR-link-2 failed; see %s\n" "$omics" "$trait" "$exposure" "$outcome_file" "$out_prefix" "$err_file" >>"$status_file"
			echo "MR-link-2 [$job_index/$total_jobs]: $safe_trait failed; see $err_file" >&2
		fi
	done

	# Combine only the exact 11-column MR-link-2 result files.
	rebuild_file="${all_file}.tmp"
	: >"$rebuild_file"
	first_result=true
	while IFS=$'\t' read -r omics trait out_prefix; do
		[[ -z "${out_prefix:-}" || ! -s "$out_prefix" ]] && continue
		if [[ "$first_result" == true ]]; then
			awk -v OFS='\t' -v om="$omics" -v tr="$trait" 'NR==1{print "omics","trait",$0; next} NR>1{print om,tr,$0}' "$out_prefix" >>"$rebuild_file"
			first_result=false
		else
			awk -v OFS='\t' -v om="$omics" -v tr="$trait" 'NR>1{print om,tr,$0}' "$out_prefix" >>"$rebuild_file"
		fi
	done < <(awk -F '\t' 'NR>1 && $5=="ok" {key=$1 FS $2; if(!seen[key]++) print $1 FS $2 FS $6}' "$status_file")
	mv -f "$rebuild_file" "$all_file"

	if ! awk -F '\t' -v expected="$total_jobs" 'NR>1 {n++;if($5!="ok" && $5!="no_estimate")bad=1} END{exit (bad || n!=expected)}' "$status_file"; then
		rm -f "$complete_file"
		echo "MR-link-2 incomplete: failed/not_run tasks must be retried" >&2
		exit 2
	fi
	printf "completed\t%s\n" "$(date -Is)" >"$complete_file"

	echo "MR-link-2 status: $status_file"
	echo "MR-link-2 combined results: $all_file"
)

case "${1:---help}" in
	index-cis)
		shift
		index_cis_main "$@"
		;;
	mr-link2)
		shift
		mr_link2_main "$@"
		;;
	-h | --help)
		echo "Usage: c2.cause.sh {mr-link2|index-cis} [options]"
		echo "Run either subcommand with --help for its options."
		;;
	*)
		echo "Unknown C2 command: $1" >&2
		exit 2
		;;
esac

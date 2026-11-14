#!/bin/bash -l

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    liftOver  Convert genotype coordinates
    pgs       Calculate polygenic scores
    snp       Extract selected variants
    apoe      Derive APOE genotypes

Main options:
    --dir-in DIR --dir-out DIR  Input/output
    --dir-gen DIR --data NAME  Genotype dataset
    --chain-file FILE  Liftover chain
    --chr CHR  Process only chrCHR.pvar (e.g. 20); default: all chromosomes
    --ref-fasta FILE  Indexed GRCh38 FASTA [${DIR_GEN}/1hgp/38/GRCH38.fasta.gz]
    liftOver preserves row count/order and all fields except POS; failures get POS=-1.
    Original .pgen/.psam are reused through relative symlinks.
    --input FILE --label NAME  PGS input and label
    --dir-pgs-pfile DIR  PGS genotypes [${DIR_GEN}/ukb/38/imp]
    --submit local|hpc|auto  Auto-detect HPC/local using shared 0f/0os.sh [auto]; optional override
    --write-cmd / --dry-run  Write the .cmd only; --run-cmd FALSE is equivalent
    --queue NAME --slots N --memory-gb N --walltime HHH:MM  LSF resources
    --threads N  Analysis threads [4]; slots round up by the shared CPU/RAM policy
    --dir-jobs DIR  Command/log directory [OUTPUT_DIR/logs]
    --wait TRUE|FALSE  Wait for verified LSF completion [FALSE]
    --background  Legacy alias for --submit hpc
    UKB_PYTHON  Python override [SHARED_PYTHON from 0f/0os.sh]
    -h, --help  Show options and examples below

Examples:
    source ../0f/0os.sh
    ./gen.sh pgs --input "$DIR_GWAS/bmi.ref" --label bmi --write-cmd
    ./gen.sh pgs --input "$DIR_GWAS/bmi.ref" --label bmi
    ./gen.sh snp --dir-gen "$DIR_GEN/ukb/37" --data typ --snp-file "$DIR_FILES/Wang_MY.bed" --snp-file-type bed --dir-out "$ANALYSIS_ROOT/ukb" --label Wang
PUBLIC_ENTRY_HELP
    exit 0
fi
# Plan every module as a .cmd before local execution or LSF submission.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/f/jobs.sh"
ukb_job_init "$@" || exit $?
set -- "${UKB_ARGS[@]}"
source "$shared_dir/cli_dirs.sh"
script_directory_args "$SCRIPT_DIR/gen.sh" "$@"
set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
worker_args=("$@")
set +H
set -eo pipefail

usage() {
	bash "${BASH_SOURCE[0]}" --help
}

init_gen_env() {
	dir0=${DIR0:-${DIR_BIG}}
	local phe_f="$shared_dir/phenotype.sh"

	if [[ -f "$phe_f" ]]; then
		PHE_F=$phe_f
	elif [[ -f "$dir0/scripts/0f/phenotype.sh" ]]; then
		PHE_F=$dir0/scripts/0f/phenotype.sh
	else
		echo "ERROR: cannot find phenotype.sh beside gen.sh or under $dir0/scripts/0f" >&2
		exit 1
	fi
	source "$PHE_F"
}

to_wsl_path() {
	local path=$1
	if [[ "$path" =~ ^[A-Za-z]:[\\/].* ]]; then
		command -v wslpath >/dev/null 2>&1 || {
			echo "ERROR: Windows path supplied but wslpath is unavailable: $path" >&2
			exit 1
		}
		wslpath -u "$path"
	else
		printf '%s\n' "$path"
	fi
}

pvar_stream() {
	cat -- "$1"
}

liftover_output_is_complete() {
	# The log is written last and acts as the completion marker.  Validate the
	# marker, dependencies, links, and PVAR rows before resuming past it.
	local input_pvar=$1 input_pgen=$2 input_psam=$3 output_pvar=$4 output_log=$5
	local output_dir output_name output_stem
	output_dir=$(dirname "$output_pvar")
	output_name=$(basename "$output_pvar")
	output_stem=${output_name%.pvar}

	[[ -s "$output_pvar" && -s "$output_log" ]] || return 1
	[[ "$output_pvar" -nt "$input_pvar" ]] || return 1
	[[ "$output_pvar" -nt "$chain_file" ]] || return 1
	[[ "$output_pvar" -nt "$ref_fasta" && "$output_pvar" -nt "$ref_fasta.fai" ]] || return 1
	[[ $ref_fasta != *.gz && $ref_fasta != *.bgz && $ref_fasta != *.bgzf ]] || [[ -s $ref_fasta.gzi && "$output_pvar" -nt "$ref_fasta.gzi" ]] || return 1
	[[ -L "$output_dir/$output_stem.pgen" && "$output_dir/$output_stem.pgen" -ef "$input_pgen" ]] || return 1
	[[ -L "$output_dir/$output_stem.psam" && "$output_dir/$output_stem.psam" -ef "$input_psam" ]] || return 1

	awk -F '\t' -v expected_input="$input_pvar" -v expected_output="$output_pvar" '
        $1=="input" {
            input_count++
            input_ok=($2==expected_input)
        }
        $1=="output" {
            output_count++
            output_ok=($2==expected_output)
        }
        $1=="total_variants" && $2 ~ /^[0-9]+$/ {
            total_count++
            total=$2+0
        }
        $1=="unique_forward_same_chr_length" && $2 ~ /^[0-9]+$/ {
            candidate_count++
            candidate=$2+0
        }
        $1=="verified_ref" && $2 ~ /^[0-9]+$/ {
            ok_count++
            ok=$2+0
        }
        $1=="failed_pos_minus_1" && $2 ~ /^[0-9]+$/ {
            fail_count++
            fail=$2+0
        }
        END {
            valid=(input_count==1 && input_ok && output_count==1 && output_ok &&
                   total_count==1 && candidate_count==1 && ok_count==1 && fail_count==1 &&
                   total>0 && candidate<=total && ok<=candidate && fail==total-ok)
            exit !valid
        }
    ' "$output_log" || return 1

	awk -v expected="$(awk -F '\t' '$1=="total_variants" {print $2}' "$output_log")" '
        /^##/ {next}
        /^#CHROM[\t ]/ {header=1; next}
        !/^#/ {if (NF<5) bad=1; n++}
        END {exit !(header && !bad && n==expected)}
    ' "$output_pvar"
}

run_liftover() {
	if [[ -n "${bim_file:-}" ]]; then
		python3 "$SCRIPT_DIR/f/gen.py" lift-bim "$(to_wsl_path "$bim_file")" \
			"$(to_wsl_path "${dir_out:?ERROR: --dir-out required}")" \
			"$(to_wsl_path "${chain_file:?ERROR: --chain-file required}")" \
			"$(to_wsl_path "${ref_fasta:-${DIR_GEN}/1hgp/38/GRCH38.fasta.gz}")"
		return
	fi
	# Keep every PVAR row in its original position.  Only POS changes; rows that
	# cannot be mapped and verified against the target reference receive POS=-1.
	set -o pipefail
	: "${dir_in:?ERROR: liftOver requires --dir-in}"
	: "${dir_out:?ERROR: liftOver requires --dir-out}"
	: "${chain_file:?ERROR: liftOver requires --chain-file}"

	dir_in=$(to_wsl_path "$dir_in")
	dir_out=$(to_wsl_path "$dir_out")
	chain_file=$(to_wsl_path "$chain_file")
	ref_fasta=$(to_wsl_path "${ref_fasta:-${DIR_GEN}/1hgp/38/GRCH38.fasta.gz}")

	[[ -d "$dir_in" ]] || {
		echo "ERROR: input directory not found: $dir_in" >&2
		exit 1
	}
	[[ -f "$chain_file" ]] || {
		echo "ERROR: chain file not found: $chain_file" >&2
		exit 1
	}
	[[ -f "$ref_fasta" ]] || {
		echo "ERROR: GRCh38 FASTA not found: $ref_fasta" >&2
		exit 1
	}
	[[ -f "$ref_fasta.fai" ]] || {
		echo "ERROR: FASTA index not found: $ref_fasta.fai" >&2
		exit 1
	}
	[[ $ref_fasta != *.gz && $ref_fasta != *.bgz && $ref_fasta != *.bgzf ]] || [[ -s $ref_fasta.gzi ]] || {
		echo "ERROR: BGZF index not found: $ref_fasta.gzi" >&2
		exit 1
	}

	local lift_over=${LIFTOVER_BIN:-}
	if [[ -z "$lift_over" ]]; then
		if command -v liftOver >/dev/null 2>&1; then
			lift_over=$(command -v liftOver)
		elif [[ -x ${DIR_SOFTWARE}/bin/liftOver ]]; then
			lift_over=${DIR_SOFTWARE}/bin/liftOver
		else
			echo "ERROR: liftOver not found (set LIFTOVER_BIN if needed)" >&2
			exit 1
		fi
	fi
	if ! command -v bedtools >/dev/null 2>&1 && [[ -x $DIR_WORK/.conda/envs/gu/bin/bedtools ]]; then
		export PATH="$DIR_WORK/.conda/envs/gu/bin:$PATH"
	fi
	command -v bedtools >/dev/null 2>&1 || {
		echo "ERROR: bedtools not found" >&2
		exit 1
	}
	command -v realpath >/dev/null 2>&1 || {
		echo "ERROR: realpath not found" >&2
		exit 1
	}

	mkdir -p "$dir_out"
	dir_in=$(cd "$dir_in" && pwd -P)
	dir_out=$(cd "$dir_out" && pwd -P)
	chain_file=$(cd "$(dirname "$chain_file")" && pwd -P)/$(basename "$chain_file")
	ref_fasta=$(cd "$(dirname "$ref_fasta")" && pwd -P)/$(basename "$ref_fasta")
	[[ "$dir_in" != "$dir_out" ]] || {
		echo "ERROR: --dir-in and --dir-out must differ" >&2
		exit 1
	}

	local -a pvars
	mapfile -t pvars < <(find "$dir_in" -maxdepth 1 -type f -name 'chr*.pvar' -printf '%p\n' | sort -V)
	((${#pvars[@]} > 0)) || {
		echo "ERROR: no chr*.pvar files found in $dir_in" >&2
		exit 1
	}

	if [[ -n ${lift_chr:-} ]]; then
		pvars=("$dir_in/chr$lift_chr.pvar")
		[[ -f ${pvars[0]} ]] || { echo "ERROR: missing ${pvars[0]}" >&2; exit 1; }
	fi

	local lift_work_dir=
	cleanup_liftover_work() {
		if [[ -n "$lift_work_dir" && -d "$lift_work_dir" ]]; then
			find "$lift_work_dir" -maxdepth 1 -type f -delete
			rmdir "$lift_work_dir" 2>/dev/null || true
		fi
	}
	trap cleanup_liftover_work EXIT
	trap 'exit 130' INT TERM

	local pvar name stem pgen psam out_pvar out_log
	local source_bed mapped_bed unmapped_bed sorted_bed candidate_bed verified_map out_tmp out_log_tmp
	local n_total n_candidate n_ok n_fail rel target
	for pvar in "${pvars[@]}"; do
		name=$(basename "$pvar")
		stem=${name%.pvar}
		pgen="$dir_in/$stem.pgen"
		psam="$dir_in/$stem.psam"
		[[ -e "$pgen" ]] || {
			echo "ERROR: missing companion file: $pgen" >&2
			exit 1
		}
		[[ -e "$psam" ]] || {
			echo "ERROR: missing companion file: $psam" >&2
			exit 1
		}

		out_pvar="$dir_out/$name"
		mkdir -p "$dir_out/logs" "$dir_out/tmp"
		out_log="$dir_out/logs/$name.liftOver.log"
		# Completed legacy outputs remain reusable without rewriting their receipt.
		if [[ ! -e $out_log && -e $dir_out/$name.liftOver.log ]]; then
			out_log="$dir_out/$name.liftOver.log"
		fi
		if liftover_output_is_complete "$pvar" "$pgen" "$psam" "$out_pvar" "$out_log"; then
			echo "SKIP $name (complete)"
			continue
		fi
		if [[ -e "$out_pvar" || -e "$out_log" ]]; then
			echo "REBUILD $name (existing result is incomplete or stale)"
		fi

		lift_work_dir="$dir_out/tmp/.liftOver.${stem}.$$"
		mkdir "$lift_work_dir"
		source_bed="$lift_work_dir/source.bed"
		mapped_bed="$lift_work_dir/mapped.bed"
		unmapped_bed="$lift_work_dir/unmapped.bed"
		sorted_bed="$lift_work_dir/mapped.sorted.bed"
		candidate_bed="$lift_work_dir/candidate.bed"
		verified_map="$lift_work_dir/verified.tsv"
		out_tmp="$lift_work_dir/$name"
		out_log_tmp="$lift_work_dir/$name.liftOver.log"

		echo "LIFT $name"
		# The BED name carries row number, expected chromosome, and original REF.
		# It survives liftOver and avoids loading millions of variants into memory.
		n_total=$(pvar_stream "$pvar" | awk -v bed="$source_bed" '
            BEGIN { OFS="\t" }
            /^#/ { next }
            {
                row++
                chr=$1
                sub(/^chr/, "", chr)
                if (chr==23) chr="X"
                else if (chr==24) chr="Y"
                else if (chr==25) chr="XY"
                else if (chr==26 || chr=="M") chr="MT"
                target_chr="chr" chr
                ref=toupper($4)
                if ($2 ~ /^[0-9]+$/ && $2>0 && ref ~ /^[ACGTN]+$/) {
                    start=$2-1
                    print target_chr,start,start+length(ref),row "|" target_chr "|" ref,0,"+" > bed
                }
            }
            END { print row+0 }
        ')
		((n_total > 0)) || {
			echo "ERROR: no variant rows in $pvar" >&2
			exit 1
		}
		: >"$mapped_bed"
		: >"$unmapped_bed"
		: >"$sorted_bed"
		if [[ -s "$source_bed" ]]; then
			"$lift_over" "$source_bed" "$chain_file" "$mapped_bed" "$unmapped_bed" >/dev/null
			LC_ALL=C sort -k4,4n -k2,2n "$mapped_bed" >"$sorted_bed"
		fi

		# Accept exactly one forward-strand mapping to the expected chromosome,
		# with the same interval length.  REF is verified in the next pass.
		: >"$candidate_bed"
		awk -v out="$candidate_bed" '
            BEGIN { OFS="\t" }
            function load_record(    a,n) {
                n=split($4,a,"|")
                id=a[1]+0
                expected_chr=a[2]
                ref=a[3]
                tchr=$1; tstart=$2; tend=$3
                ok=(n==3 && $6=="+" && $1==expected_chr && $2>=0 && $3>$2 && $3-$2==length(ref))
            }
            function emit_group() {
                if (count==1 && ok) {
                    print tchr,tstart,tend,id "|" (tstart+1) "|" ref > out
                }
            }
            NR==1 { load_record(); count=1; next }
            { split($4,current,"|") }
            current[1]+0==id { count++; next }
            { emit_group(); load_record(); count=1 }
            END { if (NR) emit_group() }
        ' "$sorted_bed"

		: >"$verified_map"
		if [[ -s "$candidate_bed" ]]; then
			bedtools getfasta -fi "$ref_fasta" -bed "$candidate_bed" -nameOnly -tab 2>/dev/null |
				awk -F '\t' 'BEGIN{OFS="\t"} {
                    n=split($1,a,"|")
                    if (n==3 && toupper($2)==toupper(a[3])) print a[1],a[2]
                }' >"$verified_map"
		fi
		n_candidate=$(wc -l <"$candidate_bed")
		n_ok=$(wc -l <"$verified_map")
		n_fail=$((n_total - n_ok))

		# verified_map is sorted by row number.  Consume it as a stream while
		# rereading the PVAR, so memory use stays bounded for imputed chromosomes.
		pvar_stream "$pvar" | awk -v map="$verified_map" '
            BEGIN {
                OFS="\t"
                if ((getline line < map)>0) { split(line,a,"\t"); next_id=a[1]+0; next_pos=a[2]+0 }
                else next_id=0
            }
            /^#/ { print; next }
            {
                row++
                if (row==next_id) {
                    $2=next_pos
                    if ((getline line < map)>0) { split(line,a,"\t"); next_id=a[1]+0; next_pos=a[2]+0 }
                    else next_id=0
                } else {
                    $2=-1
                }
                print
            }
            END {
                close(map)
                if (row==0 || next_id!=0) exit 42
            }
        ' >"$out_tmp"
		[[ $(cat -- "$out_tmp" | awk '!/^#/{n++} END{print n+0}') -eq "$n_total" ]] || {
			echo "ERROR: output row-count check failed for $name" >&2
			exit 1
		}
		mv -f "$out_tmp" "$out_pvar"

		local companion_name
		for target in "$pgen" "$psam"; do
			companion_name=$(basename "$target")
			if [[ -e "$dir_out/$companion_name" && ! -L "$dir_out/$companion_name" ]]; then
				echo "ERROR: refusing to replace non-symlink: $dir_out/$companion_name" >&2
				exit 1
			fi
			rel=$(realpath --relative-to="$dir_out" "$target")
			ln -sfn -- "$rel" "$dir_out/$companion_name"
		done

		{
			printf 'input\t%s\n' "$pvar"
			printf 'output\t%s\n' "$out_pvar"
			printf 'total_variants\t%s\n' "$n_total"
			printf 'unique_forward_same_chr_length\t%s\n' "$n_candidate"
			printf 'verified_ref\t%s\n' "$n_ok"
			printf 'failed_pos_minus_1\t%s\n' "$n_fail"
			printf 'pgen_link\t%s\n' "$(readlink "$dir_out/$stem.pgen")"
			printf 'psam_link\t%s\n' "$(readlink "$dir_out/$stem.psam")"
		} >"$out_log_tmp"
		mv -f "$out_log_tmp" "$out_log"
		echo "OK $name: total=$n_total lifted=$n_ok POS=-1=$n_fail"

		cleanup_liftover_work
		lift_work_dir=
	done
	trap - EXIT INT TERM
	echo "OK: lifted PVARs and relative pfile links written to: $dir_out"
}

run_pgs() {
	init_gen_env
	: "${pgs_input:?ERROR: pgs requires --input FILE}"
	: "${label:?ERROR: pgs requires --label LABEL}"
	pgs_input=$(to_wsl_path "$pgs_input")
	pgs_pfile_dir=$(to_wsl_path "${pgs_pfile_dir:-${DIR_GEN}/ukb/38/imp}")
	dir_out=$(to_wsl_path "${dir_out:-$ANALYSIS_ROOT/ukb/gen/pgs/$label}")
	pgs_threads=${pgs_threads:-$UKB_THREADS}
	[[ -s "$pgs_input" ]] || {
		echo "ERROR: PGS score input is missing/empty: $pgs_input" >&2
		exit 1
	}
	[[ "$pgs_threads" =~ ^[1-9][0-9]*$ ]] || {
		echo "ERROR: --pgs-threads must be a positive integer" >&2
		exit 1
	}
	declare -F pgs_plink_calc >/dev/null 2>&1 || {
		echo "ERROR: pgs_plink_calc is missing from $PHE_F" >&2
		exit 1
	}

	mkdir -p "$dir_out/jobs" "$dir_out/tmp" "$dir_out/qc" "$dir_out/logs"
	local cmd_file="$dir_out/jobs/$label.pgs.worker.sh"
	{
		printf '%s\n' '#!/bin/bash -l' 'set +H' 'set -euo pipefail'
		printf 'score_file=%q\n' "$pgs_input"
		printf 'pfile_dir=%q\n' "$pgs_pfile_dir"
		printf 'out_dir=%q\n' "$dir_out"
		printf 'label=%q\n' "$label"
		printf 'threads=%q\n' "$pgs_threads"
		printf 'phe_f=%q\n' "$PHE_F"
		cat <<'EOF2'

source "$phe_f"
declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from $phe_f" >&2; exit 1; }
collect_pgs_logs() {
    local status=$? log
    trap - EXIT
    for log in "$out_dir/tmp/$label".chr*.log; do
        [[ ! -f $log ]] || mv -f -- "$log" "$out_dir/logs/"
    done
    exit "$status"
}
trap collect_pgs_logs EXIT
pgs_plink_calc \
    --input "$score_file" \
    --pfile-dir "$pfile_dir" \
    --output "$out_dir/$label.pgs.txt" \
    --label "$label" \
    --work-dir "$out_dir/tmp" \
    --status-file "$out_dir/qc/$label.pgs.status.tsv" \
    --threads "$threads"
echo "OK: PGS written to: $out_dir/$label.pgs.txt"
EOF2
	} >"$cmd_file"
	chmod +x "$cmd_file"

	bash "$cmd_file"
}

run_snp() {


	# Extract SNP data
	: "${dir_gen:?ERROR: snp requires --dir-gen}"
	: "${snp_data:?ERROR: snp requires --data typ|imp}"
	: "${snp_file:?ERROR: snp requires --snp-file}"
	: "${snp_file_type:?ERROR: snp requires --snp-file-type bed|snp}"
	: "${dir_out:?ERROR: snp requires --dir-out}"
	: "${label:?ERROR: snp requires --label}"
	[[ "$snp_data" == typ || "$snp_data" == imp ]] || {
		echo "ERROR: --data must be typ or imp" >&2
		exit 1
	}
	[[ "$snp_file_type" == bed || "$snp_file_type" == snp ]] || {
		echo "ERROR: --snp-file-type must be bed or snp" >&2
		exit 1
	}
	if [[ "$snp_file_type" == bed ]]; then
		flank=${flank:-10000}
		[[ "$flank" =~ ^[0-9]+$ ]] || {
			echo "ERROR: --flank must be a nonnegative integer" >&2
			exit 1
		}
	elif [[ "${flank_set:-N}" == Y ]]; then
		echo "ERROR: --flank is only valid with --snp-file-type bed" >&2
		exit 1
	fi
	[[ -d "$dir_gen/$snp_data" ]] || {
		echo "ERROR: missing data directory: $dir_gen/$snp_data" >&2
		exit 1
	}
	[[ -f "$snp_file" ]] || {
		echo "ERROR: missing SNP file: $snp_file" >&2
		exit 1
	}

	dir_gen=$(cd "$dir_gen" && pwd)
	snp_file=$(cd "$(dirname "$snp_file")" && pwd)/$(basename "$snp_file")
	mkdir -p "$dir_out"
	dir_out=$(cd "$dir_out" && pwd)
	work_dir="$dir_out/tmp/${label}.snp-work"
	mkdir -p "$work_dir"

	if [[ "$snp_file_type" == bed ]]; then
		range_file="$work_dir/$label.flank${flank}.bed"
		awk -v flank="$flank" 'BEGIN{OFS="\t"} /^[[:space:]]*#/ || NF<3 {next} {s=$2-flank; if(s<0)s=0; print $1,s,$3+flank,(NF>=4?$4:".")}' "$snp_file" >"$range_file"
		[[ -s "$range_file" ]] || {
			echo "ERROR: no BED intervals found in $snp_file" >&2
			exit 1
		}
		extract_args="--extract bed0 \"$range_file\""
		mapfile -t chromosomes < <(awk '{c=$1; sub(/^chr/,"",c); print c}' "$range_file" | sort -uV)
	else
		extract_file="$work_dir/$label.extract.snp"
		awk 'BEGIN{FS="[[:space:]]+"}
             /^[[:space:]]*#/ || /^[[:space:]]*$/ {next}
             !seen && toupper($1) ~ /^(SNP|RSID|VARIANT|ID)$/ {seen=1; next}
             {seen=1; print $1}' "$snp_file" >"$extract_file"
		[[ -s "$extract_file" ]] || {
			echo "ERROR: no variant IDs found in $snp_file" >&2
			exit 1
		}
		if [[ $(sort "$extract_file" | uniq -d | wc -l) -gt 0 ]]; then
			echo "ERROR: duplicate variant IDs found in first column of $snp_file" >&2
			sort "$extract_file" | uniq -d | head -20 >&2
			exit 1
		fi
		extract_args="--extract \"$extract_file\""
		if [[ "$snp_data" == typ ]]; then
			mapfile -t chromosomes < <(find "$dir_gen/typ" -maxdepth 1 -type f -name 'chr*.bed' -printf '%f\n' | sed -E 's/^chr(.+)\.bed$/\1/' | sort -uV)
		else
			mapfile -t chromosomes < <(find "$dir_gen/imp" -maxdepth 1 -type f -name 'chr*.pgen' -printf '%f\n' | sed -E 's/^chr(.+)\.pgen$/\1/' | sort -uV)
		fi
	fi
	((${#chromosomes[@]} > 0)) || {
		echo "ERROR: no chromosome data files found" >&2
		exit 1
	}

	for chr in "${chromosomes[@]}"; do
		if [[ "$snp_data" == typ ]]; then
			[[ -f "$dir_gen/typ/chr$chr.bed" && -f "$dir_gen/typ/chr$chr.bim" && -f "$dir_gen/typ/chr$chr.fam" ]] || {
				echo "ERROR: incomplete typ files for chr$chr" >&2
				exit 1
			}
		else
			[[ -f "$dir_gen/imp/chr$chr.pgen" && -f "$dir_gen/imp/chr$chr.pvar" && -f "$dir_gen/imp/chr$chr.psam" ]] || {
				echo "ERROR: incomplete imp files for chr$chr" >&2
				exit 1
			}
		fi
	done

	cmd_file="$work_dir/$label.worker.sh"
	{
		cat <<EOF2
#!/bin/bash -l
set +H
set -e

scratch=\$(mktemp -d "$work_dir/run.XXXXXX")
trap 'rm -rf -- "\$scratch"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP
cd "\$scratch"
rm -f merge-list.txt
for chr in ${chromosomes[*]}; do
    if [[ "$snp_data" == "typ" ]]; then
        input_args="--bfile \"$dir_gen/typ/chr\$chr\""
    else
        input_args="--pfile \"$dir_gen/imp/chr\$chr\""
    fi
    set +e
    eval plink2 --threads "$UKB_THREADS" --memory "$UKB_PLINK_MEMORY_MB" \$input_args $extract_args --make-pgen --out "chr\$chr"
    plink_status=\$?
    set -e
    if (( plink_status != 0 )); then
        if grep -Eq 'No variants remaining|No variants loaded|0 variants remaining' "chr\$chr.log"; then
            echo "WARN: no variants extracted from chr\$chr" >&2
            rm -f "chr\$chr.pgen" "chr\$chr.pvar" "chr\$chr.psam"
            continue
        fi
        echo "ERROR: PLINK extraction failed for chr\$chr (exit \$plink_status)" >&2
        exit "\$plink_status"
    fi
    if awk '!/^#/{found=1} END{exit !found}' "chr\$chr.pvar"; then
        echo "chr\$chr" >> merge-list.txt
    else
        rm -f "chr\$chr.pgen" "chr\$chr.pvar" "chr\$chr.psam"
        echo "WARN: no variants extracted from chr\$chr" >&2
    fi
done
[[ -s merge-list.txt ]] || { echo "ERROR: no variants were extracted" >&2; exit 1; }
if [[ \$(wc -l < merge-list.txt) -eq 1 ]]; then
    merged=\$(head -1 merge-list.txt)
else
    merged=merged
    plink2 --threads "$UKB_THREADS" --memory "$UKB_PLINK_MEMORY_MB" --pmerge-list merge-list.txt --make-pgen --out "\$merged"
    while read -r part; do rm -f "\$part.pgen" "\$part.pvar" "\$part.psam"; done < merge-list.txt
fi
python3 "$SCRIPT_DIR/f/gen.py" export --output "$dir_out/$label.raw.gz" -- plink2 --threads "$UKB_THREADS" --memory "$UKB_PLINK_MEMORY_MB" --pfile "\$merged" --maj-ref --export A
EOF2
	} >"$cmd_file"
	chmod +x "$cmd_file"

	bash "$cmd_file"
	echo "OK: SNP data written to: $dir_out/$label.raw.gz"
}

run_apoe() {


	# Extract APOE haplotype
	# e1: C-T; e2: T-T; e3: T-C; e4: C-C
	dir0=${DIR_WORK}
	dirout=$(realpath -m -- "${dir_out:-$ANALYSIS_ROOT/ukb/gen/apoe}")
	dir_gen=${dir_gen:-$DIR_GEN/ukb/37}
	mkdir -p "$dirout/tmp/apoe" "$dirout/jobs" "$dirout/logs"

	echo -e "rs429358\nrs7412" >"$dirout/tmp/apoe/apoe.snps"

	cat >"$dirout/jobs/apoe.worker.sh" <<EOF2
#!/bin/bash -l
set +H
set -e

cd "$dirout/tmp/apoe"
plink2 --threads "$UKB_THREADS" --memory "$UKB_PLINK_MEMORY_MB" --pfile "$dir_gen/hap/chr19" --extract apoe.snps --export vcf id-paste=iid bgz --out apoe
tabix -f apoe.vcf.gz
bcftools query -l apoe.vcf.gz | awk 'BEGIN{print "IID"}; {print \$1}' > sample.ids
bcftools query -f '%CHROM:%POS:%REF:%ALT[ %GT]\n' -r 19:45411941 apoe.vcf.gz | tr ' ' '\n' > snp1.txt
bcftools query -f '%CHROM:%POS:%REF:%ALT[ %GT]\n' -r 19:45412079 apoe.vcf.gz | tr ' ' '\n' > snp2.txt
paste sample.ids snp1.txt snp2.txt | sed 's/|/ /g' > apoe.hap.tmp
awk 'NR>1 {print \$1, \$2\$4, \$3\$5}' apoe.hap.tmp > apoe.hap.tmp2
sed -e 's/ 00/ e3/g' -e 's/ 10/ e4/g' -e 's/ 11/ e1/g' -e 's/ 01/ e2/g' apoe.hap.tmp2 > apoe.hap.tmp3
awk 'BEGIN{print "IID apoe"}{print \$1,\$2\$3}' apoe.hap.tmp3 | sed 's/e2e1/e1e2/;s/e3e1/e1e3/;s/e3e2/e2e3/;s/e4e1/e1e4/;s/e4e2/e2e4/;s/e4e3/e3e4/' > apoe.hap
awk 'NR>1{print \$2}' apoe.hap | sort | uniq -c > apoe.hap.count.txt
cat apoe.hap.count.txt
mv -f apoe.hap apoe.hap.count.txt "$dirout/"
mv -f apoe.log "$dirout/logs/"
EOF2
	chmod +x "$dirout/jobs/apoe.worker.sh"

	bash "$dirout/jobs/apoe.worker.sh"
}

pos_args=()
flank=
flank_set=N
while [[ $# -gt 0 ]]; do
	case "$1" in
		-h | --help)
			usage
			exit 0
			;;
		--dir-gen)
			dir_gen=${2:?ERROR: --dir-gen needs a value}
			shift 2
			;;
		--data)
			snp_data=${2:?ERROR: --data needs a value}
			shift 2
			;;
		--snp-file)
			snp_file=${2:?ERROR: --snp-file needs a value}
			shift 2
			;;
		--snp-file-type)
			snp_file_type=${2:?ERROR: --snp-file-type needs a value}
			shift 2
			;;
		--flank)
			flank=${2:?ERROR: --flank needs a value}
			flank_set=Y
			shift 2
			;;
		--chr)
			lift_chr=${2:?ERROR: --chr needs a value}
			lift_chr=${lift_chr#chr}
			[[ $lift_chr =~ ^([1-9]|1[0-9]|2[0-6]|X|X\.male|Y|XY|MT|M)$ ]] || { echo "ERROR: invalid --chr" >&2; exit 2; }
			shift 2
			;;
		--bim-file)
			bim_file=${2:?ERROR: --bim-file needs a value}
			shift 2
			;;
		--dir-in)
			dir_in=${2:?ERROR: --dir-in needs a value}
			shift 2
			;;
		--input | --pgs-input)
			pgs_input=${2:?ERROR: $1 needs a value}
			shift 2
			;;
		--pgs-pfile-dir)
			pgs_pfile_dir=${2:?ERROR: --pgs-pfile-dir needs a value}
			shift 2
			;;
		--pgs-threads)
			pgs_threads=${2:?ERROR: --pgs-threads needs a value}
			shift 2
			;;
		--dir-out)
			dir_out=${2:?ERROR: --dir-out needs a value}
			shift 2
			;;
		--chain-file)
			chain_file=${2:?ERROR: --chain-file needs a value}
			shift 2
			;;
		--ref-fasta)
			ref_fasta=${2:?ERROR: --ref-fasta needs a value}
			shift 2
			;;
		--label)
			label=${2:?ERROR: --label needs a value}
			shift 2
			;;
		--*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 1
			;;
		*)
			pos_args+=("$1")
			shift
			;;
	esac
done
if ((${#pos_args[@]} > 1)); then
	echo "ERROR: too many positional arguments: ${pos_args[*]}" >&2
	usage >&2
	exit 1
fi

task=${pos_args[0]:-all}
task=$(echo "$task" | tr '[:upper:]' '[:lower:]')

case "$task" in pgs|liftover|snp|apoe|all) ;; *) echo "ERROR: unknown task: $task" >&2; exit 2 ;; esac
if [[ -n ${label:-} && ! $label =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo 'ERROR: --label must be a path-free name' >&2; exit 2
fi
if [[ -n ${pgs_threads:-} ]]; then
    [[ $pgs_threads =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --pgs-threads must be positive' >&2; exit 2; }
    (( UKB_THREADS >= pgs_threads )) || UKB_THREADS=$pgs_threads
fi
# Validate required arguments before dispatch; data scans and tools run in the worker.
case "$task" in
    pgs|all) : "${pgs_input:?ERROR: pgs requires --input}" "${label:?ERROR: pgs requires --label}" ;;
esac
case "$task" in
    snp|all) : "${snp_data:?ERROR: snp requires --data}" "${snp_file:?ERROR: snp requires --snp-file}"         "${snp_file_type:?ERROR: snp requires --snp-file-type}" "${dir_out:?ERROR: snp requires --dir-out}" "${label:?ERROR: snp requires --label}" ;;
    liftover) : "${dir_out:?ERROR: liftOver requires --dir-out}" "${chain_file:?ERROR: liftOver requires --chain-file}"
        [[ -n ${bim_file:-} ]] || : "${dir_in:?ERROR: liftOver requires --dir-in}" ;;
esac
dir_gen=${dir_gen:-$DIR_GEN/ukb/37}
job_root=${dir_out:-$ANALYSIS_ROOT/ukb/gen/$task}
if [[ $task == pgs && -z ${dir_out:-} ]]; then job_root=$job_root/$label; fi
ukb_dispatch "gen.$task${label:+.$label}" "$job_root" "$SCRIPT_DIR/gen.sh" "${worker_args[@]}"
ukb_worker_runtime
if [[ $task == pgs || $task == all ]]; then
    PGS_PLINK_MEMORY_MB=${PGS_PLINK_MEMORY_MB:-8192}
    [[ $PGS_PLINK_MEMORY_MB =~ ^[1-9][0-9]*$ ]] && (( PGS_PLINK_MEMORY_MB <= UKB_MEMORY_GIB * 1024 )) || {
        echo 'ERROR: PGS_PLINK_MEMORY_MB exceeds the job RAM budget or is invalid' >&2; exit 2;
    }
    export PGS_PLINK_MEMORY_MB
fi

case "$task" in
	pgs) run_pgs ;;
	liftover) run_liftover ;;
	snp) run_snp ;;
	apoe) run_apoe ;;
	all)
		run_pgs
		run_snp
		run_apoe
		;;
	*)
		echo "ERROR: unknown task: $task" >&2
		usage >&2
		exit 1
		;;
esac

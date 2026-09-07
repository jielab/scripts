#!/bin/bash -l



set +H
set -e

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

usage() {
    cat <<'USAGE'
Usage:
  ./gen.sh [options]
  ./gen.sh [TASK]
  ./gen.sh pgs --input FILE --label LABEL [--pgs-pfile-dir DIR] [--dir-out DIR]
  ./gen.sh liftGen --dir-in DIR --dir-out DIR --chain-file FILE [--ref-fasta FILE]
  ./gen.sh snp --dir-gen DIR --data typ|imp --snp-file FILE \
      --snp-file-type bed|snp [--flank BP] --dir-out DIR --label LABEL

Options:
  -h, --h, --help        Show this help message.

liftGen options:
  --bim-file FILE        Lift one BIM instead of PVARs; preserve rows/alleles,
                         mark unsafe mappings POS=-1, link BED/FAM/NOSEX.
  --dir-in DIR           GRCh37 pfile directory containing chr*.pvar.zst.
  --dir-out DIR          GRCh38 output directory.
  --chain-file FILE      hg19-to-hg38 liftOver chain (.gz is accepted).
  --ref-fasta FILE       GRCh38 FASTA used to verify REF. Default:
                          /mnt/d/data.BIG/refGen/fasta/GRCH38.fasta

PGS options:
  --input FILE            One whitespace-delimited score file. Required headers:
                          SNP/ID, CHR, EA/A1/refA, and BETA/bJ.
  --pgs-pfile-dir DIR     Chromosome-wise UKB imputed pfiles.
                          Default: /mnt/h/ukbGen/38/imp
  --pgs-threads INT       PLINK 2 threads per chromosome. Default: 1.
  --dir-out DIR           Work/output directory. Default:
                          /mnt/g/analysis/pgs/LABEL
  --label LABEL           Output prefix; final data are LABEL.pgs.txt.
  --dry-run               Write LABEL.pgs.cmd without submitting/running it.

SNP options:
  --dir-gen DIR          Genetic data root containing typ/ and imp/.
  --data typ|imp         Genotyped PLINK bed files or imputed PLINK 2 pfiles.
  --snp-file FILE        BED file, one-ID-per-line list, or a tabular SNP file.
                          For a table (including common/snp.lst), the first
                          column is used and a SNP header row is skipped.
  --snp-file-type TYPE   bed (0-based, half-open) or snp.
  --flank BP             BED only: add BP on each side. Default: 10000.
  --dir-out DIR          Output directory.
  --label LABEL          Output prefix; final data are LABEL.raw.
  --dry-run              Write the extraction script without running it.

Tasks:
  pgs       Calculate one score input against chromosome-wise UKB pfiles.
  liftGen   Lift chr*.pvar.zst from GRCh37 to GRCh38, preserving rows.
  snp       Submit SNP extraction job.
  apoe      Submit APOE haplotype extraction job.
  all       Run all task blocks in order. Default.

Examples:
  ./gen.sh pgs --input /mnt/g/data/gwas/bmi.ref --label bmi --pgs-pfile-dir /mnt/h/ukbGen/38/imp --dir-out /mnt/g/analysis/pgs/bmi
  ./gen.sh liftGen --dir-in /mnt/h/ukbGen/37/imp --dir-out /mnt/h/ukbGen/38/imp --chain-file /mnt/d/files/liftOver/hg19ToHg38.over.chain.gz
  ./gen.sh snp --dir-gen /mnt/d/data/ukb/gen/ --data typ --snp-file /mnt/d/files/Wang_MY.bed --snp-file-type bed --dir-out /mnt/d/analysis/ukb/ --label Wang
USAGE
}

init_mnt_g_env() {
    dir0=${DIR0:-/mnt/g}
    local phe_f="$SCRIPT_DIR/../0f/0phe.f.sh"

    if [[ -f "$phe_f" ]]; then
        PHE_F=$phe_f
    elif [[ -f "$dir0/scripts/0f/0phe.f.sh" ]]; then
        PHE_F=$dir0/scripts/0f/0phe.f.sh
    else
        echo "ERROR: cannot find 0phe.f.sh beside gen.sh or under $dir0/scripts/0f" >&2
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
    local pvar=$1
    if [[ "$pvar" == *.zst ]]; then
        zstdcat -- "$pvar"
    else
        cat -- "$pvar"
    fi
}

liftgen_output_is_complete() {
    # The log is written last and acts as the completion marker.  Validate the
    # marker, dependencies, links, and compressed stream before resuming past it.
    local input_pvar=$1 input_pgen=$2 input_psam=$3 output_pvar=$4 output_log=$5
    local output_dir output_name output_stem
    output_dir=$(dirname "$output_pvar")
    output_name=$(basename "$output_pvar")
    output_stem=${output_name%.pvar.zst}

    [[ -s "$output_pvar" && -s "$output_log" ]] || return 1
    [[ "$output_pvar" -nt "$input_pvar" ]] || return 1
    [[ "$output_pvar" -nt "$chain_file" ]] || return 1
    [[ "$output_pvar" -nt "$ref_fasta" && "$output_pvar" -nt "$ref_fasta.fai" ]] || return 1
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

    zstd -q -t -- "$output_pvar" >/dev/null 2>&1
}

run_liftgen() {
    if [[ -n "${bim_file:-}" ]]; then
        python3 "$SCRIPT_DIR/f/lift_bim.py" "$(to_wsl_path "$bim_file")" \
            "$(to_wsl_path "${dir_out:?ERROR: --dir-out required}")" \
            "$(to_wsl_path "${chain_file:?ERROR: --chain-file required}")" \
            "$(to_wsl_path "${ref_fasta:-/mnt/d/data.BIG/refGen/fasta/GRCH38.fasta}")"
        return
    fi
    # Keep every PVAR row in its original position.  Only POS changes; rows that
    # cannot be mapped and verified against the target reference receive POS=-1.
    set -o pipefail
    : "${dir_in:?ERROR: liftGen requires --dir-in}"
    : "${dir_out:?ERROR: liftGen requires --dir-out}"
    : "${chain_file:?ERROR: liftGen requires --chain-file}"

    dir_in=$(to_wsl_path "$dir_in")
    dir_out=$(to_wsl_path "$dir_out")
    chain_file=$(to_wsl_path "$chain_file")
    ref_fasta=$(to_wsl_path "${ref_fasta:-/mnt/d/data.BIG/refGen/fasta/GRCH38.fasta}")

    [[ -d "$dir_in" ]] || { echo "ERROR: input directory not found: $dir_in" >&2; exit 1; }
    [[ -f "$chain_file" ]] || { echo "ERROR: chain file not found: $chain_file" >&2; exit 1; }
    [[ -f "$ref_fasta" ]] || { echo "ERROR: GRCh38 FASTA not found: $ref_fasta" >&2; exit 1; }
    [[ -f "$ref_fasta.fai" ]] || { echo "ERROR: FASTA index not found: $ref_fasta.fai" >&2; exit 1; }

    local lift_over=${LIFTOVER_BIN:-}
    if [[ -z "$lift_over" ]]; then
        if command -v liftOver >/dev/null 2>&1; then
            lift_over=$(command -v liftOver)
        elif [[ -x /mnt/d/software/bin/liftOver ]]; then
            lift_over=/mnt/d/software/bin/liftOver
        else
            echo "ERROR: liftOver not found (set LIFTOVER_BIN if needed)" >&2
            exit 1
        fi
    fi
    command -v zstd >/dev/null 2>&1 || { echo "ERROR: zstd not found" >&2; exit 1; }
    command -v zstdcat >/dev/null 2>&1 || { echo "ERROR: zstdcat not found" >&2; exit 1; }
    command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found" >&2; exit 1; }
    command -v realpath >/dev/null 2>&1 || { echo "ERROR: realpath not found" >&2; exit 1; }

    mkdir -p "$dir_out"
    dir_in=$(cd "$dir_in" && pwd -P)
    dir_out=$(cd "$dir_out" && pwd -P)
    chain_file=$(cd "$(dirname "$chain_file")" && pwd -P)/$(basename "$chain_file")
    ref_fasta=$(cd "$(dirname "$ref_fasta")" && pwd -P)/$(basename "$ref_fasta")
    [[ "$dir_in" != "$dir_out" ]] || { echo "ERROR: --dir-in and --dir-out must differ" >&2; exit 1; }

    local -a pvars
    mapfile -t pvars < <(find "$dir_in" -maxdepth 1 -type f -name 'chr*.pvar.zst' -printf '%p\n' | sort -V)
    (( ${#pvars[@]} > 0 )) || { echo "ERROR: no chr*.pvar.zst files found in $dir_in" >&2; exit 1; }

    local lift_work_dir=
    cleanup_liftgen_work() {
        if [[ -n "$lift_work_dir" && -d "$lift_work_dir" ]]; then
            find "$lift_work_dir" -maxdepth 1 -type f -delete
            rmdir "$lift_work_dir" 2>/dev/null || true
        fi
    }
    trap cleanup_liftgen_work EXIT
    trap 'exit 130' INT TERM

    local pvar name stem pgen psam out_pvar out_log
    local source_bed mapped_bed unmapped_bed sorted_bed candidate_bed verified_map out_tmp out_log_tmp
    local n_total n_candidate n_ok n_fail rel target
    for pvar in "${pvars[@]}"; do
        name=$(basename "$pvar")
        stem=${name%.pvar.zst}
        pgen="$dir_in/$stem.pgen"
        psam="$dir_in/$stem.psam"
        [[ -e "$pgen" ]] || { echo "ERROR: missing companion file: $pgen" >&2; exit 1; }
        [[ -e "$psam" ]] || { echo "ERROR: missing companion file: $psam" >&2; exit 1; }

        out_pvar="$dir_out/$name"
        out_log="$dir_out/$name.liftGen.log"
        if liftgen_output_is_complete "$pvar" "$pgen" "$psam" "$out_pvar" "$out_log"; then
            echo "SKIP $name (complete)"
            continue
        fi
        if [[ -e "$out_pvar" || -e "$out_log" ]]; then
            echo "REBUILD $name (existing result is incomplete or stale)"
        fi

        lift_work_dir="$dir_out/.liftGen.${stem}.$$"
        mkdir "$lift_work_dir"
        source_bed="$lift_work_dir/source.bed"
        mapped_bed="$lift_work_dir/mapped.bed"
        unmapped_bed="$lift_work_dir/unmapped.bed"
        sorted_bed="$lift_work_dir/mapped.sorted.bed"
        candidate_bed="$lift_work_dir/candidate.bed"
        verified_map="$lift_work_dir/verified.tsv"
        out_tmp="$lift_work_dir/$name"
        out_log_tmp="$lift_work_dir/$name.liftGen.log"

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
        (( n_total > 0 )) || { echo "ERROR: no variant rows in $pvar" >&2; exit 1; }
        : > "$mapped_bed"
        : > "$unmapped_bed"
        : > "$sorted_bed"
        if [[ -s "$source_bed" ]]; then
            "$lift_over" "$source_bed" "$chain_file" "$mapped_bed" "$unmapped_bed" >/dev/null
            LC_ALL=C sort -k4,4n -k2,2n "$mapped_bed" > "$sorted_bed"
        fi

        # Accept exactly one forward-strand mapping to the expected chromosome,
        # with the same interval length.  REF is verified in the next pass.
        : > "$candidate_bed"
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

        : > "$verified_map"
        if [[ -s "$candidate_bed" ]]; then
            bedtools getfasta -fi "$ref_fasta" -bed "$candidate_bed" -nameOnly -tab 2>/dev/null |
                awk -F '\t' 'BEGIN{OFS="\t"} {
                    n=split($1,a,"|")
                    if (n==3 && toupper($2)==toupper(a[3])) print a[1],a[2]
                }' > "$verified_map"
        fi
        n_candidate=$(wc -l < "$candidate_bed")
        n_ok=$(wc -l < "$verified_map")
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
        ' | zstd -T0 -q -f -o "$out_tmp"
        zstd -q -t "$out_tmp"
        [[ $(zstdcat -- "$out_tmp" | awk '!/^#/{n++} END{print n+0}') -eq "$n_total" ]] || {
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
        } > "$out_log_tmp"
        mv -f "$out_log_tmp" "$out_log"
        echo "OK $name: total=$n_total lifted=$n_ok POS=-1=$n_fail"

        cleanup_liftgen_work
        lift_work_dir=
    done
    trap - EXIT INT TERM
    echo "OK: lifted PVARs and relative pfile links written to: $dir_out"
}

run_pgs() {
    init_mnt_g_env
    : "${pgs_input:?ERROR: pgs requires --input FILE}"
    : "${label:?ERROR: pgs requires --label LABEL}"
    pgs_input=$(to_wsl_path "$pgs_input")
    pgs_pfile_dir=$(to_wsl_path "${pgs_pfile_dir:-/mnt/h/ukbGen/38/imp}")
    dir_out=$(to_wsl_path "${dir_out:-$dir0/analysis/pgs/$label}")
    pgs_threads=${pgs_threads:-1}
    [[ -s "$pgs_input" ]] || { echo "ERROR: PGS score input is missing/empty: $pgs_input" >&2; exit 1; }
    [[ "$pgs_threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --pgs-threads must be a positive integer" >&2; exit 1; }
    declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from $PHE_F" >&2; exit 1; }

    mkdir -p "$dir_out"
    local cmd_file="$dir_out/$label.pgs.cmd"
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
pgs_plink_calc \
    --input "$score_file" \
    --pfile-dir "$pfile_dir" \
    --output "$out_dir/$label.pgs.txt" \
    --label "$label" \
    --work-dir "$out_dir" \
    --status-file "$out_dir/$label.pgs.status.tsv" \
    --threads "$threads"
echo "OK: PGS written to: $out_dir/$label.pgs.txt"
EOF2
    } > "$cmd_file"
    chmod +x "$cmd_file"

    if [[ "${dry_run:-N}" == Y ]]; then
        echo "OK: dry run; PGS script written to: $cmd_file"
    elif command -v bsub >/dev/null 2>&1; then
        (cd "$dir_out" && bsub -q short -J "$label.pgs" -o "$label.pgs.LOG" -e "$label.pgs.ERR" < "$cmd_file")
        echo "OK: PGS job submitted: $cmd_file"
    else
        bash "$cmd_file"
    fi
}

run_snp() {


    # 🚩 Extract SNP data
    : "${dir_gen:?ERROR: snp requires --dir-gen}"
    : "${snp_data:?ERROR: snp requires --data typ|imp}"
    : "${snp_file:?ERROR: snp requires --snp-file}"
    : "${snp_file_type:?ERROR: snp requires --snp-file-type bed|snp}"
    : "${dir_out:?ERROR: snp requires --dir-out}"
    : "${label:?ERROR: snp requires --label}"
    [[ "$snp_data" == typ || "$snp_data" == imp ]] || { echo "ERROR: --data must be typ or imp" >&2; exit 1; }
    [[ "$snp_file_type" == bed || "$snp_file_type" == snp ]] || { echo "ERROR: --snp-file-type must be bed or snp" >&2; exit 1; }
    if [[ "$snp_file_type" == bed ]]; then
        flank=${flank:-10000}
        [[ "$flank" =~ ^[0-9]+$ ]] || { echo "ERROR: --flank must be a nonnegative integer" >&2; exit 1; }
    elif [[ "${flank_set:-N}" == Y ]]; then
        echo "ERROR: --flank is only valid with --snp-file-type bed" >&2
        exit 1
    fi
    [[ -d "$dir_gen/$snp_data" ]] || { echo "ERROR: missing data directory: $dir_gen/$snp_data" >&2; exit 1; }
    [[ -f "$snp_file" ]] || { echo "ERROR: missing SNP file: $snp_file" >&2; exit 1; }

    dir_gen=$(cd "$dir_gen" && pwd)
    snp_file=$(cd "$(dirname "$snp_file")" && pwd)/$(basename "$snp_file")
    mkdir -p "$dir_out"
    dir_out=$(cd "$dir_out" && pwd)
    work_dir="$dir_out/.${label}.snp-work"
    mkdir -p "$work_dir"

    if [[ "$snp_file_type" == bed ]]; then
        range_file="$work_dir/$label.flank${flank}.bed"
        awk -v flank="$flank" 'BEGIN{OFS="\t"} /^[[:space:]]*#/ || NF<3 {next} {s=$2-flank; if(s<0)s=0; print $1,s,$3+flank,(NF>=4?$4:".")}' "$snp_file" > "$range_file"
        [[ -s "$range_file" ]] || { echo "ERROR: no BED intervals found in $snp_file" >&2; exit 1; }
        extract_args="--extract bed0 \"$range_file\""
        mapfile -t chromosomes < <(awk '{c=$1; sub(/^chr/,"",c); print c}' "$range_file" | sort -uV)
    else
        extract_file="$work_dir/$label.extract.snp"
        awk 'BEGIN{FS="[[:space:]]+"}
             /^[[:space:]]*#/ || /^[[:space:]]*$/ {next}
             !seen && toupper($1) ~ /^(SNP|RSID|VARIANT|ID)$/ {seen=1; next}
             {seen=1; print $1}' "$snp_file" > "$extract_file"
        [[ -s "$extract_file" ]] || { echo "ERROR: no variant IDs found in $snp_file" >&2; exit 1; }
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
    (( ${#chromosomes[@]} > 0 )) || { echo "ERROR: no chromosome data files found" >&2; exit 1; }

    for chr in "${chromosomes[@]}"; do
        if [[ "$snp_data" == typ ]]; then
            [[ -f "$dir_gen/typ/chr$chr.bed" && -f "$dir_gen/typ/chr$chr.bim" && -f "$dir_gen/typ/chr$chr.fam" ]] || {
                echo "ERROR: incomplete typ files for chr$chr" >&2; exit 1;
            }
        else
            [[ -f "$dir_gen/imp/chr$chr.pgen" && ( -f "$dir_gen/imp/chr$chr.pvar" || -f "$dir_gen/imp/chr$chr.pvar.zst" ) && -f "$dir_gen/imp/chr$chr.psam" ]] || {
                echo "ERROR: incomplete imp files for chr$chr" >&2; exit 1;
            }
        fi
    done

    cmd_file="$work_dir/$label.cmd"
    {
        cat <<EOF2
#!/bin/bash -l
set +H
set -e

cd "$work_dir"
rm -f merge-list.txt
for chr in ${chromosomes[*]}; do
    if [[ "$snp_data" == "typ" ]]; then
        input_args="--bfile \"$dir_gen/typ/chr\$chr\""
    else
        pfile_modifier=""
        [[ ! -s "$dir_gen/imp/chr\$chr.pvar.zst" ]] || pfile_modifier="vzs"
        input_args="--pfile \"$dir_gen/imp/chr\$chr\" \$pfile_modifier"
    fi
    set +e
    eval plink2 \$input_args $extract_args --make-pgen --out "chr\$chr"
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
plink2 --pmerge-list merge-list.txt --make-pgen --out "$label"
plink2 --pfile "$label" --maj-ref --export A --out "$dir_out/$label"
EOF2
    } > "$cmd_file"
    chmod +x "$cmd_file"

    if [[ "${dry_run:-N}" == Y ]]; then
        echo "OK: dry run; extraction script written to: $cmd_file"
    elif command -v bsub >/dev/null 2>&1; then
        (cd "$work_dir" && bsub -q spec -J "$label" -o "$label.LOG" -e "$label.ERR" < "$cmd_file")
        echo "OK: SNP extraction job submitted: $cmd_file"
    else
        bash "$cmd_file"
        echo "OK: SNP data written to: $dir_out/$label.raw"
    fi
}

run_apoe() {


    # 🚩 Extract APOE haplotype
    # e1: C-T; e2: T-T; e3: T-C; e4: C-C
    dir0=/work/sph-huangj
    dirout=$dir0/data/ukb/gen/tmp
    mkdir -p "$dirout"

    echo -e "rs429358\nrs7412" > "$dirout/apoe.snps"

    cat > "$dirout/apoe.cmd" <<EOF2
#!/bin/bash -l
set +H
set -e

cd "$dirout"
plink2 --pfile $dir0/data/ukb/imp/hap/chr19 --extract apoe.snps --export vcf id-paste=iid bgz --out apoe
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
EOF2
    chmod +x "$dirout/apoe.cmd"

    (cd "$dirout" && bsub -q spec -J apoe -o apoe.LOG -e apoe.ERR < apoe.cmd)
    echo "OK: APOE extraction job submitted: $dirout/apoe.cmd"
}

pos_args=()
flank=
flank_set=N
dry_run=N
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --dir-gen) dir_gen=${2:?ERROR: --dir-gen needs a value}; shift 2 ;;
        --data) snp_data=${2:?ERROR: --data needs a value}; shift 2 ;;
        --snp-file) snp_file=${2:?ERROR: --snp-file needs a value}; shift 2 ;;
        --snp-file-type) snp_file_type=${2:?ERROR: --snp-file-type needs a value}; shift 2 ;;
        --flank) flank=${2:?ERROR: --flank needs a value}; flank_set=Y; shift 2 ;;
        --bim-file) bim_file=${2:?ERROR: --bim-file needs a value}; shift 2 ;;
        --dir-in) dir_in=${2:?ERROR: --dir-in needs a value}; shift 2 ;;
        --input|--pgs-input) pgs_input=${2:?ERROR: $1 needs a value}; shift 2 ;;
        --pgs-pfile-dir) pgs_pfile_dir=${2:?ERROR: --pgs-pfile-dir needs a value}; shift 2 ;;
        --pgs-threads) pgs_threads=${2:?ERROR: --pgs-threads needs a value}; shift 2 ;;
        --dir-out) dir_out=${2:?ERROR: --dir-out needs a value}; shift 2 ;;
        --chain-file) chain_file=${2:?ERROR: --chain-file needs a value}; shift 2 ;;
        --ref-fasta) ref_fasta=${2:?ERROR: --ref-fasta needs a value}; shift 2 ;;
        --label) label=${2:?ERROR: --label needs a value}; shift 2 ;;
        --dry-run) dry_run=Y; shift ;;
        --*) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
        *) pos_args+=("$1"); shift ;;
    esac
done
if (( ${#pos_args[@]} > 1 )); then
    echo "ERROR: too many positional arguments: ${pos_args[*]}" >&2
    usage >&2
    exit 1
fi

task=${pos_args[0]:-all}
task=$(echo "$task" | tr '[:upper:]' '[:lower:]')

case "$task" in
    pgs)    run_pgs ;;
    liftgen) run_liftgen ;;
    snp)    run_snp ;;
    apoe)   run_apoe ;;
    all)    run_pgs; run_snp; run_apoe ;;
    *)
        echo "ERROR: unknown task: $task" >&2
        usage >&2
        exit 1
        ;;
esac

#!/usr/bin/env bash
set -euo pipefail

dir0=${dir0:-/mnt/d}
dirfunc=${dirfunc:-$dir0/scripts/0f}


# 🚩 Reference-genome builder
usage() {
    cat <<'EOF'
Usage:
  ./refGen.sh make-pfile [options]
  ./refGen.sh download_archaic [options]
  ./refGen.sh lift-bed --dir-in DIR --dir-out DIR --chain-file FILE [options]

Modules:
  make-pfile
    Convert chr*.vcf.gz, chr*.bed/bim/fam, or UKB-style per-chromosome
    BGEN/SAMPLE files to PLINK 2 pfiles and/or convert existing pfiles to
    PLINK 1 bed/bim/fam files.

    --out-files pgen       Run --make-pgen vzs, retaining VCF INFO and writing
                           .pvar.zst (default).
    --out-files bed        Run --max-alleles 2 --make-bed on existing pfiles.
    --out-files pgen,bed   Generate pfiles first, then bed/bim/fam.

    The bed conversion excludes multiallelic variants because .bed/.bim cannot
    represent them, and uses classic numeric chromosome codes (X=23, PAR=25,
    MT=26). chrX always means non-PAR X, chrXY contains X PAR1/PAR2 (stored
    as XY for PLINK 1 input or PLINK2's PAR1/PAR2 codes for VCF input), and chrY always
    means non-PAR Y. When --chr is omitted, pgen discovers and processes all
    supported chromosomes in the input directory; bed-only discovers all
    complete chr*.pgen/.pvar[.zst]/.psam filesets.

  download_archaic
    Download Vindija, Altai, Chagyrskaya, Denisova, and Denisova25 VCFs and
    indexes for chromosomes 1-22 and X. Repack non-BGZF Chagyrskaya files,
    then rebuild all downloaded tabix indexes.

  lift-bed
    Lift every *.bed.gz below --dir-in with the supplied UCSC chain file.
    Write BGZF-compressed files below --dir-out with the same relative paths
    and filenames. Unmapped records are retained below --dir-out/.unmapped/.
    The chain file itself determines the source and destination builds.

Common options:
  --sample-file FILE       VCF sample metadata [dirname(input-dir)/samples.txt],
                           columns/order:
                           sample fatherID motherID sex pop super_pop
  --input-prefix PREFIX    Input basename prefix [chr]. Nonstandard inputs
                           must specify this (for example, cohort_chr).
  --vcf-prefix PREFIX      Backward-compatible alias for --input-prefix.
  --replace TRUE|FALSE     Replace complete requested outputs [FALSE].
  -h, --help               Show this help.

make-pfile options:
  --out-files LIST         pgen, bed, or pgen,bed [pgen].
  --dir-input DIR          Directory containing input files. VCF expects
                           chr*.vcf.gz; BED expects chr*.bed/bim/fam; BGEN
                           accepts <prefix><chr>.bgen/sample or a unique
                           <prefix><chr>_*.bgen plus matching *_s*.sample.
  --input-format FORMAT    auto, vcf, bgen, or bed [auto].
  --dir-vcf DIR            Backward-compatible VCF-only alias for --dir-input
                           [/mnt/d/data.BIG/refGen/1kg/<grch>/vcf].
  --dir-pfile DIR          Directory for PLINK outputs
                           [/mnt/d/data.BIG/refGen/1kg/<grch>/pfile].
                           <grch> defaults to 37 when both dirs are omitted.
  --grch 37|38             Optional chrX PAR definition. If omitted, infer the
                           build from the chrX VCF header.
  --chr LIST               Chromosomes such as 2,X or chr2,chrX. Omit for all.
                           --chr2,X and --chr=2,X are accepted aliases.
  --split-sex TRUE|FALSE   With chrX, also make chrX.male and chrX.female
                           pfiles. Male validation reports include .afreq,
                           .gcount, .smiss, and .basic-stats.txt. When bed is
                           requested, also make both sex-specific BED filesets
                           [FALSE].

download_archaic options:
  --root DIR               Archaic data root
                           [/mnt/d/data.BIG/refGen/archaic].
  --chr LIST               Chromosomes such as 2,X or autosomes [1-22,X].
  --jobs N                 Parallel tabix indexing jobs [8].
  --bgzip-threads N        Threads used to repack Chagyrskaya VCFs [4].

lift-bed options:
  --dir-in DIR             Input tree containing *.bed.gz files [required].
  --dir-out DIR            Output tree; relative paths are preserved [required].
  --chain-file FILE        UCSC liftOver chain, optionally gzip-compressed
                           [required]. Its direction determines the conversion.
  --bgzip-threads N        Threads used to compress output BED files [4].

Examples:
  # 1KG GRCh37: generate pgen/pvar/psam and bed/bim/fam for all chromosomes.
  ./refGen.sh make-pfile --out-files pgen,bed \
    --dir-vcf /mnt/d/data.BIG/refGen/1kg/37/vcf \
    --dir-pfile /mnt/d/data.BIG/refGen/1kg/37/pfile \
    --grch 37

  # 1KG GRCh38: generate pgen/pvar/psam and bed/bim/fam for all chromosomes.
  ./refGen.sh make-pfile --out-files pgen,bed \
    --dir-vcf /mnt/d/data.BIG/refGen/1kg/38/vcf \
    --dir-pfile /mnt/d/data.BIG/refGen/1kg/38/pfile \
    --grch 38

  # Convert existing GRCh37 chr2 and chrX pfiles to bed/bim/fam only.
  ./refGen.sh make-pfile --out-files bed --chr 2,X \
    --dir-vcf /mnt/d/data.BIG/refGen/1kg/37/vcf \
    --dir-pfile /mnt/d/data.BIG/refGen/1kg/37/pfile

  # 1KG GRCh37: pure non-PAR chrX, chrXY PAR1/PAR2, and sex-split chrX.
  ./refGen.sh make-pfile --out-files pgen,bed \
    --chr X --split-sex TRUE --grch 37 \
    --dir-vcf /mnt/d/data.BIG/refGen/1kg/37/vcf \
    --dir-pfile /mnt/d/data.BIG/refGen/1kg/37/pfile

  # 1KG GRCh38: the same outputs, using GRCh38 PAR definitions.
  ./refGen.sh make-pfile --out-files pgen,bed \
    --chr X --split-sex TRUE --grch 38 \
    --dir-vcf /mnt/d/data.BIG/refGen/1kg/38/vcf \
    --dir-pfile /mnt/d/data.BIG/refGen/1kg/38/pfile

  # UK Biobank phased BGEN/SAMPLE files, named raw/chr*.bgen/sample.
  ./refGen.sh make-pfile --input-format bgen \
    --dir-input /mnt/h/ukbGen/37/hap/raw --dir-pfile /mnt/h/ukbGen/37/hap

  # UK Biobank array BED/BIM/FAM files, named raw/chr*.*.
  ./refGen.sh make-pfile --input-format bed --split-sex TRUE \
    --dir-input /mnt/h/ukbGen/37/typ/raw --dir-pfile /mnt/h/ukbGen/37/typ

  ./refGen.sh download_archaic \
    --root /mnt/d/data.BIG/refGen/archaic --jobs 8

  # Archaic BED masks: GRCh37 to GRCh38 (direction comes from the chain file).
  ./refGen.sh lift-bed \
    --dir-in /mnt/d/data.BIG/refGen/archaic/37/bed \
    --dir-out /mnt/d/data.BIG/refGen/archaic/38/bed \
    --chain-file /mnt/d/files/liftOver/hg19ToHg38.over.chain.gz

Reference download locations (download manually into the corresponding vcf/):
  1KG GRCh37 Phase 3:
    https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
  1KG GRCh38 high-coverage phased (3202 samples):
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
  1KG GRCh38 phased SNV/INDEL/SV:
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
  1KG 3202-sample pedigree/sex metadata:
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
  Vindija, Altai, Denisova archaic VCFs:
    https://cdna.eva.mpg.de/neandertal/Vindija/VCF/
  Chagyrskaya archaic VCFs:
    https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/
  Denisova25 archaic VCFs:
    https://cdna.eva.mpg.de/denisova/Den25/VCF/
EOF
}


# 🚩 Shared helpers
die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_true() {
    case "${1^^}" in
        TRUE|T|YES|Y|1) return 0 ;;
        FALSE|F|NO|N|0) return 1 ;;
        *) die "expected TRUE or FALSE, got: $1" ;;
    esac
}

split_csv() {
    local input=$1
    local -n output_ref=$2
    IFS=',' read -r -a output_ref <<< "$input"
    ((${#output_ref[@]} > 0)) || die "empty comma-separated option"
}

pfile_exists() {
    local prefix=$1
    [[ -s ${prefix}.pgen && -s ${prefix}.psam && ( -s ${prefix}.pvar || -s ${prefix}.pvar.zst ) ]]
}

pfile_zst_exists() {
    local prefix=$1
    [[ -s ${prefix}.pgen && -s ${prefix}.psam && -s ${prefix}.pvar.zst ]]
}

pfile_input_args() {
    local prefix=$1
    # shellcheck disable=SC2178  # The second argument names an array.
    local -n output_ref=$2
    if [[ -s ${prefix}.pvar.zst ]]; then
        output_ref=(--pfile "$prefix" vzs)
    elif [[ -s ${prefix}.pvar ]]; then
        output_ref=(--pfile "$prefix")
    else
        die "missing pvar for: $prefix"
    fi
}

bed_exists() {
    local prefix=$1
    [[ -s ${prefix}.bed && -s ${prefix}.bim && -s ${prefix}.fam ]]
}

run_make_pfile() {
    local output=$1
    shift
    if pfile_zst_exists "$output" && ! is_true "$replace"; then
        echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
        return
    fi
    "$plink2_bin" "$@" --make-pgen vzs --out "$output"
    pfile_zst_exists "$output" || die "incomplete compressed pfile output: $output"
    # A same-prefix rebuild does not remove an obsolete uncompressed .pvar.
    # Remove it only after the new compressed fileset has been verified.
    rm -f -- "${output}.pvar"
}

run_make_pfile_replacing() {
    # Dynamic scoping makes run_make_pfile see this local override. This is
    # used only to repair an existing pfile with obsolete chromosome content.
    local replace=TRUE
    run_make_pfile "$@"
}

write_sex_update() {
    local output=$1 tmp=${1}.tmp.$$
    awk 'BEGIN{OFS="\t"}
        NR==1 {next}
        {
            s=tolower($4)
            if (s=="male" || s=="m") s=1
            else if (s=="female" || s=="f") s=2
            if (s!=1 && s!=2) {
                print "ERROR: invalid sex for sample " $1 ": " $4 > "/dev/stderr"
                bad=1; next
            }
            print $1,$1,s
        }
        END {exit bad}
    ' "$sample_file" > "$tmp"
    mv -f -- "$tmp" "$output"
}

normalize_chromosomes() {
    local input=$1 raw chr
    local -A seen=()
    local values=()
    split_csv "$input" values
    selected_chrs=()
    for raw in "${values[@]}"; do
        chr=${raw#chr}; chr=${chr#CHR}; chr=${chr^^}
        if [[ $chr == AUTOSOME || $chr == AUTOSOMES ]]; then
            for chr in {1..22}; do seen[$chr]=1; done
            continue
        fi
        [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == XY || $chr == Y || $chr == MT ]] ||
            die "invalid chromosome: $raw"
        seen[$chr]=1
    done
    mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

discover_vcf_chromosomes() {
    local file base chr
    local -A seen=()
    shopt -s nullglob
    for file in "${input_prefix}"*.vcf.gz; do
        base=${file#"$input_prefix"}
        chr=${base%.vcf.gz}; chr=${chr^^}
        [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == Y || $chr == MT ]] && seen[$chr]=1
    done
    shopt -u nullglob
    ((${#seen[@]} > 0)) || die "no chromosome VCFs found with prefix: $input_prefix"
    mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

discover_bed_chromosomes() {
    local file chr prefix
    local -A seen=()
    shopt -s nullglob
    for file in "${input_prefix}"*.bed; do
        chr=${file#"$input_prefix"}; chr=${chr%.bed}; chr=${chr^^}
        [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == XY || $chr == Y || $chr == MT ]] || continue
        prefix=${file%.bed}
        [[ -s ${prefix}.bim && -s ${prefix}.fam ]] ||
            die "incomplete BED/BIM/FAM input fileset: $prefix"
        seen[$chr]=1
    done
    shopt -u nullglob
    ((${#seen[@]} > 0)) || die "no complete BED/BIM/FAM inputs found with prefix: $input_prefix"
    mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

bgen_file_for_chr() {
    local chr=$1 exact=${input_prefix}${1}.bgen
    local candidates=()
    if [[ -s $exact ]]; then
        printf '%s\n' "$exact"
        return
    fi
    shopt -s nullglob
    candidates=("${input_prefix}${chr}"_*.bgen)
    shopt -u nullglob
    ((${#candidates[@]} == 1)) ||
        die "expected exactly one BGEN for chr$chr under $input_dir; found ${#candidates[@]}"
    printf '%s\n' "${candidates[0]}"
}

bgen_sample_file() {
    local bgen=$1 same_prefix=${1%.bgen}.sample
    local candidates=()
    if [[ -s $same_prefix ]]; then
        printf '%s\n' "$same_prefix"
        return
    fi
    shopt -s nullglob
    candidates=("${bgen%.bgen}"_s*.sample)
    shopt -u nullglob
    ((${#candidates[@]} == 1)) ||
        die "expected exactly one SAMPLE file for $bgen; found ${#candidates[@]}"
    printf '%s\n' "${candidates[0]}"
}

discover_bgen_chromosomes() {
    local chr bgen sample exact
    local candidates=()
    local -A seen=()
    for chr in {1..22} X XY Y MT; do
        exact=${input_prefix}${chr}.bgen
        if [[ -s $exact ]]; then
            seen[$chr]=1
            continue
        fi
        shopt -s nullglob
        candidates=("${input_prefix}${chr}"_*.bgen)
        shopt -u nullglob
        ((${#candidates[@]} <= 1)) ||
            die "multiple BGEN inputs match chr$chr and prefix $input_prefix"
        ((${#candidates[@]} == 1)) && seen[$chr]=1
    done
    ((${#seen[@]} > 0)) || die "no BGEN inputs found with prefix: $input_prefix"
    mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
    for chr in "${selected_chrs[@]}"; do
        bgen=$(bgen_file_for_chr "$chr")
        sample=$(bgen_sample_file "$bgen")
        [[ -s $sample ]] || die "missing SAMPLE file for chr$chr"
    done
}

detect_input_format() {
    local files=()
    case $input_format in
        vcf|bgen|bed) return ;;
        auto) ;;
        *) die "--input-format must be auto, vcf, bgen, or bed" ;;
    esac
    shopt -s nullglob
    files=("${input_prefix}"*.vcf.gz)
    if ((${#files[@]})); then input_format=vcf; shopt -u nullglob; return; fi
    files=("${input_prefix}"*.bgen)
    if ((${#files[@]})); then input_format=bgen; shopt -u nullglob; return; fi
    files=("${input_prefix}"*.bed)
    if ((${#files[@]})); then input_format=bed; shopt -u nullglob; return; fi
    shopt -u nullglob
    die "cannot detect VCF, BGEN, or BED input under $input_dir"
}

discover_input_chromosomes() {
    case $input_format in
        vcf) discover_vcf_chromosomes ;;
        bgen) discover_bgen_chromosomes ;;
        bed) discover_bed_chromosomes ;;
    esac
}

discover_pfile_chromosomes() {
    local file base chr
    local -A seen=()
    shopt -s nullglob
    for file in "$pfile_dir"/chr*.pgen; do
        base=${file##*/}; chr=${base#chr}; chr=${chr%.pgen}; chr=${chr^^}
        [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X || $chr == Y || $chr == MT ]] || continue
        pfile_exists "$pfile_dir/chr$chr" && seen[$chr]=1
    done
    shopt -u nullglob
    ((${#seen[@]} > 0)) || die "no complete chr*.pgen/.pvar/.psam filesets under $pfile_dir"
    mapfile -t selected_chrs < <(printf '%s\n' "${!seen[@]}" | sort -V)
}

detect_genome_build() {
    local vcf=$1 line detected=""
    while IFS= read -r line; do
        case $line in
            *assembly=b37*|*GRCh37*|*hs37d5*|*hg19*|*length=155270560*) detected=37; break ;;
            *assembly=b38*|*GRCh38*|*hg38*|*length=156040895*) detected=38; break ;;
            \#CHROM*) break ;;
        esac
    done < <(phe_zcat "$vcf")
    [[ $detected == 37 || $detected == 38 ]] ||
        die "cannot infer GRCh build from VCF header: $vcf (use --grch 37 or 38)"
    printf '%s\n' "$detected"
}

set_par_coordinates() {
    local build=$1
    # X PAR variables are supplied by the sourced 0phe.f.sh library.
    # shellcheck disable=SC2154
    if [[ $build == 37 ]]; then
        x_par_head_end=$X_PAR1_END_b37
        x_par_tail_start=$X_PAR2_START_b37
        y_par1_start=10001
        y_par1_end=2649520
        y_par2_start=59034050
        y_par2_end=59363566
        y_length=59373566
    else
        x_par_head_end=$X_PAR1_END_b38
        x_par_tail_start=$X_PAR2_START_b38
        y_par1_start=10001
        y_par1_end=2781479
        y_par2_start=56887903
        y_par2_end=57217415
        y_length=57227415
    fi
}

pvar_file() {
    local prefix=$1
    if [[ -s ${prefix}.pvar.zst ]]; then
        printf '%s\n' "${prefix}.pvar.zst"
    elif [[ -s ${prefix}.pvar ]]; then
        printf '%s\n' "${prefix}.pvar"
    else
        die "missing pvar for: $prefix"
    fi
}

read_pvar() {
    local pvar=$1 prefix tmp_dir
    local pfile_args=()
    case $pvar in
        *.pvar.zst)
            if command -v zstdcat >/dev/null 2>&1; then
                zstdcat -- "$pvar"
            elif command -v zstd >/dev/null 2>&1; then
                zstd -q -dc -- "$pvar"
            else
                # PLINK2 can always read its own compressed variant file.
                # Materialize a temporary uncompressed pvar for validation
                # when the standalone zstd tools are not installed.
                prefix=${pvar%.pvar.zst}
                pfile_input_args "$prefix" pfile_args
                tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/refGen.pvar.XXXXXX")
                (
                    trap 'rm -rf -- "$tmp_dir"' EXIT
                    "$plink2_bin" "${pfile_args[@]}" --make-just-pvar cols= \
                        --out "$tmp_dir/pvar" >/dev/null
                    cat -- "$tmp_dir/pvar.pvar"
                )
            fi
            ;;
        *) cat -- "$pvar" ;;
    esac
}

validate_pvar_chromosomes() {
    local prefix=$1 expected=$2 pvar
    pvar=$(pvar_file "$prefix")
    if ! pvar_chromosomes_match "$pvar" "$expected"; then
        die "$prefix failed the post-generation pure-$expected validation"
    fi
}

pvar_chromosomes_match() {
    local pvar=$1 expected=$2
    read_pvar "$pvar" | awk -v expected="$expected" '
        /^#/ {next}
        {
            count++
            if (expected == "X" && $1 != "X") bad=1
            if (expected == "Y" && $1 != "Y") bad=1
            if (expected == "PAR" && $1 != "XY" && $1 != "PAR1" && $1 != "PAR2") bad=1
        }
        END {if (!count) exit 2; exit bad}
    '
}

pvar_contains_par() {
    local pvar=$1
    read_pvar "$pvar" | awk '
        /^#/ {next}
        $1 == "PAR1" || $1 == "PAR2" {found=1}
        END {exit !found}
    '
}

validate_y_nonpar() {
    local prefix=$1 pvar
    pvar=$(pvar_file "$prefix")
    if ! read_pvar "$pvar" | awk \
        -v p1s="$y_par1_start" -v p1e="$y_par1_end" \
        -v p2s="$y_par2_start" -v p2e="$y_par2_end" '
        /^#/ {next}
        {
            count++
            if ($1 != "Y" || ($2 >= p1s && $2 <= p1e) ||
                ($2 >= p2s && $2 <= p2e)) bad=1
        }
        END {if (!count) exit 2; exit bad}
    '; then
        die "$prefix contains Y-PAR variants; use --replace TRUE after checking old outputs"
    fi
}

make_chrX_male_stats() {
    local prefix=$pfile_dir/chrX.male pvar summary_tmp
    local variant_count min_pos max_pos sample_count nonmale_count hethap_count
    local stats_complete=0 stats_build=unknown
    local pfile_args=()
    if [[ -n ${chrX_build:-} ]]; then
        stats_build=GRCh$chrX_build
    elif [[ -n $grch ]]; then
        stats_build=GRCh$grch
    fi

    if [[ -s ${prefix}.afreq && -s ${prefix}.gcount && -s ${prefix}.smiss &&
          -s ${prefix}.basic-stats.txt ]]; then
        stats_complete=1
    fi
    if ((stats_complete)) && ! is_true "$replace"; then
        echo "SKIP chrX.male statistics: reports already exist (use --replace TRUE)"
        return
    fi

    echo "STATS chrX.male (--freq, --geno-counts, --missing sample-only)"
    pfile_input_args "$prefix" pfile_args
    "$plink2_bin" "${pfile_args[@]}" --freq --out "$prefix"
    "$plink2_bin" "${pfile_args[@]}" --geno-counts --out "$prefix"
    "$plink2_bin" "${pfile_args[@]}" --missing sample-only \
        'scols=maybefid,maybesid,nmiss,nmisshh,hethap,nobs,fmiss,fmisshh' \
        --out "$prefix"
    [[ -s ${prefix}.afreq && -s ${prefix}.gcount && -s ${prefix}.smiss ]] ||
        die "incomplete chrX.male statistics under: $pfile_dir"

    pvar=$(pvar_file "$prefix")
    read -r variant_count min_pos max_pos < <(
        read_pvar "$pvar" | awk '
            /^#/ {next}
            NR == 1 || !count {min=$2}
            {count++; max=$2}
            END {print count+0, min+0, max+0}
        '
    )
    read -r sample_count nonmale_count < <(
        awk '
            /^#/ {
                for (i=1; i<=NF; i++) if ($i == "SEX") sex_col=i
                next
            }
            NF {count++; if (!sex_col || $sex_col != 1) nonmale++}
            END {print count+0, nonmale+0}
        ' "${prefix}.psam"
    )
    hethap_count=$(awk '
        NR == 1 {
            for (i=1; i<=NF; i++) if ($i == "HETHAP_CT") col=i
            next
        }
        col {sum += $col}
        END {if (!col) print "NA"; else print sum+0}
    ' "${prefix}.smiss")
    [[ $variant_count -gt 0 && $sample_count -gt 0 && $nonmale_count -eq 0 ]] ||
        die "chrX.male validation failed: variants=$variant_count samples=$sample_count nonmale=$nonmale_count"
    [[ $hethap_count == 0 ]] ||
        die "chrX.male validation failed: heterozygous-haploid count=$hethap_count"

    summary_tmp=${prefix}.basic-stats.txt.tmp.$$
    printf '%s\n' \
        "pfile=${prefix}" \
        "build=${stats_build}" \
        "chromosome=X" \
        "par_regions_present=FALSE" \
        "male_samples=${sample_count}" \
        "nonmale_samples=${nonmale_count}" \
        "variants=${variant_count}" \
        "min_position=${min_pos}" \
        "max_position=${max_pos}" \
        "plink2_encoding=male_chrX_haploid" \
        "heterozygous_haploid_calls=${hethap_count}" \
        "frequency_report=${prefix}.afreq" \
        "genotype_count_report=${prefix}.gcount" \
        "sample_missingness_report=${prefix}.smiss" \
        > "$summary_tmp"
    mv -f -- "$summary_tmp" "${prefix}.basic-stats.txt"
}

validate_psam_sex() {
    local prefix=$1 expected_code=$2 expected_label=$3
    if ! awk -v expected="$expected_code" '
        /^#/ {
            for (i=1; i<=NF; i++) if ($i == "SEX") sex_col=i
            next
        }
        NF {count++; if (!sex_col || $sex_col != expected) bad++}
        END {if (!count || bad) exit 1}
    ' "${prefix}.psam"; then
        die "$prefix failed the $expected_label-only sample validation"
    fi
}

make_chrX_sex_pfiles() {
    local input=$pfile_dir/chrX output sex keep_flag sex_code
    local pfile_args=()
    validate_pvar_chromosomes "$input" X
    pfile_input_args "$input" pfile_args
    for sex in male female; do
        output=$pfile_dir/chrX.$sex
        if [[ $sex == male ]]; then
            keep_flag=--keep-males; sex_code=1
        else
            keep_flag=--keep-females; sex_code=2
        fi
        if pfile_zst_exists "$output" && ! is_true "$replace"; then
            echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
        else
            echo "MAKE pgen chrX.$sex ($sex-only pure non-PAR X)"
            if [[ $sex == male ]]; then
                run_make_pfile "$output" "${pfile_args[@]}" "$keep_flag" \
                    --set-invalid-haploid-missing
            else
                run_make_pfile "$output" "${pfile_args[@]}" "$keep_flag"
            fi
        fi
        validate_pvar_chromosomes "$output" X
        validate_psam_sex "$output" "$sex_code" "$sex"
    done
    make_chrX_male_stats
}

make_chrX_pfiles() {
    local vcf=${input_prefix}X.vcf.gz x_output=$pfile_dir/chrX
    local xy_output=$pfile_dir/chrXY no_par_marker=$pfile_dir/chrXY.no-PAR.txt
    local build_source x_pvar xy_pvar need_import=0 x_valid=0 xy_valid=0
    [[ -s $vcf ]] || die "missing VCF: $vcf"

    if [[ -n $grch ]]; then
        chrX_build=$grch
        build_source="--grch"
    else
        chrX_build=$(detect_genome_build "$vcf")
        build_source="VCF header"
    fi
    set_par_coordinates "$chrX_build"

    if pfile_zst_exists "$x_output"; then
        x_pvar=$(pvar_file "$x_output")
        if pvar_chromosomes_match "$x_pvar" X; then
            x_valid=1
        else
            echo "REBUILD pgen: existing chrX contains PAR or non-X variants"
        fi
    fi
    if pfile_zst_exists "$xy_output"; then
        xy_pvar=$(pvar_file "$xy_output")
        if pvar_chromosomes_match "$xy_pvar" PAR; then
            xy_valid=1
        else
            echo "REBUILD pgen: existing chrXY contains non-PAR variants"
        fi
    fi

    ((x_valid)) || need_import=1
    if ((!xy_valid)) && [[ ! -s $no_par_marker ]]; then need_import=1; fi
    is_true "$replace" && need_import=1
    if ((need_import)); then
        (
            local tmp_dir raw raw_pvar raw_has_par=0
            local raw_args=()
            tmp_dir=$(mktemp -d "$pfile_dir/.chrX.make-pfile.XXXXXX")
            trap 'rm -rf -- "$tmp_dir"' EXIT
            raw=$tmp_dir/chrX.split
            echo "Using GRCh$chrX_build for chrX PAR splitting ($build_source): $vcf"
            run_make_pfile "$raw" \
                --vcf "$vcf" --double-id --update-sex "$sex_update" \
                --split-par "$x_par_head_end" "$x_par_tail_start" \
                --set-missing-var-ids '@:#:$r:$a' \
                --new-id-max-allele-len 100 missing
            raw_pvar=$(pvar_file "$raw")
            if pvar_contains_par "$raw_pvar"; then raw_has_par=1; fi
            pfile_input_args "$raw" raw_args

            if ((!x_valid)) || is_true "$replace"; then
                echo "MAKE pgen chrX (non-PAR only)"
                if pfile_exists "$x_output"; then
                    run_make_pfile_replacing "$x_output" "${raw_args[@]}" --chr X
                else
                    run_make_pfile "$x_output" "${raw_args[@]}" --chr X
                fi
            else
                echo "SKIP pgen: existing chrX is already pure non-PAR X"
            fi

            if ((raw_has_par)); then
                rm -f -- "$no_par_marker"
                if ((!xy_valid)) || is_true "$replace"; then
                    echo "MAKE pgen chrXY (PAR1 and PAR2)"
                    if pfile_exists "$xy_output"; then
                        run_make_pfile_replacing "$xy_output" "${raw_args[@]}" --chr PAR1 PAR2
                    else
                        run_make_pfile "$xy_output" "${raw_args[@]}" --chr PAR1 PAR2
                    fi
                else
                    echo "SKIP pgen: existing chrXY is already PAR1/PAR2 only"
                fi
            else
                printf 'source=%s\nbuild=GRCh%s\npar_variants=0\n' \
                    "$vcf" "$chrX_build" > "$no_par_marker"
                echo "SKIP pgen chrXY: source chrX VCF has no PAR1/PAR2 variants"
            fi
        )
    else
        echo "SKIP pgen: existing chrX/chrXY outputs are already valid"
    fi
    validate_pvar_chromosomes "$x_output" X
    if pfile_exists "$xy_output"; then
        [[ ! -s $no_par_marker ]] ||
            die "chrXY exists even though the source VCF has no PAR variants"
        validate_pvar_chromosomes "$xy_output" PAR
    else
        [[ -s $no_par_marker ]] || die "chrXY was not generated and no no-PAR marker was written"
    fi
    if is_true "$split_sex"; then make_chrX_sex_pfiles; fi
}

make_chrY_pfile() {
    local vcf=${input_prefix}Y.vcf.gz output=$pfile_dir/chrY
    local chrY_build build_source
    [[ -s $vcf ]] || die "missing VCF: $vcf"
    if [[ -n $grch ]]; then
        chrY_build=$grch
        build_source="--grch"
    else
        chrY_build=$(detect_genome_build "$vcf")
        build_source="VCF header"
    fi
    set_par_coordinates "$chrY_build"
    if pfile_zst_exists "$output" && ! is_true "$replace"; then
        echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
    else
        (
            local tmp_dir nonpar_bed
            tmp_dir=$(mktemp -d "$pfile_dir/.chrY.make-pfile.XXXXXX")
            trap 'rm -rf -- "$tmp_dir"' EXIT
            nonpar_bed=$tmp_dir/chrY.nonpar.bed1
            printf 'Y\t1\t%d\nY\t%d\t%d\nY\t%d\t%d\n' \
                "$((y_par1_start - 1))" \
                "$((y_par1_end + 1))" "$((y_par2_start - 1))" \
                "$((y_par2_end + 1))" "$y_length" > "$nonpar_bed"
            echo "Using GRCh$chrY_build Y-PAR boundaries ($build_source): $vcf"
            echo "MAKE pgen chrY (non-PAR only)"
            run_make_pfile "$output" \
                --vcf "$vcf" --double-id --update-sex "$sex_update" --chr Y \
                --extract bed1 "$nonpar_bed" --set-invalid-haploid-missing \
                --set-missing-var-ids '@:#:$r:$a' \
                --new-id-max-allele-len 100 missing
        )
    fi
    validate_pvar_chromosomes "$output" Y
    validate_y_nonpar "$output"
}

make_vcf_pgen_one() {
    local chr=$1 vcf=${input_prefix}${1}.vcf.gz output=$pfile_dir/chr${1}
    local args=(--vcf "$vcf" --double-id --update-sex "$sex_update"
        --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing)
    if [[ $chr == X ]]; then make_chrX_pfiles; return; fi
    if [[ $chr == Y ]]; then make_chrY_pfile; return; fi
    if pfile_zst_exists "$output" && ! is_true "$replace"; then
        echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
        return
    fi
    [[ -s $vcf ]] || die "missing VCF: $vcf"
    echo "MAKE pgen chr$chr"
    run_make_pfile "$output" "${args[@]}"
}

make_bgen_pgen_one() {
    local chr=$1 output=$pfile_dir/chr${1} bgen sample
    bgen=$(bgen_file_for_chr "$chr")
    sample=$(bgen_sample_file "$bgen")
    if pfile_zst_exists "$output" && ! is_true "$replace"; then
        echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
    else
        echo "MAKE pgen chr$chr from BGEN"
        run_make_pfile "$output" --bgen "$bgen" ref-first --sample "$sample" \
            --oxford-single-chr "$chr"
    fi
}

make_bed_pgen_one() {
    local chr=$1 input=${input_prefix}${1} output=$pfile_dir/chr${1}
    [[ -s ${input}.bed && -s ${input}.bim && -s ${input}.fam ]] ||
        die "incomplete BED/BIM/FAM input fileset: $input"
    if pfile_zst_exists "$output" && ! is_true "$replace"; then
        echo "SKIP pgen: $output already has .pvar.zst (use --replace TRUE)"
    else
        echo "MAKE pgen chr$chr from BED/BIM/FAM"
        run_make_pfile "$output" --bfile "$input" --indiv-sort none
    fi
}

make_pgen_one() {
    local chr=$1 output=$pfile_dir/chr${1}
    case $input_format in
        vcf) make_vcf_pgen_one "$chr"; return ;;
        bgen) make_bgen_pgen_one "$chr" ;;
        bed) make_bed_pgen_one "$chr" ;;
    esac
    case $chr in
        X) validate_pvar_chromosomes "$output" X
           if is_true "$split_sex"; then make_chrX_sex_pfiles; fi ;;
        XY) validate_pvar_chromosomes "$output" PAR ;;
        Y) validate_pvar_chromosomes "$output" Y ;;
    esac
}

make_bed_prefix() {
    local label=$1 input=$2 tmp
    local pfile_args=()
    tmp=$pfile_dir/.${label}.make-bed.$$
    pfile_exists "$input" || die "missing pfile for $label: $input (request --out-files pgen,bed to create it)"
    if bed_exists "$input" && ! is_true "$replace"; then
        echo "SKIP bed: $input already exists (use --replace TRUE)"
        return
    fi
    rm -f -- "$tmp".bed "$tmp".bim "$tmp".fam "$tmp".log "$tmp".nosex
    echo "MAKE bed  $label (multiallelic variants excluded)"
    pfile_input_args "$input" pfile_args
    "$plink2_bin" "${pfile_args[@]}" \
        --max-alleles 2 --output-chr 26 --make-bed --out "$tmp"
    [[ -s $tmp.bed && -s $tmp.bim && -s $tmp.fam ]] || die "incomplete bed output for $label"
    mv -f -- "$tmp.bed" "$input.bed"
    mv -f -- "$tmp.bim" "$input.bim"
    mv -f -- "$tmp.fam" "$input.fam"
    rm -f -- "$tmp.log" "$tmp.nosex"
}

make_bed_one() {
    local chr=$1
    if [[ $chr == X ]]; then
        validate_pvar_chromosomes "$pfile_dir/chrX" X
        make_bed_prefix chrX "$pfile_dir/chrX"
        if pfile_exists "$pfile_dir/chrXY"; then
            validate_pvar_chromosomes "$pfile_dir/chrXY" PAR
            make_bed_prefix chrXY "$pfile_dir/chrXY"
        else
            echo "SKIP bed chrXY: no chrXY pfile (source chrX VCF had no PAR variants)"
        fi
        if is_true "$split_sex"; then
            local sex
            for sex in male female; do
                validate_pvar_chromosomes "$pfile_dir/chrX.$sex" X
                make_bed_prefix "chrX.$sex" "$pfile_dir/chrX.$sex"
            done
        fi
    else
        make_bed_prefix "chr$chr" "$pfile_dir/chr$chr"
    fi
}


# 🚩 Genotype conversions
run_make_module() {
    local value
    local -A requested=()
    local outputs=()
    split_csv "$out_files_csv" outputs
    for value in "${outputs[@]}"; do
        value=${value,,}; [[ $value == pfile ]] && value=pgen
        [[ $value == pgen || $value == bed ]] || die "invalid --out-files value: $value"
        requested[$value]=1
    done

    if [[ -n $chrs_csv ]]; then
        normalize_chromosomes "$chrs_csv"
    elif [[ -n ${requested[pgen]:-} ]]; then
        discover_input_chromosomes
    else
        discover_pfile_chromosomes
    fi
    if is_true "$split_sex"; then
        local found_x=0
        for value in "${selected_chrs[@]}"; do [[ $value == X ]] && found_x=1; done
        ((found_x)) || die "--split-sex TRUE requires chrX to be selected"
        if [[ -z ${requested[pgen]:-} ]] &&
           { ! pfile_zst_exists "$pfile_dir/chrX.male" ||
             ! pfile_zst_exists "$pfile_dir/chrX.female"; }; then
            echo "sex-split chrX pfiles are missing; adding pgen generation before requested bed output"
            requested[pgen]=1
        fi
    fi
    echo "make-pfile input-format=$input_format input-prefix=$input_prefix dir-input=$input_dir dir-pfile=$pfile_dir chr=$(IFS=,; echo "${selected_chrs[*]}")"

    sex_update=""
    if [[ $input_format == vcf ]]; then
        sex_update=$pfile_dir/ALL.update-sex.txt
        write_sex_update "$sex_update"
    fi
    if [[ -n ${requested[pgen]:-} ]]; then
        for value in "${selected_chrs[@]}"; do make_pgen_one "$value"; done
    fi
    if [[ -n ${requested[bed]:-} ]]; then
        for value in "${selected_chrs[@]}"; do make_bed_one "$value"; done
    fi
    echo "DONE: make-pfile under $pfile_dir"
}

download_archaic_pair() {
    local destination=$1 url=$2
    wget -c -P "$destination" "$url" "${url}.tbi"
}


# 🚩 Archaic references
run_download_archaic_module() {
    local chr f old_file repaired_file command_name

    for command_name in wget gzip gunzip bgzip tabix find xargs; do
        command -v "$command_name" >/dev/null 2>&1 || die "missing command: $command_name"
    done
    [[ $download_jobs =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
    [[ $bgzip_threads =~ ^[1-9][0-9]*$ ]] || die "--bgzip-threads must be a positive integer"

    if [[ -n $chrs_csv ]]; then
        normalize_chromosomes "$chrs_csv"
        for chr in "${selected_chrs[@]}"; do
            [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ || $chr == X ]] ||
                die "download_archaic supports chromosomes 1-22 and X, got: $chr"
        done
    else
        selected_chrs=({1..22} X)
    fi

    mkdir -p "$archaic_root"/{Vindija,Altai,Chagyr,Denisova,Denisova25}
    echo "download_archaic root=$archaic_root chr=$(IFS=,; echo "${selected_chrs[*]}")"
    for chr in "${selected_chrs[@]}"; do
        download_archaic_pair "$archaic_root/Vindija" \
            "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${chr}_mq25_mapab100.vcf.gz"
        download_archaic_pair "$archaic_root/Altai" \
            "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${chr}_mq25_mapab100.vcf.gz"
        download_archaic_pair "$archaic_root/Chagyr" \
            "https://cdna.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${chr}.noRB.vcf.gz"
        download_archaic_pair "$archaic_root/Denisova" \
            "https://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${chr}_mq25_mapab100.vcf.gz"
        download_archaic_pair "$archaic_root/Denisova25" \
            "https://cdna.eva.mpg.de/denisova/Den25/VCF/chr${chr}.Den25.L35MQ25.B30.map35_100.vcf.gz"
    done

    while IFS= read -r -d '' f; do
        if ! gzip -t "$f" 2>/dev/null || ! tabix -l "$f" >/dev/null 2>&1; then
            old_file=${f}.old.gz
            repaired_file=${f}.repacked.$$
            mv -f -- "$f" "$old_file"
            if ! gunzip -c "$old_file" | bgzip -@ "$bgzip_threads" -c > "$repaired_file"; then
                rm -f -- "$repaired_file"
                mv -f -- "$old_file" "$f"
                die "failed to repack Chagyrskaya VCF: $f"
            fi
            mv -f -- "$repaired_file" "$f"
            echo "REPACKED: $f (original: $old_file)"
        fi
        tabix -f -p vcf "$f"
    done < <(find "$archaic_root/Chagyr" -maxdepth 1 -type f -name '*.vcf.gz' -print0)

    find "$archaic_root" -mindepth 2 -maxdepth 2 -type f -name '*.vcf.gz' -print0 |
        xargs -0 -r -P "$download_jobs" -I {} tabix -f -p vcf '{}'
    echo "DONE: download_archaic under $archaic_root"
}

run_lift_bed_module() {
    local input rel output output_dir unmapped_output unmapped_dir
    local tmp_root chain_for_liftover input_bed mapped_bed unmapped_bed
    local output_tmp unmapped_tmp chrom_style_file chrom_style
    local -a bed_files=()
    local command_name

    for command_name in awk bgzip find gzip grep liftOver mktemp sort; do
        command -v "$command_name" >/dev/null 2>&1 || die "missing command: $command_name"
    done
    [[ $bgzip_threads =~ ^[1-9][0-9]*$ ]] ||
        die "--bgzip-threads must be a positive integer"
    [[ -d $dir_in ]] || die "input directory does not exist: $dir_in"
    [[ -s $chain_file ]] || die "chain file does not exist or is empty: $chain_file"

    dir_in=${dir_in%/}
    dir_out=${dir_out%/}
    [[ $dir_out != "$dir_in" && $dir_out != "$dir_in"/* ]] ||
        die "--dir-out must not be --dir-in or a directory below it"

    mapfile -d '' bed_files < <(
        find "$dir_in" -type f -name '*.bed.gz' -print0 | sort -z
    )
    ((${#bed_files[@]} > 0)) || die "no *.bed.gz files found below: $dir_in"

    mkdir -p "$dir_out"
    tmp_root=$(mktemp -d "$dir_out/.lift-bed.tmp.XXXXXX")
    trap '[[ -n ${tmp_root:-} && -d $tmp_root ]] && rm -rf -- "$tmp_root"' EXIT
    chain_for_liftover=$chain_file
    if [[ $chain_file == *.gz ]]; then
        chain_for_liftover=$tmp_root/map.over.chain
        gzip -cd -- "$chain_file" > "$chain_for_liftover"
        [[ -s $chain_for_liftover ]] || die "failed to decompress chain file: $chain_file"
    fi

    input_bed=$tmp_root/input.bed
    mapped_bed=$tmp_root/mapped.bed
    unmapped_bed=$tmp_root/unmapped.bed
    chrom_style_file=$tmp_root/chrom.style
    echo "lift-bed files=${#bed_files[@]} dir-in=$dir_in dir-out=$dir_out chain=$chain_file"

    for input in "${bed_files[@]}"; do
        rel=${input#"$dir_in"/}
        output=$dir_out/$rel
        output_dir=$(dirname -- "$output")
        unmapped_output=$dir_out/.unmapped/${rel%.bed.gz}.unmapped.bed.gz
        unmapped_dir=$(dirname -- "$unmapped_output")

        if [[ -s $output ]] && gzip -t -- "$output" 2>/dev/null && ! is_true "$replace"; then
            echo "SKIP lift-bed: $rel already exists (use --replace TRUE)"
            continue
        fi

        mkdir -p "$output_dir" "$unmapped_dir"
        : > "$chrom_style_file"
        # UCSC chains use chr1-style names, while these archaic BED masks use
        # 1-style names. Record the input convention, add chr for liftOver,
        # then restore the original convention in both result files.
        gzip -cd -- "$input" |
            awk -v style_file="$chrom_style_file" 'BEGIN {OFS="\t"}
                /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ {print; next}
                NF == 0 {next}
                NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 > $3 {
                    print "invalid BED row at input line " NR > "/dev/stderr"
                    bad=1; next
                }
                {
                    if (!seen++) print ($1 ~ /^chr/ ? "chr" : "bare") > style_file
                    if ($1 !~ /^chr/) $1="chr" $1
                    print
                }
                END {exit bad}
            ' > "$input_bed"
        [[ -s $chrom_style_file ]] || die "no BED records found in: $input"
        chrom_style=$(<"$chrom_style_file")

        : > "$mapped_bed"
        : > "$unmapped_bed"
        echo "LIFT: $rel"
        liftOver "$input_bed" "$chain_for_liftover" "$mapped_bed" "$unmapped_bed"

        output_tmp=$output.tmp.$$
        if [[ $chrom_style == bare ]]; then
            awk 'BEGIN {OFS="\t"}
                /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ {print; next}
                {$1=substr($1, 4); print}
            ' "$mapped_bed" | bgzip -@ "$bgzip_threads" -c > "$output_tmp"
        else
            bgzip -@ "$bgzip_threads" -c "$mapped_bed" > "$output_tmp"
        fi
        gzip -t -- "$output_tmp" || die "invalid compressed output for: $rel"
        mv -f -- "$output_tmp" "$output"

        if grep -qv '^#' "$unmapped_bed"; then
            unmapped_tmp=$unmapped_output.tmp.$$
            if [[ $chrom_style == bare ]]; then
                awk 'BEGIN {OFS="\t"}
                    /^#/ {print; next}
                    {$1=substr($1, 4); print}
                ' "$unmapped_bed" | bgzip -@ "$bgzip_threads" -c > "$unmapped_tmp"
            else
                bgzip -@ "$bgzip_threads" -c "$unmapped_bed" > "$unmapped_tmp"
            fi
            gzip -t -- "$unmapped_tmp" || die "invalid unmapped output for: $rel"
            mv -f -- "$unmapped_tmp" "$unmapped_output"
        else
            rm -f -- "$unmapped_output"
        fi
    done

    rm -rf -- "$tmp_root"
    tmp_root=""
    trap - EXIT
    echo "DONE: lift-bed under $dir_out"
}

[[ $# -gt 0 ]] || { usage >&2; exit 2; }
case $1 in
    -h|--help) usage; exit 0 ;;
    make-pfile|download_archaic|lift-bed) module=$1; shift ;;
    *) die "unknown module '$1'; expected make-pfile, download_archaic, or lift-bed" ;;
esac

grch=""
input_name_prefix=chr
input_format=auto
out_files_csv=pgen
chrs_csv=""
root=""
sample_file=""
dir_input=""
dir_vcf=""
dir_pfile=""
dir_in=""
dir_out=""
chain_file=""
input_dir=""
pfile_dir=""
replace=FALSE
split_sex=FALSE
download_jobs=8
bgzip_threads=4

while (($#)); do
    case $1 in
        --grch) grch=${2:?missing value for --grch}; shift 2 ;;
        --input-prefix|--vcf-prefix) input_name_prefix=${2:?missing value for --input-prefix}; shift 2 ;;
        --input-format) input_format=${2:?missing value for --input-format}; input_format=${input_format,,}; shift 2 ;;
        --out-files) out_files_csv=${2:?missing value for --out-files}; shift 2 ;;
        --chr) chrs_csv=${2:?missing value for --chr}; shift 2 ;;
        --split-sex)
            split_sex=${2:?missing value for --split-sex}
            is_true "$split_sex" || true
            shift 2
            ;;
        --chrX-male|--chrx-male)
            die "--chrX-male was replaced by --split-sex"
            ;;
        --root) root=${2:?missing value for --root}; shift 2 ;;
        --jobs) download_jobs=${2:?missing value for --jobs}; shift 2 ;;
        --bgzip-threads) bgzip_threads=${2:?missing value for --bgzip-threads}; shift 2 ;;
        --sample-file) sample_file=${2:?missing value for --sample-file}; shift 2 ;;
        --dir-input) dir_input=${2:?missing value for --dir-input}; shift 2 ;;
        --dir-vcf) dir_vcf=${2:?missing value for --dir-vcf}; shift 2 ;;
        --dir-pfile) dir_pfile=${2:?missing value for --dir-pfile}; shift 2 ;;
        --dir-in) dir_in=${2:?missing value for --dir-in}; shift 2 ;;
        --dir-out) dir_out=${2:?missing value for --dir-out}; shift 2 ;;
        --chain-file) chain_file=${2:?missing value for --chain-file}; shift 2 ;;
        --replace) replace=${2:?missing value for --replace}; is_true "$replace" || true; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

if [[ $module == make-pfile ]]; then
    [[ -z $grch || $grch == 37 || $grch == 38 ]] ||
        die "make-pfile --grch must be 37 or 38 when specified"
    [[ -z $root ]] || die "make-pfile does not accept --root; use --dir-input and --dir-pfile"
    if [[ -n $dir_vcf ]]; then
        [[ -z $dir_input ]] || die "use only one of --dir-input and --dir-vcf"
        [[ $input_format == auto || $input_format == vcf ]] ||
            die "--dir-vcf conflicts with --input-format $input_format"
        dir_input=$dir_vcf
        input_format=vcf
    fi
    if [[ -z $dir_input && -z $dir_pfile ]]; then
        default_build=${grch:-37}
        dir_input=$dir0/data.BIG/refGen/1kg/$default_build/vcf
        dir_pfile=$dir0/data.BIG/refGen/1kg/$default_build/pfile
        input_format=vcf
    elif [[ -z $dir_input ]]; then
        dir_input=$(dirname -- "${dir_pfile%/}")/vcf
        [[ $input_format != auto ]] || input_format=vcf
    elif [[ -z $dir_pfile ]]; then
        dir_pfile=$(dirname -- "${dir_input%/}")/pfile
    fi
    input_dir=${dir_input%/}
    pfile_dir=${dir_pfile%/}
elif [[ $module == lift-bed ]]; then
    [[ -z $grch ]] || die "lift-bed does not accept --grch; the chain file determines direction"
    [[ -z $root && -z $dir_input && -z $dir_vcf && -z $dir_pfile && -z $sample_file ]] ||
        die "lift-bed accepts --dir-in, --dir-out, --chain-file, --bgzip-threads, and --replace"
    [[ -n $dir_in ]] || die "lift-bed requires --dir-in DIR"
    [[ -n $dir_out ]] || die "lift-bed requires --dir-out DIR"
    [[ -n $chain_file ]] || die "lift-bed requires --chain-file FILE"
    run_lift_bed_module
    exit 0
else
    [[ -z $grch ]] || die "download_archaic does not accept --grch"
    [[ -z $dir_input && -z $dir_vcf && -z $dir_pfile && -z $sample_file ]] ||
        die "download_archaic does not accept VCF, pfile, or sample-file options"
    archaic_root=${root:-$dir0/data.BIG/refGen/archaic}
    selected_chrs=()
    run_download_archaic_module
    exit 0
fi
phe_sh=${phe_sh:-$dirfunc/0phe.f.sh}
plink2_bin=${PLINK2_BIN:-plink2}

command -v "$plink2_bin" >/dev/null 2>&1 || die "missing PLINK2 executable: $plink2_bin"
command -v awk >/dev/null 2>&1 || die "missing command: awk"
[[ -s $phe_sh ]] || die "missing chrX definitions: $phe_sh"
# shellcheck source=/dev/null
source "$phe_sh"
mkdir -p "$pfile_dir"

if [[ $input_name_prefix == */* ]]; then input_prefix=$input_name_prefix; else input_prefix=$input_dir/$input_name_prefix; fi
detect_input_format
if [[ $input_format == vcf ]]; then
    sample_file=${sample_file:-$(dirname -- "$input_dir")/samples.txt}
    [[ -s $sample_file ]] || die "missing sample metadata: $sample_file"
    IFS=$' \t' read -r -a header < <(head -n 1 "$sample_file" | tr -d '\r')
    expected=(sample fatherID motherID sex pop super_pop)
    [[ "${header[*]}" == "${expected[*]}" ]] ||
        die "unexpected header in $sample_file; expected: ${expected[*]}"
fi

selected_chrs=()


# 🚩 Module dispatch
case $module in
    make-pfile) run_make_module ;;
esac

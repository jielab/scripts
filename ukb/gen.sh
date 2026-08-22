#!/bin/bash -l



set +H
set -e

usage() {
    cat <<'USAGE'
Usage:
  ./gen.sh [options]
  ./gen.sh [TASK]
  ./gen.sh snp --dir-gen DIR --data typ|imp --snp-file FILE \
      --snp-file-type bed|snp [--flank BP] --dir-out DIR --label LABEL

Options:
  -h, --h, --help        Show this help message.

SNP options:
  --dir-gen DIR          Genetic data root containing typ/ and imp/.
  --data typ|imp         Genotyped PLINK bed files or imputed PLINK 2 pfiles.
  --snp-file FILE        Four-column BED file or one-variant-ID-per-line list.
  --snp-file-type TYPE   bed (0-based, half-open) or snp.
  --flank BP             BED only: add BP on each side. Default: 10000.
  --dir-out DIR          Output directory.
  --label LABEL          Output prefix; final data are LABEL.raw.
  --dry-run              Write the extraction script without running it.

Tasks:
  format    Prepare format commands for chr1-22,X.
  prs       Submit PRS jobs and write step2.cmd.
  snp       Submit SNP extraction job.
  apoe      Submit APOE haplotype extraction job.
  all       Run all task blocks in order. Default.

Examples:
  ./gen.sh
  task=format ./gen.sh
  TASK=prs ./gen.sh
  ./gen.sh format
  ./gen.sh snp --dir-gen /mnt/d/data/ukb/gen/ --data typ --snp-file /mnt/d/files/Wang_MY.bed --snp-file-type bed --dir-out /mnt/d/analysis/ukb/ --label Wang
USAGE
}

init_mnt_g_env() {
    dir0=/mnt/g
    dirimp=/mnt/h/ukbGen/imp
    dirout=$dir0/data-ukb-gen/imp

    if [[ -f "$dir0/scripts/0f/0phe.f.sh" ]]; then
        source "$dir0/scripts/0f/0phe.f.sh"
    else
        echo "ERROR: cannot find $dir0/scripts/0f/0phe.f.sh"
        exit 1
    fi
}

run_format() {
    init_mnt_g_env
    mkdir -p "$dirout"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 🚩 Format raw genotype data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for chr in {1..22} X; do
        cat > "$dirout/chr${chr}.cmd" <<EOF2
#!/bin/bash -l
set +H

plink2 --threads 40 --pfile $dirimp/chr$chr --make-pgen --out chr$chr
plink2 --memory 12000 --threads 40 --pfile chr$chr --freq counts --out chr$chr
awk '/^#/{next} \$3!="." && \$3!="" {if (\$5 ~ /,/) print \$3; if (seen[\$3]++) print \$3}' chr${chr}.pvar | sort -u > chr${chr}.bad.snp

# plink2 --memory 12000 --threads 40 --pfile chr$chr --export vcf bgz id-paste=iid --out chr$chr; tabix -f -p vcf chr$chr.vcf.gz
ln -sf chr${chr}.pgen chr${chr}.b38.pgen
ln -sf chr${chr}.psam chr${chr}.b38.psam

awk 'NR>1{print "chr"\$1,\$2-1,\$2,\$3}' chr$chr.pvar > chr$chr.tolift
liftOver chr$chr.tolift $dir0/files/liftOver/hg19ToHg38.over.chain.gz chr$chr.lifted chr$chr.unmapped
awk '\$0 !~ /^#/ && NF>=4{print \$4}' chr$chr.unmapped > chr$chr.b38.pvar.excl.tmp
awk -v excl=chr$chr.b38.pvar.excl.tmp 'BEGIN{OFS="\t"} FNR==NR{c=\$1; sub(/^chr/,"",c); map[\$4]=c OFS \$3; n[\$4]++; next} FNR==1{print; next} {if(\$3 in map && n[\$3]==1){split(map[\$3],a,OFS); \$1=a[1]; \$2=a[2]} else print \$3 >> excl; print}' chr$chr.lifted chr$chr.pvar > chr$chr.b38.pvar
sort -u chr$chr.b38.pvar.excl.tmp > chr$chr.b38.pvar.excl
n_bad=\$(paste chr$chr.pvar chr$chr.b38.pvar | awk 'NR>1 && \$3!=\$8{n++} END{print n+0}')
[[ "\$n_bad" == "0" ]] || { echo ERROR_ID_order_changed_n_bad=\$n_bad; exit 1; }
rm -f chr$chr.tolift chr$chr.lifted chr$chr.unmapped chr$chr.b38.pvar.excl.tmp
EOF2
        chmod +x "$dirout/chr${chr}.cmd"
    done

    cat > "$dirout/gen.step1.cmd" <<'EOF2'
#!/bin/bash -l
set +H
printf "%s\n" {1..22} X | xargs -P 8 -I{} sh -c './chr{}.cmd'
EOF2
    chmod +x "$dirout/gen.step1.cmd"

    echo "OK: format commands written to: $dirout"
    echo "Next: cd $dirout && bash gen.step1.cmd"
}

run_prs() {
    init_mnt_g_env
    mkdir -p "$dir0/analysis/prs"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 🚩 Calculate PRS
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for label in bmi bmi.2k bald cad.2m cad.6m afib.6m; do
        echo "RUN $label"
        prs_ref=$gwasdir/$label.ref
        dirout=$dir0/analysis/prs/$label

        if [[ -d "$dirout" ]]; then
            echo "$label already run; skip existing folder: $dirout"
            continue
        fi
        mkdir -p "$dirout"

        for var in SNP SNP_col EA EA_col NEA NEA_col N N_col BETA BETA_col P P_col CHR CHR_col; do
            eval "$var=''"
        done
        header_names "$prs_ref"
        cols="$SNP_col $EA_col $BETA_col header"
        echo "score columns: $cols"

        cat > "$dirout/prs.cmd" <<EOF2
#!/bin/bash -l
set +H

for chr in {1..22} X; do
    cat $prs_ref | awk -v snp=$SNP_col -v chr=\$chr 'NR==1 {print; next} \$$CHR_col==chr && !seen[\$snp]++' > chr\$chr.ref
    nref=\$(wc -l chr\$chr.ref | awk '{printf \$1}')
    if [[ \$nref == 1 ]]; then continue; fi
    plink2 --pfile $gendir/chr\$chr --chr \$chr --score chr\$chr.ref $cols no-mean-imputation ignore-dup-ids cols=+scoresums list-variants --out chr\$chr
 done

paste -d ' ' chr*.sscore > $label.prs.tmp
awk '{print NF}' $label.prs.tmp | sort -nu
num=\$(ls -l chr*.sscore | wc -l | awk '{printf \$1}')
awk -v num=\$num '{if (NR==1) print "eid $label.allele_cnt $label.score_sum"; else {allele_cnt=score_sum=0; for (i=1;i<=num;i++) {if (\$(i*6-5) !=\$1) score_sum=score_sum""i"ERR,"; else {allele_cnt=allele_cnt+\$(i*6-3); score_sum=score_sum+\$(i*6-0)}}; print \$1, allele_cnt, score_sum} }' $label.prs.tmp > $label.prs.txt
fgrep ERR $label.prs.txt || true
EOF2
        chmod +x "$dirout/prs.cmd"

        (cd "$dirout" && bsub -q short -J "$label.prs" -o prs.LOG -e prs.ERR < prs.cmd)
    done

    cat > "$dir0/analysis/prs/step2.cmd" <<'EOF2'
#!/bin/bash -l
set +H
set -e

cd "$(dirname "$0")"
egrep -i "error|ERR" */*.LOG */*.log */*.ERR */*.err 2>/dev/null || true
paste -d ' ' */*.prs.txt > all.prs.txt
awk '{print NF}' all.prs.txt | sort -nu
num=$(ls -l */*prs.txt | wc -l | awk '{printf $1}')
awk -v num=$num '{for (i=1;i<=num;i++) {if ($(i*3-2) !=$1) print $1, $i, "ERR" }}' all.prs.txt | head
awk -v num=$num '{for (i=2;i<=num;i++) {$(i*3-2)=""} print $0}' all.prs.txt | sed 's/  / /g' > all.prs
gzip -f all.prs
EOF2
    chmod +x "$dir0/analysis/prs/step2.cmd"

    echo "OK: PRS jobs submitted. After all jobs finish: cd $dir0/analysis/prs && bash step2.cmd"
}

run_snp() {
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 🚩 Extract SNP data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        extract_args="--extract \"$snp_file\""
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
            [[ -f "$dir_gen/imp/chr$chr.pgen" && -f "$dir_gen/imp/chr$chr.pvar" && -f "$dir_gen/imp/chr$chr.psam" ]] || {
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
        input_args="--pfile \"$dir_gen/imp/chr\$chr\""
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 🚩 Extract APOE haplotype
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        -h|--h|--help) usage; exit 0 ;;
        --dir-gen) dir_gen=${2:?ERROR: --dir-gen needs a value}; shift 2 ;;
        --data) snp_data=${2:?ERROR: --data needs a value}; shift 2 ;;
        --snp-file) snp_file=${2:?ERROR: --snp-file needs a value}; shift 2 ;;
        --snp-file-type) snp_file_type=${2:?ERROR: --snp-file-type needs a value}; shift 2 ;;
        --flank) flank=${2:?ERROR: --flank needs a value}; flank_set=Y; shift 2 ;;
        --dir-out) dir_out=${2:?ERROR: --dir-out needs a value}; shift 2 ;;
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

task=${task:-${TASK:-${pos_args[0]:-all}}}
task=$(echo "$task" | tr '[:upper:]' '[:lower:]')

case "$task" in
    format) run_format ;;
    prs)    run_prs ;;
    snp)    run_snp ;;
    apoe)   run_apoe ;;
    all)    run_format; run_prs; run_snp; run_apoe ;;
    *)
        echo "ERROR: unknown task: $task" >&2
        usage >&2
        exit 1
        ;;
esac

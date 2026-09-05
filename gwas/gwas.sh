#!/usr/bin/env bash
set -eo pipefail


set +H

usage() {
    cat <<'EOF'
Usage:
  ./gwas.sh [options]

Options:
  --start-step STEP       Start from STEP. Default: prep_gwas.
  --end-step STEP         Stop after STEP.
  --steps CSV             Run only listed steps.
  --threads N             Worker threads. Default: 40.
  --data-root DIR         Data/analysis root. Default: /work/sph-huangj.
  --typed-dir DIR         PLINK typed-genotype directory.
  --imputed-dir DIR       Imputed BGEN directory.
  --saige-dir DIR         SAIGE extdata directory.
  --shared-shell FILE     Shared 0phe.f.sh library.
  --label NAME            Analysis label. Default: vip.
  --type TYPE             qt, bt, or t2e. Default: t2e.
  --phenotypes CSV        Phenotype names.
  --covariates CSV        Covariate names.
  --categorical-covariates CSV  Categorical covariates. Default: center.
  --phenotype-file FILE   Phenotype input.
  --covariate-file FILE   Covariate input.
  --output-dir DIR        GWAS output root.
  --chromosomes LIST      Comma- or space-separated chromosomes.
  --submit-prep TRUE|FALSE
  --submit-plink TRUE|FALSE
  --submit-regenie TRUE|FALSE
  --submit-saige TRUE|FALSE
  --skip-existing-regenie TRUE|FALSE
  --skip-existing-saige TRUE|FALSE
  --ldsc-gwas-dir DIR     LDSC input GWAS directory.
  --ldsc-output-dir DIR   LDSC output directory.
  --ldsc-ref-dir DIR      LDSC reference directory.
  --ldsc-software-dir DIR LDSC software directory.
  --ldsc-jobs N           Parallel LDSC jobs. Default: 8.
  --ldsc-pattern GLOB     LDSC input filename pattern.
  --ldsc-n-default N      Sample size used when N is absent.
  --run-ldsc-munge TRUE|FALSE
  --run-ldsc-rg TRUE|FALSE
  -h, --h, --help        Show this help message.

GWAS steps:
  prep_gwas plink_gwas regenie_gwas saige_gwas ldsc_gwas

Examples:
  ./gwas.sh --start-step prep_gwas --end-step regenie_gwas
  ./gwas.sh --steps prep_gwas,regenie_gwas
  ./gwas.sh --steps ldsc_gwas --run-ldsc-munge TRUE --run-ldsc-rg TRUE
EOF
}

bool01() {
    case "${1,,}" in
        true|1|yes|y) printf '1\n' ;;
        false|0|no|n) printf '0\n' ;;
        *) return 2 ;;
    esac
}

threads=40
dir0=/work/sph-huangj
dirtyp=""
dirimp=""
dirsaige=""
script_f=""
label=vip
phe_f=""
cov_f=""
type=t2e
phes=ihd,stroke,stroke_h,stroke_i,stroke_o,ihd.99,stroke.99,stroke_h.99,stroke_i.99,stroke_o.99
covs=age,sex,tdi,PC1,PC2,PC3,PC4
cat_covs=center
chr_list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
dirout=""
START_STEP=prep_gwas
END_STEP=""
RUN_STEPS=""
SUBMIT_PREP=0
SUBMIT_PLINK=0
SUBMIT_REGENIE=1
SUBMIT_SAIGE=0
SKIP_EXISTING_REGENIE=1
SKIP_EXISTING_SAIGE=0
ldsc_in=""
ldsc_out=""
ldsc_ref=""
ldsc_soft=""
ldsc_jobs=8
ldsc_n_default=100000
ldsc_pattern=""
run_ldsc_munge=1
run_ldsc_rg=1

while (( $# )); do
    case "$1" in
        --start-step) START_STEP=${2:?ERROR: --start-step requires STEP}; shift 2 ;;
        --end-step) END_STEP=${2:?ERROR: --end-step requires STEP}; shift 2 ;;
        --steps) RUN_STEPS=${2:?ERROR: --steps requires CSV}; shift 2 ;;
        --threads) threads=${2:?ERROR: --threads requires N}; shift 2 ;;
        --data-root) dir0=${2:?ERROR: --data-root requires DIR}; shift 2 ;;
        --typed-dir) dirtyp=${2:?ERROR: --typed-dir requires DIR}; shift 2 ;;
        --imputed-dir) dirimp=${2:?ERROR: --imputed-dir requires DIR}; shift 2 ;;
        --saige-dir) dirsaige=${2:?ERROR: --saige-dir requires DIR}; shift 2 ;;
        --shared-shell) script_f=${2:?ERROR: --shared-shell requires FILE}; shift 2 ;;
        --label) label=${2:?ERROR: --label requires NAME}; shift 2 ;;
        --type) type=${2:?ERROR: --type requires TYPE}; shift 2 ;;
        --phenotypes) phes=${2:?ERROR: --phenotypes requires CSV}; shift 2 ;;
        --covariates) covs=${2:?ERROR: --covariates requires CSV}; shift 2 ;;
        --categorical-covariates) cat_covs=${2:?ERROR: --categorical-covariates requires CSV}; shift 2 ;;
        --phenotype-file) phe_f=${2:?ERROR: --phenotype-file requires FILE}; shift 2 ;;
        --covariate-file) cov_f=${2:?ERROR: --covariate-file requires FILE}; shift 2 ;;
        --output-dir) dirout=${2:?ERROR: --output-dir requires DIR}; shift 2 ;;
        --chromosomes) chr_list=${2:?ERROR: --chromosomes requires LIST}; chr_list=${chr_list//,/ }; shift 2 ;;
        --submit-prep) SUBMIT_PREP=$(bool01 "${2:?ERROR: --submit-prep requires TRUE or FALSE}") || { echo "ERROR: --submit-prep requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --submit-plink) SUBMIT_PLINK=$(bool01 "${2:?ERROR: --submit-plink requires TRUE or FALSE}") || { echo "ERROR: --submit-plink requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --submit-regenie) SUBMIT_REGENIE=$(bool01 "${2:?ERROR: --submit-regenie requires TRUE or FALSE}") || { echo "ERROR: --submit-regenie requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --submit-saige) SUBMIT_SAIGE=$(bool01 "${2:?ERROR: --submit-saige requires TRUE or FALSE}") || { echo "ERROR: --submit-saige requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --skip-existing-regenie) SKIP_EXISTING_REGENIE=$(bool01 "${2:?ERROR: --skip-existing-regenie requires TRUE or FALSE}") || { echo "ERROR: --skip-existing-regenie requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --skip-existing-saige) SKIP_EXISTING_SAIGE=$(bool01 "${2:?ERROR: --skip-existing-saige requires TRUE or FALSE}") || { echo "ERROR: --skip-existing-saige requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --ldsc-gwas-dir) ldsc_in=${2:?ERROR: --ldsc-gwas-dir requires DIR}; shift 2 ;;
        --ldsc-output-dir) ldsc_out=${2:?ERROR: --ldsc-output-dir requires DIR}; shift 2 ;;
        --ldsc-ref-dir) ldsc_ref=${2:?ERROR: --ldsc-ref-dir requires DIR}; shift 2 ;;
        --ldsc-software-dir) ldsc_soft=${2:?ERROR: --ldsc-software-dir requires DIR}; shift 2 ;;
        --ldsc-jobs) ldsc_jobs=${2:?ERROR: --ldsc-jobs requires N}; shift 2 ;;
        --ldsc-pattern) ldsc_pattern=${2:?ERROR: --ldsc-pattern requires GLOB}; shift 2 ;;
        --ldsc-n-default) ldsc_n_default=${2:?ERROR: --ldsc-n-default requires N}; shift 2 ;;
        --run-ldsc-munge) run_ldsc_munge=$(bool01 "${2:?ERROR: --run-ldsc-munge requires TRUE or FALSE}") || { echo "ERROR: --run-ldsc-munge requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --run-ldsc-rg) run_ldsc_rg=$(bool01 "${2:?ERROR: --run-ldsc-rg requires TRUE or FALSE}") || { echo "ERROR: --run-ldsc-rg requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done


# 🚩 Shell options, paths, and shared inputs
dirtyp=${dirtyp:-$dir0/data/ukb/gen/typ}
dirimp=${dirimp:-$dir0/data/ukb/gen/imp}
dirsaige=${dirsaige:-$dir0/software/SAIGE/extdata}
script_f=${script_f:-$dir0/scripts/0f/0phe.f.sh}
phe_f=${phe_f:-$dir0/data/ukb/phe/common/ukb.phe}
cov_f=${cov_f:-$phe_f}
phes_space=$(echo "$phes" | sed 's/,/ /g')
dirout=${dirout:-$dir0/data/gwas/self/$label}
ldsc_in=${ldsc_in:-$dirout/regenie}
ldsc_out=${ldsc_out:-$dir0/analysis/ldsc/$label.$type}
ldsc_ref=${ldsc_ref:-$dir0/data/ldref/ldsc}
ldsc_soft=${ldsc_soft:-$dir0/software/ldsc}
ldsc_pattern=${ldsc_pattern:-*.$type.gz}
mkdir -p "$dirout"

steps=(prep_gwas plink_gwas regenie_gwas saige_gwas ldsc_gwas)

all_steps=(all "${steps[@]}")
contains() { local x="$1"; shift; for y in "$@"; do [[ "$x" == "$y" ]] && return 0; done; return 1; }

if ! contains "$START_STEP" "${all_steps[@]}"; then
    echo "ERROR: unknown START_STEP: $START_STEP" >&2
    usage >&2
    exit 1
fi
if [[ -n "$END_STEP" ]] && ! contains "$END_STEP" "${all_steps[@]}"; then
    echo "ERROR: unknown END_STEP: $END_STEP" >&2
    usage >&2
    exit 1
fi
if [[ -n "$RUN_STEPS" ]]; then
    IFS=',' read -ra _run_steps_check <<< "$RUN_STEPS"
    for _st in "${_run_steps_check[@]}"; do
        _st="${_st// /}"
        [[ -z "$_st" ]] && continue
        if ! contains "$_st" "${all_steps[@]}"; then
            echo "ERROR: unknown RUN_STEPS entry: $_st" >&2
            usage >&2
            exit 1
        fi
    done
fi

idx_of() { local x="$1"; shift; local i=0; for y in "$@"; do [[ "$x" == "$y" ]] && { echo "$i"; return 0; }; i=$((i+1)); done; echo -1; }
in_run_steps() { [[ -z "$RUN_STEPS" ]] && return 1; [[ ",${RUN_STEPS// /,}," == *",$1,"* ]]; }
should_run_step() {
    local step="$1"
    if [[ -n "$RUN_STEPS" ]]; then in_run_steps "$step"; return $?; fi
    [[ "$START_STEP" == "all" ]] && return 0
    contains "$START_STEP" "${steps[@]}" || return 1
    local i s e
    i=$(idx_of "$step" "${steps[@]}"); s=$(idx_of "$START_STEP" "${steps[@]}")
    if [[ -n "$END_STEP" ]] && contains "$END_STEP" "${steps[@]}"; then e=$(idx_of "$END_STEP" "${steps[@]}"); else e=$((${#steps[@]}-1)); fi
    [[ $i -ge $s && $i -le $e ]]
}
run_step() {
    local step="$1"; shift
    if should_run_step "$step"; then echo; echo "========== RUN STEP: $step =========="; "$@"; echo "========== DONE STEP: $step =========="; else echo "Skip step: $step"; fi
}

if [[ "$type" == bt ]]; then
    plink2_str="+a1freqcc,+a1countcc,+totallelecc"
    regenie_str="--$type"
    saige_str="--traitType=binary"
elif [[ "$type" == qt ]]; then
    plink2_str="+a1freq,+a1count,+totallelecc"
    regenie_str="--$type --apply-rint"
    saige_str="--traitType=quantitative --invNormalize=TRUE"
elif [[ "$type" == t2e ]]; then
    plink2_str=""
    regenie_str="--$type"
    saige_str=""
else
    echo "ERROR: please specify TYPE as qt, bt, or t2e" >&2
    exit 1
fi
saige_inv=""


# 🚩 Prep GWAS input and sparse GRM
prep_gwas() {
    mkdir -p "$dirtyp"
    for chr in {1..22}; do
        echo "#!/bin/bash
        plink2 --threads $threads --bfile chr$chr --maf 0.01 --mac 100 --max-alleles 2 --geno 0.1 --hwe 1e-15 --mind 0.1 --indep-pairwise 1000kb 0.1 --out chr$chr
        plink2 --threads $threads --bfile chr$chr --extract chr$chr.prune.in --make-bed --out chr$chr.prune
        " > "$dirtyp/chr$chr.prune.cmd"
        if [[ "$SUBMIT_PREP" == "1" ]]; then
            cd "$dirtyp"
            bsub -q short -n "$threads" -J typ.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.prune.cmd
        fi
    done

    echo "#!/bin/bash
	module load python/anaconda3/2022.7; source activate; conda activate R4.5.2; conda activate SAIGE
	ls -v chr*.bed | sed -e 's/\.bed//g' > merge_list.txt
	cat chr*prune.in > typ.prune.in
	plink20 --threads $threads --merge-list merge_list.txt --make-bed --out typ
	plink2 --bfile typ --extract typ.prune.in --exclude range <(echo -e \"6\\t25000000\\t35000000\") --make-bed --out typ.prune
	Rscript $dirsaige/createSparseGRM.R --nThreads=$threads --plink2File=typ.prune --numRandomMarkerforSparseKin=2000 --relatednessCutoff=0.125 --outputPrefix=typ
    " > "$dirtyp/merge.cmd"
    if [[ "$SUBMIT_PREP" == "1" ]]; then
        cd "$dirtyp"
        bsub -q short -n "$threads" -J merge -w "ended(typ.chr*)" -o merge.LOG -e merge.ERR < merge.cmd
    fi
}
run_step prep_gwas prep_gwas


# 🚩 PLINK2 GWAS
plink_gwas() {
    if [[ "$type" == "t2e" ]]; then
        echo "Skip PLINK2 GWAS: TYPE=t2e is not supported by this plink2 --glm section. Use TYPE=bt or TYPE=qt."
        return 0
    fi

    local dirout2="$dirout/plink2/$type"
    mkdir -p "$dirout2"
    for chr in $chr_list; do
        x_str=""; if [[ "$chr" == X ]]; then x_str="--xchr-model 2"; fi
        echo "#!/bin/bash
	source $script_f
	if [[ ! -f $dirimp/chr$chr.extract ]]; then
		plink2 --memory 12000 --threads $threads --pfile $dirimp/chr$chr --freq counts --out $dirimp/chr$chr
		awk '\$5>100 {print \$2}' $dirimp/chr$chr.acount | sed '1 s/ID/SNP/' > $dirimp/chr$chr.extract
	fi

	plink2 --memory 12000 --threads $threads --1 --out chr$chr\\
		--pfile $dirimp/chr$chr.b38 --extract $dirimp/chr$chr.extract \\
		--glm cols=+omitted,+nobs,+beta,$plink2_str hide-covar allow-no-covars no-x-sex \\
		--pheno $phe_f --no-input-missing-phenotype --no-psam-pheno --pheno-name $phes $x_str \\
		--covar $cov_f --covar-name $covs
		
	for t in $phes_space; do 
		file=\$(ls -1 chr$chr.\$t.glm.* | awk '{printf \$1}')
		sed -i '1 s/^#//; 1 s/?//g' \$file
		phe_header_names \$file
		cat \$file | awk -v snp_col=\$SNP_col -v chr_col=\$CHR_col -v pos_col=\$POS_col -v ea_col=\$EA_col -v nea_col=\$NEA_col -v eaf_col=\$EAF_col -v n_col=\$N_col -v beta_col=\$BETA_col -v se_col=\$SE_col -v p_col=\$P_col '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print \$snp_col,\$chr_col,\$pos_col, \$ea_col,\$nea_col,\$eaf_col,\$n_col,\$beta_col,\$se_col,\$p_col}' | grep -vwE \"nan|inf|NA\" | sed 's/ /\t/g' | sort -k 2,2n -k 3,3n | gzip -f > chr$chr.\$t.gz
	done
        " > "$dirout2/chr$chr.cmd"
        if [[ "$SUBMIT_PLINK" == "1" ]]; then
            cd "$dirout2"
            bsub -q short -n "$threads" -J $label.$type.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
        fi
    done

    echo "#!/bin/bash
	for t in $phes_space; do
		mkdir -p $dirout2/\$t.$type
		ls -v chr*.\$t.gz | awk 'NR==1{system(\"zcat \" \$0); next} {system(\"zcat \" \$0 \" | sed 1d\")}' | gzip -f > \$t.$type.gz
		zcat \$t.$type.gz | awk 'NR==1 || \$NF <=0.001' | gzip -f > \$t.$type.001.gz
		mv chr*.\$t.* \$t.$type/
	done
    " > "$dirout2/sum.cmd"
    if [[ "$SUBMIT_PLINK" == "1" ]]; then
        cd "$dirout2"
        bsub -q short -n "$threads" -J $label.$type.sum -w "done($label.$type.chr*)" -o sum.LOG -e sum.ERR < sum.cmd
    fi
}
run_step plink_gwas plink_gwas


# 🚩 REGENIE GWAS
regenie_gwas() {
    for phe in $phes_space; do
        local dirout2="$dirout/regenie/$phe.$type"
        if [[ -d "$dirout2" && "$SKIP_EXISTING_REGENIE" == "1" ]]; then
            echo "$phe.$type already processed; skip. Set SKIP_EXISTING_REGENIE=0 to overwrite command files."
            continue
        fi
        mkdir -p "$dirout2"
        if [[ "$type" == bt || "$type" == qt ]]; then phe_str="--phenoCol $phe"; elif [[ "$type" == t2e ]]; then phe_str="--phenoColList $phe.t2e --eventColList $phe.Yt2e"; fi

        echo "#!/bin/bash
	module load python/anaconda3/2020.7 
	source /share/apps/anaconda3/2020.7/bin/activate; conda activate /work/sph-huangj/.conda/envs/regenie-4.1.1
	regenie --step 1 --out step1 --lowmem --lowmem-prefix tmp --threads $threads --bsize 1000 \\
		--bed $dirtyp/typ.prune --phenoFile $phe_f $phe_str $regenie_str --covarFile $cov_f --covarColList $covs --catCovarList $cat_covs --maxCatLevels 30
        " > "$dirout2/step1.cmd"
        if [[ "$SUBMIT_REGENIE" == "1" ]]; then
            cd "$dirout2"
            bsub -q short -n "$threads" -J $phe.step1 -o step1.LOG -e step1.ERR < step1.cmd
        fi

        for chr in $chr_list; do
            echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source /share/apps/anaconda3/2020.7/bin/activate; conda activate /work/sph-huangj/.conda/envs/regenie-4.1.1
	source $script_f
	regenie --step 2 --out chr$chr --pThresh 0.01 --pred step1_pred.list --firth --approx --bsize 1000 \\
		--pgen $dirimp/chr$chr.b38 --extract $dirimp/chr$chr.extract \\
		--phenoFile $phe_f $phe_str $regenie_str --covarFile $cov_f --covarColList $covs --catCovarList $cat_covs --maxCatLevels 30
	file=\$(ls -1 chr${chr}_$phe.*regenie | head -1)
	phe_header_names \$file
	cat \$file | awk -v snp_col=\$SNP_col -v chr_col=\$CHR_col -v pos_col=\$POS_col -v ea_col=\$EA_col -v nea_col=\$NEA_col -v eaf_col=\$EAF_col -v n_col=\$N_col -v beta_col=\$BETA_col -v se_col=\$SE_col -v p_col=\$P_col -v log10p_col=\$LOG10P_col '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else { if (p_col==\"\") pval = 10^(-\$log10p_col); else pval=\$p_col; print \$snp_col, \$chr_col, \$pos_col, \$ea_col, \$nea_col, \$eaf_col, \$n_col, \$beta_col, \$se_col, pval }}' | grep -vwE \"nan|inf|NA\" | sed 's/ /\t/g' | sort -k 2,2n -k 3,3n > chr$chr.$phe
	plink2 --memory 12000 --threads $threads --pfile $dirimp/chr$chr --chr $chr --clump chr$chr.$phe --clump-p1 5e-08 --clump-p2 5e-08 --clump-r2 0.01 --clump-kb 1000 --clump-id-field \$SNP --clump-p-field \$P --out chr$chr.$phe
	gzip -f chr$chr.$phe
            " > "$dirout2/chr$chr.cmd"
            if [[ "$SUBMIT_REGENIE" == "1" ]]; then
                cd "$dirout2"
                bsub -q short -n "$threads" -J $phe.chr$chr -w "done($phe.step1)" -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
            fi
        done

        echo "#!/bin/bash
	ls -v chr*.$phe.gz | awk 'NR==1{system(\"zcat \" \$0); next} {system(\"zcat \" \$0 \" | sed 1d\")}' | gzip -f > $phe.$type.gz
	zcat $phe.$type.gz | awk 'NR==1 || \$NF <=0.001' | gzip -f > $phe.$type.001.gz
        " > "$dirout2/sum.cmd"
        if [[ "$SUBMIT_REGENIE" == "1" ]]; then
            cd "$dirout2"
            bsub -q short -n "$threads" -J $phe.sum -w "done($phe.chr*)" -o sum.LOG -e sum.ERR < sum.cmd
        fi
    done
}
run_step regenie_gwas regenie_gwas


# 🚩 SAIGE GWAS
saige_gwas() {
    if [[ "$type" == "t2e" ]]; then
        echo "Skip SAIGE GWAS: SAIGE section is configured for TYPE=bt or TYPE=qt, not TYPE=t2e."
        return 0
    fi

    for phe in $phes_space; do
        local dirout2="$dirout/saige/$phe.$type"
        if [[ -d "$dirout2" && "$SKIP_EXISTING_SAIGE" == "1" ]]; then
            echo "$phe.$type SAIGE folder exists; skip. Set SKIP_EXISTING_SAIGE=0 to overwrite command files."
            continue
        fi
        mkdir -p "$dirout2"

        echo "#!/bin/bash
module load python/anaconda3/2020.7; source activate; conda activate R4.5.2; conda activate SAIGE
Rscript $dirsaige/step1_fitNULLGLMM.R --nThreads=$threads --outputPrefix=step1 \\
	--useSparseGRMtoFitNULL=TRUE --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE \\
	--sparseGRMFile=$dirtyp/typ.sparseGRM.mtx --sparseGRMSampleIDFile=$dirtyp/typ.sparseGRM.mtx.sampleIDs.txt \\
	--plink2File=$dirtyp/typ.prune --phenoFile=$phe_f $saige_str --phenoCol=$phe $saige_inv --covarColList=$covs --qCovarColList=sex
        " > "$dirout2/step1.cmd"
        if [[ "$SUBMIT_SAIGE" == "1" ]]; then
            cd "$dirout2"
            bsub -q spec -n "$threads" -J $phe.saige.step1 -o step1.LOG -e step1.ERR < step1.cmd
        fi

        for chr in $chr_list; do
            echo "#!/bin/bash
	module load python/anaconda3/2020.7; source activate; conda activate R4.5.2; conda activate SAIGE
	Rscript $dirsaige/step2_SPAtests.R --nThreads=$threads \\
		--bgen=$dirimp/chr$chr.bgen --bgenFileIndex $dirimp/chr$chr.bgen.bgi --sampleFile $dirimp/chr$chr.sample --chrom=$chr \\
		--is_imputed_data=TRUE --minMAC=$( [[ $type == bt ]] && echo 100 || echo 200 ) \\
		--GMMATmodelFile=step1.rda --varianceRatioFile=step1.varianceRatio.txt \\
		--SAIGEOutputFile=chr$chr.txt
            " > "$dirout2/chr$chr.cmd"
            if [[ "$SUBMIT_SAIGE" == "1" ]]; then
                cd "$dirout2"
                bsub -q short -n "$threads" -J $phe.saige.chr$chr -w "done($phe.saige.step1)" -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
            fi
        done
    done
}
run_step saige_gwas saige_gwas


# 🚩 LDSC GWAS
ldsc_gwas() {
    [[ -s "$script_f" ]] || { echo "ERROR: cannot find $script_f for phe_header_names" >&2; return 1; }
    source "$script_f"

    local dirref=$ldsc_ref
    local dirsof=$ldsc_soft
    local n_default=$ldsc_n_default
    local pattern=$ldsc_pattern
    local run_munge=$run_ldsc_munge
    local run_rg=$run_ldsc_rg
    mkdir -p "$ldsc_out/input"

    if [[ ! -d "$ldsc_in" ]]; then
        echo "ERROR: --ldsc-gwas-dir not found: $ldsc_in" >&2
        return 1
    fi

    find "$ldsc_in" -type f -name "$pattern" ! -name "*.001.gz" ! -name "chr*.gz" | sort -V > "$ldsc_out/gwas.files"
    if [[ ! -s "$ldsc_out/gwas.files" ]]; then
        echo "ERROR: no GWAS .gz files found under $ldsc_in with pattern $pattern" >&2
        echo "Use --ldsc-gwas-dir and/or --ldsc-pattern if needed." >&2
        return 1
    fi

    while read -r datf; do
        [[ -n "$datf" ]] || continue
        dat=$(basename "$datf" .gz)
        if [[ -f "$ldsc_out/$dat.sumstats.gz" ]]; then echo "$dat already munged"; continue; fi
        echo RUN LDSC munge: "$dat"
        phe_header_names "$datf"
        if [[ -z "${SNP:-}" || -z "${EA:-}" || -z "${NEA:-}" || -z "${BETA:-}" || -z "${P:-}" ]]; then
            echo "ERROR: missing required header in $datf. Need SNP, EA, NEA, BETA, P." >&2
            return 1
        fi
        if [[ -z "${N:-}" ]]; then n_str="--N $n_default"; else n_str="--N-col $N"; fi
        datf2="$datf"
        bad_p=$(zcat "$datf" | awk -v p="$P_col" 'NR>1 && $p<=0 {n++} END{print n+0}')
        if [[ "$bad_p" -gt 0 ]]; then
            datf2="$ldsc_out/input/$dat.fixp.gz"
            zcat "$datf" | awk -v p="$P_col" '{if(NR>1 && $p<=1e-300) $p=1e-300; print $0}' | sed 's/ /\t/g' | gzip -f > "$datf2"
        fi
        echo "#!/bin/bash
	module load python/anaconda3/2020.7
	cd $ldsc_out
	python2 $dirsof/munge_sumstats.py --chunksize 10000 --sumstats $datf2 --merge-alleles $dirref/hm3/w_hm3.snplist --out $ldsc_out/$dat --snp $SNP --a1 $EA --a2 $NEA $n_str --signed-sumstats $BETA,0 --p $P --ignore SNPID,OR
        " > "$ldsc_out/$dat.munge.cmd"
        chmod +x "$ldsc_out/$dat.munge.cmd"
    done < "$ldsc_out/gwas.files"

    if [[ "$run_munge" == "1" ]]; then
        cd "$ldsc_out"
        ls -1 *.munge.cmd 2>/dev/null | xargs -r -P "$ldsc_jobs" -I{} sh -c './{}'
    fi

    sed 's#.*/##; s/\.gz$//' "$ldsc_out/gwas.files" > "$ldsc_out/gwas.names"
    mapfile -t dats_ARR < "$ldsc_out/gwas.names"
    local ndat=${#dats_ARR[@]}
    if [[ "$ndat" -lt 2 ]]; then
        echo "LDSC rg requires at least 2 GWAS sumstats; found $ndat. Munge commands were still generated."
        return 0
    fi

    for ((i=0; i<ndat-1; i++)); do
        first="${dats_ARR[i]}"
        follower=("${dats_ARR[@]:i+1}")
        echo "#!/bin/bash
	cd $ldsc_out
	echo $first ${follower[*]} | sed -e 's/ /.sumstats.gz,/g' -e 's/\$/.sumstats.gz/' | xargs -I % \\
	python2 $dirsof/ldsc.py --rg % --out $ldsc_out/$first.rg --ref-ld-chr $dirref/hm3/ --w-ld-chr $dirref/hm3/
	beginl=\$(awk '\$1==\"Summary\" {printf NR}' $ldsc_out/$first.rg.log)
	if [[ -n \"\$beginl\" ]]; then
		awk -v s=\$beginl 'FNR > s' $ldsc_out/$first.rg.log | head -n -3 | sed 's/.sumstats.gz//g' > $ldsc_out/$first.rg.txt
	fi
        " > "$ldsc_out/${first}.ldsc.cmd"
        chmod +x "$ldsc_out/${first}.ldsc.cmd"
    done

    if [[ "$run_rg" == "1" ]]; then
        cd "$ldsc_out"
        ls -1 *.ldsc.cmd 2>/dev/null | xargs -r -P "$ldsc_jobs" -I{} sh -c './{}'
        if ls *.rg.txt >/dev/null 2>&1; then
            awk 'NR==1 || FNR>1' *.rg.txt | sed 's/  */ /g; s/^ //' > all.rg.res
        fi
        if ls *.rg.log >/dev/null 2>&1; then
            awk 'BEGIN { OFS="\t" } /^Heritability of phenotype 1$/ {cap=1; n=0; next} /^[-]+$/ && cap { next } cap {buf[++n] = $0; if (n == 6) {print FILENAME, buf[1], buf[2], buf[3], buf[4], buf[5], buf[6]; cap = 0}}' *.rg.log > all.h2.res
        fi
    fi
}
run_step ldsc_gwas ldsc_gwas

#!/usr/bin/env bash
set -euo pipefail
set +H

# Independent GWAS modules. Statistical/data conversion helpers live in f/;
# parameter parsing, genotype selection and analysis commands are defined here.
usage() {
    cat <<'EOF'
Usage:
  ./gwas.sh prep_gwas [options]
  ./gwas.sh run_plink2 [options]
  ./gwas.sh run_regenie [options]
  ./gwas.sh run_saige [options]

Modules:
  prep_gwas     Prepare shared pruned BED and optional sparse GRM once.
  run_plink2    Linear/logistic GWAS; quantitative, ordinal-score and binary traits.
  run_regenie   REGENIE model fitting and GWAS; also supports time-to-event.
  run_saige     SAIGE sparse null model and single-variant GWAS (qt/bt).

Options:
  --grch 37|38             Imputed genome build. Default: 37.
  --threads N             Threads per job. Default: 16.
  --jobs N                Concurrent chromosome jobs per trait. Default: 1.
  --memory N              PLINK memory in MiB. Default: 24000.
  --typ-dir DIR           Shared prep directory. Default: /mnt/h/ukbGen/<grch>/typ.
  --imputed-dir DIR       Default: /mnt/h/ukbGen/<grch>/imp.
  --pheno FILE            Default: /mnt/d/data/ukb/phe/common/ukb.phe.
  --pheno-name CSV        Default: height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e,cvd_cad.adu.
  --type TYPE             auto, qt, ordinal, bt, t2e. Default: auto.
  --event-col NAME        T2E event column. Default: <trait stem>.Yt2e.
  --covariate-file FILE   Default: the --pheno file.
  --covariates CSV        Default: age,sex,tdi,PC1,PC2,PC3,PC4.
  --categorical-covariates CSV  Default: center. Use "" for none.
  --output-dir DIR        Default: /mnt/d/data.BIG/gwas/self/common.
  --chr LIST              Comma/space-separated, e.g. 1,X. Default: 1-22,X
                          (genome-wide); prep_gwas defaults to 1-22.
  --extract FILE          Variant IDs for association subset/testing.
  --keep FILE             FID IID subset for PLINK2/REGENIE; not shared prep.
  --min-mac N             Association minimum MAC. Default: 100.
  --sparse-grm TRUE|FALSE Create sparse GRM in prep_gwas. Default: TRUE.
  --run-cmd TRUE|FALSE    Write run.cmd; TRUE also executes it. Default: TRUE.
  --replace TRUE|FALSE    Re-run completed tasks. Default: FALSE.
  --plink2 PATH           PLINK2 executable (PATH or /mnt/d/software/bin).
  --plink PATH            PLINK1 executable for merging pruned BEDs.
  --regenie PATH          REGENIE executable (version >=4 for T2E).
  --rscript PATH          Rscript for phenotype preparation.
  --saige-rscript PATH    Default: isolated f/saige_Rscript.sh when installed.
  --saige-dir DIR         Directory with SAIGE R scripts; auto-detected.
  -h, --help              Show this help.

Aliases: --phenotype-file = --pheno; --phenotypes/--traits = --pheno-name.

Examples:
  # Generate phenotypes using common/adu_cip.tsv (currently provisional UKB CIP):
  ADULT_ENABLE=TRUE ADULT_TRAITS=cvd_cad ADULT_CIP_FILE=/mnt/d/data/ukb/phe/common/adu_cip.tsv \
    bash /mnt/d/scripts/ukb/phe.sh --steps phe4gwas

  # Prepare once; successful prune tasks resume when adding the GRM later:
  ./gwas.sh prep_gwas --sparse-grm TRUE --threads 16

  # Each trait has its own samples and Step 1 model; Step 2 runs by chromosome.
  # No --chr: genome-wide. At most jobs * threads worker threads are used.
  ./gwas.sh run_regenie --grch 38 --jobs 4 --threads 4 --run-cmd TRUE \
    --pheno /mnt/d/data/ukb/phe/common/ukb.phe \
    --pheno-name height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e,cvd_cad.adu

  # Run only selected chromosomes (Step 1 still uses the shared pruned panel):
  ./gwas.sh run_regenie --chr 1,X --jobs 2 --threads 4 --run-cmd TRUE

  ./gwas.sh run_regenie --grch 37 --threads 16 --pheno /mnt/d/data/ukb/phe/common/ukb.phe --pheno-name cvd_cad.adu --type qt

  ./gwas.sh run_plink2 --pheno-name height,bald,bald12,cvd_cad.Yt2e
  ./gwas.sh run_saige --pheno-name height,bald,bald12,cvd_cad.Yt2e

  # Inspect commands without starting GWAS; separate output for another build:
  ./gwas.sh run_regenie --grch 38 --run-cmd FALSE \
    --output-dir /mnt/d/data.BIG/gwas/self/grch38/common

  # QC comparison after GWAS:
  /mnt/d/scripts/gwas/gwas_compare.sh compare \
    --gwas-files /mnt/d/data.BIG/gwas/main/common/bald/gwas/bald.gz,/mnt/d/data.BIG/gwas/self/common/bald/gwas/bald.regenie.gz \
    --mplot TRUE --compare-beta TRUE --compare-EAF TRUE


* Ordinal bald uses a quantitative score without rank normalization. 
* Binary input uses 0/1; cvd_cad.t2e pairs with cvd_cad.Yt2e. PLINK2/SAIGE skip T2E explicitly.
* All *.adu traits use --qt in auto mode, with a separate REGENIE model fit.
* ADuLT is standardized once; no --apply-rint, --bt, --t2e or event column is used.
* ADuLT BETA is on its transformed quantitative scale, not log-OR or log-HR.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }
bool_word() {
    case "${1^^}" in TRUE) echo TRUE ;; FALSE) echo FALSE ;; *) die 'Use TRUE or FALSE' ;; esac
}
tool_path() {
    if command -v "$1" >/dev/null 2>&1; then command -v "$1"
    elif [[ -x /mnt/d/software/bin/$1 ]]; then echo "/mnt/d/software/bin/$1"
    else echo "$1"; fi
}

# 🚩 Defaults and argument parsing
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
helper="$here/f/gwas_runner.py"
module=${1:-}
case "$module" in
    -h|--help|'') usage; exit 0 ;;
    prep_gwas|run_plink2|run_regenie|run_saige) shift ;;
    *) die "Unknown module: $module (see --help)" ;;
esac
grch=37 threads=16 memory=24000 min_mac=100 jobs=1
typed_dir=''
imputed_dir='' phenotype_file=/mnt/d/data/ukb/phe/common/ukb.phe covariate_file=''
phenotypes=height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e,cvd_cad.adu
type=auto event_col='' covariates=age,sex,tdi,PC1,PC2,PC3,PC4 categorical_covariates=center
output_dir=/mnt/d/data.BIG/gwas/self/common
chromosomes=''
extract='' keep='' sparse_grm=TRUE run=TRUE replace=FALSE
plink2=$(tool_path plink2); plink=$(tool_path plink); regenie=$(tool_path regenie); rscript=$(tool_path Rscript)
saige_dir=/mnt/d/software/SAIGE/extdata; saige_rscript=$rscript
for saige_env in "$HOME/anaconda3/envs/saige" "$HOME/miniforge3/envs/saige"; do
    if [[ -f "$saige_env/bin/createSparseGRM.R" ]]; then
        saige_dir="$saige_env/bin"; saige_rscript="$here/f/saige_Rscript.sh"; break
    fi
done
while (( $# )); do
    case "$1" in -h|--help) usage; exit 0 ;; esac
    (( $# >= 2 )) || die "$1 requires a value"
    # run/replace are read through indirect expansion in save_jobs.
    # shellcheck disable=SC2034
    case "$1" in
        --grch) grch=$2 ;; --threads) threads=$2 ;; --memory) memory=$2 ;;
        --typ-dir) typed_dir=$2 ;; --imputed-dir) imputed_dir=$2 ;;
        --jobs) jobs=$2 ;;
        --pheno|--phenotype-file) phenotype_file=$2 ;;
        --pheno-name|--phenotypes|--traits) phenotypes=$2 ;;
        --covariate-file) covariate_file=$2 ;; --type) type=$2 ;; --event-col) event_col=$2 ;;
        --covariates) covariates=$2 ;; --categorical-covariates) categorical_covariates=$2 ;;
        --output-dir) output_dir=$2 ;; --chr) chromosomes=$2 ;;
        --extract) extract=$2 ;; --keep) keep=$2 ;; --min-mac) min_mac=$2 ;;
        --plink2) plink2=$2 ;; --plink) plink=$2 ;; --regenie) regenie=$2 ;; --rscript) rscript=$2 ;;
        --saige-rscript) saige_rscript=$2 ;; --saige-dir) saige_dir=$2 ;;
        --sparse-grm) sparse_grm=$(bool_word "$2") ;;
        --run-cmd) run=$(bool_word "$2") ;;
        --replace) replace=$(bool_word "$2") ;;
        *) die "Unknown option: $1" ;;
    esac
    shift 2
done
[[ "$grch" == 37 || "$grch" == 38 ]] || die '--grch must be 37 or 38'
case "$type" in auto|qt|ordinal|bt|t2e) ;; *) die 'Invalid --type' ;; esac
for number in threads memory min_mac jobs; do
    [[ ${!number} =~ ^[0-9]+$ ]] || die "Invalid --${number//_/-}"
    printf -v "$number" '%d' "$((10#${!number}))"
done
(( threads > 0 && jobs > 0 && memory >= 640 && min_mac > 0 )) || die 'Invalid threads/jobs/memory/min-mac'
if [[ -z "$chromosomes" ]]; then
    chromosomes=$(seq -s, 1 22)
    [[ "$module" == prep_gwas ]] || chromosomes+=,X
fi
read -r -a chrs <<< "${chromosomes//,/ }"
(( ${#chrs[@]} )) || die 'Empty chromosome list'
declare -A seen_chr=()
for chr in "${chrs[@]}"; do
    [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2]|X)$ ]] || die "Invalid chromosome: $chr"
    [[ ! ${seen_chr[$chr]+yes} ]] || die "Duplicate chromosome: $chr"
    [[ "$module" != prep_gwas || "$chr" != X ]] || die 'Prep uses autosomes only'
    seen_chr[$chr]=1
done
mapfile -t chrs < <(printf '%s\n' "${chrs[@]}" | sort -V)
chromosomes=$(IFS=,; echo "${chrs[*]}")
typed_dir=${typed_dir:-/mnt/h/ukbGen/$grch/typ}
imputed_dir=${imputed_dir:-/mnt/h/ukbGen/$grch/imp}
covariate_file=${covariate_file:-$phenotype_file}
[[ -n "$phenotypes" ]] || die 'Empty phenotype list'
[[ "$phenotypes" != ,* && "$phenotypes" != *, && "$phenotypes" != *,,* ]] || die 'Empty trait name in --pheno-name'

# 🚩 Command-plan utilities. The Python helper only handles JSON/checkpoints,
# logged execution, and PGEN file views; it does not define GWAS commands.
jobs_file=$(mktemp)
trap 'rm -f -- "$jobs_file"' EXIT
cmd=() inputs=() outputs=() pfile_args=() geno_inputs=()
add_job() {
    python3 "$helper" plan-node "$1" "${2:-}" "${3:-}" \
        --inputs "${inputs[@]}" --outputs "${outputs[@]}" -- "${cmd[@]}" >> "$jobs_file"
}
save_jobs() {
    local key config=()
    for key in module grch threads memory jobs typed_dir imputed_dir phenotype_file covariate_file \
        phenotypes type event_col covariates categorical_covariates output_dir chromosomes \
        extract keep min_mac plink2 plink regenie rscript saige_rscript saige_dir sparse_grm run replace; do
        config+=("$key=${!key}")
    done
    python3 "$helper" finalize-plan "$1" "$jobs_file" "${config[@]}"
}
genotype() {
    pfile_args=(--pfile "$1")
    if [[ -f $1.pvar.zst ]]; then
        pfile_args+=(vzs); geno_inputs=("$1.pgen" "$1.pvar.zst" "$1.psam")
    else geno_inputs=("$1.pgen" "$1.pvar" "$1.psam"); fi
}

# 🚩 Shared typed-genotype pruning and sparse GRM
prep_gwas() {
    local root cache chr src dst merged prefix actual
    local beds=() merge_inputs=() base=("$plink2" --threads "$threads" --memory "$memory")
    root=$(realpath -m -- "$typed_dir"); cache="$root/.prep"; mkdir -p "$cache"
    : > "$jobs_file"
    printf '6\t25000000\t35000000\tMHC\n' > "$cache/mhc.range"
    for chr in "${chrs[@]}"; do
        src="$root/chr$chr"; dst="$cache/chr$chr.prune"; genotype "$src"
        # k=0 explicitly preserves the fixed HWE cutoff for the shared GRM panel.
        cmd=("${base[@]}" "${pfile_args[@]}" --maf .01 --mac 100 --max-alleles 2 \
            --geno .1 --hwe 1e-15 0 --indep-pairwise 1000kb 0.1 --out "$dst")
        [[ "$chr" != 6 ]] || cmd+=(--exclude range "$cache/mhc.range")
        inputs=("${geno_inputs[@]}"); outputs=("$dst.prune.in"); add_job "chr$chr.prune"
        cmd=("${base[@]}" "${pfile_args[@]}" --extract "$dst.prune.in" --make-bed --out "$dst")
        inputs=("${geno_inputs[@]}" "$dst.prune.in"); outputs=("$dst.bed" "$dst.bim" "$dst.fam")
        add_job "chr$chr.bed"; beds+=("$dst"); merge_inputs+=("${outputs[@]}")
    done
    merged="$root/typ.prune"
    printf '%s\n' "${beds[@]:1}" > "$cache/merge.list"
    cmd=("$plink" --bfile "${beds[0]}" --make-bed --out "$merged" --threads "$threads" --memory "$memory")
    (( ${#beds[@]} <= 1 )) || cmd+=(--merge-list "$cache/merge.list")
    inputs=("${merge_inputs[@]}"); outputs=("$merged.bed" "$merged.bim" "$merged.fam"); add_job merge
    if [[ "$sparse_grm" == TRUE ]]; then
        prefix="$root/typ"; actual="${prefix}_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
        cmd=("$saige_rscript" "$saige_dir/createSparseGRM.R" "--plinkFile=$merged" "--nThreads=$threads" \
            --numRandomMarkerforSparseKin=2000 --relatednessCutoff=0.125 "--outputPrefix=$prefix")
        cmd=(python3 "$helper" create-grm "$actual" "$prefix" "${cmd[@]}")
        inputs=("$merged.bed" "$merged.bim" "$merged.fam"); outputs=("$prefix.sparseGRM.mtx" "$prefix.sparseGRM.mtx.sampleIDs.txt"); add_job sparseGRM
    fi
    save_jobs "$cache"
}

# 🚩 Shared phenotype setup and output collection
setup_trait() {
    skip_trait=FALSE
    trait=$1; trait_type=$type
    [[ -n "$trait" && "$trait" != . && "$trait" != .. && "$trait" != *[/\\]* && "$trait" != *$'\n'* && "$trait" != *$'\r'* ]] || die 'Invalid trait name'
    if [[ "$trait_type" == auto ]]; then
        case "$trait" in
            height|*.adu) trait_type=qt ;; bald) trait_type=ordinal ;;
            bald12|cvd_cad.Yt2e) trait_type=bt ;; cvd_cad.t2e) trait_type=t2e ;;
            *) die "Specify --type for $trait" ;;
        esac
    fi
    root="$(realpath -m -- "$output_dir")/$trait/gwas"; work="$root/$method"; mkdir -p "$work"
    event=${event_col:-${trait%.t2e}.Yt2e}
    phe="$work/analysis.phe"; covars="$phe.covars"; pruned="$typed_dir/typ.prune"; grm="$typed_dir/typ.sparseGRM.mtx"
    keep_args=(); extract_args=(); raw=(); : > "$jobs_file"
    [[ -z "$keep" ]] || keep_args=(--keep "$keep")
    [[ -z "$extract" ]] || extract_args=(--extract "$extract")
    if [[ "$trait_type" == t2e && "$method" != regenie ]]; then
        printf '%s\n' "$method: time-to-event is unsupported; use run_regenie --type t2e" | tee "$work/UNSUPPORTED.txt"
        skip_trait=TRUE; return 0
    fi
    cmd=("$rscript" "$here/f/gwas_inputs.R" "$phenotype_file" "$covariate_file" "$trait" "$trait_type" "$event" "$covariates" "$categorical_covariates" "$phe")
    inputs=("$phenotype_file" "$covariate_file" "$here/f/gwas_inputs.R"); outputs=("$phe" "$covars"); add_job phenotype
    common_inputs=("$phe" "$covars"); [[ -z "$keep" ]] || common_inputs+=("$keep")
}
association_inputs() {
    genotype "$prefix"; inputs=("${geno_inputs[@]}" "${common_inputs[@]}")
    [[ -z "$extract" ]] || inputs+=("$extract")
}
finish_trait() {
    local final="$root/$trait.$method.gz"
    cmd=(python3 "$here/f/gwas_io.py" "$method" "$final" "${raw[@]}")
    inputs=("${raw[@]}" "$here/f/gwas_io.py"); outputs=("$final"); add_job merge-results
    cmd=(python3 "$helper" build-metadata "$final" "$grch")
    inputs=("$final"); outputs=("$final.grch"); add_job build-metadata
    save_jobs "$work"
}

# 🚩 PLINK2 linear / logistic association
run_plink2() {
    method=plink2
    for trait in "${trait_list[@]}"; do
        setup_trait "$trait"
        [[ "$skip_trait" != TRUE ]] || continue
        for chr in "${chrs[@]}"; do
            prefix="$imputed_dir/chr$chr"; dest="$work/chr$chr"; association_inputs
            glm=(hide-covar allow-no-covars omit-ref no-x-sex 'cols=+omitted,+a1freq,+nobs')
            result="$dest.$trait.glm.linear"
            if [[ "$trait_type" == bt ]]; then glm+=(firth-fallback); result="$dest.$trait.glm.logistic.hybrid"; fi
            cmd=("$plink2" --threads "$threads" --memory "$memory" "${pfile_args[@]}" --glm "${glm[@]}" \
                --pheno "$phe" --pheno-name "$trait" --no-psam-pheno --mac "$min_mac" --max-alleles 2 --out "$dest" "${keep_args[@]}" "${extract_args[@]}")
            [[ "$trait_type" != bt ]] || cmd+=(--1)
            [[ "$chr" != X ]] || cmd+=(--xchr-model 2)
            outputs=("$result"); add_job "chr$chr.association" "$covars" "$method"; raw+=("$result")
        done
        finish_trait
    done
}

# 🚩 REGENIE model fitting and per-chromosome association
run_regenie() {
    method=regenie
    for trait in "${trait_list[@]}"; do
        setup_trait "$trait"
        [[ "$skip_trait" != TRUE ]] || continue
        flags=(--qt); phflags=(--phenoCol "$trait")
        case "$trait_type" in
            bt) flags=(--bt) ;;
            t2e) flags=(--t2e); phflags=(--phenoColList "$trait" --eventColList "$event") ;;
        esac
        common=("$regenie" --threads "$threads" --bsize 1000 --phenoFile "$phe" "${phflags[@]}" "${flags[@]}" "${keep_args[@]}")
        step1="$work/step1"
        cmd=("${common[@]}" --step 1 --bed "$pruned" --lowmem --lowmem-prefix "$work/tmp" --out "$step1")
        inputs=("${common_inputs[@]}" "$pruned.bed" "$pruned.bim" "$pruned.fam")
        outputs=("${step1}_pred.list" "${step1}_1.loco"); add_job null-model "$covars" "$method"
        for chr in "${chrs[@]}"; do
            prefix="$imputed_dir/chr$chr"; dest="$work/chr$chr"; view="$work/chr$chr.input"; genotype "$prefix"
            # Plain PVAR view for REGENIE; large PGEN stays in place via a symlink.
            cmd=(python3 "$helper" pgen-view "$prefix" "$view")
            inputs=("${geno_inputs[@]}"); outputs=("$view.pgen" "$view.pvar" "$view.psam"); add_job "chr$chr.pgen-view"
            association_inputs
            inputs+=("${step1}_pred.list" "${step1}_1.loco" "$view.pgen" "$view.pvar" "$view.psam")
            cmd=("${common[@]}" --step 2 --pgen "$view" --pred "${step1}_pred.list" --minMAC "$min_mac" --out "$dest" "${extract_args[@]}")
            [[ "$trait_type" != bt ]] || cmd+=(--firth --approx --pThresh 0.01)
            [[ "$chr" != X ]] || cmd+=(--par-region "b$grch")
            result="${dest}_$trait.regenie"; outputs=("$result")
            add_job "chr$chr.association" "$covars" "$method"; raw+=("$result")
        done
        finish_trait
    done
}

# 🚩 SAIGE sparse null model and dosage-based association
run_saige() {
    method=saige
    [[ -z "$keep" ]] || die 'SAIGE --keep: provide a subset --pheno file instead'
    for trait in "${trait_list[@]}"; do
        setup_trait "$trait"
        [[ "$skip_trait" != TRUE ]] || continue
        step1="$work/step1"; saige_type=quantitative
        [[ "$trait_type" != bt ]] || saige_type=binary
        cmd=("$saige_rscript" "$saige_dir/step1_fitNULLGLMM.R" "--plinkFile=$pruned" \
            "--phenoFile=$phe" "--phenoCol=$trait" --sampleIDColinphenoFile=IID "--traitType=$saige_type" \
            --invNormalize=FALSE --useSparseGRMtoFitNULL=TRUE --skipVarianceRatioEstimation=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE --LOCO=FALSE "--nThreads=$threads" "--outputPrefix=$step1" \
            "--sparseGRMFile=$grm" "--sparseGRMSampleIDFile=$grm.sampleIDs.txt")
        inputs=("${common_inputs[@]}" "$pruned.bed" "$pruned.bim" "$pruned.fam" "$grm" "$grm.sampleIDs.txt")
        outputs=("$step1.rda" "$step1.varianceRatio.txt"); add_job null-model "$covars" "$method"
        for chr in "${chrs[@]}"; do
            prefix="$imputed_dir/chr$chr"; dest="$work/chr$chr"; vcf="$dest.dosage.vcf.gz"; genotype "$prefix"
            # Force DS output, retaining dosage for SAIGE instead of hard-call BED.
            cmd=("$plink2" --threads "$threads" --memory "$memory" "${pfile_args[@]}" --mac "$min_mac" \
                --max-alleles 2 --export vcf bgz vcf-dosage=DS-force id-paste=iid --out "$dest.dosage" "${extract_args[@]}")
            inputs=("${geno_inputs[@]}"); [[ -z "$extract" ]] || inputs+=("$extract")
            outputs=("$vcf"); add_job "chr$chr.export"
            cmd=("$(tool_path tabix)" -f -C -p vcf "$vcf")
            inputs=("$vcf"); outputs=("$vcf.csi"); add_job "chr$chr.index"
            result="$dest.saige.txt"; association_inputs
            inputs+=("$vcf" "$vcf.csi" "$step1.rda" "$step1.varianceRatio.txt")
            cmd=("$saige_rscript" "$saige_dir/step2_SPAtests.R" "--vcfFile=$vcf" "--vcfFileIndex=$vcf.csi" --vcfField=DS "--chrom=$chr" \
                "--GMMATmodelFile=$step1.rda" "--varianceRatioFile=$step1.varianceRatio.txt" \
                "--SAIGEOutputFile=$result" --LOCO=FALSE "--minMAC=$min_mac" --is_output_moreDetails=TRUE --is_overwrite_output=TRUE)
            [[ "$trait_type" != bt ]] || cmd+=(--is_Firth_beta=TRUE --pCutoffforFirth=0.01)
            outputs=("$result"); add_job "chr$chr.association"; raw+=("$result")
        done
        finish_trait
    done
}

# 🚩 Dispatch one independent module, never a sequence of GWAS methods.
IFS=, read -r -a trait_list <<< "$phenotypes"
if [[ "$module" != prep_gwas && "$run" == TRUE ]]; then
    # Check every requested trait before starting any expensive model fit.
    python3 - "$phenotype_file" "$phenotypes" <<'PY'
import gzip
import sys
path, requested = sys.argv[1:]
opener = gzip.open if path.endswith('.gz') else open
with opener(path, 'rt') as handle:
    columns = handle.readline().split()
missing = [name for name in requested.split(',') if name not in columns]
if missing:
    sys.exit('ERROR: Missing phenotype columns in ' + path + ': ' + ', '.join(missing)
             + '. Generate the requested phenotypes before running GWAS.')
PY
fi
"$module"

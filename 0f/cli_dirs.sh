#!/usr/bin/env bash
# Translate public directory options at shell boundaries; retain legacy backend argv.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    case "${1:-}" in
        -h|--help)
            cat <<'HELP'
Modules:
    Library  Normalize --dir-XXX options without changing backend parameter names

Main options:
    Legacy directory option names remain accepted for saved commands.
    -h, --help  Show help; conventions: ../0.README.md

Example:
    source ./cli_dirs.sh
HELP
            exit 0;;
    esac
fi

script_directory_args() {
    local script=${1:?script path required}
    shift
    local group pair arg option value stopped=0
    local -A aliases=()
    case "$script" in
        /mnt/d/scripts/0f/*) group="--dir-pfile=--pfile-dir --dir-tmp=--tmp-dir --dir-work=--work-dir";;
        /mnt/d/scripts/allele/*|/mnt/d/scripts/allele3/*) group="--dir-atlas=--atlas-dir --dir-cache=--cache-dir --dir-out=--out --dir-root=--root --dir-vcf=--vcf-dir";;
        /mnt/d/scripts/bmonth/*) group="--dir-in=--in --dir-out=--out --dir-temp=--temp-dir";;
        /mnt/d/scripts/ems120/*) group="--dir-analysis=--analysis-dir";;
        /mnt/d/scripts/grid/*) group="--dir-analysis-root=--analysis-root --dir-arg=--arg-dir --dir-arg-hap=--arg-hap-dir --dir-arg-hap-root=--arg-hap-root --dir-arg-map=--arg-map-dir --dir-arg-trees=--arg-trees-dir --dir-cache=--cache-dir --dir-csx-ref=--csx-ref-dir --dir-data-root=--data-root --dir-frequency=--frequency-dir --dir-ldscore=--ldscore-dir --dir-out=--out-dir --dir-out-root=--out-root --dir-output=--output-dir --dir-output-root=--output-root --dir-pca=--pca-dir --dir-posterior-frequency=--posterior-frequency-dir --dir-prsformer-cache=--prsformer-cache-dir --dir-prsformer-root=--prsformer-root --dir-prsformer-run=--prsformer-run-dir --dir-raw=--raw-dir --dir-ref=--ref-dir --dir-roster=--roster-dir --dir-run=--run-dir --dir-score=--score-dir --dir-sumstats=--sumstats-dir --dir-target=--target-dir --dir-ukb-hap-root=--ukb-hap-root --dir-upstream=--upstream-dir";;
        /mnt/d/scripts/gu/*) group="--dir-analysis-root=--analysis-root --dir-archaic-root=--archaic-root --dir-arg=--arg-dir --dir-arg-hap=--arg-hap-dir --dir-arg-map=--arg-map-dir --dir-data-root=--data-root --dir-gen4arg=--gen4arg-dir --dir-hap=--hap-dir --dir-ibdmix=--ibdmix-dir --dir-map=--map-dir --dir-model=--model-dir --dir-output=--output-dir --dir-raw=--raw-dir --dir-reference-panel=--reference-panel-dir --dir-reference-root=--reference-root --dir-root=--root --dir-software=--software-dir --dir-target-root=--target-root --dir-target-vcf=--target-vcf-dir --dir-tmp-root=--tmp-root --dir-trace=--trace-dir --dir-ukb-1kg-vcf=--ukb-1kg-vcf-dir --dir-ukb-archaic-root=--ukb-archaic-root --dir-ukb-hap-root=--ukb-hap-root --dir-ukb-mask-root=--ukb-mask-root --dir-ukb-pgen-root=--ukb-pgen-root --dir-ukb-results-root=--ukb-results-root --dir-ukb-run=--ukb-run-dir --dir-ukb-typed-root=--ukb-typed-root";;
        /mnt/d/scripts/gwas/*) group="--dir-block=--block-dir --dir-chain=--chain-dir --dir-data-root=--data-root --dir-gene=--gene-dir --dir-imputed=--imputed-dir --dir-ld=--ld-dir --dir-ldsc-software=--ldsc-software-dir --dir-output=--output-dir --dir-pfile=--pfile-dir --dir-pgs-pfile=--pgs-pfile-dir --dir-project=--project-dir --dir-reference=--reference-dir --dir-refgen-id=--refgen-id-dir --dir-root=--root --dir-saige=--saige-dir --dir-typ=--typ-dir --dir-work=--work-dir";;
        /mnt/d/scripts/le8/*) group="--dir-ukb-phe=--ukb-phe --dir-analysis-root=--analysis-root --dir-cache-root=--cache-root --dir-data-root=--data-root --dir-gwas=--gwas-dir --dir-match-tmp=--match-tmp-dir --dir-mqtl=--mqtl-dir --dir-mrlink2-ref-id=--mrlink2-ref-id-dir --dir-mrlink2-ref-pfile=--mrlink2-ref-pfile-dir --dir-out=--outdir --dir-pqtl=--pqtl-dir --dir-reference-id=--reference-id-dir --dir-reference-pfile=--reference-pfile-dir --dir-refgen-root=--refgen-root --dir-tmp=--tmp-dir";;
        /mnt/d/scripts/maha/*) group="--dir-chns-project=--chns-project-dir --dir-helper=--helper-dir --dir-output=--output-dir";;
        /mnt/d/scripts/metag/*) group="--dir-cache=--cache-dir --dir-data-root=--data-root --dir-db=--db --dir-out=--out --dir-reference-root=--reference-root --dir-root=--root";;
        /mnt/d/scripts/nextstrain/*) group="--dir-data-root=--data-root --dir-gisaid=--gisaid-dir --dir-input=--input-dir --dir-out=--outdir --dir-output=--output-dir";;
        /mnt/d/scripts/pfi/*) group="--dir-analysis-root=--analysis-root --dir-in=--in_dir --dir-model=--model_dir --dir-out=--out_dir --dir-raw-root=--raw-root --dir-summary=--summary_dir --dir-temp=--temp_dir --dir-work=--work_dir";;
        /mnt/d/scripts/phole/*) group="--dir-data-root=--data-root --dir-gwas-root=--gwas-root --dir-matrix=--matrix-dir --dir-model=--model-dir --dir-out=--out --dir-validation=--validation-dir";;
        /mnt/d/scripts/ukb/*) group="--dir-data-root=--data-root --dir-genotype=--genotype-dir --dir-output=--output-dir --dir-pfile=--pfile-dir --dir-pgs-pfile=--pgs-pfile-dir --dir-phenotype=--phenotype-dir --dir-script=--script-dir --dir-work=--work-dir";;
        *) group="--dir-root=--root";;
    esac
    for pair in $group; do aliases[${pair%%=*}]=${pair#*=}; done
    SCRIPT_DIRECTORY_ARGS=()
    for arg in "$@"; do
        if ((stopped)); then
            SCRIPT_DIRECTORY_ARGS+=("$arg")
        elif [[ "$arg" == -- ]]; then
            stopped=1
            SCRIPT_DIRECTORY_ARGS+=("$arg")
        else
            option=${arg%%=*}
            if [[ -n "$option" && -n ${aliases[$option]:-} ]]; then
                SCRIPT_DIRECTORY_ARGS+=("${aliases[$option]}")
                if [[ "$arg" == *=* ]]; then SCRIPT_DIRECTORY_ARGS+=("${arg#*=}"); fi
            else
                SCRIPT_DIRECTORY_ARGS+=("$arg")
            fi
        fi
    done
}

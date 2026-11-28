#!/usr/bin/env bash
# Normalize named permanent-output options before existing planners parse argv.
# No directories are created here; existing defaults and execution remain unchanged.
script_output_directory_args() {
    local script=$1; shift
    local key="$(basename -- "$(dirname -- "$script")")/$(basename -- "$script")"
    local arg option target value name tail
    local -A aliases=()
    case "$key" in
        allele3/0.ancient_split.sh)
            aliases[--dir-out-ancient37]='--out37'
            aliases[--dir-out-ancient38]='--out38'
            ;;
        allele3/1.ancient_pca.sh)
            aliases[--dir-out-ancient_pca]='--out'
            ;;
        grid/1.pca.sh)
            aliases[--dir-out-pca_projection]='@pca'
            ;;
        ygb/gen.sh)
            aliases[--dir-out-pfile]='--dir-pfile'
            aliases[--dir-out-lift_bed]='--dir-out'
            aliases[--dir-out-liftover]='--dir-out'
            aliases[--dir-out-pgs]='env:YGB_PGS_OUTPUT_DIR'
            aliases[--dir-out-snp]='env:YGB_SNP_OUTPUT_DIR'
            ;;
        ygb/phe.sh)
            aliases[--dir-out-phenotype]='--phenotype-dir'
            ;;
        gwas/format.sh)
            aliases[--dir-out-gwas_format]='--dir-out'
            ;;
        gu/1.arg.sh)
            aliases[--dir-out-arg]='--dir-arg'
            aliases[--dir-out-gen4arg]='--dir-gen4arg'
            ;;
        gu/gu.paper.sh)
            aliases[--dir-out-ibdmix_report]='--dir-out'
            ;;
        gu/gu.sh)
            aliases[--dir-out-ygb_vcf]='--ygb-vcf-out'
            ;;
        flt/flt.sh)
            aliases[--dir-out-phewas_cache]='env:FLT_PHEWAS_CACHE'
            ;;
        metag/1.lib.sh)
            aliases[--dir-out-metag_library]='env:METAG_DIR_OUT_LIBRARY_DIR'
            aliases[--dir-out-metag_metadata]='env:METAG_DIR_OUT_METADATA_DIR'
            aliases[--dir-out-metag_reference]='env:METAG_DIR_OUT_REFERENCE_RAW'
            aliases[--dir-out-human_reference]='env:METAG_DIR_OUT_HUMAN_REFERENCE_DIR'
            ;;
        metag/metag.sh)
            aliases[--dir-out-metag_species]='env:METAG_DIR_OUT_OUTPUT_DIR'
            aliases[--dir-out-metag_reads]='env:METAG_DIR_OUT_READS_RAW'
            ;;
        metag/metag.paper.sh)
            aliases[--dir-out-metag_species]='env:METAG_DIR_OUT_OUTPUT_DIR'
            ;;
        pfi/pfi.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            aliases[--dir-out-journal_rank]='file:PFI_JOURNAL_RANK_REAL:journal_rank_real.csv'
            ;;
        pfi/00_run_pfi_all.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            aliases[--dir-out-journal_rank]='file:PFI_JOURNAL_RANK_REAL:journal_rank_real.csv'
            ;;
        pfi/02_download_pubmed.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            ;;
        pfi/03_parse_pubmed.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            ;;
        pfi/04_download_pmc_oa_bulk_xml.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            ;;
        pfi/05_parse_pmc_oa_bulk_xml.sh)
            aliases[--dir-out-pfi_raw]='env:PFI_RAW_ROOT'
            ;;
        pfi/00_run_pfi_after08_auto.sh)
            aliases[--dir-out-journal_rank]='file:PFI_JOURNAL_RANK_REAL:journal_rank_real.csv'
            ;;
        pfi/08b_normalize_journal_rank.sh)
            aliases[--dir-out-journal_rank]='file:PFI_JOURNAL_RANK_REAL:journal_rank_real.csv'
            ;;
    esac
    SCRIPT_OUTPUT_DIRECTORY_ARGS=()
    while (($#)); do
        arg=$1; shift
        if [[ $arg == -- ]]; then
            SCRIPT_OUTPUT_DIRECTORY_ARGS+=("$arg" "$@"); break
        fi
        option=${arg%%=*}
        target=
        if [[ -n $option ]]; then target=${aliases[$option]:-}; fi
        if [[ -z $target ]]; then
            SCRIPT_OUTPUT_DIRECTORY_ARGS+=("$arg"); continue
        fi
        if [[ $arg == *=* ]]; then value=${arg#*=}; else
            if (($# == 0)) || [[ $1 == --* ]]; then
                printf 'ERROR: %s requires DIR\n' "$option" >&2; return 2
            fi
            value=$1; shift
        fi
        [[ -n $value ]] || { printf 'ERROR: %s requires nonempty DIR\n' "$option" >&2; return 2; }
        case "$target" in
            env:*)
                [[ $value == /* ]] || value="$PWD/$value"
                export "${target#env:}=$value"
                if [[ $target == env:METAG_DIR_OUT_HUMAN_REFERENCE_DIR ]]; then export HUMAN_DIR="$value"; fi ;;
            file:*)
                [[ $value == /* ]] || value="$PWD/$value"
                tail=${target#file:}; name=${tail%%:*}
                export "$name=${value%/}/${tail#*:}" ;;
            @pca)
                SCRIPT_OUTPUT_DIRECTORY_ARGS+=(--pca-file "${value%/}/ygb.discodivas.pca.tsv.gz" --ancestry-file "${value%/}/ygb.ancestry.auto.tsv.gz") ;;
            *) SCRIPT_OUTPUT_DIRECTORY_ARGS+=("$target" "$value") ;;
        esac
    done
}

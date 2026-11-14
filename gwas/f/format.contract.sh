#!/usr/bin/env bash
# Variables are passed explicitly; no eval or command-text serialization.
gwas_contract() {
    local action=$1 stage=$2 key
    local -a settings=()
    for key in GRCH COL_MAP CHR_PREFIX HM3_MODE P_HM3 RSID_MODE RSID_UNMATCHED FILL_N FILL_EAF N_TOTAL DO_LIFTOVER MAGMA_N GWAS_N MAGMA_WINDOW CHRS P_LEAD LEAD_WINDOW CLUMP_KB REFGEN_POP RAW PHEF PERF_F DBSNP HM3 HM3_POS REFGEN_CLUMP REFGEN_COJO REFGEN_KEEP CHAIN FINAL MAGMA_REF GENE_LOC SYNONYMS MAGMA_PREFIX MAGMA_DIR MERGED GWAS_DIR GWAS GWAS_POST_MAGMA_ROWS; do
        settings+=("$key=${!key-}")
    done
    env "${settings[@]}" python3 "${PERF_F%/*}/format.contract.py" "$action" "$stage" "${QC_PREFIX}.${stage}.dependencies.json"
}
# Preserve the optimized worker and wrap its existing artifact checks.
eval "$(declare -f gwas_post_magma | sed '1s/gwas_post_magma/gwas_post_magma_without_contract/')"
gwas_post_magma() {
    [[ ",${DO_STEP}," == *,magma,* ]] || return 0
    # Resolve per-SNP mode from the same mapped rows the worker will consume.
    # Explicit/sidecar N can reuse a receipt without remapping the GWAS.
    if [[ -z "${MAGMA_N-}" && ! -s "$GWAS_DIR/$GWAS.magma.N" ]]; then
        gwas_post_prepare_magma_rows || return $?
    fi
    if [[ "$REPLACE" != TRUE ]] && gwas_contract check magma; then
        gwas_post_log "MAGMA dependencies and outputs verified: $MAGMA_PREFIX"
        return 0
    fi
    local REPLACE=TRUE magma_identity
    magma_identity=$(gwas_contract identity magma) || return $?
    gwas_post_magma_without_contract || return $?
    if [[ "$(gwas_contract identity magma)" != "$magma_identity" ]]; then
        echo "ERROR: MAGMA dependencies changed during execution; refusing receipt" >&2
        return 1
    fi
    gwas_contract seal magma
}

# One pure resolver is used by both receipt snapshot and execution.
gwas_magma_n_argument() {
    env "MAGMA_N=${MAGMA_N-}" "GWAS_N=${GWAS_N-}" "GWAS_DIR=$GWAS_DIR" "GWAS=$GWAS" \
        "GWAS_POST_MAGMA_ROWS=${GWAS_POST_MAGMA_ROWS-}" \
        python3 "${PERF_F%/*}/format.contract.py" resolve-n
}

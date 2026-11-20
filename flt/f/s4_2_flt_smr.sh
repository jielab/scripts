#!/usr/bin/env bash
set -euo pipefail
root="$FLT_OUT/04smr"
run_dir=$(mktemp -d "$root/run.XXXXXXXX")
for tissue in eQTLGen Whole_Blood; do
    ref=$FLT_EQTLGEN
    if [[ $tissue == Whole_Blood ]]; then ref=$FLT_GTEX; fi
    for trait in G732 mdd_2025; do
        prefix="$run_dir/${tissue}_${trait}"
        "$FLT_SMR" --bfile "$MR_LD_PREFIX" --gwas-summary "$FLT_SMR_INPUT/$trait.ma" \
            --beqtl-summary "$ref" --out "$prefix" --cis-wind 2000 --peqtl-smr 5e-8 \
            --heidi-max-m 20 --diff-freq 0.20 --diff-freq-prop 0.05 --thread-num "$WORKFLOW_CORES" > "$prefix.stdout" 2>&1
        [[ -s "$prefix.smr" ]] || { echo "Missing SMR output: $prefix.smr" >&2; exit 1; }
    done
done
"$FLT_PYTHON" "$FLT_CODE/collect_smr.py" "$run_dir"

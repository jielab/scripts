#!/usr/bin/env bash
# TRACE action controls. Backend identity also selects the saved output tree;
# extraction and summarization never need a CUDA device.

gu_trace_needs_gpu() {
    [[ $1 == cuda && ( $2 == run || $2 == infer ) ]]
}

gu_trace_configure_device() {
    case "$TRACE_BACKEND" in
        cpu|cuda) ;;
        *) echo "ERROR: unknown TRACE backend: $TRACE_BACKEND" >&2; return 2 ;;
    esac
    if gu_trace_needs_gpu "$TRACE_BACKEND" "$ACTION"; then
        export WORKFLOW_DEVICE=gpu
        JIN=1
        python3 "$F/trace.py" infer --gpu-check
    else
        export WORKFLOW_DEVICE=cpu CUDA_VISIBLE_DEVICES=''
    fi
}

gu_trace_replace_outputs() {
    local action=$1 replace=$2 output=$3 sample_map=$4 folder
    local -a folders=()
    [[ $replace == 1 ]] || return 0
    case "$action" in
        run|extract) folders=(extract datafiles infer calls final) ;;
        infer) folders=(infer calls final) ;;
        segments) folders=(calls final) ;;
        check) return 0 ;;
        *) echo "ERROR: invalid TRACE action: $action" >&2; return 2 ;;
    esac
    for folder in "${folders[@]}"; do
        rm -rf -- "$output/$folder"
        mkdir -p -- "$output/$folder"
    done
    if [[ ( $action == run || $action == extract ) && $sample_map == "$output/samples/trace_sample_map.tsv" ]]; then
        rm -f -- "$sample_map"
    fi
}

gu_trace_prepare_partial() (
    # prepare_run refreshes only small manifests/maps/sample lists. Preserve
    # them as well as the authoritative record if input validation fails.
    # Numerical extraction/inference arrays are never copied or hashed here.
    set -euo pipefail
    local backup ready=0 name status
    local -a metadata=(manifest regions samples maps run.meta.tsv run.meta.next.tsv)
    backup=$(mktemp -d "$OUT/.partial-prepare.XXXXXX")
    restore_partial_metadata() {
        status=$?
        trap - EXIT
        if ((ready)) && { ((status != 0)) || [[ ${GU_TRACE_VALIDATE_ONLY:-0} == 1 ]]; }; then
            for name in "${metadata[@]}"; do
                rm -rf -- "$OUT/$name"
                [[ ! -e $backup/$name && ! -L $backup/$name ]] || cp -a -- "$backup/$name" "$OUT/$name"
            done
        fi
        rm -rf -- "$backup"
        exit "$status"
    }
    trap restore_partial_metadata EXIT
    for name in "${metadata[@]}"; do
        [[ ! -e $OUT/$name && ! -L $OUT/$name ]] || cp -a -- "$OUT/$name" "$backup/$name"
    done
    ready=1
    prepare_run
)

gu_trace_validate_classes() {
    local class class_chrs
    for class in autosomes chrX; do
        if [[ $class == autosomes ]]; then class_chrs="${autosomes[*]}"; else class_chrs=X; fi
        env GU_TRACE_VALIDATE_ONLY=1 TRACE_CHRS="$class_chrs" TRACE_OUT="$OUT/chromosome_classes/$class" TRACE_SAMPLE_MAP= \
            TRACE_FINAL_NAME=trace_haplotype_segments.internal.tsv.gz bash "$F/trace.sh" "$ACTION" || return $?
    done
}

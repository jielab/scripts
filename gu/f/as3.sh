#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
cat <<'HELP'
ArchaicSeeker3 analysis variables:
  AS3_DATA_IN=DIR              prepared input directory; must contain manifest.tsv
  AS3_RUNTIME=DIR              bundled AS3 runtime, default f/as3_runtime
  AS3_OUT=DIR                  default /mnt/d/analysis/gu/as3
  AS3_GPUS="0 1"               one serial worker per GPU
  AS3_CHRS="1 2 ... 22 X"      optional filter; X includes XPAR/XNONPAR_F/XNONPAR_M
  AS3_TARGET_CHUNK_SIZE=N      target-sample chunk size inside AS3
  AS3_MERGE=5000
  AS3_FORCE=1                  rerun completed tasks
  AS3_ALLOW_CPU=1              permit an intentional CPU fallback; default 0

Run through ./gu.sh as3 [run|check] --chr LIST; preparation is automatic.
The manifest may contain multiple cross-fit
panels. Every target sample is assigned to one panel per analysis unit. chrX rows
are marked experimental because the public pretrained model is not documented as
chrX-calibrated.
HELP
exit 0
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
action=${GU_ACTION:-as3_run}
runtime=${AS3_RUNTIME:-$script_dir/as3_runtime}
data=${AS3_DATA_IN:-${dir_archaic:-${dir_ref:-/mnt/d/data.BIG/refGen}/archaic/GRCH${GRCH:-38}}/as3/preprocessed}
out=${AS3_OUT:-/mnt/d/analysis/gu/as3}
gpus=${AS3_GPUS:-0}
merge=${AS3_MERGE:-5000}
force=${AS3_FORCE:-0}
allow_cpu=${AS3_ALLOW_CPU:-0}
if [[ -n ${AS3_EXEC:-} ]]; then
  exe=$AS3_EXEC
else
  exe=""
  for cand in "$runtime/ArchaicSeeker3.1-mamba" "$runtime/ArchaicSeeker3.0-mamba" "$runtime/ArchaicSeeker3-mamba"; do
    [[ -f $cand ]] && { exe=$cand; break; }
  done
  [[ -n $exe ]] || exe=$runtime/ArchaicSeeker3.1-mamba
fi
chr_filter=${AS3_CHRS:-}
manifest=$data/manifest.tsv

[[ $action == as3_run || $action == as3_check ]] || { echo "ERROR unknown AS3 action: $action" >&2; exit 1; }
[[ $force == 0 || $force == 1 ]] || { echo "ERROR AS3_FORCE must be 0 or 1" >&2; exit 1; }
[[ $allow_cpu == 0 || $allow_cpu == 1 ]] || { echo "ERROR AS3_ALLOW_CPU must be 0 or 1" >&2; exit 1; }
[[ $merge =~ ^[1-9][0-9]*$ ]] || { echo "ERROR AS3_MERGE must be a positive integer" >&2; exit 1; }
if [[ -n ${AS3_TARGET_CHUNK_SIZE:-} && ! ${AS3_TARGET_CHUNK_SIZE} =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR AS3_TARGET_CHUNK_SIZE must be a positive integer" >&2
    exit 1
fi
for cmd in python bcftools sort comm awk; do command -v "$cmd" >/dev/null || { echo "ERROR missing command: $cmd" >&2; exit 1; }; done
[[ -s $manifest ]] || { echo "ERROR AS3 manifest missing after automatic preparation: $manifest" >&2; exit 1; }
[[ -f $exe ]] || { echo "ERROR AS3 executable missing: $exe" >&2; exit 1; }
for required in \
    "$runtime/src/dataloaders/__init__.py" \
    "$runtime/src/models/__init__.py" \
    "$runtime/src/stepsagnostic/__init__.py" \
    "$runtime/exp/Basemodel_Mamba_4096/args.pckl" \
    "$runtime/exp/Basemodel_Mamba_4096/models/best_model.pth" \
    "$runtime/exp/Smoother_512_Kernel_8192/args.pckl" \
    "$runtime/exp/Smoother_512_Kernel_8192/models/best_model.pth"; do
    [[ -s $required ]] || { echo "ERROR bundled AS3 runtime file missing: $required" >&2; exit 1; }
done
python -c 'import torch,mamba_ssm,pandas,pysam' 2>/dev/null || { echo "ERROR AS3 Python dependencies are incomplete" >&2; exit 1; }
read -r -a gpu_array <<< "$gpus"
(( ${#gpu_array[@]} > 0 )) || { echo "ERROR AS3_GPUS is empty" >&2; exit 1; }
for gpu in "${gpu_array[@]}"; do [[ $gpu =~ ^[0-9]+$ ]] || { echo "ERROR invalid GPU ID: $gpu" >&2; exit 1; }; done
mkdir -p "$out"/{log,results}

unit_selected() {
    local unit=$1 c
    [[ -z $chr_filter ]] && return 0
    for c in $chr_filter; do
        c=${c#chr}
        if [[ $c == X || $c == 23 ]]; then
            [[ $unit == X* ]] && return 0
        elif [[ $unit == "$c" ]]; then
            return 0
        fi
    done
    return 1
}

# Load selected manifest rows. Field order is defined by corrected as3_prep.sh.
tasks=()
declare -A seen_keys=()
while IFS=$'\t' read -r panel unit target ref map target_list afr_list experimental; do
    [[ $panel == panel_id ]] && continue
    unit_selected "$unit" || continue
    key=${panel}.chr${unit}
    [[ -z ${seen_keys[$key]:-} ]] || { echo "ERROR duplicate AS3 manifest key: $key" >&2; exit 1; }
    seen_keys[$key]=1
    tasks+=("$panel"$'\t'"$unit"$'\t'"$target"$'\t'"$ref"$'\t'"$map"$'\t'"$target_list"$'\t'"$afr_list"$'\t'"$experimental")
done < "$manifest"
(( ${#tasks[@]} > 0 )) || { echo "ERROR no manifest tasks selected" >&2; exit 1; }

validate_task() {
    local panel=$1 unit=$2 target=$3 ref=$4 map=$5 target_list=$6 afr_list=$7 experimental=$8 key work nt nr
    key=${panel}.chr${unit}
    for f in "$target" "$ref" "$map" "$target_list" "$afr_list"; do [[ -s $f ]] || { echo "ERROR $key missing input: $f" >&2; return 1; }; done
    bcftools view -h "$target" >/dev/null 2>&1 || { echo "ERROR unreadable target VCF: $target" >&2; return 1; }
    bcftools view -h "$ref" >/dev/null 2>&1 || { echo "ERROR unreadable reference VCF: $ref" >&2; return 1; }
    nt=$(bcftools index -n "$target" 2>/dev/null || echo 0)
    nr=$(bcftools index -n "$ref" 2>/dev/null || echo 0)
    (( nt > 0 && nr > 0 )) || { echo "ERROR empty/unindexed VCF for $key" >&2; return 1; }
    [[ $nt -eq $nr ]] || { echo "ERROR target/reference variant-count mismatch for $key: $nt vs $nr" >&2; return 1; }
    work=$(mktemp -d "${TMPDIR:-/tmp}/as3_check_${key}.XXXXXX")
    bcftools query -l "$target" | sort -u > "$work/target.samples"
    bcftools query -l "$ref" | sort -u > "$work/ref.samples"
    sort -u "$target_list" > "$work/target.expected"
    sort -u "$afr_list" > "$work/afr.expected"
    awk 'NF{print $1}' "$map" | sort -u > "$work/map.samples"
    comm -12 "$work/target.samples" "$work/ref.samples" > "$work/overlap"
    if [[ -s $work/overlap ]]; then echo "ERROR target/reference overlap for $key" >&2; head -10 "$work/overlap" >&2; rm -rf "$work"; return 1; fi
    cmp -s "$work/target.samples" "$work/target.expected" || { echo "ERROR target sample-list mismatch for $key" >&2; rm -rf "$work"; return 1; }
    cmp -s "$work/ref.samples" "$work/map.samples" || { echo "ERROR reference.map mismatch for $key" >&2; rm -rf "$work"; return 1; }
    comm -23 "$work/afr.expected" "$work/ref.samples" > "$work/missing_afr"
    if [[ -s $work/missing_afr ]]; then echo "ERROR African-reference samples missing from reference VCF for $key" >&2; head -10 "$work/missing_afr" >&2; rm -rf "$work"; return 1; fi
    rm -rf "$work"
    if [[ $experimental == 1 ]]; then echo "[AS3 $key] WARNING: chrX prediction is experimental; pretrained-model calibration on chrX is not documented" >&2; fi
}

status=0
for task in "${tasks[@]}"; do
    IFS=$'\t' read -r panel unit target ref map target_list afr_list experimental <<< "$task"
    validate_task "$panel" "$unit" "$target" "$ref" "$map" "$target_list" "$afr_list" "$experimental" || status=1
done
(( status == 0 )) || { echo "ERROR AS3 manifest/input validation failed" >&2; exit 1; }

if [[ $action == as3_check ]]; then
    python - <<'PY'
import torch
print("CUDA available:", torch.cuda.is_available())
print("CUDA GPU count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("CUDA GPUs:", [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())])
PY
    echo "CHECK PASSED: manifest=$manifest tasks=${#tasks[@]} requested_GPUs=${gpu_array[*]}"
    printf '%s\n' "${tasks[@]}" | awk -F'\t' '{print "  panel=" $1 " unit=" $2 " target=" $3 " reference=" $4}'
    exit 0
fi

if [[ $allow_cpu == 0 ]]; then
    python - "${gpu_array[@]}" <<'PY'
import sys, torch
if not torch.cuda.is_available():
    raise SystemExit("ERROR: CUDA is unavailable; set AS3_ALLOW_CPU=1 only for an intentional CPU test")
n = torch.cuda.device_count()
ids = [int(x) for x in sys.argv[1:]]
bad = [x for x in ids if x >= n]
if bad:
    raise SystemExit(f"ERROR: requested GPU IDs {bad} but only {n} CUDA devices are visible")
print("CUDA GPUs:", n, [torch.cuda.get_device_name(i) for i in range(n)])
PY
fi

completion_record() {
    local target=$1 ref=$2 map=$3
    printf 'target\t%s\nreference\t%s\nmap\t%s\nexecutable\t%s\nbase_model\t%s\nsmoother_model\t%s\nmerge\t%s\ntarget_chunk_size\t%s\n' \
        "$(stat -c '%n:%s:%Y' "$target")" \
        "$(stat -c '%n:%s:%Y' "$ref")" \
        "$(stat -c '%n:%s:%Y' "$map")" \
        "$(stat -c '%n:%s:%Y' "$exe")" \
        "$(stat -c '%n:%s:%Y' "$runtime/exp/Basemodel_Mamba_4096/models/best_model.pth")" \
        "$(stat -c '%n:%s:%Y' "$runtime/exp/Smoother_512_Kernel_8192/models/best_model.pth")" \
        "$merge" "${AS3_TARGET_CHUNK_SIZE:-}"
}

run_task() {
    local panel=$1 unit=$2 target=$3 ref=$4 map=$5 experimental=$6 gpu=$7 key o run_dir old marker n
    key=${panel}.chr${unit}
    o=$out/results/$panel/chr$unit
    mkdir -p "$out/results/$panel"
    marker=$o/.complete
    if [[ $force == 0 && -s $marker && -e $o/introgression_prediction.txt && -e $o/introgression_prediction.bed ]] && \
       cmp -s <(completion_record "$target" "$ref" "$map") "$marker"; then
        echo "[AS3 $key] SKIP"
        return 0
    fi
    run_dir=$out/results/$panel/.chr${unit}.run.$BASHPID
    old=$out/results/$panel/.chr${unit}.old.$BASHPID
    rm -rf "$run_dir" "$old"
    mkdir -p "$run_dir"
    local -a extra=()
    [[ -n ${AS3_TARGET_CHUNK_SIZE:-} ]] && extra+=(--target-chunk-size "$AS3_TARGET_CHUNK_SIZE")
    n=$(bcftools query -l "$target" | wc -l)
    echo "[AS3 $key] RUN gpu=$gpu samples=$n experimental_chrX=$experimental"
    if ! CUDA_VISIBLE_DEVICES=$gpu python "$exe" -t "$target" -r "$ref" -m "$map" --merge "$merge" -o "$run_dir" "${extra[@]}" > "$out/log/${key}.log" 2>&1; then
        echo "ERROR $key AS3 process failed; see $out/log/${key}.log" >&2
        rm -rf "$run_dir"
        return 1
    fi
    [[ -e $run_dir/introgression_prediction.txt && -e $run_dir/introgression_prediction.bed ]] || {
        echo "ERROR $key missing prediction outputs" >&2
        rm -rf "$run_dir"
        return 1
    }
    completion_record "$target" "$ref" "$map" > "$run_dir/.complete"
    [[ -e $o ]] && mv "$o" "$old"
    if ! mv "$run_dir" "$o"; then
        [[ -e $old ]] && mv "$old" "$o"
        return 1
    fi
    rm -rf "$old"
    echo "[AS3 $key] DONE output=$o"
}

status=0
for gpu_index in "${!gpu_array[@]}"; do
(
    worker_status=0
    for i in "${!tasks[@]}"; do
        (( i % ${#gpu_array[@]} == gpu_index )) || continue
        IFS=$'\t' read -r panel unit target ref map target_list afr_list experimental <<< "${tasks[$i]}"
        run_task "$panel" "$unit" "$target" "$ref" "$map" "$experimental" "${gpu_array[$gpu_index]}" || worker_status=1
    done
    exit "$worker_status"
) &
done
for pid in $(jobs -p); do wait "$pid" || status=1; done
(( status == 0 )) || { echo "ERROR one or more AS3 tasks failed; see $out/log" >&2; exit 1; }
cp "$manifest" "$out/input.manifest.tsv"
echo "ALL DONE: AS3 output=$out tasks=${#tasks[@]}"

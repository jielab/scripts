#!/bin/bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
write_only="${1:-0}"
root="$OUT"
read -r -a chroms <<< "$("$PYTHON" -c 'import os,sys; sys.path.insert(0,sys.argv[1]); from ld_pipeline import parse_chroms; print(*parse_chroms(os.environ["CHRS"]))' "$script_dir/f")"
mkdir -p "$root"
slots=$((THREADS + 1))

write_header() {
    printf '#!/bin/bash -l\nset -euo pipefail\ncd %q\n' "$script_dir"
    for name in GRCH INPUT_FORMAT PFILE_DIR VCF_DIR PANEL PLINK2 THREADS MEMORY_MB D_GRID_KB MAX_D_KB STABLE_CONFIRMATIONS REFERENCE_COUNT_POLICY REFERENCE_MAX_REL_DEVIATION COINCIDENT_POLICY LD_STORAGE MIN_FREE_GB LD_TARGETS; do
        if [[ -v "$name" ]]; then printf 'export %s=%q\n' "$name" "${!name}"; fi
    done
    printf 'export PARALLEL_CHROMS=1\nexport LD_RESUME_DIR=\n'
    printf 'if [[ -z "${LSB_JOBID:-}" ]]; then echo "Submit with bsub" >&2; exit 2; fi\n'
}

for chr in "${chroms[@]}"; do
    {
        write_header
        printf 'export OUT=%q\nexport LD_TMP_DIR=%q\n' "$root/chr$chr" "$root/chr$chr/tmp"
        printf 'if (( ${LSB_DJOB_NUMPROC:-0} < %d )); then echo "Insufficient CPU slots" >&2; exit 2; fi\n' "$slots"
        printf 'exec bash %q --execute --grch %q --chr %q\n' "$script_dir/ld.sh" "$GRCH" "$chr"
    } > "$root/chr$chr.cmd.tmp.$$"
    mv -- "$root/chr$chr.cmd.tmp.$$" "$root/chr$chr.cmd"
    echo "[LD] Written $root/chr$chr.cmd"
done
{
    write_header
    printf 'export OUT=%q\nexport CHRS=%q\n' "$root" "1-22"
    printf 'source %q\n' "$script_dir/f/env.sh"
    printf 'exec "$PYTHON" %q\n' "$script_dir/f/merge.py"
} > "$root/merge.cmd.tmp.$$"
mv -- "$root/merge.cmd.tmp.$$" "$root/merge.cmd"

if (( write_only )); then
    echo "[LD] Command files written only; no jobs submitted."
    echo "[LD] Run ./ld.sh --grch $GRCH --chr $CHRS to submit with merge dependency."
    exit 0
fi
bash "$script_dir/f/check.sh"
batch="ld${GRCH}"
record="$root/submitted.$batch.$(date +%Y%m%d_%H%M%S).$$.tsv"
printf 'task\tjob_id\tjob_name\n' > "$record"
submit() {
    local task="$1" dependency="$2" cores="$3" memory="$4" response name
    if [[ "$task" == merge ]]; then name="$batch.step2.merge"; else name="$batch.step1.$task"; fi
    local options=(-q ser -n "$cores" -R "span[hosts=1] rusage[mem=$memory]" -W 168:00
                   -J "$name" -o "$root/$task.%J.LOG" -e "$root/$task.%J.ERR")
    [[ -z "$dependency" ]] || options+=(-w "$dependency")
    response=$(bsub "${options[@]}" < "$root/$task.cmd")
    echo "$response"
    if [[ "$response" =~ \<([0-9]+)\> ]]; then
        job_id="${BASH_REMATCH[1]}"
    else
        echo '[LD] Cannot parse submitted job ID. Check bjobs before resubmitting.' >&2
        exit 2
    fi
    printf '%s\t%s\t%s\n' "$task" "$job_id" "$name" >> "$record"
}

for chr in "${chroms[@]}"; do
    submit "chr$chr" '' "$slots" 16384
done
# All chromosome jobs are submitted before this wildcard is resolved by LSF.
# Users must not submit overlapping runs with the same names.
submit merge "done(\"$batch.step1.chr*\")" 1 4096
echo "[LD] Submission IDs: $record"

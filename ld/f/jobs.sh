#!/bin/bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
write_only="${1:-0}"
root="$OUT"
read -r -a chroms <<< "$("$PYTHON" -c 'import os,sys; sys.path.insert(0,sys.argv[1]); from ld_pipeline import parse_chroms; print(*parse_chroms(os.environ["CHRS"]))' "$script_dir/f")"
mkdir -p "$root"
queue="${LD_QUEUE:-${HPC_CPU_QUEUE:-33}}"
slots=$((THREADS + 1))

write_header() {
    printf '#!/bin/bash -l\n'
    printf '#BSUB -q %s\n#BSUB -n %s\n#BSUB -W 144:00\n' "$queue" "$slots"
    printf '#BSUB -R "span[hosts=1] rusage[mem=%sGB/host]"\n' "$MEMORY_DEFAULT_GIB"
    printf '#BSUB -J ld%s.%s\n' "$GRCH" "${task_label:-merge}"
    printf '#BSUB -o "%s/%s.%%J.LOG"\n#BSUB -e "%s/%s.%%J.ERR"\n' "$root" "${task_label:-merge}" "$root" "${task_label:-merge}"
    printf 'set -euo pipefail\ncd %q\n' "$script_dir"
    for name in DIR_DAT DIR_BIG DIR_TMP DIR_DATA DIR_SCRATCH DIR_GEN ANALYSIS_ROOT DATE_SCRATCH GRCH INPUT_FORMAT PFILE_DIR VCF_DIR PANEL PLINK2 THREADS MEMORY_MB D_GRID_KB MAX_D_KB AUTO_EXTEND_D STABLE_CONFIRMATIONS REFERENCE_COUNT_POLICY REFERENCE_MAX_REL_DEVIATION COINCIDENT_POLICY LD_STORAGE MIN_FREE_GB LD_TARGETS; do
        if [[ -v "$name" ]]; then printf 'export %s=%q\n' "$name" "${!name}"; fi
    done
    printf 'export PARALLEL_CHROMS=1\nexport LD_RESUME_DIR=\n'
    printf 'export submit=%q SUBMIT=%q WORKFLOW_WORKER=1\n' "$submit" "$submit"
    if [[ $submit == hpc ]]; then
        printf 'if [[ -z "${LSB_JOBID:-}" ]]; then echo "Submit with bsub" >&2; exit 2; fi\n'
    fi
}

for chr in "${chroms[@]}"; do
    {
        task_label=chr$chr write_header
        printf 'export OUT=%q\nexport LD_TMP_DIR=%q\n' "$root/chr$chr" "$root/chr$chr/tmp"
        [[ $submit != hpc ]] || printf 'if (( ${LSB_DJOB_NUMPROC:-0} < %d )); then echo "Insufficient CPU slots" >&2; exit 2; fi\n' "$slots"
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
commands=()
for chr in "${chroms[@]}"; do commands+=("$root/chr$chr.cmd"); done
runner=("$PYTHON" "$script_dir/../0f/submit_jobs.py" --submit "$submit" --jobs-N "$jobs_N")
[[ $SUBMIT_WAIT != TRUE ]] || runner+=(--wait)
exec "${runner[@]}" --merge "$root/merge.cmd" "${commands[@]}"

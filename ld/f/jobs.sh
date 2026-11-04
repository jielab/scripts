#!/bin/bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
write_only="${1:-0}"
root="$OUT"
read -r -a chroms <<< "$("$PYTHON" -c 'import os,sys; sys.path.insert(0,sys.argv[1]); from ld_pipeline import parse_chroms; print(*parse_chroms(os.environ["CHRS"]))' "$script_dir/f")"
mkdir -p "$root"
queue="${LD_QUEUE:-${HPC_CPU_QUEUE:-33}}"
slots=$((THREADS + 1))
[[ $LD_MEMORY_GIB =~ ^[1-9][0-9]*$ && $LD_WALLTIME =~ ^[1-9][0-9]*:[0-5][0-9]$ ]] || {
    echo 'LD_MEMORY_GIB must be positive integer GiB; LD_WALLTIME must be HHH:MM.' >&2; exit 2;
}
(( MEMORY_MB <= LD_MEMORY_GIB * 1024 )) || {
    echo 'PLINK MEMORY_MB exceeds the total LD_MEMORY_GIB allocation.' >&2; exit 2;
}
if [[ $submit == hpc ]]; then
    queue=$("$PYTHON" - "$script_dir/../0f" "$queue" "$slots" "$LD_MEMORY_GIB" "$LD_WALLTIME" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from qiming_resources import cpu_queue
hours, minutes = map(int, sys.argv[5].split(':'))
try:
    print(cpu_queue(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), hours + minutes / 60))
except ValueError as exc:
    raise SystemExit(str(exc))
PY
    )
fi

write_header() {
    local header_slots=$slots header_memory=$LD_MEMORY_GIB header_wall=$LD_WALLTIME
    # Merge reads the completed BED/JSON outputs; it does not compute LD.
    if [[ ${task_label:-merge} == merge ]]; then
        header_slots=1; header_memory=4; header_wall=1:00
    fi
    printf '#!/bin/bash -l\n'
    printf '#BSUB -q %s\n#BSUB -n %s\n#BSUB -W %s\n' "$queue" "$header_slots" "$header_wall"
    printf '#BSUB -M %sGB\n' "$header_memory"
    printf '#BSUB -R "span[hosts=1] rusage[mem=%sGB/host]"\n' "$header_memory"
    printf '#BSUB -J ld%s.%s\n' "$GRCH" "${task_label:-merge}"
    printf '#BSUB -o "%s/%s.%%J.LOG"\n#BSUB -e "%s/%s.%%J.ERR"\n' "$root" "${task_label:-merge}" "$root" "${task_label:-merge}"
    printf 'set -euo pipefail\ncd %q\n' "$script_dir"
    for name in DIR_WORK DIR_DAT DIR_BIG DIR_TMP DIR_DATA DIR_SCRATCH DIR_GEN DIR_SOFTWARE ANALYSIS_ROOT DATE_SCRATCH GRCH INPUT_FORMAT PFILE_DIR VCF_DIR PANEL PLINK2 THREADS MEMORY_MB LD_MEMORY_GIB LD_WALLTIME LD_PYTHON LD_CONDA_SH LD_INPUT_MANIFEST D_GRID_KB MAX_D_KB AUTO_EXTEND_D STABLE_CONFIRMATIONS REFERENCE_COUNT_POLICY REFERENCE_MAX_REL_DEVIATION COINCIDENT_POLICY LD_STORAGE MIN_FREE_GB LD_TARGETS; do
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

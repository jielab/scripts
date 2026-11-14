#!/usr/bin/env bash
# Called only by the recorded workflow worker, never during planning.
set -euo pipefail
umask 077

fail() { echo "[EMS120] ERROR: $*" >&2; exit 2; }
[[ ${WORKFLOW_WORKER:-} == 1 ]] || fail "Use ems120.sh; this is an internal worker."
# bsub uses LC_ALL=C for its machine-readable receipt. Chinese R identifiers
# and workbook columns require UTF-8 after entering the allocated worker.
if [[ $(locale charmap 2>/dev/null || true) != UTF-8 ]]; then
    for utf8_locale in C.UTF-8 en_US.UTF-8; do
        if [[ $(LC_ALL=$utf8_locale locale charmap 2>/dev/null || true) == UTF-8 ]]; then
            export LANG=$utf8_locale LC_ALL=$utf8_locale
            break
        fi
    done
fi
[[ $(locale charmap 2>/dev/null || true) == UTF-8 ]] || fail "A UTF-8 locale is required for Chinese EMS columns and R source."
export EMS120_PHASE=${EMS120_PHASE:?}
case $EMS120_PHASE in prepare|phone|geo|dx|report) ;; *) fail "Unknown EMS120_PHASE: $EMS120_PHASE" ;; esac
export EMS120_DEVICE=${EMS120_DEVICE:-cpu}
case $EMS120_DEVICE in cpu|gpu) ;; *) fail "EMS120_DEVICE must be cpu or gpu" ;; esac
if [[ $EMS120_DEVICE == gpu ]]; then
    [[ $EMS120_PHASE == phone || $EMS120_PHASE == geo || $EMS120_PHASE == dx ]] || fail "Only ML stages can request GPUs."
    [[ ${WORKFLOW_GPUS:-0} == 1 ]] || fail "GPU execution requires exactly one recorded GPU allocation."
else
    [[ ${WORKFLOW_GPUS:-0} == 0 ]] || fail "CPU stages must request zero GPUs."
    export CUDA_VISIBLE_DEVICES=
fi
if [[ ${submit:-hpc} == hpc ]]; then
    [[ -n ${LSB_JOBID:-} ]] || fail "HPC analysis requires an LSF allocation."
    [[ ${LSB_DJOB_NUMPROC:-0} -ge ${WORKFLOW_CORES:-1} ]] || fail "Insufficient allocated CPU slots."
fi

: "${EMS120_SCRIPT_DIR:?}" "${EMS120_SHARED_DIR:?}" "${EMS120_ANALYSIS_DIR:?}" "${EMS120_MODEL_DIR:?}"
for file in f/ems120.R f/ems120.py f/session.R; do
    [[ -s $EMS120_SCRIPT_DIR/$file ]] || fail "Missing project file: $file"
done

# Lock canonical output roots across distinct plans and symlink spellings.
# flock is nonblocking: an accidental duplicate exits instead of occupying a
# GPU while waiting. A shared external --dir-model has its own lock as well.
if [[ ${EMS120_PREFLIGHT:-0} != 1 ]]; then
    command -v flock >/dev/null 2>&1 || fail "flock (util-linux) is required."
    mkdir -p -- "$EMS120_ANALYSIS_DIR" "$EMS120_MODEL_DIR"
    analysis_real=$(cd -- "$EMS120_ANALYSIS_DIR" && pwd -P)
    model_real=$(cd -- "$EMS120_MODEL_DIR" && pwd -P)
    exec 8>"$analysis_real/.ems120.lock"
    flock -n 8 || fail "Another EMS120 job is using analysis directory $analysis_real"
    if [[ $model_real != "$analysis_real" ]]; then
        exec 9>"$model_real/.ems120.lock"
        flock -n 9 || fail "Another EMS120 job is using model directory $model_real"
    fi
    for helper in phenotype.R plotting.R png_output.R png_output.py; do
        [[ -s $EMS120_SHARED_DIR/$helper ]] || fail "Missing shared helper: $EMS120_SHARED_DIR/$helper"
    done
    # A cache-only redraw may not need the original Excel directory. Retain
    # that existing behavior; the R stages still validate inputs they need.
    if [[ $EMS120_PHASE == report ]]; then
        [[ -s $EMS120_ANALYSIS_DIR/dat/dat.list.rds ]] || fail "CPU report requires dat/dat.list.rds; it cannot train models."
    elif [[ ! -d ${EMS120_DATA_DIR:?}/clean && ! -s $EMS120_ANALYSIS_DIR/dat/dat0.list.rds ]]; then
        fail "Missing raw input directory $EMS120_DATA_DIR/clean and no dat/dat0.list.rds cache is available."
    fi
fi

find_conda() {
    local c
    if [[ -n ${EMS120_CONDA_EXE:-} ]]; then
        command -v "$EMS120_CONDA_EXE" || return 1
        return 0
    fi
    for c in "${CONDA_EXE:-}" "${DIR_WORK:-}/software/gu/conda/bin/conda" "${DIR_WORK:-}/software/anaconda3/bin/conda" \
             "${DIR_WORK:-}/software/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
             "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda" /opt/conda/bin/conda; do
        if [[ -n $c && -x $c ]]; then printf '%s\n' "$c"; return 0; fi
    done
    command -v conda
}
conda_exe=$(find_conda) || fail "Conda not found; supply --conda-exe /absolute/path/bin/conda."
conda_base=$("$conda_exe" info --base) || fail "Could not resolve Conda base."
set +u
if [[ -f $conda_base/etc/profile.d/conda.sh ]]; then
    source "$conda_base/etc/profile.d/conda.sh"
else
    conda_hook=$("$conda_exe" shell.bash hook) || fail "Could not initialize Conda."
    eval "$conda_hook"
fi
conda activate "${EMS120_CONDA_ENV:-ai}" || fail "Cannot activate existing Conda environment ${EMS120_CONDA_ENV:-ai}."
set -u
[[ -n ${CONDA_PREFIX:-} && -x $CONDA_PREFIX/bin/python ]] || fail "Activated environment has no Python."
export RETICULATE_PYTHON=$CONDA_PREFIX/bin/python
export RETICULATE_AUTOCONFIGURE=FALSE PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
unset PYTHONHOME
export PNG_PYTHON=$RETICULATE_PYTHON PYTHON_BIN=$RETICULATE_PYTHON
if [[ -n ${EMS120_PYTHON_LIB:-} ]]; then
    [[ -d $EMS120_PYTHON_LIB ]] || fail "Missing project Python library: $EMS120_PYTHON_LIB"
    export PYTHONPATH="$EMS120_PYTHON_LIB${PYTHONPATH:+:$PYTHONPATH}"
fi
rscript=${EMS120_RSCRIPT:-}
if [[ -z $rscript ]]; then rscript=$(command -v Rscript || true); fi
[[ -n $rscript ]] && command -v "$rscript" >/dev/null 2>&1 || fail "Rscript not found in the selected environment; use --rscript."
rscript=$(command -v "$rscript")
rprefix=$(cd -- "$(dirname -- "$rscript")/.." && pwd -P)
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# A separate R environment needs its own newer OpenSSL/C++ runtime first.
ssl_prefix=$CONDA_PREFIX
if [[ -d $rprefix/lib/R && $rprefix != "$CONDA_PREFIX" ]]; then
    export LD_LIBRARY_PATH="$rprefix/lib:$LD_LIBRARY_PATH"
    export PATH="$rprefix/bin:$PATH"
    ssl_prefix=$rprefix
fi
export EMS120_SSL_PREFIX=$ssl_prefix
# Preserve the original local R/OpenSSL compatibility workaround, scoped to
# the chosen environment. WSL driver libraries are used only where present.
preload=()
for lib in "$ssl_prefix/lib/libcrypto.so.3" "$ssl_prefix/lib/libssl.so.3"; do
    [[ ! -f $lib ]] || preload+=("$lib")
done
if [[ -f /usr/lib/wsl/lib/libcuda.so.1 && -f /usr/lib/wsl/lib/libnvidia-ml.so.1 ]]; then
    export LD_LIBRARY_PATH="/usr/lib/wsl/lib:$LD_LIBRARY_PATH"
    preload+=(/usr/lib/wsl/lib/libcuda.so.1 /usr/lib/wsl/lib/libnvidia-ml.so.1)
fi
if ((${#preload[@]})); then
    preload_join=$(IFS=:; printf '%s' "${preload[*]}")
    export LD_PRELOAD="$preload_join${LD_PRELOAD:+:$LD_PRELOAD}"
fi

export EMS120_CORES=${EMS120_CORES:-${WORKFLOW_CORES:-8}}
[[ $EMS120_CORES =~ ^[1-9][0-9]*$ && $EMS120_CORES -le ${WORKFLOW_CORES:-8} ]] || fail "Compute threads exceed the allocated CPU slots."
export OMP_NUM_THREADS=$EMS120_CORES OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TOKENIZERS_PARALLELISM=false MPLBACKEND=Agg
if [[ $EMS120_DEVICE == cpu ]]; then export CUDA_VISIBLE_DEVICES=; fi
if [[ ${submit:-local} == local && $EMS120_DEVICE == gpu && ! -v CUDA_VISIBLE_DEVICES ]]; then
    # Select one local card only when the caller has not selected a device.
    # Never replace the scheduler's device assignment on HPC.
    export CUDA_VISIBLE_DEVICES=0
fi
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
printf '[EMS120] mode=%s job=%s cores=%s device=%s\n' "${submit:-local}" "${LSB_JOBID:-local}" "$EMS120_CORES" "$EMS120_DEVICE"
printf '[EMS120] Python=%s Rscript=%s\n' "$RETICULATE_PYTHON" "$rscript"
printf '[EMS120] data=%s analysis=%s models=%s\n' "$EMS120_DATA_DIR" "$EMS120_ANALYSIS_DIR" "$EMS120_MODEL_DIR"

# Compile without writing into the source checkout (which may be read-only).
"$RETICULATE_PYTHON" -c 'import pathlib,sys; p=pathlib.Path(sys.argv[1]); compile(p.read_bytes(),str(p),"exec")' "$EMS120_SCRIPT_DIR/f/ems120.py"
if [[ $EMS120_PHASE == phone || $EMS120_PHASE == geo || $EMS120_PHASE == dx ]]; then
    "$rscript" --vanilla "$EMS120_SCRIPT_DIR/f/session.R" check
else
    echo "[EMS120] Data/report phase=$EMS120_PHASE; CUDA hidden, ML runtime check skipped."
fi
if [[ ${EMS120_PREFLIGHT:-0} == 1 ]]; then
    "$rscript" --vanilla "$EMS120_SCRIPT_DIR/f/preflight.R"
    echo "[EMS120] Preflight completed; production analysis and models were not generated."
    exit 0
fi


cd -- "$EMS120_ANALYSIS_DIR"
args=(--start_step "${START_STEP:-all}")
[[ -z ${END_STEP:-} ]] || args+=(--end_step "$END_STEP")
"$rscript" --vanilla "$EMS120_SCRIPT_DIR/f/ems120.R" "${args[@]}"

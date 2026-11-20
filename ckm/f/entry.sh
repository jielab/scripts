# Shared CKM worker options.
config=$CKM_ROOT/config/ckm.json
out_dir=${CKM_OUT_DIR:-$ANALYSIS_ROOT/ckm}
threads=${CKM_THREADS:-$WORKFLOW_CORES}
rscript=${CKM_RSCRIPT:-Rscript}
if [[ $submit == hpc && -z ${CKM_RSCRIPT:-} ]]; then
    rscript=Rscript
fi
extra=()
module=$1; shift
while (($#)); do
    case "$1" in
        -h|--help|help) usage; exit 0 ;;
        --config) config=${2:?--config requires a file}; shift 2 ;;
        --out-dir|--output-dir|--analysis-root) out_dir=${2:?--out-dir requires a directory}; shift 2 ;;
        --threads) threads=${2:?--threads requires N}; shift 2 ;;
        --rscript) rscript=${2:?--rscript requires an interpreter}; shift 2 ;;
        --phenotype-dir|--stage3-method|--definition|--domains)
            extra+=("$1" "${2:?option requires a value}"); shift 2 ;;
        *) echo "Unknown CKM option: $1" >&2; exit 2 ;;
    esac
done
(( threads <= WORKFLOW_CORES )) || { echo '--threads exceeds --cores' >&2; exit 2; }
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo '--threads must be a positive integer' >&2; exit 2; }
mkdir -p -- "$out_dir/tmp"
out_dir=$(cd -- "$out_dir" && pwd -P)
export CKM_THREADS=$threads OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=$threads TMPDIR=$out_dir/tmp PYTHONDONTWRITEBYTECODE=1
if [[ -z ${CKM_PYTHON:-} && -x $(dirname -- "$rscript")/python3 ]]; then
    export CKM_PYTHON=$(dirname -- "$rscript")/python3
fi
export CKM_PYTHON=${CKM_PYTHON:-${SHARED_PYTHON:-python3}}
command -v "$rscript" >/dev/null || { echo "Rscript not found: $rscript. Run install.sh or use --rscript." >&2; exit 2; }
exec "$rscript" "$CKM_ROOT/f/run.R" --module "$module" --config "$config" --out-dir "$out_dir" --threads "$threads" "${extra[@]}"

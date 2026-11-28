#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
Build the allele3 manuscript from completed analyses.

  ./allele3.paper.sh                      Show help.
  ./allele3.paper.sh run                  Build the manuscript; auto-detect local/LSF.
  ./allele3.paper.sh run --plan-only      Write jobs/paper.cmd and display paths.
  ./allele3.paper.sh run --wait TRUE      Wait for the LSF job to finish.
  ./allele3.paper.sh report               Alias of run.

Input and output:
  --analysis-root DIR    Project directory [ANALYSIS_ROOT/allele3].
  --result-root DIR      Parent of avi_window_enrichment/, alt_impact_similarity/, hotspots/, etc.
                         Default lookup: main/ only.
  --out DIR              Dedicated report directory [PROJECT/paper].
  --scan-root DIR        3.allele3_scan; contains gnomad/ and allofus/.
  --avi-root DIR         Completed 4.avi_scan.
  --herv-root DIR        HERV analysis directory (default: analysis/allele3/2.herv_scan).
  --pca-root DIR         Permanent data/gen/ancient/38/pca3d.
  --targets CSV          Cohorts [gnomad,allofus]. One cohort per estimate.
  --near BP,BP           Paired proximal distance band [0,1000].
  --far BP,BP            Paired distal distance band [100000,1000000].
  --bootstrap N --seed N Resampling for the paired comparison [1000, 20261115].
  --python PATH          Python with NumPy and Matplotlib [allele3-evolution].
  --queue NAME           Optional LSF queue. Job prefix: allele3.paper.
  --cores N              LSF CPU slots [1]. Memory and runtime follow queue policy.
  --submit auto|local|hpc  Optional override. Environment comes from ../0f/0os.sh.
  --plan-only            Prepare and print the job command without running it.

Files: paper/main/results/*.tsv, paper/main/reports/{index.html,*.png},
       paper/{jobs,logs,qc,tmp}/. The report reads existing scores and scans.
       To compute new analyses use allele3.sh avi_window_enrichment or allele3.sh alt_impact_similarity.
HELP
}
action=${1:-help}
case $action in
    help|-h|--help) usage; exit 0 ;;
    run|report|paper) shift ;;
    *) printf 'Unknown command: %s\n' "$action" >&2; exit 2 ;;
esac
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
args=()
mode=${submit:-${SUBMIT:-auto}}; plan=0; slots=1; wait_for_job=FALSE; queue_name=; output=; python_bin=${ALLELE_PYTHON:-}
while (($#)); do
    if [[ $1 == --*=* ]]; then option=$1; shift; set -- "${option%%=*}" "${option#*=}" "$@"; fi
    case $1 in
        --help|-h|--help-full) usage; exit 0 ;;
        --plan-only) plan=1; shift; continue ;;
        --submit|--python|--cores|--target|--queue|--wait|--analysis-root|--out|--result-root|--pca-root|--avi-root|--herv-root|--scan-root|--near|--far|--targets|--seed|--grch|--bootstrap)
            [[ $# -ge 2 && -n $2 && $2 != --* ]] || { echo "Missing value for $1" >&2; exit 2; } ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    case $1 in
        --cores)
            [[ $2 =~ ^[1-9][0-9]{0,3}$ ]] || { echo "$1 needs a positive integer" >&2; exit 2; }
            slots=$2 ;;
        --target) args+=(--targets "$2") ;;
        --submit) mode=$2 ;;
        --python) python_bin=$2 ;;
        --queue) queue_name=$2 ;;
        --wait) wait_for_job=${2^^} ;;
        --analysis-root) export ALLELE_ANALYSIS_ROOT="$(realpath -m -- "$2")" ;;
        --out) output=$(realpath -m -- "$2") ;;
        *) args+=("$1" "$2") ;;
    esac
    shift 2
done
[[ $wait_for_job == TRUE || $wait_for_job == FALSE ]] || { echo '--wait takes TRUE or FALSE' >&2; exit 2; }
export submit="$mode"
source "$root/../0f/0os.sh"
set -a
source "${ALLELE_CONFIG:-$root/allele3.env}"
set +a
project=$(realpath -m -- "$ALLELE_ANALYSIS_ROOT")
output=${output:-"$project/paper"}
if [[ -z $python_bin && -x ${ALLELE_PYTHON:-/nonexistent} ]]; then python_bin=$ALLELE_PYTHON; fi
if [[ -z $python_bin ]]; then
    for candidate in "${CONDA_PREFIX:-/nonexistent}/bin/python" "$DIR_WORK/.conda/envs/allele3-evolution/bin/python" "$DIR_HOME/anaconda3/envs/allele3-evolution/bin/python"; do
        if [[ -x $candidate && ( $candidate != "${CONDA_PREFIX:-}/bin/python" || ${CONDA_DEFAULT_ENV:-} != base ) ]]; then python_bin=$candidate; break; fi
    done
fi
python_bin=$(command -v "${python_bin:-python3}")
case "$project/" in "$output/"*) echo 'Choose a dedicated manuscript output directory' >&2; exit 2 ;; esac
[[ $output != "$project/main" ]] || { echo 'main/ contains scientific products' >&2; exit 2; }
mkdir -p -- "$output"/{logs,jobs,qc,tmp}
job=$(mktemp "$output/tmp/command-XXXXXXXX")
trap 'rm -f -- "$job"' EXIT
{
    printf '#!/usr/bin/env bash\nset -euo pipefail\ncd -- %q\n' "$PWD"
    printf 'export MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=%q\n' "$output/tmp/matplotlib"
    printf 'exec %q ' "$python_bin"
    printf '%q ' "$root/f/paper.py" --analysis-root "$project" --out "$output" "${args[@]}"
    printf '\n'
} > "$job"
record=$(mktemp "$output/jobs/.command-XXXXXXXX")
cp -- "$job" "$record"
mv -f -- "$record" "$output/jobs/paper.cmd"
echo "$output/jobs/paper.cmd"
if ((plan)); then cat -- "$job"; exit 0; fi
if [[ $submit == hpc ]]; then
    queue_name=${queue_name:-${ALLELE_JOB_QUEUE:-}}
    options=(-J "allele3.paper.${job##*-}" -n "$slots" -oo "$output/logs/paper.%J.log" -eo "$output/logs/paper.%J.err")
    [[ -z $queue_name ]] || options+=(-q "$queue_name")
    [[ $wait_for_job != TRUE ]] || options+=(-K)
    bsub "${options[@]}" < "$job"
else
    bash "$job" 2>&1 | tee -a "$output/logs/paper.log"
fi

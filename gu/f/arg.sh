#!/usr/bin/env bash
# Internal ARG workflow; public options live in ../arg.sh.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# Wait for every PID explicitly; short successful jobs must not turn into wait -n=127.
refgen_children=()
refgen_stop_tree(){
  local pid=$1 child
  while read -r child; do [[ -z $child ]] || refgen_stop_tree "$child"; done < <(pgrep -P "$pid" || true)
  kill -TERM "$pid" 2>/dev/null || true
}
refgen_stop_jobs(){
  local pid
  for pid in "${refgen_children[@]}"; do refgen_stop_tree "$pid"; done
  for pid in "${refgen_children[@]}"; do wait "$pid" 2>/dev/null || true; done
  refgen_children=()
}
run_parallel(){
  local fn=$1 c pid failed=0 limit=${GRID_ARG_JOBS:-${JOBS:-1}}; shift
  refgen_children=()
  for c in "$@"; do
    "$fn" "$c" & refgen_children+=("$!")
    if ((${#refgen_children[@]}>=limit)); then
      for pid in "${refgen_children[@]}"; do wait "$pid" || failed=1; done
      refgen_children=()
      ((failed==0)) || { echo "ERROR: $fn failed; no further jobs started" >&2; return 1; }
    fi
  done
  for pid in "${refgen_children[@]}"; do wait "$pid" || failed=1; done
  refgen_children=()
  ((failed==0)) || { echo "ERROR: $fn failed" >&2; return 1; }
}
# AA-aware VCF preparation shared by prep_gen and build.
prepared_record(){
  local c=$1 in=$IN_DIR/chr$1.vcf.gz index
  [[ -s $in ]] || { echo "ERROR: missing phased VCF: $in" >&2; return 1; }
  if [[ -s $in.tbi ]]; then index=$in.tbi; elif [[ -s $in.csi ]]; then index=$in.csi; else echo "ERROR: missing VCF index: $in" >&2; return 1; fi
  python3 "$HERE/arg.py" cache --file "$in" --file "$index" --text-file "$HERE/arg.sh" --text-file "$HERE/arg.py" --value "missing_max=$MISS;mac_min=$MAC;build=$REFGEN_GRCH"
  if [[ $c == X ]]; then
    python3 "$HERE/arg.py" cache --text-file "$HERE/arg.py"
    [[ ! -s ${REFGEN_PFILE_DIR}/chrX.psam ]] || python3 "$HERE/arg.py" cache --text-file "${REFGEN_PFILE_DIR}/chrX.psam"
  fi
}
prep_one(){
  local c=$1 in=$IN_DIR/chr$1.vcf.gz out meta expected replace_word=FALSE
  out=$(prepared_vcf "$c"); meta=$out.input.meta.tsv; expected=$WORK/chr$c.prepared.meta.tsv
  prepared_record "$c" > "$expected"
  if [[ $REPLACE == FALSE && -s $out && -s $out.tbi && -s $meta ]] && cmp -s "$expected" "$meta"; then
    python3 "$HERE/arg.py" cache --file "$out" --file "$out.tbi" > "$expected.outputs"
    if cmp -s "$expected.outputs" "$out.outputs.json"; then echo "SKIP verified ARG VCF: $out"; return; fi
  fi
  rm -f "$meta" "$out.outputs.json"
  [[ $REPLACE == FALSE && ! -e $out ]] || replace_word=TRUE
  local -a x_args=()
  if [[ $c == X && -s ${REFGEN_PFILE_DIR}/chrX.psam ]]; then x_args+=(--x-psam "${REFGEN_PFILE_DIR}/chrX.psam"); fi
  bash "$HERE/arg.sh" prep_vcf --in "$in" --out "$out" --threads "$THREADS" --missing-max "$MISS" --mac-min "$MAC" --tmp-root "$(tmp_for "$c")" --replace "$replace_word" --chr "$c" --build "$REFGEN_GRCH" "${x_args[@]}"
  python3 "$HERE/arg.py" cache --file "$out" --file "$out.tbi" > "$out.outputs.json.next"
  mv -f "$out.outputs.json.next" "$out.outputs.json"
  mv -f "$expected" "$meta"
}

arg_cli() (
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

usage(){ exec bash "$HERE/../arg.sh" --help; }

die(){ echo "ERROR: $*" >&2; exit 2; }
bool(){ case "${1^^}" in TRUE|T|YES|Y|1|ON) echo TRUE;; FALSE|F|NO|N|0|OFF) echo FALSE;; *) die "expected TRUE or FALSE, got: $1";; esac; }
need(){ [[ $# -ge 2 && -n ${2:-} && ${2:-} != --* ]] || die "$1 requires a value"; }

method=needle
format=native
action=
dir_gen=
arg_dir=
dir_pfile=
dir_vcf=
map_dir=
map_pattern=
ancestry_file=
sample_file=
scratch=
custom_keep=
arg_needle_home=
chrs=
grch=37
replace=FALSE
full=FALSE
affinity=FALSE
max_individuals=20000
seed_haplotypes=4000
anchors_per_pop=1000
threads=8
jobs=1
missing_max=0.05
missing_set=0
mac_min=2
sample_match_min_work=50000000
gen4arg_dir=
extract=
needle_mode=auto
maf_min=0.001
normalize=TRUE
normalize_demography=
seed=20260904
features=auto
x_odd_male=drop-last
threads_mode=wgs
threads_demography=
threads_fit=FALSE
mutation_rate=1.25e-8

while (($#)); do
  case "$1" in
    --method) need "$@"; method=${2,,}; shift 2;;
    --format) need "$@"; format=${2,,}; shift 2;;
    --action) need "$@"; action=${2,,}; shift 2;;
    --dir-gen) need "$@"; dir_gen=${2%/}; shift 2;;
    --arg-dir) need "$@"; arg_dir=${2%/}; shift 2;;
    --dir-pfile|--hap-dir|--arg-hap-dir) need "$@"; dir_pfile=${2%/}; shift 2;;
    --dir-vcf|--target-vcf-dir) need "$@"; dir_vcf=${2%/}; shift 2;;
    --map-dir|--arg-map-dir) need "$@"; map_dir=${2%/}; shift 2;;
    --map-pattern|--arg-map-pattern) need "$@"; map_pattern=$2; shift 2;;
    --ancestry-file) need "$@"; ancestry_file=$2; shift 2;;
    --sample-file|--sample-panel) need "$@"; sample_file=$2; shift 2;;
    --scratch|--arg-scratch) need "$@"; scratch=${2%/}; shift 2;;
    --keep|--arg-keep) need "$@"; custom_keep=$2; shift 2;;
    --arg-needle-home) need "$@"; arg_needle_home=${2%/}; shift 2;;
    --chr|--chrs) need "$@"; chrs=$2; shift 2;;
    --grch) need "$@"; grch=${2,,}; shift 2;;
    --replace|--replace-arg) need "$@"; replace=$(bool "$2"); shift 2;;
    --full|--arg-full) need "$@"; full=$(bool "$2"); shift 2;;
    --affinity|--arg-affinity) need "$@"; affinity=$(bool "$2"); shift 2;;
    --max-individuals|--arg-max-individuals) need "$@"; max_individuals=$2; shift 2;;
    --seed-haplotypes|--arg-seed-haplotypes) need "$@"; seed_haplotypes=$2; shift 2;;
    --anchors-per-pop|--arg-anchors-per-pop) need "$@"; anchors_per_pop=$2; shift 2;;
    --threads|--arg-threads) need "$@"; threads=$2; shift 2;;
    --jobs|--arg-jobs) need "$@"; jobs=$2; shift 2;;
    --missing-max) need "$@"; missing_max=$2; missing_set=1; shift 2;;
    --mac-min) need "$@"; mac_min=$2; shift 2;;
    --sample-match-min-work) need "$@"; sample_match_min_work=$2; shift 2;;
    --gen4arg-dir) need "$@"; gen4arg_dir=${2%/}; shift 2;;
    --keep-snv|--extract) need "$@"; extract=$2; shift 2;;
    --needle-mode) need "$@"; needle_mode=$2; shift 2;;
    --maf-min) need "$@"; maf_min=$2; shift 2;;
    --normalize) need "$@"; normalize=$(bool "$2"); shift 2;;
    --normalize-demography) need "$@"; normalize_demography=$2; shift 2;;
    --seed) need "$@"; seed=$2; shift 2;;
    --features) need "$@"; features=$(bool "$2"); shift 2;;
    --x-odd-male) need "$@"; x_odd_male=$2; shift 2;;
    --threads-mode) need "$@"; threads_mode=$2; shift 2;;
    --demography) need "$@"; threads_demography=$2; shift 2;;
    --threads-fit-to-data) need "$@"; threads_fit=$(bool "$2"); shift 2;;
    --mutation-rate) need "$@"; mutation_rate=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "unknown build option: $1";;
  esac
done

case "$method" in needle|tsinfer|threads) ;; *) die "--method must be needle, tsinfer or threads";; esac
case "$x_odd_male" in error|drop-last) ;; *) die '--x-odd-male must be error or drop-last';; esac
if [[ $method != tsinfer && $missing_set == 0 ]]; then missing_max=0; fi
if [[ $method == tsinfer && -n $extract ]]; then die '--keep-snv is currently a Needle option; tsinfer uses AA-aware sequence variants'; fi
case "$format" in native|trace) ;; *) die "--format must be native or trace";; esac
case "$grch" in 37|b37|grch37|hg19) grch=37;; 38|b38|grch38|hg38) grch=38;; *) die "--grch must be 37 or 38";; esac
for pair in "threads:$threads" "jobs:$jobs" "max-individuals:$max_individuals" "seed-haplotypes:$seed_haplotypes" "anchors-per-pop:$anchors_per_pop"; do
  name=${pair%%:*}; value=${pair#*:}; [[ $value =~ ^[0-9]+$ ]] || die "--$name must be a non-negative integer"
done
((threads>0 && jobs>0 && seed_haplotypes>0)) || die "--threads, --jobs and --seed-haplotypes must be positive"

dir_gen=${dir_gen:-/mnt/d/data.BIG/refGen/1kg/$grch}
[[ -d $dir_gen ]] || die "--dir-gen does not exist: $dir_gen"
base=${dir_gen##*/}
case "$base" in pfile|hap|typ|vcf|imp|gen) gen_root=$(dirname -- "$dir_gen");; *) gen_root=$dir_gen;; esac
if [[ -z $arg_dir ]]; then
  if [[ $method == threads ]]; then arg_dir=$gen_root/arg.threads; else arg_dir=$gen_root/arg; fi
fi
map_dir=${map_dir:-/mnt/d/data.BIG/refGen/maps/GRCh$grch}
sample_file=${sample_file:-$gen_root/samples.txt}
if [[ -z $ancestry_file ]]; then
  if [[ -s $gen_root/samples.txt ]]; then ancestry_file=$gen_root/samples.txt
  else ancestry_file=/mnt/d/data/ukb/phe/pca/ukb.ancestry.auto.tsv.gz
  fi
fi
if [[ -z $dir_pfile ]]; then
  if find "$dir_gen" -maxdepth 1 -type f -name 'chr*.pgen' -print -quit 2>/dev/null | grep -q .; then dir_pfile=$dir_gen
  elif find "$gen_root/hap" -maxdepth 1 -type f -name 'chr*.pgen' -print -quit 2>/dev/null | grep -q .; then dir_pfile=$gen_root/hap
  else dir_pfile=$gen_root/pfile
  fi
fi
dir_vcf=${dir_vcf:-$gen_root/vcf}
gen_root=$(realpath "$gen_root")
dir_pfile=$(realpath -m "$dir_pfile")
dir_vcf=$(realpath -m "$dir_vcf")
arg_dir=$(realpath -m "$arg_dir")
dataset=$(basename -- "$(dirname -- "$gen_root")")
identity=$(printf '%s\n' "$gen_root" "$dir_pfile" "$arg_dir" | sha256sum | cut -c1-16)
scratch=${scratch:-$HOME/refgen_arg_scratch/${dataset}.b${grch}.${identity}}
scratch=$(realpath -m "$scratch")
if [[ -z $gen4arg_dir ]]; then
  if [[ $method == threads ]]; then gen4arg_dir=${dir_pfile}.4arg.threads
  elif [[ $method == needle ]]; then gen4arg_dir=${dir_pfile}.4arg; else gen4arg_dir=${dir_vcf}.4arg; fi
fi
gen4arg_dir=$(realpath -m "$gen4arg_dir")
[[ $gen4arg_dir != "$dir_pfile" && $gen4arg_dir != "$dir_vcf" && $gen4arg_dir != "$arg_dir" ]] || die 'prepared, source and ARG directories must be distinct'
[[ -z $extract || -s $extract ]] || die "missing --keep-snv file: $extract (use /mnt/d/... in WSL)"
case "$needle_mode" in auto|array|sequence) ;; *) die '--needle-mode must be auto, array or sequence';; esac
if [[ $needle_mode == auto ]]; then
  if [[ ${dataset,,} == 1kg && -z $extract ]]; then needle_mode=sequence; else needle_mode=array; fi
fi
[[ -z $normalize_demography || -s $normalize_demography ]] || die 'missing --normalize-demography file'
[[ $seed =~ ^[0-9]+$ ]] || die '--seed must be a nonnegative integer'
python3 - "$maf_min" "$missing_max" "$mac_min" <<'PY'
import sys
assert 0 <= float(sys.argv[1]) < .5, 'maf-min must be in [0,.5)'
assert 0 <= float(sys.argv[2]) < 1, 'missing-max must be in [0,1)'
assert int(sys.argv[3]) >= 1, 'mac-min must be positive'
PY
chrs=${chrs:-1-22,X}
chrs=${chrs//,/ }; chrs=${chrs//;/ }
chrs=$(python3 - "$chrs" "$method" <<'PY'
import re,sys
out=[]
for token in sys.argv[1].split():
    token=re.sub('^chr','',token,flags=re.I).upper()
    if '-' in token:
        a,b=map(int,token.split('-'))
        if not 1 <= a <= b <= 22: raise SystemExit('Invalid chromosome range')
        out.extend(map(str,range(a,b+1)))
    else: out.append(token)
out=['X' if x.upper() in ('X','23') else x for x in out]
allowed=set(map(str,range(1,23))) | {'X'}
if not out or not set(out)<=allowed: raise SystemExit('Unsupported chromosome selection')
print(' '.join(dict.fromkeys(out)))
PY
)
if [[ ${REFGEN_GEN4ARG_ONLY:-0} == 1 ]]; then action=prepare; fi

if [[ $method == threads ]]; then
  action=${action:-all}; [[ $action != build ]] || action=all
  case "$action" in all|prepare|infer|convert|check) ;; *) die "bad threads --action=$action";; esac
  case "$threads_mode" in wgs|array) ;; *) die '--threads-mode must be wgs or array';; esac
  [[ -z $threads_demography || -s $threads_demography ]] || die 'missing --demography file'
  [[ $features != TRUE && $affinity != TRUE ]] || die 'Threads supports ARG and TRACE output; --features/--affinity are Needle-only'
  [[ -d $dir_pfile ]] || die "phased PGEN directory is missing: $dir_pfile"
elif [[ $method == needle ]]; then
  action=${action:-all}
  case "$action" in all|prepare|infer|convert|features|affinity|check) ;; *) die "bad needle --action=$action";; esac
  [[ -d $dir_pfile ]] || die "phased PGEN/BGEN directory is missing: $dir_pfile"
  [[ -d $map_dir || -n $map_pattern ]] || die "genetic maps are missing: $map_dir"
  [[ -s $ancestry_file ]] || die "ancestry table is missing: $ancestry_file"
else
  action=${action:-build}; [[ $action == all ]] && action=build
  case "$action" in build|prepare|check) ;; *) die "bad tsinfer --action=$action";; esac
  [[ -d $dir_vcf ]] || die "phased VCF directory is missing: $dir_vcf"
fi

export REFGEN_ARG_METHOD=$method REFGEN_ARG_ACTION=$action REFGEN_ARG_FORMAT=$format REFGEN_GEN_ROOT=$gen_root
export REFGEN_ARG_DIR=$arg_dir REFGEN_PFILE_DIR=$dir_pfile REFGEN_VCF_DIR=$dir_vcf
export REFGEN_MAP_DIR=$map_dir REFGEN_MAP_PATTERN=$map_pattern REFGEN_ANCESTRY_FILE=$ancestry_file
export REFGEN_SAMPLE_FILE=$sample_file REFGEN_ARG_SCRATCH=$scratch REFGEN_ARG_KEEP=$custom_keep
export REFGEN_ARG_NEEDLE_HOME=$arg_needle_home REFGEN_ARG_CHRS=$chrs REFGEN_GRCH=$grch
export REFGEN_REPLACE=$replace REFGEN_ARG_FULL=$full REFGEN_ARG_AFFINITY=$affinity
export REFGEN_MAX_INDIVIDUALS=$max_individuals REFGEN_SEED_HAPLOTYPES=$seed_haplotypes
export REFGEN_ANCHORS_PER_POP=$anchors_per_pop REFGEN_THREADS=$threads REFGEN_JOBS=$jobs
export REFGEN_MISSING_MAX=$missing_max REFGEN_MAC_MIN=$mac_min
export REFGEN_SAMPLE_MATCH_MIN_WORK=$sample_match_min_work
export REFGEN_GEN4ARG_DIR=$gen4arg_dir REFGEN_EXTRACT=$extract REFGEN_NEEDLE_MODE=$needle_mode
export REFGEN_ARG_MAF_MIN=$maf_min REFGEN_ARG_GENO_MAX=$missing_max REFGEN_SEED=$seed
export REFGEN_ARG_NORMALIZE=$normalize REFGEN_NORMALIZE_DEMOGRAPHY=$normalize_demography
export REFGEN_X_ODD_MALE=$x_odd_male
export REFGEN_THREADS_MODE=$threads_mode REFGEN_THREADS_DEMOGRAPHY=$threads_demography
export REFGEN_THREADS_FIT=$threads_fit REFGEN_THREADS_MUTATION_RATE=$mutation_rate
if [[ $features == auto ]]; then [[ $format == trace || $method == threads ]] && features=FALSE || features=TRUE; fi
export REFGEN_ARG_FEATURES=$features

if [[ $action == check ]]; then
  [[ -d $arg_dir ]] || die "ARG directory is missing: $arg_dir"
else
  mkdir -p "$gen4arg_dir"
  # Keep lock files permanently: unlinking a flock file permits two lock inodes.
  exec 8>"$gen4arg_dir/.refgen.lock"
  flock -n 8 || die "another build/prep_gen is still using prepared data: $gen4arg_dir. Wait for that task to finish, then rerun this command; verified inputs and inference checkpoints will be reused."
  if [[ $action != prepare ]]; then
    mkdir -p "$arg_dir"
    exec 9>"$arg_dir/.refgen.lock"
    flock -n 9 || die "another build is still using ARG output: $arg_dir. Wait for that task to finish, then rerun this command."
  fi
fi
if [[ $action == prepare ]]; then
  echo "prep_gen method=$method build=GRCh$grch chr=$chrs dir-gen=$gen_root output=$gen4arg_dir"
else
  echo "build method=$method format=$format action=$action build=GRCh$grch chr=$chrs dir-gen=$gen_root input=$gen4arg_dir output=$arg_dir"
fi
if [[ $method == needle ]]; then echo "Needle input mode=$needle_mode SNPs=${extract:-all-QC-passing} maf-min=$maf_min missing-max=$missing_max"; fi
bash "$HERE/arg.sh" "$method"
if [[ $action != prepare && $action != check && ! ( $method == threads && $action == infer ) ]]; then
  # Retain the input/output locks until validation finishes. A failed check
  # must propagate to the caller instead of reporting a successful build.
  export REFGEN_ARG_ACTION=check
  bash "$HERE/arg.sh" "$method"
  echo "BUILD PASS: output validation complete: $arg_dir"
fi

)

arg_threads() (
  local runtime=${GU_THREADS_ENV:-$HOME/.venvs/gu-threads}
  [[ -x $runtime/bin/python ]] || { echo "ERROR: Threads environment missing; run bash $HERE/install_threads.sh" >&2; exit 2; }
  export PATH="$runtime/bin:$HOME/miniforge3/envs/grid/bin:$HOME/anaconda3/envs/gu/bin:$PATH"
  export PYTHONNOUSERSITE=1
  exec "$runtime/bin/python" -u "$HERE/arg_threads.py"
)

arg_needle() (
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
SCRIPTS_ROOT=$(cd -- "$HERE/../.." && pwd -P)
GRID_ROOT=${REFGEN_GRID_ROOT:-$SCRIPTS_ROOT/grid}
GRID_F=$GRID_ROOT/f
[[ -s $GRID_F/common.sh ]] || { echo "ERROR: GRID helpers are missing: $GRID_F" >&2; exit 2; }

# Use the validated ARG-Needle environment without requiring an interactive
# `conda activate` in data-preparation jobs.
GRID_CONDA_ENV=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid-argneedle}
if [[ ! -x $GRID_CONDA_ENV/bin/python3 && -x $HOME/miniforge3/envs/grid/bin/python3 ]]; then
  GRID_CONDA_ENV=$HOME/miniforge3/envs/grid
fi
if [[ -x $GRID_CONDA_ENV/bin/python3 ]]; then export PATH="$GRID_CONDA_ENV/bin:$PATH"; fi
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null
[[ ! -d $GRID_CONDA_ENV/lib/R/library ]] || export R_LIBS_USER=$GRID_CONDA_ENV/lib/R/library

export GRID_ARG_ACTION=${REFGEN_ARG_ACTION:?}
export GRID_ARG_FORMAT=${REFGEN_ARG_FORMAT:-native}
export GRID_ARG_HAP_DIR=${REFGEN_PFILE_DIR:?}
export GRID_ARG_MAP_DIR=${REFGEN_MAP_DIR:-}
export GRID_ARG_MAP_PATTERN=${REFGEN_MAP_PATTERN:-}
export GRID_ARG_OUT=${REFGEN_ARG_DIR:?}
export GRID_ARG_TREES_DIR=$GRID_ARG_OUT/trees
export GRID_ARG_SCRATCH=${REFGEN_ARG_SCRATCH:-}
export GRID_ARG_FULL=${REFGEN_ARG_FULL:-FALSE}
export GRID_ARG_MAX_INDIVIDUALS=${REFGEN_MAX_INDIVIDUALS:-20000}
export GRID_ARG_SEED_HAPLOTYPES=${REFGEN_SEED_HAPLOTYPES:-4000}
export GRID_ARG_THREADS=${REFGEN_THREADS:-8}
export GRID_ARG_JOBS=${REFGEN_JOBS:-1}
export GRID_ARG_ANCHORS_PER_POP=${REFGEN_ANCHORS_PER_POP:-1000}
export GRID_ARG_KEEP=${REFGEN_ARG_KEEP:-}
export GRID_ARG_AFFINITY=${REFGEN_ARG_AFFINITY:-FALSE}
export GRID_ARG_NEEDLE_HOME=${REFGEN_ARG_NEEDLE_HOME:-}
export GRID_ANCESTRY_FILE=${REFGEN_ANCESTRY_FILE:?}
export GRID_CHRS=${REFGEN_ARG_CHRS:-1-22}
export GRID_REPLACE=${REFGEN_REPLACE:-FALSE}
export GRID_DRY_RUN=FALSE
export GRID_OUTPUT_ROOT=$GRID_ARG_OUT
export GRID_SEED=${REFGEN_SEED:-20260904}
export GRID_ANCHOR_PROB_MIN=${REFGEN_ANCHOR_PROB_MIN:-0.999}
export GRID_ARG_MAF_MIN=${REFGEN_ARG_MAF_MIN:-0.001}
export GRID_ARG_GENO_MAX=${REFGEN_ARG_GENO_MAX:-0.05}
export GRID_ARG_NORMALIZE=${REFGEN_ARG_NORMALIZE:-TRUE}
export GRID_ARG_WINDOW_BP=${REFGEN_ARG_WINDOW_BP:-1000000}

# shellcheck source=/dev/null
source "$GRID_F/common.sh"
grid_run(){ printf '[REFGEN ARG]'; printf ' %q' "$@"; printf '\n' >&2; "$@"; }
grid_run_logged(){
  local log=$1; shift
  mkdir -p "$(dirname -- "$log")"
  printf '[REFGEN ARG]'; printf ' %q' "$@"; printf '\n' | tee -a "$log" >&2
  "$@" >>"$log" 2>&1
}
if [[ -z $GRID_ARG_NEEDLE_HOME ]]; then
  GRID_ARG_NEEDLE_HOME=$(python3 - <<'PY' 2>/dev/null || true
import importlib.util
from pathlib import Path
s=importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
if s is not None and s.origin is not None: print(Path(s.origin).resolve().parents[2])
PY
)
fi
if [[ -z $GRID_ARG_NEEDLE_HOME && -s $HOME/software/arg-needle/arg-needle-scripts/arg_needle/infer_args_advanced.py ]]; then
  GRID_ARG_NEEDLE_HOME=$HOME/software/arg-needle/arg-needle-scripts
fi
read -r -a CHRS <<< "$(grid_expand_chrs "$GRID_CHRS")"
(( ${#CHRS[@]} )) || _grid_die 'No chromosomes selected'

WORK=${GRID_ARG_SCRATCH:+$GRID_ARG_SCRATCH/refgen_argneedle}
PREP=${REFGEN_GEN4ARG_DIR:?}
INF=${WORK:+$WORK/infer}
KEEP=$PREP/panel/arg.keep.txt
PANEL=$PREP/panel/arg.panel.tsv
ORDER=$PREP/panel/arg.order.txt
TRACE_DIR=$GRID_ARG_OUT/trace/needle

check_chr(){
  local c=$1 f
  for f in "$GRID_ARG_OUT/argn/chr$c.argn" "$GRID_ARG_TREES_DIR/chr$c.trees" \
           "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"; do
    [[ -s $f ]] || { echo "ERROR: missing ARG-Needle output for chr$c: $f" >&2; return 1; }
  done
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE ]]; then gzip -t "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"; fi
  python3 - "$GRID_ARG_TREES_DIR/chr$c.trees" "$c" <<'PY'
import sys,tskit
ts=tskit.load(sys.argv[1])
if min(ts.num_trees,ts.num_samples,ts.num_sites,ts.num_mutations) < 1:
    raise SystemExit(f"ERROR: empty ARG-Needle tree for chr{sys.argv[2]}: trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
print(f"ARG-Needle CHECK PASS chr{sys.argv[2]} trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}")
PY
  if [[ $c == X ]]; then
    python3 "$HERE/arg.py" needle-x check --tree "$GRID_ARG_TREES_DIR/chrX.trees" --identities "$PREP/chrX.haplotypes.tsv" --output "$GRID_ARG_TREES_DIR/chrX.sample_map.tsv" --psam "$(grid_sample_for X)" --keep "$PREP/panel/X/arg.keep.txt"
  fi
}

trace_chr(){
  local c=$1
  local -a cmd=(python3 "$HERE/arg.py" trace --method needle
    --tree "$GRID_ARG_TREES_DIR/chr$c.trees"
    --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
    --out-dir "$TRACE_DIR" --chr "$c")
  [[ $GRID_ARG_ACTION != check ]] || cmd+=(--check)
  [[ $GRID_REPLACE == FALSE ]] || cmd+=(--replace)
  grid_run "${cmd[@]}"
}

if [[ $GRID_ARG_ACTION == check ]]; then
  for c in "${CHRS[@]}"; do check_chr "$c"; done
  if [[ $GRID_ARG_FORMAT == trace ]]; then
    for c in "${CHRS[@]}"; do trace_chr "$c"; done
  fi
  exit 0
fi

[[ -n $GRID_ARG_SCRATCH ]] || _grid_die '--scratch is required for ARG-Needle construction'
case "$GRID_ARG_SCRATCH" in /mnt/[a-zA-Z]|/mnt/[a-zA-Z]/*) _grid_die '--scratch must use the native Linux filesystem, not /mnt/*';; esac
mkdir -p "$PREP/panel" "$PREP/log" "$INF"
exec 7>"$WORK/.lock"
flock -n 7 || _grid_die "another task uses scratch: $WORK"
owner=$(printf '%s\n' "$(realpath "$PREP")" "$(realpath -m "$GRID_ARG_OUT")")
if [[ -e $WORK/owner.txt ]]; then
  [[ $(cat "$WORK/owner.txt") == "$owner" ]] || _grid_die "scratch belongs to a different dataset/output; choose another --scratch: $WORK"
elif find "$INF" -type f -print -quit | grep -q . || [[ -d $WORK/prepare ]]; then
  _grid_die "unverified legacy scratch; choose a NEW --scratch directory: $WORK"
fi
printf '%s\n' "$owner" > "$WORK/owner.txt"
if [[ $GRID_ARG_ACTION != prepare ]]; then mkdir -p "$GRID_ARG_TREES_DIR" "$GRID_ARG_OUT/log"; fi

select_panel(){
  local dir=$PREP/panel
  [[ $1 != X ]] || dir=$dir/X
  KEEP=$dir/arg.keep.txt; PANEL=$dir/arg.panel.tsv; ORDER=$dir/arg.order.txt
}

make_panel_one(){
  local c=$1
  select_panel "$c"
  mkdir -p "$(dirname -- "$KEEP")"
  local sample; sample=$(grid_sample_for "$c" || true)
  [[ -s $sample ]] || _grid_die "Missing phased sample/PSAM for chr${CHRS[0]} under $GRID_ARG_HAP_DIR"
  local -a cmd=(python3 "$GRID_F/make_arg_keep.py" --sample "$sample" --ancestry "$GRID_ANCESTRY_FILE" --anchors-per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN" --seed "$GRID_SEED")
  if grid_is_true "$GRID_ARG_FULL"; then cmd+=(--full --max-individuals 0); else cmd+=(--max-individuals "$GRID_ARG_MAX_INDIVIDUALS"); fi
  [[ -z $GRID_ARG_KEEP ]] || cmd+=(--custom-keep "$GRID_ARG_KEEP")
  # Regenerate deterministically, but preserve timestamps when contents match.
  cmd+=(--out-keep "$KEEP.next" --out-panel "$PANEL.next")
  grid_run "${cmd[@]}"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.py" needle-x panel --keep "$KEEP.next" --psam "$sample" --policy "${REFGEN_X_ODD_MALE:-drop-last}" --qc "$PREP/panel/X/selection.qc.json"
    python3 - "$KEEP.next" "$PANEL.next" <<'PY'
import sys,pandas as pd
from pathlib import Path
wanted={r.split()[-1] for r in Path(sys.argv[1]).read_text().splitlines()[1:] if r.strip()}
p=pd.read_csv(sys.argv[2],sep='\t',dtype={'eid':str}); p[p.eid.isin(wanted)].to_csv(sys.argv[2],sep='\t',index=False)
PY
  fi
  for f in "$KEEP" "$PANEL"; do
    if cmp -s "$f.next" "$f"; then rm -f "$f.next"; else mv -f "$f.next" "$f"; fi
  done
  tail -n +2 "$KEEP" > "$ORDER.next"
  if cmp -s "$ORDER.next" "$ORDER"; then rm -f "$ORDER.next"; else mv -f "$ORDER.next" "$ORDER"; fi
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE || $GRID_ARG_AFFINITY == TRUE ]]; then
    python3 - "$PANEL" "$GRID_ARG_ANCHORS_PER_POP" "$GRID_ANCHOR_PROB_MIN" <<'PY'
import sys,pandas as pd
p=pd.read_csv(sys.argv[1],sep='\t'); prob=pd.to_numeric(p['prob'],errors='coerce')
for pop in ('AFR','EAS','EUR','SAS'):
    n=min(int(sys.argv[2]),int(((p['pop']==pop)&(prob.isna()|(prob>=float(sys.argv[3])))).sum()))
    if n<10: raise SystemExit(f'GRID features need >=10 anchors for {pop} (have {n}); use --features FALSE for ARG/TRACE only')
PY
  fi
}

make_panel(){
  local c made=0
  for c in "${CHRS[@]}"; do
    if [[ $c == X ]]; then make_panel_one X
    elif ((made==0)); then make_panel_one "$c"; made=1; fi
  done
}

prepare_record(){
  local c=$1 pf s m
  select_panel "$c"
  pf=$(grid_pfile_for "$c" || true); s=$(grid_sample_for "$c" || true); m=$(grid_map_for "$c" || true)
  local -a args=(--value "build=${REFGEN_GRCH};maf=$GRID_ARG_MAF_MIN;geno=$GRID_ARG_GENO_MAX;mode=${REFGEN_NEEDLE_MODE};seed=$GRID_SEED"
    --text-file "$KEEP" --text-file "$PANEL" --file "$m" --file "$s"
    --text-file "$HERE/arg.sh" --text-file "$HERE/arg.py" --text-file "$GRID_F/build_argneedle_map.py" --text-file "$GRID_F/validate_haps.py")
  if [[ -n $pf ]]; then
    args+=(--file "$pf.pgen")
    if [[ -s $pf.pvar.zst ]]; then args+=(--file "$pf.pvar.zst"); else args+=(--file "$pf.pvar"); fi
  else args+=(--file "$(grid_bgen_for "$c")"); fi
  [[ -z ${REFGEN_EXTRACT:-} ]] || args+=(--text-file "$REFGEN_EXTRACT")
  if [[ $c == X ]]; then
    args+=(--text-file "$PREP/panel/X/selection.qc.json" --value "x_odd_male=${REFGEN_X_ODD_MALE:-drop-last}")
  fi
  args+=(--file "$(command -v plink2)")
  python3 "$HERE/arg.py" cache "${args[@]}"
}
prepared_outputs(){
  local p=$PREP/chr$1
  local -a args=(--file "$p.haps.gz" --text-file "$p.sample" --text-file "$p.map" --text-file "$p.qc.json" --text-file "$p.positions.qc.json")
  if [[ $1 == X ]]; then args+=(--text-file "$p.haplotypes.tsv" --text-file "$p.packed.keep"); fi
  python3 "$HERE/arg.py" cache "${args[@]}"
}
require_prepared(){
  local c=$1 p=$PREP/chr$1
  prepare_record "$c" > "$p.request.next"
  [[ -s $p.request.json ]] && cmp -s "$p.request.next" "$p.request.json" || _grid_die "prepared chr$c does not match this request; run prep_gen with the SAME input/panel/filter options"
  prepared_outputs "$c" > "$p.outputs.next"
  cmp -s "$p.outputs.next" "$p.outputs.json" || _grid_die "prepared chr$c was modified or incomplete; run prep_gen"
  rm -f "$p.request.next" "$p.outputs.next"
}

export_chr(){
  local c=$1 b pf s srcmap prefix subset haps sampleout gz
  select_panel "$c"
  b=$(grid_bgen_for "$c" || true); pf=$(grid_pfile_for "$c" || true)
  s=$(grid_sample_for "$c" || true); srcmap=$(grid_map_for "$c" || true)
  [[ -s $b || -n $pf ]] || _grid_die "Missing phased BGEN or PGEN for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $s ]] || _grid_die "Missing Oxford sample or PSAM for chr$c in $GRID_ARG_HAP_DIR"
  [[ -s $srcmap ]] || _grid_die "Missing GRCh${REFGEN_GRCH:-37} genetic map for chr$c"
  prefix=$PREP/chr$c; subset=$prefix.subset
  prepare_record "$c" > "$prefix.request.next"
  if [[ $GRID_REPLACE == FALSE && -s $prefix.outputs.json ]] && cmp -s "$prefix.request.next" "$prefix.request.json"; then
    if prepared_outputs "$c" > "$prefix.outputs.next" && cmp -s "$prefix.outputs.next" "$prefix.outputs.json"; then
      echo "SKIP verified genetic input: $prefix"; rm -f "$prefix.request.next" "$prefix.outputs.next"; return
    fi
  fi
  rm -f "$prefix.request.json" "$prefix.outputs.json"
  rm -f "$prefix.haps" "$prefix.haps.gz" "$prefix.sample" "$prefix.map" "$prefix.log"
  local -a cmd
  if [[ -n $pf ]]; then
    cmd=(plink2 --pfile "$pf"); [[ ! -s $pf.pvar.zst ]] || cmd+=(vzs)
  else cmd=(plink2 --bgen "$b" ref-first --sample "$s"); fi
  cmd+=(--threads "$GRID_ARG_THREADS")
  if [[ $c == X ]]; then
    [[ -n $pf ]] || _grid_die 'X preparation requires a PGEN/PSAM input with known sex'
    local lo=2699521 hi=154931043
    if [[ $REFGEN_GRCH == 38 ]]; then lo=2781480; hi=155701382; fi
    cmd+=(--chr X --from-bp "$lo" --to-bp "$hi")
  fi
  [[ -z ${REFGEN_EXTRACT:-} ]] || cmd+=(--extract "$REFGEN_EXTRACT")
  cmd+=(--keep "$KEEP" --indiv-sort f "$ORDER" --snps-only just-acgt --max-alleles 2 --maf "$GRID_ARG_MAF_MIN" --geno "$GRID_ARG_GENO_MAX" --sort-vars --make-pgen vzs --out "$subset")
  grid_run_logged "$PREP/log/arg.prepare.chr$c.log" "${cmd[@]}"
  grid_run_logged "$PREP/log/arg.prepare.chr$c.log" plink2 --pfile "$subset" vzs --threads "$GRID_ARG_THREADS" --export haps --out "$prefix"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.py" needle-x pack --prefix "$prefix" --psam "$subset.psam" --keep "$KEEP" --build "$REFGEN_GRCH"
  fi
  grid_run python3 "$HERE/arg.py" unique-haps --haps "$prefix.haps" --chr "$c" --out "$prefix.positions.qc.json"
  haps=$(find "$PREP" -maxdepth 1 -type f \( -name "chr$c.haps" -o -name "chr$c.haps.gz" \) -print -quit)
  sampleout=$(find "$PREP" -maxdepth 1 -type f -name "chr$c.sample" -print -quit)
  [[ -s $haps && -s $sampleout ]] || _grid_die "PLINK2 did not create chr$c Oxford HAPS/sample"
  if [[ $haps != *.gz ]]; then
    if command -v pigz >/dev/null 2>&1; then grid_run pigz -f -p "$GRID_ARG_THREADS" "$haps"; else grid_run gzip -f "$haps"; fi
    gz=$haps.gz
  else gz=$haps; fi
  [[ $gz == "$prefix.haps.gz" ]] || mv -f "$gz" "$prefix.haps.gz"
  rm -f "$subset.pgen" "$subset.pvar" "$subset.pvar.zst" "$subset.psam" "$subset.log"
  grid_run python3 "$GRID_F/build_argneedle_map.py" --source "$srcmap" --haps "$prefix.haps.gz" --chr "$c" --out "$prefix.map"
  local validation_keep=$KEEP
  [[ $c != X ]] || validation_keep=$prefix.packed.keep
  grid_run python3 "$GRID_F/validate_haps.py" --haps "$prefix.haps.gz" --sample "$prefix.sample" --keep "$validation_keep" --map "$prefix.map" --out "$prefix.qc.json"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.py" needle-x qc --qc "$prefix.qc.json" --identities "$prefix.haplotypes.tsv"
  fi
  prepared_outputs "$c" > "$prefix.outputs.next"
  mv -f "$prefix.outputs.next" "$prefix.outputs.json"
  mv -f "$prefix.request.next" "$prefix.request.json"
}

infer_chr(){
  local c=$1 prefix nind seed final
  select_panel "$c"
  prefix=$PREP/chr$c
  require_prepared "$c"
  mkdir -p "$WORK/input"
  # Stage verified inputs on native Linux storage without changing identical files.
  for ext in haps.gz sample map; do
    if [[ ! -s $WORK/input/chr$c.$ext ]] || ! cmp -s "$prefix.$ext" "$WORK/input/chr$c.$ext"; then
      cp "$prefix.$ext" "$WORK/input/chr$c.$ext.next"; mv -f "$WORK/input/chr$c.$ext.next" "$WORK/input/chr$c.$ext"
    fi
  done
  prefix=$WORK/input/chr$c
  [[ -s $prefix.haps.gz && -s $prefix.sample && -s $prefix.map ]] || _grid_die "Run arg.sh prep_gen first for chr$c"
  nind=$(( $(wc -l < "$prefix.sample") - 2 )); seed=$GRID_ARG_SEED_HAPLOTYPES
  ((seed<=2*nind)) || seed=$((2*nind))
  ((2*nind>=300)) || _grid_die "ARG-Needle ASMC decoding requires at least 300 total haplotypes (selected $((2*nind)))"
  if [[ ${REFGEN_NEEDLE_MODE} == sequence && $seed -ge $((2*nind)) ]]; then seed=$((2*nind-1)); fi
  ((seed>=2)) || _grid_die "Too few scaffold haplotypes for chr$c"
  local chromosome=$c
  [[ $c != X ]] || chromosome=23
  if [[ $c == X ]]; then echo "Needle X: $(wc -l < "$PREP/chrX.haplotypes.tsv") identity rows (including header); normalization=$(chr_normalize "$c")"; fi
  local -a cmd=(python3 "$GRID_F/run_argneedle_advanced.py" --home "$GRID_ARG_NEEDLE_HOME" --haps "$prefix.haps.gz" --map "$prefix.map" --out "$INF/chr$c" --chr "$chromosome" --seed-haplotypes "$seed" --threads "$GRID_ARG_THREADS" --normalize "$(chr_normalize "$c")" --random-seed "$GRID_SEED" --log "$GRID_ARG_OUT/log/arg.infer.chr$c.log")
  [[ $GRID_REPLACE == FALSE ]] || cmd+=(--replace)
  cmd+=(--mode "$REFGEN_NEEDLE_MODE")
  [[ -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ]] || cmd+=(--normalize-demography "$REFGEN_NORMALIZE_DEMOGRAPHY")
  grid_run "${cmd[@]}"
  final=$INF/chr$c.argn; [[ -s $final ]] || _grid_die "Missing final ARG-Needle file $final"
  mkdir -p "$GRID_ARG_OUT/argn"
  if ! cmp -s "$final" "$GRID_ARG_OUT/argn/chr$c.argn"; then
    cp "$final" "$GRID_ARG_OUT/argn/chr$c.argn.next"; mv -f "$GRID_ARG_OUT/argn/chr$c.argn.next" "$GRID_ARG_OUT/argn/chr$c.argn"
  fi
  inferred_record "$c" > "$GRID_ARG_OUT/argn/chr$c.input.json.next"
  mv -f "$GRID_ARG_OUT/argn/chr$c.input.json.next" "$GRID_ARG_OUT/argn/chr$c.input.json"
}

chr_normalize(){
  # The bundled CEU normalization is not an X-specific demography.
  if [[ $GRID_ARG_NORMALIZE != TRUE || ( $1 == X && -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ) ]]; then echo 0; else echo 1; fi
}

inferred_record(){
  local c=$1
  local -a args=(--file "$GRID_ARG_OUT/argn/chr$c.argn" --text-file "$PREP/chr$c.request.json" --text-file "$PREP/chr$c.outputs.json"
    --value "mode=$REFGEN_NEEDLE_MODE;normalize=$(chr_normalize "$c");scaffold=$GRID_ARG_SEED_HAPLOTYPES;seed=$GRID_SEED")
  [[ -z ${REFGEN_NORMALIZE_DEMOGRAPHY:-} ]] || args+=(--file "$REFGEN_NORMALIZE_DEMOGRAPHY")
  python3 "$HERE/arg.py" cache "${args[@]}"
}
require_inferred(){
  local c=$1
  require_prepared "$c"
  inferred_record "$c" > "$WORK/chr$c.inferred.next"
  cmp -s "$WORK/chr$c.inferred.next" "$GRID_ARG_OUT/argn/chr$c.input.json" || _grid_die "ARG chr$c does not match prepared input/options; run infer first"
  rm -f "$WORK/chr$c.inferred.next"
}

convert_chr(){
  local c=$1 argn=$GRID_ARG_OUT/argn/chr$1.argn trees=$GRID_ARG_TREES_DIR/chr$1.trees sampleout=$PREP/chr$1.sample
  select_panel "$c"
  [[ -s $argn && -s $sampleout ]] || _grid_die "Missing chr$c ARG/HAPS preparation output"
  require_inferred "$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.convert.chr$c.log" python3 "$GRID_F/argn_to_trees.py" --argn "$argn" --haps "$PREP/chr$c.haps.gz" --out "$trees" --home "$GRID_ARG_NEEDLE_HOME"
  if [[ $c == X ]]; then
    grid_run python3 "$HERE/arg.py" needle-x map --tree "$trees" --identities "$PREP/chrX.haplotypes.tsv" --output "$GRID_ARG_TREES_DIR/chrX.sample_map.tsv"
  else
    grid_run python3 "$GRID_F/make_sample_map.py" --trees "$trees" --sample "$sampleout" --out "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
  fi
  if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE || $GRID_ARG_AFFINITY == TRUE ]]; then
    grid_run python3 "$GRID_F/make_anchors.py" --panel "$PANEL" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --out "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --per-pop "$GRID_ARG_ANCHORS_PER_POP" --prob-min "$GRID_ANCHOR_PROB_MIN"
  fi
}

features_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.variants.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $PREP/chr$c.haps.gz && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing tree/HAPS/anchors for chr$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.features.chr$c.log" python3 "$GRID_F/arg_features.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --haps "$PREP/chr$c.haps.gz" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --chr "$c" --window-bp "$GRID_ARG_WINDOW_BP" --out "$out"
}

affinity_chr(){
  local c=$1 out=$GRID_ARG_TREES_DIR/chr$1.affinity.tsv.gz
  [[ -s $GRID_ARG_TREES_DIR/chr$c.trees && -s $GRID_ARG_TREES_DIR/chr$c.sample_map.tsv && -s $GRID_ARG_TREES_DIR/chr$c.anchors.tsv ]] || _grid_die "Missing ARG conversion files for chr$c"
  grid_run_logged "$GRID_ARG_OUT/log/arg.affinity.chr$c.log" python3 "$GRID_F/arg_affinity.py" --trees "$GRID_ARG_TREES_DIR/chr$c.trees" --sample-map "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" --anchors "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" --out "$out"
}

trap 'refgen_stop_jobs; exit 130' INT
trap 'refgen_stop_jobs; exit 143' TERM
trap 'refgen_stop_jobs; exit 129' HUP

case "$GRID_ARG_ACTION" in
  prepare) make_panel; run_parallel export_chr "${CHRS[@]}";;
  infer) make_panel; run_parallel infer_chr "${CHRS[@]}";;
  convert) make_panel; run_parallel convert_chr "${CHRS[@]}";;
  features) make_panel; run_parallel features_chr "${CHRS[@]}";;
  affinity) run_parallel affinity_chr "${CHRS[@]}";;
  all)
    make_panel
    run_parallel export_chr "${CHRS[@]}"
    run_parallel infer_chr "${CHRS[@]}"
    run_parallel convert_chr "${CHRS[@]}"
    if [[ ${REFGEN_ARG_FEATURES:-TRUE} == TRUE ]]; then run_parallel features_chr "${CHRS[@]}"; fi
    if grid_is_true "$GRID_ARG_AFFINITY"; then run_parallel affinity_chr "${CHRS[@]}"; fi
    ;;
esac

if [[ $GRID_ARG_FORMAT == trace && ( $GRID_ARG_ACTION == all || $GRID_ARG_ACTION == convert ) ]]; then
  run_parallel trace_chr "${CHRS[@]}"
  printf 'method\tneedle\nformat\ttrace\nbuild\tGRCh%s\nchromosomes\t%s\n' \
    "${REFGEN_GRCH:-37}" "${CHRS[*]}" > "$TRACE_DIR/ARG_TRACE_BUILD.tsv"
fi

if [[ $GRID_ARG_ACTION == prepare ]]; then echo "Genetic input ready: $PREP"; exit 0; fi
cat > "$GRID_ARG_OUT/ARG_NEEDLE_BUILD.txt" <<META
created=$(date -Is)
producer=/mnt/d/scripts/gu/arg.sh build
method=needle
build=GRCh${REFGEN_GRCH:-37}
genotype_dir=$GRID_ARG_HAP_DIR
chromosomes=${CHRS[*]}
full=$GRID_ARG_FULL
max_individuals=$GRID_ARG_MAX_INDIVIDUALS
seed_haplotypes=$GRID_ARG_SEED_HAPLOTYPES
scratch=$GRID_ARG_SCRATCH
META

echo "ARG-Needle data preparation complete: $GRID_ARG_OUT"

)

arg_tsinfer() (
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
GU_F=$HERE
for f in arg.py; do
  [[ -s $GU_F/$f ]] || { echo "ERROR: ARG helper is missing: $GU_F/$f" >&2; exit 2; }
done

GU_ENV=${GU_ENV_PREFIX:-$HOME/anaconda3/envs/gu}
[[ ! -x $GU_ENV/bin/python ]] || export PATH="$GU_ENV/bin:$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ACTION=${REFGEN_ARG_ACTION:-build}
IN_DIR=${REFGEN_VCF_DIR:?}
OUT_DIR=${REFGEN_ARG_DIR:?}
ARG_VCF_DIR=${REFGEN_GEN4ARG_DIR:?}
THREADS=${REFGEN_THREADS:-8}
export JOBS=${REFGEN_JOBS:-1}
SAMPLE_MATCH_MIN_WORK=${REFGEN_SAMPLE_MATCH_MIN_WORK:-50000000}
export MISS=${REFGEN_MISSING_MAX:-0.05}
export MAC=${REFGEN_MAC_MIN:-2}
REPLACE=${REFGEN_REPLACE:-FALSE}
FORMAT=${REFGEN_ARG_FORMAT:-native}
TRACE_DIR=$OUT_DIR/trace/tsinfer
TMP_PARENT=${REFGEN_TSINFER_TMP_ROOT:-/tmp/refgen/arg/tsinfer/$(basename -- "${REFGEN_GEN_ROOT:?}").b${REFGEN_GRCH:-37}}
WORK=
if [[ $ACTION == check ]]; then
  [[ -d $OUT_DIR ]] || { echo "ERROR: ARG directory is missing: $OUT_DIR" >&2; exit 2; }
else
  mkdir -p "$ARG_VCF_DIR" "$TMP_PARENT"
  [[ $ACTION == prepare ]] || mkdir -p "$OUT_DIR"
  TMP_PARENT=$(cd -- "$TMP_PARENT" && pwd -P)
  WORK=$(mktemp -d "$TMP_PARENT/run.XXXXXX")
  cleanup(){ case "$WORK" in "$TMP_PARENT"/run.*) rm -rf -- "$WORK";; *) echo "WARNING: refusing to remove unexpected temporary path: $WORK" >&2;; esac; }
  trap cleanup EXIT
fi

CHRS=$(python3 - "${REFGEN_ARG_CHRS:-1-22,X}" <<'PY'
import re,sys
out=[]
for x in re.split(r"[,;\s]+",sys.argv[1].strip()):
    x=re.sub(r"^chr","",x,flags=re.I)
    m=re.fullmatch(r"(\d+)-(\d+)",x)
    if m: out += [str(i) for i in range(int(m[1]),int(m[2])+1)]
    elif x: out.append(x.upper())
print(" ".join(dict.fromkeys(out)))
PY
)

prepared_vcf(){ printf '%s/chr%s.vcf.gz\n' "$ARG_VCF_DIR" "$1"; }
tmp_for(){ printf '%s/chr%s\n' "$WORK" "$1"; }

build_record(){
  local c=$1
  prepared_record "$c"
  python3 "$HERE/arg.py" cache --text-file "$HERE/arg.sh" --text-file "$GU_F/arg.py" --text-file "$HERE/comm.py" --file "$(prepared_vcf "$c")" --file "$(prepared_vcf "$c").input.meta.tsv"
  printf 'sample_match_min_work\t%s\nrecombination_rate\t%s\nmutation_rate\t%s\none_bit\t%s\n' \
    "$SAMPLE_MATCH_MIN_WORK" "${REFGEN_RECOMB_RATE:-1e-8}" "${REFGEN_MUT_RATE:-1.25e-8}" "${REFGEN_ARG_ONE_BIT:-0}"
}
build_one(){
  local c=$1 in out meta expected
  local -a extra=()
  in=$(prepared_vcf "$c"); out=$OUT_DIR/chr$c.trees; meta=$out.input.meta.tsv; expected=$WORK/chr$c.arg.meta.tsv
  build_record "$c" > "$expected"
  [[ ${REFGEN_ARG_ONE_BIT:-0} != 1 ]] || extra+=(--one-bit)
  if [[ $REPLACE == FALSE && -s $out && -s $meta && -s $OUT_DIR/chr$c.sample_map.tsv && -s $OUT_DIR/chr$c.arg_qc.json ]] && cmp -s "$expected" "$meta"; then
    python3 "$HERE/arg.py" cache --file "$out" --file "$OUT_DIR/chr$c.sample_map.tsv" --file "$OUT_DIR/chr$c.arg_qc.json" > "$expected.outputs"
    if cmp -s "$expected.outputs" "$out.outputs.json"; then check_one "$c"; echo "SKIP verified tsinfer ARG: $out"; return; fi
  fi
  rm -f "$meta" "$out.outputs.json"
  python3 "$GU_F/arg.py" infer-tsinfer --chr "$c" --vcf "$in" --out "$out" --sample-map "$OUT_DIR/chr$c.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr$c.arg_qc.json" --threads "$THREADS" --tmp-root "$(tmp_for "$c")" \
    --sample-match-min-work "$SAMPLE_MATCH_MIN_WORK" --recombination-rate "${REFGEN_RECOMB_RATE:-1e-8}" \
    --mutation-rate "${REFGEN_MUT_RATE:-1.25e-8}" "${extra[@]}"
  python3 "$GU_F/arg.py" check-tsinfer --trees "$out" --sample-map "$OUT_DIR/chr$c.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr$c.arg_qc.json" --input-vcf "$in" --chr "$c"
  python3 "$HERE/arg.py" cache --file "$out" --file "$OUT_DIR/chr$c.sample_map.tsv" --file "$OUT_DIR/chr$c.arg_qc.json" > "$out.outputs.json.next"
  mv -f "$out.outputs.json.next" "$out.outputs.json"
  mv -f "$expected" "$meta"
}
check_one(){
  local c=$1 tree=$OUT_DIR/chr$1.trees map=$OUT_DIR/chr$1.sample_map.tsv qc=$OUT_DIR/chr$1.arg_qc.json
  for f in "$tree" "$map" "$qc"; do [[ -s $f ]] || { echo "ERROR: missing GU TRACE ARG output for chr$c: $f" >&2; return 1; }; done
  if [[ $c == X ]]; then
    python3 "$GU_F/arg.py" check-tsinfer --trees "$tree" --sample-map "$map" --qc-json "$qc" --input-vcf "$(prepared_vcf "$c")" --chr "$c"
    return
  fi
  python3 - "$tree" "$map" "$qc" "$c" <<'PY'
import csv,json,sys,tskit
tree,sample_map,qc,chrom=sys.argv[1:]
ts=tskit.load(tree)
with open(sample_map,newline="") as f: n=sum(1 for _ in csv.DictReader(f,delimiter="\t"))
with open(qc) as f: json.load(f)
if min(ts.num_trees,ts.num_sites,ts.num_samples,n)<1: raise SystemExit(f"ERROR: empty GU TRACE ARG for chr{chrom}")
if n != ts.num_samples: raise SystemExit(f"ERROR: sample map rows {n} != tree samples {ts.num_samples} for chr{chrom}")
print(f"tsinfer CHECK PASS chr{chrom} trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites}")
PY
}
trace_one(){
  local c=$1
  local -a cmd=(python3 "$HERE/arg.py" trace --method tsinfer
    --tree "$OUT_DIR/chr$c.trees" --sample-map "$OUT_DIR/chr$c.sample_map.tsv"
    --out-dir "$TRACE_DIR" --chr "$c")
  [[ $ACTION != check ]] || cmd+=(--check)
  [[ $REPLACE == FALSE ]] || cmd+=(--replace)
  "${cmd[@]}"
}

trap 'refgen_stop_jobs; exit 130' INT
trap 'refgen_stop_jobs; exit 143' TERM
trap 'refgen_stop_jobs; exit 129' HUP
build_chr(){ prep_one "$1"; build_one "$1"; }
read -r -a CHR_ARRAY <<< "$CHRS"

case "$ACTION" in
  prepare) run_parallel prep_one "${CHR_ARRAY[@]}";;
  build)
    run_parallel build_chr "${CHR_ARRAY[@]}"
    cat > "$OUT_DIR/ARG_TSINFER_BUILD.txt" <<META
created=$(date -Is)
producer=/mnt/d/scripts/gu/arg.sh build
method=tsinfer
build=GRCh${REFGEN_GRCH:-37}
vcf_dir=$IN_DIR
chromosomes=$CHRS
META
    ;;
  check) run_parallel check_one "${CHR_ARRAY[@]}";;
esac
if [[ $FORMAT == trace && $ACTION != prepare ]]; then
  run_parallel trace_one "${CHR_ARRAY[@]}"
  if [[ $ACTION != check ]]; then
    printf 'method\ttsinfer\nformat\ttrace\nbuild\tGRCh%s\nchromosomes\t%s\n' \
      "${REFGEN_GRCH:-37}" "$CHRS" > "$TRACE_DIR/ARG_TRACE_BUILD.tsv"
  fi
fi
echo "tsinfer ARG data preparation $ACTION complete: $OUT_DIR"

)

arg_prep_vcf() (
# Prepare phased VCF input for ARG inference.
set -euo pipefail

usage() {
	cat <<'USAGE'
Usage:
  arg.sh --in chrN.vcf.gz --out chrN.vcf.gz [options]

Options:
  --threads N           bcftools threads (default: 4)
  --missing-max X       require F_MISSING < X (default: 0.05)
  --mac-min N           require minor allele count >= N (default: 2)
  --tmp-root DIR        temporary work root (default: /tmp/gu-arg-vcf-$USER)
  --ref-fasta FILE      optional reference FASTA for REF validation/normalization
  --replace TRUE|FALSE  replace an existing output (default: FALSE)

The input must contain INFO/AA and phased GT data. INFO/AA values such as
"A|||" are normalized to "A". Sites whose AA is unknown or does not equal
REF or ALT are removed.
USAGE
}

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
ensure_fill_tags() {
	if bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then return 0; fi
	if env -u BCFTOOLS_PLUGINS bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then
		unset BCFTOOLS_PLUGINS
		return 0
	fi
	return 1
}

in="" out="" threads=4 missing_max=0.05 mac_min=2
tmp_root=${ARG_TMP_ROOT:-/tmp/gu-arg-vcf-${USER:-user}}
ref_fasta="" replace=0 chrom="" build=37 x_psam=""
while (( $# )); do
	case "$1" in
		--in) [[ $# -ge 2 ]] || die "--in requires a file"; in=$2; shift 2 ;;
		--out) [[ $# -ge 2 ]] || die "--out requires a file"; out=$2; shift 2 ;;
		--threads) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || die "--threads requires a positive integer"; threads=$2; shift 2 ;;
		--missing-max) [[ $# -ge 2 ]] || die "--missing-max requires a number"; missing_max=$2; shift 2 ;;
		--mac-min) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || die "--mac-min requires a positive integer"; mac_min=$2; shift 2 ;;
		--tmp-root) [[ $# -ge 2 ]] || die "--tmp-root requires a directory"; tmp_root=$2; shift 2 ;;
		--ref-fasta) [[ $# -ge 2 ]] || die "--ref-fasta requires a file"; ref_fasta=$2; shift 2 ;;
		--chr) chrom=$2; shift 2 ;;
		--build) build=$2; shift 2 ;;
		--x-psam) x_psam=$2; shift 2 ;;
		--replace)
			[[ $# -ge 2 ]] || die "--replace requires TRUE or FALSE"
			case "${2,,}" in true|1|yes) replace=1 ;; false|0|no) replace=0 ;; *) die "--replace requires TRUE or FALSE" ;; esac
			shift 2 ;;
		-h|--help) usage; exit 0 ;;
		*) die "unknown option: $1" ;;
	esac
done

[[ -n $in && -n $out ]] || { usage >&2; exit 1; }
[[ $out == *.vcf.gz ]] || die "output must end in .vcf.gz: $out"
if [[ -s $out && -s $out.tbi && $replace -eq 0 ]]; then
	echo "SKIP: prepared ARG VCF already exists: $out"
	exit 0
fi
[[ -s $in ]] || die "input VCF is missing or empty: $in"
[[ -z $ref_fasta || -s $ref_fasta ]] || die "reference FASTA is missing: $ref_fasta"
for cmd in bcftools bgzip tabix awk sed stat; do need "$cmd"; done
ensure_fill_tags || die "bcftools fill-tags plugin is unavailable"
ARG_BCFTOOLS_BIN=$(type -P bcftools)
bcftools() {
	"$ARG_BCFTOOLS_BIN" "$@" 2> >(sed '/^Warning: cannot determine contig names given the \.csi index alone$/d' >&2)
}
bcftools view -h "$in" | grep '^##INFO=<ID=AA,' >/dev/null || die "input VCF has no INFO/AA header: $in"
python3 - "$missing_max" <<'PY' || die "--missing-max must be in [0,1)"
import sys
x=float(sys.argv[1]); assert 0 <= x < 1
PY

mkdir -p "$(dirname "$out")" "$tmp_root"
meta=$out.input.meta.tsv
completion_record() {
	printf 'input\t%s\nmissing_max\t%s\nmac_min\t%s\nref_fasta\t%s\n' \
		"$(stat -c '%n:%s:%Y' "$in")" "$missing_max" "$mac_min" \
		"$([[ -n $ref_fasta ]] && stat -c '%n:%s:%Y' "$ref_fasta" || printf none)"
}
lock=$out.lock
exec 9>"$lock"
if command -v flock >/dev/null 2>&1; then flock -n 9 || die "another process is preparing $out"; fi
work=$(mktemp -d "$tmp_root/argvcf.XXXXXX")
cleanup() { rm -rf -- "$work"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

raw_bcf=$work/pass.bcf
tagged_bcf=$work/tagged.bcf
filtered_bcf=$work/filtered.bcf
aa_tsv=$work/aa.tsv
aa_tsv_gz=$work/aa.tsv.gz
aa_hdr=$work/aa.hdr
info_only_bcf=$work/info-only.bcf
# Keep the standard .vcf.gz suffix so htslib/bcftools can discover the
# temporary file's .tbi index during validation.
out_part=${out%.vcf.gz}.part.$$.vcf.gz
rm -f "$out_part" "$out_part.tbi"

echo "[ARG-VCF] PASS + biallelic SNP selection: $in"
source_vcf=$in
view_args=(view --threads "$threads" -f PASS -m2 -M2 -v snps -Ob -o "$raw_bcf")
if [[ $chrom == X ]]; then
	# Keep filtering in bcftools; validate biological ploidy after variant QC,
	# avoiding a Python pass over billions of raw, mostly rare genotypes.
	x_contig=$(bcftools index -s "$in" | awk '$1=="X" || $1=="chrX" || $1=="23" {print $1}')
	[[ -n $x_contig && $x_contig != *$'\n'* ]] || die 'X input must have exactly one X contig'
	x_lo=2699521; x_hi=154931043
	if [[ $build == 38 ]]; then x_lo=2781480; x_hi=155701382; fi
	view_args+=(-t "$x_contig:$x_lo-$x_hi")
fi
if [[ -n $ref_fasta ]]; then
	# SNP-only data need no left alignment, but this excludes REF mismatches and
	# removes duplicate records at the same position before downstream filtering.
	bcftools "${view_args[@]}" "$source_vcf"
	bcftools index -f "$raw_bcf"
	bcftools norm --threads "$threads" -f "$ref_fasta" -c x -d all -Ob -o "$work/norm.bcf" "$raw_bcf"
	mv "$work/norm.bcf" "$raw_bcf"
	bcftools index -f "$raw_bcf"
else
	bcftools "${view_args[@]}" "$source_vcf"
	bcftools index -f "$raw_bcf"
	bcftools norm --threads "$threads" -d all -Ob -o "$work/dedup.bcf" "$raw_bcf"
	mv "$work/dedup.bcf" "$raw_bcf"
	bcftools index -f "$raw_bcf"
fi

echo "[ARG-VCF] Recalculate allele counts and missingness from GT"
# Some bcftools fill-tags versions do not provide MAC. For the biallelic sites
# selected above, MAC >= N is equivalent to AC >= N and AN - AC >= N.
bcftools +fill-tags "$raw_bcf" --threads "$threads" -Ob -o "$tagged_bcf" -- -t AC,AN,F_MISSING
bcftools index -f "$tagged_bcf"
bcftools view --threads "$threads" \
	-i "INFO/F_MISSING<=${missing_max} && INFO/AC>=${mac_min} && INFO/AN-INFO/AC>=${mac_min} && INFO/AA!='.'" \
	-Ob -o "$filtered_bcf" "$tagged_bcf"
bcftools index -f "$filtered_bcf"

echo "[ARG-VCF] Normalize INFO/AA and retain only AA matching REF or ALT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' "$filtered_bcf" | \
awk 'BEGIN{FS=OFS="\t"}
{
	aa=toupper($5); sub(/\|.*/,"",aa); sub(/,.*/,"",aa)
	ref=toupper($3); alt=toupper($4)
	if (aa ~ /^[ACGT]$/ && (aa==ref || aa==alt)) print $1,$2,$3,$4,aa
}' > "$aa_tsv"
[[ -s $aa_tsv ]] || die "no variants retained after INFO/AA normalization: $in"
bgzip -c "$aa_tsv" > "$aa_tsv_gz"
tabix -f -s 1 -b 2 -e 2 "$aa_tsv_gz"
printf '##INFO=<ID=AA,Number=1,Type=String,Description="Normalized ancestral allele parsed from original INFO/AA">\n' > "$aa_hdr"

# Apply one documented complement operation per annotation family. This avoids
# relying on ambiguous parsing of multiple caret-prefixed families in one list.
bcftools annotate --threads "$threads" -x INFO/AA "$filtered_bcf" -Ou | \
	bcftools annotate --threads "$threads" -a "$aa_tsv_gz" -h "$aa_hdr" \
		-c CHROM,POS,REF,ALT,INFO/AA -Ou | \
	bcftools view --threads "$threads" -i 'INFO/AA!="."' -Ou | \
	bcftools annotate --threads "$threads" -x '^INFO/AA' -Ob -o "$info_only_bcf"
format_remove=$(bcftools view -h "$info_only_bcf" | sed -n 's/^##FORMAT=<ID=\([^,>]*\).*/FORMAT\/\1/p' | grep -Fvx FORMAT/GT | paste -sd, - || true)
if [[ -n $format_remove ]]; then
	bcftools annotate --threads "$threads" -x "$format_remove" -Oz -o "$out_part" "$info_only_bcf"
else
	bcftools view --threads "$threads" -Oz -o "$out_part" "$info_only_bcf"
fi
tabix -f -p vcf "$out_part"

[[ $(bcftools query -l "$out_part" | wc -l) -gt 0 ]] || die "output has no samples: $out_part"
final_n=$(bcftools index -n "$out_part")
[[ $final_n -gt 0 ]] || die "output has no variants: $out_part"
unexpected_info=$(bcftools view -h "$out_part" | sed -n 's/^##INFO=<ID=\([^,>]*\).*/\1/p' | grep -Fvx AA || true)
unexpected_format=$(bcftools view -h "$out_part" | sed -n 's/^##FORMAT=<ID=\([^,>]*\).*/\1/p' | grep -Fvx GT || true)
[[ -z $unexpected_info ]] || die "output retained unexpected INFO annotations: $(tr '\n' ',' <<< "$unexpected_info")"
[[ -z $unexpected_format ]] || die "output retained unexpected FORMAT annotations: $(tr '\n' ',' <<< "$unexpected_format")"
bad_aa=$(bcftools query -f '%REF\t%ALT\t%INFO/AA\n' "$out_part" | \
	awk 'BEGIN{n=0} $3!=$1 && $3!=$2{n++} END{print n+0}')
[[ $bad_aa -eq 0 ]] || die "output contains $bad_aa sites with AA not matching REF/ALT"

validation_args=(--chr "$chrom" --build "$build")
[[ -z $x_psam ]] || validation_args+=(--x-psam "$x_psam")
python3 "$HERE/arg.py" validate-vcf "$out_part" "${validation_args[@]}"
mv -f "$out_part" "$out"
rm -f "$out.csi"
mv -f "$out_part.tbi" "$out.tbi"
completion_record > "$meta"
raw_n=$(bcftools index -n "$raw_bcf")
filtered_n=$(bcftools index -n "$filtered_bcf")
sample_n=$(bcftools query -l "$out" | wc -l)
printf 'metric\tvalue\ninput_pass_biallelic_snp\t%s\npost_mac_missing\t%s\npost_valid_AA\t%s\nsamples\t%s\nmissing_max\t%s\nmac_min\t%s\n' \
	"$raw_n" "$filtered_n" "$final_n" "$sample_n" "$missing_max" "$mac_min" > "${out%.vcf.gz}.qc.tsv"
echo "[ARG-VCF] DONE: $out variants=$final_n samples=$sample_n"

)

arg_trace_check() (
set -euo pipefail

# ARG construction belongs to reference-data preparation. This private helper
# lets TRACE validate ARG outputs before starting its analysis scheduler.
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
ARG_SH=${GU_ARG_SH:-$ROOT/arg.sh}
ACTION=${1:-check}
[[ $ACTION == check ]] || { echo "ERROR: internal ARG validation is check-only; build ARG data with /mnt/d/scripts/gu/arg.sh build" >&2; exit 2; }
[[ -s $ARG_SH ]] || { echo "ERROR: arg.sh is missing: $ARG_SH" >&2; exit 2; }

backend=${GU_ARG_METHOD:-${GU_ARG_BACKEND:-tsinfer}}
case "${backend,,}" in
  tsinfer) method=tsinfer; ref_action=check;;
  external) method=tsinfer; ref_action=check;;
  needle) method=needle; ref_action=check;;
  threads) method=threads; ref_action=check;;
  *) echo "ERROR: GU ARG method must be needle, tsinfer or threads" >&2; exit 2;;
esac
format=${GU_ARG_FORMAT:-native}
case "$format" in native|trace) ;; *) echo "ERROR: GU_ARG_FORMAT must be native or trace" >&2; exit 2;; esac
replace=FALSE
chrs=${GU_CHRS:-1-22,X}

args=(check --method "$method" --format "$format" --action "$ref_action"
      --dir-gen "${GU_TARGET_ROOT:?}" --arg-dir "${GU_ARG_DIR:-$GU_TARGET_ROOT/arg}"
      --chr "$chrs" --grch "${GU_BUILD:-37}" --replace "$replace"
      --threads "${GU_ARG_THREADS:-8}" --jobs "${GU_UNIT_JOBS:-1}"
      --sample-match-min-work "${GU_ARG_SAMPLE_MATCH_MIN_WORK:-50000000}")
[[ -z ${GU_TARGET_VCF_DIR:-} ]] || args+=(--dir-vcf "$GU_TARGET_VCF_DIR")
[[ -z ${GU_SAMPLE_PANEL:-} ]] || args+=(--sample-file "$GU_SAMPLE_PANEL")
if [[ $method == needle && -n ${GU_TARGET_GEN_PREFIX:-} ]]; then
  args+=(--dir-pfile "$(dirname -- "$GU_TARGET_GEN_PREFIX")")
fi
export REFGEN_RECOMB_RATE=${GU_RECOMB_RATE:-1e-8}
export REFGEN_MUT_RATE=${GU_MUT_RATE:-1.25e-8}
export REFGEN_ARG_ONE_BIT=${GU_ARG_ONE_BIT:-0}

echo "[GU CHECK] internal read-only ARG validation through:"
printf '  bash %q' "$ARG_SH"; printf ' %q' "${args[@]}"; printf '\n'
exec bash "$ARG_SH" "${args[@]}"

)

arg_logged() (
# Persist the whole public invocation, including startup errors and worker output.
set -euo pipefail
export PYTHONUNBUFFERED=1
action=$1; shift
[[ $action != prep_gen ]] || export REFGEN_GEN4ARG_ONLY=1
dir_gen= arg_dir= grch=37 method=needle
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case "${args[i]}" in
    --dir-gen) dir_gen=${args[i+1]:-};;
    --arg-dir) arg_dir=${args[i+1]:-};;
    --grch) grch=${args[i+1]:-37};;
    --method) method=${args[i+1]:-needle};;
  esac
done
case "${grch,,}" in 38|b38|grch38|hg38) grch=38;; *) grch=37;; esac
dir_gen=${dir_gen:-/mnt/d/data.BIG/refGen/1kg/$grch}
dir_gen=${dir_gen%/}
case "${dir_gen##*/}" in pfile|hap|typ|vcf|imp|gen) gen_root=$(dirname -- "$dir_gen");; *) gen_root=$dir_gen;; esac
if [[ -z $arg_dir ]]; then
  if [[ ${method,,} == threads ]]; then arg_dir=$gen_root/arg.threads; else arg_dir=$gen_root/arg; fi
fi
log_dir=$arg_dir/log
mkdir -p "$log_dir"
log_dir=$(cd -- "$log_dir" && pwd -P)
log=$(mktemp "$log_dir/arg.$action.$(date +%Y%m%dT%H%M%S).XXXXXX.log")
finish(){
  local status=$?
  printf '\n[ARG RUN END] time=%s exit_status=%s\n' "$(date -Is)" "$status" >> "$log"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
{
  printf '[ARG RUN START] time=%s pid=%s cwd=%q\n' "$(date -Is)" "$$" "$PWD"
  printf '[ARG COMMAND] ./arg.sh %q' "$action"
  printf ' %q' "$@"
  printf '\n'
} >> "$log"
printf 'Full run log: %s\n' "$log"
set +e
bash "$HERE/arg.sh" cli "$@" 2>&1 | (trap '' HUP; exec tee --output-error=warn-nopipe -a "$log")
statuses=("${PIPESTATUS[@]}")
status=${statuses[0]}
if ((status==0 && statuses[1]!=0)); then status=${statuses[1]}; fi
printf 'Run exit status: %s; log: %s\n' "$status" "$log"
exit "$status"

)

command=${1:---help}; shift || true
case "$command" in
  cli|needle|tsinfer|threads|prep_vcf|trace_check|logged) "arg_$command" "$@";;
  -h|--help) exec bash "$HERE/../arg.sh" --help;;
  *) echo "ERROR: unknown internal ARG stage: $command" >&2; exit 2;;
esac

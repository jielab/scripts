#!/usr/bin/env bash
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
    -h|--help) usage; exit 0;;
    *) die "unknown build option: $1";;
  esac
done

case "$method" in needle|tsinfer) ;; *) die "--method must be needle or tsinfer";; esac
case "$x_odd_male" in error|drop-last) ;; *) die '--x-odd-male must be error or drop-last';; esac
if [[ $method == needle && $missing_set == 0 ]]; then missing_max=0; fi
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
arg_dir=${arg_dir:-$gen_root/arg}
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
  if [[ $method == needle ]]; then gen4arg_dir=${dir_pfile}.4arg; else gen4arg_dir=${dir_vcf}.4arg; fi
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

if [[ $method == needle ]]; then
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
if [[ $features == auto ]]; then [[ $format == trace ]] && features=FALSE || features=TRUE; fi
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
bash "$HERE/arg.${method}.sh"
if [[ $action != prepare && $action != check ]]; then
  # Retain the input/output locks until validation finishes. A failed check
  # must propagate to the caller instead of reporting a successful build.
  export REFGEN_ARG_ACTION=check
  bash "$HERE/arg.${method}.sh"
  echo "BUILD PASS: output validation complete: $arg_dir"
fi

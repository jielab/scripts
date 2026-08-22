#!/bin/bash -l
# PRS-CSx -> 20-PC projection -> automatic ancestry -> PRS-CSx/DiscoDivas evaluation.

set -euo pipefail
export LC_ALL=C
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)

usage() {
  cat <<'USAGE'
Usage: ./grid.sh csx|pca|ancestry|eval|all [options]

Core options:
  --trait STR                phenotype [height]
  --jobs INT                 simultaneous chromosome jobs [12]
  --dir0 PATH                D-drive root under WSL [/mnt/d]
  --dir-gen PATH             target PLINK bfiles chr1..chr22 [/mnt/d/data/ukb/gen/typ]
  --dir-imp PATH             target PLINK 2 pfiles chr1..chr22 [/mnt/h/ukbGen/imp]
  --dir-gwas PATH            cleaned ancestry-specific GWAS files
  --outdir PATH              trait output [/mnt/d/analysis/grid/TRAIT]
  --phe-file PATH            RDS with eid, TRAIT, age, sex
  --n-gwas SPEC              AFR=...,EAS=...,EUR=...,SAS=...
  --phi NUM                  fixed PRS-CSx global shrinkage [1e-2]
  --mcmc-iter INT            PRS-CSx iterations [4000]
  --mcmc-burnin INT          PRS-CSx burn-in [2000]
  --mcmc-thin INT            PRS-CSx thinning [5]

PCA + ancestry options:
  --med-file PATH            released four-population PC medians
  --cov-pcs INT              projected PCs retained and used as regression covariates [20]
  --distance-pcs INT         PCs used for DiscoDivas distances [10]

  --auto-ancestry BOOL       build ancestry automatically at the end of PCA [TRUE]
  --auto-ancestry-file PATH  automatic output [D:/data/ukb/phe/pca/ukb.ancestry.auto.tsv.gz]
  --ethnicity-col STR        self-reported UKB ethnicity column in phe.rds [auto]
  --ancestry-prob-threshold NUM assign AFR/EAS/EUR/SAS only above this probability [0.90]
  --anchor-max-per-group INT balanced pure anchors per group for automatic LDA [10000]
  --anchor-margin-quantile NUM discard the least separated anchor fraction [0.25]

Optional external ancestry override:
  --ancestry-file PATH       custom TSV/TSV.GZ; overrides automatic ancestry
  --ancestry-id-col STR      ID column in custom ancestry file [eid]
  --ancestry-col STR         genetic ancestry column in custom file
  --ancestry-prob-col STR    assignment probability column in custom file
  --self-report-col STR      self-reported ancestry column in custom file
  --target-ids PATH          optional one-column analysis IDs
  --exclude-ids PATH         discovery-overlap/related IDs to exclude

Evaluation options:
  --prscsx-repeats INT       PRS-CSx 50:50 validation/testing repetitions [100]
  --prscsx-prob-min NUM      minimum ancestry probability for PRS-CSx analysis [0.9]
  --validation-frac NUM      validation fraction for PRS-CSx repeated splits [0.5]
  --disco-repeats INT        DiscoDivas fine-tuning/testing repetitions [20]
  --fine-tune-n INT          participants per AFR/EAS/EUR/SAS fine-tuning cohort [1300]
  --fine-tune-prob-min NUM   ancestry probability for DiscoDivas fine-tuning [0.999999]
  --require-self-match BOOL  require genetic/self-report ancestry agreement [TRUE]
  --min-test-group INT       minimum testing N reported per group [200]
  --disco-a STR              A values in AFR,EAS,EUR,SAS order [1,1,1,1]
  --run-prscsx BOOL          run PRS-CSx repeated-split analysis [TRUE]
  --run-disco BOOL           run DiscoDivas empirical-design analysis [TRUE]
  --write-predictions BOOL   save first-repeat participant audit data [FALSE]
  --force BOOL               recreate completed outputs [FALSE]

Compatibility aliases:
  --eval-repeats = --prscsx-repeats
  --eval-pcs = --cov-pcs
  --min-group = --min-test-group

Important:
  * PCA is trait-independent; `./grid.sh pca` is sufficient.
  * `grid.sh pca` runs PCA and then ancestry. A complete PCA cache is skipped, after
    which any missing or stale ancestry outputs are built automatically.
  * Core PCA/ancestry results and QC files are written under
    /mnt/d/data/ukb/phe/pca. Chromosome intermediates, commands, and logs are under
    /mnt/d/analysis/grid/pca.
  * Ancestry creates ukb.ancestry.auto.tsv.gz using projected PCs, released
    1KG centers, and an auto-detected UKB self-reported ethnicity column from phe.rds.
  * Mixed/uncertain participants are assigned OTH; AMR is not inferred from the
    released four-population centers unless a custom ancestry file is supplied.
  * There is no separate disco module. DiscoDivas is reconstructed inside eval.

Examples:
  ./grid.sh pca
  ./grid.sh eval --trait height
  ./grid.sh all --trait height --jobs 12
  ./grid.sh ancestry --force TRUE  # optional ancestry-only rebuild
USAGE
}

dir0=/mnt/d
dir_gen=/mnt/d/data/ukb/gen/typ
dir_imp=/mnt/h/ukbGen/imp
dir_gwas=""
outdir=""
trait=height
step=""
jobs=12
force=FALSE
phe_file=""
n_gwas_spec="AFR=300000,EAS=500000,EUR=4000000,SAS=80000"
phi="1e-2"
mcmc_iter=4000
mcmc_burnin=2000
mcmc_thin=5
med_file=""
disco_a="1,1,1,1"
cov_pcs=20
distance_pcs=10
ancestry_file=""
ancestry_id_col=eid
ancestry_col=""
ancestry_prob_col=""
self_report_col=""
auto_ancestry=TRUE
auto_ancestry_file=""
ethnicity_col=auto
ancestry_prob_threshold=0.90
anchor_max_per_group=10000
anchor_margin_quantile=0.25
target_ids=""
exclude_ids=""
prscsx_repeats=100
prscsx_prob_min=0.9
validation_frac=0.5
disco_repeats=20
fine_tune_n=1300
fine_tune_prob_min=0.999999
require_self_match=TRUE
min_test_group=200
run_prscsx=TRUE
run_disco=TRUE
write_predictions=FALSE

need_value() { [[ -n "${2-}" && "${2-}" != --* ]] || { echo "ERROR: missing value for $1" >&2; exit 2; }; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    csx|pca|ancestry|eval|all) [[ -z "$step" ]] || { echo "ERROR: multiple modules specified: $step and $1" >&2; exit 2; }; step=${1,,}; shift ;;
    --jobs) need_value "$1" "${2-}"; jobs=$2; shift 2 ;;
    --trait) need_value "$1" "${2-}"; trait=$2; shift 2 ;;
    --dir0) need_value "$1" "${2-}"; dir0=${2%/}; shift 2 ;;
    --dir-gen|--dir_gen) need_value "$1" "${2-}"; dir_gen=${2%/}; shift 2 ;;
    --dir-imp|--dir_imp) need_value "$1" "${2-}"; dir_imp=${2%/}; shift 2 ;;
    --dir-gwas|--dir_gwas) need_value "$1" "${2-}"; dir_gwas=${2%/}; shift 2 ;;
    --outdir) need_value "$1" "${2-}"; outdir=${2%/}; shift 2 ;;
    --phe-file|--phe_file) need_value "$1" "${2-}"; phe_file=$2; shift 2 ;;
    --n-gwas|--n_gwas) need_value "$1" "${2-}"; n_gwas_spec=$2; shift 2 ;;
    --phi) need_value "$1" "${2-}"; phi=$2; shift 2 ;;
    --mcmc-iter) need_value "$1" "${2-}"; mcmc_iter=$2; shift 2 ;;
    --mcmc-burnin) need_value "$1" "${2-}"; mcmc_burnin=$2; shift 2 ;;
    --mcmc-thin) need_value "$1" "${2-}"; mcmc_thin=$2; shift 2 ;;
    --med-file|--med_file) need_value "$1" "${2-}"; med_file=$2; shift 2 ;;
    --disco-a|--disco_a) need_value "$1" "${2-}"; disco_a=$2; shift 2 ;;
    --cov-pcs|--cov_pcs|--eval-pcs|--eval_pcs) need_value "$1" "${2-}"; cov_pcs=$2; shift 2 ;;
    --distance-pcs|--distance_pcs) need_value "$1" "${2-}"; distance_pcs=$2; shift 2 ;;
    --ancestry-file|--ancestry_file) need_value "$1" "${2-}"; ancestry_file=$2; shift 2 ;;
    --ancestry-id-col|--ancestry_id_col) need_value "$1" "${2-}"; ancestry_id_col=$2; shift 2 ;;
    --ancestry-col|--ancestry_col) need_value "$1" "${2-}"; ancestry_col=$2; shift 2 ;;
    --ancestry-prob-col|--ancestry_prob_col) need_value "$1" "${2-}"; ancestry_prob_col=$2; shift 2 ;;
    --self-report-col|--self_report_col) need_value "$1" "${2-}"; self_report_col=$2; shift 2 ;;
    --auto-ancestry|--auto_ancestry) need_value "$1" "${2-}"; auto_ancestry=${2^^}; shift 2 ;;
    --auto-ancestry-file|--auto_ancestry_file) need_value "$1" "${2-}"; auto_ancestry_file=$2; shift 2 ;;
    --ethnicity-col|--ethnicity_col) need_value "$1" "${2-}"; ethnicity_col=$2; shift 2 ;;
    --ancestry-prob-threshold|--ancestry_prob_threshold) need_value "$1" "${2-}"; ancestry_prob_threshold=$2; shift 2 ;;
    --anchor-max-per-group|--anchor_max_per_group) need_value "$1" "${2-}"; anchor_max_per_group=$2; shift 2 ;;
    --anchor-margin-quantile|--anchor_margin_quantile) need_value "$1" "${2-}"; anchor_margin_quantile=$2; shift 2 ;;
    --target-ids|--target_ids) need_value "$1" "${2-}"; target_ids=$2; shift 2 ;;
    --exclude-ids|--exclude_ids) need_value "$1" "${2-}"; exclude_ids=$2; shift 2 ;;
    --prscsx-repeats|--prscsx_repeats|--eval-repeats|--eval_repeats) need_value "$1" "${2-}"; prscsx_repeats=$2; shift 2 ;;
    --prscsx-prob-min|--prscsx_prob_min) need_value "$1" "${2-}"; prscsx_prob_min=$2; shift 2 ;;
    --validation-frac|--validation_frac) need_value "$1" "${2-}"; validation_frac=$2; shift 2 ;;
    --disco-repeats|--disco_repeats) need_value "$1" "${2-}"; disco_repeats=$2; shift 2 ;;
    --fine-tune-n|--fine_tune_n) need_value "$1" "${2-}"; fine_tune_n=$2; shift 2 ;;
    --fine-tune-prob-min|--fine_tune_prob_min) need_value "$1" "${2-}"; fine_tune_prob_min=$2; shift 2 ;;
    --require-self-match|--require_self_match) need_value "$1" "${2-}"; require_self_match=${2^^}; shift 2 ;;
    --min-test-group|--min_test_group|--min-group|--min_group) need_value "$1" "${2-}"; min_test_group=$2; shift 2 ;;
    --run-prscsx|--run_prscsx) need_value "$1" "${2-}"; run_prscsx=${2^^}; shift 2 ;;
    --run-disco|--run_disco) need_value "$1" "${2-}"; run_disco=${2^^}; shift 2 ;;
    --write-predictions|--write_predictions) need_value "$1" "${2-}"; write_predictions=${2^^}; shift 2 ;;
    --force) need_value "$1" "${2-}"; force=${2^^}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$step" ]] || { echo "ERROR: a module is required: csx, pca, ancestry, eval, or all" >&2; exit 2; }
case "$step" in csx|pca|ancestry|eval|all) ;; *) echo "ERROR: invalid module $step" >&2; exit 2 ;; esac
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs must be positive" >&2; exit 2; }
[[ "$prscsx_repeats" =~ ^[1-9][0-9]*$ && "$disco_repeats" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: repeat counts must be positive integers" >&2; exit 2; }
[[ "$fine_tune_n" =~ ^[1-9][0-9]*$ && "$cov_pcs" =~ ^[1-9][0-9]*$ && "$distance_pcs" =~ ^[1-9][0-9]*$ && "$anchor_max_per_group" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid fine-tuning/PCA/ancestry settings" >&2; exit 2; }
awk -v x="$ancestry_prob_threshold" 'BEGIN { exit !(x > 0.5 && x < 1) }' || { echo "ERROR: --ancestry-prob-threshold must be in (0.5,1)" >&2; exit 2; }
awk -v x="$anchor_margin_quantile" 'BEGIN { exit !(x >= 0 && x < 0.9) }' || { echo "ERROR: --anchor-margin-quantile must be in [0,0.9)" >&2; exit 2; }
(( distance_pcs <= cov_pcs )) || { echo "ERROR: --distance-pcs cannot exceed --cov-pcs" >&2; exit 2; }
awk -v x="$validation_frac" 'BEGIN { exit !(x > 0.1 && x < 0.9) }' || { echo "ERROR: --validation-frac must be between 0.1 and 0.9" >&2; exit 2; }
for v in force require_self_match run_prscsx run_disco write_predictions auto_ancestry; do
  eval "x=\${$v}"; case "$x" in TRUE|FALSE) ;; *) echo "ERROR: --${v//_/-} must be TRUE/FALSE" >&2; exit 2 ;; esac
done
[[ -z "$ancestry_file" || -n "$ancestry_col" ]] || { echo "ERROR: --ancestry-file requires --ancestry-col" >&2; exit 2; }

[[ -n "$dir_gwas" ]] || dir_gwas="$dir0/data.BIG/gwas/main/clean"
[[ -n "$outdir" ]] || outdir="$dir0/analysis/grid/$trait"

grid_privacy_root="$dir0/analysis/grid"
privacy_cleanup_completed_run=FALSE
cleanup_ukb_individual_outputs() {
  local status=$? f
  trap - EXIT
  if (( status != 0 )) || [[ "$privacy_cleanup_completed_run" != TRUE ]]; then
    echo "GRID did not reach full completion (exit=$status); preserving participant-level score outputs for resume." >&2
    exit "$status"
  fi
  case "$grid_privacy_root" in ""|/|/mnt|/mnt/|/mnt/d|/mnt/d/|/mnt/d/analysis|/mnt/d/analysis/)
    echo "PRIVACY ERROR: unsafe grid output root: '$grid_privacy_root'" >&2; exit 70 ;;
  esac
  if [[ -d "$grid_privacy_root" ]]; then
    echo "Privacy cleanup: removing GRID UKB participant-level score outputs"
    while IFS= read -r -d '' f; do rm -f -- "$f"; done < <(find "$grid_privacy_root" -type f \( \
      -name '*.sscore' -o -name 'first_repeat_predictions.tsv.gz' \) -print0)
  fi
  exit "$status"
}
trap cleanup_ukb_individual_outputs EXIT
[[ -n "$phe_file" ]] || phe_file="$dir0/data/ukb/phe/Rdata/phe.rds"
[[ -d "$dir0/software/bin" ]] && export PATH="$dir0/software/bin:$PATH"
dir_csx="$outdir/csx"
dir_pca_cache="$dir0/data/ukb/phe/pca"
dir_pca="$dir_pca_cache"
dir_pca_raw="$dir0/analysis/grid/pca"
pca_file="$dir_pca_cache/ukb.discodivas.pca.tsv.gz"
[[ -n "$auto_ancestry_file" ]] || auto_ancestry_file="$dir_pca_cache/ukb.ancestry.auto.tsv.gz"
dir_eval="$outdir/evaluate"
dir_prscsx="$SCRIPT_DIR/csx"
dir_discodivas="$SCRIPT_DIR/DiscoDivas"
code_dir="$SCRIPT_DIR/f"
[[ -d "$code_dir" ]] || code_dir="$SCRIPT_DIR"
combine_csx_r="$code_dir/combine_csx_prs.R"
prepare_sumstats_r="$code_dir/prepare_prscsx_sumstats.R"
combine_pca_r="$code_dir/combine_disco_pca.R"
evaluate_r="$code_dir/evaluate_prs.R"
prepare_ancestry_auto_r="$code_dir/prepare_ancestry_auto.R"
[[ -n "$med_file" ]] || med_file="$dir_discodivas/files/med.g1000.4pop.tsv"
races=(AFR EAS EUR SAS)
model_id=$(printf 'phi_%s' "$phi" | sed 's/[^A-Za-z0-9]/_/g')
case "$step" in
  pca|ancestry) audit_file="$dir_pca_raw/pipeline_runtime.tsv" ;;
  *) audit_file="$outdir/pipeline_runtime.tsv" ;;
esac
privacy_cleanup_completed_run=TRUE
mkdir -p "$(dirname "$audit_file")"
[[ -e "$audit_file" ]] || printf 'step\tstarted\tfinished\tseconds\tstatus\tnote\n' > "$audit_file"
log() { echo "[$(date '+%F %T')] $*" >&2; }
file_stamp() { [[ -e "$1" ]] && stat -c '%s:%Y' "$1" || printf 'MISSING'; }
run_timed_step() {
  local name=$1; shift
  local start_epoch end_epoch started finished rc
  start_epoch=$(date +%s); started=$(date '+%F %T')
  log "START step: $name"
  if "$@"; then rc=0; else rc=$?; fi
  end_epoch=$(date +%s); finished=$(date '+%F %T')
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$started" "$finished" "$((end_epoch-start_epoch))" "$([[ $rc -eq 0 ]] && echo OK || echo FAILED)" "step=$step force=$force" >> "$audit_file"
  [[ $rc -eq 0 ]] || return "$rc"
  log "DONE step: $name in $((end_epoch-start_epoch)) seconds"
}
need_file() { [[ -s "$1" ]] || { echo "ERROR: missing/empty file: $1" >&2; exit 1; }; }
run_parallel_cmds() {
  local label=$1; shift
  local -a cmds=("$@")
  ((${#cmds[@]})) || { log "$label: nothing to run"; return; }
  log "$label: ${#cmds[@]} jobs with concurrency $jobs"
  printf '%s\0' "${cmds[@]}" | xargs -0 -r -P "$jobs" -n 1 bash -c
}
need_chr_bfiles() { local c e; for c in {1..22}; do for e in bed bim fam; do need_file "$dir_gen/chr$c.$e"; done; done; }
need_chr_pfiles() { local c e; for c in {1..22}; do for e in pgen pvar psam; do need_file "$dir_imp/chr$c.$e"; done; done; }

parse_n_gwas() {
  declare -gA N_GWAS=()
  local item key val
  IFS=',' read -r -a items <<< "$n_gwas_spec"
  for item in "${items[@]}"; do
    key=${item%%=*}; val=${item#*=}; key=${key^^}
    [[ "$key" != "$item" && "$val" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid --n-gwas item: $item" >&2; exit 2; }
    N_GWAS[$key]=$val
  done
  for key in "${races[@]}"; do [[ -n "${N_GWAS[$key]-}" ]] || { echo "ERROR: missing N for $key" >&2; exit 2; }; done
}

step_csx() {
  parse_n_gwas
  need_file "$dir_prscsx/PRScsx.py"; need_file "$prepare_sumstats_r"; need_file "$combine_csx_r"
  command -v Rscript >/dev/null || { echo "ERROR: Rscript not found" >&2; exit 1; }
  command -v plink2 >/dev/null || { echo "ERROR: plink2 not found" >&2; exit 1; }
  local python_bin; python_bin=$(command -v python3 || command -v python || true); [[ -n "$python_bin" ]] || { echo "ERROR: python not found" >&2; exit 1; }
  need_chr_bfiles
  mkdir -p "$dir_csx/input/log" "$dir_csx/$model_id/log"
  local race chr sst_str n_str race_str model_dir="$dir_csx/$model_id"
  race_str=$(IFS=,; echo "${races[*]}")
  n_str=""
  for race in "${races[@]}"; do
    n_str+="${n_str:+,}${N_GWAS[$race]}"
    need_file "$dir_gwas/$trait.$race.gz"
    mkdir -p "$dir_csx/input/$race" "$model_dir/$race"
  done
  log "PRS-CSx: phi=$phi; n_gwas=$n_str; populations=$race_str"

  # Cache safety: posterior files are reusable only when all model-defining inputs match.
  local config_file="$model_dir/model_config.tsv" config_tmp="$model_dir/model_config.current.tsv"
  {
    printf 'key\tvalue\n'
    printf 'trait\t%s\nphi\t%s\nn_gwas\t%s\npopulations\t%s\n' "$trait" "$phi" "$n_gwas_spec" "$race_str"
    printf 'n_iter\t%s\nburnin\t%s\nthin\t%s\n' "$mcmc_iter" "$mcmc_burnin" "$mcmc_thin"
    printf 'dir_gen\t%s\nref_dir\t%s\n' "$dir_gen" "$dir0/data.BIG/refLD/csx"
    printf 'prscsx_py\t%s|%s\n' "$dir_prscsx/PRScsx.py" "$(file_stamp "$dir_prscsx/PRScsx.py")"
    printf 'prepare_sumstats\t%s|%s\n' "$prepare_sumstats_r" "$(file_stamp "$prepare_sumstats_r")"
    for chr in {1..22}; do printf 'target_bim_chr%s\t%s|%s\n' "$chr" "$dir_gen/chr$chr.bim" "$(file_stamp "$dir_gen/chr$chr.bim")"; done
    for race in "${races[@]}"; do
      printf 'gwas_%s\t%s|%s\n' "$race" "$dir_gwas/$trait.$race.gz" "$(file_stamp "$dir_gwas/$trait.$race.gz")"
    done
  } > "$config_tmp"
  local legacy_cache=FALSE
  [[ -e "$model_dir/chr1.done" || -s "$model_dir/AFR/AFR.chr1.pst_eff.txt" ]] && legacy_cache=TRUE
  if [[ -s "$config_file" ]] && ! cmp -s "$config_file" "$config_tmp" && [[ "$force" != TRUE ]]; then
    echo "ERROR: PRS-CSx cached outputs were created with a different configuration. Re-run once with --force TRUE." >&2
    diff -u "$config_file" "$config_tmp" >&2 || true
    exit 1
  fi
  if [[ ! -s "$config_file" && "$legacy_cache" == TRUE && "$force" != TRUE ]]; then
    echo "ERROR: legacy PRS-CSx cache has no model_config.tsv, so its n_gwas/MCMC inputs cannot be verified. Re-run once with --force TRUE." >&2
    exit 1
  fi

  # Stream each ancestry-specific GWAS once and split it into chr1-chr22 files.
  # This avoids reading a large GWAS 22 separate times.
  local -a prep_cmds=()
  for race in "${races[@]}"; do
    local prep_complete=TRUE
    [[ -s "$dir_csx/input/$race/$race.sumstats_qc.tsv" ]] || prep_complete=FALSE
    for chr in {1..22}; do
      [[ -s "$dir_csx/input/$race/$race.chr$chr.dat" ]] || prep_complete=FALSE
    done
    if [[ "$force" == TRUE || "$prep_complete" != TRUE ]]; then
      prep_cmds+=("Rscript '$prepare_sumstats_r' --in '$dir_gwas/$trait.$race.gz' --out-dir '$dir_csx/input/$race' --pop '$race' > '$dir_csx/input/log/$race.prepare.log' 2>&1")
    fi
  done
  run_parallel_cmds "Prepare PRS-CSx summary statistics" "${prep_cmds[@]}"

  local -a model_cmds=()
  for chr in {1..22}; do
    sst_str=""
    for race in "${races[@]}"; do sst_str+="${sst_str:+,}$dir_csx/input/$race/$race.chr$chr.dat"; done
    local cmd="$model_dir/chr$chr.cmd"
    cat > "$cmd" <<EOF2
#!/bin/bash -l
set -euo pipefail
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1
exec > >(tee "$model_dir/log/chr$chr.log") 2>&1
complete=TRUE
for race in ${races[*]}; do
  w="$model_dir/\$race/\$race.chr$chr.pst_eff.txt"
  [[ -s "\$w" ]] || complete=FALSE
done
if [[ "$force" == TRUE || "\$complete" != TRUE ]]; then
  "$python_bin" "$dir_prscsx/PRScsx.py" --ref_dir="$dir0/data.BIG/refLD/csx" \
    --bim_prefix="$dir_gen/chr$chr" --chrom="$chr" --pop="$race_str" \
    --sst_file="$sst_str" --n_gwas="$n_str" --phi="$phi" \
    --n_iter="$mcmc_iter" --n_burnin="$mcmc_burnin" --thin="$mcmc_thin" \
    --out_dir="$model_dir" --out_name="chr$chr" --meta=FALSE --seed=$((20260801 + chr))
fi
for race in ${races[*]}; do
  w="$model_dir/\$race/\$race.chr$chr.pst_eff.txt"
  [[ -s "\$w" ]] || { echo "ERROR: missing \$w" >&2; exit 1; }
  awk 'NF != 6 {exit 1} END {if (NR < 100) exit 1}' "\$w" || { echo "ERROR: malformed \$w" >&2; exit 1; }
  out="$model_dir/\$race/\$race.chr$chr"
  if [[ "$force" == TRUE || ! -s "\$out.sscore" || ! -e "\$out.score.done" ]]; then
    plink2 --bfile "$dir_gen/chr$chr" --score "\$w" 2 4 6 cols=nallele,scoresums list-variants --out "\$out"
    touch "\$out.score.done"
  fi
done
touch "$model_dir/chr$chr.done"
EOF2
    chmod +x "$cmd"
    local complete=TRUE
    [[ -e "$model_dir/chr$chr.done" ]] || complete=FALSE
    for race in "${races[@]}"; do [[ -s "$model_dir/$race/$race.chr$chr.sscore" && -e "$model_dir/$race/$race.chr$chr.score.done" ]] || complete=FALSE; done
    [[ "$force" == TRUE || "$complete" != TRUE ]] && model_cmds+=("$cmd")
  done
  log "PRS-CSx chromosomes requiring work: ${#model_cmds[@]}/22"
  run_parallel_cmds "PRS-CSx MCMC and scoring" "${model_cmds[@]}"

  local chr_qc="$model_dir/chromosome_qc.tsv"
  printf 'chr\tpopulation\tweight_rows\tweight_bytes\tscore_rows\tscore_bytes\tweight_modified\tscore_modified\n' > "$chr_qc"
  for chr in {1..22}; do
    for race in "${races[@]}"; do
      local w="$model_dir/$race/$race.chr$chr.pst_eff.txt" sc="$model_dir/$race/$race.chr$chr.sscore"
      need_file "$w"; need_file "$sc"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$chr" "$race" "$(wc -l < "$w")" "$(stat -c %s "$w")" "$(wc -l < "$sc")" "$(stat -c %s "$sc")" "$(stat -c %y "$w" | cut -d. -f1)" "$(stat -c %y "$sc" | cut -d. -f1)" >> "$chr_qc"
    done
  done

  rm -f "$dir_csx/score_manifest.tsv"
  Rscript "$combine_csx_r" --dir "$model_dir" --pops "$race_str" --model-id "$model_id" --phi "$phi" \
    --manifest "$dir_csx/score_manifest.tsv" --qc "$model_dir/score_qc.tsv"
  need_file "$dir_csx/score_manifest.tsv"
  mv -f "$config_tmp" "$config_file"
  log "PRS-CSx complete: $dir_csx/score_manifest.tsv"
}


pca_cache_has_required_pcs() {
  [[ -s "$pca_file" && -s "$dir_pca/ukb_reference_distances.tsv.gz" \
     && -s "$dir_pca/pca_qc.xlsx" && -s "$dir_pca/Fig1.PCA_QC.png" ]] || return 1
  local header
  if [[ "$pca_file" == *.gz ]]; then
    gzip -t "$pca_file" 2>/dev/null || return 1
    header=$(gzip -cd "$pca_file" 2>/dev/null | head -n 1 || true)
  else
    header=$(head -n 1 "$pca_file" || true)
  fi
  [[ $'\t'"$header"$'\t' == *$'\tPC'"$cov_pcs"$'\t'* ]]
}

step_ancestry() {
  [[ "$auto_ancestry" == TRUE ]] || { log "Automatic ancestry construction is disabled"; return; }
  need_file "$phe_file"; need_file "$pca_file"; need_file "$med_file"; need_file "$prepare_ancestry_auto_r"
  mkdir -p "$dir_pca" "$dir_pca_cache" "$dir_pca_raw"
  local config="$auto_ancestry_file.config.tsv" current="$dir_pca_raw/ancestry.config.current.tsv"
  {
    printf 'key\tvalue\n'
    printf 'phe\t%s|%s\n' "$phe_file" "$(file_stamp "$phe_file")"
    printf 'pca\t%s|%s\n' "$pca_file" "$(file_stamp "$pca_file")"
    printf 'med\t%s|%s\n' "$med_file" "$(file_stamp "$med_file")"
    printf 'script\t%s|%s\n' "$prepare_ancestry_auto_r" "$(file_stamp "$prepare_ancestry_auto_r")"
    printf 'ethnicity_col\t%s\nprob_threshold\t%s\nfine_prob_threshold\t%s\n' "$ethnicity_col" "$ancestry_prob_threshold" "$fine_tune_prob_min"
    printf 'cov_pcs\t%s\ndistance_pcs\t%s\nanchor_max_per_group\t%s\nanchor_margin_quantile\t%s\n' "$cov_pcs" "$distance_pcs" "$anchor_max_per_group" "$anchor_margin_quantile"
  } > "$current"
  if [[ "$force" != TRUE && -s "$auto_ancestry_file" && -s "$config" \
        && -s "$dir_pca/ancestry_auto_qc.xlsx" && -s "$dir_pca/Fig2.Ancestry_QC.png" ]] \
        && cmp -s "$config" "$current"; then
    log "Automatic ancestry cache is current: $auto_ancestry_file"
    return
  fi
  Rscript "$prepare_ancestry_auto_r" --phe "$phe_file" --pca "$pca_file" --med "$med_file" \
    --out "$auto_ancestry_file" --outdir "$dir_pca" --ethnicity-col "$ethnicity_col" \
    --n-pc "$cov_pcs" --distance-pcs "$distance_pcs" --prob-threshold "$ancestry_prob_threshold" \
    --fine-prob-threshold "$fine_tune_prob_min" --anchor-max-per-group "$anchor_max_per_group" \
    --anchor-margin-quantile "$anchor_margin_quantile"
  need_file "$auto_ancestry_file"; need_file "$dir_pca/ancestry_auto_qc.xlsx"; need_file "$dir_pca/Fig2.Ancestry_QC.png"
  mv -f "$current" "$config"
  log "Automatic ancestry complete: $auto_ancestry_file"
}

step_pca() {
  local weight="$dir_discodivas/files/g1k_hm3_maf5_woamb_wolr.pca.weight"
  need_file "$weight"; need_file "$med_file"; need_file "$combine_pca_r"
  mkdir -p "$dir_pca_raw/log" "$dir_pca_cache"

  if [[ "$force" != TRUE ]] && pca_cache_has_required_pcs; then
    log "Complete PCA results and QC already exist; all PCA commands are skipped"
    step_ancestry
    return
  fi

  if [[ -s "$pca_file" && "$force" != TRUE ]]; then
    log "Existing PCA cache lacks PC$cov_pcs and will be rebuilt"
  fi
  command -v plink2 >/dev/null || { echo "ERROR: plink2 not found" >&2; exit 1; }
  need_chr_pfiles
  local available_mb pca_jobs plink_memory=4096
  available_mb=$(awk '/^MemAvailable:/ {print int($2/1024)}' /proc/meminfo)
  # Large imputed chromosomes need substantially more than 1 GiB for PLINK's
  # variant metadata and scoring workspace.  Cap PCA at four 4-GiB processes
  # so other programs retain predictable memory headroom; reduce it further
  # automatically when less memory is available.
  pca_jobs=$((available_mb * 80 / 100 / plink_memory))
  ((pca_jobs < 1)) && pca_jobs=1
  ((pca_jobs > 4)) && pca_jobs=4
  ((pca_jobs > jobs)) && pca_jobs=$jobs
  log "PCA memory plan: ${available_mb} MiB available, concurrency $pca_jobs, ${plink_memory} MiB per PLINK job"

  local -a cmds=(); local chr cmd
  for chr in {1..22}; do
    cmd="$dir_pca_raw/chr$chr.cmd"
    cat > "$cmd" <<EOF2
#!/bin/bash -l
set -euo pipefail
exec > >(tee "$dir_pca_raw/log/chr$chr.log") 2>&1
plink2 --pfile "$dir_imp/chr$chr" \\
  --score "$weight" 2 6 header-read no-mean-imputation list-variants cols=nallele,scoresums \\
  --score-col-nums 7-$((6 + cov_pcs)) --rm-dup exclude-mismatch --memory "$plink_memory" --threads 1 \\
  --out "$dir_pca_raw/chr$chr"
touch "$dir_pca_raw/chr$chr.pca${cov_pcs}.done"
EOF2
    chmod +x "$cmd"
    # A previous 10-PC cache must not be reused for the 20-PC projection.
    if [[ "$force" == TRUE || ! -s "$dir_pca_raw/chr$chr.sscore" || ! -s "$dir_pca_raw/chr$chr.sscore.vars" || ! -e "$dir_pca_raw/chr$chr.pca${cov_pcs}.done" ]]; then
      cmds+=("$cmd")
    fi
  done
  local configured_jobs=$jobs
  jobs=$pca_jobs
  run_parallel_cmds "DiscoDivas PCA projection (PC1-PC$cov_pcs)" "${cmds[@]}"
  jobs=$configured_jobs
  Rscript "$combine_pca_r" --dir "$dir_pca_raw" --out "$pca_file" --med "$med_file" --outdir "$dir_pca" \
    --n-pc "$cov_pcs" --distance-pcs "$distance_pcs"
  need_file "$pca_file"; need_file "$dir_pca/ukb_reference_distances.tsv.gz"; need_file "$dir_pca/Fig1.PCA_QC.png"
  step_ancestry
}

step_eval() {
  need_file "$phe_file"; need_file "$pca_file"; need_file "$med_file"
  need_file "$dir_csx/score_manifest.tsv"; need_file "$evaluate_r"
  mkdir -p "$dir_eval"
  local -a cmd=(Rscript "$evaluate_r"
    --trait "$trait" --phe "$phe_file" --pca "$pca_file" --med "$med_file"
    --csx-manifest "$dir_csx/score_manifest.tsv" --out "$dir_eval"
    --prscsx-repeats "$prscsx_repeats" --prscsx-prob-min "$prscsx_prob_min"
    --validation-frac "$validation_frac" --disco-repeats "$disco_repeats"
    --fine-tune-n "$fine_tune_n" --fine-tune-prob-min "$fine_tune_prob_min"
    --require-self-match "$require_self_match" --cov-pcs "$cov_pcs"
    --distance-pcs "$distance_pcs" --min-test-group "$min_test_group"
    --disco-a "$disco_a" --run-prscsx "$run_prscsx" --run-disco "$run_disco"
    --write-predictions "$write_predictions")
  local anc_file="$ancestry_file" anc_col="$ancestry_col" anc_prob="$ancestry_prob_col" anc_self="$self_report_col" anc_id="$ancestry_id_col"
  if [[ -z "$anc_file" && "$auto_ancestry" == TRUE ]]; then
    step_ancestry
    anc_file="$auto_ancestry_file"; anc_id=eid; anc_col=genetic_ancestry
    anc_prob=ancestry_probability; anc_self=self_report_ancestry
  fi
  [[ -n "$anc_file" ]] && cmd+=(--ancestry-file "$anc_file" --ancestry-id-col "$anc_id")
  [[ -n "$anc_col" ]] && cmd+=(--ancestry-col "$anc_col")
  [[ -n "$anc_prob" ]] && cmd+=(--ancestry-prob-col "$anc_prob")
  [[ -n "$anc_self" ]] && cmd+=(--self-report-col "$anc_self")
  [[ -n "$target_ids" ]] && cmd+=(--target-ids "$target_ids")
  [[ -n "$exclude_ids" ]] && cmd+=(--exclude-ids "$exclude_ids")
  "${cmd[@]}"

  need_file "$dir_eval/evaluation_results.xlsx"
  [[ "$run_prscsx" == TRUE ]] && need_file "$dir_eval/Fig1.PRSCSx_evaluation.png"
  [[ "$run_disco" == TRUE ]] && {
    need_file "$dir_eval/Fig2.DiscoDivas_primary.png"
    need_file "$dir_eval/Fig3.model_transfer.png"
    need_file "$dir_eval/Fig4.interpolation_diagnostics.png"
    need_file "$dir_eval/Fig5.continuum.png"
    need_file "$dir_eval/Fig6.calibration_stability.png"
  }
  need_file "$dir_eval/Fig7.ancestry_data_QC.png"
  need_file "$dir_eval/deployment_models.rds"
  need_file "$dir_eval/individual_input_template.tsv"
  log "Evaluation complete: $dir_eval"
}

case "$step" in
  csx) run_timed_step csx step_csx ;;
  pca) run_timed_step pca step_pca ;;
  ancestry) run_timed_step ancestry step_ancestry ;;
  eval) run_timed_step eval step_eval ;;
  all) run_timed_step csx step_csx; run_timed_step pca step_pca; run_timed_step eval step_eval ;;
esac

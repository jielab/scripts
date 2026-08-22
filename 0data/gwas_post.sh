#!/bin/bash
# Clean very large GWAS Catalog summary statistics.
# Version: gwas_post.sh
#
# Logic:
#   1) Format raw GWAS; optionally retain only HM3/P <= p_small/cis variants
#   2) MAGMA gene analysis on the full standardized GWAS
#   3) optional liftOver / small-GWAS and other post-processing
#
# Main outputs:
#   clean/<GWAS>/<GWAS>.gz       final standardized GWAS used downstream
#   clean/<GWAS>/<GWAS>.cis.gz   optional cis subset, when --cis-bed is provided
#   clean/<GWAS>/<GWAS>.awk.snp  distance-based lead SNPs from awk
#   clean/<GWAS>/<GWAS>.<suffix> concatenated clump/COJO outputs
#   clean/<GWAS>/clump, cojo     temporary chr-specific outputs, removed after successful lead step
#   cmd/<GWAS>.cmd               one runnable command file per GWAS

set -euo pipefail
export LC_ALL=C

SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
ORIGINAL_ARGS=("$@")

usage() {
cat <<'USAGE'
Usage:
  ./gwas_post.sh [step] [--dir-raw dir_raw] [--dir-clean dir_clean] [options]

Core options:
  step                   comma-separated format|magma|liftover|cis|lead|mh_plot, or all [all]
  --dir0 PATH            root path [/mnt/d]
  --dir-raw PATH         raw GWAS input dir; if PATH/raw exists, that subdir is used [data.BIG/gwas/[label]/raw]
  --dir-out PATH         output dir; clean/cmd/qc are created here [data.BIG/gwas/[label]]
  --dir-clean PATH       standardized input folder containing <GWAS>/<GWAS>.gz [<dir-out>/clean]
  --dir-magma PATH       MAGMA output root [sibling magma folder of --dir-clean]
  --gwas LIST            comma-separated GWAS names to process [all discovered]
  --label STR            GWAS label under data.BIG/gwas/[label]; default is basename of --dir-raw, e.g. met or ppp
  --small TRUE/FALSE     during format, retain HM3/P <= p-small/cis variants [TRUE]
  --jobs INT             local parallel jobs [4]
  --force TRUE/FALSE     recreate existing outputs [FALSE]
  --run-cmd TRUE/FALSE   execute generated cmd files after writing them [FALSE]
  --is.bsub TRUE/FALSE   submit cmd files to LSF bsub instead of local parallel [FALSE]
Step inputs and defaults:
  format:    raw *.gz from [dir-raw], --phef [/mnt/d/scripts/0f/0phe.f.sh],
             --hm3 [/mnt/d/data.BIG/refLD/ldsc/hm3/w_hm3.snp], --p-small [1e-4],
             retaining all variants in --cis-bed regions when provided,
             excluding chr*.bad.snp from --refGen_clump
  magma:     --dir-clean/<GWAS>/<GWAS>.gz -> --dir-magma/<GWAS>/*;
             defaults: magma from PATH, GRCh37, /mnt/d/files/NCBI37.3.gene.loc,
             /mnt/d/data.BIG/refLD/magma/g1000_eur and dbSNP151 synonyms
  liftover:  clean/<GWAS>/<GWAS>.small.gz, --chain [/mnt/d/files/liftOver/hg19ToHg38.over.chain.gz],
             --liftover-bin [liftOver], --liftOver [FALSE]
  cis:       clean/<GWAS>/<GWAS>.gz,
             --cis-bed [cis_bed], --cis-flank [100000]
  lead:      clean/<GWAS>/<GWAS>.gz,
             --refGen_clump [/mnt/d/data.BIG/refGen/1kg/GRCH37/pfile/] containing chr*.psam/pgen/pvar,
             --refGen_cojo [/mnt/d/data.BIG/refGen/1kg/GRCH37/pfile/EUR.] containing EUR.chr*.bed/bim/fam,
             --p-lead [5e-8], --lead-window [1000000],
             temporary chr-specific outputs in clean/<GWAS>/clump and clean/<GWAS>/cojo;
             concatenated outputs in clean/<GWAS>; clump/cojo dirs and *.cma.cojo are removed after success
  mh_plot:   compare clean/<GWAS>/<GWAS>.gz with mh_plot/<GWAS>.png and plot missing PNGs;
             uses --mh-plot-bed first for other loci, falling back to --cis-bed;
             with --small FALSE, plots only HM3 or P < 0.001 without changing clean/<GWAS>/<GWAS>.gz
  all:       for raw inputs missing small output: small + optional liftover + mh_plot + cis + lead
             Existing small outputs are skipped even when downstream outputs are missing;
             use separate cis, lead, and mh_plot steps to fill those missing outputs.
             Before dispatch, deletes clean/<GWAS> folders having a non-empty .err
             but no valid <GWAS>.gz, then processes raw inputs missing completed output.

Format options:
  --phef FILE             source file containing header_names()
  --data-f FILE           shared gwas_post functions [/mnt/d/scripts/0f/0data.f.sh]
  --hm3 FILE              HapMap3 rsID list
  --p-small FLOAT         also keep variants with P <= this threshold
  --delete-raw TRUE/FALSE delete raw file after Small GWAS is gzip-valid [FALSE]

Manhattan plot options:
  --plot-f FILE           source file containing mh_plot() [/mnt/d/scripts/0f/plot.f.R]
  --mh-plot-bed FILE      primary gene BED for significant-locus annotation
                          [/mnt/d/files/glist.b38.bed]; unmatched SNPs fall back to --cis-bed

liftOver options:
  --liftOver TRUE/FALSE   after Small GWAS, liftOver clean/<GWAS>/<GWAS>.small.gz into clean/<GWAS>/<GWAS>.gz
  --chain FILE            liftOver chain
  --liftover-bin PATH     liftOver binary

cis options:
  --cis-bed FILE          cis BED-like file: CHR START END GWAS_NAME
                          [/mnt/d/data.BIG/gwas/ppp/ppp_3k.b38.bed]
  --cis-flank INT         cis flank added to START/END [100000]

Lead SNP options:
  --p-lead FLOAT          lead SNP genome-wide threshold [5e-8]
  --lead-window INT       awk top-SNP window/bin in bp [1000000]
  --chrs STR              comma-separated chromosomes for clump/cojo: autosome, autosome,X,Y, X,Y [autosome]
  --refGen_clump PATH     plink2 pfile folder/prefix; chr* is appended when PATH ends with / or . [/mnt/d/data.BIG/refGen/1kg/GRCH37/pfile/]
  --refGen_cojo PATH      gcta bfile folder/prefix; chr* is appended when PATH ends with / or . [/mnt/d/data.BIG/refGen/1kg/GRCH37/pfile/EUR.]

MAGMA options:
  --grch 37|38            GWAS genome build [37]
  --magma-ref PREFIX      LD reference bfile prefix [/mnt/d/data.BIG/refLD/magma/g1000_eur]
  --gene-loc FILE         NCBI gene location file [NCBI37.3.gene.loc or NCBI38.gene.loc]
  --synonyms FILE         dbSNP alias mapping [/mnt/d/data.BIG/annot/dbsnp/dbsnp151.synonyms]
  --window U,D            upstream/downstream gene window in kb [0,0]
  --N INT                 fixed sample size; otherwise reads <GWAS>.magma.N, uses a
                          usable per-SNP N column, or falls back to gwas_N=100000

Examples:
  ./gwas_post.sh all --dir-raw /mnt/h/gwas/met --dir-out /mnt/d/data.BIG/gwas/met --liftOver FALSE --run-cmd TRUE --jobs 12
  ./gwas_post.sh all --dir-raw /mnt/h/gwas/ppp --dir-out /mnt/d/data.BIG/gwas/ppp --liftOver FALSE --run-cmd TRUE --jobs 12
  ./gwas_post.sh format,lead,mh_plot --dir-raw /mnt/d/Downloads --dir-out /mnt/d/data.BIG/gwas/main --small FALSE --liftOver FALSE --run-cmd TRUE --jobs 3
  ./gwas_post.sh magma --dir-clean /mnt/d/data.BIG/gwas/main/clean --gwas cvd_cad,masld
  tail -f /mnt/d/data.BIG/gwas/ppp/clean/gwas_post.background.log
  pgrep -af 'gwas_post.sh'
  pstree -ap "$(pgrep -f 'gwas_post.sh --dir-raw')"

  # Re-run only lead SNP discovery from existing clean/<GWAS>/<GWAS>.gz
  ./gwas_post.sh lead --label met --jobs 12
USAGE
}

# ----------------------------- defaults -----------------------------
dir0=/mnt/d
dir_raw_arg=""
dir_out_arg=""
dir_clean_arg=""
dir_magma_arg=""
gwas_arg=""
label=""
step=all
step_set=0
jobs=4
force=FALSE
run_cmd=FALSE
run_cmd_set=0
is_bsub=FALSE
foreground=FALSE

phef=""
data_f=""
hm3=""
p_small=1e-4
small_mode=TRUE
plot_f=""
mh_plot_bed=""
delete_raw=FALSE

liftOver=FALSE
chain=""
liftover_bin=liftOver

cis_bed=""
cis_flank=100000

p_lead=5e-8
lead_window=1000000
chrs=autosome
refGen_clump=""
refGen_cojo=""

grch=37
magma_ref=""
gene_loc=""
synonyms=""
magma_window="0,0"
magma_N=""
gwas_N=100000

need_arg_value() {
  local opt="$1" val="${2-}"
  if [[ -z "$val" || "$val" == --* ]]; then
    echo "ERROR: missing value for $opt" >&2
    usage >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    format|magma|liftover|cis|lead|mh_plot|all|*,*)
      (( step_set == 0 )) || { echo "ERROR: multiple step modules supplied: $step and $1" >&2; exit 2; }
      step="$1"; step_set=1; shift;;
    --dir0) need_arg_value "$1" "${2-}"; dir0="$2"; shift 2;;
    --dir-raw|--dir_raw) need_arg_value "$1" "${2-}"; dir_raw_arg="$2"; shift 2;;
    --dir-out|--dir_out) need_arg_value "$1" "${2-}"; dir_out_arg="$2"; shift 2;;
    --dir-clean|--dir_clean) need_arg_value "$1" "${2-}"; dir_clean_arg="$2"; shift 2;;
    --dir-magma|--dir_magma) need_arg_value "$1" "${2-}"; dir_magma_arg="$2"; shift 2;;
    --gwas) need_arg_value "$1" "${2-}"; gwas_arg="$2"; shift 2;;
    --dir-gwas|--dir_gwas) need_arg_value "$1" "${2-}"; dir_raw_arg="$2"; [[ -z "$dir_out_arg" ]] && dir_out_arg="$2"; shift 2;;
    --label) need_arg_value "$1" "${2-}"; label="$2"; shift 2;;
    --jobs) need_arg_value "$1" "${2-}"; jobs="$2"; shift 2;;
    --force) need_arg_value "$1" "${2-}"; force="$2"; shift 2;;
    --run-cmd|--run_cmd) need_arg_value "$1" "${2-}"; run_cmd="$2"; run_cmd_set=1; shift 2;;
    --is.bsub|--is_bsub) need_arg_value "$1" "${2-}"; is_bsub="$2"; shift 2;;
    --foreground) need_arg_value "$1" "${2-}"; foreground="$2"; shift 2;;

    --phef) need_arg_value "$1" "${2-}"; phef="$2"; shift 2;;
    --data-f|--data_f) need_arg_value "$1" "${2-}"; data_f="$2"; shift 2;;
    --hm3) need_arg_value "$1" "${2-}"; hm3="$2"; shift 2;;
    --p-small|--p_small) need_arg_value "$1" "${2-}"; p_small="$2"; shift 2;;
    --small) need_arg_value "$1" "${2-}"; small_mode="$2"; shift 2;;
    --plot-f|--plot_f) need_arg_value "$1" "${2-}"; plot_f="$2"; shift 2;;
    --mh-plot-bed|--mh_plot_bed) need_arg_value "$1" "${2-}"; mh_plot_bed="$2"; shift 2;;
    --delete-raw|--delete_raw) need_arg_value "$1" "${2-}"; delete_raw="$2"; shift 2;;

    --liftOver|--liftover) need_arg_value "$1" "${2-}"; liftOver="$2"; shift 2;;
    --chain) need_arg_value "$1" "${2-}"; chain="$2"; shift 2;;
    --liftover-bin|--liftOver-bin|--liftover_bin) need_arg_value "$1" "${2-}"; liftover_bin="$2"; shift 2;;

    --cis-bed|--cis_bed) need_arg_value "$1" "${2-}"; cis_bed="$2"; shift 2;;
    --cis-flank|--cis_flank) need_arg_value "$1" "${2-}"; cis_flank="$2"; shift 2;;

    --p-lead|--p_lead) need_arg_value "$1" "${2-}"; p_lead="$2"; shift 2;;
    --lead-window|--lead_window) need_arg_value "$1" "${2-}"; lead_window="$2"; shift 2;;
    --chrs) need_arg_value "$1" "${2-}"; chrs="$2"; shift 2;;
    --refGen_clump|--refGen-clump|--refgen-clump|--ref-gen-clump|--refgen_clump) need_arg_value "$1" "${2-}"; refGen_clump="$2"; shift 2;;
    --refGen_cojo|--refGen-cojo|--refgen-cojo|--ref-gen-cojo|--refgen_cojo) need_arg_value "$1" "${2-}"; refGen_cojo="$2"; shift 2;;

    --grch|--GRCH) need_arg_value "$1" "${2-}"; grch="$2"; shift 2;;
    --magma-ref|--magma_ref) need_arg_value "$1" "${2-}"; magma_ref="$2"; shift 2;;
    --gene-loc|--gene_loc) need_arg_value "$1" "${2-}"; gene_loc="$2"; shift 2;;
    --synonyms) need_arg_value "$1" "${2-}"; synonyms="$2"; shift 2;;
    --window) need_arg_value "$1" "${2-}"; magma_window="$2"; shift 2;;
    --N|--n) need_arg_value "$1" "${2-}"; magma_N="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "ERROR: unknown step/module or option: $1" >&2; usage; exit 2;;
  esac
done

upper() { echo "$1" | tr '[:lower:]' '[:upper:]'; }
force=$(upper "$force")
run_cmd=$(upper "$run_cmd")
is_bsub=$(upper "$is_bsub")
foreground=$(upper "$foreground")
liftOver=$(upper "$liftOver")
delete_raw=$(upper "$delete_raw")
small_mode=$(upper "$small_mode")
step=$(echo "$step" | tr '[:upper:]' '[:lower:]')

[[ "$small_mode" == "TRUE" || "$small_mode" == "FALSE" ]] || { echo "ERROR: --small must be TRUE or FALSE" >&2; exit 2; }
# Keep `small` as a deprecated step alias, but expose the operation as `format`.
step=${step//small/format}
IFS=',' read -r -a requested_steps <<< "$step"
for requested_step in "${requested_steps[@]}"; do
  case "$requested_step" in
    format|magma|liftover|cis|lead|mh_plot|all) ;;
    *) echo "ERROR: step/module must contain format|magma|liftover|cis|lead|mh_plot, or all" >&2; exit 2;;
  esac
done
has_step(){ [[ ",$step," == *",$1,"* || "$step" == "all" ]]; }
wants_magma(){ [[ ",$step," == *",magma,"* ]]; }

# A standalone MAGMA request is an execution command, not just a command-file generator.
if wants_magma && (( run_cmd_set == 0 )); then run_cmd=TRUE; fi

if ! [[ "$jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --jobs must be a positive integer: $jobs" >&2
  exit 2
fi

[[ -z "$phef" ]] && phef="$dir0/scripts/0f/0phe.f.sh"
[[ -z "$data_f" ]] && data_f="$dir0/scripts/0f/0data.f.sh"
[[ -z "$plot_f" ]] && plot_f="$dir0/scripts/0f/plot.f.R"
[[ -z "$mh_plot_bed" ]] && mh_plot_bed="$dir0/files/glist.b38.bed"
[[ -z "$hm3" ]] && hm3="$dir0/data.BIG/refLD/ldsc/hm3/w_hm3.snp"
[[ -z "$cis_bed" ]] && cis_bed="$dir0/data.BIG/gwas/ppp/ppp_3k.b38.bed"
[[ -z "$chain" ]] && chain="$dir0/files/liftOver/hg19ToHg38.over.chain.gz"
[[ -z "$refGen_clump" ]] && refGen_clump="$dir0/data.BIG/refGen/1kg/GRCH37/pfile/"
[[ -z "$refGen_cojo" ]] && refGen_cojo="$dir0/data.BIG/refGen/1kg/GRCH37/pfile/EUR."
[[ -z "$magma_ref" ]] && magma_ref="$dir0/data.BIG/refLD/magma/g1000_eur"
[[ -z "$synonyms" ]] && synonyms="$dir0/data.BIG/annot/dbsnp/dbsnp151.synonyms"
if [[ -z "$gene_loc" ]]; then
  [[ "$grch" == 37 ]] && gene_loc="$dir0/files/NCBI37.3.gene.loc" || gene_loc="$dir0/files/NCBI38.gene.loc"
fi
[[ "$grch" == 37 || "$grch" == 38 ]] || { echo "ERROR: --grch must be 37 or 38" >&2; exit 2; }

# Keep non-interactive/background runs consistent with the normal workstation setup.
[[ -d "$dir0/software/bin" ]] && export PATH="$dir0/software/bin:$PATH"

label_from_dir_raw() {
  local p="$1" b
  p="${p%/}"
  b=$(basename "$p")
  if [[ "$b" == "raw" ]]; then
    p=$(dirname "$p")
    b=$(basename "$p")
  fi
  echo "$b"
}

if [[ -z "$label" ]]; then
  if [[ -n "$dir_raw_arg" ]]; then
    label=$(label_from_dir_raw "$dir_raw_arg")
  elif [[ -n "$dir_clean_arg" ]]; then
    label=$(basename "$(dirname "${dir_clean_arg%/}")")
  else
    label=met
  fi
fi

if [[ -n "$dir_out_arg" ]]; then
  dir_out="$dir_out_arg"
else
  dir_out="$dir0/data.BIG/gwas/$label"
fi

if [[ -n "$dir_raw_arg" ]]; then
  if [[ -d "$dir_raw_arg/raw" ]]; then
    dir_raw="$dir_raw_arg/raw"
  else
    dir_raw="$dir_raw_arg"
  fi
else
  dir_raw="$dir0/data.BIG/gwas/$label/raw"
fi

if [[ -n "$dir_clean_arg" ]]; then dir_clean="${dir_clean_arg%/}"; else dir_clean="$dir_out/clean"; fi
if [[ -n "$dir_magma_arg" ]]; then dir_magma="${dir_magma_arg%/}"; else dir_magma="$(dirname "$dir_clean")/magma"; fi
dir_cmd="$dir_out/cmd"
dir_qc="$dir_out/qc"
dir_log="$dir_clean"
mkdir -p "$dir_clean" "$dir_cmd" "$dir_qc"

if wants_magma; then
  if [[ "$grch" == "38" ]]; then
    echo "❗ GRCh build 38" >&2
  else
    echo "GRCh build 37" >&2
  fi
  echo "   gene location : $gene_loc" >&2
  echo "   LD reference  : $magma_ref" >&2
  echo "   MAGMA output  : $dir_magma" >&2
fi

if [[ "$run_cmd" == "TRUE" && "$is_bsub" != "TRUE" && "$foreground" != "TRUE" ]]; then
  background_log="$dir_clean/gwas_post.background.log"
  [[ -s "$phef" ]] || { echo "ERROR: missing or empty file: $phef" >&2; exit 1; }
  # shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
  source "$phef"
  if declare -F phe_check_existing_tasks >/dev/null 2>&1; then
    phe_check_existing_tasks "gwas_post.sh"
  fi
  phe_run_background --title "background gwas_post" "$background_log" bash "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" --foreground TRUE
  exit 0
fi

log() { echo "[$(date '+%F %T')] $*" >&2; }
need_file() { [[ -s "$1" ]] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
need_dir() { [[ -d "$1" ]] || { echo "ERROR: missing directory: $1" >&2; exit 1; }; }
need_refgen_clump() {
  local d="$1"
  if [[ -f "${d}.pgen" && -f "${d}.pvar" && -f "${d}.psam" ]]; then
    return 0
  fi
  if [[ -d "$d" ]]; then
    compgen -G "$d/chr*.pgen" >/dev/null || { echo "ERROR: no chr*.pgen files found in refGen_clump: $d" >&2; exit 1; }
    compgen -G "$d/chr*.pvar" >/dev/null || { echo "ERROR: no chr*.pvar files found in refGen_clump: $d" >&2; exit 1; }
    compgen -G "$d/chr*.psam" >/dev/null || { echo "ERROR: no chr*.psam files found in refGen_clump: $d" >&2; exit 1; }
    return 0
  fi
  compgen -G "${d}chr*.pgen" >/dev/null || { echo "ERROR: no chr*.pgen files found in refGen_clump prefix: $d" >&2; exit 1; }
  compgen -G "${d}chr*.pvar" >/dev/null || { echo "ERROR: no chr*.pvar files found in refGen_clump prefix: $d" >&2; exit 1; }
  compgen -G "${d}chr*.psam" >/dev/null || { echo "ERROR: no chr*.psam files found in refGen_clump prefix: $d" >&2; exit 1; }
}

need_refgen_cojo() {
  local d="$1"
  if [[ -f "${d}.bed" && -f "${d}.bim" && -f "${d}.fam" ]]; then
    return 0
  fi
  if [[ -d "$d" ]]; then
    compgen -G "$d/chr*.bed" >/dev/null || { echo "ERROR: no chr*.bed files found in refGen_cojo: $d" >&2; exit 1; }
    compgen -G "$d/chr*.bim" >/dev/null || { echo "ERROR: no chr*.bim files found in refGen_cojo: $d" >&2; exit 1; }
    compgen -G "$d/chr*.fam" >/dev/null || { echo "ERROR: no chr*.fam files found in refGen_cojo: $d" >&2; exit 1; }
    return 0
  fi
  compgen -G "${d}chr*.bed" >/dev/null || { echo "ERROR: no chr*.bed files found in refGen_cojo prefix: $d" >&2; exit 1; }
  compgen -G "${d}chr*.bim" >/dev/null || { echo "ERROR: no chr*.bim files found in refGen_cojo prefix: $d" >&2; exit 1; }
  compgen -G "${d}chr*.fam" >/dev/null || { echo "ERROR: no chr*.fam files found in refGen_cojo prefix: $d" >&2; exit 1; }
}

ensure_magma_resources() {
  for ext in bed bim fam; do need_file "${magma_ref}.${ext}"; done
  need_file "$gene_loc"
  need_file "$synonyms"
}
q() { printf '%q' "$1"; }

gzip_ok() { [[ -s "$1" ]] && gzip -t "$1" >/dev/null 2>&1; }

magma_meta_grch() {
  local meta="$1"
  [[ -s "$meta" ]] || return 1
  awk -F '\t' '$1=="grch" {gsub(/\r/,"",$2); print $2; exit}' "$meta"
}

magma_output_matches_build() {
  local magma_dir="$1" magma_prefix="$2" old_grch
  [[ -s "$magma_dir/magma.done" && -s "$magma_prefix.genes.out" && -s "$magma_dir/magma.meta.tsv" ]] || return 1
  old_grch=$(magma_meta_grch "$magma_dir/magma.meta.tsv" || true)
  [[ "$old_grch" == "$grch" ]]
}

cleanup_failed_output_dirs() {
  [[ "$step" == "all" ]] || return 0
  local d g final clean_real d_real
  local -a failed_dirs=()
  clean_real=$(realpath -m -- "$dir_clean")

  while IFS= read -r -d '' d; do
    g=$(basename "$d")
    final="$d/$g.gz"
    if find "$d" -type f -name '*.err' -size +0c -print -quit 2>/dev/null | grep -q . && ! gzip_ok "$final"; then
      failed_dirs+=("$d")
    fi
  done < <(find "$dir_clean" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null)

  (( ${#failed_dirs[@]} > 0 )) || { log "No failed output folders to delete."; return 0; }
  log "Failed output folders with non-empty .err and no valid <GWAS>.gz:"
  printf '  %s\n' "${failed_dirs[@]}" >&2

  for d in "${failed_dirs[@]}"; do
    d_real=$(realpath -m -- "$d")
    if [[ "$(dirname "$d_real")" != "$clean_real" || "$d_real" == "$clean_real" ]]; then
      log "ERROR: refusing to delete output folder outside clean directory: $d_real"
      exit 1
    fi
    rm -rf -- "$d_real"
  done
  log "these folders are deleted:"
  printf '  %s\n' "${failed_dirs[@]}" >&2
}

# Strip common GWAS file extensions without destroying phenotype names containing dots.
gwas_name_from_file() {
  local b
  b=$(basename "$1")
  b=${b%.gz}
  b=${b%.bgz}
  b=${b%.tsv}
  b=${b%.txt}
  b=${b%.sumstats}
  b=${b%.assoc}
  echo "$b"
}

list_raw_files() {
  [[ -d "$dir_raw" ]] || return 0
  # Bash globbing is more reliable than find on /mnt/* (DrvFS) immediately
  # after a Windows-side download/rename becomes visible to WSL.
  (
    shopt -s nullglob
    local f
    for f in "$dir_raw"/*.gz "$dir_raw"/*.bgz "$dir_raw"/*.tsv "$dir_raw"/*.txt "$dir_raw"/*.sumstats "$dir_raw"/*.assoc; do
      [[ -f "$f" && -s "$f" && "$f" != *.aria2 ]] && printf '%s\n' "$f"
    done
  ) | sort -u -V
}

validate_cis_bed_for_raw_gz() {
  local invalid_tmp bed_names_tmp f b n_missing=0
  local -a gz_files=() missing=()

  invalid_tmp=$(mktemp)
  bed_names_tmp=$(mktemp)

  awk '
    BEGIN { FS="[[:space:]]+" }
    /^[[:space:]]*($|#)/ { next }
    {
      sub(/\r$/, "", $NF)
      if (NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2) {
        print NR ":" $0
      } else {
        print $4 > names
      }
    }
  ' names="$bed_names_tmp" "$cis_bed" > "$invalid_tmp"

  if [[ -s "$invalid_tmp" ]]; then
    log "ERROR: invalid [--cis-bed] $cis_bed; every data row must have at least 4 columns and column 3 must be greater than column 2."
    sed -n '1,20p' "$invalid_tmp" >&2
    rm -f "$invalid_tmp" "$bed_names_tmp"
    exit 1
  fi

  sort -u -o "$bed_names_tmp" "$bed_names_tmp"
  mapfile -d '' -t gz_files < <(find "$dir_raw" -maxdepth 1 -type f -name '*.gz' ! -name '*.aria2' -size +0c -print0)
  if (( ${#gz_files[@]} == 0 )); then
    log "ERROR: no .gz files found in [--dir-raw] ($dir_raw)."
    rm -f "$invalid_tmp" "$bed_names_tmp"
    exit 1
  fi

  for f in "${gz_files[@]}"; do
    b=$(basename "$f")
    b=${b%.gz}
    if ! grep -Fqx -- "$b" "$bed_names_tmp"; then
      missing+=("$b")
      ((n_missing+=1))
    fi
  done

  if (( n_missing == 0 )); then
    log "${#gz_files[@]} .gz files in [--dir-raw] ($dir_raw), all have entries in [--cis-bed] ($cis_bed)."
  else
    local missing_csv
    missing_csv=$(IFS=', '; echo "${missing[*]}")
    log "ERROR: ${#gz_files[@]} .gz files in [--dir-raw] ($dir_raw), $n_missing ($missing_csv) do not have entries in [--cis-bed] ($cis_bed). Please fix this and then re-run."
    rm -f "$invalid_tmp" "$bed_names_tmp"
    exit 1
  fi

  rm -f "$invalid_tmp" "$bed_names_tmp"
}

list_names_from_dir() {
  local d="$1"
  [[ -d "$d" ]] || return 0
  {
    find "$d" -maxdepth 1 -type f \( -name '*.gz' -o -name '*.bgz' \) \
      ! -name '*.small.gz' ! -name '*.cis.gz' ! -name '*.sig.tsv.gz' ! -name '*.lead.tsv.gz' \
      -size +0c 2>/dev/null
    find "$d" -mindepth 2 -maxdepth 2 -type f \( -name '*.gz' -o -name '*.bgz' \) \
      ! -name '*.small.gz' ! -name '*.cis.gz' ! -name '*.sig.tsv.gz' ! -name '*.lead.tsv.gz' \
      -size +0c 2>/dev/null
  } | while read -r f; do gwas_name_from_file "$f"; done | sort -u -V
}

list_liftover_names() {
  [[ -d "$dir_clean" ]] || return 0
  {
    find "$dir_clean" -maxdepth 1 -type f -name '*.small.gz' -size +0c 2>/dev/null
    find "$dir_clean" -mindepth 2 -maxdepth 2 -type f -name '*.small.gz' -size +0c 2>/dev/null
  } | while read -r f; do
        b=$(basename "$f")
        b=${b%.gz}
        b=${b%.small}
        echo "$b"
      done | sort -u -V
}

raw_file_for_name() {
  local g="$1" f
  for f in "$dir_raw/$g.gz" "$dir_raw/$g.bgz" "$dir_raw/$g.tsv.gz" "$dir_raw/$g.txt.gz" "$dir_raw/$g.sumstats.gz" "$dir_raw/$g.assoc.gz" "$dir_raw/$g.tsv" "$dir_raw/$g.txt" "$dir_raw/$g.sumstats" "$dir_raw/$g.assoc" "$dir_raw/$g"; do
    [[ -s "$f" ]] && { echo "$f"; return 0; }
  done
  list_raw_files | while read -r f; do [[ "$(gwas_name_from_file "$f")" == "$g" ]] && { echo "$f"; break; }; done
}

collect_gwas_names() {
  if [[ -n "$gwas_arg" ]]; then
    tr ',' '\n' <<< "$gwas_arg" | sed '/^[[:space:]]*$/d' | sort -u -V
  elif has_step format; then
    list_raw_files | while read -r f; do gwas_name_from_file "$f"; done | sort -u -V
  elif [[ "$step" == "liftover" ]]; then
      list_liftover_names
  else
    list_names_from_dir "$dir_clean"
  fi
}

run_cmds() {
  local list="$1" n rc=0 base
  [[ -s "$list" ]] || { log "No command files in $list"; return 0; }
  n=$(wc -l < "$list" | tr -d ' ')

  if [[ "$is_bsub" == "TRUE" ]]; then
    log "Submitting $n command files to bsub; per-GWAS logs: $dir_log/<GWAS>/<GWAS>.log"
    while read -r cmd; do
      [[ -s "$cmd" ]] || continue
      base=$(basename "$cmd" .cmd)
      bsub -J "gwas_post_$base" -oo "$dir_cmd/$base.bsub.log" -eo "$dir_cmd/$base.bsub.err" \
        "mkdir -p '$dir_log/$base'; bash '$cmd' > '$dir_log/$base/$base.log' 2>&1"
    done < "$list"
    return 0
  fi

  log "Running $n command files locally (jobs=$jobs); per-GWAS logs: $dir_log/<GWAS>/<GWAS>.log"
  run_one_cmd(){
    local cmd="$1" base log_dir log_file err_file tool_err rc
    [[ -s "$cmd" ]] || return 0
    base=$(basename "$cmd" .cmd)
    log_dir="$dir_log/$base"
    log_file="$log_dir/$base.log"
    err_file="$log_dir/$base.err"
    mkdir -p "$log_dir"
    echo "[$(date '+%F %T')] START $base" >&2
    rm -f "$err_file"
    if bash "$cmd" > "$log_file" 2>&1; then
      echo "[$(date '+%F %T')] DONE  $base" >&2
    else
      rc=$?
      {
        echo "ERROR: $base failed with exit=$rc"
        echo "ERROR: log=$log_file"
        while IFS= read -r tool_err; do
          [[ -s "$tool_err" ]] || continue
          echo "ERROR: detail=$tool_err"
          grep -Ei 'error|failed|invalid' "$tool_err" || true
        done < <(find "$log_dir" -mindepth 2 -type f -name '*.err' -size +0c 2>/dev/null | sort -V)
      } > "$err_file"
      echo "[$(date '+%F %T')] FAIL  $base exit=$rc log=$err_file" >&2
      return "$rc"
    fi
  }
  export -f run_one_cmd
  export dir_log
  if command -v parallel >/dev/null 2>&1; then
    rm -f "$list.joblog"
    parallel --line-buffer -j "$jobs" --joblog "$list.joblog" run_one_cmd {} :::: "$list" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
      log "ERROR: one or more command files failed. First failed jobs from $list.joblog:"
      awk 'NR>1 && $7 != 0 {print "  exit="$7" cmd="$9; n++; if(n>=10) exit}' "$list.joblog" >&2 || true
      return "$rc"
    fi
  else
    xargs -I{} -P "$jobs" bash -c 'run_one_cmd "$1"' _ {} < "$list" || rc=$?
    return "$rc"
  fi
}

write_gwas_cmd() {
  local gwas="$1" raw gwas_dir small final cmd awk_snp cis_out qc_prefix clump_dir cojo_dir merged_prefix lead_done mh_png magma_dir magma_prefix clump_kb
  raw=$(raw_file_for_name "$gwas" || true)
  gwas_dir="$dir_clean/$gwas"
  final="$gwas_dir/$gwas.gz"
  if [[ "$liftOver" == "TRUE" ]]; then
    small="$gwas_dir/$gwas.small.gz"
  else
    small="$final"
  fi
  awk_snp="$gwas_dir/$gwas.awk.snp"
  cis_out="$gwas_dir/$gwas.cis.gz"
  qc_prefix="$dir_qc/$gwas"
  clump_dir="$gwas_dir/clump"
  cojo_dir="$gwas_dir/cojo"
  merged_prefix="$gwas_dir/$gwas"
  lead_done="$gwas_dir/$gwas.lead.done"
  mh_png="$dir_out/mh_plot/$gwas.png"
  magma_dir="$dir_magma/$gwas"
  magma_prefix="$magma_dir/$gwas"
  clump_kb=$(( (lead_window + 999) / 1000 ))
  cmd="$dir_cmd/$gwas.cmd"

  if [[ "$force" != "TRUE" && ! -s "$lead_done" && -s "$gwas_dir/$gwas.log" ]] &&
     grep -q 'lead SNP discovery done' "$gwas_dir/$gwas.log" &&
     [[ ! -d "$clump_dir" && ! -d "$cojo_dir" ]]; then
    {
      printf 'GWAS\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n'
      printf '%s\t%s\t%s\t%s\tlegacy-log-complete\n' "$gwas" "$p_lead" "$lead_window" "$chrs"
    } > "$lead_done"
  fi

  if [[ "$force" != "TRUE" ]]; then
    case "$step" in
      small)
        if gzip_ok "$small"; then
          log "SKIP completed GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      mh_plot)
        if [[ -s "$mh_png" ]]; then
          log "SKIP completed Manhattan plot: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      magma)
        if magma_output_matches_build "$magma_dir" "$magma_prefix"; then
          log "SKIP completed MAGMA: $gwas (GRCh$grch)"
          rm -f "$cmd"
          return 0
        elif [[ -s "$magma_dir/magma.done" && -s "$magma_prefix.genes.out" ]]; then
          old_grch=$(magma_meta_grch "$magma_dir/magma.meta.tsv" 2>/dev/null || true)
          log "MAGMA build mismatch/incomplete metadata for $gwas: existing GRCh${old_grch:-unknown}, requested GRCh$grch; rerun."
        fi
        ;;
      liftover)
        if gzip_ok "$final"; then
          log "SKIP completed GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      cis)
        if gzip_ok "$cis_out"; then
          log "SKIP completed GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      lead)
        if gzip_ok "$final" && [[ -s "$lead_done" ]]; then
          log "SKIP completed GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      all)
        if gzip_ok "$small"; then
          log "SKIP GWAS with completed small output: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
    esac
  fi

  cat > "$cmd" <<CMD_TOP
#!/bin/bash
set -euo pipefail
export LC_ALL=C

GWAS=$(q "$gwas")
GWAS_DIR=$(q "$gwas_dir")
RAW=$(q "$raw")
SMALL=$(q "$small")
FINAL=$(q "$final")
AWK_SNP=$(q "$awk_snp")
CIS_OUT=$(q "$cis_out")
QC_PREFIX=$(q "$qc_prefix")
CLUMP=$(q "$clump_dir")
COJO=$(q "$cojo_dir")
MERGED=$(q "$merged_prefix")
LEAD_DONE=$(q "$lead_done")

PHEF=$(q "$phef")
DATA_F=$(q "$data_f")
PLOT_F=$(q "$plot_f")
MH_PLOT_BED=$(q "$mh_plot_bed")
HM3=$(q "$hm3")
P_SMALL=$(q "$p_small")
SMALL_MODE=$(q "$small_mode")
MH_PNG=$(q "$mh_png")
MAGMA_DIR=$(q "$magma_dir")
MAGMA_PREFIX=$(q "$magma_prefix")
DELETE_RAW=$(q "$delete_raw")

DO_STEP=$(q "$step")
FORCE=$(q "$force")
DO_LIFTOVER=$(q "$liftOver")
CHAIN=$(q "$chain")
LIFTOVER_BIN=$(q "$liftover_bin")

CIS_BED=$(q "$cis_bed")
CIS_FLANK=$(q "$cis_flank")

P_LEAD=$(q "$p_lead")
LEAD_WINDOW=$(q "$lead_window")
CLUMP_KB=$(q "$clump_kb")
CHRS=$(q "$chrs")
REFGEN_CLUMP=$(q "$refGen_clump")
REFGEN_COJO=$(q "$refGen_cojo")

GRCH=$(q "$grch")
MAGMA_REF=$(q "$magma_ref")
GENE_LOC=$(q "$gene_loc")
SYNONYMS=$(q "$synonyms")
MAGMA_WINDOW=$(q "$magma_window")
MAGMA_N=$(q "$magma_N")
GWAS_N=$(q "$gwas_N")

if [[ ",\$DO_STEP," == *,format,* || "\$DO_STEP" == "all" ]]; then
  if [[ ! -s "\$RAW" ]] || ! gzip -t "\$RAW" >/dev/null 2>&1; then
    echo "ERROR: broken raw gzip for \$GWAS: \$RAW" >&2
    echo "ERROR: removing incomplete outputs for \$GWAS before re-run." >&2
    find "\$GWAS_DIR" -mindepth 1 -depth ! -path "\$GWAS_DIR/\$GWAS.log" -delete 2>/dev/null || true
    rm -f -- "\${QC_PREFIX}"* "\$MH_PNG"
    exit 1
  fi
fi

source "\$DATA_F"

# Format every raw data row into the standard 11-column schema.  Unlike
# std_small() from DATA_F, this intentionally performs no variant filtering.
std_format(){
  local src="\$1" out="\$2" header_out="\$3" tmp input_fs
  gwas_post_need_file "\$src"
  gwas_post_log "Format full GWAS: \$src -> \$out"
  gwas_post_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_post_col "\${SNP_col:-}"); CHR_col=\$(gwas_post_col "\${CHR_col:-}")
  POS_col=\$(gwas_post_col "\${POS_col:-}"); EA_col=\$(gwas_post_col "\${EA_col:-}")
  NEA_col=\$(gwas_post_col "\${NEA_col:-}"); EAF_col=\$(gwas_post_col "\${EAF_col:-}")
  N_col=\$(gwas_post_col "\${N_col:-}"); BETA_col=\$(gwas_post_col "\${BETA_col:-}")
  SE_col=\$(gwas_post_col "\${SE_col:-}"); P_col=\$(gwas_post_col "\${P_col:-}")
  LOG10P_col=\$(gwas_post_col "\${LOG10P_col:-}")
  [[ "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 ]] || { echo "ERROR: CHR/POS not detected: \$src" >&2; exit 1; }
  [[ "\$P_col" -gt 0 || "\$LOG10P_col" -gt 0 ]] || { echo "ERROR: no P or LOG10P detected: \$src" >&2; exit 1; }
  input_fs=\$(gwas_post_detect_fs "\$src"); tmp="\${out}.tmp.\$\$"
  gwas_post_zcat "\$src" | awk -v FS="\$input_fs" -v OFS='\t' \
    -v snp_col="\$SNP_col" -v chr_col="\$CHR_col" -v pos_col="\$POS_col" \
    -v ea_col="\$EA_col" -v nea_col="\$NEA_col" -v eaf_col="\$EAF_col" -v n_col="\$N_col" \
    -v beta_col="\$BETA_col" -v se_col="\$SE_col" -v p_col="\$P_col" -v logp_col="\$LOG10P_col" '
    function get(c){return c>0 ? \$c : ""}
    function val(c, x){x=get(c); return x=="" ? "NA" : x}
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function normchr(x){gsub(/^chr/,"",x); if(x=="X")x="23"; if(x=="Y")x="24"; if(x=="MT"||x=="M")x="25"; return x}
    NR==1{print "SNP","CHR","POS","EA","NEA","EAF","N","BETA","SE","P","LOG10P"; next}
    {chr=normchr(get(chr_col)); pos=get(pos_col); snp=get(snp_col); ea=get(ea_col); nea=get(nea_col)
     if(snp==""||snp=="NA"||snp=="."){snp=chr":"pos; if(ea!="")snp=snp":"ea; if(nea!="")snp=snp":"nea}
     p=(p_col>0 ? get(p_col) : ""); lp=(logp_col>0 ? get(logp_col) : "")
     if(p=="" && isnum(lp))p=10^(-lp); if(lp=="" && isnum(p) && p>0)lp=-log(p)/log(10)
     if(snp=="")snp="NA"; if(chr=="")chr="NA"; if(pos=="")pos="NA"; if(ea=="")ea="NA"; if(nea=="")nea="NA"
     if(p=="")p="NA"; if(lp=="")lp="NA"
     print snp,chr,pos,ea,nea,val(eaf_col),val(n_col),val(beta_col),val(se_col),p,lp}' | \
    { IFS= read -r format_header; printf '%s\n' "\$format_header"; sort -t \$'\t' -k2,2n -k3,3n; } | \
    gwas_post_compress > "\$tmp"
  gzip -t "\$tmp"; mv -f "\$tmp" "\$out"
  gwas_post_zcat "\$out" | awk -v g="\$GWAS" 'NR==1{next} END{print g"\t"NR-1}' > "\${QC_PREFIX}.format.nrow.tsv"
}

if [[ "\$SMALL_MODE" == "FALSE" ]]; then
  std_small(){ std_format "\$@"; }
fi

gwas_post_magma(){
  [[ ",\$DO_STEP," == *,magma,* ]] || return 0
  [[ -s "\$FINAL" ]] || { echo "ERROR: missing or empty file: \$FINAL" >&2; exit 1; }

  if [[ "\$GRCH" == "38" ]]; then
    echo "❗ [\$GWAS] GRCh build 38" >&2
  else
    echo "[\$GWAS] GRCh build 37" >&2
  fi
  echo "   gene location : \$GENE_LOC" >&2
  echo "   LD reference  : \$MAGMA_REF" >&2

  old_grch=""
  if [[ -s "\$MAGMA_DIR/magma.meta.tsv" ]]; then
    old_grch=\$(awk -F '\t' '\$1=="grch" {gsub(/\\r/,"",\$2); print \$2; exit}' "\$MAGMA_DIR/magma.meta.tsv" || true)
  fi
  if [[ "\$FORCE" != TRUE && -s "\$MAGMA_DIR/magma.done" && -s "\$MAGMA_PREFIX.genes.out" && "\$old_grch" == "\$GRCH" ]]; then
    echo "[\$(date '+%F %T')] MAGMA exists: \$MAGMA_PREFIX.genes.out (GRCh\$GRCH)" >&2
    return 0
  fi
  if [[ "\$FORCE" != TRUE && -s "\$MAGMA_DIR/magma.done" && -s "\$MAGMA_PREFIX.genes.out" && "\$old_grch" != "\$GRCH" ]]; then
    echo "[\$(date '+%F %T')] MAGMA build mismatch: existing GRCh\${old_grch:-unknown}, requested GRCh\$GRCH; rerun." >&2
  fi
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH" >&2; exit 1; }
  for ext in bed bim fam; do [[ -s "\$MAGMA_REF.\$ext" ]] || { echo "ERROR: missing or empty file: \$MAGMA_REF.\$ext" >&2; exit 1; }; done
  [[ -s "\$GENE_LOC" ]] || { echo "ERROR: missing or empty file: \$GENE_LOC" >&2; exit 1; }
  [[ -s "\$SYNONYMS" ]] || { echo "ERROR: missing or empty file: \$SYNONYMS" >&2; exit 1; }
  mkdir -p "\$MAGMA_DIR"
  # Never leave an old completion marker/meta file behind during a rerun.
  rm -f "\$MAGMA_DIR/magma.done" "\$MAGMA_DIR/magma.meta.tsv"

  if [[ -z "\$MAGMA_N" && -s "\$GWAS_DIR/\$GWAS.magma.N" ]]; then
    MAGMA_N=\$(awk 'NF{print \$1; exit}' "\$GWAS_DIR/\$GWAS.magma.N")
  fi
  if [[ -z "\$MAGMA_N" && ! "\$GWAS_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid default gwas_N: \$GWAS_N" >&2; exit 1
  fi
  if [[ -n "\$MAGMA_N" && ! "\$MAGMA_N" =~ ^[1-9][0-9]*\$ ]]; then
    echo "ERROR: invalid MAGMA sample size for \$GWAS: \$MAGMA_N" >&2; exit 1
  fi

  header=\$(set +o pipefail; gwas_clean_zcat "\$FINAL" | head -1)
  for col in SNP CHR POS P; do
    awk -F '\t' -v c="\$col" '{for(i=1;i<=NF;i++)if(\$i==c)exit 0;exit 1}' <<< "\$header" ||
      { echo "ERROR: \$FINAL lacks required column \$col" >&2; exit 1; }
  done
  snploc="\$MAGMA_PREFIX.snp.loc"; pval="\$MAGMA_PREFIX.pval"
  gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
    NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
    {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
     if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-2])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
    sort -k1,1 -k2,2n -k3,3n -u > "\$snploc"

  has_n=FALSE
  awk -F '\t' '{for(i=1;i<=NF;i++)if(\$i=="N")exit 0;exit 1}' <<< "\$header" && has_n=TRUE || true
  if [[ -n "\$MAGMA_N" ]]; then
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  elif [[ "\$has_n" == TRUE ]]; then
    { printf 'SNP\tP\tN\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);n=\$(c["N"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1&&n+0>=50)print s,p,n}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg='ncol=N'
  else
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"
  fi
  nloc=\$(wc -l < "\$snploc" | tr -d ' '); npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  if (( npval <= 1000 )) && [[ "\$narg" == ncol=N ]]; then
    echo "WARNING: no usable per-SNP N for \$GWAS; falling back to gwas_N=\$GWAS_N" >&2
    MAGMA_N="\$GWAS_N"
    { printf 'SNP\tP\n'; gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}{s=\$(c["SNP"]);p=\$(c["P"]);if(s!=""&&s!="NA"&&s!="."&&p+0>0&&p+0<=1)print s,p}' |
      sort -k1,1 -k2,2g | awk -F '\t' '!seen[\$1]++'; } > "\$pval"
    narg="N=\$MAGMA_N"; npval=\$((\$(wc -l < "\$pval" | tr -d ' ')-1))
  fi
  (( nloc > 1000 && npval > 1000 )) || { echo "ERROR: too few MAGMA SNPs: snploc=\$nloc pval=\$npval" >&2; exit 1; }
  if [[ "\$MAGMA_WINDOW" == "0,0" || "\$MAGMA_WINDOW" == 0 ]]; then
    magma --annotate --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$MAGMA_PREFIX"
  else
    magma --annotate window="\$MAGMA_WINDOW" --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$MAGMA_PREFIX"
  fi
  magma --bfile "\$MAGMA_REF" synonyms="\$SYNONYMS" --pval "\$pval" "\$narg" --gene-annot "\$MAGMA_PREFIX.genes.annot" --out "\$MAGMA_PREFIX"
  [[ -s "\$MAGMA_PREFIX.genes.out" ]] || { echo "ERROR: MAGMA output missing: \$MAGMA_PREFIX.genes.out" >&2; exit 1; }
  { printf 'key\tvalue\n'; printf 'gwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\n' \
      "\$GWAS" "\$GRCH" "\$GENE_LOC" "\$MAGMA_REF" "\$SYNONYMS" "\$MAGMA_WINDOW" "\$nloc" "\$npval"; } > "\$MAGMA_DIR/magma.meta.tsv"
  date '+%F %T' > "\$MAGMA_DIR/magma.done"
  echo "[\$(date '+%F %T')] MAGMA done: \$MAGMA_PREFIX.genes.out" >&2
}

gwas_post_mh_plot(){
  [[ ",\$DO_STEP," == *,mh_plot,* || "\$DO_STEP" == "all" ]] || return 0
  gwas_post_need_file "\$FINAL"
  if [[ "\$FORCE" != "TRUE" && -s "\$MH_PNG" ]]; then
    gwas_post_log "Manhattan plot exists: \$MH_PNG"
    return 0
  fi
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found; required for Manhattan plotting" >&2; exit 1; }
  mkdir -p "\$(dirname "\$MH_PNG")"
  plot_input="\$FINAL"
  if [[ "\$SMALL_MODE" == "FALSE" ]]; then
    plot_input="\${QC_PREFIX}.mh_plot.small.gz"
    gwas_post_log "Subset Manhattan input to HM3 or P < 0.001: \$plot_input"
    gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v hm3_file="\$HM3" -v pthr='0.001' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{while((getline x<hm3_file)>0){split(x,a,/[ \t]+/); if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1} close(hm3_file)}
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i; print; next}
      (("SNP" in c)&&\$(c["SNP"]) in hm3) || (("P" in c)&&isnum(\$(c["P"]))&&\$(c["P"])+0<pthr){print}' | gwas_post_compress > "\$plot_input"
    gzip -t "\$plot_input"
  fi
  gwas_post_log "Manhattan plot: \$plot_input -> \$MH_PNG"
  Rscript - "\$plot_input" "\$MH_PNG" "\$GWAS" "\$PLOT_F" "\$CIS_BED" "\$MH_PLOT_BED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
output <- args[[2]]
gwas <- args[[3]]
plot_f <- args[[4]]
cis_bed <- args[[5]]
mh_plot_bed <- args[[6]]
source(plot_f)

cis_gene_arg <- NULL
if (nzchar(cis_bed) && file.exists(cis_bed)) {
  cis_rows <- tryCatch(utils::read.table(cis_bed, header=FALSE, comment.char="#", stringsAsFactors=FALSE),
                       error=function(e) NULL)
  if (!is.null(cis_rows) && ncol(cis_rows) >= 4 && any(as.character(cis_rows[[4]]) == gwas))
    cis_gene_arg <- gwas
}

png_args <- list(filename = output, width = 13.333, height = 7.5, units = "in", res = 300)
if (capabilities("cairo")) png_args[["type"]] <- "cairo"
do.call(grDevices::png, png_args)
on.exit(grDevices::dev.off(), add = TRUE)
graphics::par(mar = c(5, 5, 3, 1) + 0.1)
print(mh_plot(input, col=c("gray", "darkgray"), cis_gene=cis_gene_arg, cis_bed=cis_bed, mh_plot_bed=mh_plot_bed,
              cis_color="red", other_top_color="green",
              other_top_max_per_chr=2, locus_size=1e6, main=gwas))
RSCRIPT
  [[ "\$plot_input" == "\$FINAL" ]] || rm -f "\$plot_input"
  [[ -s "\$MH_PNG" ]] || { echo "ERROR: Manhattan plot was not created: \$MH_PNG" >&2; exit 1; }
}

if [[ "\$DO_STEP" == "all" ]]; then
  # all is selected only for GWAS without small output. Defer cis so the
  # execution order is small/liftOver -> mh_plot -> cis -> lead.
  DEFER_CIS=TRUE gwas_post_run_core
  gwas_post_mh_plot
  gwas_post_make_cis
else
  REQUESTED_STEP="\$DO_STEP"
  if [[ ",\$REQUESTED_STEP," == *,format,* ]]; then DO_STEP=small; gwas_post_run_core; fi
  DO_STEP="\$REQUESTED_STEP"; gwas_post_magma
  if [[ ",\$REQUESTED_STEP," == *,liftover,* ]]; then DO_STEP=liftover; gwas_post_run_core; fi
  DO_STEP="\$REQUESTED_STEP"; gwas_post_mh_plot
  if [[ ",\$REQUESTED_STEP," == *,cis,* ]]; then DO_STEP=cis; gwas_post_run_core; fi
  if [[ ",\$REQUESTED_STEP," == *,lead,* ]]; then DO_STEP=lead; else DO_STEP=none; fi
fi

gwas_post_lead_outputs_exist(){
  [[ -s "\${MERGED}.clumps" || -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]]
}

gwas_post_lead_awk_has_no_rows(){
  [[ -s "\$AWK_SNP" ]] || return 1
  awk 'NR>1{found=1} END{exit found ? 1 : 0}' "\$AWK_SNP"
}

gwas_post_lead_done(){
  [[ "\$FORCE" != "TRUE" ]] || return 1
  [[ -s "\$AWK_SNP" ]] || return 1
  [[ ! -d "\$CLUMP" && ! -d "\$COJO" ]] || return 1
  [[ -s "\$LEAD_DONE" ]] && return 0
  gwas_post_lead_outputs_exist && return 0
  gwas_post_lead_awk_has_no_rows && return 0
  return 1
}

gwas_post_mark_lead_done(){
  local tmp
  tmp="\${LEAD_DONE}.tmp.\$\$"
  mkdir -p "\$(dirname "\$LEAD_DONE")"
  {
    printf 'GWAS\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n'
    printf '%s\t%s\t%s\t%s\tcomplete\n' "\$GWAS" "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS"
  } > "\$tmp"
  mv -f "\$tmp" "\$LEAD_DONE"
}

# PLINK 2 requires an A1 column when an ID refers to a multiallelic reference
# variant.  Clump only needs SNP/P here, so remove those ambiguous IDs instead.
# The reference-derived list is cached under qc/ref_bad and reused by all GWAS.
gwas_post_filter_clump_multiallelic(){
  local ref="\$1" tag="\$2" assoc="\$3" pvar bad_dir bad tmp filtered n_before n_after n_removed
  pvar="\${ref}.pvar"
  gwas_post_need_file "\$pvar"
  bad_dir="\$(dirname "\$QC_PREFIX")/ref_bad"
  bad="\$bad_dir/\${tag}.multiallelic.snp"
  mkdir -p "\$bad_dir"

  if [[ ! -e "\$bad" || "\$pvar" -nt "\$bad" ]]; then
    tmp="\${bad}.tmp.\$\$"
    awk -v FS='\t' '
      NR==1 {for(i=1;i<=NF;i++){x=\$i; sub(/^#/,"",x); if(x=="ID")id=i; if(x=="ALT")alt=i} next}
      id>0 && alt>0 && index(\$alt,",")>0 && \$id!="" && \$id!="." {print \$id}
    ' "\$pvar" | sort -u > "\$tmp"
    mv -f "\$tmp" "\$bad"
  fi

  n_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  filtered="\${assoc}.filtered.\$\$"
  awk -v FS='\t' -v OFS='\t' '
    NR==FNR {bad[\$1]=1; next}
    FNR==1 {print; next}
    !(\$1 in bad) {print}
  ' "\$bad" "\$assoc" > "\$filtered"
  mv -f "\$filtered" "\$assoc"
  n_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  n_removed=\$((n_before-n_after))
  printf 'GWAS\tCHR\tN_BEFORE\tN_REMOVED_MULTIALLELIC\tN_AFTER\n%s\t%s\t%s\t%s\t%s\n' \
    "\$GWAS" "\$tag" "\$n_before" "\$n_removed" "\$n_after" > "\${QC_PREFIX}.\${tag}.clump_multiallelic.tsv"
  gwas_post_log "clump multiallelic filter \$tag: before=\$n_before removed=\$n_removed after=\$n_after"
}

if run_lead; then
  if gwas_post_lead_done; then
    gwas_post_log "lead SNP discovery exists: \$AWK_SNP"
  else
    lead_t0=\$(date +%s)
    rm -rf -- "\$CLUMP" "\$COJO"
    mkdir -p "\$CLUMP" "\$COJO"
    gwas_post_need_file "\$FINAL"
    gwas_post_log "awk distance lead SNPs: \$FINAL -> \$AWK_SNP"
    # 1) awk: no-LD top SNP per chromosome/window among P <= P_LEAD
    gwas_post_zcat "\$FINAL" | awk -v pthr="\$P_LEAD" -v win="\$LEAD_WINDOW" -f <(awk_lead) > "\${AWK_SNP}.tmp.\$\$"
    mv -f "\${AWK_SNP}.tmp.\$\$" "\$AWK_SNP"
    awk -v g="\$GWAS" 'NR==1{next} END{print g"\\t"NR-1}' "\$AWK_SNP" > "\${QC_PREFIX}.awk.nrow.tsv"

    labels="\${CLUMP}/labels.\$\$.txt"
    refs="\${CLUMP}/refs.\$\$.tsv"
    : > "\$labels"
    : > "\$refs"

  # 2) Prepare all chromosome inputs in one pass through FINAL.
  while read -r ref lab chr unused; do
    [[ -n "\$ref" ]] || continue
    want_chr "\$chr" "\$CHRS" || continue
    tag=\$(chr_label "\$chr")
    echo "\$tag" >> "\$labels"
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    assoc="\${cp}.assoc"
    ma="\${jp}.ma"
    cref=\$(cojo_bfile "\$REFGEN_COJO" "\$chr")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$ref" "\$tag" "\$chr" "\$assoc" "\$ma" "\$cref" >> "\$refs"
  done < <(ref_clump_pfiles "\$REFGEN_CLUMP")
  gwas_post_log "prepare chromosome inputs in one pass: \$FINAL"
  prep_inputs_all_chr "\$refs"

  # 3) plink2 --clump and gcta --cojo per chromosome.
  while IFS=\$'\t' read -r ref tag chr assoc ma cref; do
    [[ -n "\$ref" ]] || continue
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    gwas_post_filter_clump_multiallelic "\$ref" "\$tag" "\$assoc"
    if gwas_post_has_data_rows "\$assoc"; then
      run_tool "\$cp" plink2 --pfile "\$ref" --clump "\$assoc" --clump-snp-field SNP --clump-p-field P --clump-p1 "\$P_LEAD" --clump-kb "\$CLUMP_KB" --out "\$cp"
    else
      gwas_post_log "skip plink2 chr\$chr: no SNPs in \$assoc"
    fi

    # 3) gcta --cojo: COJO selection with the chromosome-specific bed/bim/fam reference
    if gwas_post_has_data_rows "\$ma"; then
      run_tool "\$jp" gcta --bfile "\$cref" --cojo-file "\$ma" --cojo-slct --cojo-p "\$P_LEAD" --out "\$jp"
    else
      gwas_post_log "skip gcta chr\$chr: no SNPs in \$ma"
    fi
  done < "\$refs"

    concat_chr_outputs "\$CLUMP" "\$MERGED" "\$labels" clumps
    concat_chr_outputs "\$COJO" "\$MERGED" "\$labels" jma.cojo cma.cojo ldr.cojo
    rm -f "\$labels" "\$refs"
    gwas_post_log "delete lead intermediates: \$CLUMP \$COJO \${MERGED}.cma.cojo"
    rm -rf -- "\$CLUMP" "\$COJO"
    rm -f -- "\${MERGED}.cma.cojo"
    gwas_post_mark_lead_done
    gwas_post_log "lead SNP discovery done in \$((\$(date +%s)-lead_t0)) sec"
  fi
fi

CMD_TOP

  chmod +x "$cmd"
  echo "$cmd"
}

# ----------------------------- checks and dispatch -----------------------------
need_file "$phef"
need_file "$data_f"
if wants_magma; then
  ensure_magma_resources
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH ($PATH)" >&2; exit 1; }
fi
if [[ "$step" == "all" ]]; then
  cleanup_failed_output_dirs
fi
if has_step mh_plot; then
  need_file "$plot_f"
fi
if { has_step format && [[ "$small_mode" == "TRUE" ]]; } || has_step mh_plot; then
  need_file "$hm3"
fi
if [[ "$liftOver" == "TRUE" ]] && has_step liftover; then
  need_file "$chain"
fi
if [[ "$step" == "cis" && -z "$cis_bed" ]]; then
  echo "ERROR: --cis-bed is required for the cis module" >&2
  exit 2
fi
if [[ -n "$cis_bed" ]] && { { has_step format && [[ "$small_mode" == "TRUE" ]]; } || has_step cis; }; then
  need_file "$cis_bed"
  need_dir "$dir_raw"
  validate_cis_bed_for_raw_gz
fi
if [[ -n "$cis_bed" ]] && has_step mh_plot; then
  need_file "$cis_bed"
  need_file "$mh_plot_bed"
fi
if { has_step format && [[ "$small_mode" == "TRUE" ]]; } || has_step lead; then
  need_refgen_clump "$refGen_clump"
fi
if has_step lead; then
  need_refgen_cojo "$refGen_cojo"
fi

log "label=$label step=$step small=$small_mode liftOver=$liftOver cis-bed=$cis_bed mh-plot-bed=$mh_plot_bed cis-flank=$cis_flank chrs=$chrs jobs=$jobs run-cmd=$run_cmd is.bsub=$is_bsub"
log "raw=$dir_raw out=$dir_out clean=$dir_clean cmd=$dir_cmd refGen_clump=$refGen_clump refGen_cojo=$refGen_cojo"
if wants_magma; then
  if [[ "$grch" == "38" ]]; then
    log "❗ GRCh build 38 | magma=$dir_magma ref=$magma_ref gene-loc=$gene_loc synonyms=$synonyms"
  else
    log "GRCh build 37 | magma=$dir_magma ref=$magma_ref gene-loc=$gene_loc synonyms=$synonyms"
  fi
fi

cmd_list="$dir_cmd/gwas_post.cmd.list"
: > "$cmd_list"

names_tmp="$dir_cmd/gwas_post.names.tmp"
collect_gwas_names > "$names_tmp"

if [[ ! -s "$names_tmp" ]]; then
  echo "ERROR: no GWAS files found for step=$step" >&2
  echo "  raw:   $dir_raw/*.gz" >&2
  echo "  clean: $dir_clean/*.gz" >&2
  exit 1
fi

log "Discovered $(wc -l < "$names_tmp" | tr -d ' ') GWAS: $(paste -sd, "$names_tmp")"

while read -r gwas; do
  [[ -n "$gwas" ]] || continue
  write_gwas_cmd "$gwas" >> "$cmd_list"
done < "$names_tmp"
rm -f "$names_tmp"

log "Created $(wc -l < "$cmd_list") per-GWAS command files."
if [[ "$run_cmd" == "TRUE" ]]; then
  run_cmds "$cmd_list"
else
  log "run-cmd=FALSE; generated command files only: $cmd_list"
fi
log "DONE"

#!/usr/bin/env bash
source /mnt/d/scripts/0f/console.sh
# Clean GWAS Catalog statistics and prepare MAGMA/clump/COJO/PGS outputs.
# Defaults, CLI and project workflow live here; implementations live in f/.

set -euo pipefail
export LC_ALL=C

SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
ORIGINAL_ARGS=("$@")


# 🚩 Command-line help
usage(){ cat <<'HELP'
Usage: gwas_format.sh format|thin|magma|liftover|cis|lead|mplot|h2|pgs|all [options]
Combine modules with commas, for example: format,liftover.

Project and execution:
  --dir-raw DIR --raw-file FILE --gwas NAME[,NAME] --label LABEL
  --dir-out DIR                Default: /mnt/i/gwas/<label>
  --category common           Output: <project>/<category>/<trait>/{gwas,magma,pgs,qc}
  --dir-clean DIR             Restrict to one trait's gwas folder
  --grch auto|37|38           Detect the build per GWAS by default
  --jobs 4 --replace FALSE --run-cmd FALSE --foreground TRUE --submit-bsub FALSE
  Explicit magma requests run by default unless --run-cmd is supplied.

Key analysis settings (all defaults are editable below in this script):
  --hm3 FALSE --p-hm3 1e-3 --fill-eaf FALSE --fill-n [NUMBER] --n-total NUMBER
  --thin TRUE --thin-chr-max 10000
    --hm3 FALSE: [gwas].gz keeps all standardized SNPs.
    --hm3 TRUE: [gwas].gz keeps HM3 SNPs or P < --p-hm3 (default 1e-3).
    --thin TRUE independently creates [gwas].thin.gz, even with --hm3 FALSE:
      HM3 or P < 1e-3 -> thinP0 on P > 1e-3 -> strict per-chromosome cap.
    With --hm3 TRUE --thin TRUE, --p-hm3 must be at least 1e-3.
    If P<=1e-3 alone exceeds the cap, retain the strongest signals and record counts.
    Thin outputs only [gwas].thin.gz and its .tbi index; each thin run regenerates them.
    The standalone thin module processes existing standardized GWAS, without raw files.
  --hm3-file FILE --hm3-pos FILE   HM3 rsID list and build-specific positions
  --liftover FALSE --chain FILE --liftover-bin liftOver
  --cis-bed FILE --cis-flank 100000
  --p-lead 5e-8 --lead-window 1000000 --chr all --refgen-pop EUR
  --refgen-id-dir DIR --refgen-clump PREFIX --refgen-cojo PREFIX
  --magma-ref PREFIX --gene-loc FILE --synonyms FILE --window 0,0 --sample-size NUMBER
  --add-panel none|magma --plot-width 13.333333 --plot-height auto --plot-res 180
  --write-sig FALSE --add-signal FILE --match-col Protein --match-value '[trait]'
  --pgs-pfile-dir DIR --pgs-threads 1
  --h2-sex unknown|male|female|mixed --h2-conda-env ldsc --h2-python PATH
  --h2-ref-ld-chr PREFIX --h2-w-ld-chr PREFIX --h2-merge-alleles FILE

Examples:
cd /mnt/d/scripts/gwas

# Main GWAS: raw downloads -> standardized tables -> MAGMA/lead SNPs -> plots/PGS.
./gwas_format.sh format --dir-raw /mnt/d/Downloads --dir-out /mnt/i/gwas/main \
  --grch auto --hm3 FALSE --run-cmd TRUE --foreground TRUE --jobs 4
./gwas_format.sh magma --label main --grch auto --run-cmd TRUE --foreground TRUE --jobs 4
./gwas_format.sh lead --label main --grch auto --run-cmd TRUE --foreground TRUE --jobs 4
./gwas_format.sh mplot --label main --grch auto --add-panel magma --write-sig TRUE --run-cmd TRUE --foreground TRUE --jobs 4
./gwas_format.sh pgs --label main --grch auto --run-cmd TRUE --foreground TRUE --jobs 4
bash /mnt/i/gwas/main/pgs/pgs.step2.cmd

# Optional protein cis extraction after formatting the prot project.
./gwas_format.sh cis --label prot --grch 38 --cis-bed /mnt/d/files/ppp_3k.38.bed --run-cmd TRUE --foreground TRUE --jobs 4

# 为已有 GWAS 生成 .thin.gz，并用它重画 mplot；不重新 format 原始 .gz。
# 在同一个 Bash / WSL 终端复制运行完整一段。
cd /mnt/d/scripts/gwas
for project in 4grid main met prot; do
  ./gwas_format.sh thin,mplot \
    --dir-out "/mnt/i/gwas/$project" --label "$project" --category common \
    --grch auto --hm3 FALSE --thin TRUE --thin-chr-max 10000 \
    --add-panel none --replace FALSE --run-cmd TRUE --foreground TRUE --jobs 4 || break
done
HELP
}


# 🚩 Defaults: inputs, execution, plots, reference panels and sample sizes
dir0=/mnt/d
dir_raw_arg=""
raw_file_arg=""
n_total=""
h2_sex=unknown
h2_ref_ld="/mnt/i/refLD/ldsc/1000G/1000G_Phase3_ldscores/LDscore."
h2_w_ld="/mnt/i/refLD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
h2_merge_alleles="/mnt/i/refLD/ldsc/hm3/w_hm3.snplist"
h2_python=""
h2_conda_env=ldsc
dir_out_arg=""
dir_clean_arg=""
dir_magma_arg=""
gwas_arg=""
label=""
category=common
step=all
step_set=0
jobs=4
replace=FALSE
fill_eaf=FALSE
fill_n=""
run_cmd=FALSE
run_cmd_set=0
is_bsub=FALSE
foreground=TRUE

phef=""
data_f=""
hm3_file=""
hm3_pos=""
p_hm3=1e-3
hm3_mode=FALSE
thin=TRUE
thin_chr_max=10000
thin_r="${SCRIPT_PATH%/*}/f/thin.R"
plot_f=""
mplot_r=""
mh_plot_bed=""
add_panel=none
plot_width=13.333333
plot_height=auto
plot_res=180
write_sig=FALSE
add_signal=""
signal_match_col=Protein
signal_match_value='[trait]'
signal_locus_pos=Gene_chr,Gene_Start,Gene_end
signal_display_col=Gene,beta,p-value
delete_raw=FALSE

liftOver=FALSE
chain=""
liftover_bin=liftOver

cis_bed=""
cis_flank=100000

p_lead=5e-8
lead_window=1000000
chrs=all
refGen_pop=EUR
refGen_id_dir=""
refGen_clump=""
refGen_cojo=""
pgs_pfile_dir=
pgs_threads=1

grch=auto
magma_ref=""
gene_loc=""
synonyms=""
magma_window="0,0"
magma_annot_cache_arg=""
magma_N=""
gwas_N=100000

# Helpers define functions only; keep user-facing configuration above.
# shellcheck source=f/gwas_format.f.sh
source "${SCRIPT_PATH%/*}/f/gwas_format.f.sh"
# shellcheck source=f/gwas_format_cmd.f.sh
source "${SCRIPT_PATH%/*}/f/gwas_format_cmd.f.sh"


# 🚩 Command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    format|thin|magma|liftover|cis|lead|mplot|h2|pgs|all|*,*)
      (( step_set == 0 )) || { echo "ERROR: multiple step modules supplied: $step and $1" >&2; exit 2; }
      step="$1"; step_set=1; shift;;
    --data-root) need_arg_value "$1" "${2-}"; dir0="$2"; shift 2;;
    --dir-raw) need_arg_value "$1" "${2-}"; dir_raw_arg="$2"; shift 2;;
    --raw-file) need_arg_value "$1" "${2-}"; raw_file_arg="$2"; shift 2;;
    --n-total) need_arg_value "$1" "${2-}"; n_total="$2"; shift 2;;
    --h2-sex) need_arg_value "$1" "${2-}"; h2_sex="$2"; shift 2;;
    --h2-ref-ld-chr) need_arg_value "$1" "${2-}"; h2_ref_ld="$2"; shift 2;;
    --h2-w-ld-chr) need_arg_value "$1" "${2-}"; h2_w_ld="$2"; shift 2;;
    --h2-merge-alleles) need_arg_value "$1" "${2-}"; h2_merge_alleles="$2"; shift 2;;
    --h2-python) need_arg_value "$1" "${2-}"; h2_python="$2"; shift 2;;
    --h2-conda-env) need_arg_value "$1" "${2-}"; h2_conda_env="$2"; shift 2;;
    --dir-out) need_arg_value "$1" "${2-}"; dir_out_arg="$2"; shift 2;;
    --category) need_arg_value "$1" "${2-}"; category="$2"; shift 2;;
    --dir-clean) need_arg_value "$1" "${2-}"; dir_clean_arg="$2"; shift 2;;
    --dir-magma) need_arg_value "$1" "${2-}"; dir_magma_arg="$2"; shift 2;;
    --gwas) need_arg_value "$1" "${2-}"; gwas_arg="$2"; shift 2;;
    --label) need_arg_value "$1" "${2-}"; label="$2"; shift 2;;
    --jobs) need_arg_value "$1" "${2-}"; jobs="$2"; shift 2;;
    --replace) need_arg_value "$1" "${2-}"; replace="$2"; shift 2;;
    --fill-eaf) need_arg_value "$1" "${2-}"; fill_eaf="$2"; shift 2;;
    --fill-n)
      if [[ -n "${2-}" && "${2-}" != --* ]]; then fill_n="$2"; shift 2; else fill_n=100000; shift; fi;;
    --run-cmd) need_arg_value "$1" "${2-}"; run_cmd="$2"; run_cmd_set=1; shift 2;;
    --submit-bsub) need_arg_value "$1" "${2-}"; is_bsub="$2"; shift 2;;
    --foreground) need_arg_value "$1" "${2-}"; foreground="$2"; shift 2;;

    --phef) need_arg_value "$1" "${2-}"; phef="$2"; shift 2;;
    --data-f) need_arg_value "$1" "${2-}"; data_f="$2"; shift 2;;
    --hm3-file) need_arg_value "$1" "${2-}"; hm3_file="$2"; shift 2;;
    --hm3-pos) need_arg_value "$1" "${2-}"; hm3_pos="$2"; shift 2;;
    --p-hm3) need_arg_value "$1" "${2-}"; p_hm3="$2"; shift 2;;
    --hm3) need_arg_value "$1" "${2-}"; hm3_mode="$2"; shift 2;;
    --small) echo 'ERROR: --small has been renamed to --hm3 TRUE/FALSE.' >&2; exit 2;;
    --p-small) echo 'ERROR: --p-small has been renamed to --p-hm3 (default: 1e-3).' >&2; exit 2;;
    --thin) need_arg_value "$1" "${2-}"; thin="$2"; shift 2;;
    --thin-chr-max) need_arg_value "$1" "${2-}"; thin_chr_max="$2"; shift 2;;
    --add-panel) need_arg_value "$1" "${2-}"; add_panel="$2"; shift 2;;
    --plot-width) need_arg_value "$1" "${2-}"; plot_width="$2"; shift 2;;
    --plot-height) need_arg_value "$1" "${2-}"; plot_height="$2"; shift 2;;
    --plot-res) need_arg_value "$1" "${2-}"; plot_res="$2"; shift 2;;
    --write-sig) need_arg_value "$1" "${2-}"; write_sig="$2"; shift 2;;
    --plot-f) need_arg_value "$1" "${2-}"; plot_f="$2"; shift 2;;
    --mh-plot-bed) need_arg_value "$1" "${2-}"; mh_plot_bed="$2"; shift 2;;
    --add-signal) need_arg_value "$1" "${2-}"; add_signal="$2"; shift 2;;
    --match-col) need_arg_value "$1" "${2-}"; signal_match_col="$2"; shift 2;;
    --match-value) need_arg_value "$1" "${2-}"; signal_match_value="$2"; shift 2;;
    --locus-pos) need_arg_value "$1" "${2-}"; signal_locus_pos="$2"; shift 2;;
    --display-col) need_arg_value "$1" "${2-}"; signal_display_col="$2"; shift 2;;
    --delete-raw) need_arg_value "$1" "${2-}"; delete_raw="$2"; shift 2;;

    --liftover) need_arg_value "$1" "${2-}"; liftOver="$2"; shift 2;;
    --chain) need_arg_value "$1" "${2-}"; chain="$2"; shift 2;;
    --liftover-bin) need_arg_value "$1" "${2-}"; liftover_bin="$2"; shift 2;;

    --cis-bed) need_arg_value "$1" "${2-}"; cis_bed="$2"; shift 2;;
    --cis-flank) need_arg_value "$1" "${2-}"; cis_flank="$2"; shift 2;;

    --p-lead) need_arg_value "$1" "${2-}"; p_lead="$2"; shift 2;;
    --lead-window) need_arg_value "$1" "${2-}"; lead_window="$2"; shift 2;;
    --chr) need_arg_value "$1" "${2-}"; chrs="$2"; shift 2;;
    --refgen-pop) need_arg_value "$1" "${2-}"; refGen_pop="$2"; shift 2;;
    --refgen-id-dir) need_arg_value "$1" "${2-}"; refGen_id_dir="$2"; shift 2;;
    --refgen-clump) need_arg_value "$1" "${2-}"; refGen_clump="$2"; shift 2;;
    --refgen-cojo) need_arg_value "$1" "${2-}"; refGen_cojo="$2"; shift 2;;
    --pgs-pfile-dir) need_arg_value "$1" "${2-}"; pgs_pfile_dir="$2"; shift 2;;
    --pgs-threads) need_arg_value "$1" "${2-}"; pgs_threads="$2"; shift 2;;

    --grch) need_arg_value "$1" "${2-}"; grch="$2"; shift 2;;
    --magma-ref) need_arg_value "$1" "${2-}"; magma_ref="$2"; shift 2;;
    --gene-loc) need_arg_value "$1" "${2-}"; gene_loc="$2"; shift 2;;
    --synonyms) need_arg_value "$1" "${2-}"; synonyms="$2"; shift 2;;
    --window) need_arg_value "$1" "${2-}"; magma_window="$2"; shift 2;;
    --magma-annot-cache) need_arg_value "$1" "${2-}"; magma_annot_cache_arg="$2"; shift 2;;
    --sample-size) need_arg_value "$1" "${2-}"; magma_N="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "ERROR: unknown step/module or option: $1" >&2; usage; exit 2;;
  esac
done

gwas_format_validate_options

# A standalone MAGMA request is an execution command, not just a command-file generator.
if wants_magma && (( run_cmd_set == 0 )); then run_cmd=TRUE; fi

# 🚩 Resolve shared resources and project paths
[[ -z "$phef" ]] && phef="$dir0/scripts/0f/0phe.f.sh"
phe_r="${phef%.sh}.R"
[[ -z "$data_f" ]] && data_f="$dir0/scripts/0f/0data.f.sh"
index_f="$dir0/scripts/0f/gwas_index.f.sh"
perf_f="${SCRIPT_PATH%/*}/f/gwas_post_perf.f.sh"
[[ -z "$plot_f" ]] && plot_f="$dir0/scripts/0f/mplot.f.R"
[[ -z "$mplot_r" ]] && mplot_r="${SCRIPT_PATH%/*}/f/gwas_post_mplot.R"
if [[ -z "$mh_plot_bed" && "$grch" != auto ]]; then mh_plot_bed="$dir0/files/glist.${grch}.bed"; fi
[[ -z "$hm3_file" ]] && hm3_file="/mnt/i/refGen/hm3/hapmap3_r3.snp"
[[ -z "$hm3_pos" ]] && hm3_pos="/mnt/i/refGen/hm3/hapmap3_r3_grch{grch}.snplist"
[[ -z "$chain" ]] && chain="$dir0/files/liftOver/hg19ToHg38.over.chain.gz"
if [[ "$grch" != auto ]]; then
  [[ -z "$refGen_clump" ]] && refGen_clump="/mnt/i/refGen/1kg/${grch}/pfile/"
  [[ -z "$refGen_cojo" ]] && refGen_cojo="/mnt/i/refGen/1kg/${grch}/pfile/${refGen_pop}/"
fi
[[ -z "$magma_ref" ]] && magma_ref="/mnt/i/refLD/magma/g1000_eur"
[[ -z "$synonyms" ]] && synonyms="/mnt/i/annot/dbsnp/dbsnp151.synonyms"
if [[ -z "$gene_loc" && "$grch" != auto ]]; then
  [[ "$grch" == 37 ]] && gene_loc="$dir0/files/NCBI.37.gene.loc" || gene_loc="$dir0/files/NCBI.38.gene.loc"
fi
[[ "$grch" == 37 || "$grch" == 38 || "$grch" == auto ]] || { echo "ERROR: --grch must be 37, 38, or auto" >&2; exit 2; }
# Keep non-interactive/background runs consistent with the normal workstation setup.
[[ -d "$dir0/software/bin" ]] && export PATH="$dir0/software/bin:$PATH"

if [[ -z "$label" ]]; then
  if [[ -n "$dir_raw_arg" ]]; then
    label=$(label_from_dir_raw "$dir_raw_arg")
  elif [[ -n "$dir_clean_arg" ]]; then
    label=$(label_from_dir_clean "$dir_clean_arg")
  else
    label=met
  fi
fi
[[ -z "$cis_bed" && "$label" == prot ]] && cis_bed="$dir0/files/ppp_3k.38.bed"

if [[ -n "$dir_out_arg" ]]; then
  dir_out="$dir_out_arg"
else
  dir_out="/mnt/i/gwas/$label"
fi

if [[ -n "$dir_raw_arg" ]]; then
  if [[ -d "$dir_raw_arg/raw" ]]; then
    dir_raw="$dir_raw_arg/raw"
  else
    dir_raw="$dir_raw_arg"
  fi
else
  dir_raw="$dir_out"
fi

if [[ -n "$dir_clean_arg" ]]; then
  dir_clean="${dir_clean_arg%/}"
  [[ "$(basename "$dir_clean")" == "gwas" ]] || {
    echo "ERROR: --dir-clean must point to one trait's gwas folder in the new layout: $dir_clean" >&2
    exit 2
  }
else
  dir_clean="$dir_out/$category"
fi
dir_magma="${dir_magma_arg%/}"
# Keep coordinator state isolated by requested module set.  This prevents a
# magma run from overwriting a concurrent pgs run's name list or command files.
step_key=${step//,/_}
dir_cmd="$dir_out/.project/$category/cmd/$step_key"
dir_log="$dir_out"
if [[ -n "$magma_annot_cache_arg" ]]; then
  magma_annot_cache="${magma_annot_cache_arg%/}"
else
  magma_annot_cache="$dir_out/.project/magma/annotation"
fi
mkdir -p "$dir_cmd"

if wants_magma; then
  if [[ "$grch" == auto ]]; then
    echo "GRCh build auto (per GWAS, 39 rsID sentinels)" >&2
  elif [[ "$grch" == "38" ]]; then
    echo " GRCh build 38" >&2
  else
    echo "GRCh build 37" >&2
  fi
  echo "   gene location : ${gene_loc:-<auto:NCBI37.3/NCBI38>}" >&2
  echo "   LD reference  : $magma_ref" >&2
  echo "   annot cache   : $magma_annot_cache" >&2
  echo "   MAGMA output  : $dir_magma" >&2
fi

gwas_format_maybe_background


# 🚩 Check resources, discover GWAS, generate commands and dispatch
need_file "$phef"
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$phef"
gwas_format_check_resources

log "label=$label step=$step hm3=$hm3_mode thin=$thin thin-chr-max=$thin_chr_max liftOver=$liftOver add-panel=$add_panel plot=${plot_width}x${plot_height}in@${plot_res}dpi write-sig=$write_sig add-signal=${add_signal:-none} match=$signal_match_col:$signal_match_value locus-pos=$signal_locus_pos display-col=$signal_display_col cis-bed=$cis_bed mh-plot-bed=${mh_plot_bed:-<auto:glist.37/glist.38>} cis-flank=$cis_flank chrs=$chrs jobs=$jobs run-cmd=$run_cmd is.bsub=$is_bsub"
log "raw=$dir_raw project=$dir_out category=$category layout=<project>/<category>/<trait>/{gwas,magma,pgs,qc} mplot=<project>/mplot coordinator-cmd=$dir_cmd refGen_clump=$refGen_clump refGen_cojo=$refGen_cojo pgs_pfile_dir=${pgs_pfile_dir:-<auto:/mnt/i/ukbGen/GRCh/imp>}"
if has_step pgs; then
  write_pgs_step2_cmd
fi
if wants_magma; then
  if [[ "$grch" == auto ]]; then
    log "GRCh build auto (per GWAS, 39 rsID sentinels) | magma=$dir_magma ref=$magma_ref gene-loc=<auto:NCBI37.3/NCBI38> annot-cache=$magma_annot_cache synonyms=$synonyms"
  elif [[ "$grch" == "38" ]]; then
    log " GRCh build 38 | magma=$dir_magma ref=$magma_ref gene-loc=$gene_loc annot-cache=$magma_annot_cache synonyms=$synonyms"
  else
    log "GRCh build 37 | magma=$dir_magma ref=$magma_ref gene-loc=$gene_loc annot-cache=$magma_annot_cache synonyms=$synonyms"
  fi
fi

cmd_list="$dir_cmd/gwas_post.cmd.list"
: > "$cmd_list"

names_tmp="$dir_cmd/gwas_post.names.tmp"
collect_gwas_names > "$names_tmp"

if [[ ! -s "$names_tmp" ]]; then
  echo "ERROR: no GWAS files found for step=$step" >&2
  echo "  raw:   $dir_raw/*.gz" >&2
  echo "  output: $dir_out/$category/<GWAS>/gwas/<GWAS>.gz" >&2
  echo "  pgs input: $dir_out/$category/<GWAS>/gwas/<GWAS>.jma.cojo" >&2
  exit 1
fi

n_discovered=$(wc -l < "$names_tmp" | tr -d ' ')
name_preview=$(head -n 20 "$names_tmp" | paste -sd, -)
(( n_discovered <= 20 )) || name_preview="$name_preview,..."
log "Discovered $n_discovered GWAS: $name_preview"

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
log "DONE [$step_key]"

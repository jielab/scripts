#!/bin/bash
# Clean GWAS Catalog statistics and prepare MAGMA/clump/COJO/PGS outputs.

set -euo pipefail
export LC_ALL=C

SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
ORIGINAL_ARGS=("$@")

usage() {
cat <<'USAGE'
Usage:
  ./gwas_post.sh [step] [--dir-raw dir_raw] [--dir-clean dir_clean] [options]

Core options:
  step                   comma-separated format|magma|liftover|cis|lead|mplot|pgs, or all [all]
  --data-root PATH       root path [/mnt/d]
  --dir-raw PATH         raw GWAS input dir; if PATH/raw exists, that subdir is used [data.BIG/gwas/[label]/raw]
  --dir-out PATH         project root [data.BIG/gwas/[label]]
  --category NAME        output category below the project, e.g. common or rare [common]
  --dir-clean PATH       one trait's gwas folder override (optional; project layout is auto-discovered)
  --dir-magma PATH       one trait's MAGMA folder override (single-trait use)
  --gwas LIST            comma-separated GWAS names to process [all discovered]
  --label STR            GWAS label under data.BIG/gwas/[label]; default is basename of --dir-raw, e.g. met or prot
  --small TRUE/FALSE     format output mode: TRUE keeps HM3 or P <= --p-small;
                          FALSE keeps all valid variants [FALSE]
  --jobs INT             local parallel jobs [4]
  --replace TRUE/FALSE   recreate existing outputs [FALSE]
  --fill-eaf TRUE/FALSE  fill missing/invalid EAF from build-matched 1KG pfiles [FALSE]
  --fill-n [VALUE]       fill missing/invalid N with VALUE; VALUE defaults to 100000
  --run-cmd TRUE/FALSE   execute generated cmd files after writing them [FALSE]
  --foreground TRUE/FALSE Run local commands in the current terminal; FALSE starts
                          a managed background run [FALSE]
  --submit-bsub TRUE/FALSE Submit cmd files to LSF instead of local parallel [FALSE]
Step inputs and defaults:
  format:    raw *.gz from [dir-raw], --phef [/mnt/d/scripts/0f/0phe.f.sh];
             --small FALSE standardizes all valid variants, while --small TRUE keeps
              the union of --hm3 rsIDs, --hm3-pos coordinates, and variants
              with P <= --p-small [1e-4]
  magma:     <project>/<category>/<GWAS>/gwas/<GWAS>.gz ->
             <project>/<category>/<GWAS>/magma/*;
             defaults: magma from PATH, per-GWAS auto-detected GRCh/gene locations,
             /mnt/d/data.BIG/refLD/magma/g1000_eur and dbSNP151 synonyms
  liftover:  <project>/<category>/<GWAS>/gwas/<GWAS>.small.gz,
             --chain [/mnt/d/files/liftOver/hg19ToHg38.over.chain.gz],
             --liftover-bin [liftOver], --liftover [FALSE]
  cis:       <project>/<category>/<GWAS>/gwas/<GWAS>.gz -> <GWAS>.cis.gz,
             retaining variants in --cis-bed regions with --cis-flank [100000]
  lead:      <project>/<category>/<GWAS>/gwas/<GWAS>.gz,
              required columns: SNP CHR POS EA NEA BETA SE P; EAF is optional,
              but a missing/invalid EAF skips GCTA-COJO and leaves cojo.done absent,
             --refgen-clump [/mnt/d/data.BIG/refGen/1kg/37/pfile/] containing chr*.psam/pgen/pvar[.zst],
             --refgen-cojo [/mnt/d/data.BIG/refGen/1kg/37/pfile/EUR/] containing chr*.bed/bim/fam,
             default population EUR uses <1kg>/<build>/id/EUR.id.2col for PLINK2,
             while COJO BED defaults to <1kg>/<build>/bfile/EUR/chr*;
             --p-lead [5e-8], --lead-window [1000000],
              temporary chr-specific outputs in <project>/<category>/<GWAS>/gwas/clump and cojo;
              independent <GWAS>.clump.done and <GWAS>.cojo.done markers;
              concatenated outputs in <project>/<category>/<GWAS>/gwas; intermediates are removed after success
  mplot:     plot <project>/<category>/<GWAS>/gwas/<GWAS>.gz as
              <project>/mplot/<GWAS>.png using the built-in self method;
              retains --mh-plot-bed/--cis-bed locus annotations. --add-panel magma mirrors
              the matching MAGMA gene-level Manhattan below Y=0 on the same chromosome x scale,
              with fixed 60% GWAS / 40% MAGMA height allocation, thinP1 compression when the
              GWAS maximum exceeds -log10(P)=50, unlabeled orange dashed GWAS P=5e-8 and
              MAGMA P=2.5e-6 threshold lines, cis red, trans green, and bold orange chromosome
              labels at Y=0. --add-signal adds matched blue vertical locus lines and bottom labels.
              Each run also appends new rows to <project>/mplot/0flag.tsv when a cis locus
              has no genome-wide significant SNP but another chromosome does.
              The display copy is thinned with thinP0(P=1e-3), so most HM3 background
              dots above that P value are omitted without changing the full GWAS or hits.
              --write-sig TRUE writes matching common lead loci and all significant
              MAGMA genes to <project>/mplot/<GWAS>.sig.txt.
  pgs:       score COJO-selected variants from
               <project>/<category>/<GWAS>/gwas/<GWAS>.jma.cojo against chromosome-wise
               UKB imputed pfiles. Uses refA as the effect allele and bJ as its joint-effect
               weight, then sums dosage-weighted chromosome scores into
               <project>/<category>/<GWAS>/pgs/<GWAS>.pgs.gz. After successful completion,
               per-chromosome <GWAS>.chr* working files are removed; only .pgs.gz,
               .pgs.done, and .pgs.meta.tsv outputs are retained. If no COJO variant
               matches the pfiles, a valid header-only .pgs.gz and .pgs.done are written
               and diagnostic PLINK logs are retained. With --replace FALSE, non-empty
               .pgs.gz and .pgs.done files alone mark a PGS as complete. For every label,
               <project>/pgs/pgs.step2.cmd is generated to merge completed trait PGS files.
               The merged <label>.pgs.gz contains eid plus one three-decimal
               <trait>.SCORE_SUM column per trait; constant PLINK2 allele counts are stored once
               per trait in <label>.ALLELE_CT.tsv instead of repeated for every participant.
  all:       for raw inputs missing formatted output: format according to --small,
               then optional liftover + mplot + cis + lead + pgs.
             Existing formatted outputs are skipped even when downstream outputs are missing;
             use separate cis, lead, and mplot steps to fill those missing outputs.
             Before dispatch, deletes clean/<GWAS> folders having a non-empty .err
             but no valid <GWAS>.gz, then processes raw inputs missing completed output.

Format options:
  --phef FILE             source file containing phe_header_names()
  --data-f FILE           shared gwas_post functions [/mnt/d/scripts/0f/0data.f.sh]
  --small TRUE/FALSE      TRUE writes only HM3 (rsID/coordinate) or P <= --p-small;
                          FALSE writes all [FALSE]
  --hm3 FILE              HapMap3 rsID list used by --small TRUE
                          [/mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp]
  --hm3-pos FILE          optional CHR/SNP/CM/POS coordinate fallback; {grch}
                          resolves per trait [hapmap3_r3_grch{grch}.snplist]
  --p-small FLOAT         P threshold used by --small TRUE [1e-4]
  --delete-raw TRUE/FALSE delete raw file after the format output is gzip-valid [FALSE]
  --fill-eaf TRUE/FALSE   preserve valid raw EAF; fill missing values from 1KG [FALSE]
  --fill-n [VALUE]        preserve valid raw N; fill missing values with VALUE [100000]

Manhattan plot options:
  --add-panel PANEL       add MAGMA below the SNP plot: magma or none [none]
  --plot-width FLOAT      PNG width in inches [13.333333]
  --plot-height FLOAT     PNG height in inches; omitted/auto preserves the current
                          4000:2600 MAGMA or 4000:1900 single-panel ratio [auto]
  --plot-res INT          PNG resolution in pixels per inch [180]
  --write-sig TRUE/FALSE  write <project>/mplot/<GWAS>.sig.txt with the plotted
                          significant common and MAGMA hits [FALSE]
  --plot-f FILE           source file containing reusable mh_plot() [/mnt/d/scripts/0f/mplot.f.R]
  --mh-plot-bed FILE      primary gene BED for self significant-locus annotation;
                          omitted selects glist.b37/b38.bed per detected build
  --add-signal FILE       tabular signal file for blue vertical mplot lines (optional)
  --match-col COLUMN      column in --add-signal matched to the current trait [Protein]
  --match-value VALUE     match value; [trait] substitutes the current GWAS trait [[trait]]
  --locus-pos C,S,E       signal chromosome,start,end columns [Gene_chr,Gene_Start,Gene_end]
  --display-col COLUMNS   comma-separated signal fields printed at plot bottom [Gene,beta,p-value]

liftOver options:
  --liftover TRUE/FALSE   after Small GWAS, liftOver clean/<GWAS>/<GWAS>.small.gz into clean/<GWAS>/<GWAS>.gz
  --chain FILE            liftOver chain
  --liftover-bin PATH     liftOver binary

cis options:
  --cis-bed FILE          cis BED-like file: CHR START END GWAS_NAME
                          [/mnt/d/files/ppp_3k.38.bed for --label prot]
  --cis-flank INT         cis flank added to START/END [100000]

Lead SNP options:
  --p-lead FLOAT          lead SNP genome-wide threshold [5e-8]
  --lead-window INT       awk top-SNP window/bin in bp [1000000]
  --chr STR               comma-separated chromosomes for clump/cojo: autosome, autosome,X,Y, X,Y [autosome,X]
  --refgen-clump PATH     plink2 pfile folder/prefix; chr* is appended when PATH ends with / or . [/mnt/d/data.BIG/refGen/1kg/37/pfile/]
  --refgen-cojo PATH      gcta bfile folder/prefix; chr* is appended when PATH ends with / or . [/mnt/d/data.BIG/refGen/1kg/37/pfile/EUR/]
  --refgen-pop POP        PLINK2 reference super-population [EUR]; use ALL to disable --keep
  --refgen-id-dir PATH    directory containing <POP>.id.2col [build-matched 1KG id/]
                           GWAS IDs are automatically matched to each reference's unchanged IDs by CHR/POS/alleles

PGS options:
  --pgs-pfile-dir PATH    chromosome-wise imputed pfile directory; expects chr<CHR>.pgen,
                           chr<CHR>.pvar[.zst], and chr<CHR>.psam
                           [/mnt/h/ukbGen/<GWAS_GRCh>/imp]
  --pgs-threads INT       plink2 threads per GWAS job [1]

MAGMA options:
  --grch 37|38|auto       GWAS genome build; fixed 37/38 bypasses check_GRCH;
                          omitted/auto detects each GWAS from 39 rsID sentinels [auto]
  --magma-ref PREFIX      LD reference bfile prefix [/mnt/d/data.BIG/refLD/magma/g1000_eur]
  --gene-loc FILE         NCBI gene location file [NCBI.37.gene.loc or NCBI.38.gene.loc]
  --synonyms FILE         dbSNP alias mapping [/mnt/d/data.BIG/annot/dbsnp/dbsnp151.synonyms]
  --window U,D            upstream/downstream gene window in kb [0,0]
  --magma-annot-cache PATH Shared exact-match SNP-to-gene annotation cache
                          [<project>/.project/magma/annotation]
  --sample-size INT       fixed sample size; otherwise reads <GWAS>.magma.N, uses a
                          usable per-SNP N column, or falls back to gwas_N=100000

Examples:
  ./gwas_post.sh format --dir-raw /mnt/d/Downloads --dir-out /mnt/d/data.BIG/gwas/main --grch auto --small TRUE --p-small 1e-4 --run-cmd TRUE
  ./gwas_post.sh format --dir-raw /mnt/d/Downloads --dir-out /mnt/d/data.BIG/gwas/main --grch auto --small FALSE --run-cmd TRUE
  ./gwas_post.sh all --dir-raw /mnt/h/gwas/met --dir-out /mnt/d/data.BIG/gwas/met --grch 38 --liftover FALSE --run-cmd TRUE --jobs 12
  ./gwas_post.sh all --dir-raw /mnt/h/gwas/prot --dir-out /mnt/d/data.BIG/gwas/prot --grch 38 --liftover FALSE --run-cmd TRUE --jobs 12
  ./gwas_post.sh format,lead,mplot --dir-raw /mnt/d/Downloads --dir-out /mnt/d/data.BIG/gwas/main --grch auto --small FALSE --liftover FALSE --run-cmd TRUE --jobs 3
  ./gwas_post.sh magma --label prot --category common --grch 38 --jobs 6
  ./gwas_post.sh magma --label met --category common --grch 38 --jobs 6
  ./gwas_post.sh magma --label main --category common --grch auto --jobs 6
  ./gwas_post.sh magma,mplot --label main --add-panel magma --grch auto --jobs 6
  ./gwas_post.sh mplot --label prot --add-panel magma --grch 38 --write-sig TRUE --add-signal /mnt/d/data.BIG/gwas/prot/rare/rare_collapse.tsv --match-col Protein --match-value '[trait]' --locus-pos Gene_chr,Gene_Start,Gene_end --display-col Gene,beta,p-value --replace TRUE --run-cmd TRUE --jobs 6
  ./gwas_post.sh mplot --label met --add-panel magma --grch 38 --replace TRUE --run-cmd TRUE --jobs 6
  ./gwas_post.sh mplot --label main --add-panel magma --grch auto --replace TRUE --run-cmd TRUE --jobs 6
  ./gwas_post.sh pgs --label prot --category common --grch 38 --pgs-pfile-dir /mnt/h/ukbGen/38/imp --run-cmd TRUE --jobs 6 --foreground FALSE
  ./gwas_post.sh pgs --label met --category common --grch 38 --pgs-pfile-dir /mnt/h/ukbGen/38/imp --run-cmd TRUE --jobs 6 --foreground TRUE
  tail -f /mnt/d/data.BIG/gwas/prot/.project/common/cmd/pgs/gwas_post.background.log
  bash /mnt/d/data.BIG/gwas/prot/pgs/pgs.step2.cmd
  bash /mnt/d/data.BIG/gwas/met/pgs/pgs.step2.cmd
  pgrep -af 'gwas_post.sh'
  pstree -ap "$(pgrep -f 'gwas_post.sh --dir-raw')"

  # Re-run only lead SNP discovery from existing common/<GWAS>/gwas/<GWAS>.gz
  ./gwas_post.sh lead --label met --grch 38 --jobs 12
USAGE
}


# 🚩 defaults
dir0=/mnt/d
dir_raw_arg=""
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
foreground=FALSE

phef=""
data_f=""
hm3=""
hm3_pos=""
p_small=1e-4
small_mode=FALSE
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
chrs=autosome,X
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
    format|magma|liftover|cis|lead|mplot|pgs|all|*,*)
      (( step_set == 0 )) || { echo "ERROR: multiple step modules supplied: $step and $1" >&2; exit 2; }
      step="$1"; step_set=1; shift;;
    --data-root) need_arg_value "$1" "${2-}"; dir0="$2"; shift 2;;
    --dir-raw) need_arg_value "$1" "${2-}"; dir_raw_arg="$2"; shift 2;;
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
    --hm3) need_arg_value "$1" "${2-}"; hm3="$2"; shift 2;;
    --hm3-pos) need_arg_value "$1" "${2-}"; hm3_pos="$2"; shift 2;;
    --p-small) need_arg_value "$1" "${2-}"; p_small="$2"; shift 2;;
    --small) need_arg_value "$1" "${2-}"; small_mode="$2"; shift 2;;
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

upper() { echo "$1" | tr '[:lower:]' '[:upper:]'; }
replace=$(upper "$replace")
fill_eaf=$(upper "$fill_eaf")
run_cmd=$(upper "$run_cmd")
is_bsub=$(upper "$is_bsub")
refGen_pop=$(upper "$refGen_pop")
if [[ -z "$refGen_pop" || "$refGen_pop" == *[!A-Z0-9_-]* ]]; then
  echo "ERROR: --refgen-pop must be a simple population label: $refGen_pop" >&2
  exit 2
fi
foreground=$(upper "$foreground")
liftOver=$(upper "$liftOver")
delete_raw=$(upper "$delete_raw")
small_mode=$(upper "$small_mode")
write_sig=$(upper "$write_sig")
step=$(echo "$step" | tr '[:upper:]' '[:lower:]')
category=$(echo "$category" | tr '[:upper:]' '[:lower:]')
add_panel=$(echo "$add_panel" | tr '[:upper:]' '[:lower:]')
plot_height=$(echo "$plot_height" | tr '[:upper:]' '[:lower:]')
[[ "$add_panel" == none || "$add_panel" == magma ]] || { echo "ERROR: --add-panel must be none or magma: $add_panel" >&2; exit 2; }

[[ "$small_mode" == "TRUE" || "$small_mode" == "FALSE" ]] || { echo "ERROR: --small must be TRUE or FALSE" >&2; exit 2; }
[[ "$replace" == "TRUE" || "$replace" == "FALSE" ]] || { echo "ERROR: --replace must be TRUE or FALSE" >&2; exit 2; }
[[ "$fill_eaf" == "TRUE" || "$fill_eaf" == "FALSE" ]] || { echo "ERROR: --fill-eaf must be TRUE or FALSE" >&2; exit 2; }
[[ "$run_cmd" == "TRUE" || "$run_cmd" == "FALSE" ]] || { echo "ERROR: --run-cmd must be TRUE or FALSE" >&2; exit 2; }
[[ "$is_bsub" == "TRUE" || "$is_bsub" == "FALSE" ]] || { echo "ERROR: --submit-bsub must be TRUE or FALSE" >&2; exit 2; }
[[ "$foreground" == "TRUE" || "$foreground" == "FALSE" ]] || { echo "ERROR: --foreground must be TRUE or FALSE" >&2; exit 2; }
[[ "$liftOver" == "TRUE" || "$liftOver" == "FALSE" ]] || { echo "ERROR: --liftover must be TRUE or FALSE" >&2; exit 2; }
[[ "$delete_raw" == "TRUE" || "$delete_raw" == "FALSE" ]] || { echo "ERROR: --delete-raw must be TRUE or FALSE" >&2; exit 2; }
[[ "$write_sig" == "TRUE" || "$write_sig" == "FALSE" ]] || { echo "ERROR: --write-sig must be TRUE or FALSE" >&2; exit 2; }
[[ -n "$category" && "$category" != "." && "$category" != ".." && "$category" != *[!a-z0-9._-]* ]] || {
  echo "ERROR: --category must be a simple folder name such as common or rare: $category" >&2
  exit 2
}
if [[ -n "$fill_n" ]]; then
  awk -v n="$fill_n" 'BEGIN{exit !(n ~ /^[0-9]+([.][0-9]+)?$/ && n+0>0)}' || { echo "ERROR: --fill-n must be a positive number: $fill_n" >&2; exit 2; }
fi
step=${step//small/format}
IFS=',' read -r -a requested_steps <<< "$step"
for requested_step in "${requested_steps[@]}"; do
  case "$requested_step" in
    format|magma|liftover|cis|lead|mplot|pgs|all) ;;
    *) echo "ERROR: step/module must contain format|magma|liftover|cis|lead|mplot|pgs, or all" >&2; exit 2;;
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
if ! [[ "$pgs_threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --pgs-threads must be a positive integer: $pgs_threads" >&2
  exit 2
fi
awk -v x="$plot_width" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+$/ && x+0>0)}' || {
  echo "ERROR: --plot-width must be a positive number of inches: $plot_width" >&2
  exit 2
}
if [[ "$plot_height" != auto ]]; then
  awk -v x="$plot_height" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+$/ && x+0>0)}' || {
    echo "ERROR: --plot-height must be auto or a positive number of inches: $plot_height" >&2
    exit 2
  }
fi
if ! [[ "$plot_res" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --plot-res must be a positive integer: $plot_res" >&2
  exit 2
fi
awk -v p="$p_lead" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' || {
  echo "ERROR: --p-lead must be a number in (0,1]: $p_lead" >&2
  exit 2
}
awk -v p="$p_small" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' || {
  echo "ERROR: --p-small must be a number in (0,1]: $p_small" >&2
  exit 2
}
if ! [[ "$lead_window" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --lead-window must be a positive integer: $lead_window" >&2
  exit 2
fi

[[ -z "$phef" ]] && phef="$dir0/scripts/0f/0phe.f.sh"
[[ -z "$data_f" ]] && data_f="$dir0/scripts/0f/0data.f.sh"
index_f="$dir0/scripts/0f/gwas_index.f.sh"
perf_f="$dir0/scripts/0data/f/gwas_post_perf.f.sh"
[[ -z "$plot_f" ]] && plot_f="$dir0/scripts/0f/mplot.f.R"
[[ -z "$mplot_r" ]] && mplot_r="$dir0/scripts/0data/f/gwas_post_mplot.R"
if [[ -z "$mh_plot_bed" && "$grch" != auto ]]; then mh_plot_bed="$dir0/files/glist.${grch}.bed"; fi
[[ -z "$hm3" ]] && hm3="$dir0/data.BIG/refGen/hm3/hapmap3_r3.snp"
[[ -z "$hm3_pos" ]] && hm3_pos="$dir0/data.BIG/refGen/hm3/hapmap3_r3_grch{grch}.snplist"
[[ -z "$chain" ]] && chain="$dir0/files/liftOver/hg19ToHg38.over.chain.gz"
if [[ "$grch" != auto ]]; then
  [[ -z "$refGen_clump" ]] && refGen_clump="$dir0/data.BIG/refGen/1kg/${grch}/pfile/"
  [[ -z "$refGen_cojo" ]] && refGen_cojo="$dir0/data.BIG/refGen/1kg/${grch}/pfile/${refGen_pop}/"
fi
[[ -z "$magma_ref" ]] && magma_ref="$dir0/data.BIG/refLD/magma/g1000_eur"
[[ -z "$synonyms" ]] && synonyms="$dir0/data.BIG/annot/dbsnp/dbsnp151.synonyms"
if [[ -z "$gene_loc" && "$grch" != auto ]]; then
  [[ "$grch" == 37 ]] && gene_loc="$dir0/files/NCBI.37.gene.loc" || gene_loc="$dir0/files/NCBI.38.gene.loc"
fi
[[ "$grch" == 37 || "$grch" == 38 || "$grch" == auto ]] || { echo "ERROR: --grch must be 37, 38, or auto" >&2; exit 2; }
# Keep non-interactive/background runs consistent with the normal workstation setup.
[[ -d "$dir0/software/bin" ]] && export PATH="$dir0/software/bin:$PATH"

label_from_dir_raw() {
  local p="$1" b
  p="${p%/}"
  b=$(basename "$p")
  [[ "$b" != "raw" ]] || b=$(basename "$(dirname "$p")")
  echo "$b"
}

label_from_dir_clean() {
  local p="$1"
  p="${p%/}"
  basename "$(dirname "$(dirname "$(dirname "$p")")")"
}

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
  dir_out="$dir0/data.BIG/gwas/$label"
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

if [[ "$run_cmd" == "TRUE" && "$is_bsub" != "TRUE" && "$foreground" != "TRUE" ]]; then
  background_log="$dir_cmd/gwas_post.background.log"
  [[ -s "$phef" ]] || { echo "ERROR: missing or empty file: $phef" >&2; exit 1; }
  # shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
  source "$phef"
  if declare -F phe_check_existing_tasks >/dev/null 2>&1; then
    # Only reject another run of the same module set.  Different modules use
    # separate coordinator files, generated commands, and per-module logs.
    phe_check_existing_tasks "gwas_post.sh $step"
  fi
  phe_run_background --title "background gwas_post $label/$category/$step" "$background_log" bash "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" --foreground TRUE
  exit 0
fi

log() { echo "[$(date '+%F %T')] $*" >&2; }
need_file() { [[ -s "$1" ]] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
need_dir() { [[ -d "$1" ]] || { echo "ERROR: missing directory: $1" >&2; exit 1; }; }
need_refgen_clump() {
  local d="$1"
  if [[ -f "${d}.pgen" && -f "${d}.psam" ]] &&
     [[ -f "${d}.pvar" || -f "${d}.pvar.zst" ]]; then
    return 0
  fi
  if [[ -d "$d" ]]; then
    compgen -G "$d/chr*.pgen" >/dev/null || { echo "ERROR: no chr*.pgen files found in refGen_clump: $d" >&2; exit 1; }
    { compgen -G "$d/chr*.pvar" >/dev/null || compgen -G "$d/chr*.pvar.zst" >/dev/null; } || { echo "ERROR: no chr*.pvar[.zst] files found in refGen_clump: $d" >&2; exit 1; }
    compgen -G "$d/chr*.psam" >/dev/null || { echo "ERROR: no chr*.psam files found in refGen_clump: $d" >&2; exit 1; }
    return 0
  fi
  compgen -G "${d}chr*.pgen" >/dev/null || { echo "ERROR: no chr*.pgen files found in refGen_clump prefix: $d" >&2; exit 1; }
  { compgen -G "${d}chr*.pvar" >/dev/null || compgen -G "${d}chr*.pvar.zst" >/dev/null; } || { echo "ERROR: no chr*.pvar[.zst] files found in refGen_clump prefix: $d" >&2; exit 1; }
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

need_pgs_pfiles() {
  local d="${1%/}"
  [[ -d "$d" ]] || { echo "ERROR: missing PGS pfile directory: $d" >&2; exit 1; }
  compgen -G "$d/chr*.pgen" >/dev/null || { echo "ERROR: no chr*.pgen files found in PGS pfile directory: $d" >&2; exit 1; }
  { compgen -G "$d/chr*.pvar" >/dev/null || compgen -G "$d/chr*.pvar.zst" >/dev/null; } || {
    echo "ERROR: no chr*.pvar[.zst] files found in PGS pfile directory: $d" >&2
    exit 1
  }
  compgen -G "$d/chr*.psam" >/dev/null || { echo "ERROR: no chr*.psam files found in PGS pfile directory: $d" >&2; exit 1; }
}

ensure_magma_resources() {
  for ext in bed bim fam; do need_file "${magma_ref}.${ext}"; done
  if [[ -n "$gene_loc" ]]; then
    need_file "$gene_loc"
  elif [[ "$grch" == auto ]]; then
    need_file "$dir0/files/NCBI.37.gene.loc"
    need_file "$dir0/files/NCBI.38.gene.loc"
  else
    echo "ERROR: no MAGMA gene-location file resolved for GRCh$grch" >&2
    exit 1
  fi
  need_file "$synonyms"
}
q() { printf '%q' "$1"; }

write_pgs_step2_cmd() {
  local root="$dir_out/pgs" script="$dir_out/pgs/pgs.step2.cmd"
  mkdir -p "$root"
  rm -f -- "$root/pgs.step.cmd"
  {
    cat <<'PGS_STEP_HEADER'
#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
PGS_STEP_HEADER
    printf 'DEFAULT_PROJECT=%q\nDEFAULT_CATEGORY=%q\nDEFAULT_LABEL=%q\n' "$dir_out" "$category" "$label"
    cat <<'PGS_STEP_BODY'
PROJECT=${1:-$DEFAULT_PROJECT}
CATEGORY=${2:-$DEFAULT_CATEGORY}
LABEL=${3:-$DEFAULT_LABEL}
BATCH_SIZE=${4:-32}
[[ "$BATCH_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: batch size must be a positive integer: $BATCH_SIZE" >&2; exit 2; }
ROOT="$PROJECT/pgs"
OUTPUT="$ROOT/$LABEL.pgs.gz"
ALLELE_OUTPUT="$ROOT/$LABEL.ALLELE_CT.tsv"
DONE="$ROOT/$LABEL.pgs.done"
MANIFEST="$ROOT/$LABEL.pgs.files.tsv"
TMP_ROOT="$ROOT/.tmp"

mkdir -p "$ROOT" "$TMP_ROOT"
TMP_DIR=$(mktemp -d "$TMP_ROOT/merge.XXXXXX")
case "$TMP_DIR" in "$TMP_ROOT"/merge.*) ;; *) echo "ERROR: unsafe temporary directory: $TMP_DIR" >&2; exit 1;; esac
declare -a ACTIVE_STREAM_DIRS=()
cleanup(){
  local d
  for d in "${ACTIVE_STREAM_DIRS[@]}"; do
    case "$d" in /tmp/gwas-post-pgs-streams.*)
      [[ ! -d "$d" ]] || { rm -f -- "$d"/*; rmdir -- "$d" 2>/dev/null || true; }
      ;;
    esac
  done
  rm -rf -- "$TMP_DIR"
  rmdir -- "$TMP_ROOT" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

for tool in awk cp find gzip head mkfifo mktemp mv paste rm rmdir sort; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required command not found: $tool" >&2; exit 1; }
done
compress=(gzip -c)
if command -v pigz >/dev/null 2>&1; then compress=(pigz -c -p 4); fi

manifest_tmp="$TMP_DIR/files.tsv"
allele_tmp="$TMP_DIR/$LABEL.ALLELE_CT.tsv"
printf 'GWAS\tSTATUS\tPGS\n' > "$manifest_tmp"
printf 'trait\tALLELE_CT\n' > "$allele_tmp"
declare -a inputs=() preview=()
missing=0
empty=0
expected=0

while IFS= read -r -d '' score_file; do
  trait=${score_file##*/}
  trait=${trait%.jma.cojo}
  trait_dir="$PROJECT/$CATEGORY/$trait"
  pgs_file="$trait_dir/pgs/$trait.pgs.gz"
  done_file="$trait_dir/pgs/$trait.pgs.done"
  err_file="$trait_dir/$trait.pgs.err"
  ((expected+=1))

  if [[ ! -s "$pgs_file" || ! -s "$done_file" || -s "$err_file" ]]; then
    printf '%s\tincomplete\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    echo "ERROR: incomplete PGS for $trait: output=$pgs_file done=$done_file err=$err_file" >&2
    ((missing+=1))
    continue
  fi
  gzip -t -- "$pgs_file"
  preview=()
  mapfile -t preview < <(set +o pipefail; gzip -cd -- "$pgs_file" | head -n 2)
  expected_header=$(printf '#IID\t%s.ALLELE_CT\t%s.SCORE_SUM' "$trait" "$trait")
  [[ "${preview[0]:-}" == "$expected_header" ]] || {
    echo "ERROR: unexpected PGS header for $trait: ${preview[0]:-<empty>}" >&2
    exit 1
  }
  if [[ -z "${preview[1]:-}" ]]; then
    printf '%s\tempty_no_matched_variants\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    ((empty+=1))
  else
    IFS=$'\t' read -r first_iid allele_ct first_score extra <<< "${preview[1]}"
    [[ -n "$first_iid" && "$allele_ct" =~ ^[0-9]+$ && -n "$first_score" && -z "${extra:-}" ]] || {
      echo "ERROR: invalid first PGS data row for $trait: ${preview[1]}" >&2
      exit 1
    }
    printf '%s\t%s\n' "$trait" "$allele_ct" >> "$allele_tmp"
    printf '%s\tincluded\t%s\n' "$trait" "$pgs_file" >> "$manifest_tmp"
    inputs+=("$pgs_file")
  fi
done < <(find "$PROJECT/$CATEGORY" -mindepth 3 -maxdepth 3 -type f -path '*/gwas/*.jma.cojo' -print0 | sort -z)

(( expected > 0 )) || { echo "ERROR: no .jma.cojo inputs found under $PROJECT/$CATEGORY" >&2; exit 1; }
if (( missing > 0 )); then
  echo "ERROR: $missing of $expected expected PGS files are incomplete; merge not written." >&2
  exit 1
fi
(( ${#inputs[@]} > 0 )) || { echo "ERROR: every completed PGS is empty; merge not written." >&2; exit 1; }

merge_trait_batch(){
  local out="$1"; shift
  local n=$# f fifo pid rc producer_rc=0 i=0
  local stream_dir
  local -a fifos=() producer_pids=()
  stream_dir=$(mktemp -d /tmp/gwas-post-pgs-streams.XXXXXX)
  case "$stream_dir" in /tmp/gwas-post-pgs-streams.*) ;; *) echo "ERROR: unsafe stream directory: $stream_dir" >&2; return 1;; esac
  ACTIVE_STREAM_DIRS+=("$stream_dir")
  for f in "$@"; do
    fifo="$stream_dir/in.$(printf '%04d' "$i")"
    mkfifo -- "$fifo"
    fifos+=("$fifo")
    gzip -cd -- "$f" > "$fifo" &
    producer_pids+=("$!")
    ((i+=1))
  done
  set +e
  paste "${fifos[@]}" | awk -v FS='\t' -v OFS='\t' -v n="$n" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?$/}
    {
      expected=3*n
      if(NF!=expected){print "ERROR: batch field count mismatch at row " NR ": expected " expected ", got " NF > "/dev/stderr";exit 2}
      if(NR==1){
        printf "eid"
        for(i=0;i<n;i++){
          base=3*i
          if($(base+1)!="#IID" || $(base+3)!~ /[.]SCORE_SUM$/){print "ERROR: invalid PGS header in batch input " i+1 > "/dev/stderr";exit 2}
          printf "%s%s",OFS,$(base+3)
        }
        printf "\n"
        next
      }
      eid=$1
      printf "%s",eid
      for(i=0;i<n;i++){
        base=3*i
        if($(base+1)!=eid){print "ERROR: IID mismatch in batch at row " NR ", input " i+1 > "/dev/stderr";exit 2}
        score=$(base+3)
        if(!isnum(score)){print "ERROR: non-numeric SCORE_SUM in batch at row " NR ", input " i+1 ": " score > "/dev/stderr";exit 2}
        score+=0
        if(score>-0.0005 && score<0.0005)score=0
        printf "%s%.3f",OFS,score
      }
      printf "\n"
    }
  ' | "${compress[@]}" > "$out"
  rc=$?
  for pid in "${producer_pids[@]}"; do wait "$pid" || producer_rc=1; done
  set -e
  rm -f -- "${fifos[@]}"
  rmdir -- "$stream_dir"
  (( rc == 0 && producer_rc == 0 )) || return 1
  gzip -t -- "$out"
}

declare -a chunks=() chunk_counts=() batch=()
chunk_index=0
for f in "${inputs[@]}"; do
  batch+=("$f")
  if (( ${#batch[@]} == BATCH_SIZE )); then
    chunk="$TMP_DIR/chunk.$(printf '%04d' "$chunk_index").gz"
    merge_trait_batch "$chunk" "${batch[@]}"
    chunks+=("$chunk"); chunk_counts+=("${#batch[@]}")
    batch=(); ((chunk_index+=1))
  fi
done
if (( ${#batch[@]} > 0 )); then
  chunk="$TMP_DIR/chunk.$(printf '%04d' "$chunk_index").gz"
  merge_trait_batch "$chunk" "${batch[@]}"
  chunks+=("$chunk"); chunk_counts+=("${#batch[@]}")
fi

final_tmp="$TMP_DIR/$LABEL.pgs.gz"
if (( ${#chunks[@]} == 1 )); then
  cp -- "${chunks[0]}" "$final_tmp"
else
  stream_dir=$(mktemp -d /tmp/gwas-post-pgs-streams.XXXXXX)
  case "$stream_dir" in /tmp/gwas-post-pgs-streams.*) ;; *) echo "ERROR: unsafe stream directory: $stream_dir" >&2; exit 1;; esac
  ACTIVE_STREAM_DIRS+=("$stream_dir")
  fifos=(); producer_pids=(); i=0; producer_rc=0
  for f in "${chunks[@]}"; do
    fifo="$stream_dir/in.$(printf '%04d' "$i")"
    mkfifo -- "$fifo"
    fifos+=("$fifo")
    gzip -cd -- "$f" > "$fifo" &
    producer_pids+=("$!")
    ((i+=1))
  done
  counts=$(IFS=,; echo "${chunk_counts[*]}")
  set +e
  paste "${fifos[@]}" | awk -v FS='\t' -v OFS='\t' -v counts="$counts" '
    BEGIN{nchunk=split(counts,n,","); expected=0; for(j=1;j<=nchunk;j++)expected+=1+n[j]}
    {
      if(NF!=expected){print "ERROR: final field count mismatch at row " NR ": expected " expected ", got " NF > "/dev/stderr";exit 2}
      eid=$1; pos=1
      printf "%s",eid
      for(j=1;j<=nchunk;j++){
        if($(pos)!=eid){print "ERROR: IID mismatch between chunks at row " NR ", chunk " j > "/dev/stderr";exit 2}
        for(k=1;k<=n[j];k++)printf "%s%s",OFS,$(pos+k)
        pos+=1+n[j]
      }
      printf "\n"
    }
  ' | "${compress[@]}" > "$final_tmp"
  rc=$?
  for pid in "${producer_pids[@]}"; do wait "$pid" || producer_rc=1; done
  set -e
  rm -f -- "${fifos[@]}"
  rmdir -- "$stream_dir"
  (( rc == 0 && producer_rc == 0 )) || { echo "ERROR: final PGS merge pipeline failed" >&2; exit 1; }
fi
gzip -t -- "$final_tmp"
mv -f -- "$final_tmp" "$OUTPUT"
mv -f -- "$allele_tmp" "$ALLELE_OUTPUT"
mv -f -- "$manifest_tmp" "$MANIFEST"
done_tmp="$TMP_DIR/$LABEL.pgs.done"
printf 'LABEL\tSTATUS\tEXPECTED\tINCLUDED\tEMPTY\tTIME\n%s\tcomplete\t%s\t%s\t%s\t%s\n' \
  "$LABEL" "$expected" "${#inputs[@]}" "$empty" "$(date '+%F %T')" > "$done_tmp"
mv -f -- "$done_tmp" "$DONE"
echo "PGS merge done: $OUTPUT; allele counts: $ALLELE_OUTPUT (expected=$expected included=${#inputs[@]} empty=$empty)"
PGS_STEP_BODY
  } > "$script"
  chmod +x "$script"
  log "Wrote $label PGS merge command: $script"
}

# A tabix sidecar makes completion checks O(1). Legacy gzip files retain the
# full integrity fallback until the one-time migration creates their index.
gzip_ok() {
  local file="$1"
  [[ -s "$file" ]] || return 1
  if { [[ -s "${file}.tbi" ]] || [[ -s "${file}.csi" ]]; } && command -v tabix >/dev/null 2>&1; then
    tabix -l "$file" >/dev/null 2>&1
  else
    gzip -t "$file" >/dev/null 2>&1
  fi
}

cis_output_complete() {
  local file="$1" index=""
  [[ -s "$file" ]] || return 1
  [[ -s "${file}.tbi" ]] && index="${file}.tbi"
  [[ -z "$index" && -s "${file}.csi" ]] && index="${file}.csi"
  [[ -n "$index" && ! "$file" -nt "$index" ]] || return 1
  tabix -l "$file" >/dev/null 2>&1
}


magma_output_complete() {
  local magma_dir="$1" magma_prefix="$2"
  [[ -s "$magma_dir/magma.done" && -s "$magma_prefix.genes.out" && -s "$magma_prefix.genes.raw" ]]
}

mplot_output_complete() {
  local png="$1" final="$2" genes="$3" meta="$4" expected_grch="$5" flag="$6" aggregate="$7" sig="$8" cojo="$9"
  local signal_input=none
  [[ -z "$add_signal" ]] || signal_input="$add_signal"
  [[ -s "$png" && -s "$meta" && -s "$flag" && -s "$aggregate" && "$png" -nt "$final" ]] || return 1
  awk -F '\t' -v panel="$add_panel" -v grch="$expected_grch" -v signal="$signal_input" \
    -v match_col="$signal_match_col" -v match_value="$signal_match_value" \
    -v locus_pos="$signal_locus_pos" -v display_col="$signal_display_col" -v write_sig="$write_sig" \
    -v plot_width="$plot_width" -v plot_height="$plot_height" -v plot_res="$plot_res" '
    $1=="add_panel"&&$2==panel{p=1}
    $1=="grch"&&$2==grch{g=1}
    $1=="magma_threshold"&&$2+0==2.5e-6{t=1}
    $1=="mplot_style"&&$2==8{s=1}
    $1=="plot_width"&&$2+0==plot_width+0{x=1}
    $1=="plot_height"&&$2==plot_height{y=1}
    $1=="plot_res"&&$2+0==plot_res+0{r=1}
    $1=="write_sig"&&$2==write_sig{w=1}
    $1=="signal_input"&&$2==signal{a=1}
    $1=="signal_match_col"&&$2==match_col{b=1}
    $1=="signal_match_value"&&$2==match_value{c=1}
    $1=="signal_locus_pos"&&$2==locus_pos{d=1}
    $1=="signal_display_col"&&$2==display_col{e=1}
    END{exit !(p&&g&&t&&s&&x&&y&&r&&w&&a&&b&&c&&d&&e)}' "$meta" || return 1
  if [[ "$add_panel" == magma ]]; then
    [[ -s "$genes" && "$png" -nt "$genes" ]] || return 1
  fi
  if [[ -n "$add_signal" ]]; then [[ -s "$add_signal" && "$png" -nt "$add_signal" ]] || return 1; fi
  if [[ "$write_sig" == TRUE ]]; then
    [[ -s "$sig" && "$sig" -nt "$final" ]] || return 1
    [[ ! -s "$cojo" || "$sig" -nt "$cojo" ]] || return 1
    [[ "$add_panel" != magma || "$sig" -nt "$genes" ]] || return 1
  fi
}

pgs_output_complete() {
  local output="$1" done_file="$2"
  [[ -s "$output" && -s "$done_file" ]]
}

prune_pgs_dir() {
  local d="$1" gwas="$2" f keep_logs=FALSE meta="$1/$2.pgs.meta.tsv"
  [[ -d "$d" ]] || return 0
  if [[ -s "$meta" ]] && awk -F '\t' '$1=="matched_variants"&&$2=="none"{found=1} END{exit !found}' "$meta"; then
    keep_logs=TRUE
  fi
  for f in "$d/${gwas}.chr"*; do
    [[ -f "$f" ]] || continue
    if [[ "$keep_logs" == TRUE && "$f" == *.log ]]; then continue; fi
    rm -f -- "$f"
  done
}

prune_magma_dir() {
  local d="$1"
  [[ -d "$d" ]] || return 0
  find "$d" -mindepth 1 -maxdepth 1 -type f \
    ! -name '*.genes.out' ! -name '*.genes.raw' ! -name '*.log' \
    ! -name 'magma.meta.tsv' ! -name 'magma.done' -delete
}

cleanup_failed_output_dirs() {
  [[ "$step" == "all" ]] || return 0
  log "Automatic failed-folder deletion is disabled for the category-first layout."
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
    find "$dir_raw" -mindepth 4 -maxdepth 4 -type f -path "*/$category/*/raw/*" \
      \( -name '*.gz' -o -name '*.bgz' -o -name '*.tsv' -o -name '*.txt' -o -name '*.sumstats' -o -name '*.assoc' \) \
      -size +0c 2>/dev/null
    for f in "$dir_raw"/*.gz "$dir_raw"/*.bgz "$dir_raw"/*.tsv "$dir_raw"/*.txt "$dir_raw"/*.sumstats "$dir_raw"/*.assoc; do
      [[ -f "$f" && -s "$f" && "$f" != *.aria2 ]] && printf '%s\n' "$f"
    done
  ) | sort -u -V
}

list_names_from_dir() {
  local d="$1"
  [[ -d "$d" ]] || return 0
  find "$d" -type f \( -name '*.gz' -o -name '*.bgz' \) \
    ! -name '*.small.gz' ! -name '*.cis.gz' ! -name '*.sig.tsv.gz' ! -name '*.lead.tsv.gz' \
    -size +0c 2>/dev/null |
    while read -r f; do
      [[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
      gwas_name_from_file "$f"
    done | sort -u -V
}

list_liftover_names() {
  [[ -d "$dir_clean" ]] || return 0
  find "$dir_clean" -type f -name '*.small.gz' -size +0c 2>/dev/null |
    while read -r f; do
        [[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
        b=$(basename "$f")
        b=${b%.gz}
        b=${b%.small}
        echo "$b"
      done | sort -u -V
}

list_pgs_names() {
  [[ -d "$dir_clean" ]] || return 0
  find "$dir_clean" -type f -name '*.jma.cojo' -size +0c 2>/dev/null |
    while read -r f; do
      [[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
      b=$(basename "$f")
      echo "${b%.jma.cojo}"
    done | sort -u -V
}

raw_file_for_name() {
  local g="$1" f
  for f in \
    "$dir_out/$category/$g/raw/$g.gz" "$dir_out/$category/$g/raw/$g.bgz" \
    "$dir_out/$category/$g/raw/$g.tsv.gz" "$dir_out/$category/$g/raw/$g.txt.gz" \
    "$dir_raw/$g.gz" "$dir_raw/$g.bgz" "$dir_raw/$g.tsv.gz" "$dir_raw/$g.txt.gz" \
    "$dir_raw/$g.sumstats.gz" "$dir_raw/$g.assoc.gz" "$dir_raw/$g.tsv" "$dir_raw/$g.txt" \
    "$dir_raw/$g.sumstats" "$dir_raw/$g.assoc" "$dir_raw/$g"; do
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
  elif [[ "$step" == "pgs" ]]; then
      list_pgs_names
  else
    list_names_from_dir "$dir_clean"
  fi
}

run_cmds() {
  local list="$1" n rc=0 base
  [[ -s "$list" ]] || { log "No command files in $list"; return 0; }
  n=$(wc -l < "$list" | tr -d ' ')

  if [[ "$is_bsub" == "TRUE" ]]; then
    log "Submitting $n command files to bsub; per-GWAS logs: $dir_log/$category/<GWAS>/<GWAS>.$step_key.log"
    while read -r cmd; do
      [[ -s "$cmd" ]] || continue
      base=$(basename "$cmd" .cmd)
      bsub -J "gwas_post_${step_key}_$base" -oo "$dir_log/$category/$base/$base.$step_key.bsub.log" -eo "$dir_log/$category/$base/$base.$step_key.bsub.err" \
        "mkdir -p '$dir_log/$category/$base'; bash '$cmd' > '$dir_log/$category/$base/$base.$step_key.log' 2>&1"
    done < "$list"
    return 0
  fi

  log "Running $n command files locally (jobs=$jobs); per-GWAS logs: $dir_log/$category/<GWAS>/<GWAS>.$step_key.log"
  run_one_cmd(){
    local cmd="$1" base log_dir log_file err_file tool_err rc
    [[ -s "$cmd" ]] || return 0
    base=$(basename "$cmd" .cmd)
    log_dir="$dir_log/$category/$base"
    log_file="$log_dir/$base.$step_key.log"
    err_file="$log_dir/$base.$step_key.err"
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
  export category
  export step_key
  if command -v parallel >/dev/null 2>&1; then
    rm -f "$list.joblog"
    parallel --line-buffer -j "$jobs" --joblog "$list.joblog" run_one_cmd {} :::: "$list" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
      log "ERROR: one or more command files failed. First failed jobs from $list.joblog:"
      awk -F '\t' 'NR>1 && $7 != 0 {print "  exit="$7" cmd="$9; n++; if(n>=10) exit}' "$list.joblog" >&2 || true
      return "$rc"
    fi
  else
    xargs -I{} -P "$jobs" bash -c 'run_one_cmd "$1"' _ {} < "$list" || rc=$?
    return "$rc"
  fi
}

write_gwas_cmd() {
  local gwas="$1" raw trait_dir gwas_dir small final cmd awk_snp cis_out qc_prefix clump_dir cojo_dir merged_prefix clump_done cojo_done mh_png mh_meta mh_flag mh_sig mplot_flag_file magma_dir magma_prefix clump_kb
  local pgs_dir pgs_score_file pgs_output pgs_done pgs_meta gwas_pgs_pfile_dir
  local gwas_grch gwas_grch_cache gwas_refGen_clump gwas_refGen_cojo gwas_refGen_id_dir gwas_refGen_keep gwas_gene_loc gwas_mh_plot_bed gwas_hm3_pos detection_input cached_grch rc
  raw=""
  # Only formatting consumes the raw GWAS.  Downstream-only steps use files in
  # the trait's gwas directory, so resolving raw here would rescan the entire
  # project once per trait when no matching raw file exists.
  if has_step format; then
    raw=$(raw_file_for_name "$gwas" || true)
  fi
  if [[ -n "$dir_clean_arg" ]]; then
    gwas_dir="$dir_clean"
    trait_dir=$(dirname "$gwas_dir")
  else
    trait_dir="$dir_out/$category/$gwas"
    gwas_dir="$trait_dir/gwas"
  fi
  final="$gwas_dir/$gwas.gz"
  if [[ "$liftOver" == "TRUE" ]]; then
    small="$gwas_dir/$gwas.small.gz"
  else
    small="$final"
  fi
  awk_snp="$gwas_dir/$gwas.awk.snp"
  cis_out="$gwas_dir/$gwas.cis.gz"
  qc_prefix="$trait_dir/qc/$gwas"
  gwas_grch_cache="${qc_prefix}.grch"
  clump_dir="$gwas_dir/clump"
  cojo_dir="$gwas_dir/cojo"
  merged_prefix="$gwas_dir/$gwas"
  clump_done="$gwas_dir/$gwas.clump.done"
  cojo_done="$gwas_dir/$gwas.cojo.done"
  mh_png="$dir_out/mplot/$gwas.png"
  mh_meta="$dir_out/.project/$category/mplot/$gwas.state"
  mh_flag="$dir_out/.project/$category/mplot/flag/$gwas.tsv"
  mh_sig="$dir_out/mplot/$gwas.sig.txt"
  mplot_flag_file="$dir_out/mplot/0flag.tsv"
  pgs_dir="$trait_dir/pgs"
  pgs_score_file="${merged_prefix}.jma.cojo"
  pgs_output="$pgs_dir/$gwas.pgs.gz"
  pgs_done="$pgs_dir/$gwas.pgs.done"
  pgs_meta="$pgs_dir/$gwas.pgs.meta.tsv"
  if [[ -n "$dir_magma" ]]; then magma_dir="$dir_magma"; else magma_dir="$trait_dir/magma"; fi
  magma_prefix="$magma_dir/$gwas"
  clump_kb=$(( (lead_window + 999) / 1000 ))
  cmd="$dir_cmd/$gwas.cmd"
  mkdir -p "$gwas_dir" "$trait_dir/qc" "$(dirname "$mh_meta")"

  # Keep the standalone PGS fast path intentionally minimal: the two non-empty
  # result/marker files alone define completion.  Skip build detection, pfile
  # validation, gzip tests, metadata reads, and timestamp comparisons.
  if [[ "$replace" != "TRUE" && "$step" == pgs ]] &&
     pgs_output_complete "$pgs_output" "$pgs_done"; then
    prune_pgs_dir "$pgs_dir" "$gwas"
    log "SKIP completed PGS: $gwas"
    rm -f "$cmd"
    return 0
  fi

  # A completed lead result does not need GWAS-build detection or reference
  # validation.  Require independent clump/COJO markers plus their artifacts,
  # except for the explicit no-significant-row case.
  if [[ "$replace" != "TRUE" && "$step" == lead ]] && gzip_ok "$final" &&
     [[ -s "$awk_snp" && ! -d "$clump_dir" && ! -d "$cojo_dir" ]] &&
     { { [[ -s "$clump_done" && -s "${merged_prefix}.clumps" && -s "$cojo_done" ]] &&
           [[ -s "${merged_prefix}.jma.cojo" || -s "${merged_prefix}.ldr.cojo" ]]; } ||
       { [[ -s "$clump_done" && -s "$cojo_done" ]] && awk 'NR>1{exit 1}' "$awk_snp"; }; }; then
    log "SKIP completed lead GWAS: $gwas"
    rm -f "$cmd"
    return 0
  fi

  gwas_grch="$grch"
  # With no fixed build (omitted or --grch auto), resolve every GWAS independently
  # from the 39 sentinel rsIDs.  Use RAW before format, SMALL before a standalone
  # liftOver, and the standardized FINAL for downstream-only requests.
  if [[ "$gwas_grch" == auto ]]; then
    if has_step format; then
      detection_input="$raw"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing raw GWAS for GRCh detection: $gwas" >&2; exit 1; }
    elif [[ "$step" == liftover ]]; then
      detection_input="$small"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing small GWAS for GRCh detection: $gwas" >&2; exit 1; }
    else
      detection_input="$final"
      [[ -s "$detection_input" ]] || { echo "ERROR: missing standardized GWAS for GRCh detection: $gwas" >&2; exit 1; }
    fi
    cached_grch=""
    if [[ -s "$gwas_grch_cache" && "$gwas_grch_cache" -nt "$detection_input" ]]; then
      cached_grch=$(awk 'NR==1 && ($1==37 || $1==38){print $1}' "$gwas_grch_cache")
    fi
    if [[ -n "$cached_grch" ]]; then
      gwas_grch="$cached_grch"
      log "Reuse cached GRCh$gwas_grch: $gwas"
    elif check_GRCH "$detection_input" >&2; then
      gwas_grch="$CHECK_GRCH_RESULT"
      printf '%s\n' "$gwas_grch" > "${gwas_grch_cache}.tmp.$$"
      mv -f "${gwas_grch_cache}.tmp.$$" "$gwas_grch_cache"
    else
      rc=$?
      echo "ERROR: automatic GRCh detection failed for $gwas; the input must contain usable rsIDs from ${CHECK_GRCH_SNP_LIST:-$dir0/data/ukb/phe/common/snp.lst}, or specify --grch 37/38." >&2
      return "$rc"
    fi
  fi
  gwas_hm3_pos=${hm3_pos//\{grch\}/$gwas_grch}
  gwas_pgs_pfile_dir="${pgs_pfile_dir:-/mnt/h/ukbGen/${gwas_grch}/imp}"
  if has_step pgs; then
    need_pgs_pfiles "$gwas_pgs_pfile_dir"
  fi
  gwas_gene_loc="$gene_loc"
  if [[ -z "$gwas_gene_loc" ]]; then
    [[ "$gwas_grch" == 37 ]] && gwas_gene_loc="$dir0/files/NCBI.37.gene.loc" || gwas_gene_loc="$dir0/files/NCBI.38.gene.loc"
  fi
  gwas_mh_plot_bed=""
  gwas_mh_plot_bed="$mh_plot_bed"
  [[ -n "$gwas_mh_plot_bed" ]] || gwas_mh_plot_bed="$dir0/files/glist.${gwas_grch}.bed"
  if has_step format && [[ "$liftOver" == TRUE && "$gwas_grch" != 37 ]]; then
    echo "ERROR: --liftover TRUE uses the GRCh37-to-GRCh38 chain, but $gwas was detected/configured as GRCh$gwas_grch" >&2
    exit 1
  fi
  gwas_refGen_clump="${refGen_clump:-$dir0/data.BIG/refGen/1kg/${gwas_grch}/pfile/}"
  gwas_refGen_cojo="${refGen_cojo:-$dir0/data.BIG/refGen/1kg/${gwas_grch}/pfile/${refGen_pop}/}"
  gwas_refGen_id_dir="${refGen_id_dir:-$dir0/data.BIG/refGen/1kg/${gwas_grch}/id}"
  gwas_refGen_keep=""
  if has_step lead || wants_magma || { has_step format && [[ "$fill_eaf" == TRUE ]]; }; then
    if [[ "$refGen_pop" != ALL ]]; then
      gwas_refGen_keep="${gwas_refGen_id_dir%/}/${refGen_pop}.id.2col"
      if ! awk 'BEGIN{FS="[ \t]+"}{a=$1;b=$2;gsub(/\r/,"",a);gsub(/\r/,"",b);if(NF!=2||a!=b)bad=1;n++}END{exit bad||n==0?2:0}' "$gwas_refGen_keep"; then
        echo "ERROR: invalid PLINK two-column keep file: $gwas_refGen_keep" >&2
        exit 1
      fi
    fi
    need_refgen_clump "$gwas_refGen_clump"
  fi
  if has_step lead; then
    need_refgen_cojo "$gwas_refGen_cojo"
  fi

  if [[ "$replace" != "TRUE" ]]; then
    case "$step" in
      small)
        if gzip_ok "$small"; then
          log "SKIP completed GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      mplot)
        if mplot_output_complete "$mh_png" "$final" "$magma_prefix.genes.out" "$mh_meta" "$gwas_grch" "$mh_flag" "$mplot_flag_file" "$mh_sig" "${merged_prefix}.jma.cojo"; then
          log "SKIP completed Manhattan plot: $gwas panel=$add_panel"
          rm -f "$cmd"
          return 0
        fi
        ;;
      magma)
        if magma_output_complete "$magma_dir" "$magma_prefix"; then
          prune_magma_dir "$magma_dir"
          log "SKIP completed MAGMA: $gwas"
          rm -f "$cmd"
          return 0
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
        if cis_output_complete "$cis_out"; then
          log "SKIP completed indexed cis GWAS: $gwas"
          rm -f "$cmd"
          return 0
        fi
        ;;
      lead)
        if gzip_ok "$final" && [[ -s "$awk_snp" && ! -d "$clump_dir" && ! -d "$cojo_dir" ]] &&
           { { [[ -s "$clump_done" && -s "${merged_prefix}.clumps" && -s "$cojo_done" ]] &&
                 [[ -s "${merged_prefix}.jma.cojo" || -s "${merged_prefix}.ldr.cojo" ]]; } ||
             { [[ -s "$clump_done" && -s "$cojo_done" ]] && awk 'NR>1{exit 1}' "$awk_snp"; }; }; then
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
export PATH=$(q "$dir0/software/bin"):\$PATH

gwas_post_need_file(){
  [[ -s "\$1" ]] || { echo "ERROR: missing or empty file: \$1" >&2; exit 1; }
}
gwas_post_log(){
  printf '[%s] %s\n' "\$(date '+%F %T')" "\$*"
}
gwas_post_zcat(){
  case "\$1" in *.gz|*.bgz) gzip -cd -- "\$1";; *) cat -- "\$1";; esac
}
gwas_post_has_data_rows(){
  [[ -s "\$1" ]] && awk 'NR>1{found=1; exit} END{exit found ? 0 : 1}' "\$1"
}

GWAS=$(q "$gwas")
GWAS_DIR=$(q "$gwas_dir")
GWAS_POST_TMP_ROOT="\$GWAS_DIR/.tmp"
mkdir -p "\$GWAS_POST_TMP_ROOT"
# Isolate temporary files per process and remove them on exit.
for stale_tmp in "\$GWAS_POST_TMP_ROOT"/run.*; do
  [[ -d "\$stale_tmp" ]] || continue
  stale_pid=\${stale_tmp##*.}
  if [[ ! "\$stale_pid" =~ ^[1-9][0-9]*\$ ]] || ! kill -0 "\$stale_pid" 2>/dev/null; then
    rm -rf -- "\$stale_tmp"
  fi
done
GWAS_POST_TMP="\$GWAS_POST_TMP_ROOT/run.\$\$"
mkdir -p "\$GWAS_POST_TMP"
export TMPDIR="\$GWAS_POST_TMP"
gwas_post_cleanup_tmp(){
  rm -rf -- "\$GWAS_POST_TMP"
  rmdir -- "\$GWAS_POST_TMP_ROOT" 2>/dev/null || true
}
trap gwas_post_cleanup_tmp EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
RAW=$(q "$raw")
SMALL=$(q "$small")
FINAL=$(q "$final")
AWK_SNP=$(q "$awk_snp")
CIS_OUT=$(q "$cis_out")
QC_PREFIX=$(q "$qc_prefix")
CLUMP=$(q "$clump_dir")
COJO=$(q "$cojo_dir")
MERGED=$(q "$merged_prefix")
CLUMP_DONE=$(q "$clump_done")
COJO_DONE=$(q "$cojo_done")
REF_BAD_DIR=$(q "$dir_out/.project/$category/qc/ref_bad")

PHEF=$(q "$phef")
DATA_F=$(q "$data_f")
INDEX_F=$(q "$index_f")
PERF_F=$(q "$perf_f")
PLOT_F=$(q "$plot_f")
MPLOT_R=$(q "$mplot_r")
PLOT_METHOD=self
ADD_PANEL=$(q "$add_panel")
PLOT_WIDTH=$(q "$plot_width")
PLOT_HEIGHT=$(q "$plot_height")
PLOT_RES=$(q "$plot_res")
WRITE_SIG=$(q "$write_sig")
ADD_SIGNAL=$(q "$add_signal")
SIGNAL_MATCH_COL=$(q "$signal_match_col")
SIGNAL_MATCH_VALUE=$(q "$signal_match_value")
SIGNAL_LOCUS_POS=$(q "$signal_locus_pos")
SIGNAL_DISPLAY_COL=$(q "$signal_display_col")
PGS_STEPS=$(q "$step")
LEAD_REF_CACHE=$(q "$dir_out/.project/$category/lead_reference")
MH_PLOT_BED=$(q "$gwas_mh_plot_bed")
HM3=$(q "$hm3")
HM3_POS=$(q "$gwas_hm3_pos")
P_SMALL=$(q "$p_small")
SMALL_MODE=$(q "$small_mode")
MH_PNG=$(q "$mh_png")
MH_META=$(q "$mh_meta")
MH_FLAG=$(q "$mh_flag")
MH_SIG=$(q "$mh_sig")
COJO_FILE=$(q "${merged_prefix}.jma.cojo")
MPLOT_FLAG_FILE=$(q "$mplot_flag_file")
MAGMA_DIR=$(q "$magma_dir")
MAGMA_PREFIX=$(q "$magma_prefix")
MAGMA_ANNOT_CACHE=$(q "$magma_annot_cache")
DELETE_RAW=$(q "$delete_raw")

DO_STEP=$(q "$step")
REPLACE=$(q "$replace")
FILL_EAF=$(q "$fill_eaf")
FILL_N=$(q "$fill_n")
DO_LIFTOVER=$(q "$liftOver")
CHAIN=$(q "$chain")
LIFTOVER_BIN=$(q "$liftover_bin")

CIS_BED=$(q "$cis_bed")
CIS_FLANK=$(q "$cis_flank")

P_LEAD=$(q "$p_lead")
LEAD_WINDOW=$(q "$lead_window")
CLUMP_KB=$(q "$clump_kb")
CHRS=$(q "$chrs")
REFGEN_CLUMP=$(q "$gwas_refGen_clump")
REFGEN_POP=$(q "$refGen_pop")
REFGEN_KEEP=$(q "$gwas_refGen_keep")
REFGEN_COJO=$(q "$gwas_refGen_cojo")

GRCH=$(q "$gwas_grch")
MAGMA_REF=$(q "$magma_ref")
GENE_LOC=$(q "$gwas_gene_loc")
SYNONYMS=$(q "$synonyms")
MAGMA_WINDOW=$(q "$magma_window")
MAGMA_N=$(q "$magma_N")
GWAS_N=$(q "$gwas_N")
PGS_DIR=$(q "$pgs_dir")
PGS_SCORE_FILE=$(q "$pgs_score_file")
PGS_OUTPUT=$(q "$pgs_output")
PGS_DONE=$(q "$pgs_done")
PGS_META=$(q "$pgs_meta")
PGS_PFILE_DIR=$(q "${gwas_pgs_pfile_dir%/}")
PGS_THREADS=$(q "$pgs_threads")

# Do not delete trait-wide error files here: another module may be using the
# same trait concurrently.  Each coordinator now owns a module-specific error.

if [[ ",\$DO_STEP," == *,format,* || "\$DO_STEP" == "all" ]]; then
  if [[ ! -s "\$RAW" ]]; then
    echo "ERROR: missing/empty raw GWAS for \$GWAS: \$RAW" >&2
    echo "ERROR: removing incomplete outputs for \$GWAS before re-run." >&2
    find "\$GWAS_DIR" -mindepth 1 -depth ! -path "\$GWAS_DIR/\$GWAS.log" -delete 2>/dev/null || true
    rm -f -- "\${QC_PREFIX}"* "\$MH_PNG" "\$MH_META" "\${MH_PNG}.meta.tsv" "\$MH_SIG"
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
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); EAF_col=\$(gwas_clean_col "\${EAF_col:-}")
  N_col=\$(gwas_clean_col "\${N_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  LOG10P_col=\$(gwas_clean_col "\${LOG10P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  input_fs=\$(gwas_clean_detect_fs "\$src"); tmp="\${out}.tmp.\$\$"
  gwas_post_zcat "\$src" | awk -v FS="\$input_fs" -v OFS='\t' \
    -v snp_col="\$SNP_col" -v chr_col="\$CHR_col" -v pos_col="\$POS_col" \
    -v ea_col="\$EA_col" -v nea_col="\$NEA_col" -v eaf_col="\$EAF_col" -v n_col="\$N_col" \
    -v beta_col="\$BETA_col" -v se_col="\$SE_col" -v p_col="\$P_col" -v logp_col="\$LOG10P_col" '
    function get(c, x){x=(c>0 ? \$c : "");gsub(/\r/,"",x);return x}
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
    gwas_clean_compress > "\$tmp"
  gzip -t "\$tmp"; mv -f "\$tmp" "\$out"
  gwas_post_zcat "\$out" | awk -v g="\$GWAS" 'NR==1{next} END{print g"\t"NR-1}' > "\${QC_PREFIX}.format.nrow.tsv"
}

gwas_post_validate_format_columns(){
  local src="\$1" header_out="\${QC_PREFIX}.header.required.txt"
  gwas_clean_header_names "\$src" "\$header_out"
  SNP_col=\$(gwas_clean_col "\${SNP_col:-}"); CHR_col=\$(gwas_clean_col "\${CHR_col:-}")
  POS_col=\$(gwas_clean_col "\${POS_col:-}"); EA_col=\$(gwas_clean_col "\${EA_col:-}")
  NEA_col=\$(gwas_clean_col "\${NEA_col:-}"); BETA_col=\$(gwas_clean_col "\${BETA_col:-}")
  SE_col=\$(gwas_clean_col "\${SE_col:-}"); P_col=\$(gwas_clean_col "\${P_col:-}")
  [[ "\$SNP_col" -gt 0 && "\$CHR_col" -gt 0 && "\$POS_col" -gt 0 && "\$EA_col" -gt 0 && "\$NEA_col" -gt 0 && "\$BETA_col" -gt 0 && "\$SE_col" -gt 0 && "\$P_col" -gt 0 ]] || {
    echo "ERROR: required raw GWAS columns are SNP, CHR, POS, EA, NEA, BETA, SE, and P: \$src" >&2
    exit 1
  }
  gwas_post_log "format column map: SNP=\$SNP_col CHR=\$CHR_col POS=\$POS_col EA=\$EA_col NEA=\$NEA_col BETA=\$BETA_col SE=\$SE_col P=\$P_col"
}

gwas_post_fill_missing_fields(){
  local target="\$1" tmp eaf_tmp
  local -a eaf_keep_args=()
  [[ -s "\$target" ]] || { echo "ERROR: formatted GWAS is missing: \$target" >&2; return 1; }

  if [[ -n "\$FILL_N" ]]; then
    tmp="\${target}.fill_n.tmp.\$\$"
    gwas_post_log "fill missing/invalid N with \$FILL_N: \$target"
    gwas_post_zcat "\$target" | awk -v FS='\t' -v OFS='\t' -v fill_n="\$FILL_N" -v audit="\${QC_PREFIX}.fill_n.tsv" '
      function isnum(x){return x ~ /^[0-9]+([.][0-9]+)?\$/ && x+0>0}
      NR==1{for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}if(!("N" in c)){print "ERROR: standardized GWAS lacks N" > "/dev/stderr";exit 2}print;next}
      {if(isnum(\$(c["N"])))kept++;else{\$(c["N"])=fill_n;filled++}print}
      END{print "STATUS\tN" > audit;print "existing\t" kept+0 >> audit;print "filled\t" filled+0 >> audit}
    ' | gzip -c > "\$tmp"
    gzip -t "\$tmp"; mv -f "\$tmp" "\$target"
  fi

  if [[ "\$FILL_EAF" == TRUE ]]; then
    declare -F match_EAF >/dev/null 2>&1 || gwas_clean_load_phef
    eaf_tmp="\${target}.fill_eaf.tmp.\$\$.gz"
    gwas_post_log "fill missing/invalid EAF from 1KG GRCh\$GRCH: \$target"
    [[ -z "\$REFGEN_KEEP" ]] || eaf_keep_args=(--keep "\$REFGEN_KEEP")
    match_EAF --reference "\$REFGEN_CLUMP" "\${eaf_keep_args[@]}" --output "\$eaf_tmp" \
      --audit "\${QC_PREFIX}.fill_eaf.tsv" "\$target"
    gzip -t "\$eaf_tmp"; mv -f "\$eaf_tmp" "\$target"
  fi
}

if [[ "\$SMALL_MODE" == "FALSE" ]]; then
  std_small(){ std_format "\$@"; }
fi

gwas_post_prune_magma_dir(){
  [[ -d "\$MAGMA_DIR" ]] || return 0
  find "\$MAGMA_DIR" -mindepth 1 -maxdepth 1 -type f \
    ! -name '*.genes.out' ! -name '*.genes.raw' ! -name '*.log' \
    ! -name 'magma.meta.tsv' ! -name 'magma.done' -delete
}

gwas_post_magma_annotation(){
  local snploc="\${1:-}" cache_key window_tag resource_tag cache_dir annot meta
  local lock_file lock_fd tmp_prefix annot_tmp meta_tmp nloc snploc_hash
  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock not found; required for the shared MAGMA annotation cache" >&2; exit 1; }

  window_tag=\$(printf '%s' "\$MAGMA_WINDOW" | tr -c 'A-Za-z0-9._-' '_')
  resource_tag=\$(printf '%s\n%s\n' "\$MAGMA_REF" "\$GENE_LOC" | sha256sum | awk '{print substr(\$1,1,16)}')
  cache_key="v2.GRCh\${GRCH}.window_\${window_tag}.magma_\${resource_tag}"
  cache_dir="\$MAGMA_ANNOT_CACHE/v2/GRCh\$GRCH/window_\$window_tag/magma_\$resource_tag"
  annot="\$cache_dir/genes.annot"
  meta="\$cache_dir/annotation.meta.tsv"
  mkdir -p "\$(dirname "\$cache_dir")"

  lock_file="\${cache_dir}.lock"
  exec {lock_fd}> "\$lock_file"
  flock "\$lock_fd"
  if [[ ! -s "\$annot" || ! -s "\$meta" ]]; then
    mkdir -p "\$cache_dir"
    if [[ -z "\$snploc" ]]; then
      snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"
      gwas_post_log "Prepare MAGMA annotation coordinates from \$FINAL"
      gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
        NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
        {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
         if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
        sort -T "\$GWAS_POST_TMP" -k1,1 -k2,2n -k3,3n -u > "\$snploc"
    fi
    gwas_post_need_file "\$snploc"
    nloc=\$(wc -l < "\$snploc" | tr -d ' ')
    (( nloc > 1000 )) || { echo "ERROR: too few MAGMA annotation SNPs: \$nloc" >&2; return 1; }
    snploc_hash=\$(sha256sum "\$snploc" | awk '{print \$1}')
    tmp_prefix="\$GWAS_POST_TMP/annotation.\$\$"
    rm -f -- "\${tmp_prefix}"*
    gwas_post_log "Build shared MAGMA annotation from GWAS rsID coordinates: \$annot"
    if [[ "\$MAGMA_WINDOW" == "0,0" || "\$MAGMA_WINDOW" == 0 ]]; then
      magma --annotate --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    else
      magma --annotate window="\$MAGMA_WINDOW" --snp-loc "\$snploc" --gene-loc "\$GENE_LOC" --out "\$tmp_prefix"
    fi
    gwas_post_need_file "\$tmp_prefix.genes.annot"
    annot_tmp="\${annot}.tmp.\$\$"
    mv -f "\$tmp_prefix.genes.annot" "\$annot_tmp"
    mv -f "\$annot_tmp" "\$annot"
    meta_tmp="\${meta}.tmp.\$\$"
    {
      printf 'key\tvalue\n'
      printf 'grch\t%s\nwindow_kb\t%s\ngene_loc\t%s\nld_reference\t%s\nsource_gwas\t%s\nsnp_loc_n\t%s\nsnp_loc_sha256\t%s\n' \
        "\$GRCH" "\$MAGMA_WINDOW" "\$GENE_LOC" "\$MAGMA_REF" "\$GWAS" "\$nloc" "\$snploc_hash"
    } > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$meta"
  else
    gwas_post_log "Reuse shared MAGMA annotation: \$annot"
  fi
  flock -u "\$lock_fd"
  exec {lock_fd}>&-

  MAGMA_ANNOT="\$annot"
  MAGMA_ANNOT_KEY="\$cache_key"
  MAGMA_SNPLOC_N=\$(awk -F '\t' '\$1=="snp_loc_n"{print \$2;exit}' "\$meta")
  MAGMA_SNPLOC_HASH=\$(awk -F '\t' '\$1=="snp_loc_sha256"{print \$2;exit}' "\$meta")
}

gwas_post_magma(){
  [[ ",\$DO_STEP," == *,magma,* ]] || return 0
  [[ -s "\$FINAL" ]] || { echo "ERROR: missing or empty file: \$FINAL" >&2; exit 1; }

  if [[ "\$GRCH" == "38" ]]; then
    echo " [\$GWAS] GRCh build 38" >&2
  else
    echo "[\$GWAS] GRCh build 37" >&2
  fi
  echo "   gene location : \$GENE_LOC" >&2
  echo "   LD reference  : \$MAGMA_REF" >&2

  if [[ "\$REPLACE" != TRUE && -s "\$MAGMA_DIR/magma.done" && -s "\$MAGMA_PREFIX.genes.out" && -s "\$MAGMA_PREFIX.genes.raw" ]]; then
    gwas_post_prune_magma_dir
    echo "[\$(date '+%F %T')] MAGMA exists: \$MAGMA_PREFIX.genes.out (GRCh\$GRCH)" >&2
    return 0
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
  snploc="\$GWAS_POST_TMP/\$GWAS.snp.loc"; pval="\$GWAS_POST_TMP/\$GWAS.pval"
  gwas_clean_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' '
    NR==1{for(i=1;i<=NF;i++)c[\$i]=i;next}
    {s=\$(c["SNP"]);ch=\$(c["CHR"]);pos=\$(c["POS"]);gsub(/^chr/,"",ch)
     if(s!=""&&s!="NA"&&s!="."&&ch~/^([1-9]|1[0-9]|2[0-3])\$/&&pos~/^[0-9]+\$/)print s,ch,pos}' |
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
  gwas_post_magma_annotation "\$snploc"
  magma --bfile "\$MAGMA_REF" synonyms="\$SYNONYMS" --pval "\$pval" "\$narg" --gene-annot "\$MAGMA_ANNOT" --out "\$MAGMA_PREFIX"
  [[ -s "\$MAGMA_PREFIX.genes.out" ]] || { echo "ERROR: MAGMA output missing: \$MAGMA_PREFIX.genes.out" >&2; exit 1; }
  [[ -s "\$MAGMA_PREFIX.genes.raw" ]] || { echo "ERROR: MAGMA intermediate gene result missing: \$MAGMA_PREFIX.genes.raw" >&2; exit 1; }
  meta_tmp="\$MAGMA_DIR/magma.meta.tsv.tmp.\$\$"
  { printf 'key\tvalue\n'; printf 'gwas\t%s\ngrch\t%s\ngene_loc\t%s\nld_reference\t%s\nsynonyms\t%s\nwindow_kb\t%s\nsnp_loc_n\t%s\npval_n\t%s\nannotation_cache\t%s\nannotation_key\t%s\nsnp_loc_sha256\t%s\n' \
      "\$GWAS" "\$GRCH" "\$GENE_LOC" "\$MAGMA_REF" "\$SYNONYMS" "\$MAGMA_WINDOW" "\$nloc" "\$npval" "\$MAGMA_ANNOT" "\$MAGMA_ANNOT_KEY" "\$MAGMA_SNPLOC_HASH"; } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$MAGMA_DIR/magma.meta.tsv"
  gwas_post_prune_magma_dir
  date '+%F %T' > "\$MAGMA_DIR/magma.done"
  echo "[\$(date '+%F %T')] MAGMA done: \$MAGMA_PREFIX.genes.out" >&2
}

gwas_post_mplot(){
  [[ ",\$DO_STEP," == *,mplot,* || "\$DO_STEP" == "all" ]] || return 0
  gwas_post_need_file "\$FINAL"
  if [[ "\$REPLACE" != "TRUE" && -s "\$MH_PNG" ]]; then
    gwas_post_log "Manhattan plot exists: \$MH_PNG"
    return 0
  fi
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found; required for Manhattan plotting" >&2; exit 1; }
  mkdir -p "\$(dirname "\$MH_PNG")"
  plot_input="\$FINAL"
  if [[ "\$SMALL_MODE" == "FALSE" ]]; then
    plot_input="\${QC_PREFIX}.mplot.small.gz"
    gwas_post_log "Subset Manhattan input to HM3 or P < 0.001: \$plot_input"
    gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v hm3_file="\$HM3" -v pthr='0.001' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{while((getline x<hm3_file)>0){split(x,a,/[ \t]+/); if(a[1]!=""&&a[1]!="SNP")hm3[a[1]]=1} close(hm3_file)}
      NR==1{for(i=1;i<=NF;i++)c[\$i]=i; print; next}
      (("SNP" in c)&&\$(c["SNP"]) in hm3) || (("P" in c)&&isnum(\$(c["P"]))&&\$(c["P"])+0<pthr){print}' | gwas_clean_compress > "\$plot_input"
    gzip -t "\$plot_input"
  fi
  gwas_post_log "Manhattan plot: \$plot_input -> \$MH_PNG"
  Rscript - "\$plot_input" "\$MH_PNG" "\$GWAS" "\$PLOT_F" "\$CIS_BED" "\$MH_PLOT_BED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
output <- args[[2]]
gwas <- args[[3]]
plot_f <- args[[4]]
cis_bed <- if (nzchar(args[[5]])) args[[5]] else NULL
mh_plot_bed <- if (nzchar(args[[6]])) args[[6]] else NULL
source(plot_f)

cis_gene_arg <- NULL
if (!is.null(cis_bed) && file.exists(cis_bed)) {
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

gwas_post_pgs_current(){
  [[ -s "\$PGS_OUTPUT" && -s "\$PGS_DONE" ]]
}

gwas_post_prune_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    rm -f -- "\$f"
  done
}

gwas_post_prune_empty_pgs_dir(){
  local f
  [[ -d "\$PGS_DIR" ]] || return 0
  for f in "\$PGS_DIR/\${GWAS}.chr"*; do
    [[ -f "\$f" ]] || continue
    case "\$f" in *.log) continue;; esac
    rm -f -- "\$f"
  done
}

gwas_post_pgs(){
  [[ ",\$PGS_STEPS," == *,pgs,* || "\$PGS_STEPS" == "all" ]] || return 0
  if [[ "\$REPLACE" != "TRUE" ]] && gwas_post_pgs_current; then
    if [[ -s "\$PGS_META" ]] && awk -F '\t' '\$1=="matched_variants"&&\$2=="none"{found=1} END{exit !found}' "\$PGS_META"; then
      gwas_post_prune_empty_pgs_dir
    else
      gwas_post_prune_pgs_dir
    fi
    gwas_post_log "PGS exists and is current: \$PGS_OUTPUT"
    return 0
  fi
  gwas_post_need_file "\$PGS_SCORE_FILE"
  gwas_clean_load_phef
  declare -F pgs_plink_calc >/dev/null 2>&1 || { echo "ERROR: pgs_plink_calc is missing from \$PHEF" >&2; exit 1; }
  mkdir -p "\$PGS_DIR"
  rm -f -- "\$PGS_DONE"

  local status_file="\$GWAS_POST_TMP/\$GWAS.pgs.status.tsv"
  local pgs_match_status score_chrs intermediate_status meta_tmp done_tmp done_status=complete
  pgs_plink_calc \
    --input "\$PGS_SCORE_FILE" \
    --pfile-dir "\$PGS_PFILE_DIR" \
    --output "\$PGS_OUTPUT" \
    --label "\$GWAS" \
    --work-dir "\$PGS_DIR" \
    --status-file "\$status_file" \
    --threads "\$PGS_THREADS"
  pgs_match_status=\$(awk -F '\t' '\$1=="matched_variants"{print \$2}' "\$status_file")
  score_chrs=\$(awk -F '\t' '\$1=="chromosomes"{print \$2}' "\$status_file")
  intermediate_status=\$(awk -F '\t' '\$1=="chromosome_intermediates"{print \$2}' "\$status_file")
  [[ "\$pgs_match_status" == scored || "\$pgs_match_status" == none ]] || { echo "ERROR: invalid PGS status: \$pgs_match_status" >&2; exit 1; }
  [[ -n "\$score_chrs" && -n "\$intermediate_status" ]] || { echo "ERROR: incomplete PGS status file: \$status_file" >&2; exit 1; }
  [[ "\$pgs_match_status" != none ]] || done_status=complete_no_matched_variants

  meta_tmp="\${PGS_META}.tmp.\$\$"
  {
    printf 'key\tvalue\n'
    printf 'gwas\t%s\nsource\t%s\neffect_allele\trefA\nscore_weight\tbJ\nscore_stat\tdosage_weighted_sum\nmissing_genotype\tno_mean_imputation\nmatched_variants\t%s\npfile_dir\t%s\nchromosomes\t%s\noutput\t%s\nchromosome_intermediates\t%s\n' \
      "\$GWAS" "\$PGS_SCORE_FILE" "\$pgs_match_status" "\$PGS_PFILE_DIR" "\$score_chrs" "\$PGS_OUTPUT" "\$intermediate_status"
  } > "\$meta_tmp"
  mv -f "\$meta_tmp" "\$PGS_META"
  done_tmp="\${PGS_DONE}.tmp.\$\$"
  printf 'GWAS\tSTATUS\tTIME\n%s\t%s\t%s\n' "\$GWAS" "\$done_status" "\$(date '+%F %T')" > "\$done_tmp"
  mv -f "\$done_tmp" "\$PGS_DONE"
  if [[ "\$pgs_match_status" == none ]]; then
    gwas_post_log "PGS done with no matched variants: \$PGS_OUTPUT (header only; PLINK2 logs retained)"
  else
    gwas_post_log "PGS done: \$PGS_OUTPUT (per-chromosome working files removed from \$PGS_DIR)"
  fi
}

# Load the index and performance overrides after all legacy definitions so
# the optimized implementations replace them deliberately.
gwas_post_magma_annotation_impl=\$(declare -f gwas_post_magma_annotation)
source "\$INDEX_F"
source "\$PERF_F"
# The performance helper's MAGMA override derives annotation IDs from
# REFGEN_CLUMP.  Those IDs are not guaranteed to match MAGMA_REF (for example,
# GRCh38 pvars use CHR:POS:REF:ALT while g1000_eur.bim uses rsIDs).  Keep the
# optimized p-value preparation, but restore the compatible annotation builder.
eval "\$gwas_post_magma_annotation_impl"
unset gwas_post_magma_annotation_impl

if [[ "\$DO_STEP" == "all" ]]; then
  gwas_post_validate_format_columns "\$RAW"
  DO_STEP=small; gwas_clean_run_core
  gwas_post_fill_missing_fields "\$SMALL"
  gwas_post_ensure_index "\$SMALL"
  if [[ "\$DO_LIFTOVER" == TRUE ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP=all
  gwas_post_ensure_index "\$FINAL"
  gwas_post_prepare_views
  gwas_post_mplot
  [[ -z "\$CIS_BED" ]] || gwas_clean_make_cis
else
  REQUESTED_STEP="\$DO_STEP"
  if [[ ",\$REQUESTED_STEP," == *,format,* ]]; then
    gwas_post_validate_format_columns "\$RAW"
    DO_STEP=small; gwas_clean_run_core
    gwas_post_fill_missing_fields "\$SMALL"
    gwas_post_ensure_index "\$SMALL"
  fi
  if [[ ",\$REQUESTED_STEP," == *,liftover,* ]]; then
    DO_STEP=liftover; gwas_clean_run_core
  fi
  DO_STEP="\$REQUESTED_STEP"
  gwas_post_ensure_index "\$FINAL"
  gwas_post_prepare_views
  gwas_post_magma
  gwas_post_mplot
  if [[ ",\$REQUESTED_STEP," == *,cis,* ]]; then gwas_clean_make_cis; fi
  if [[ ",\$REQUESTED_STEP," == *,lead,* ]]; then DO_STEP=lead; else DO_STEP=none; fi
fi

gwas_post_clump_complete(){
  [[ -s "\$CLUMP_DONE" ]] && { [[ -s "\${MERGED}.clumps" ]] || gwas_post_lead_awk_has_no_rows; }
}

gwas_post_cojo_complete(){
  [[ -s "\$COJO_DONE" ]] && {
    [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]] || gwas_post_lead_awk_has_no_rows
  }
}

gwas_post_lead_awk_has_no_rows(){
  [[ -s "\$AWK_SNP" ]] || return 1
  awk 'NR>1{found=1} END{exit found ? 1 : 0}' "\$AWK_SNP"
}

gwas_post_lead_done(){
  [[ "\$REPLACE" != "TRUE" ]] || return 1
  [[ -s "\$AWK_SNP" ]] || return 1
  [[ ! -d "\$CLUMP" && ! -d "\$COJO" ]] || return 1
  gwas_post_clump_complete && gwas_post_cojo_complete
}

gwas_post_mark_phase_done(){
  local phase="\$1" status="\${2:-complete}" marker tmp
  case "\$phase" in
    clump) marker="\$CLUMP_DONE" ;;
    cojo) marker="\$COJO_DONE" ;;
    *) echo "ERROR: unknown lead phase: \$phase" >&2; return 2 ;;
  esac
  tmp="\${marker}.tmp.\$\$"
  mkdir -p "\$(dirname "\$marker")"
  {
    printf 'GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "\$GWAS" "\$phase" "\$P_LEAD" "\$LEAD_WINDOW" "\$CHRS" "\$status"
  } > "\$tmp"
  mv -f "\$tmp" "\$marker"
}

# Exclude ambiguous PLINK IDs. Reuse the build/chromosome output until deleted.
gwas_post_filter_clump_multiallelic(){
  local ref="\$1" tag="\$2" assoc="\$3" pvar pvar_real bad_dir bad ref_meta tmp meta_tmp all_ids multi_ids dup_ids filtered n_before n_after n_removed
  if [[ -s "\${ref}.pvar.zst" ]]; then
    pvar="\${ref}.pvar.zst"
  else
    pvar="\${ref}.pvar"
  fi
  gwas_post_need_file "\$pvar"
  pvar_real=\$(readlink -f -- "\$pvar")
  [[ -n "\$pvar_real" ]] || { echo "ERROR: cannot resolve reference pvar: \$pvar" >&2; exit 1; }
  bad_dir="\$REF_BAD_DIR/GRCh\${GRCH}"
  bad="\$bad_dir/\${tag}.ambiguous.snp"
  ref_meta="\${bad}.reference.tsv"
  mkdir -p "\$bad_dir"

  command -v flock >/dev/null 2>&1 || { echo "ERROR: flock is required for shared lead reference caches" >&2; exit 1; }
  lock_file="\${bad}.lock"
  exec {bad_lock_fd}> "\$lock_file"
  flock "\$bad_lock_fd"
  if [[ ! -e "\$bad" ]]; then
    tmp="\${bad}.tmp.\$\$"
    meta_tmp="\${ref_meta}.tmp.\$\$"
    all_ids="\${tmp}.all"
    multi_ids="\${tmp}.multi"
    dup_ids="\${tmp}.dup"
    : > "\$multi_ids"
    : > "\$dup_ids"
    gwas_clean_zcat "\$pvar" | awk -v FS='\t' -v multi="\$multi_ids" '
      \$1=="#CHROM" {for(i=1;i<=NF;i++){x=\$i; sub(/^#/,"",x); if(x=="ID")id=i; if(x=="ALT")alt=i} next}
      /^##/ {next}
      id>0 && \$id!="" && \$id!="." {print \$id; if(alt>0 && index(\$alt,",")>0) print \$id > multi}
    ' > "\$all_ids"
    sort "\$all_ids" | uniq -d > "\$dup_ids"
    cat "\$multi_ids" "\$dup_ids" 2>/dev/null | sort -u > "\$tmp"
    mv -f "\$tmp" "\$bad"
    printf 'GRCH\tPVAR\n%s\t%s\n' "\$GRCH" "\$pvar_real" > "\$meta_tmp"
    mv -f "\$meta_tmp" "\$ref_meta"
    rm -f "\$all_ids" "\$multi_ids" "\$dup_ids"
  fi
  flock -u "\$bad_lock_fd"
  exec {bad_lock_fd}>&-

  n_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  filtered="\${assoc}.filtered.\$\$"
  awk -v FS='\t' -v OFS='\t' '
    FILENAME==ARGV[1] {bad[\$1]=1; next}
    FNR==1 {print; next}
    !(\$1 in bad) {print}
  ' "\$bad" "\$assoc" > "\$filtered"
  mv -f "\$filtered" "\$assoc"
  n_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  n_removed=\$((n_before-n_after))
  printf 'GWAS\tCHR\tN_BEFORE\tN_REMOVED_AMBIGUOUS_REF_ID\tN_AFTER\n%s\t%s\t%s\t%s\t%s\n' \
    "\$GWAS" "\$tag" "\$n_before" "\$n_removed" "\$n_after" > "\${QC_PREFIX}.\${tag}.clump_multiallelic.tsv"
  gwas_post_log "clump ambiguous-reference-ID filter \$tag: before=\$n_before removed=\$n_removed after=\$n_after cache=GRCh\${GRCH}/\${tag}.ambiguous.snp"
}

gwas_post_prep_lead_inputs(){
  local refs="\$1" suffix=".tmp.\$\$" tag assoc ma
  rm -f -- "\${QC_PREFIX}".*.cojo_skip.log
  while IFS=\$'\t' read -r _ tag _ assoc ma _; do
    rm -f "\${assoc}\${suffix}" "\${ma}\${suffix}"
  done < "\$refs"
  gwas_post_zcat "\$FINAL" | awk -v FS='\t' -v OFS='\t' -v refs="\$refs" \
    -v suffix="\$suffix" -v default_n="\$GWAS_N" -v p_lead="\$P_LEAD" '
    function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
    function valid(x){return x!="" && x!="NA" && x!="."}
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;return x}
    function get(name){return name in c ? \$(c[name]) : "NA"}
    BEGIN{
      while((getline line < refs)>0){split(line,a,"\t");chr=a[3];if(chr=="")continue;
        assoc[chr]=a[4] suffix;ma[chr]=a[5] suffix
        print "SNP","CHR","POS","EA","NEA","P" > assoc[chr]
        print "SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","N" > ma[chr]}
      close(refs)
    }
    NR==1{
      for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);c[h]=i}
      required[1]="SNP";required[2]="CHR";required[3]="POS";required[4]="EA"
      required[5]="NEA";required[6]="BETA";required[7]="SE";required[8]="P"
      for(i=1;i<=8;i++)if(!(required[i] in c)){print "ERROR: required GWAS column is missing: " required[i] > "/dev/stderr";fatal=1}
      if(fatal)exit 2
      next
    }
    {
      snp=get("SNP");chr=normchr(get("CHR"));pos=get("POS");ea=get("EA");nea=get("NEA")
      beta=get("BETA");se=get("SE");p=get("P");eaf=get("EAF");n=get("N")
      # Filter significant variants before the in-memory ID match.
      if(!(chr in assoc)||pos!~/^[0-9]+\$/||pos+0<=0||!valid(snp)||!valid(ea)||!valid(nea)||!isnum(p)||p+0>p_lead)next
      print snp,chr,pos,ea,nea,p > assoc[chr]
      if(!isnum(beta)||!isnum(se))next
      if(!isnum(n))n=default_n
      print snp,chr,pos,ea,nea,eaf,beta,se,p,n > ma[chr]
    }
  '
  while IFS=\$'\t' read -r _ _ _ assoc ma _; do
    mv -f "\${assoc}\${suffix}" "\$assoc";mv -f "\${ma}\${suffix}" "\$ma"
  done < "\$refs"
}

gwas_post_subset_match_reference(){
  local ref="\$1" query="\$2" out="\$3" kind="\$4" ref_file="" ext
  case "\$kind" in
    pvar)
      case "\$ref" in
        *.pvar|*.pvar.gz|*.pvar.bgz|*.pvar.zst) ref_file="\$ref" ;;
        *)
          for ext in .pvar .pvar.gz .pvar.bgz .pvar.zst; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    bim)
      case "\$ref" in
        *.bim|*.bim.gz|*.bim.bgz) ref_file="\$ref" ;;
        *)
          for ext in .bim .bim.gz .bim.bgz; do
            [[ -s "\${ref}\${ext}" ]] && { ref_file="\${ref}\${ext}"; break; }
          done
          ;;
      esac
      ;;
    *) echo "ERROR: invalid reference subset kind: \$kind" >&2; return 2 ;;
  esac
  gwas_post_need_file "\$ref_file"

  # Restrict the reference to query coordinates before match_SNP.
  awk -v FS='[ \t]+' -v kind="\$kind" '
    function normchr(x){gsub(/^chr/,"",x);x=toupper(x);if(x=="X")x=23;else if(x=="Y")x=24;else if(x=="MT"||x=="M")x=25;sub(/^0+/,"",x);return x=="" ? "0" : x}
    NR==FNR{
      if(FNR==1){for(i=1;i<=NF;i++){h=toupper(\$i);sub(/^#/,"",h);if(h=="CHR")qchr=i;else if(h=="POS")qpos=i}next}
      if(!qchr||!qpos){fatal=1;next}
      if(\$(qpos)~/^[0-9]+\$/)wanted[normchr(\$(qchr)) SUBSEP \$(qpos)]=1
      next
    }
    kind=="pvar" && /^##/{print;next}
    kind=="pvar" && /^#CHROM/{print;header=1;next}
    kind=="pvar"{
      if((normchr(\$1) SUBSEP \$2) in wanted){print;kept++}
      next
    }
    kind=="bim"{
      if(NF>=4 && (normchr(\$1) SUBSEP \$4) in wanted){print;kept++}
      next
    }
    END{
      if(fatal){print "ERROR: lead query requires CHR and POS columns" > "/dev/stderr";exit 2}
      if(kind=="pvar"&&!header){print "ERROR: invalid .pvar header" > "/dev/stderr";exit 2}
      # Keep an empty BIM subset non-empty so match_SNP can report all query
      # rows as unmatched instead of rejecting the reference argument.
      if(kind=="bim"&&kept==0)print "0\t.\t0\t0\tN\tN"
    }
  ' "\$query" <(gwas_clean_zcat "\$ref_file") > "\$out"
}

# GCTA treats only numeric chromosomes up to --autosome-num as usable.  PLINK
# references commonly encode chrX/chrY as X/Y, which makes GCTA silently load
# just one unusable record from an otherwise complete sex-chromosome BIM.  For
# sex chromosomes, make a small per-trait BED containing only the already
# matched COJO candidates and emit numeric chromosome codes (X=23, Y=24).
gwas_post_prepare_gcta_bfile(){
  local source="\$1" chr="\$2" ma="\$3" out="\$4"
  GCTA_BFILE="\$source"
  GCTA_CHR_ARGS=()
  (( chr > 22 )) || return 0
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 is required to prepare chr\$chr for GCTA" >&2; return 1; }
  gwas_post_need_file "\${source}.bed"
  gwas_post_need_file "\${source}.bim"
  gwas_post_need_file "\${source}.fam"
  gwas_post_need_file "\$ma"
  rm -f -- "\${out}.bed" "\${out}.bim" "\${out}.fam" "\${out}.log" "\${out}.err"
  run_tool "\$out" plink2 --bfile "\$source" --extract "\$ma" --make-bed --output-chr 26 --out "\$out"
  gwas_post_need_file "\${out}.bed"
  gwas_post_need_file "\${out}.bim"
  gwas_post_need_file "\${out}.fam"
  GCTA_BFILE="\$out"
  GCTA_CHR_ARGS=(--autosome-num "\$chr")
}

gwas_post_match_lead_inputs(){
  local ref="\$1" cref="\$2" tag="\$3" assoc="\$4" ma="\$5"
  local matched ref_subset n_assoc_before n_assoc_after n_ma_before n_ma_after
  declare -F match_SNP >/dev/null 2>&1 || gwas_clean_load_phef

  n_assoc_before=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  if (( n_assoc_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.clump.reference.pvar"
    gwas_post_subset_match_reference "\$ref" "\$assoc" "\$ref_subset" pvar
    matched="\${assoc}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.clump_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv" "\$assoc"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$6+0)<best[\$1]{best[\$1]=\$6+0;line[\$1]=\$1 OFS \$6}
      END{print "SNP","P";for(snp in line)print line[snp]}' "\$matched" > "\${assoc}.tmp.\$\$"
    mv -f "\${assoc}.tmp.\$\$" "\$assoc"; rm -f "\$matched"
    n_assoc_after=\$(awk 'NR>1{n++} END{print n+0}' "\$assoc")
  else
    rm -f "\${QC_PREFIX}.\${tag}.clump_match.tsv" "\${QC_PREFIX}.\${tag}.clump_unmatched.tsv"
    n_assoc_after=0
  fi

  n_ma_before=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  if (( n_ma_before > 0 )); then
    ref_subset="\$GWAS_POST_TMP/\${tag}.cojo.reference.bim"
    gwas_post_subset_match_reference "\$cref" "\$ma" "\$ref_subset" bim
    matched="\${ma}.matched.\$\$"
    match_SNP --reference "\$ref_subset" --output "\$matched" \
      --audit "\${QC_PREFIX}.\${tag}.cojo_match.tsv" \
      --unmatched "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv" "\$ma"
    awk -v FS='\t' -v OFS='\t' '
      NR==1{next}
      !((\$1) in best) || (\$9+0)<best[\$1]{best[\$1]=\$9+0;line[\$1]=\$1 OFS \$4 OFS \$5 OFS \$6 OFS \$7 OFS \$8 OFS \$9 OFS \$10}
      END{print "SNP","A1","A2","freq","b","se","p","N";for(snp in line)print line[snp]}' "\$matched" > "\${ma}.tmp.\$\$"
    mv -f "\${ma}.tmp.\$\$" "\$ma"; rm -f "\$matched"
    n_ma_after=\$(awk 'NR>1{n++} END{print n+0}' "\$ma")
  else
    rm -f "\${QC_PREFIX}.\${tag}.cojo_match.tsv" "\${QC_PREFIX}.\${tag}.cojo_unmatched.tsv"
    n_ma_after=0
  fi

  gwas_post_log "reference-ID match \$tag: clump=\$n_assoc_after/\$n_assoc_before cojo=\$n_ma_after/\$n_ma_before"
}

if run_lead; then
  if gwas_post_lead_done; then
    gwas_post_log "lead SNP discovery exists: \$AWK_SNP"
  else
    lead_t0=\$(date +%s)
    clump_was_done=FALSE
    cojo_was_done=FALSE
    cojo_skipped=FALSE
    if [[ "\$REPLACE" != "TRUE" ]]; then
      [[ -s "\$CLUMP_DONE" && -s "\${MERGED}.clumps" ]] && clump_was_done=TRUE
      if [[ -s "\$COJO_DONE" ]] && [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]]; then
        cojo_was_done=TRUE
      fi
    fi
    rm -f -- "\${MERGED}.lead.done"
    if [[ "\$clump_was_done" != "TRUE" ]]; then
      rm -f -- "\$CLUMP_DONE" "\${MERGED}.clumps"
    fi
    if [[ "\$cojo_was_done" != "TRUE" ]]; then
      rm -f -- "\$COJO_DONE" "\${MERGED}.jma.cojo" "\${MERGED}.cma.cojo" "\${MERGED}.ldr.cojo"
    fi
    rm -rf -- "\$CLUMP" "\$COJO"
    mkdir -p "\$CLUMP" "\$COJO"
    gwas_post_need_file "\$FINAL"
    labels="\${CLUMP}/labels.\$\$.txt"
    refs="\${CLUMP}/refs.\$\$.tsv"
    : > "\$labels"
    : > "\$refs"

  # 1) Resolve selected reference chromosomes. GWAS coordinates, alleles, and
  # frequencies are never filled from the reference panel.
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
  gwas_post_log "prepare P<=\$P_LEAD lead inputs from required SNP/CHR/POS/EA/NEA/BETA/SE/P fields: \$FINAL"
  gwas_post_prep_lead_inputs "\$refs"

  # 2) No-LD top SNP per chromosome/window and per-chromosome tool inputs.
  gwas_post_log "awk distance lead SNPs: \${GWAS_POST_LEAD_VIEW:-\$FINAL} -> \$AWK_SNP"
  gwas_post_zcat "\${GWAS_POST_LEAD_VIEW:-\$FINAL}" | awk -v FS='\t' -v OFS='\t' -v pthr="\$P_LEAD" -v win="\$LEAD_WINDOW" -f <(awk_lead) > "\${AWK_SNP}.tmp.\$\$"
  mv -f "\${AWK_SNP}.tmp.\$\$" "\$AWK_SNP"
  awk -v g="\$GWAS" 'NR==1{next} END{print g"\\t"NR-1}' "\$AWK_SNP" > "\${QC_PREFIX}.awk.nrow.tsv"

  # 3) plink2 --clump and gcta --cojo per chromosome.
  while IFS=\$'\t' read -r ref tag chr assoc ma cref; do
    [[ -n "\$ref" ]] || continue
    cp="\${CLUMP}/\${tag}"
    jp="\${COJO}/\${tag}"
    gwas_post_match_lead_inputs "\$ref" "\$cref" "\$tag" "\$assoc" "\$ma"
    # Inputs were already restricted to P<=P_LEAD before reference matching,
    # which bounds memory independently of the full GWAS size.
    cojo_skip_log="\${QC_PREFIX}.\${tag}.cojo_skip.log"
    awk -v FS='\t' -v OFS='\t' '
      function isnum(x){return x ~ /^[-+]?([0-9]*[.])?[0-9]+([eE][-+]?[0-9]+)?\$/}
      BEGIN{print "STATUS","REASON","SNP","A1","A2","EAF","P"}
      NR>1 && (!isnum(\$4)||\$4+0<0||\$4+0>1){print "SKIP_GCTA","missing_or_invalid_EAF_after_reference_match",\$1,\$2,\$3,\$4,\$7}
    ' "\$ma" > "\${cojo_skip_log}.tmp.\$\$"
    if gwas_post_has_data_rows "\${cojo_skip_log}.tmp.\$\$"; then
      mv -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    else
      rm -f "\${cojo_skip_log}.tmp.\$\$" "\$cojo_skip_log"
    fi
    if [[ "\$clump_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP plink2 chr\$chr: clump phase already complete"
    elif gwas_post_has_data_rows "\$assoc"; then
      gwas_post_filter_clump_multiallelic "\$ref" "\$tag" "\$assoc"
      if gwas_post_has_data_rows "\$assoc"; then
        pfile_modifier=()
        [[ ! -s "\${ref}.pvar.zst" ]] || pfile_modifier=(vzs)
        pfile_keep=()
        [[ -z "\$REFGEN_KEEP" ]] || pfile_keep=(--keep "\$REFGEN_KEEP")
        run_tool "\$cp" plink2 --pfile "\$ref" "\${pfile_modifier[@]}" "\${pfile_keep[@]}" --clump "\$assoc" --clump-snp-field SNP --clump-p-field P --clump-p1 "\$P_LEAD" --clump-kb "\$CLUMP_KB" --out "\$cp"
      else
        gwas_post_log "skip plink2 chr\$chr: all significant IDs are ambiguous in the reference"
      fi
    else
      gwas_post_log "skip plink2 chr\$chr: no SNPs in \$assoc"
    fi

    # 3) gcta --cojo: COJO selection with the chromosome-specific bed/bim/fam reference
    if [[ "\$cojo_was_done" == "TRUE" ]]; then
      gwas_post_log "SKIP gcta chr\$chr: COJO phase already complete"
    elif [[ -s "\${QC_PREFIX}.\${tag}.cojo_skip.log" ]]; then
      cojo_skipped=TRUE
      gwas_post_log "SKIP gcta chr\$chr: significant variant has missing/invalid EAF; see \${QC_PREFIX}.\${tag}.cojo_skip.log"
    elif gwas_post_has_data_rows "\$ma"; then
      gwas_post_prepare_gcta_bfile "\$cref" "\$chr" "\$ma" "\${jp}.gcta_ref"
      run_tool "\$jp" gcta --bfile "\$GCTA_BFILE" "\${GCTA_CHR_ARGS[@]}" \
        --cojo-file "\$ma" --cojo-slct --cojo-p "\$P_LEAD" --out "\$jp"
    else
      gwas_post_log "skip gcta chr\$chr: no SNPs in \$ma"
    fi
  done < "\$refs"

    if gwas_post_lead_awk_has_no_rows; then
      rm -f "\$labels" "\$refs"
      rm -rf -- "\$CLUMP" "\$COJO"
      gwas_post_mark_phase_done clump no_significant_variants
      gwas_post_mark_phase_done cojo no_significant_variants
    else
      if [[ "\$clump_was_done" != "TRUE" ]]; then
        concat_chr_outputs "\$CLUMP" "\$MERGED" "\$labels" clumps
        if [[ -s "\${MERGED}.clumps" ]]; then
          gwas_post_mark_phase_done clump
        fi
      fi
      if [[ "\$cojo_was_done" != "TRUE" && "\$cojo_skipped" != "TRUE" ]]; then
        concat_chr_outputs "\$COJO" "\$MERGED" "\$labels" jma.cojo cma.cojo ldr.cojo
        if [[ -s "\${MERGED}.jma.cojo" || -s "\${MERGED}.ldr.cojo" ]]; then
          gwas_post_mark_phase_done cojo
        fi
      fi
      rm -f "\$labels" "\$refs"
      if ! gwas_post_clump_complete; then
        echo "ERROR: clump phase did not create \${MERGED}.clumps; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      if ! gwas_post_cojo_complete; then
        echo "ERROR: COJO phase is incomplete for \$GWAS; GCTA skip details are under \${QC_PREFIX}.*.cojo_skip.log; retaining \$CLUMP and \$COJO" >&2
        exit 1
      fi
      gwas_post_log "delete lead intermediates: \$CLUMP \$COJO \${MERGED}.cma.cojo"
      rm -rf -- "\$CLUMP" "\$COJO"
      rm -f -- "\${MERGED}.cma.cojo"
    fi
    gwas_post_log "lead SNP discovery done in \$((\$(date +%s)-lead_t0)) sec"
  fi
fi

# PGS depends on the genome-wide .jma.cojo created by the lead/COJO phase, so
# it is deliberately invoked after lead even though it is declared after mplot.
gwas_post_pgs

CMD_TOP

  chmod +x "$cmd"
  echo "$cmd"
}


# 🚩 checks and dispatch
need_file "$phef"
need_file "$data_f"
need_file "$index_f"
need_file "$perf_f"
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$phef"
command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip not found; install htslib" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found; install htslib" >&2; exit 1; }
if wants_magma; then
  ensure_magma_resources
  command -v magma >/dev/null 2>&1 || { echo "ERROR: magma not found in PATH ($PATH)" >&2; exit 1; }
fi
if has_step pgs; then
  command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 not found in PATH ($PATH)" >&2; exit 1; }
fi
if [[ "$step" == "all" ]]; then
  cleanup_failed_output_dirs
fi
if has_step mplot; then
  need_file "$mplot_r"
  need_file "$plot_f"
  [[ -z "$add_signal" ]] || need_file "$add_signal"
fi
if { has_step format && [[ "$small_mode" == "TRUE" ]]; } || has_step mplot; then
  need_file "$hm3"
fi
if [[ "$liftOver" == "TRUE" ]] && has_step liftover; then
  need_file "$chain"
fi
if has_step cis && [[ -z "$cis_bed" ]]; then
  echo "ERROR: --cis-bed is required for the cis module" >&2
  exit 2
fi
if [[ -n "$cis_bed" ]] && has_step cis; then
  need_file "$cis_bed"
fi
if has_step mplot; then
  # Create project-level plot destinations before the potentially long
  # per-GWAS command-generation pass, so background progress is visible at once.
  mkdir -p "$dir_out/mplot" "$dir_out/.project/$category/mplot/flag"
  [[ -z "$cis_bed" ]] || need_file "$cis_bed"
  if [[ -n "$mh_plot_bed" ]]; then
    need_file "$mh_plot_bed"
  else
    need_file "$dir0/files/glist.37.bed"
    need_file "$dir0/files/glist.38.bed"
  fi
fi
if { has_step lead || wants_magma; } && [[ "$grch" != auto ]]; then
  need_refgen_clump "$refGen_clump"
fi
if has_step lead && [[ "$grch" != auto ]]; then
  need_refgen_cojo "$refGen_cojo"
fi

log "label=$label step=$step small=$small_mode liftOver=$liftOver add-panel=$add_panel plot=${plot_width}x${plot_height}in@${plot_res}dpi write-sig=$write_sig add-signal=${add_signal:-none} match=$signal_match_col:$signal_match_value locus-pos=$signal_locus_pos display-col=$signal_display_col cis-bed=$cis_bed mh-plot-bed=${mh_plot_bed:-<auto:glist.37/glist.38>} cis-flank=$cis_flank chrs=$chrs jobs=$jobs run-cmd=$run_cmd is.bsub=$is_bsub"
log "raw=$dir_raw project=$dir_out category=$category layout=<project>/<category>/<trait>/{gwas,magma,pgs,qc} mplot=<project>/mplot coordinator-cmd=$dir_cmd refGen_clump=$refGen_clump refGen_cojo=$refGen_cojo pgs_pfile_dir=${pgs_pfile_dir:-<auto:/mnt/h/ukbGen/GRCh/imp>}"
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
log "DONE"

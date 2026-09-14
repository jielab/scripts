#!/usr/bin/env bash
# Sourced by gwas_format.sh: validation, completion checks, discovery and execution.
# Configuration comes from the entry point; sourcing only defines functions.

source "${BASH_SOURCE[0]%/*}/gwas_thin_state.f.sh"

need_arg_value() {
  local opt="$1" val="${2-}"
  if [[ -z "$val" || "$val" == --* ]]; then
    echo "ERROR: missing value for $opt" >&2
    usage >&2
    exit 2
  fi
}

upper() { echo "$1" | tr '[:lower:]' '[:upper:]'; }
has_step(){ [[ ",$step," == *",$1,"* || "$step" == "all" ]]; }
wants_magma(){ [[ ",$step," == *",magma,"* ]]; }

gwas_format_validate_options() {
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
  hm3_mode=$(upper "$hm3_mode")
  thin=$(upper "$thin")
  write_sig=$(upper "$write_sig")
  step=$(echo "$step" | tr '[:upper:]' '[:lower:]')
  category=$(echo "$category" | tr '[:upper:]' '[:lower:]')
  add_panel=$(echo "$add_panel" | tr '[:upper:]' '[:lower:]')
  plot_height=$(echo "$plot_height" | tr '[:upper:]' '[:lower:]')
  [[ "$add_panel" == none || "$add_panel" == magma ]] || { echo "ERROR: --add-panel must be none or magma: $add_panel" >&2; exit 2; }

  [[ "$hm3_mode" == "TRUE" || "$hm3_mode" == "FALSE" ]] || { echo "ERROR: --hm3 must be TRUE or FALSE; use --hm3-file FILE for the reference list" >&2; exit 2; }
  [[ "$thin" == TRUE || "$thin" == FALSE ]] || { echo "ERROR: --thin must be TRUE or FALSE" >&2; exit 2; }
  [[ "$thin_chr_max" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --thin-chr-max must be a positive integer" >&2; exit 2; }
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
  IFS=',' read -r -a requested_steps <<< "$step"
  for requested_step in "${requested_steps[@]}"; do
    case "$requested_step" in
      format|thin|magma|liftover|cis|lead|mplot|h2|pgs|all) ;;
      *) echo "ERROR: step/module must contain format|thin|magma|liftover|cis|lead|mplot|h2|pgs, or all" >&2; exit 2;;
    esac
  done

  [[ -z "$raw_file_arg" || ( -f "$raw_file_arg" && -s "$raw_file_arg" && -n "$gwas_arg" && "$gwas_arg" != *,* ) ]] || { echo 'ERROR: --raw-file requires an existing file and one --gwas NAME' >&2; exit 2; }
  [[ -z "$n_total" || "$n_total" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --n-total must be a positive integer' >&2; exit 2; }
  [[ -z "$n_total" || "$hm3_mode" == FALSE ]] || { echo 'ERROR: --n-total currently requires --hm3 FALSE' >&2; exit 2; }
  case "$h2_sex" in unknown|male|female|mixed) ;; *) echo 'ERROR: invalid --h2-sex' >&2; exit 2;; esac
  if has_step liftover && [[ "$liftOver" == TRUE ]]; then
    [[ "$step" == liftover || "$step" == format,liftover ]] || { echo 'ERROR: run format,liftover first, then downstream modules in a separate invocation with --grch auto' >&2; exit 2; }
  fi

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
  awk -v p="$p_hm3" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' || {
    echo "ERROR: --p-hm3 must be a number in (0,1]: $p_hm3" >&2
    exit 2
  }
  if has_step format && [[ "$hm3_mode" == TRUE && "$thin" == TRUE ]]; then
    awk -v p="$p_hm3" 'BEGIN{exit !(p+0>=0.001)}' || {
      echo 'ERROR: --hm3 TRUE --thin TRUE requires --p-hm3 >= 1e-3 so formatting retains all thin candidates.' >&2
      exit 2
    }
  fi
  if ! [[ "$lead_window" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --lead-window must be a positive integer: $lead_window" >&2
    exit 2
  fi

}

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

gwas_format_maybe_background() {
  if [[ "$run_cmd" == "TRUE" && "$is_bsub" != "TRUE" && "$foreground" != "TRUE" ]]; then
    background_log="$dir_cmd/gwas_post.background.log"
    [[ -s "$phef" ]] || { echo "ERROR: missing or empty file: $phef" >&2; exit 1; }
    # shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
    source "$phef"
    if declare -F phe_check_existing_tasks >/dev/null 2>&1; then
      # Only reject another run of the same module set.  Different modules use
      # separate coordinator files, generated commands, and per-module logs.
      phe_check_existing_tasks "format_gwas.sh $step"
    fi
    phe_run_background --title "background gwas_post $label/$category/$step" "$background_log" bash "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" --foreground TRUE
    exit 0
  fi

}

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

# A result is reusable only for the requested trait, phase and lead settings.
# In particular, an autosome-only marker must not satisfy all chromosomes.
gwas_lead_marker_matches() {
  local marker="$1" trait="$2" phase="$3" p="$4" window="$5" chromosomes="$6"
  [[ -s "$marker" ]] || return 1
  awk -F '\t' -v trait="$trait" -v phase="$phase" -v p="$p" \
    -v window="$window" -v chromosomes="$chromosomes" '
    NR==1 {header=($0=="GWAS\tPHASE\tP_LEAD\tLEAD_WINDOW\tCHRS\tSTATUS")}
    NR==2 {ok=(header && NF==6 && $1==trait && $2==phase &&
      $3==p && $4==window && $5==chromosomes &&
      ($6=="complete" || $6=="no_significant_variants" || $6=="no_reference_matched_variants"))}
    END {exit !(NR==2 && ok)}
  ' "$marker"
}

# Embedded in each worker. GCTA reports an empty reference-QC intersection as
# exit 1; retain that evidence without treating an empty chromosome as a crash.
gwas_post_run_cojo() {
  local prefix="$1" audit="$2" rc ext
  shift 2
  GCTA_COJO_HAS_MATCH=FALSE
  mkdir -p "$(dirname "$audit")"
  rm -f -- "${audit}.tsv" "${audit}.log" "${audit}.freq.badsnps" "${audit}.badsnps"
  rm -f -- "${prefix}.freq.badsnps" "${prefix}.badsnps"
  if run_tool "$prefix" "$@"; then
    GCTA_COJO_HAS_MATCH=TRUE
    return 0
  else
    rc=$?
  fi

  # Match only GCTA's explicit empty-intersection error after it has read the
  # summary statistics and reached reference matching. Other failures propagate.
  if [[ "$rc" == 1 ]] &&
     grep -Fxq 'Error: none of the SNPs in the GWAS summary data can be found in the genotype data.' "${prefix}.log" &&
     grep -Fxq 'Matching the GWAS meta-analysis results to the genotype data ...' "${prefix}.log" &&
     grep -Eq '^GWAS summary statistics of [1-9][0-9]* SNPs read from ' "${prefix}.log"; then
    cp -- "${prefix}.log" "${audit}.log"
    for ext in freq.badsnps badsnps; do
      if [[ -s "${prefix}.${ext}" ]]; then cp -- "${prefix}.${ext}" "${audit}.${ext}"; fi
    done
    printf 'STATUS\tREASON\tUSABLE_SNPS\nEMPTY\tno_reference_matched_variants_after_gcta_qc\t0\n' > "${audit}.tsv"
    rm -f -- "${prefix}.err" "${prefix}.jma.cojo" "${prefix}.cma.cojo" "${prefix}.ldr.cojo"
    gwas_post_log "COJO has no usable reference SNPs after GCTA QC: $prefix; evidence: ${audit}.tsv"
    return 0
  fi
  return "$rc"
}

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
  [[ -s "$magma_dir/magma.done" && -s "$magma_prefix.genes.out" && -s "$magma_prefix.genes.raw" &&
     -s "${magma_dir%/*}/qc/$(basename "$magma_prefix").magma.chromosomes.log" ]]
}

mplot_output_complete() {
  local png="$1" final="$2" genes="$3" meta="$4" expected_grch="$5" flag="$6" aggregate="$7" sig="$8" cojo="$9"
  local signal_input=none thin_output="${final%.gz}.thin.gz"
  if [[ "${thin:-FALSE}" == TRUE ]]; then
    gwas_thin_complete "$final" "$thin_output" "$expected_grch" "$thin_chr_max" "$hm3_file" "${hm3_pos//\{grch\}/$expected_grch}" "$thin_r" "$phe_r" "$meta" || return 1
    [[ "$png" -nt "$thin_output" ]] || return 1
  fi
  [[ -z "$add_signal" ]] || signal_input="$add_signal"
  [[ -s "$png" && -s "$meta" && -s "$flag" && -s "$aggregate" && "$png" -nt "$final" ]] || return 1
  awk -F '\t' -v panel="$add_panel" -v grch="$expected_grch" -v signal="$signal_input" \
    -v match_col="$signal_match_col" -v match_value="$signal_match_value" \
    -v locus_pos="$signal_locus_pos" -v display_col="$signal_display_col" -v write_sig="$write_sig" \
    -v thin="${thin:-FALSE}" -v cap="${thin_chr_max:-10000}" \
    -v plot_width="$plot_width" -v plot_height="$plot_height" -v plot_res="$plot_res" '
    $1=="plot_method"&&$2=="self"{m=1}
    $1=="add_panel"&&$2==panel{p=1}
    $1=="grch"&&$2==grch{g=1}
    $1=="magma_threshold"&&$2+0==2.5e-6{t=1}
    $1=="mplot_style"&&$2==9{s=1}
    $1=="thin_mode"&&$2==thin{tn=1}
    $1=="thin_chr_max"&&$2==cap{tc=1}
    $1=="plot_width"&&$2+0==plot_width+0{x=1}
    $1=="plot_height"&&$2==plot_height{y=1}
    $1=="plot_res"&&$2+0==plot_res+0{r=1}
    $1=="write_sig"&&$2==write_sig{w=1}
    $1=="signal_input"&&$2==signal{a=1}
    $1=="signal_match_col"&&$2==match_col{b=1}
    $1=="signal_match_value"&&$2==match_value{c=1}
    $1=="signal_locus_pos"&&$2==locus_pos{d=1}
    $1=="signal_display_col"&&$2==display_col{e=1}
    END{exit !(m&&p&&g&&t&&s&&x&&y&&r&&w&&a&&b&&c&&d&&e&&tn&&(thin!="TRUE"||tc))}' "$meta" || return 1
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
    ! -name '*.small.gz' ! -name '*.hm3.gz' ! -name '*.thin.gz' ! -name '*.cis.gz' ! -name '*.sig.tsv.gz' ! -name '*.lead.tsv.gz' \
    -size +0c 2>/dev/null |
    while read -r f; do
      [[ "$(basename "$(dirname "$f")")" == "gwas" ]] || continue
      gwas_name_from_file "$f"
    done | sort -u -V
}

list_liftover_names() {
  [[ -d "$dir_clean" ]] || return 0
  find "$dir_clean" -type f \( -name '*.hm3.gz' -o -name '*.small.gz' -o -name '*.source.grch37.gz' \) -size +0c 2>/dev/null |
    while read -r f; do
        [[ "$(basename "$(dirname "$f")")" == "gwas" || "$(basename "$(dirname "$f")")" == qc ]] || continue
        b=$(basename "$f")
        b=${b%.gz}
        b=${b%.hm3}
        b=${b%.small}
        b=${b%.source.grch37}
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
    echo "[$(date '+%F %T')] START [$step_key] $base" >&2
    rm -f "$err_file"
    if bash "$cmd" > "$log_file" 2>&1; then
      echo "[$(date '+%F %T')] DONE  [$step_key] $base" >&2
    else
      rc=$?
      {
        echo "ERROR: [$step_key] $base failed with exit=$rc"
        echo "ERROR: log=$log_file"
        while IFS= read -r tool_err; do
          [[ -s "$tool_err" ]] || continue
          echo "ERROR: detail=$tool_err"
          grep -Ei 'error|failed|invalid' "$tool_err" || true
        done < <(find "$log_dir" -mindepth 2 -type f -name '*.err' -size +0c 2>/dev/null | sort -V)
      } > "$err_file"
      echo "[$(date '+%F %T')] FAIL  [$step_key] $base exit=$rc log=$err_file" >&2
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
      log "ERROR: [$step_key] one or more command files failed. First failed jobs from $list.joblog:"
      awk -F '\t' 'NR>1 && $7 != 0 {print "  exit="$7" cmd="$9; n++; if(n>=10) exit}' "$list.joblog" >&2 || true
      return "$rc"
    fi
  else
    xargs -I{} -P "$jobs" bash -c 'run_one_cmd "$1"' _ {} < "$list" || rc=$?
    return "$rc"
  fi
}

gwas_format_check_resources() {
  need_file "$phef"
  need_file "$data_f"
  need_file "$index_f"
  need_file "$perf_f"
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
  if [[ "$thin" == TRUE ]] && { has_step format || has_step thin || has_step mplot || has_step liftover; }; then
    need_file "$thin_r"
    need_file "$phe_r"
    need_file "$hm3_file"
    command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript is required for --thin TRUE" >&2; exit 1; }
  fi
  if { has_step format && [[ "$hm3_mode" == "TRUE" ]]; } || has_step mplot; then
    need_file "$hm3_file"
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

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X_LEN_b37=155270560; X_LEN_b38=156040895

X_PAR1_START_b37=60001; X_PAR1_END_b37=2699520
X_PAR2_START_b37=154931044; X_PAR2_END_b37=155260560
X_NONPAR1_START_b37=1; X_NONPAR1_END_b37=60000
X_NONPAR2_START_b37=2699521; X_NONPAR2_END_b37=154931043
X_NONPAR3_START_b37=155260561; X_NONPAR3_END_b37=155270560

X_PAR1_START_b38=10001; X_PAR1_END_b38=2781479
X_PAR2_START_b38=155701383; X_PAR2_END_b38=156030895
X_NONPAR1_START_b38=1; X_NONPAR1_END_b38=10000
X_NONPAR2_START_b38=2781480; X_NONPAR2_END_b38=155701382
X_NONPAR3_START_b38=156030896; X_NONPAR3_END_b38=156040895

MHC_START_b37=28477897; MHC_END_b37=33448354
MHC_START_b38=28510120; MHC_END_b38=33480577

X_PAR_b37="X:${X_PAR1_START_b37}-${X_PAR1_END_b37},X:${X_PAR2_START_b37}-${X_PAR2_END_b37}"
X_NONPAR_b37="X:${X_NONPAR1_START_b37}-${X_NONPAR1_END_b37},X:${X_NONPAR2_START_b37}-${X_NONPAR2_END_b37},X:${X_NONPAR3_START_b37}-${X_NONPAR3_END_b37}"

X_PAR_b38="X:${X_PAR1_START_b38}-${X_PAR1_END_b38},X:${X_PAR2_START_b38}-${X_PAR2_END_b38}"
X_NONPAR_b38="X:${X_NONPAR1_START_b38}-${X_NONPAR1_END_b38},X:${X_NONPAR2_START_b38}-${X_NONPAR2_END_b38},X:${X_NONPAR3_START_b38}-${X_NONPAR3_END_b38}"

MHC_b37="6:${MHC_START_b37}-${MHC_END_b37}"
MHC_b38="6:${MHC_START_b38}-${MHC_END_b38}"

log() { echo "[$(date '+%F %T')] $*"; }

count_lines() { [[ -s $1 ]] && awk 'NF{n++} END{print n+0}' "$1" || echo 0; }

data_rows() { [[ -s $1 ]] && awk 'NR>1{n++} END{print n+0}' "$1" || echo 0; }

phe_background_job_pattern() {
  local x b
  for x in "$@"; do
    b="$(basename -- "$x")"
    case "$b" in
      *.sh|*.py) printf '%s\n' "$b"; return 0 ;;
    esac
  done
  b="$(basename -- "${1:-background}")"
  printf '%s\n' "$b"
}

phe_running_tasks() {
  local pattern="$1" self="${2:-$$}" bash_self="${3:-${BASHPID:-$$}}" parent="${4:-${PPID:-}}"
  ps -eo pid=,ppid=,args= 2>/dev/null | PHE_TASK_PATTERN="$pattern" awk -v self="$self" -v bash_self="$bash_self" -v parent="$parent" '
    BEGIN { pattern = ENVIRON["PHE_TASK_PATTERN"] }
    {
      pid = $1
      ppid = $2
      $1 = ""
      $2 = ""
      sub(/^[[:space:]]+/, "")
      cmd = $0
    }
    pid == self { next }
    pid == bash_self { next }
    pid == parent { next }
    ppid == self { next }
    ppid == bash_self { next }
    index(cmd, pattern) > 0 { print pid "\t" cmd }
  '
}

phe_check_existing_tasks() {
  local pattern="$1" old pids self bash_self parent
  self="$$"
  bash_self="${BASHPID:-$$}"
  parent="${PPID:-}"
  old="$(phe_running_tasks "$pattern" "$self" "$bash_self" "$parent" || true)"
  [[ -z "$old" ]] && return 0
  pids="$(printf '%s\n' "$old" | awk '{print $1}' | paste -sd ' ' -)"
  {
    echo "ERROR: existing $pattern task(s) are still running:"
    printf '%s\n' "$old" | sed 's/^/  /'
    echo
    echo "Stop them first, then rerun this command:"
    echo "  kill $pids"
    echo
    echo "If they do not exit after a few seconds, use:"
    echo "  kill -9 $pids"
  } >&2
  return 1
}

phe_run_background() {
  local title="background job" job_pattern
  if [[ "${1-}" == "--title" ]]; then
    title="$2"
    shift 2
  fi
  local log_file="$1"
  shift || {
    echo "ERROR: phe_run_background requires: LOG_FILE COMMAND [ARGS...]" >&2
    return 2
  }
  [[ $# -gt 0 ]] || {
    echo "ERROR: phe_run_background missing command" >&2
    return 2
  }
  mkdir -p -- "$(dirname -- "$log_file")"
  echo "Starting $title. Log: $log_file"
  nohup "$@" > "$log_file" 2>&1 < /dev/null &
  echo "PID: $!"
  echo "Progress: tail -f $log_file"
  job_pattern="$(phe_background_job_pattern "$@")"
  echo "Job: pgrep -af '$job_pattern'"
}

phe_zcat() {
  local f="$1"
  if declare -F gwas_clean_zcat >/dev/null 2>&1; then
    gwas_clean_zcat "$f"
    return
  fi
  case "$f" in
    *.gz|*.bgz)
      if command -v pigz >/dev/null 2>&1; then
        pigz -dc -p "${GWAS_CLEAN_DECOMP_THREADS:-1}" "$f"
      else
        gzip -dc "$f"
      fi
      ;;
    *)
      cat "$f"
      ;;
  esac
}

n_data_rows() { [[ -s $1 ]] && phe_zcat "$1" 2>/dev/null | awk 'NR>1{n++} END{print n+0}' || echo 0; }

gzip_ok() { [[ -s $1 ]] && gzip -t "$1" >/dev/null 2>&1; }

join_by() { local IFS="$1"; shift; echo "$*"; }

join_comma() { local IFS=,; echo "$*"; }

limit_jobs() { local max=$1; while (( $(jobs -pr | wc -l) >= max )); do wait -n || true; done; }

wait_all() { local rc=0; while (( $(jobs -pr | wc -l) > 0 )); do wait -n || rc=1; done; return $rc; }

header_names() {
    Arr1=("SNP" "CHR" "POS" "EA" "NEA" "EAF" "N" "BETA" "SE" "P" "LOG10P")
    Arr2=("^variant_id$|^snp$|^id$|^rsid$|^variant_ids$" "^chr$|^chrom$|^chromosome$|^chr_name$" "^pos$|^bp$|^base_pair$|^base_pair_location$|^genpos$|^hm_pos$" "^a1$|^ea$|^eff.allele$|^effect_allele$|^allele1$" "^OMITTED$|^nea$|^non_effect_allele$|^other_allele$|^allele0$|^allele2$|^ref.allele$|^reference_allele$|^ref$" "^eaf$|^a1freq$|^a1_freq$|^effect_allele_frequency" "^n$|^obs_ct$|^Neff$" "^beta$|^effect_weight$" "^se$|^standard_error" "^p$|^pval$|^p_value$|^p_bolt_lmm$" "^log10p$|^neg.log.10.p.value$|^neg.log10.p.value$|^negative.log.10.p.value$|^negative.log10.p.value$|^minus.log10.p$|^mlog10p$")    
    datf=$1
    head_row=$(phe_zcat "$datf" | head -1 | sed 's/\t/ /g')
    head_row=$(echo "$head_row" | dos2unix)
    head_list=$(echo "$head_row" | tr ' ' '\n')
    for i in ${!Arr1[@]}; do
        matches=""
        IFS='|' read -r -a prefer_patterns <<< "${Arr2[$i]}"
        for pat in "${prefer_patterns[@]}"; do
            matches=$(echo "$head_list" | grep -Einw "$pat" | sed 's/:.*//')
            match_count=$(echo "$matches" | grep -c -v "^$")
            [[ $match_count -gt 0 ]] && break
        done
        match_count=$(echo "$matches" | grep -c -v "^$")
        if [[ $match_count -gt 0 ]]; then
            idx=$(echo "$matches" | head -n 1)
            eval "${Arr1[$i]}=$(echo "$head_row" | cut -d' ' -f $idx)"
            eval "${Arr1[$i]}_col=$idx"
        else
            eval "${Arr1[$i]}="
            eval "${Arr1[$i]}_col="
        fi
    done
    echo "dat $datf, snp $SNP $SNP_col, chr $CHR $CHR_col, pos $POS $POS_col, ea $EA $EA_col, nea $NEA $NEA_col, eaf $EAF $EAF_col, n $N $N_col, beta $BETA $BETA_col, se $SE $SE_col, p $P $P_col, log10p $LOG10P $LOG10P_col"
}

#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/common.sh" "$@"; grid_parse_args "$@"
trait=${GRID_TRAIT,,}; POPS=($(grid_csv_words "$GRID_POPS")); CHRS=($(grid_expand_chrs "$GRID_CHRS"))
out="$GRID_OUTPUT_ROOT/$trait"; prep="$out/sumstats"; raw="$out/csx/raw"; weights="$out/csx/weights"; score="$out/scores"
mkdir -p "$prep" "$raw" "$weights" "$score" "$out/log"
prscx=$(find "$ROOT/f/csx" -maxdepth 4 -type f \( -iname 'PRScsx.py' -o -iname 'prscsx.py' \) -print -quit 2>/dev/null || true)
[[ -s $prscx ]] || _grid_die "Bundled PRS-CSx program not found below $ROOT/f/csx"
[[ -s $GRID_CSX_SNPINFO ]] || _grid_die "Missing PRS-CSx SNP info: $GRID_CSX_SNPINFO"
[[ -s $GRID_CSX_BIM_PREFIX.bim ]] || _grid_die "Missing target validation BIM: $GRID_CSX_BIM_PREFIX.bim"

n_override(){
  python3 - "$1" "$GRID_N_GWAS" <<'PY'
import sys,re
p=sys.argv[1].upper(); s=sys.argv[2].strip()
if not s: raise SystemExit(1)
if re.fullmatch(r'[0-9.eE+-]+',s): print(s); raise SystemExit
for x in re.split(r'[,; ]+',s):
 if '=' in x:
  k,v=x.split('=',1)
  if k.upper()==p: print(v); raise SystemExit
raise SystemExit(1)
PY
}

sst=(); ng=(); poplow=()
for p in "${POPS[@]}"; do
  p=${p^^}; pl=${p,,}; g=$(grid_find_gwas "$trait" "$p" || true); [[ -s $g ]] || _grid_die "Missing $GRID_GWAS_DIR/$trait.$p.gz"
  std="$prep/$trait.$p.hm3.tsv.gz"; meta="$prep/$trait.$p.meta.json"
  if [[ ! -s $std || $GRID_REPLACE == TRUE ]]; then
    grid_run_logged "$out/log/csx.prepare.$p.log" python3 "$ROOT/f/prepare_sumstats.py" --input "$g" --output "$std" --metadata "$meta" --snpinfo "$GRID_CSX_SNPINFO" --trait "$trait" --pop "$p" --chunk "$GRID_SUMSTATS_CHUNK"
  fi
  grid_run python3 "$ROOT/f/split_sumstats.py" --input "$std" --out-dir "$prep/bychr" --prefix "$trait.$p" --chrs "${CHRS[*]}"
  if n=$(n_override "$p" 2>/dev/null); then :
  else n=$(python3 - "$meta" <<'PY'
import json,sys
x=json.load(open(sys.argv[1])).get('n_gwas_median')
if x is None: raise SystemExit(1)
print(int(round(float(x))))
PY
    ) || _grid_die "GWAS sample size unavailable for $trait.$p; provide --n-gwas '$p=N'"
  fi
  sst+=("$prep/bychr/$trait.$p.chrCHR.tsv.gz"); ng+=("$n"); poplow+=("$p")
done

join_comma(){ local IFS=,; echo "$*"; }
for c in "${CHRS[@]}"; do
  files=(); for x in "${sst[@]}"; do files+=("${x/CHR/$c}"); done
  outname="$trait.chr$c"
  # Run all discovery populations jointly so the continuous shrinkage prior is shared.
  marker="$raw/$outname.done"
  if [[ ! -s $marker || $GRID_REPLACE == TRUE ]]; then
    cmd=(python3 "$prscx" --ref_dir="$GRID_CSX_REF_DIR" --bim_prefix="$GRID_CSX_BIM_PREFIX" --sst_file="$(join_comma "${files[@]}")" --n_gwas="$(join_comma "${ng[@]}")" --pop="$(join_comma "${poplow[@]}")" --chrom="$c" --phi="$GRID_PHI" --n_iter="$GRID_MCMC_ITER" --n_burnin="$GRID_MCMC_BURNIN" --thin="$GRID_MCMC_THIN" --out_dir="$raw" --out_name="$outname")
    grid_run_logged "$out/log/csx.chr$c.log" "${cmd[@]}"
    date -Is > "$marker"
  fi
  for p in "${POPS[@]}"; do
    p=${p^^}; pl=${p,,}
    # Prefer the repository's normalized output layout, then standard upstream names.
    cand=$(find "$raw" -type f \( -path "*/$p/$p.chr$c.pst_eff.txt" -o -iname "*chr${c}*${p}*pst_eff*.txt" -o -iname "*chr${c}*${pl}*pst_eff*.txt" \) -print | sort | head -n1)
    [[ -s $cand ]] || _grid_die "Cannot locate PRS-CSx posterior weights for $trait $p chr$c under $raw"
    w="$weights/$p.chr$c.tsv"
    [[ -s $w && $GRID_REPLACE == FALSE ]] || grid_run python3 "$ROOT/f/normalize_csx_weights.py" --input "$cand" --output "$w"
  done
done

score_one(){
  local p=$1 c=$2 mode prefix w o
  mode=$(grid_target_mode "$c" || true); [[ -n $mode ]] || _grid_die "Missing target chr$c bfile/pfile under $GRID_TARGET_DIR"
  prefix="$GRID_TARGET_DIR/chr$c"; w="$weights/$p.chr$c.tsv"; o="$score/tmp/CSX_${p}.chr$c"; mkdir -p "$score/tmp"
  [[ -s $o.sscore && $GRID_REPLACE == FALSE ]] && return
  cmd=(plink2 "--$mode" "$prefix" --score "$w" 1 2 3 header-read no-mean-imputation cols=+scoresums --threads "$GRID_THREADS" --out "$o")
  [[ -z $GRID_KEEP ]] || cmd+=(--keep "$GRID_KEEP")
  [[ -z $GRID_REMOVE ]] || cmd+=(--remove "$GRID_REMOVE")
  grid_run_logged "$out/log/csx.score.$p.chr$c.log" "${cmd[@]}"
}
for p in "${POPS[@]}"; do
  p=${p^^}; running=0; status=0
  for c in "${CHRS[@]}"; do score_one "$p" "$c" & ((++running)); if ((running>=GRID_JOBS)); then wait -n || status=1; running=$((running-1)); fi; done
  while ((running)); do wait -n || status=1; running=$((running-1)); done
  ((status==0)) || _grid_die "CSX scoring failed for $p"
  inputs=(); for c in "${CHRS[@]}"; do inputs+=("$score/tmp/CSX_${p}.chr$c.sscore"); done
  grid_run python3 "$ROOT/f/combine_scores.py" --inputs "${inputs[@]}" --name "CSX_$p" --output "$score/CSX_$p.tsv.gz"
done
merge=(); for p in "${POPS[@]}"; do merge+=("$score/CSX_${p^^}.tsv.gz"); done
grid_run python3 "$ROOT/f/merge_scores.py" --inputs "${merge[@]}" --output "$score/csx.tsv.gz"
{
  echo -e 'trait\tpop\tgwas\tn_gwas\tweights_dir\tscore'
  i=0; for p in "${POPS[@]}"; do g=$(grid_find_gwas "$trait" "$p"); echo -e "$trait\t${p^^}\t$g\t${ng[$i]}\t$weights\t$score/CSX_${p^^}.tsv.gz"; ((++i)); done
} > "$out/csx/manifest.tsv"
echo "PRS-CSx completed: $score/csx.tsv.gz"

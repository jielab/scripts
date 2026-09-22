#!/usr/bin/env bash
# Joint four-population inference -> permanent weights -> UKB scores.
if [[ -z ${GRID_COMMAND_FILE:-} ]]; then
  exec bash "$(dirname -- "${BASH_SOURCE[0]}")/pipeline.sh" csx "$@"
fi
set -euo pipefail
trait=$GRID_TRAIT
score_home="$GRID_SCORE_DIR/$trait${suffix:+/${suffix#.}}"
io="$ROOT/f/pipeline_io.py"
prscx="$ROOT/f/csx/PRScsx.py"
gwas=(); finals=()
for p in "${POPS[@]}"; do
  g=$(grid_find_gwas "$trait" "$p") || _grid_die "Missing GWAS for $trait.$p below $GRID_GWAS_DIR"
  finals+=("${g%.gz}$suffix.csx.gz")
  gwas+=("$g")
done
python3 - "$GRID_MCMC_ITER" "$GRID_MCMC_BURNIN" "$GRID_MCMC_THIN" "$GRID_SEED" "$GRID_PHI" <<'PY'
import sys,math
n,b,t,s=map(int,sys.argv[1:5]); phi=sys.argv[5]
assert n>b>=0 and t>0 and n//t-b//t>=2 and s>=0, 'Need valid MCMC settings and at least two retained samples'
assert phi=='auto' or (math.isfinite(float(phi)) and float(phi)>0), 'phi must be auto or positive'
PY
need "$GRID_CSX_SNPINFO"; need "$GRID_CSX_BIM_PREFIX.bim"
# PRS-CSx chooses reference type by SNPINFO filename; explicitly stage only the selected one.
case $(basename -- "$GRID_CSX_SNPINFO") in
  snpinfo_mult_1kg_hm3) ref_type=1kg;;
  snpinfo_mult_ukbb_hm3) ref_type=ukbb;;
  *) _grid_die 'Use snpinfo_mult_1kg_hm3 or snpinfo_mult_ukbb_hm3';;
esac
ref_dirs=(); ld_files=()
for p in "${POPS[@]}"; do
  ref="$GRID_CSX_REF_DIR/ldblk_${ref_type}_${p,,}"
  [[ -d $ref ]] || ref="$GRID_CSX_REF_DIR/ldblk_${ref_type}_$p"
  for c in "${CHRS[@]}"; do need "$ref/ldblk_${ref_type}_chr$c.hdf5"; ld_files+=("$ref/ldblk_${ref_type}_chr$c.hdf5"); done
  ref_dirs+=("$ref")
done
if [[ $GRID_STAGE != weights ]]; then
  command -v plink2 >/dev/null || _grid_die 'plink2 is missing'
  for c in "${CHRS[@]}"; do grid_target_mode "$c" >/dev/null || _grid_die "Missing target chr$c under $GRID_TARGET_DIR"; done
  [[ -z $GRID_KEEP ]] || need "$GRID_KEEP"
  [[ -z $GRID_REMOVE || -f $GRID_REMOVE ]] || _grid_die "Missing withdrawal file: $GRID_REMOVE"
fi
python3 "$io" inspect "$GRID_CSX_SNPINFO" "${gwas[@]}"
python3 "$io" coverage "${CHRS[*]}" "${gwas[@]}"
for f in "${finals[@]}"; do echo "output SNP weights: $f"; done
echo "output score files: $score_home/csx.pgs.gz"
echo "CSx: joint AFR,EAS,EUR,SAS; phi=$GRID_PHI; iter/burnin/thin=$GRID_MCMC_ITER/$GRID_MCMC_BURNIN/$GRID_MCMC_THIN; seed=$GRID_SEED; chromosomes=${CHRS[*]}"
if [[ $GRID_STAGE == score ]]; then
  for f in "${finals[@]}"; do need "$f"; need "$f.signature"; done
fi
[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || { echo 'CHECK/PLAN complete; no inference or scoring executed'; return 0; }
# Include source code and input metadata so interrupted runs cannot reuse incompatible results.
sig=$(python3 "$ROOT/f/csx_config.py" "$GRID_PHI" "$GRID_MCMC_ITER" "$GRID_MCMC_BURNIN" "$GRID_MCMC_THIN" "$GRID_SEED" "$GRID_N_GWAS" "${CHRS[*]}" "$GRID_CSX_SNPINFO" "$GRID_CSX_BIM_PREFIX" "${gwas[@]}" "${ld_files[@]}")
run="$work/$trait/$sig"; mkdir -p "$run" "$work/log/$trait"
exec {lock}>"$work/$trait/run.lock"
flock -n "$lock" || _grid_die "Another CSx run is active for $trait"
ref="$run/reference"; mkdir -p "$ref"
ln -sfn "$GRID_CSX_SNPINFO" "$ref/$(basename -- "$GRID_CSX_SNPINFO")"
for i in "${!POPS[@]}"; do ln -sfn "${ref_dirs[$i]}" "$ref/ldblk_${ref_type}_${POPS[$i],,}"; done
if [[ $GRID_STAGE != score ]]; then
  complete=TRUE
  for f in "${finals[@]}"; do [[ -s $f && -s $f.signature && $(cat "$f.signature") == "$sig" ]] || complete=FALSE; done
  if [[ $complete == TRUE && $GRID_REPLACE == FALSE ]]; then
    echo "SKIP $trait inference: matching permanent weights"
  else
    mkdir -p "$run/sumstats" "$run/raw" "$run/weights"
    ng=()
    for i in "${!POPS[@]}"; do
      p=${POPS[$i]}; std="$run/sumstats/$p.tsv.gz"; meta="$run/sumstats/$p.json"
      grid_run_logged "$work/log/$trait/prepare.$p.log" python3 "$ROOT/f/sumstats_cache.py" --input "${gwas[$i]}" --output "$std" --metadata "$meta" --snpinfo "$GRID_CSX_SNPINFO" --trait "$trait" --pop "$p" --chunk "$GRID_SUMSTATS_CHUNK" --work "$work" --replace "$GRID_REPLACE"
      tail -n 1 "$work/log/$trait/prepare.$p.log"
      grid_run_logged "$work/log/$trait/split.$p.log" python3 "$ROOT/f/split_sumstats.py" --input "$std" --out-dir "$run/sumstats" --prefix "$p" --chrs "${CHRS[*]}"
      n=$(python3 - "$meta" "$GRID_N_GWAS" "$p" <<'PY'
import sys,json,re,math
meta,override,pop=sys.argv[1:]; value=None
if override:
 if '=' not in override: value=float(override)
 else:
  opts=dict(x.split('=',1) for x in re.split('[,; ]+',override) if x)
  value=float(opts[pop]) if pop in opts else None
if value is None: value=json.load(open(meta))['n_gwas_median']
if value is None or not math.isfinite(value) or value<=0: raise SystemExit('Missing/invalid GWAS N; provide --n-gwas')
print(round(value))
PY
      )
      python3 - "$meta" "$n" "$GRID_PHI" <<'PYMETA'
import json,sys
p,n,phi=sys.argv[1:];d=json.load(open(p));d.update(n_gwas_used=int(n),inference_phi=phi);open(p,'w').write(json.dumps(d,indent=2)+'\n')
PYMETA
      ng+=("$n"); echo "  $trait.$p: n_gwas=$n"
    done
    infer_chr(){
      local c=$1 raw="$run/raw/chr$1" marker="$run/raw/chr$1/done" p
      mkdir -p "$raw"
      if [[ -f $marker && $GRID_REPLACE == FALSE ]]; then
        local valid=TRUE
        for p in "${POPS[@]}"; do [[ -s $raw/$p/$p.chr$c.pst_eff.txt ]] || valid=FALSE; done
        compgen -G "$raw/*META*pst_eff*.txt" >/dev/null || valid=FALSE
        [[ $valid != TRUE ]] || { echo "SKIP $trait chr$c: completed posterior"; return; }
      fi
      rm -f -- "$marker"
      local files=() cmd=()
      for p in "${POPS[@]}"; do files+=("$run/sumstats/$p.chr$c.tsv"); done
      cmd=(env OMP_NUM_THREADS="$GRID_THREADS" OPENBLAS_NUM_THREADS="$GRID_THREADS" MKL_NUM_THREADS="$GRID_THREADS" python3 "$prscx" --ref_dir="$ref" --bim_prefix="$GRID_CSX_BIM_PREFIX" --sst_file="$(join_comma "${files[@]}")" --n_gwas="$(join_comma "${ng[@]}")" --pop="$(join_comma "${POPS[@]}")" --chrom="$c" --n_iter="$GRID_MCMC_ITER" --n_burnin="$GRID_MCMC_BURNIN" --thin="$GRID_MCMC_THIN" --seed="$((GRID_SEED+c))" --out_dir="$raw" --out_name="$trait" --meta=TRUE)
      [[ $GRID_PHI == auto ]] || cmd+=(--phi="$GRID_PHI")
      echo "RUN $trait chr$c: joint PRS-CSx"
      grid_run_logged "$work/log/$trait/csx.chr$c.log" "${cmd[@]}" || return $?
      for p in "${POPS[@]}"; do need "$raw/$p/$p.chr$c.pst_eff.txt"; done
      touch "$marker"
    }
    # Wait for each bounded batch explicitly: no lost child failures with wait -n.
    pids=()
    for c in "${CHRS[@]}"; do
      infer_chr "$c" & pids+=("$!")
      if ((${#pids[@]} >= GRID_JOBS)); then status=0; for pid in "${pids[@]}"; do wait "$pid" || status=1; done; ((status==0)) || _grid_die 'CSx inference failed'; pids=(); fi
    done
    status=0; for pid in "${pids[@]}"; do wait "$pid" || status=1; done; ((status==0)) || _grid_die 'CSx inference failed'
    for i in "${!POPS[@]}"; do
      p=${POPS[$i]}; inputs=()
      for c in "${CHRS[@]}"; do
        w="$run/weights/$p.chr$c.tsv"
        grid_run python3 "$ROOT/f/normalize_csx_weights.py" --input "$run/raw/chr$c/$p/$p.chr$c.pst_eff.txt" --output "$w"
        inputs+=("$w")
      done
      grid_run python3 "$io" weights "$run/weights/$p.csx.gz" "${inputs[@]}"
      publish "$run/weights/$p.csx.gz" "${finals[$i]}"
      printf '%s\n' "$sig" > "$run/signature"; publish "$run/signature" "${finals[$i]}.signature"
      publish "$run/sumstats/$p.json" "${finals[$i]}.metadata.json"
    done
  fi
fi
if [[ $GRID_STAGE == weights ]]; then echo "DONE $trait: permanent SNP weights saved"; return 0; fi
source "$ROOT/f/csx_score.sh"

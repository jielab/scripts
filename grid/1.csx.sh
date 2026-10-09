#!/usr/bin/env bash
# PRS-CSx: GWAS preparation -> joint MCMC -> permanent SNP weights -> UKB PRS.
# Inference, population/combined scores and resource caps are orchestrated here.
# Python helpers: f/1.csx.py; official PRS-CSx code: f/csx/.
# Reference: https://github.com/getian107/PRScsx
set -euo pipefail
usage() {
	cat <<'HELP'
PRS-CSx — height / ldl / t2dm, jointly using AFR,EAS,EUR,SAS GWAS

Usage examples (WSL):
  cd /mnt/d/scripts/grid
  ./1.csx.sh --traits height,ldl,t2dm --check
  ./1.csx.sh --traits height,ldl,t2dm --jobs 4 --threads 4

  # The full run above enables --posterior TRUE by default.
  # Optional: append only auto/meta scores (does not generate individual posterior scores).
  ./1.csx.sh --traits height,ldl,t2dm --models auto,meta --jobs 4 --threads 4
  # Optional: run inference and scoring separately; retain the inference directory.
  ./1.csx.sh --trait height --stage weights --jobs 4 --threads 4
  ./1.csx.sh --trait height --stage score --jobs 4 --threads 4

Modules:
  weights  Normalize BETA+SE GWAS -> joint four-population PRS-CSx by chromosome.
  score    Apply saved SNP weights to UKB; four individual PRS files.
  all      Both modules (default). No PCA is needed for CSx.

Important parameters:
  --trait NAME / --traits LIST   Default: height,ldl,t2dm; traits run sequentially.
  --stage all|weights|score      Default: all.
  --models LIST                Default: populations,auto,meta.
                               Use auto,meta to append only the two combined scores.
                               auto: learned phi + posterior meta; meta: fixed --phi + posterior meta.
  --dir-gwas DIR                /mnt/f/gwas/4grid/common
  --dir-gen DIR                 /mnt/f/gen/ukb/37/hap (chr*.pgen/pvar/psam)
  --csx-bim-prefix PREFIX       Target variant list; default: above DIR/ukb_array
  --csx-ref-dir DIR             /mnt/f/refLD/csx (1000 Genomes LD)
  --csx-snpinfo FILE            Default: snpinfo_mult_1kg_hm3 in default LD folder
  --phi VALUE|auto              Default: 1e-2; fixed value, NOT phenotype-tuned.
  --posterior TRUE|FALSE       TRUE; save synchronized population draws and score them.
  --posterior-frequency-dir DIR Optional AFR/EAS/EUR/SAS.tsv.gz discovery EAF tables.
                               Defaults to normalized GWAS EAF; MAF is not accepted.
  --memory-cap-gb N             16 GiB RAM for the entire invocation, including all parallel jobs.
  --swap-cap-gb N               2 GiB swap for the entire invocation; 0 disables task swap.
                               Requires systemd user manager + cgroup v2; fails closed if unavailable.
                               Env defaults: GRID_MEMORY_CAP_GB / GRID_SWAP_CAP_GB.
  --score-memory MB             2048 MiB workspace per ordinary PLINK scoring call.
  --posterior-memory MB         8192 per PLINK posterior scoring call; chromosomes scored sequentially.
  --mcmc-iter N                 4000; --mcmc-burnin 2000; --mcmc-thin 5
  --n-gwas N|AFR=N,EAS=N,...     Default: median N over retained HM3 SNPs.
                                For t2dm, verify N is effective sample size.
  --jobs N / --threads N        Parallel chromosomes / threads per job (4 / 8).
  --chrs 1-22                  Autosomes; subsets get separate filenames/directories.
  --seed N                     20260904 (chromosome number added for each chain).
  --keep FILE / --remove FILE  PLINK sample lists; default remove: /mnt/d/files/ukb.exclude.id.
  --score-dir DIR              /mnt/d/data/ukb/pgs
  --output-root DIR            /mnt/d/analysis/grid (temporary).
  --replace TRUE|FALSE         Default FALSE; reuse only matching saved results.
  --check                      Validate inputs only; no inference/scoring.
  --dry-run TRUE               Same preflight/plan behavior; no completion markers.

Permanent outputs:
  <GWAS folder>/<original-name>.csx.sumstats.gz (normalized GWAS BETA,SE,N,...)
  Matching .csx.sumstats.json enables reuse after deleting the temporary work directory.
  <GWAS folder>/<original-name>.csx.gz  (SNP,A1,BETA,CHR,BP,A2; NOT individual PRS)
  Example: /mnt/f/gwas/4grid/common/height.AFR/gwas/height.AFR.csx.gz
  /mnt/d/data/ukb/pgs/<trait>/1csx.scores.rds (eid, csx.AFR/EAS/EUR/SAS, csx.auto, csx.meta)
  1csx.scores.provenance.json connects every published column to its source weights.
  Rerun --stage score to attach this provenance to older scores without rerunning MCMC.
  1csx.posterior.rds: four centred posterior means + ten covariance terms per person.
  RDS attributes retain discovery EAF centering and the fitted-input identity.
  Keep joint_posterior.h5 files under temporary inference directories to rescore without MCMC.
  Inference directories: csx/<trait>/phi-1e-2/ or auto/; config.json records settings.
  Different inputs/settings use numbered suffixes; no hash-based directory names.
  Combined SNP weights: <score-dir>/<trait>/.weights/csx.{auto,meta}.gz
  Per-population score exchanges and signatures: /tmp/grid-cache/.
  t2dm.AFA is mapped to AFR LD, but its GWAS output retains the name t2dm.AFA.

Temporary files/commands/logs: /mnt/d/analysis/grid/csx/
Input coordinates must be GRCh37; the preflight checks rsID-position agreement.
For another reference/target, set --csx-snpinfo and --csx-bim-prefix explicitly.
HELP
}
case "${1:-}" in -h | --help | help)
	usage
	exit 0
	;;
esac
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/0.common.sh"
memory_cap_enter csx "${GRID_MEMORY_CAP_GB:-16}" "${GRID_SWAP_CAP_GB:-2}" "$ROOT/1.csx.sh" "$@"
set -- "${MEMORY_CAP_ARGS[@]}"
export GRID_MEMORY_CAP_GB=$MEMORY_CAP_GB


# 🚩 csx_score_run
csx_score_run() {
	# Score every UKB individual with each population's permanent SNP weights.
	for f in "${finals[@]}"; do
		need "$f"
		need "$f.signature"
		[[ $(cat "$f.signature") == "$sig" ]] || _grid_die "Weights/settings mismatch: $f; rerun --stage weights with these settings"
	done
	score_inputs=("${finals[@]}")
	for c in "${CHRS[@]}"; do
		mode=$(grid_target_mode "$c")
		prefix="$GRID_TARGET_DIR/chr$c"
		if [[ $mode == pfile ]]; then
			score_inputs+=("$prefix.pgen" "$prefix.psam")
			if [[ -s $prefix.pvar ]]; then score_inputs+=("$prefix.pvar"); else score_inputs+=("$prefix.pvar.zst"); fi
		else score_inputs+=("$prefix.bed" "$prefix.bim" "$prefix.fam"); fi
	done
	[[ -z $GRID_KEEP ]] || score_inputs+=("$GRID_KEEP")
	[[ -z $GRID_REMOVE ]] || score_inputs+=("$GRID_REMOVE")
	# Bind newly attested scores to weight CONTENT, including an initial rescore
	# of older caches whose keys contained only path/size/mtime information.
	weight_hashes=$(sha256sum -- "${finals[@]}")
	score_sig=$(python3 "$ROOT/f/1.csx.py" score-config "$sig" "$GRID_KEEP" "$GRID_REMOVE" "$weight_hashes" --files "${score_inputs[@]}")
	score_run=$(python3 "$ROOT/f/1.csx.py" workspace "$run/scores" run "$score_sig")
	mkdir -p "$score_home"
	score_cache="$(python3 "$ROOT/f/0.common.py" cache-path "$work/scores")/$trait${suffix:+/${suffix#.}}"
	mkdir -p "$score_cache"
	for i in "${!POPS[@]}"; do
		p=${POPS[$i]}
		dest="$score_cache/csx.$p.tsv.gz"
		if [[ -s $dest && -s $dest.signature && $(cat "$dest.signature") == "$score_sig" && $GRID_REPLACE == FALSE ]]; then
			echo "SKIP $trait $p: matching cached scores"
			continue
		fi
		score_chr() {
			local c=$1 mode prefix o
			mode=$(grid_target_mode "$c")
			prefix="$GRID_TARGET_DIR/chr$c"
			o="$score_run/$p.chr$c"
			if [[ -s $o.sscore && -f $o.done && $GRID_REPLACE == FALSE ]]; then return; fi
			rm -f -- "$o.done"
			local cmd=(plink2 "--$mode" "$prefix")
			if [[ $mode == pfile && ! -s $prefix.pvar ]]; then cmd+=(vzs); fi
			cmd+=(--score "${finals[$i]}" 1 2 3 header-read no-mean-imputation list-variants cols=+scoresums --threads "$GRID_THREADS" --memory "$GRID_SCORE_MEMORY" --out "$o")
			[[ -z $GRID_KEEP ]] || cmd+=(--keep "$GRID_KEEP")
			[[ ! -s $GRID_REMOVE ]] || cmd+=(--remove "$GRID_REMOVE")
			grid_run_logged "$logdir/$trait/score.$p.chr$c.log" "${cmd[@]}" || return $?
			need "$o.sscore"
			touch "$o.done"
		}
		echo "RUN $trait $p: UKB scoring"
		pids=()
		for c in "${CHRS[@]}"; do
			score_chr "$c" &
			pids+=("$!")
			if ((${#pids[@]} >= GRID_SCORE_JOBS)); then
				status=0
				for pid in "${pids[@]}"; do wait "$pid" || status=1; done
				((status == 0)) || _grid_die 'CSx scoring failed'
				pids=()
			fi
		done
		status=0
		for pid in "${pids[@]}"; do wait "$pid" || status=1; done
		((status == 0)) || _grid_die 'CSx scoring failed'
		inputs=()
		for c in "${CHRS[@]}"; do inputs+=("$score_run/$p.chr$c.sscore"); done
		grid_run_logged "$logdir/$trait/combine.$p.log" python3 "$ROOT/f/0.common.py" combine-scores --inputs "${inputs[@]}" --name "CSX_$p" --output "$score_run/$p.tsv.gz"
		publish "$score_run/$p.tsv.gz" "$dest"
		printf '%s\n' "$score_sig" >"$score_run/signature"
		publish "$score_run/signature" "$dest.signature"
	done
	merge=()
	for p in "${POPS[@]}"; do merge+=("$score_cache/csx.$p.tsv.gz"); done
	grid_run python3 "$ROOT/f/0.common.py" merge-scores --inputs "${merge[@]}" --output "$score_run/csx.tsv.gz"
	provenance_weights=()
	for i in "${!POPS[@]}"; do provenance_weights+=("csx.${POPS[$i]}=${finals[$i]}"); done
	grid_run python3 "$ROOT/f/0.common.py" score-provenance --weights "${provenance_weights[@]}" --signature "$sig" --chrs "${CHRS[*]}" --output "$score_run/provenance.json"
	grid_run python3 "$ROOT/f/0.common.py" publish csx "$score_run/csx.tsv.gz" "$score_home/1csx.scores.rds" --remove "$GRID_REMOVE" --provenance "$score_run/provenance.json"
	{
		printf 'trait\tpopulation\tinput_gwas\tweights\tukb_score\n'
		for i in "${!POPS[@]}"; do printf '%s\t%s\t%s\t%s\t%s\n' "$trait" "${POPS[$i]}" "${gwas[$i]}" "${finals[$i]}" "$score_home/1csx.scores.rds"; done
	} >"$score_run/manifest.tsv"
	echo "DONE $trait: $score_home/1csx.scores.rds"
}


# 🚩 csx_combined_run
csx_combined_run() {
	set -euo pipefail
	args=(--trait "$GRID_TRAIT" --mode "$csx_model" --gwas-dir "$GRID_GWAS_DIR"
		--target-dir "$GRID_TARGET_DIR" --snpinfo "$GRID_CSX_SNPINFO" --ref-dir "$GRID_CSX_REF_DIR"
		--bim "$GRID_CSX_BIM_PREFIX" --work "$work" --score-home "$GRID_SCORE_DIR/$GRID_TRAIT${suffix:+/${suffix#.}}"
		--chrs "${CHRS[*]}" --jobs "$GRID_JOBS" --threads "$GRID_THREADS" --iterations "$GRID_MCMC_ITER"
		--score-jobs "$GRID_SCORE_JOBS" --score-memory "$GRID_SCORE_MEMORY"
		--burnin "$GRID_MCMC_BURNIN" --thin "$GRID_MCMC_THIN" --seed "$GRID_SEED" --phi "$GRID_PHI"
		--n-gwas "$GRID_N_GWAS" --remove "$GRID_REMOVE" --keep "$GRID_KEEP" --stage "$GRID_STAGE" --replace "$GRID_REPLACE")
	[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || args+=(--check)
	printf 'python3 %q combined ' "$ROOT/f/1.csx.py" >>"$GRID_COMMAND_FILE"
	printf '%q ' "${args[@]}" >>"$GRID_COMMAND_FILE"
	printf '\n' >>"$GRID_COMMAND_FILE"
	python3 "$ROOT/f/1.csx.py" combined "${args[@]}"
}


# 🚩 csx_population_run
csx_population_run() {
	# Joint four-population inference -> permanent weights -> UKB scores.
	set -euo pipefail
	trait=$GRID_TRAIT
	score_home="$GRID_SCORE_DIR/$trait${suffix:+/${suffix#.}}"
	io="$ROOT/f/0.common.py"
	prscx="$ROOT/f/csx/PRScsx.py"
	gwas=()
	finals=()
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
	if [[ $GRID_POSTERIOR == TRUE ]]; then
		((GRID_MCMC_ITER / GRID_MCMC_THIN - GRID_MCMC_BURNIN / GRID_MCMC_THIN >= 100)) || _grid_die 'Posterior scoring requires at least 100 retained draws; increase iterations or reduce thinning'
	fi
	need "$GRID_CSX_SNPINFO"
	need "$GRID_CSX_BIM_PREFIX.bim"
	# PRS-CSx chooses reference type by SNPINFO filename; explicitly stage only the selected one.
	case $(basename -- "$GRID_CSX_SNPINFO") in
		snpinfo_mult_1kg_hm3) ref_type=1kg ;;
		snpinfo_mult_ukbb_hm3) ref_type=ukbb ;;
		*) _grid_die 'Use snpinfo_mult_1kg_hm3 or snpinfo_mult_ukbb_hm3' ;;
	esac
	ref_dirs=()
	ld_files=()
	for p in "${POPS[@]}"; do
		ref="$GRID_CSX_REF_DIR/ldblk_${ref_type}_${p,,}"
		[[ -d $ref ]] || ref="$GRID_CSX_REF_DIR/ldblk_${ref_type}_$p"
		for c in "${CHRS[@]}"; do
			need "$ref/ldblk_${ref_type}_chr$c.hdf5"
			ld_files+=("$ref/ldblk_${ref_type}_chr$c.hdf5")
		done
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
	echo "output score files: $score_home/1csx.scores.rds"
	echo "CSx: joint AFR,EAS,EUR,SAS; phi=$GRID_PHI; iter/burnin/thin=$GRID_MCMC_ITER/$GRID_MCMC_BURNIN/$GRID_MCMC_THIN; seed=$GRID_SEED; chromosomes=${CHRS[*]}"
	if [[ $GRID_STAGE == score ]]; then
		for f in "${finals[@]}"; do
			need "$f"
			need "$f.signature"
		done
	fi
	[[ $GRID_CHECK == FALSE && $GRID_DRY_RUN == FALSE ]] || {
		echo 'CHECK/PLAN complete; no inference or scoring executed'
		return 0
	}
	# Keep settings/input metadata in JSON; choose readable directories independently.
	sig=$(python3 "$ROOT/f/1.csx.py" signature "$GRID_PHI" "$GRID_MCMC_ITER" "$GRID_MCMC_BURNIN" "$GRID_MCMC_THIN" "$GRID_SEED" "$GRID_N_GWAS" "${CHRS[*]}" "$GRID_CSX_SNPINFO" "$GRID_CSX_BIM_PREFIX" "${gwas[@]}" "${ld_files[@]}")
	mkdir -p "$work/$trait" "$logdir/$trait"
	local lock_dir
	lock_dir=$(python3 "$io" cache-path "$work/$trait")
	mkdir -p "$lock_dir"
	exec {lock}>"$lock_dir/run.lock"
	flock -n "$lock" || _grid_die "Another CSx run is active for $trait"
	run=$(python3 "$ROOT/f/1.csx.py" workspace "$work/$trait" inference "$sig")
	sumstats=$(python3 "$io" cache-path "$run/sumstats")
	ref=$(python3 "$io" cache-path "$run/reference")
	mkdir -p "$ref"
	ln -sfn "$GRID_CSX_SNPINFO" "$ref/$(basename -- "$GRID_CSX_SNPINFO")"
	for i in "${!POPS[@]}"; do ln -sfn "${ref_dirs[$i]}" "$ref/ldblk_${ref_type}_${POPS[$i],,}"; done
	if [[ $GRID_STAGE != score ]]; then
		complete=TRUE
		for f in "${finals[@]}"; do [[ -s $f && -s $f.signature && $(cat "$f.signature") == "$sig" ]] || complete=FALSE; done
		if [[ $GRID_POSTERIOR == TRUE ]]; then
			for c in "${CHRS[@]}"; do [[ -s $run/raw/chr$c/joint_posterior.h5 ]] || complete=FALSE; done
		fi
		if [[ $complete == TRUE && $GRID_REPLACE == FALSE ]]; then
			echo "SKIP $trait inference: matching permanent weights"
		else
			mkdir -p "$sumstats" "$run/raw" "$run/weights"
			ng=()
			for i in "${!POPS[@]}"; do
				p=${POPS[$i]}
				std="$sumstats/$p.tsv.gz"
				meta="$sumstats/$p.json"
				grid_run_logged "$logdir/$trait/prepare.$p.log" python3 "$ROOT/f/0.common.py" sumstats-cache --input "${gwas[$i]}" --output "$std" --metadata "$meta" --snpinfo "$GRID_CSX_SNPINFO" --trait "$trait" --pop "$p" --chunk "$GRID_SUMSTATS_CHUNK" --work "$work" --replace "$GRID_REPLACE"
				tail -n 1 "$logdir/$trait/prepare.$p.log"
				grid_run_logged "$logdir/$trait/split.$p.log" python3 "$ROOT/f/0.common.py" split-sumstats --input "$std" --out-dir "$sumstats" --prefix "$p" --chrs "${CHRS[*]}"
				n=$(
					python3 - "$meta" "$GRID_N_GWAS" "$p" <<'PY'
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
				ng+=("$n")
				echo "  $trait.$p: n_gwas=$n"
			done
			infer_chr() {
				local c=$1 raw="$run/raw/chr$1" marker="$run/raw/chr$1/done" p
				mkdir -p "$raw"
				if [[ -f $marker && $GRID_REPLACE == FALSE ]]; then
					local valid=TRUE
					for p in "${POPS[@]}"; do [[ -s $raw/$p/$p.chr$c.pst_eff.txt ]] || valid=FALSE; done
					compgen -G "$raw/*META*pst_eff*.txt" >/dev/null || valid=FALSE
					[[ $GRID_POSTERIOR != TRUE || -s $raw/joint_posterior.h5 ]] || valid=FALSE
					[[ $valid != TRUE ]] || {
						echo "SKIP $trait chr$c: completed posterior"
						return
					}
				fi
				rm -f -- "$marker"
				local files=() cmd=()
				for p in "${POPS[@]}"; do files+=("$sumstats/$p.chr$c.tsv"); done
				cmd=(env OMP_NUM_THREADS="$GRID_THREADS" OPENBLAS_NUM_THREADS="$GRID_THREADS" MKL_NUM_THREADS="$GRID_THREADS" python3 "$prscx" --ref_dir="$ref" --bim_prefix="$GRID_CSX_BIM_PREFIX" --sst_file="$(join_comma "${files[@]}")" --n_gwas="$(join_comma "${ng[@]}")" --pop="$(join_comma "${POPS[@]}")" --chrom="$c" --n_iter="$GRID_MCMC_ITER" --n_burnin="$GRID_MCMC_BURNIN" --thin="$GRID_MCMC_THIN" --seed="$((GRID_SEED + c))" --out_dir="$raw" --out_name="$trait" --meta=TRUE)
				cmd+=(--write_pst="$GRID_POSTERIOR")
				[[ $GRID_PHI == auto ]] || cmd+=(--phi="$GRID_PHI")
				echo "RUN $trait chr$c: joint PRS-CSx"
				grid_run_logged "$logdir/$trait/csx.chr$c.log" "${cmd[@]}" || return $?
				for p in "${POPS[@]}"; do need "$raw/$p/$p.chr$c.pst_eff.txt"; done
				touch "$marker"
			}
			# Wait for each bounded batch explicitly: no lost child failures with wait -n.
			pids=()
			for c in "${CHRS[@]}"; do
				infer_chr "$c" &
				pids+=("$!")
				if ((${#pids[@]} >= GRID_JOBS)); then
					status=0
					for pid in "${pids[@]}"; do wait "$pid" || status=1; done
					((status == 0)) || _grid_die 'CSx inference failed'
					pids=()
				fi
			done
			status=0
			for pid in "${pids[@]}"; do wait "$pid" || status=1; done
			((status == 0)) || _grid_die 'CSx inference failed'
			for i in "${!POPS[@]}"; do
				p=${POPS[$i]}
				inputs=()
				for c in "${CHRS[@]}"; do
					w="$run/weights/$p.chr$c.tsv"
					grid_run python3 "$ROOT/f/1.csx.py" normalize-weights --input "$run/raw/chr$c/$p/$p.chr$c.pst_eff.txt" --output "$w"
					inputs+=("$w")
				done
				grid_run python3 "$io" weights "$run/weights/$p.csx.gz" "${inputs[@]}"
				publish "$run/weights/$p.csx.gz" "${finals[$i]}"
				printf '%s\n' "$sig" >"$run/signature"
				publish "$run/signature" "${finals[$i]}.signature"
				publish "$sumstats/$p.json" "${finals[$i]}.metadata.json"
			done
		fi
	fi
	if [[ $GRID_STAGE == weights ]]; then
		echo "DONE $trait: permanent SNP weights saved"
		return 0
	fi
	csx_score_run
	if [[ $GRID_POSTERIOR == TRUE ]]; then
		# Harmonized inputs are disposable; rebuild them without rerunning MCMC.
		mkdir -p "$sumstats"
		for i in "${!POPS[@]}"; do
			p=${POPS[$i]}
			# Validate the preparation signature even when a working table exists;
			# corrected EAF orientation must reach posterior centering on --stage score.
			grid_run_logged "$logdir/$trait/prepare.$p.log" python3 "$io" sumstats-cache --input "${gwas[$i]}" --output "$sumstats/$p.tsv.gz" --metadata "$sumstats/$p.json" --snpinfo "$GRID_CSX_SNPINFO" --trait "$trait" --pop "$p" --chunk "$GRID_SUMSTATS_CHUNK" --work "$work" --replace FALSE
			grid_run_logged "$logdir/$trait/split.$p.log" python3 "$io" split-sumstats --input "$sumstats/$p.tsv.gz" --out-dir "$sumstats" --prefix "$p" --chrs "${CHRS[*]}"
		done
		for c in "${CHRS[@]}"; do need "$run/raw/chr$c/joint_posterior.h5"; done
		posterior_args=(--raw-dir "$run/raw" --sumstats-dir "$sumstats" --target-dir "$GRID_TARGET_DIR"
			--output "$score_home/1csx.posterior.rds" --chrs "${CHRS[*]}" --threads "$GRID_THREADS"
			--memory "$GRID_POSTERIOR_MEMORY" --keep "$GRID_KEEP" --remove "$GRID_REMOVE" --replace "$GRID_REPLACE")
		[[ -z $GRID_POSTERIOR_FREQ_DIR ]] || posterior_args+=(--frequency-dir "$GRID_POSTERIOR_FREQ_DIR")
		grid_run_logged "$logdir/$trait/posterior.log" python3 "$ROOT/f/1.csx.py" posterior "${posterior_args[@]}"
		echo "DONE $trait: individual posterior covariance in $score_home/1csx.posterior.rds"
	fi
}

grid_pipeline csx "$@"

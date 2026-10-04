#!/usr/bin/env bash
# GRID step implementations; the public entry is ../3grid.sh.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
source "$ROOT/f/0.common.sh"
grid_activate_environment


# 🚩 grid_setup
grid_setup() {
	set -euo pipefail
	grid_configure "$@"
	grid_parse_args "$@"
	trait=${GRID_TRAIT,,}
	POPS=($(grid_csv_words "$GRID_POPS"))
	CHRS=($(grid_expand_chrs "$GRID_CHRS"))
	suffix=''
	[[ ${CHRS[*]} == "$(seq -s ' ' 1 22)" ]] || suffix="/chr$(
		IFS=,
		echo "${CHRS[*]}"
	)"
	[[ $trait == height || $trait == ldl || $trait == t2dm ]] || _grid_die 'Unknown trait'
	for c in "${CHRS[@]}"; do [[ $c =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || _grid_die "Invalid autosome: $c"; done
	out="$GRID_OUTPUT_ROOT/$trait$suffix"
	gdir="$out/grid"
	mkdir -p "$gdir" "$out/log" "$out/scores" "$gdir/weights" "$gdir/model"
}


# 🚩 grid_arg
grid_arg() (
	set -euo pipefail
	# shellcheck source=f/0.common.sh
	grid_configure "$@"
	grid_parse_args "$@"

	read -r -a CHRS <<<"$(grid_expand_chrs "$GRID_CHRS")"
	((${#CHRS[@]})) || _grid_die 'No chromosomes selected'
	case "$GRID_ARG_ACTION" in
		all | check) ;;
		prepare | infer | convert | features | affinity)
			echo "ERROR: ARG construction moved to refGen.sh data preparation." >&2
			echo "Run: bash /mnt/d/scripts/gu/arg.sh build --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --dir-pfile $GRID_ARG_HAP_DIR --map-dir ${GRID_ARG_MAP_DIR:-MAP_DIR} --ancestry-file $GRID_ANCESTRY_FILE --chr $GRID_CHRS" >&2
			exit 2
			;;
		*) _grid_die "bad --arg-action=$GRID_ARG_ACTION (GRID now supports check/all only)" ;;
	esac

	missing=0
	for c in "${CHRS[@]}"; do
		files=(
			"$GRID_ARG_OUT/argn/chr$c.argn"
			"$GRID_ARG_TREES_DIR/chr$c.trees"
			"$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv"
			"$GRID_ARG_TREES_DIR/chr$c.anchors.tsv"
			"$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
		)
		for f in "${files[@]}"; do
			if [[ ! -s $f ]]; then
				echo "ERROR: missing GRID ARG-Needle data for chr$c: $f" >&2
				missing=1
			fi
		done
	done
	if ((missing)); then
		echo "Prepare the reusable ARG data with:" >&2
		echo "  bash /mnt/d/scripts/gu/arg.sh build --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --dir-pfile $GRID_ARG_HAP_DIR --map-dir ${GRID_ARG_MAP_DIR:-MAP_DIR} --ancestry-file $GRID_ANCESTRY_FILE --chr $GRID_CHRS" >&2
		exit 2
	fi

	for c in "${CHRS[@]}"; do
		gzip -t "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
		python3 - "$GRID_ARG_TREES_DIR/chr$c.trees" "$c" <<'PY'
import sys,tskit
ts=tskit.load(sys.argv[1])
if min(ts.num_trees,ts.num_samples,ts.num_sites,ts.num_mutations) < 1:
    raise SystemExit(
        f"ERROR: invalid GRID ARG-Needle tree for chr{sys.argv[2]}: "
        f"trees={ts.num_trees} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}"
    )
print(
    f"GRID ARG CHECK PASS chr{sys.argv[2]} trees={ts.num_trees} "
    f"samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations}"
)
PY
	done
	echo "GRID ARG-Needle data ready: $GRID_ARG_OUT"
)


# 🚩 grid_ld
grid_ld() (
	set -euo pipefail
	grid_setup "$@"
	ld_hdf5_for() {
		local p=${1,,} c=$2
		local ref_type=1kg d f
		[[ $(basename -- "$GRID_CSX_SNPINFO") != snpinfo_mult_ukbb_hm3 ]] || ref_type=ukbb
		d="$GRID_CSX_REF_DIR/ldblk_${ref_type}_$p"
		[[ -d $d ]] || d="$GRID_CSX_REF_DIR/ldblk_${ref_type}_${p^^}"
		f="$d/ldblk_${ref_type}_chr$c.hdf5"
		[[ -s $f ]] || return 1
		printf '%s\n' "$f"
	}
	make_ld() {
		local p c f dest
		for p in "${POPS[@]}"; do
			p=${p^^}
			mkdir -p "$GRID_LDSCORE_DIR/$p"
			for c in "${CHRS[@]}"; do
				dest="$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz"
				f=$(ld_hdf5_for "$p" "$c")
				[[ -s $f ]] || _grid_die "No PRS-CSx HDF5 LD file for $p chr$c below $GRID_CSX_REF_DIR"
				key=$(python3 "$ROOT/f/0.common.py" signature "$p" "$c" --files "$f" "$ROOT/f/3grid.py")
				[[ -s $dest && -f $dest.signature && $(cat "$dest.signature") == "$key" && $GRID_REPLACE == FALSE ]] && continue
				grid_run_logged "$out/log/grid.ld.$p.chr$c.log" python3 "$ROOT/f/3grid.py" ld --hdf5 "$f" --pop "$p" --chr "$c" --out "$dest.tmp.gz"
				mv -f -- "$dest.tmp.gz" "$dest"
				printf '%s\n' "$key" >"$dest.signature"
			done
		done
	}

	make_ld
	echo "GRID ld completed: $gdir"
)


# 🚩 grid_transport
grid_transport() (
	set -euo pipefail
	grid_setup "$@"
	make_transport() {
		for c in "${CHRS[@]}"; do
			[[ -s $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz ]] || _grid_die "Missing ARG feature file $GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
			if grid_is_true "$GRID_REQUIRE_LD"; then for p in "${POPS[@]}"; do [[ -s $GRID_LDSCORE_DIR/${p^^}/chr$c.ldscore.tsv.gz ]] || _grid_die "Missing LD score for ${p^^} chr$c"; done; fi
		done
		t="$gdir/transport.tsv.gz"
		cmd=(python3 "$ROOT/f/3grid.py" transport --trait "$trait" --pops "$GRID_POPS" --chrs "${CHRS[*]}" --sumstats-dir "$out/sumstats/bychr" --arg-dir "$GRID_ARG_TREES_DIR" --ldscore-dir "$GRID_LDSCORE_DIR" --centers "$GRID_MED_FILE" --max-snps-per-chr "$GRID_MAX_SNPS_PER_CHR" --out "$t")
		[[ -z $GRID_EXTERNAL_AGE ]] || cmd+=(--external-age "$GRID_EXTERNAL_AGE")
		inputs=("$out/grid/inputs.signature" "$GRID_MED_FILE" "$ROOT/f/3grid.py")
		[[ -z $GRID_EXTERNAL_AGE ]] || inputs+=("$GRID_EXTERNAL_AGE")
		for c in "${CHRS[@]}"; do
			inputs+=("$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz")
			for p in "${POPS[@]}"; do
				inputs+=("$out/sumstats/bychr/$trait.$p.chr$c.tsv.gz")
				[[ ! -f $GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz ]] || inputs+=("$GRID_LDSCORE_DIR/$p/chr$c.ldscore.tsv.gz")
			done
		done
		key=$(python3 "$ROOT/f/0.common.py" signature "$GRID_MAX_SNPS_PER_CHR" "$GRID_REQUIRE_LD" --files "${inputs[@]}")
		if [[ -s $t && -f $t.signature && $(cat "$t.signature") == "$key" && $GRID_REPLACE == FALSE ]]; then
			echo 'SKIP matching GRID transport'
			return
		fi
		rm -f -- "$t.signature"
		grid_run_logged "$out/log/grid.transport.log" "${cmd[@]}"
		printf '%s\n' "$key" >"$t.signature"
	}

	make_transport
	echo "GRID transport completed: $gdir"
)


# 🚩 grid_fit
grid_fit() (
	set -euo pipefail
	grid_setup "$@"
	make_fit() {
		local t="$gdir/transport.tsv.gz"
		[[ -s $t ]] || _grid_die 'GRID transport table missing; rerun ./3grid.sh grid'
		grid_run_logged "$out/log/grid.transport.fit.log" python3 "$ROOT/f/3grid.py" fit --input "$t" --out-dir "$gdir/model" --ridge-alpha "$GRID_RIDGE_ALPHA" --model "$GRID_TRANSPORT_MODEL" --conservation-min "$GRID_CONSERVATION_MIN" --conservation-max "$GRID_CONSERVATION_MAX" --seed "$GRID_SEED"
	}

	make_fit
	echo "GRID fit completed: $gdir"
)


# 🚩 grid_weights
grid_weights() (
	set -euo pipefail
	grid_setup "$@"
	make_weights() {
		[[ -s $gdir/model/variant_conservation.tsv.gz ]] || _grid_die 'GRID model output missing; rerun ./3grid.sh grid'
		[[ -s $out/csx/manifest.tsv ]] || _grid_die 'Run ./1csx.sh, then ./3grid.sh grid first'
		grid_run_logged "$out/log/grid.weights.log" python3 "$ROOT/f/3grid.py" weights --weights-dir "$out/csx/weights" --conservation "$gdir/model/variant_conservation.tsv.gz" --manifest "$out/csx/manifest.tsv" --out-dir "$gdir/weights" --pops "$GRID_POPS" --chrs "${CHRS[*]}"
	}

	make_weights
	echo "GRID weights completed: $gdir"
)


# 🚩 grid_score
grid_score() (
	set -euo pipefail
	grid_setup "$@"
	score_one() {
		local label=$1 c=$2 mode prefix w o
		mode=$(grid_target_mode "$c" || true)
		[[ -n $mode ]] || _grid_die "Missing target genotype chr$c"
		prefix="$GRID_TARGET_DIR/chr$c"
		w="$gdir/weights/$label.chr$c.tsv"
		o="$out/scores/tmp/$label.chr$c"
		mkdir -p "$out/scores/tmp"
		[[ -s $w ]] || _grid_die "Missing $w"
		input_files=("$w" "$ROOT/f/3grid.f.sh" "$ROOT/f/0.common.py")
		if [[ $mode == pfile ]]; then
			input_files+=("$prefix.pgen" "$prefix.psam")
			if [[ -s $prefix.pvar ]]; then input_files+=("$prefix.pvar"); else input_files+=("$prefix.pvar.zst"); fi
		else input_files+=("$prefix.bed" "$prefix.bim" "$prefix.fam"); fi
		[[ -z $GRID_KEEP ]] || input_files+=("$GRID_KEEP")
		[[ -z $GRID_REMOVE ]] || input_files+=("$GRID_REMOVE")
		score_sig=$(python3 "$ROOT/f/0.common.py" signature --files "${input_files[@]}")
		[[ -s $o.sscore && -f $o.signature && $(cat "$o.signature") == "$score_sig" && $GRID_REPLACE == FALSE ]] && return
		rm -f -- "$o.signature"
		cmd=(plink2 "--$mode" "$prefix")
		[[ $mode != pfile || -s $prefix.pvar ]] || cmd+=(vzs)
		cmd+=(--score "$w" 1 2 3 header-read no-mean-imputation cols=+scoresums --threads "$GRID_THREADS" --out "$o")
		[[ -z $GRID_KEEP ]] || cmd+=(--keep "$GRID_KEEP")
		[[ ! -s $GRID_REMOVE ]] || cmd+=(--remove "$GRID_REMOVE")
		grid_run_logged "$out/log/grid.score.$label.chr$c.log" "${cmd[@]}"
		[[ -s $o.sscore ]] || _grid_die "Missing $o.sscore"
		printf '%s\n' "$score_sig" >"$o.signature"
	}
	make_scores() {
		labels=(GRID_shared)
		for p in "${POPS[@]}"; do labels+=("GRID_${p^^}"); done
		files=()
		for label in "${labels[@]}"; do
			pids=()
			status=0
			for c in "${CHRS[@]}"; do
				score_one "$label" "$c" &
				pids+=("$!")
				if ((${#pids[@]} >= GRID_JOBS)); then
					for pid in "${pids[@]}"; do wait "$pid" || status=1; done
					((status == 0)) || _grid_die "Scoring failed for $label"
					pids=()
				fi
			done
			for pid in "${pids[@]}"; do wait "$pid" || status=1; done
			((status == 0)) || _grid_die "Scoring failed for $label"
			inp=()
			for c in "${CHRS[@]}"; do inp+=("$out/scores/tmp/$label.chr$c.sscore"); done
			grid_run python3 "$ROOT/f/0.common.py" combine-scores --inputs "${inp[@]}" --name "$label" --output "$out/scores/$label.tsv.gz"
			files+=("$out/scores/$label.tsv.gz")
		done
		grid_run python3 "$ROOT/f/0.common.py" merge-scores --inputs "${files[@]}" --output "$out/scores/grid_population.tsv.gz"
		grid_run python3 "$ROOT/f/3grid.py" mix --scores "$out/scores/grid_population.tsv.gz" --ancestry "$GRID_ANCESTRY_FILE" --prefix GRID --out "$out/scores/grid_mixed.tsv.gz"
		grid_run python3 "$ROOT/f/0.common.py" merge-scores --inputs "$out/scores/grid_mixed.tsv.gz" "$out/scores/GRID_shared.tsv.gz" --output "$out/scores/grid.tsv.gz"
		# Matched/posterior PRS-CSx baselines use the same genotype-only ancestry probabilities.
		if [[ -s $out/scores/csx.tsv.gz ]]; then grid_run python3 "$ROOT/f/3grid.py" mix --scores "$out/scores/csx.tsv.gz" --ancestry "$GRID_ANCESTRY_FILE" --prefix CSX --out "$out/scores/csx_mixed.tsv.gz"; fi
	}

	make_scores
	echo "GRID score completed: $gdir"
)


# 🚩 grid_grid
grid_grid() (
	set -euo pipefail
	grid_setup "$@"
	if [[ $GRID_DRY_RUN == TRUE ]]; then
		echo 'PLAN GRID: arg check -> bridge 1csx inputs -> LD -> transport -> blocked fit -> weights -> scores'
		exit 0
	fi
	# Internal GRID steps run together; evaluation uses the separate Yeval.sh entry.
	exec {grid_lock}>"$gdir/run.lock"
	flock -n "$grid_lock" || _grid_die "Another GRID run is active: $gdir"
	rm -f "$gdir/GRID_RUN.txt"
	grid_arg "$@"
	grid_run python3 "$ROOT/f/3grid.py" inputs --trait "$trait" --gwas-dir "$GRID_GWAS_DIR" --snpinfo "$GRID_CSX_SNPINFO" --chrs "${CHRS[*]}" --out "$out"
	steps=(ld transport fit weights score)
	for step in "${steps[@]}"; do "grid_$step" "$@"; done
	cat >"$gdir/GRID_RUN.txt" <<META
created=$(date -Is)
trait=$trait
chromosomes=${CHRS[*]}
method=genealogy-informed shrinkage of PRS-CSx population effects toward a shared effect
validation=blocked out-of-fold transportability model
local_genealogy=$GRID_ARG_TREES_DIR/chrCHR.variants.tsv.gz
ld_source=PRS-CSx HDF5 reference panels
participant_phenotype_used_for_weights=FALSE
META
	echo "GRID completed: $out/scores/grid.tsv.gz"
)

command=${1:-help}
shift || true
case "$command" in
	run) grid_grid "$@" ;;
	arg | ld | transport | fit | weights | score) "grid_$command" "$@" ;;
	-h | --help | help) echo 'Usage: bash f/3grid.f.sh run|arg|ld|transport|fit|weights|score [shared options]' ;;
	*)
		echo "Unknown GRID step: $command" >&2
		exit 2
		;;
esac

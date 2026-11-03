#!/usr/bin/env bash

# Conda activate/deactivate hooks in common geospatial packages are not nounset-safe.
set -eo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=f/0.common.sh
source "$ROOT/f/0.common.sh"
source "$ROOT/../0f/0os.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
gu_load_config
ENV_NAME=${GU_ENV_NAME:-gu}
SOFT=${GU_SOFT:-$HOME/.local/share/gu}
STATE=${GU_INSTALL_STATE:-$SOFT/.gu-install.state}
YML=$ROOT/environment.yml
IBDMIX_DIR=${IBDMIX_RUNTIME:-$SOFT/ibdmix}
TRACE_DIR=${TRACE_RUNTIME:-$SOFT/trace}
DATA_ROOT=$GU_DATA_ROOT
REF_ROOT=$GU_REF_ROOT
TARGET_ROOT=${GU_TARGET_ROOT:-$REF_ROOT/1kg/37}
TARGET_VCF_DIR=${GU_TARGET_VCF_DIR:-$TARGET_ROOT/vcf}
SAMPLE_PANEL=${GU_SAMPLE_PANEL:-$REF_ROOT/1kg/38/samples.txt}
ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-${GU_ARCHAIC37_ROOT:-$REF_ROOT/archaic/37/vcf}}


# 🚩 Command-line help
usage() {
	cat <<'HELP'
Modules:
    Default  Install GU software and environment
    --arg  Install Threads and SINGER environments

Main options:
    --GPU TRUE|FALSE  Build optional IBDmix CUDA backend (default FALSE)
    --check  Check GU software
    --check-ibdmix  Check only the CPU IBDmix runtime (no installation)
    --check-references  Check reference inputs
    --dir-software DIR  Software installation root (see f/runtime.env)
    -h, --help  Show options and examples below

Example:
    ./install.sh                 # GU software and environment
    ./install.sh --arg           # Threads + SINGER, independent environments
    ./install.sh --check         # Check GU software
    ./install.sh --check-references

    ./install.sh --arg
    ./install.sh
HELP
}

install_threads() (
	# Independent of the Needle and GU environments (Threads pins arg-needle-lib).
	set -euo pipefail
	runtime=${GU_THREADS_ENV:-$HOME/.venvs/gu-threads}
	python=${GU_THREADS_PYTHON:-python3}
	[[ ! -x $runtime/bin/python ]] || python=$runtime/bin/python
	"$python" -c 'import sys; assert sys.version_info >= (3, 12), "Threads tested environment requires Python >=3.12; set GU_THREADS_PYTHON"'
	if [[ ! -d $runtime/conda-meta ]]; then "$python" -m venv "$runtime"; fi
	# Pinned threads environment; keep independent of the other ARG backend.
	"$runtime/bin/python" -m pip install -r /dev/stdin <<'GU_THREADS_REQUIREMENTS'
threads-arg==0.2.1
arg-needle-lib==1.1.3
pysam>=0.23,<0.24
tszip>=0.2,<1

GU_THREADS_REQUIREMENTS
	"$runtime/bin/python" -m pip freeze >"$runtime/requirements.lock.txt"
	"$runtime/bin/threads" infer --help
)

install_singer() (
	set -euo pipefail
	runtime=${GU_SINGER_ENV:-$HOME/.venvs/gu-singer}
	prefix=${GU_SINGER_HOME:-$HOME/software/singer-0.1.9}
	archive=${GU_SINGER_ARCHIVE:-$SOFT/sources/singer-0.1.9-beta-linux-x86_64.tar.gz}
    if [[ ! -s $archive ]]; then
        [[ -z ${GU_SINGER_ARCHIVE:-} ]] || { echo "missing $GU_SINGER_ARCHIVE" >&2; exit 2; }
        mkdir -p "$(dirname -- "$archive")"
        pending=$(mktemp --suffix=.tar.gz)
        trap 'rm -f "$pending"' EXIT
        /usr/bin/curl -fL --connect-timeout 20 --max-time 600 --retry 2 -o "$pending" https://raw.githubusercontent.com/popgenmethods/SINGER/eb8e39b1a15be4a9a4df4fdaab61847bf73515d7/releases/singer-0.1.9-beta-linux-x86_64.tar.gz
        printf '%s  %s\n' accd2cf432092aec64f27ddc88a25098974060f25d54e477c0603afbad8e91b8 "$pending" | sha256sum -c -
        cp "$pending" "$archive"
    fi
    printf '%s  %s\n' accd2cf432092aec64f27ddc88a25098974060f25d54e477c0603afbad8e91b8 "$archive" | sha256sum -c -
	mkdir -p "$prefix"
	tar -xzf "$archive" -C "$prefix" --strip-components=1
	if [[ ! -d $runtime/conda-meta ]]; then "${GU_SINGER_PYTHON:-python3}" -m venv "$runtime"; fi
	# Pinned singer environment; keep independent of the other ARG backend.
	"$runtime/bin/python" -m pip install -r /dev/stdin <<'GU_SINGER_REQUIREMENTS'
tskit>=1.0,<2
numpy>=2,<3
pysam==0.23.3
tszip==0.3.1
msprime>=1.3,<2
GU_SINGER_REQUIREMENTS
	"$runtime/bin/python" -m pip freeze >"$runtime/requirements.lock.txt"
	"$runtime/bin/python" "$prefix/singer_master" --help >/dev/null
	"$runtime/bin/python" "$prefix/convert_to_tskit" --help >/dev/null
	printf 'SINGER 0.1.9 installed: %s\nPython runtime: %s\n' "$prefix" "$runtime"
)

for option in "$@"; do
	case "$option" in -h | --help)
		usage
		exit 0
		;;
	esac
done
case "${1:-}" in
	--arg | --threads | --singer)
		[[ $# == 1 ]] || {
			echo "ERROR: $1 must be used alone; configure backend paths with GU_THREADS_ENV / GU_SINGER_ENV / GU_SINGER_HOME." >&2
			exit 2
		}
		case "$1" in
			--arg)
				install_threads
				install_singer
				;;
			--threads) install_threads ;;
			--singer) install_singer ;;
		esac
		exit 0
		;;
esac

need_value() { [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || {
	echo "ERROR: $1 requires a value" >&2
	exit 2
}; }

INSTALL_ACTION=install
CONDA_INPUT=""
BUILD_GPU=FALSE
IBDMIX_COMMIT=6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef
TRACE_COMMIT=70cab6c6687900826eb0b2ded69c39554b1b4146


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in
		--GPU)
			need_value "$@"
			case "${2^^}" in TRUE|FALSE) BUILD_GPU=${2^^} ;; *) echo "ERROR: --GPU requires TRUE or FALSE" >&2; exit 2 ;; esac
			shift 2 ;;
		--check-ibdmix)
			INSTALL_ACTION=check-ibdmix
			shift
			;;
		--check)
			INSTALL_ACTION=check
			shift
			;;
		--check-references)
			INSTALL_ACTION=check-references
			shift
			;;
		--env-name)
			need_value "$@"
			ENV_NAME=$2
			shift 2
			;;
		--conda)
			need_value "$@"
			CONDA_INPUT=$2
			shift 2
			;;
		--dir-software)
			need_value "$@"
			SOFT=$2
			STATE=$SOFT/.gu-install.state
			IBDMIX_DIR=$SOFT/ibdmix
			TRACE_DIR=$SOFT/trace
			shift 2
			;;
		--dir-ibdmix)
			need_value "$@"
			IBDMIX_DIR=$2
			shift 2
			;;
		--dir-trace)
			need_value "$@"
			TRACE_DIR=$2
			shift 2
			;;
		--dir-data-root)
			need_value "$@"
			DATA_ROOT=$2
			REF_ROOT=$DATA_ROOT/gen
			TARGET_ROOT=$REF_ROOT/1kg/37
			TARGET_VCF_DIR=$TARGET_ROOT/vcf
			SAMPLE_PANEL=$REF_ROOT/1kg/38/samples.txt
			ARCHAIC_ROOT=$REF_ROOT/archaic/37/vcf
			shift 2
			;;
		--dir-reference-root)
			need_value "$@"
			REF_ROOT=$2
			TARGET_ROOT=$REF_ROOT/1kg/37
			TARGET_VCF_DIR=$TARGET_ROOT/vcf
			SAMPLE_PANEL=$REF_ROOT/1kg/38/samples.txt
			ARCHAIC_ROOT=$REF_ROOT/archaic/37/vcf
			shift 2
			;;
		--dir-target-root)
			need_value "$@"
			TARGET_ROOT=$2
			TARGET_VCF_DIR=$TARGET_ROOT/vcf
			SAMPLE_PANEL=$TARGET_ROOT/samples.txt
			shift 2
			;;
		--dir-target-vcf)
			need_value "$@"
			TARGET_VCF_DIR=$2
			shift 2
			;;
		--sample-panel)
			need_value "$@"
			SAMPLE_PANEL=$2
			shift 2
			;;
		--dir-archaic-root)
			need_value "$@"
			ARCHAIC_ROOT=$2
			shift 2
			;;
		-h | --help)
			usage
			exit 0
			;;
		*)
			echo "ERROR: unknown option: $1" >&2
			usage >&2
			exit 2
			;;
	esac
done

CONDA_EXE_PATH=""
for candidate in "$CONDA_INPUT" "${GU_CONDA_EXE:-}" "${GU_CONDA_BASE:+$GU_CONDA_BASE/bin/conda}" "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" /opt/conda/bin/conda "$(command -v conda 2>/dev/null || true)"; do
	[[ -n $candidate && -x $candidate ]] || continue
	CONDA_EXE_PATH=$candidate
	break
done
[[ -n $CONDA_EXE_PATH ]] || {
	echo "ERROR: conda is unavailable; install Miniconda/Anaconda or use --conda FILE" >&2
	exit 1
}
CONDA_BASE_PATH=$("$CONDA_EXE_PATH" info --base 2>/dev/null) || {
	echo "ERROR: cannot determine Conda base from $CONDA_EXE_PATH" >&2
	exit 1
}
[[ -d $CONDA_BASE_PATH ]] || {
	echo "ERROR: invalid Conda base from $CONDA_EXE_PATH: $CONDA_BASE_PATH" >&2
	exit 1
}

die() {
	echo "ERROR: $*" >&2
	exit 1
}
have_env() { [[ -d "$GU_ENV_PREFIX/conda-meta" ]]; }
activate() {
	source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"
	conda activate "$GU_ENV_PREFIX"
	# The shared diagnostic wrapper was sourced before conda activation.
	GU_BCFTOOLS_BIN=$(type -P bcftools)
}
check() {
	local bad=0 cmd plugin_help
	local F=$ROOT/f
	[[ ${CONDA_PREFIX:-} == "$GU_ENV_PREFIX" ]] || activate
	export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"
	gu_check_publication || bad=1
	for cmd in python python3 bcftools bgzip tabix bedtools samtools plink2 bgenix phyml cmake; do
		command -v "$cmd" >/dev/null 2>&1 || {
			echo "MISSING command: $cmd" >&2
			bad=1
		}
	done
	plugin_help=$(bcftools plugin -l 2>/dev/null || true)
	grep -Fx fill-tags <<<"$plugin_help" >/dev/null || {
		echo "MISSING bcftools plugin: fill-tags" >&2
		bad=1
	}
	python - <<'PY' || bad=1
import importlib
for name in ('numpy','pandas','pysam','cyvcf2','pybedtools','zarr','tskit','tszip',
             'tsinfer','tsdate','bio2zarr'):
    importlib.import_module(name)
PY
	for cmd in trace-extract trace-infer trace-summarize; do command -v "$cmd" >/dev/null 2>&1 || {
		echo "MISSING command: $cmd" >&2
		bad=1
	}; done
	for cmd in ibdmix generate_gt gt_lods; do
		[[ -x "$IBDMIX_DIR/build/src/$cmd" ]] || {
			echo "MISSING IBDmix executable: $IBDMIX_DIR/build/src/$cmd" >&2
			bad=1
		}
	done
	((bad == 0))
}

find_archaic_vcf() {
	local ref=$1 chr=$2 d found
	for d in "$ARCHAIC_ROOT/$ref" "$ARCHAIC_ROOT/avcf/$ref"; do
		[[ -d $d ]] || continue
		found=$(find "$d" -maxdepth 1 -type f \( -name "*chr${chr}_*.vcf.gz" -o -name "*chr${chr}.*.vcf.gz" -o -name "*chr${chr}.vcf.gz" \) -print -quit 2>/dev/null)
		[[ -n $found ]] && {
			printf '%s\n' "$found"
			return 0
		}
	done
	return 1
}

check_references() {
	local bad=0 chr ref vcf vcf_base pbase
	[[ -s $SAMPLE_PANEL ]] || {
		echo "MISSING 1KG sample panel: $SAMPLE_PANEL" >&2
		bad=1
	}
	if [[ -s $SAMPLE_PANEL ]]; then
		awk 'BEGIN{FS="[[:space:]]+"}NR==1{for(i=1;i<=NF;i++)h[tolower($i)]=1;exit !(h["sample"]&&h["sex"])}' "$SAMPLE_PANEL" || {
			echo "INVALID samples.txt: sample and sex columns are required" >&2
			bad=1
		}
	fi
	for chr in {1..22} X; do
		vcf_base=$(gu_target_genotype_base "$TARGET_VCF_DIR/chr" "$chr")
		pbase=$(gu_target_genotype_base "$TARGET_ROOT/pfile/chr" "$chr")
		vcf=$vcf_base.vcf.gz
		if [[ ! -s $vcf || (! -s $vcf.tbi && ! -s $vcf.csi) ]]; then
			[[ -s $pbase.pgen && (-s $pbase.pvar || -s $pbase.pvar.zst) && -s $pbase.psam ]] || {
				echo "MISSING target chr$chr: expected indexed $vcf or complete pfile $pbase" >&2
				bad=1
			}
		fi
	done
	for ref in Altai Chagyr Denisova Denisova25 Vindija; do
		for chr in {1..22} X; do
			vcf=$(find_archaic_vcf "$ref" "$chr" || true)
			[[ -n $vcf && -s $vcf && (-s $vcf.tbi || -s $vcf.csi) ]] || {
				echo "MISSING indexed archaic VCF: $ref chr$chr under $ARCHAIC_ROOT" >&2
				bad=1
			}
		done
	done
	((bad == 0))
}

check_ibdmix() {
    local bad=0 cmd
    [[ ${CONDA_PREFIX:-} == "$GU_ENV_PREFIX" ]] || activate
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
    export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"
    for cmd in python3 bcftools bgzip tabix bedtools samtools plink2; do
        command -v "$cmd" >/dev/null 2>&1 || { echo "MISSING command: $cmd" >&2; bad=1; }
    done
    bcftools plugin -l | grep -Fx fill-tags >/dev/null || { echo "MISSING bcftools plugin: fill-tags" >&2; bad=1; }
    python3 -c 'import numpy, pandas, pysam' || bad=1
    for cmd in "$IBDMIX_DIR/build/src/ibdmix" "$IBDMIX_DIR/build/src/generate_gt" "$IBDMIX_DIR/build/src/gt_lods" "${IBDMIX_PREPARE_BIN:-$SOFT/ibdmix-prepare/prepare}"; do
        [[ -x $cmd ]] || { echo "MISSING executable: $cmd" >&2; bad=1; continue; }
        if ldd "$cmd" 2>/dev/null | grep -q 'not found'; then
            echo "MISSING shared library: $cmd" >&2; bad=1
        fi
    done
    ((bad == 0))
}

if [[ $INSTALL_ACTION == check-ibdmix ]]; then
    check_ibdmix && { echo 'OK (IBDmix CPU software; reference inputs checked separately)'; exit 0; }
    exit 1
fi

mkdir -p "$SOFT"

if [[ $INSTALL_ACTION == check ]]; then
	check && echo OK
	exit
fi
if [[ $INSTALL_ACTION == check-references ]]; then
	check_references && echo OK
	exit
fi

# Reuse only this pinned release; old marker files do not certify the source version.
if [[ -s $STATE ]] && grep -Fxq "trace=$TRACE_COMMIT" "$STATE" && grep -Fxq "ibdmix=$IBDMIX_COMMIT" "$STATE" && [[ $BUILD_GPU == FALSE ]] && have_env && check; then
	echo OK
	exit 0
fi

# HPC compute nodes have no external DNS. A prepared project environment can
# be reused explicitly for the offline build; the health check still runs.
if [[ ${GU_INSTALL_ENV_READY:-0} != 1 ]]; then
	if have_env; then
		"$CONDA_EXE_PATH" env update -p "$GU_ENV_PREFIX" -f "$YML"
	else
		"$CONDA_EXE_PATH" env create -p "$GU_ENV_PREFIX" -f "$YML"
	fi
fi
activate

# Immutable upstream source snapshots. Never overwrite local upstream edits.
checkout_pinned() {
	local repo=$1 destination=$2 commit=$3
	if [[ ! -d $destination/.git ]]; then
		[[ ! -e $destination ]] || die "non-git directory already exists: $destination"
		git clone "$repo" "$destination"
	fi
	[[ -z $(git -C "$destination" status --porcelain --untracked-files=no) ]] || die "local edits exist in $destination; preserve them before installing the pinned version"
	if ! git -C "$destination" cat-file -e "$commit^{commit}" 2>/dev/null; then
		git -C "$destination" fetch origin "$commit"
	fi
	git -C "$destination" checkout --detach "$commit"
}
checkout_pinned https://github.com/YulinZhang9806/trace.git "$TRACE_DIR" "$TRACE_COMMIT"
# Reinstall from the verified source, including its declared CPU dependencies.
if [[ ${GU_INSTALL_ENV_READY:-0} == 1 ]] && python -c 'from importlib.metadata import version; assert version("tracehmm") == "0.0.1a0"' >/dev/null 2>&1; then
    echo "Reusing TRACE installed during pinned online preparation"
else
    python -m pip install --no-build-isolation "$TRACE_DIR"
fi
checkout_pinned https://github.com/PrincetonUniversity/IBDmix.git "$IBDMIX_DIR" "$IBDMIX_COMMIT"
ibdmix_cmake_args=()
if [[ -d $SOFT/sources/deps/cli11-2.2.0 ]]; then
    ibdmix_cmake_args+=("-DFETCHCONTENT_SOURCE_DIR_CLI11=$SOFT/sources/deps/cli11-2.2.0")
fi
cmake --fresh -S "$IBDMIX_DIR" -B "$IBDMIX_DIR/build" "${ibdmix_cmake_args[@]}"
cmake --build "$IBDMIX_DIR/build" --parallel "${GU_BUILD_JOBS:-4}"

IBDMIX_PREPARE_BIN="${IBDMIX_PREPARE_BIN:-$SOFT/ibdmix-prepare/prepare}" bash "$ROOT/f/ibdmix.sh" build-prepare --if-needed

if [[ $BUILD_GPU == TRUE ]]; then
	GU_SOFT="$SOFT" IBDMIX_GPU_BIN="${IBDMIX_GPU_BIN:-$SOFT/ibdmix-gpu/ibdmix-gpu}" bash "$ROOT/f/ibdmix.sh" build-gpu
fi
check || die "installation finished but the software health check still reports missing components"
printf 'ibdmix=%s\ntrace=%s\n' "$IBDMIX_COMMIT" "$TRACE_COMMIT" >"$STATE"
"$CONDA_EXE_PATH" list -n "$ENV_NAME" --explicit >"$SOFT/gu-conda-explicit.txt"
python -m pip freeze >"$SOFT/gu-pip-freeze.txt"
echo "OK (software). Run ./install.sh --check-references separately for the full reference inventory."

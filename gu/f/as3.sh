#!/usr/bin/env bash


# 🚩 as3_prepare
as3_prepare() (
	# Validate direct ArchaicSeeker3 Ref1028 inputs and write a lightweight manifest.
	set -euo pipefail

	need() { command -v "$1" >/dev/null 2>&1 || {
		echo "FATAL ERROR: missing command: $1" >&2
		exit 1
	}; }
	for cmd in bcftools awk sort comm cmp mktemp readlink stat; do need "$cmd"; done
	AS3_BCFTOOLS_BIN=$(type -P bcftools)
	bcftools() {
		"$AS3_BCFTOOLS_BIN" "$@" 2> >(sed '/^\[W::bcf_hdr_check_sanity\] AC should be declared as Number=A$/d' >&2)
	}

	dir0=${dir0:-/mnt/d}
	dir_ref=${dir_ref:-/mnt/f/gen}
	reference_dir=${AS3_REFERENCE_PANEL_DIR:-$dir_ref/archaic/38/vcf}
	reference_map=${AS3_REFERENCE_MAP:-}
	out=${AS3_DATA_OUT:?AS3_DATA_OUT is required; run this internal helper through gu.sh}
	target_dir=${AS3_TARGET_VCF_DIR:-${AS3_MODERN_HUMAN_VCF_DIR:-$dir_ref/1kg/${GRCH:-38}/vcf}}
	target_template=${AS3_TARGET_VCF_TEMPLATE:-}
	reference_template=${AS3_REFERENCE_VCF_TEMPLATE:-}
	[[ -n $target_template ]] || target_template='chr{CHR}.vcf.gz'
	[[ -n $reference_template ]] || reference_template='Ref_Panel.chr{CHR}.vcf.gz'
	panel_id=${AS3_REFERENCE_PANEL_ID:-Ref1028}
	genome_build=${AS3_GENOME_BUILD:-b${GRCH:-38}}
	requested=${AS3_CHRS:-"$(seq 1 22)"}
	check_only=${AS3_PREP_CHECK_ONLY:-0}

	while (($#)); do
		case "$1" in
			--reference-panel-dir)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --reference-panel-dir requires DIR" >&2
					exit 2
				}
				reference_dir=$2
				shift 2
				;;
			--reference-map)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --reference-map requires FILE" >&2
					exit 2
				}
				reference_map=$2
				shift 2
				;;
			--out)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --out requires DIR" >&2
					exit 2
				}
				out=$2
				shift 2
				;;
			--chr)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --chr requires LIST" >&2
					exit 2
				}
				requested=${2//,/ }
				shift 2
				;;
			--target-vcf-dir)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --target-vcf-dir requires DIR" >&2
					exit 2
				}
				target_dir=$2
				shift 2
				;;
			--target-vcf-template)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --target-vcf-template requires TEMPLATE" >&2
					exit 2
				}
				target_template=$2
				shift 2
				;;
			--reference-vcf-template)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --reference-vcf-template requires TEMPLATE" >&2
					exit 2
				}
				reference_template=$2
				shift 2
				;;
			--genome-build)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: --genome-build requires b38" >&2
					exit 2
				}
				genome_build=$2
				shift 2
				;;
			--check-only)
				check_only=1
				shift
				;;
			--threads | --max-procs | --tmp-root)
				[[ $# -ge 2 ]] || {
					echo "FATAL ERROR: $1 requires a value" >&2
					exit 2
				}
				shift 2
				;;
			--allow-chrx)
				echo "FATAL ERROR: official Ref1028 supports GRCh38 chr1-22 only" >&2
				exit 2
				;;
			*)
				echo "FATAL ERROR: unknown AS3 input option: $1" >&2
				exit 2
				;;
		esac
	done

	case "${genome_build,,}" in
		b38 | 38 | grch38 | hg38) genome_build=b38 ;;
		*)
			echo "FATAL ERROR: official Ref1028 is GRCh38 and supports chr1-22 only; got $genome_build" >&2
			exit 2
			;;
	esac
	[[ $check_only == 0 || $check_only == 1 ]] || {
		echo "FATAL ERROR: AS3_PREP_CHECK_ONLY must be 0 or 1" >&2
		exit 2
	}
	[[ -d $reference_dir ]] || {
		echo "FATAL ERROR: Ref1028 directory missing: $reference_dir" >&2
		exit 1
	}
	[[ -d $target_dir ]] || {
		echo "FATAL ERROR: target VCF directory missing: $target_dir" >&2
		exit 1
	}
	reference_map=${reference_map:-$reference_dir/Ref_Panel.map.txt}
	[[ -s $reference_map ]] || {
		echo "FATAL ERROR: Ref1028 map missing or empty: $reference_map" >&2
		exit 1
	}

	render() {
		local dir=$1 template=$2 chr=$3 path
		path=$(awk -v value="$template" -v chr="$chr" 'BEGIN{
      token="{CHR}";p=index(value,token)
      if(!p){token="{chr}";p=index(value,token)}
      if(p)value=substr(value,1,p-1) chr substr(value,p+length(token))
      print value
    }')
		[[ $path == /* ]] && printf '%s\n' "$path" || printf '%s/%s\n' "$dir" "$path"
	}

	valid_vcf() {
		local file=$1 count
		[[ -s $file && (-s $file.tbi || -s $file.csi) ]] || return 1
		bcftools view -h "$file" >/dev/null 2>&1 || return 1
		count=$(bcftools index -n "$file" 2>/dev/null || true)
		[[ $count =~ ^[1-9][0-9]*$ ]]
	}

	validate_single_record_contig() {
		local file=$1 chr=$2 role=$3
		bcftools index -s "$file" | awk -v expected="$chr" -v role="$role" -v file="$file" '
      $3+0>0 {
        c=$1; sub(/^chr/,"",c); n++
        if(c!=expected){print "FATAL ERROR: " role " VCF has records on unexpected contig " $1 ": " file > "/dev/stderr"; bad=1}
      }
      END{
        if(n!=1){print "FATAL ERROR: " role " VCF must contain records on exactly one chromosome (found " n "): " file > "/dev/stderr"; bad=1}
        exit bad
      }'
	}

	tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/gu-as3-ref1028.XXXXXX")
	trap 'rm -rf -- "$tmp_root"' EXIT
	awk 'BEGIN{FS="[[:space:]]+"}
  NF {
    if(NF!=2){print "FATAL ERROR: Ref1028 map must have exactly two columns at line " NR > "/dev/stderr"; bad=1}
    if($2!="AFR"&&$2!="DEN"&&$2!="NEAN"){print "FATAL ERROR: invalid Ref1028 ancestry label " $2 " at line " NR > "/dev/stderr"; bad=1}
    print $1
    n[$2]++
  }
  END{
    if(n["AFR"]!=146||n["DEN"]!=1||n["NEAN"]!=3){print "FATAL ERROR: Ref1028 map requires exactly AFR=146, DEN=1, and NEAN=3" > "/dev/stderr"; bad=1}
    exit bad
  }' "$reference_map" >"$tmp_root/map.samples.unsorted"
	sort "$tmp_root/map.samples.unsorted" >"$tmp_root/map.samples"
	duplicate=$(uniq -d "$tmp_root/map.samples" | sed -n '1p')
	[[ -z $duplicate ]] || {
		echo "FATAL ERROR: duplicate Ref1028 map sample: $duplicate" >&2
		exit 1
	}
	read -r n_afr n_den n_nean < <(awk 'BEGIN{FS="[[:space:]]+"}NF{n[$2]++}END{print n["AFR"]+0,n["DEN"]+0,n["NEAN"]+0}' "$reference_map")

	mapfile -t chromosomes < <(tr ', ' '\n\n' <<<"$requested" | awk 'NF{c=$1;sub(/^chr/,"",c);sub(/^CHR/,"",c);if(c==23)c="X";if(!seen[c]++)print c}')
	((${#chromosomes[@]} > 0)) || {
		echo "FATAL ERROR: AS3 chromosome list is empty" >&2
		exit 2
	}
	for chr in "${chromosomes[@]}"; do
		[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || {
			echo "FATAL ERROR: official Ref1028 supports GRCh38 chr1-22 only; got chr$chr" >&2
			exit 2
		}
	done

	mkdir -p "$out"
	manifest_tmp=$out/.manifest.$$.tsv
	release_tmp=$out/.release.$$.tsv
	printf 'panel_id\tunit\ttarget_vcf\treference_vcf\treference_map\treference_target_overlap\n' >"$manifest_tmp"
	printf 'key\tvalue\npanel_id\t%s\ngenome_build\tGRCh38\nreference_dir\t%s\nreference_map\t%s\nmap_AFR\t%s\nmap_DEN\t%s\nmap_NEAN\t%s\n' \
		"$panel_id" "$reference_dir" "$reference_map" "$n_afr" "$n_den" "$n_nean" >"$release_tmp"

	for chr in "${chromosomes[@]}"; do
		target=$(render "$target_dir" "$target_template" "$chr")
		reference=$(render "$reference_dir" "$reference_template" "$chr")
		valid_vcf "$target" || {
			echo "FATAL ERROR: target chr$chr VCF/index is missing, empty, or unreadable: $target" >&2
			exit 1
		}
		valid_vcf "$reference" || {
			echo "FATAL ERROR: Ref1028 chr$chr VCF/index is missing, empty, or unreadable: $reference" >&2
			exit 1
		}
		# A run directory may expose a native VCF through a short-lived symlink.
		# Persist only the canonical source in the manifest so the AS3 worker never
		# depends on that directory after input preparation.
		target=$(readlink -f -- "$target")
		reference=$(readlink -f -- "$reference")
		validate_single_record_contig "$target" "$chr" target
		validate_single_record_contig "$reference" "$chr" Ref1028

		bcftools query -l "$reference" | sort >"$tmp_root/reference.chr${chr}.samples"
		cmp -s "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" || {
			echo "FATAL ERROR: Ref1028 chr$chr VCF samples do not exactly match $reference_map" >&2
			echo "Only in VCF:" >&2
			comm -23 "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" | sed -n '1,10p' >&2
			echo "Only in map:" >&2
			comm -13 "$tmp_root/reference.chr${chr}.samples" "$tmp_root/map.samples" | sed -n '1,10p' >&2
			exit 1
		}
		target_samples_tmp=$tmp_root/target.chr${chr}.samples
		bcftools query -l "$target" | sort >"$target_samples_tmp"
		[[ -s $target_samples_tmp ]] || {
			echo "FATAL ERROR: target chr$chr VCF has no samples: $target" >&2
			exit 1
		}
		target_duplicate=$(uniq -d "$target_samples_tmp" | sed -n '1p')
		[[ -z $target_duplicate ]] || {
			echo "FATAL ERROR: duplicate target chr$chr sample: $target_duplicate" >&2
			exit 1
		}
		overlap=$(comm -12 "$target_samples_tmp" "$tmp_root/map.samples" | awk 'END{print NR+0}')
		printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$panel_id" "$chr" "$target" "$reference" "$reference_map" "$overlap" >>"$manifest_tmp"
		target_n=$(awk 'END{print NR+0}' "$target_samples_tmp")
		ref_n=$(awk 'END{print NR+0}' "$tmp_root/map.samples")
		target_variants=$(bcftools index -n "$target")
		ref_variants=$(bcftools index -n "$reference")
		printf '[AS3 INPUT chr%s] target_samples=%s reference_samples=%s target_variants=%s reference_variants=%s target_reference_sample_overlap=%s\n' \
			"$chr" "$target_n" "$ref_n" "$target_variants" "$ref_variants" "$overlap"
	done

	if [[ ! -s $out/manifest.tsv ]] || ! cmp -s "$manifest_tmp" "$out/manifest.tsv"; then mv -f "$manifest_tmp" "$out/manifest.tsv"; else rm -f "$manifest_tmp"; fi
	if [[ ! -s $out/release.tsv ]] || ! cmp -s "$release_tmp" "$out/release.tsv"; then mv -f "$release_tmp" "$out/release.tsv"; else rm -f "$release_tmp"; fi
	printf 'AS3 REF1028 INPUT CHECK PASSED: build=GRCh38 panel=%s chromosomes=%s map(AFR=%s,DEN=%s,NEAN=%s) manifest=%s mode=%s\n' \
		"$panel_id" "${chromosomes[*]}" "$n_afr" "$n_den" "$n_nean" "$out/manifest.tsv" "$([[ $check_only == 1 ]] && printf check || printf prepare)"
)

if [[ ${1:-} == prepare ]]; then
	shift
	as3_prepare "$@"
	exit
fi
set -euo pipefail

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
	cat <<'HELP'
Internal ArchaicSeeker3 runner. Use ./gu.sh as3 [check] --chr LIST; preparation
is automatic. AS3 requires CUDA and supports GRCh38 autosomes 1-22 only.
Use ./gu.sh as3 --replace-as3 TRUE to overwrite completed tasks.
Target samples are processed in memory-bounded chunks (default 64); use
--as3-target-chunk-size through gu.sh or AS3_TARGET_CHUNK_SIZE to override.
If inference completes but canonical post-processing fails or is interrupted,
the raw inference checkpoint is retained and the next identical run resumes at
post-processing instead of repeating GPU inference.
Each inference child is isolated in a cgroup with AS3_MEMORY_HIGH,
AS3_MEMORY_MAX, and AS3_MEMORY_SWAP_MAX (defaults: 20G, 24G, and 8G).
HELP
	exit 0
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
action=${GU_ACTION:-as3_run}
runtime=${AS3_RUNTIME:-$script_dir/as3}
model_dir=${AS3_MODEL_DIR:-${dir_ref:-/mnt/f/gen}/archaic/38/models}
data=${AS3_DATA_IN:?AS3_DATA_IN is required; run this internal helper through gu.sh}
out=${AS3_OUT:-${GU_ANALYSIS_ROOT:-/mnt/d/analysis/gu/1kg}/as3}
log_root=/tmp/gu-logs/as3/$(printf '%s' "$out" | sha256sum | cut -c1-16)
runtime_work=${AS3_RUNTIME_WORK_DIR:?AS3_RUNTIME_WORK_DIR is required; run this internal helper through gu.sh}
gpus=${AS3_GPUS:-0}
merge=${AS3_MERGE:-10000}
replace=${AS3_REPLACE:-0}
target_chunk_size=${AS3_TARGET_CHUNK_SIZE:-64}
memory_high=${AS3_MEMORY_HIGH:-20G}
memory_max=${AS3_MEMORY_MAX:-24G}
memory_swap_max=${AS3_MEMORY_SWAP_MAX:-8G}
as3_python=${AS3_PYTHON:-}
as3_env_name=${AS3_ENV_NAME:-as3_mamba}
if [[ -z $as3_python ]]; then
	conda_base=$(conda info --base 2>/dev/null || true)
	for cand in \
		"${conda_base:+$conda_base/envs/$as3_env_name/bin/python}" \
		"${GU_SOFT:-/mnt/d/software/gu}/mambaforge/envs/$as3_env_name/bin/python" \
		"${CONDA_PREFIX:-}/../$as3_env_name/bin/python" \
		"$(command -v python3 2>/dev/null || true)"; do
		[[ -x $cand ]] && {
			as3_python=$cand
			break
		}
	done
fi
exe=$runtime/ArchaicSeeker3.1-mamba
base_model_args=$model_dir/Basemodel_Mamba_4096/args.pckl
base_model_cp=$model_dir/Basemodel_Mamba_4096/best_model.pth
smoother_model_args=$model_dir/Smoother_512_Kernel_8192/args.pckl
smoother_model_cp=$model_dir/Smoother_512_Kernel_8192/best_model.pth
chr_filter=${AS3_CHRS:-}
manifest=$data/manifest.tsv

[[ $action == as3_run || $action == as3_check || $action == as3_env_check ]] || {
	echo "ERROR unknown AS3 action: $action" >&2
	exit 1
}
[[ $replace == 0 || $replace == 1 ]] || {
	echo "ERROR: internal replace-as3 state must be 0 or 1" >&2
	exit 1
}
[[ $merge =~ ^[1-9][0-9]*$ ]] || {
	echo "ERROR AS3_MERGE must be a positive integer" >&2
	exit 1
}
[[ ${AS3_STRIDE:-512} =~ ^[1-9][0-9]*$ ]] || {
	echo "ERROR AS3_STRIDE must be a positive integer" >&2
	exit 1
}
[[ ${AS3_ANC:-0} =~ ^[0-9]+$ ]] || {
	echo "ERROR AS3_ANC must be a non-negative integer" >&2
	exit 1
}
if [[ ! $target_chunk_size =~ ^[1-9][0-9]*$ ]]; then
	echo "ERROR AS3_TARGET_CHUNK_SIZE must be a positive integer" >&2
	exit 1
fi
for memory_setting in "$memory_high" "$memory_max" "$memory_swap_max"; do
	[[ $memory_setting =~ ^[1-9][0-9]*[KMGTPE]$ ]] || {
		echo "ERROR AS3 cgroup memory values must use an integer plus K/M/G/T/P/E (for example 24G): $memory_setting" >&2
		exit 1
	}
done
[[ -x $as3_python ]] || {
	echo "ERROR: AS3 Python was not found; run ./install.sh" >&2
	exit 1
}
for cmd in bcftools sort awk systemctl systemd-run stat; do command -v "$cmd" >/dev/null || {
	echo "ERROR missing command: $cmd" >&2
	exit 1
}; done
AS3_BCFTOOLS_BIN=$(type -P bcftools)
bcftools() {
	"$AS3_BCFTOOLS_BIN" "$@" 2> >(sed '/^\[W::bcf_hdr_check_sanity\] AC should be declared as Number=A$/d' >&2)
}
[[ $(stat -fc %T /sys/fs/cgroup 2>/dev/null) == cgroup2fs ]] || {
	echo "ERROR AS3 requires cgroup v2 for process-level memory isolation" >&2
	exit 1
}
systemctl --user is-system-running >/dev/null 2>&1 || {
	echo "ERROR AS3 requires a running per-user systemd manager for memory isolation" >&2
	exit 1
}
[[ -f $exe ]] || {
	echo "ERROR AS3 executable missing: $exe" >&2
	exit 1
}
for required in \
	"$runtime/src/dataloaders/__init__.py" \
	"$runtime/src/models/__init__.py" \
	"$runtime/src/stepsagnostic/__init__.py" \
	"$runtime/src/standard_postprocess.py" \
	"$base_model_args" \
	"$base_model_cp" \
	"$smoother_model_args" \
	"$smoother_model_cp"; do
	[[ -s $required ]] || {
		echo "ERROR bundled AS3 runtime file missing: $required" >&2
		exit 1
	}
done
AS3_EXPECT_PYTHON=${AS3_PYTHON_VERSION:-3.9} "$as3_python" "$script_dir/as3.py" health || {
	echo "ERROR AS3 Python dependency/Torch/CUDA health check failed: $as3_python; run ./install.sh --repair-as3" >&2
	exit 1
}
"$as3_python" "$script_dir/as3.py" model-check --runtime "$runtime" --model-dir "$model_dir" || {
	echo "ERROR AS3 bundled model definitions/checkpoints failed to load" >&2
	exit 1
}
read -r -a gpu_array <<<"$gpus"
((${#gpu_array[@]} > 0)) || {
	echo "ERROR AS3_GPUS is empty" >&2
	exit 1
}
for gpu in "${gpu_array[@]}"; do [[ $gpu =~ ^[0-9]+$ ]] || {
	echo "ERROR invalid GPU ID: $gpu" >&2
	exit 1
}; done
if [[ $action == as3_env_check ]]; then
	"$as3_python" - "${gpu_array[@]}" <<'PY'
import sys, torch
ids = [int(x) for x in sys.argv[1:]]
available = torch.cuda.is_available()
n = torch.cuda.device_count()
print("CUDA available:", available)
print("CUDA GPU count:", n)
if available:
    print("CUDA GPUs:", [torch.cuda.get_device_name(i) for i in range(n)])
bad = [x for x in ids if x >= n]
if bad:
    raise SystemExit(f"ERROR: requested GPU IDs {bad} but only {n} CUDA devices are visible")
if not available:
    raise SystemExit("ERROR: CUDA is unavailable; AS3 requires a CUDA GPU")
PY
	echo "AS3 ENV CHECK PASSED: python=$as3_python runtime=$runtime requested_GPUs=${gpu_array[*]} target_chunk_size=$target_chunk_size cgroup_memory_high=$memory_high cgroup_memory_max=$memory_max cgroup_swap_max=$memory_swap_max"
	exit 0
fi
[[ -s $manifest ]] || {
	echo "ERROR AS3 manifest missing after automatic preparation: $manifest" >&2
	exit 1
}
mkdir -p "$log_root" "$out/results"

unit_selected() {
	local unit=$1 c
	[[ -z $chr_filter ]] && return 0
	for c in $chr_filter; do
		c=${c#chr}
		if [[ $c == X || $c == 23 ]]; then
			[[ $unit == X* ]] && return 0
		elif [[ $unit == "$c" ]]; then
			return 0
		fi
	done
	return 1
}

# Load selected direct Ref1028 manifest rows.
tasks=()
declare -A seen_keys=()
while IFS=$'\t' read -r panel unit target ref map sample_overlap; do
	[[ $panel == panel_id ]] && continue
	unit_selected "$unit" || continue
	key=${panel}.chr${unit}
	[[ -z ${seen_keys[$key]:-} ]] || {
		echo "ERROR duplicate AS3 manifest key: $key" >&2
		exit 1
	}
	seen_keys[$key]=1
	tasks+=("$panel"$'\t'"$unit"$'\t'"$target"$'\t'"$ref"$'\t'"$map"$'\t'"$sample_overlap")
done <"$manifest"
((${#tasks[@]} > 0)) || {
	echo "ERROR no manifest tasks selected" >&2
	exit 1
}

validate_task() {
	local panel=$1 unit=$2 target=$3 ref=$4 map=$5 sample_overlap=$6 key work nt nr duplicate
	key=${panel}.chr${unit}
	for f in "$target" "$ref" "$map"; do [[ -s $f ]] || {
		echo "ERROR $key missing input: $f" >&2
		return 1
	}; done
	bcftools view -h "$target" >/dev/null 2>&1 || {
		echo "ERROR unreadable target VCF: $target" >&2
		return 1
	}
	bcftools view -h "$ref" >/dev/null 2>&1 || {
		echo "ERROR unreadable reference VCF: $ref" >&2
		return 1
	}
	nt=$(bcftools index -n "$target" 2>/dev/null || echo 0)
	nr=$(bcftools index -n "$ref" 2>/dev/null || echo 0)
	((nt > 0 && nr > 0)) || {
		echo "ERROR empty/unindexed VCF for $key" >&2
		return 1
	}
	work=$(mktemp -d "${TMPDIR:-/tmp}/as3_check_${key}.XXXXXX")
	bcftools query -l "$target" | sort >"$work/target.samples"
	bcftools query -l "$ref" | sort -u >"$work/ref.samples"
	awk 'NF{print $1}' "$map" | sort -u >"$work/map.samples"
	[[ -s $work/target.samples ]] || {
		echo "ERROR target VCF has no samples for $key" >&2
		rm -rf "$work"
		return 1
	}
	duplicate=$(uniq -d "$work/target.samples" | sed -n '1p')
	[[ -z $duplicate ]] || {
		echo "ERROR duplicate target sample for $key: $duplicate" >&2
		rm -rf "$work"
		return 1
	}
	cmp -s "$work/ref.samples" "$work/map.samples" || {
		echo "ERROR reference.map mismatch for $key" >&2
		rm -rf "$work"
		return 1
	}
	rm -rf "$work"
	[[ $sample_overlap =~ ^[0-9]+$ ]] || {
		echo "ERROR invalid target/reference sample-overlap count for $key: $sample_overlap" >&2
		return 1
	}
	echo "[AS3 $key] validated direct inputs: target_variants=$nt reference_variants=$nr target_reference_sample_overlap=$sample_overlap"
}

status=0
for task in "${tasks[@]}"; do
	IFS=$'\t' read -r panel unit target ref map sample_overlap <<<"$task"
	validate_task "$panel" "$unit" "$target" "$ref" "$map" "$sample_overlap" || status=1
done
((status == 0)) || {
	echo "ERROR AS3 manifest/input validation failed" >&2
	exit 1
}

if [[ $action == as3_check ]]; then
	"$as3_python" - <<'PY'
import torch
print("CUDA available:", torch.cuda.is_available())
print("CUDA GPU count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("CUDA GPUs:", [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())])
PY
	echo "CHECK PASSED: manifest=$manifest tasks=${#tasks[@]} requested_GPUs=${gpu_array[*]}"
	printf '%s\n' "${tasks[@]}" | awk -F'\t' '{print "  panel=" $1 " unit=" $2 " target=" $3 " reference=" $4}'
	exit 0
fi

"$as3_python" - "${gpu_array[@]}" <<'PY'
import sys, torch
if not torch.cuda.is_available():
    raise SystemExit("ERROR: CUDA is unavailable; AS3 requires a CUDA GPU")
n = torch.cuda.device_count()
ids = [int(x) for x in sys.argv[1:]]
bad = [x for x in ids if x >= n]
if bad:
    raise SystemExit(f"ERROR: requested GPU IDs {bad} but only {n} CUDA devices are visible")
print("CUDA GPUs:", n, [torch.cuda.get_device_name(i) for i in range(n)])
PY

completion_record() {
	local target=$1 ref=$2 map=$3
	printf 'schema\t2\ntarget\t%s\nreference\t%s\nmap\t%s\nexecutable\t%s\nbase_model_args\t%s\nbase_model\t%s\nsmoother_model_args\t%s\nsmoother_model\t%s\nmerge\t%s\nstride\t%s\nanc\t%s\ntarget_chunk_size\t%s\n' \
		"$(stat -c '%n:%s:%Y' "$target")" \
		"$(stat -c '%n:%s:%Y' "$ref")" \
		"$(stat -c '%n:%s:%Y' "$map")" \
		"$(stat -c '%n:%s:%Y' "$exe")" \
		"$(stat -c '%n:%s:%Y' "$base_model_args")" \
		"$(stat -c '%n:%s:%Y' "$base_model_cp")" \
		"$(stat -c '%n:%s:%Y' "$smoother_model_args")" \
		"$(stat -c '%n:%s:%Y' "$smoother_model_cp")" \
		"$merge" "${AS3_STRIDE:-512}" "${AS3_ANC:-0}" "$target_chunk_size"
}

inference_record() {
	local target=$1 ref=$2 map=$3
	printf 'schema\t1\ntarget\t%s\nreference\t%s\nmap\t%s\nexecutable\t%s\nbase_model_args\t%s\nbase_model\t%s\nsmoother_model_args\t%s\nsmoother_model\t%s\nstride\t%s\nanc\t%s\ntarget_chunk_size\t%s\n' \
		"$(stat -c '%n:%s:%Y' "$target")" \
		"$(stat -c '%n:%s:%Y' "$ref")" \
		"$(stat -c '%n:%s:%Y' "$map")" \
		"$(stat -c '%n:%s:%Y' "$exe")" \
		"$(stat -c '%n:%s:%Y' "$base_model_args")" \
		"$(stat -c '%n:%s:%Y' "$base_model_cp")" \
		"$(stat -c '%n:%s:%Y' "$smoother_model_args")" \
		"$(stat -c '%n:%s:%Y' "$smoother_model_cp")" \
		"${AS3_STRIDE:-512}" "${AS3_ANC:-0}" "$target_chunk_size"
}

prediction_output_present() {
	local dir=$1
	[[ -e $dir/introgression_prediction.txt ]] || return 1
	[[ -e $dir/introgression.bed ]]
}

resumable_output_present() {
	local dir=$1 target=$2 ref=$3 map=$4
	[[ -s $dir/.inference.complete ]] || return 1
	[[ -s $dir/.inference.provenance ]] || return 1
	[[ -e $dir/introgression_prediction.txt ]] || return 1
	[[ -e $dir/introgression.raw.bed ]] || return 1
	[[ -e $dir/introgression.raw.snps.gz ]] || return 1
	cmp -s <(inference_record "$target" "$ref" "$map") "$dir/.inference.provenance"
}

write_task_owner() {
	local dir=$1 key=$2 start_ticks
	start_ticks=$(awk '{print $22}' "/proc/$BASHPID/stat")
	printf 'schema\t1\npid\t%s\nstart_ticks\t%s\nkey\t%s\n' \
		"$BASHPID" "$start_ticks" "$key" >"$dir/.owner"
}

task_dir_owner_live() {
	local dir=$1 pid start_ticks actual cmdline
	if [[ -s $dir/.owner ]]; then
		pid=$(awk -F '\t' '$1=="pid"{print $2}' "$dir/.owner")
		start_ticks=$(awk -F '\t' '$1=="start_ticks"{print $2}' "$dir/.owner")
		[[ $pid =~ ^[0-9]+$ && $start_ticks =~ ^[0-9]+$ ]] || return 1
		[[ -r /proc/$pid/stat ]] || return 1
		actual=$(awk '{print $22}' "/proc/$pid/stat" 2>/dev/null || true)
		[[ $actual == "$start_ticks" ]]
		return
	fi

	# Legacy task directories encoded only a PID. Require both a live PID and
	# an AS3-looking command line so an unrelated process that reused the PID
	# cannot block the task forever.
	pid=${dir##*.}
	[[ $pid =~ ^[0-9]+$ && -r /proc/$pid/cmdline ]] || return 1
	cmdline=$(tr '\0' ' ' <"/proc/$pid/cmdline" 2>/dev/null || true)
	[[ $cmdline == *'/as3.sh'* || $cmdline == *'ArchaicSeeker3.1-mamba'* ]]
}

resume_dir=
prune_stale_task_dirs() {
	local panel=$1 unit=$2 target=$3 ref=$4 map=$5 parent d selected_resume=
	parent=$out/results/$panel
	for d in "$parent"/.chr"${unit}".run.* "$parent"/.chr"${unit}".old.*; do
		[[ -e $d ]] || continue
		if task_dir_owner_live "$d"; then
			echo "ERROR AS3 temporary directory belongs to a live process: $d" >&2
			echo "ERROR refusing a concurrent run for panel=$panel unit=$unit" >&2
			return 1
		fi
		if [[ $d == *.run.* ]] && resumable_output_present "$d" "$target" "$ref" "$map"; then
			if [[ -z $selected_resume || $d -nt $selected_resume ]]; then
				[[ -z $selected_resume ]] || {
					echo "[AS3 ${panel}.chr${unit}] REMOVE superseded inference checkpoint: $selected_resume"
					rm -rf -- "$selected_resume"
				}
				selected_resume=$d
			else
				echo "[AS3 ${panel}.chr${unit}] REMOVE superseded inference checkpoint: $d"
				rm -rf -- "$d"
			fi
			continue
		fi
		echo "[AS3 ${panel}.chr${unit}] REMOVE stale temporary directory: $d"
		rm -rf -- "$d"
	done
	resume_dir=$selected_resume
}

run_task() (
	local panel=$1 unit=$2 target=$3 ref=$4 map=$5 sample_overlap=$6 gpu=$7 key o run_dir old marker n rc log_file resume=0 preserve=0
	key=${panel}.chr${unit}
	o=$out/results/$panel/chr$unit
	mkdir -p "$out/results/$panel"
	resume_dir=
	prune_stale_task_dirs "$panel" "$unit" "$target" "$ref" "$map" || return 1
	marker=$o/.complete
	if [[ $replace == 0 && -s $marker ]] && prediction_output_present "$o"; then
		echo "[AS3 $key] SKIP"
		return 0
	fi
	run_dir=$out/results/$panel/.chr${unit}.run.$BASHPID
	old=$out/results/$panel/.chr${unit}.old.$BASHPID
	rm -rf "$old"
	if [[ -n $resume_dir ]]; then
		if [[ $resume_dir != "$run_dir" ]]; then
			rm -rf "$run_dir"
			mv -- "$resume_dir" "$run_dir"
		fi
		resume=1
		echo "[AS3 $key] RESUME checkpoint=$run_dir stage=postprocess"
	else
		rm -rf "$run_dir"
		mkdir -p "$run_dir"
		inference_record "$target" "$ref" "$map" >"$run_dir/.inference.provenance"
	fi
	write_task_owner "$run_dir" "$key"
	cleanup_task() {
		local rc=$?
		trap - EXIT HUP INT TERM
		set +e
		if ((preserve == 1)) || resumable_output_present "$run_dir" "$target" "$ref" "$map"; then
			rm -f -- "$run_dir/.owner"
			echo "[AS3 $key] PRESERVE inference checkpoint=$run_dir stage=postprocess exit=$rc" >&2
		else
			rm -rf -- "$run_dir"
		fi
		if [[ -e $old ]]; then
			if [[ -e $o ]]; then rm -rf -- "$old"; else mv -- "$old" "$o"; fi
		fi
		exit "$rc"
	}
	trap cleanup_task EXIT
	trap 'exit 130' HUP INT TERM
	local -a extra=(--stride "${AS3_STRIDE:-512}" --anc "${AS3_ANC:-0}" --target-chunk-size "$target_chunk_size")
	n=$(bcftools query -l "$target" | wc -l)
	# CUDA_VISIBLE_DEVICES selects a physical GPU ID.  Inside the masked child
	# process, that selected card is exposed to PyTorch as logical cuda:0.
	echo "[AS3 $key] RUN cuda_required=true physical_gpu_id=$gpu CUDA_VISIBLE_DEVICES=$gpu process_device=cuda:0 samples=$n target_chunk_size=$target_chunk_size cgroup_memory_high=$memory_high cgroup_memory_max=$memory_max cgroup_swap_max=$memory_swap_max Ref1028_overlap=$sample_overlap"
	mkdir -p "$runtime_work"
	local -a resource_scope=(systemd-run --user --scope --quiet --collect
		-p "MemoryHigh=$memory_high" -p "MemoryMax=$memory_max" -p "MemorySwapMax=$memory_swap_max")
	log_file=$log_root/${key}.log
	# The tool writes run.log inside its output folder; keep that log disposable.
	if [[ -f $run_dir/run.log && ! -L $run_dir/run.log ]]; then
		mv -- "$run_dir/run.log" "$log_root/${key}.tool.log"
	fi
	ln -sfn -- "$log_root/${key}.tool.log" "$run_dir/run.log"
	set +e
	if ((resume == 1)); then
		(cd "$runtime_work" && CUDA_VISIBLE_DEVICES=$gpu "${resource_scope[@]}" "$as3_python" "$exe" -t "$target" -r "$ref" -m "$map" \
			--base-model-cp "$base_model_cp" --base-model-args "$base_model_args" \
			--smoother-model-cp "$smoother_model_cp" --smoother-model-args "$smoother_model_args" \
			--merge "$merge" -o "$run_dir" "${extra[@]}" --postprocess-only) >>"$log_file" 2>&1
	else
		(cd "$runtime_work" && CUDA_VISIBLE_DEVICES=$gpu "${resource_scope[@]}" "$as3_python" "$exe" -t "$target" -r "$ref" -m "$map" \
			--base-model-cp "$base_model_cp" --base-model-args "$base_model_args" \
			--smoother-model-cp "$smoother_model_cp" --smoother-model-args "$smoother_model_args" \
			--merge "$merge" -o "$run_dir" "${extra[@]}") >"$log_file" 2>&1
	fi
	rc=$?
	set -e
	if ((rc != 0)); then
		if ((rc == 137)) || grep -Eqi 'CUDA out of memory|OutOfMemoryError|oom-kill|Killed process' "$log_file"; then
			echo "ERROR $key AS3 ran out of memory (exit=$rc, MemoryMax=$memory_max); see $log_file" >&2
		else
			echo "ERROR $key AS3 failed (exit=$rc); see $log_file" >&2
		fi
		if resumable_output_present "$run_dir" "$target" "$ref" "$map"; then
			preserve=1
			echo "ERROR $key raw inference is complete and was retained; the next identical run will resume post-processing" >&2
		fi
		return 1
	fi
	prediction_output_present "$run_dir" || {
		echo "ERROR $key missing prediction outputs" >&2
		rm -rf "$run_dir"
		return 1
	}
	completion_record "$target" "$ref" "$map" >"$run_dir/.complete"
	rm -f -- "$run_dir/.owner"
	[[ -e $o ]] && mv "$o" "$old"
	if ! mv "$run_dir" "$o"; then
		[[ -e $old ]] && mv "$old" "$o"
		return 1
	fi
	rm -rf "$old"
	echo "[AS3 $key] DONE output=$o"
)

status=0
for gpu_index in "${!gpu_array[@]}"; do
	(
		worker_status=0
		for i in "${!tasks[@]}"; do
			((i % ${#gpu_array[@]} == gpu_index)) || continue
			IFS=$'\t' read -r panel unit target ref map sample_overlap <<<"${tasks[$i]}"
			run_task "$panel" "$unit" "$target" "$ref" "$map" "$sample_overlap" "${gpu_array[$gpu_index]}" || worker_status=1
		done
		exit "$worker_status"
	) &
done
for pid in $(jobs -p); do wait "$pid" || status=1; done
((status == 0)) || {
	echo "ERROR one or more AS3 tasks failed; see $log_root" >&2
	exit 1
}
cp "$manifest" "$out/input.manifest.tsv"
echo "ALL DONE: AS3 output=$out tasks=${#tasks[@]}"

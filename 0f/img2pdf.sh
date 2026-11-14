#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Default  Convert input images to PDF

Main options:
    -o FILE  Output PDF [merged.pdf]
    -h, --help  Show options and examples below

Examples:
    ./img2pdf.sh *.png -o figures.pdf
    ls -1v *.png | ./img2pdf.sh -o figures.pdf
PUBLIC_ENTRY_HELP
    exit 0
fi
# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
set -euo pipefail


# 🚩 Command-line help
usage() {
	bash "${BASH_SOURCE[0]}" --help
}

OUT="merged.pdf"
declare -a ARG_FILES=()


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in
		-o | --output)
			if (($# < 2)); then
				echo "ERROR: $1 requires a PDF filename." >&2
				usage >&2
				exit 2
			fi
			OUT="$2"
			shift 2
			;;
		-h | --help)
			usage
			exit 0
			;;
		--)
			shift
			while (($#)); do
				ARG_FILES+=("$1")
				shift
			done
			;;
		-*)
			echo "ERROR: Unknown option: $1" >&2
			usage >&2
			exit 2
			;;
		*)
			ARG_FILES+=("$1")
			shift
			;;
	esac
done

# Dependencies
if command -v magick >/dev/null 2>&1; then
	IM=(magick)
elif command -v convert >/dev/null 2>&1; then
	IM=(convert)
else
	echo "ERROR: ImageMagick not found (need 'magick' or 'convert')." >&2
	exit 127
fi

if ! command -v img2pdf >/dev/null 2>&1; then
	echo "ERROR: img2pdf not found. Install it, e.g.:" >&2
	echo "  sudo apt install img2pdf" >&2
	exit 127
fi

# Collect files.
# If stdin is piped/redirected, read one filename per line first.
# Then append any filenames provided as positional arguments.
declare -a FILES=()

if [[ ! -t 0 ]]; then
	while IFS= read -r line || [[ -n "$line" ]]; do
		[[ -z "$line" ]] && continue
		FILES+=("$line")
	done
fi

if ((${#ARG_FILES[@]})); then
	FILES+=("${ARG_FILES[@]}")
fi

if ((${#FILES[@]} == 0)); then
	echo "ERROR: No input images were provided." >&2
	echo >&2
	usage >&2
	exit 2
fi

# Validate before doing any work.
for f in "${FILES[@]}"; do
	if [[ ! -f "$f" ]]; then
		echo "ERROR: File not found: $f" >&2
		exit 1
	fi
done

TMPDIR_WORK="$(mktemp -d "${TMPDIR:-$DIR_RUNTIME_TMP}/img2pdf.XXXXXX")"
cleanup() {
	rm -rf "$TMPDIR_WORK"
}
trap cleanup EXIT

declare -a PAGES=()
i=0
total=${#FILES[@]}

for f in "${FILES[@]}"; do
	i=$((i + 1))
	printf 'Processing %d/%d: %s\n' "$i" "$total" "$f" >&2

	# Get source width so header text size can be adapted for long names.
	if command -v identify >/dev/null 2>&1; then
		width="$(identify -format '%w' "$f" 2>/dev/null || echo 1200)"
	else
		width="$("${IM[@]}" identify -format '%w' "$f" 2>/dev/null || echo 1200)"
	fi

	[[ "$width" =~ ^[0-9]+$ ]] || width=1200

	# Header height and font size scale with image width, but stay sensible.
	pointsize=$((width / 10))
	((pointsize < 36)) && pointsize=36
	((pointsize > 90)) && pointsize=90

	header=$((pointsize * 2 + 40))
	((header < 120)) && header=120

	page="$TMPDIR_WORK/page_$(printf '%06d' "$i").png"

	# Add a white strip above the original image and write the original
	# filename/path there.  The original image is otherwise unchanged.
	"${IM[@]}" "$f" \
		-background white \
		-alpha remove \
		-alpha off \
		-gravity north \
		-splice "0x${header}" \
		-fill black \
		-pointsize "$pointsize" \
		-gravity north \
		-annotate "+0+$((pointsize / 2))" "$f" \
		-alpha off \
		"PNG24:$page"

	PAGES+=("$page")
done

# img2pdf preserves argument order, hence the PDF follows stdin/argument order.
img2pdf "${PAGES[@]}" -o "$OUT"

printf 'Created: %s (%d pages)\n' "$OUT" "$total" >&2

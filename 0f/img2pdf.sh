#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  # Recommended: preserve exactly the order produced by ls -v
  ls -1v *.png *.jpeg abc/*.png | img2pdf.sh [-o OUTPUT.pdf]

  # Or pass files directly (order is whatever the shell supplies)
  img2pdf.sh [-o OUTPUT.pdf] IMAGE [IMAGE ...]

Examples:
  ls -1v *.png | img2pdf.sh -o figures.pdf
  ls -1v *.png *.jpeg abc/*.png | img2pdf.sh -o all.pdf
  img2pdf.sh -o figures.pdf Fig1.png Fig2.png abc/Fig3.png

Default output:
  merged.pdf

Notes:
  - In pipe mode, input order is preserved exactly. The script does NOT sort.
  - The page header uses the filename/path exactly as supplied, e.g. abc/XXX.png.
  - One input image becomes one PDF page.
EOF
}

OUT="merged.pdf"
declare -a ARG_FILES=()

while (($#)); do
    case "$1" in
        -o|--output)
            if (($# < 2)); then
                echo "ERROR: $1 requires a PDF filename." >&2
                usage >&2
                exit 2
            fi
            OUT="$2"
            shift 2
            ;;
        -h|--help)
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

TMPDIR_WORK="$(mktemp -d "${TMPDIR:-/tmp}/img2pdf.XXXXXX")"
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
    pointsize=$(( width / 10 ))
    (( pointsize < 36 )) && pointsize=36
    (( pointsize > 90 )) && pointsize=90

    header=$(( pointsize * 2 + 40 ))
    (( header < 120 )) && header=120

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

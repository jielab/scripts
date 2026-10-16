#!/usr/bin/env python3
"""Summarize PASS SNP and indel sites separately by chromosome with bcftools."""

import argparse
import csv
import io
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from vcf_inputs import chrom_key, validate_inputs

CATEGORIES = (
    "singleton", "non_singleton_biallelic", "tri_allelic", "four_plus_allelic"
)
VARIANT_TYPES = ("SNP", "INDEL")
CHROM = r"(?:[1-9][0-9]*|X|Y|MT|M)"

# Paths are positional shell arguments, never interpolated into shell source.
# Discard unrelated INFO fields before joining (e.g. allele-specific annotations).
# Split first so mixed records lose only failing/zero-count/non-SNP/indel ALTs.
# Join with +both to keep SNP and indel groups distinct at the same coordinate.
PIPELINE = r"""
set -euo pipefail
bcftools view --no-update -G -f PASS -Ou "$1" |
bcftools annotate -x 'ID,QUAL,^INFO/AC' -Ou |
bcftools norm -m -any -Ou |
bcftools view -v snps,indels -i 'INFO/AC>0' -Ou |
bcftools norm -m +both -Ou |
bcftools query -f '%CHROM\t%POS\t%TYPE\t%ALT\t%INFO/AC\n' |
awk -F '\t' -v expected="$2" '
function fail(message) {
    print "check_ac: " message > "/dev/stderr"
    bad=1
    exit 2
}
{
    chrom=$1
    sub(/^chr/, "", chrom)
    if (chrom != expected) fail("filename/CHROM mismatch at " $1 ":" $2)
    if (NR > 1 && $2 < previous) fail("VCF must be sorted at " $1 ":" $2)
    if ($2 != previous) delete seen
    previous=$2
    type=$3
    if (type != "SNP" && type != "INDEL") fail("unexpected variant type at " $1 ":" $2 ": " type)
    if (seen[type]++) fail("unmerged " type " positions at " $1 ":" $2)
    nalt=split($4, alt, ",")
    nac=split($5, ac, ",")
    if (nalt != nac) fail("AC/ALT lengths differ at " $1 ":" $2)
    for (i=1; i<=nac; i++)
        if (ac[i] !~ /^[0-9]+$/ || ac[i]+0 < 1)
            fail("invalid AC at " $1 ":" $2)
    if      (nalt == 1 && ac[1]+0 == 1) k=1
    else if (nalt == 1)                k=2
    else if (nalt == 2)                k=3
    else                              k=4
    count[type,k]++
}
END {
    if (!bad) {
        printf "SNP\t%.0f\t%.0f\t%.0f\t%.0f\n", count["SNP",1], count["SNP",2], count["SNP",3], count["SNP",4]
        printf "INDEL\t%.0f\t%.0f\t%.0f\t%.0f\n", count["INDEL",1], count["INDEL",2], count["INDEL",3], count["INDEL",4]
    }
}'
"""


def discover(directory, prefix, selected=None):
    folder = Path(directory).expanduser().resolve()
    if not folder.is_dir():
        raise ValueError(f"VCF directory does not exist: {folder}")
    if not prefix or "/" in prefix or "\\" in prefix:
        raise ValueError("--prefix must be a nonempty filename prefix, without a directory")
    pattern = re.compile(
        re.escape(prefix) + rf"(?P<chrom>{CHROM})(?:[._].*)?\.vcf(?:\.bgz|\.gz)?$"
    )
    requested = None
    if selected and selected != "all":
        requested = {part.strip().removeprefix("chr") for part in selected.split(",")}
        if any(not re.fullmatch(CHROM, chrom) for chrom in requested):
            raise ValueError("--chr must be all or a comma list such as 18,19,21,X,Y")
    matches = {}
    for path in folder.iterdir():
        match = pattern.fullmatch(path.name)
        if not match or not path.is_file():
            continue
        chrom = match["chrom"]
        if requested is not None and chrom not in requested:
            continue
        if chrom in matches:
            raise ValueError(f"Multiple VCFs match chr{chrom}; narrow --prefix: {matches[chrom]}, {path}")
        matches[chrom] = path
    if requested is not None and requested - matches.keys():
        absent = sorted(requested - matches.keys(), key=chrom_key)
        raise ValueError("Missing requested chromosomes: " + ",".join(absent))
    if not matches:
        raise ValueError(f"No chromosome VCFs match prefix {prefix!r} in {folder}")
    return [(chrom, matches[chrom]) for chrom in sorted(matches, key=chrom_key)]


def check_ac_header(path):
    result = subprocess.run(
        ["bcftools", "view", "--no-version", "-h", str(path)],
        text=True, capture_output=True,
    )
    if result.returncode:
        raise ValueError(f"Cannot read VCF header: {path}\n{result.stderr.strip()}")
    if not re.search(r"^##INFO=<ID=AC,Number=A,Type=Integer,", result.stdout, re.M):
        raise ValueError(f"{path}: INFO/AC must be Number=A, Type=Integer for safe allele merging")


def summarize(path, chrom):
    result = subprocess.run(
        ["bash", "-c", PIPELINE, "check_ac", str(path), chrom],
        text=True, capture_output=True,
    )
    if result.returncode:
        raise ValueError(f"Statistics failed for {path}\n{result.stderr.strip()}")
    rows = [line.split("\t") for line in result.stdout.strip().splitlines()]
    if ([row[0] for row in rows] != list(VARIANT_TYPES)
            or any(len(row) != len(CATEGORIES) + 1 for row in rows)):
        raise ValueError(f"Unexpected statistics output for {path}: {result.stdout!r}")
    return {row[0]: [int(value) for value in row[1:]] for row in rows}


def output_path(value, files):
    """Validate the optional destination before a potentially long VCF scan."""
    if value is None:
        return None
    path = Path(value).expanduser().resolve()
    if not path.parent.is_dir():
        raise ValueError(f"Output directory does not exist: {path.parent}")
    if path.exists() and not path.is_file():
        raise ValueError(f"--out must name a TSV file, not a directory: {path}")
    protected = {
        Path(str(vcf) + suffix).resolve()
        for _, vcf in files
        for suffix in ("", ".tbi", ".csi", ".aria2", ".tbi.aria2", ".csi.aria2")
    }
    if path in protected:
        raise ValueError(f"--out must not overwrite an input VCF or index: {path}")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="./allele.sh check_ac",
        description="按染色体分别统计 PASS、AC>0 的 SNP 和 indel；显示各类别数量及比例。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例：\n  ./allele.sh check_ac --dir /mnt/f/annot/gnomad "
               "--prefix gnomad.genomes.v4.1.sites.chr --chr 18,19,21 \\\n"
               "    --out /mnt/f/annot/gnomad/check_ac.chr18_19_21.tsv\n\n"
               "每条染色体完成后显示 SNP、INDEL 各自的分类计数和比例。\n"
               "全部成功后保存完整 TSV；省略 --out 时只输出到终端。",
    )
    parser.add_argument("--dir", required=True, help="VCF 文件所在目录（不递归）")
    parser.add_argument("--prefix", required=True, help="染色体编号之前的文件名前缀")
    parser.add_argument("--chr", help="如 18,19,21；省略或 all 时选择全部匹配染色体，含 X/Y")
    parser.add_argument("--out", help="另存完整 TSV 的文件路径（父目录须已存在）；仍输出到终端")
    args = parser.parse_args(argv)
    for command in ("bcftools", "bash", "awk"):
        if not shutil.which(command):
            parser.error(f"missing command: {command}")
    try:
        files = discover(args.dir, args.prefix, args.chr)
        output = output_path(args.out, files)
        # Validate the whole selection before scanning any chromosome; no partial
        # download is silently skipped. No index is needed for generic full VCFs.
        validate_inputs([path for _, path in files])
        for _, path in files:
            check_ac_header(path)
        print(f"[check_ac] {len(files)} chromosome(s); PASS SNP and INDEL counted separately; ALT AC>0", file=sys.stderr)
        destination = f"stdout + {output} (saved after all chromosomes succeed)" if output else "stdout only (use --out FILE to save)"
        print(f"[check_ac] output: {destination}", file=sys.stderr, flush=True)
        summaries = []
        for chrom, path in files:
            print(f"[check_ac] chr{chrom}: {path}", file=sys.stderr, flush=True)
            counts_by_type = summarize(path, chrom)
            for variant_type, counts in counts_by_type.items():
                total = sum(counts)
                print(f"[check_ac] chr{chrom} {variant_type}: {total:,} merged sites", file=sys.stderr, flush=True)
                for category, count in zip(CATEGORIES, counts):
                    pct = f"{100 * count / total:.4f}" if total else "NA"
                    summaries.append((f"chr{chrom}", variant_type, category, count, total, pct))
                    percent_label = f"{pct}%" if total else "NA"
                    print(f"[check_ac]   {category}: {count:,} ({percent_label})", file=sys.stderr, flush=True)
        # Publish only the small final table, after every input succeeds.
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(("CHR", "TYPE", "CATEGORY", "COUNT", "TOTAL_SITES", "PERCENT"))
        writer.writerows(summaries)
        table = buffer.getvalue()
        if output:
            # Replace only a complete table; keep any previous file on scan/write failure.
            with tempfile.TemporaryDirectory(prefix=".check_ac-", dir=output.parent) as staging:
                temporary = Path(staging) / "summary.tsv"
                temporary.write_text(table, encoding="utf-8")
                temporary.replace(output)
            print(f"[check_ac] saved: {output}", file=sys.stderr, flush=True)
        sys.stdout.write(table)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

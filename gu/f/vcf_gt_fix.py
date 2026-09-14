#!/usr/bin/env python3
"""Stream a VCF while validating/fixing GT ploidy for X-chromosome units.

This does not invent phase. With --require-phased, unphased heterozygotes fail.
Haploid calls can be duplicated as homozygous diploid calls when a downstream
program requires two alleles per sample.
"""
from __future__ import annotations

import argparse
import sys


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", help="replace CHROM with this value (e.g. 23)")
    p.add_argument("--duplicate-haploid", choices=["slash", "pipe"])
    p.add_argument("--require-diploid", action="store_true")
    p.add_argument("--require-phased", action="store_true")
    args = p.parse_args()

    n_records = n_calls = n_haploid = n_unphased_het = 0
    for line_no, raw in enumerate(sys.stdin, start=1):
        if raw.startswith("#"):
            sys.stdout.write(raw)
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 10:
            raise SystemExit(f"Malformed VCF line {line_no}: expected >=10 columns")
        n_records += 1
        if args.chrom:
            fields[0] = args.chrom
        fmt = fields[8].split(":")
        try:
            gt_idx = fmt.index("GT")
        except ValueError:
            raise SystemExit(f"VCF line {line_no} has no GT field")
        for i in range(9, len(fields)):
            vals = fields[i].split(":")
            if gt_idx >= len(vals):
                raise SystemExit(f"VCF line {line_no}, sample column {i+1}: missing GT value")
            gt = vals[gt_idx]
            n_calls += 1
            if gt in {".", ""}:
                if args.duplicate_haploid:
                    sep = "|" if args.duplicate_haploid == "pipe" else "/"
                    vals[gt_idx] = f".{sep}."
                elif args.require_diploid:
                    raise SystemExit(f"Haploid/missing-width GT {gt!r} at VCF line {line_no}")
            elif "/" not in gt and "|" not in gt:
                n_haploid += 1
                if args.duplicate_haploid:
                    sep = "|" if args.duplicate_haploid == "pipe" else "/"
                    vals[gt_idx] = f"{gt}{sep}{gt}"
                elif args.require_diploid:
                    raise SystemExit(f"Haploid GT {gt!r} at VCF line {line_no}; duplication was not enabled")
            else:
                sep = "|" if "|" in gt else "/"
                alleles = gt.split(sep)
                if len(alleles) != 2:
                    raise SystemExit(f"Non-diploid GT {gt!r} at VCF line {line_no}")
                if args.require_phased and sep == "/" and len(set(alleles) - {"."}) > 1:
                    n_unphased_het += 1
                    if n_unphased_het <= 5:
                        print(f"UNPHASED_HET line={line_no} gt={gt}", file=sys.stderr)
            fields[i] = ":".join(vals)
        sys.stdout.write("\t".join(fields) + "\n")
    if args.require_phased and n_unphased_het:
        raise SystemExit(f"Found {n_unphased_het} unphased heterozygous GT calls; AS3 requires phased target haplotypes")
    print(f"records={n_records} calls={n_calls} haploid_duplicated_or_seen={n_haploid}", file=sys.stderr)


if __name__ == "__main__":
    main()

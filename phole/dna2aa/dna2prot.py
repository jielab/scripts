# pip install alphagenome biopython
# 启明 module load python/anaconda3/2020.7; source activate; conda activate alphagenome

#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


COLS = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--chr", required=True)          # 无 chr 前缀，例如 9
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    ap.add_argument("--hap-fa", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-fa", required=True)
    return ap.parse_args()


def parse_attr(attr: str) -> dict:
    return {k: v for k, v in ATTR_RE.findall(attr)}


def load_gtf(gtf_path: str) -> pd.DataFrame:
    return pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=COLS,
        compression="infer",
        low_memory=False
    )


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def wrap60(s: str):
    for i in range(0, len(s), 60):
        yield s[i:i+60]


def main():
    args = parse_args()

    gene = args.gene
    chrom = args.chr
    gstart = args.start
    gend = args.end
    out_fa = Path(args.out_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    rec = next(SeqIO.parse(args.hap_fa, "fasta"))
    hap = str(rec.seq).upper().replace("U", "T")
    expected_len = gend - gstart + 1

    print(f"gene={gene}")
    print(f"region={chrom}:{gstart}-{gend}")
    print(f"hap_fa={args.hap_fa}")
    print(f"hap_len={len(hap)}")
    if len(hap) != expected_len:
        print(f"[warn] haplotype length {len(hap)} != expected gene span {expected_len}")

    gtf = load_gtf(args.gtf)

    # 兼容 GTF 里 seqname 可能是 9 或 chr9
    valid_chr = {chrom, f"chr{chrom}"}
    gtf = gtf[(gtf["feature"] == "CDS") & (gtf["seqname"].isin(valid_chr))].copy()

    attrs = gtf["attributes"].apply(parse_attr)
    gtf["gene_name"] = attrs.apply(lambda d: d.get("gene_name"))
    gtf["transcript_id"] = attrs.apply(lambda d: d.get("transcript_id"))
    gtf["transcript_biotype"] = attrs.apply(lambda d: d.get("transcript_biotype") or d.get("transcript_type"))
    gtf["protein_id"] = attrs.apply(lambda d: d.get("protein_id"))
    gtf["frame"] = pd.to_numeric(gtf["frame"], errors="coerce")

    cds = gtf[gtf["gene_name"] == gene].copy()
    if cds.empty:
        raise ValueError(f"No CDS entries found for gene={gene}")

    cds = cds[(cds["end"] >= gstart) & (cds["start"] <= gend)].copy()
    if cds.empty:
        raise ValueError(f"No CDS overlaps found in region {chrom}:{gstart}-{gend}")

    def extract_hap(seg_start: int, seg_end: int) -> str:
        # gene region is 1-based inclusive
        i0 = max(seg_start, gstart) - gstart
        i1 = min(seg_end, gend) - gstart + 1
        if i0 < 0 or i1 > len(hap) or i0 >= i1:
            return ""
        return hap[i0:i1]

    def clip_by_frame(piece: str, frame, strand: str) -> str:
        if pd.isna(frame) or int(frame) == 0:
            return piece
        frame = int(frame)
        if strand == "+":
            return piece[frame:]
        return piece[:-frame] if frame <= len(piece) else ""

    aa_records = []
    for tid, sub in cds.groupby("transcript_id", sort=False):
        strand = sub["strand"].iloc[0]
        biotype = sub["transcript_biotype"].iloc[0]
        protein_id = sub["protein_id"].dropna().iloc[0] if sub["protein_id"].notna().any() else ""

        sub = sub.sort_values("start", ascending=(strand == "+")).copy()
        pieces = []

        for _, row in sub.iterrows():
            piece = extract_hap(int(row["start"]), int(row["end"]))
            piece = clip_by_frame(piece, row["frame"], strand)
            pieces.append(piece)

        cds_seq = "".join(pieces)

        if strand == "-":
            cds_seq = revcomp(cds_seq)

        if not cds_seq:
            continue

        if len(cds_seq) % 3 != 0:
            drop = len(cds_seq) % 3
            print(f"[warn] {tid}: CDS length not multiple of 3; trimming {drop} base(s).")
            cds_seq = cds_seq[:len(cds_seq) - drop]

        aa_seq = str(Seq(cds_seq).translate(to_stop=False))
        aa_records.append((tid, protein_id, strand, biotype, len(cds_seq), len(aa_seq), aa_seq))

    if not aa_records:
        raise ValueError(f"No protein sequences generated for gene={gene}")

    aa_df = pd.DataFrame(
        aa_records,
        columns=["transcript_id","protein_id","strand","biotype","cds_len","aa_len","aa_seq"]
    ).sort_values(["aa_len", "transcript_id"], ascending=[False, True])

    print(aa_df[["transcript_id","protein_id","strand","biotype","cds_len","aa_len"]].head(20).to_string(index=False))

    with out_fa.open("w") as fh:
        for _, row in aa_df.iterrows():
            header = f">{row.transcript_id}|protein_id={row.protein_id}|gene={gene}|aa_len={row.aa_len}|strand={row.strand}"
            fh.write(header + "\n")
            for chunk in wrap60(row["aa_seq"]):
                fh.write(chunk + "\n")

    print(f"Wrote: {out_fa}")


if __name__ == "__main__":
    main()
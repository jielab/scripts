#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from alphagenome.models import dna_client
from alphagenome.data import genome


# 这里并不是从 DNA 唯一推出“蛋白质”，
# 而是从 DNA + 注释模型，推出某些候选转录本对应的预测蛋白序列。
# AlphaGenome 在本流程中主要用于提供 RNA / splicing 层面的辅助证据，
# 而不是直接输出 protein FASTA 或“唯一正确 isoform”。

COLS = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--chr", required=True)          # 这里传 chr9 这种格式，供 AlphaGenome interval 使用
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    ap.add_argument("--hap-fa", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--api-key", default=os.getenv("ALPHA_GENOME_API_KEY", ""))
    ap.add_argument("--ontology-terms", default="")  # 逗号分隔，可留空
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

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not args.api_key:
        raise ValueError("AlphaGenome API key not found. Set ALPHA_GENOME_API_KEY or use --api-key.")

    model = dna_client.create(args.api_key)

    gene = args.gene
    chr_with_prefix = args.chr
    chr_no_prefix = chr_with_prefix.replace("chr", "")
    gstart = args.start
    gend = args.end

    rec = next(SeqIO.parse(args.hap_fa, "fasta"))
    hap = str(rec.seq).upper().replace("U", "T")
    hap_len = len(hap)
    expected_len = gend - gstart + 1

    print(f"gene={gene}")
    print(f"hap_fa={args.hap_fa}")
    print(f"interval gene span={chr_with_prefix}:{gstart}-{gend}")
    print(f"hap_len={hap_len}")
    if hap_len != expected_len:
        print(f"[warn] hap_len={hap_len} != expected_len={expected_len}")

    # 建立 1 Mb AlphaGenome 输入窗口
    L = dna_client.SEQUENCE_LENGTH_1MB
    center = (gstart + gend) // 2
    iv_start = center - (L // 2)
    iv_end = iv_start + L
    interval = genome.Interval(chr_with_prefix, iv_start, iv_end)

    left_pad_len = gstart - iv_start
    right_pad_len = L - left_pad_len - hap_len
    if left_pad_len < 0 or right_pad_len < 0:
        raise ValueError("Gene span does not fit into 1 Mb window. Check coordinates.")

    seq_for_model = "N" * left_pad_len + hap + "N" * right_pad_len
    assert len(seq_for_model) == L

    print(f"AlphaGenome interval={interval.chromosome}:{interval.start}-{interval.end}")
    print(f"Haplotype placed at {chr_with_prefix}:{gstart}-{gend}")

    requested_outputs = [
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.SPLICE_SITES,
        dna_client.OutputType.SPLICE_SITE_USAGE,
        dna_client.OutputType.SPLICE_JUNCTIONS,
    ]

    ontology_terms = None
    if args.ontology_terms.strip():
        ontology_terms = [x.strip() for x in args.ontology_terms.split(",") if x.strip()]

    out = model.predict_sequences(
        sequences=[seq_for_model],
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms,
        intervals=[interval],
        progress_bar=True,
    )[0]

    # 1) 写 metadata summaries
    if out.rna_seq is not None:
        out.rna_seq.metadata.to_csv(outdir / "rna_seq.metadata.tsv", sep="\t", index=False)
        print(f"RNA-seq tracks: {out.rna_seq.values.shape}")

    if out.splice_sites is not None:
        out.splice_sites.metadata.to_csv(outdir / "splice_sites.metadata.tsv", sep="\t", index=False)
        print(f"Splice-sites tracks: {out.splice_sites.values.shape}")

    if out.splice_site_usage is not None:
        out.splice_site_usage.metadata.to_csv(outdir / "splice_site_usage.metadata.tsv", sep="\t", index=False)
        print(f"Splice-site-usage tracks: {out.splice_site_usage.values.shape}")

    # JunctionData 往往不是普通 TrackData，下面尽量做稳妥输出
    if out.splice_junctions is not None:
        try:
            if hasattr(out.splice_junctions, "metadata") and out.splice_junctions.metadata is not None:
                out.splice_junctions.metadata.to_csv(outdir / "splice_junctions.metadata.tsv", sep="\t", index=False)
        except Exception as e:
            print(f"[warn] failed to save splice_junctions metadata: {e}")

    # 2) 同时再按固定 GTF 翻译 AA，便于和主流程对照
    gtf = load_gtf(args.gtf)
    valid_chr = {chr_no_prefix, chr_with_prefix}
    gtf = gtf[(gtf["feature"] == "CDS") & (gtf["seqname"].isin(valid_chr))].copy()

    attrs = gtf["attributes"].apply(parse_attr)
    gtf["gene_name"] = attrs.apply(lambda d: d.get("gene_name"))
    gtf["transcript_id"] = attrs.apply(lambda d: d.get("transcript_id"))
    gtf["transcript_biotype"] = attrs.apply(lambda d: d.get("transcript_biotype") or d.get("transcript_type"))
    gtf["protein_id"] = attrs.apply(lambda d: d.get("protein_id"))
    gtf["frame"] = pd.to_numeric(gtf["frame"], errors="coerce")

    cds = gtf[gtf["gene_name"] == gene].copy()
    cds = cds[(cds["end"] >= gstart) & (cds["start"] <= gend)].copy()
    if cds.empty:
        raise ValueError(f"No CDS entries found for {gene}")

    def extract_hap(seg_start: int, seg_end: int) -> str:
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

    aa_df = pd.DataFrame(
        aa_records,
        columns=["transcript_id","protein_id","strand","biotype","cds_len","aa_len","aa_seq"]
    ).sort_values(["aa_len", "transcript_id"], ascending=[False, True])

    aa_df.drop(columns=["aa_seq"]).to_csv(outdir / "predicted_proteins.summary.tsv", sep="\t", index=False)

    with open(outdir / "predicted_proteins.fa", "w") as fh:
        for _, row in aa_df.iterrows():
            header = f">{row.transcript_id}|protein_id={row.protein_id}|gene={gene}|aa_len={row.aa_len}|strand={row.strand}"
            fh.write(header + "\n")
            for chunk in wrap60(row["aa_seq"]):
                fh.write(chunk + "\n")

    print(f"Wrote: {outdir / 'predicted_proteins.fa'}")


if __name__ == "__main__":
    main()
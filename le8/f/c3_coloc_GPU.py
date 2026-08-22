#!/usr/bin/env python3
import csv
import gzip
import math
import os
import re
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import pyarrow.feather as feather
except Exception as e:
    raise SystemExit("This script needs pandas, numpy and pyarrow in the active Python environment. Install with: python -m pip install pandas numpy pyarrow\n" + str(e))


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def norm(x):
    return re.sub(r"[^A-Z0-9]+", "", str(x).upper())


def pick(cols, names):
    mapping = {norm(c): c for c in cols}
    for n in names:
        if norm(n) in mapping:
            return mapping[norm(n)]
    return None


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def read_header(path):
    with opener(path) as f:
        line = f.readline()
    delim = "\t" if line.count("\t") >= line.count(",") else ","
    cols = line.rstrip("\n\r").split(delim)
    return cols, delim


def load_region(path, chrom=None, start=None, end=None):
    cols, delim = read_header(path)
    c_snp = pick(cols, ["SNP", "RSID", "RS_NUMBER", "VARIANT_ID", "ID", "MARKERNAME"])
    c_chr = pick(cols, ["CHR", "CHROM", "CHROMOSOME", "chromosome"])
    c_pos = pick(cols, ["POS", "BP", "POSITION", "BASE_PAIR_LOCATION", "position"])
    c_ea = pick(cols, ["EA", "A1", "ALT", "ALLELE1", "EFFECT_ALLELE", "effect_allele"])
    c_nea = pick(cols, ["NEA", "A2", "REF", "ALLELE0", "OTHER_ALLELE", "reference_allele"])
    # GPU-coloc receives the full regional QTL/GWAS files, so prefer marginal effects.
    c_beta = pick(cols, ["BETA", "beta", "B", "EFFECT", "LOG_ODDS", "BJ", "bJ"])
    c_se = pick(cols, ["SE", "se", "SEBETA", "STDERR", "BJ_SE", "bJ_se"])
    c_p = pick(cols, ["P", "pval", "PVALUE", "P_VALUE", "PJ", "pJ"])
    required = [c_chr, c_pos, c_ea, c_nea, c_beta, c_se]
    if any(x is None for x in required):
        raise RuntimeError(f"Missing required columns in {path}. Need CHR POS EA NEA BETA SE.")
    usecols = [x for x in [c_snp, c_chr, c_pos, c_ea, c_nea, c_beta, c_se, c_p] if x]
    chunks = []
    for ch in pd.read_csv(path, sep=delim, usecols=lambda c: c in set(usecols), chunksize=500000, compression="infer"):
        ch = ch.rename(columns={c_chr: "CHR", c_pos: "POS", c_ea: "EA", c_nea: "NEA", c_beta: "BETA", c_se: "SE"})
        if c_snp:
            ch = ch.rename(columns={c_snp: "SNP"})
        else:
            ch["SNP"] = ""
        if c_p:
            ch = ch.rename(columns={c_p: "P"})
        else:
            ch["P"] = pd.NA
        ch["CHR"] = ch["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
        ch["POS"] = pd.to_numeric(ch["POS"], errors="coerce")
        ch["BETA"] = pd.to_numeric(ch["BETA"], errors="coerce")
        ch["SE"] = pd.to_numeric(ch["SE"], errors="coerce")
        ch["P"] = pd.to_numeric(ch["P"], errors="coerce")
        ch = ch.dropna(subset=["CHR", "POS", "BETA", "SE"])
        ch = ch[ch["SE"] > 0]
        if chrom is not None:
            ch = ch[(ch["CHR"].astype(str) == str(chrom)) & (ch["POS"] >= start) & (ch["POS"] <= end)]
        if len(ch):
            chunks.append(ch)
    if not chunks:
        return pd.DataFrame(columns=["SNP", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P"])
    d = pd.concat(chunks, ignore_index=True)
    d["EA"] = d["EA"].astype(str).str.upper()
    d["NEA"] = d["NEA"].astype(str).str.upper()
    d = d[(d["EA"] != "") & (d["NEA"] != "")]
    return d


def find_lead_and_region(path, region, window_kb):
    if region and re.match(r"^[A-Za-z0-9]+:\d+-\d+$", region):
        chrom, rest = region.split(":", 1)
        start, end = rest.split("-", 1)
        return chrom.replace("chr", ""), int(float(start)), int(float(end))
    d = load_region(path)
    if d.empty:
        raise RuntimeError("No usable QTL rows in " + path)
    if d["P"].notna().any():
        lead = d.loc[d["P"].idxmin()]
    else:
        z = (d["BETA"] / d["SE"]).abs()
        lead = d.loc[z.idxmax()]
    chrom = str(lead["CHR"]).replace("chr", "")
    pos = int(lead["POS"])
    w = window_kb * 1000
    return chrom, max(1, pos - w), pos + w


def to_signal(d, signal, typ):
    """Return the one-row, variant-as-column matrix required by gpu-coloc --format."""
    if d.empty:
        raise RuntimeError("Empty region for signal " + signal)
    d = d.copy()
    d["z"] = d["BETA"] / d["SE"]
    sd_prior = 0.15 if typ == "quant" else 0.20
    v = d["SE"] ** 2
    r = (sd_prior ** 2) / ((sd_prior ** 2) + v)
    d["lbf"] = 0.5 * (np.log1p(-r) + (r * (d["z"] ** 2)))
    # Log-BFs are invariant to effect-allele direction, so canonicalise the allele pair.
    # This prevents QTL/outcome files with swapped EA/NEA from creating different column names.
    def canonical_pair(row):
        pair = "_".join(sorted([str(row["NEA"]).upper(), str(row["EA"]).upper()]))
        comp = str.maketrans("ATCG", "TAGC")
        pair_comp = "_".join(sorted([str(row["NEA"]).upper().translate(comp), str(row["EA"]).upper().translate(comp)]))
        return min(pair, pair_comp)
    allele_pair = d.apply(canonical_pair, axis=1)
    d["variant"] = ("chr" + d["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
                    + "_" + d["POS"].astype(int).astype(str) + "_" + allele_pair)
    d = d[["variant", "lbf"]].replace([np.inf, -np.inf], np.nan).dropna().drop_duplicates("variant")
    if d.empty:
        raise RuntimeError("No variants after formatting for signal " + signal)
    lead_idx = d["lbf"].idxmax()
    lead = str(d.loc[lead_idx, "variant"])
    strength = float(d.loc[lead_idx, "lbf"])
    # gpu-coloc/format.py reads each signal file and treats its columns as variants,
    # taking df.iloc[0] as the signal's log-BF vector.
    wide = pd.DataFrame([d["lbf"].to_numpy()], columns=d["variant"].tolist())
    return wide, lead, strength


def prepare_signals(manifest, cad_gwas, outdir, window_kb, outcome_type="cc"):
    outdir = Path(outdir)
    qtl_dir = outdir / "qtl_signals"
    cad_dir = outdir / "cad_signals"
    qtl_dir.mkdir(parents=True, exist_ok=True)
    cad_dir.mkdir(parents=True, exist_ok=True)
    for old_file in list(qtl_dir.glob("*.feather")) + list(cad_dir.glob("*.feather")):
        old_file.unlink()

    with open(manifest, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    qtl_summary = []
    cad_summary = []
    status_rows = []
    for row in rows:
        omics = row.get("omics", "protein") or "protein"
        trait = row.get("trait", "")
        qtl_file = row.get("file", row.get("exposure", ""))
        typ = row.get("type", "quant") or "quant"
        region = row.get("region", "") or ""
        if not trait or not qtl_file or not os.path.exists(qtl_file):
            status_rows.append({"omics": omics, "trait": trait, "status": "not_run", "message": "missing qtl file"})
            continue
        try:
            chrom, start, end = find_lead_and_region(qtl_file, region, window_kb)
            qtl = load_region(qtl_file, chrom, start, end)
            cad = load_region(cad_gwas, chrom, start, end)
            qtl_signal = safe_name(f"QTL__{omics}__{trait}__chr{chrom}_{start}_{end}")
            cad_signal = safe_name(f"CAD__{omics}__{trait}__chr{chrom}_{start}_{end}")
            qtl_s, qtl_lead, qtl_strength = to_signal(qtl, qtl_signal, typ)
            cad_s, cad_lead, cad_strength = to_signal(cad, cad_signal, outcome_type)
            feather.write_feather(qtl_s, qtl_dir / f"{qtl_signal}.feather")
            feather.write_feather(cad_s, cad_dir / f"{cad_signal}.feather")
            qtl_summary.append({"signal": qtl_signal, "chromosome": chrom, "location_min": start, "location_max": end, "signal_strength": qtl_strength, "lead_variant": qtl_lead, "omics": omics, "trait": trait, "source_file": qtl_file})
            cad_summary.append({"signal": cad_signal, "chromosome": chrom, "location_min": start, "location_max": end, "signal_strength": cad_strength, "lead_variant": cad_lead, "omics": omics, "trait": trait, "source_file": cad_gwas})
            status_rows.append({"omics": omics, "trait": trait, "status": "ok", "message": f"chr{chrom}:{start}-{end}; qtl_n={qtl_s.shape[1]}; cad_n={cad_s.shape[1]}"})
        except Exception as e:
            status_rows.append({"omics": omics, "trait": trait, "status": "failed", "message": str(e)})

    pd.DataFrame(qtl_summary).to_csv(outdir / "qtl_summary.tsv", sep="\t", index=False)
    pd.DataFrame(cad_summary).to_csv(outdir / "cad_summary.tsv", sep="\t", index=False)
    pd.DataFrame(status_rows).to_csv(outdir / "signal_preparation_status.tsv", sep="\t", index=False)
    if len(qtl_summary) == 0 or len(cad_summary) == 0:
        raise SystemExit("No GPU-coloc signal files were created. See signal_preparation_status.tsv")


def main():
    if len(sys.argv) not in (5, 6):
        raise SystemExit("Usage: c3_coloc_GPU.py <manifest.tsv> <outcome_gwas> <outdir> <window_kb> [cc|quant]")
    outcome_type = sys.argv[5].lower() if len(sys.argv) == 6 else "cc"
    if outcome_type not in {"cc", "quant"}:
        raise SystemExit("outcome_type must be cc or quant")
    prepare_signals(sys.argv[1], sys.argv[2], sys.argv[3], int(float(sys.argv[4])), outcome_type)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 TRACE helper: sample map normalization and report building
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import json
from pathlib import Path


def read_panel(panel_file):
    import pandas as pd
    panel = pd.read_csv(panel_file, sep="\t", dtype=str)
    if len(panel.columns) < 3:
        raise SystemExit(f"sample panel must have at least three columns: {panel_file}")
    panel.columns = ["sample", "pop", "super_pop"] + list(panel.columns[3:])
    return panel


def normalize_sample_map(src, panel_file, out):
    import pandas as pd
    x = pd.read_csv(src, sep="\t", dtype=str)
    cols = {c.lower(): c for c in x.columns}
    node = cols.get("tree_node_id") or cols.get("node") or cols.get("haplotype_node")
    sample = cols.get("sample") or cols.get("sample_id") or cols.get("individual") or cols.get("id")
    hap = cols.get("haplotype") or cols.get("hap")
    if node is None:
        raise SystemExit("TRACE_SAMPLE_MAP needs a tree_node_id/node column")
    if sample is None:
        x["sample"] = x[node].map(lambda v: f"node_{v}")
        sample = "sample"
    if hap is None:
        x["haplotype"] = x.groupby(sample).cumcount() + 1
        hap = "haplotype"
    panel = read_panel(panel_file)
    y = x[[node, sample, hap]].copy()
    y.columns = ["tree_node_id", "sample", "haplotype"]
    y["tree_node_id"] = y["tree_node_id"].astype(str)
    y["sample"] = y["sample"].astype(str)
    y = y.merge(panel[["sample", "pop", "super_pop"]], on="sample", how="left")
    y.to_csv(out, sep="\t", index=False)


def metadata_dict(obj):
    metadata = getattr(obj, "metadata", None)
    if isinstance(metadata, dict):
        return metadata
    if isinstance(metadata, bytes):
        try:
            return json.loads(metadata.decode())
        except Exception:
            return {}
    return {}


def auto_sample_map(tree_file, panel_file, out):
    import pandas as pd
    import tskit
    import tszip

    tree_file = str(tree_file)
    ts = tszip.decompress(tree_file) if tree_file.endswith(".tsz") else tskit.load(tree_file)
    panel = read_panel(panel_file)
    panel_samples = set(panel["sample"].astype(str))
    rows = []
    seen = {}
    for node_id in ts.samples():
        node = ts.node(int(node_id))
        candidates = []
        node_metadata = metadata_dict(node)
        for key in ("sample", "sample_id", "individual", "individual_id", "name", "id"):
            if key in node_metadata:
                candidates.append(str(node_metadata[key]))
        if node.individual != tskit.NULL:
            individual_metadata = metadata_dict(ts.individual(node.individual))
            for key in ("sample", "sample_id", "individual", "individual_id", "name", "id"):
                if key in individual_metadata:
                    candidates.append(str(individual_metadata[key]))
        sample = next((value for value in candidates if value in panel_samples), None)
        if sample is None:
            sample = candidates[0] if candidates else f"node_{node_id}"
        seen[sample] = seen.get(sample, 0) + 1
        rows.append({"tree_node_id": str(int(node_id)), "sample": sample, "haplotype": seen[sample]})
    out_df = pd.DataFrame(rows).merge(panel[["sample", "pop", "super_pop"]], on="sample", how="left")
    out_df.to_csv(out, sep="\t", index=False)


def _first_existing(cols, candidates):
    low = {c.lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None


def normalize_trace_segments(hap):
    import pandas as pd
    if hap.empty:
        return hap
    cols = list(hap.columns)
    chrom_col = _first_existing(cols, ["chrom", "chromosome", "chr"])
    start_col = _first_existing(cols, ["start", "start_bp", "start(bp)", "left", "begin"])
    end_col = _first_existing(cols, ["end", "end_bp", "end(bp)", "right", "stop"])
    post_col = _first_existing(cols, ["mean_posterior", "posterior", "posterior_mean", "prob", "probability"])
    len_bp_col = _first_existing(cols, ["length_bp", "length(bp)", "length", "physical_length", "physical_length_bp"])
    len_cm_col = _first_existing(cols, ["length_cm", "length(cM)", "genetic_length", "genetic_length_cm"])

    if chrom_col is None or start_col is None or end_col is None:
        raise SystemExit("TRACE summary files must contain chromosome/chrom, start, and end columns")

    if chrom_col != "chrom":
        hap["chrom"] = hap[chrom_col].astype(str)
    else:
        hap["chrom"] = hap["chrom"].astype(str)
    hap["chrom"] = hap["chrom"].str.replace(r"^chr", "", regex=True)
    hap["chrom"] = hap["chrom"].str.replace(r"_part[0-9]+$", "", regex=True)
    hap.loc[hap["chrom"] == "23", "chrom"] = "X"
    hap["start"] = pd.to_numeric(hap[start_col], errors="coerce").astype("Int64")
    hap["end"] = pd.to_numeric(hap[end_col], errors="coerce").astype("Int64")
    if len_bp_col is not None:
        hap["length_bp"] = pd.to_numeric(hap[len_bp_col], errors="coerce")
    else:
        hap["length_bp"] = pd.to_numeric(hap["end"], errors="coerce") - pd.to_numeric(hap["start"], errors="coerce")
    if len_cm_col is not None:
        hap["length_cM"] = pd.to_numeric(hap[len_cm_col], errors="coerce")
    if post_col is not None:
        hap["mean_posterior"] = pd.to_numeric(hap[post_col], errors="coerce")
    hap = hap.dropna(subset=["chrom", "start", "end", "length_bp"])
    hap = hap[hap["end"].astype(float) > hap["start"].astype(float)].copy()
    hap["start"] = hap["start"].astype(int)
    hap["end"] = hap["end"].astype(int)
    hap["length_bp"] = hap["length_bp"].astype(float)
    hap["anc"] = "TRACE_Archaic"
    return hap


def chr_lengths_grch37(chroms):
    lengths = {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "23": 155270560,
        "Y": 59373566, "MT": 16569, "M": 16569,
    }
    return {c: lengths[c] for c in chroms if c in lengths}


def build_windows(chroms, window_size):
    import pandas as pd
    chroms = [str(c).replace("chr", "") for c in chroms]
    rank_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]
    lengths = chr_lengths_grch37(chroms)
    rows = []
    genome_pos = 1
    for chrom in sorted(lengths, key=lambda x: rank_order.index(x) if x in rank_order else 999):
        length = lengths[chrom]
        nwin = (length + window_size - 1) // window_size
        for win_id in range(1, nwin + 1):
            start = (win_id - 1) * window_size + 1
            end = min(win_id * window_size, length)
            rows.append({"chrom": chrom, "win_id": win_id, "window": f"{chrom}:{start}-{end}",
                         "win_start": start, "win_end": end, "chr_rank": rank_order.index(chrom) + 1 if chrom in rank_order else 999,
                         "genome_pos": genome_pos})
            genome_pos += end - start + 1
    return pd.DataFrame(rows)


def write_excel_sheets(path, sheets):
    try:
        import pandas as pd  # noqa
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for name, df in sheets.items():
                safe = name[:31]
                df.to_excel(writer, sheet_name=safe, index=False)
    except Exception:
        out = Path(str(path) + ".csv_sheets")
        out.mkdir(parents=True, exist_ok=True)
        for name, df in sheets.items():
            df.to_csv(out / f"{name}.csv", index=False)


def plot_percent(report_dir, group_percent):
    if group_percent.empty:
        return
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except Exception:
        return
    for group_col in [c for c in ("super_pop", "pop") if c in group_percent.columns]:
        d = group_percent.dropna(subset=[group_col]).copy()
        if d.empty:
            continue
        piv = d.pivot_table(index=group_col, columns="anc", values="segment_percent_diploid", fill_value=0)
        ax = piv.plot(kind="bar", figsize=(8, 4.8))
        ax.set_xlabel(group_col)
        ax.set_ylabel("% of diploid genome")
        ax.set_title("TRACE archaic segment burden")
        plt.tight_layout()
        plt.savefig(Path(report_dir) / f"trace_all_percent_by_{group_col}.png", dpi=220)
        plt.close()


def combine(root, sample_file=None, window_size=1000000, stats_by_group="super_pop", genome_build="GRCh37"):
    import pandas as pd
    root = Path(root)
    sample_map = pd.read_csv(root / "samples" / "trace_sample_map.tsv", sep="\t", dtype=str)
    sample_map["tree_node_id"] = sample_map["tree_node_id"].astype(str)
    rows = []
    for summary_file in sorted((root / "summary").glob("hap*.summary.txt")):
        node = summary_file.name.split(".")[0].replace("hap", "")
        try:
            x = pd.read_csv(summary_file, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            continue
        if x.empty:
            continue
        x.insert(0, "tree_node_id", str(node))
        rows.append(x)
    if rows:
        hap = pd.concat(rows, ignore_index=True)
    else:
        hap = pd.DataFrame(columns=["tree_node_id", "chrom", "start", "end", "mean_posterior", "length_bp", "length_cM"])
    hap = normalize_trace_segments(hap)
    hap = hap.merge(sample_map, on="tree_node_id", how="left")
    for c in ["sample", "haplotype", "tree_node_id", "pop", "super_pop", "anc", "chrom", "start", "end", "mean_posterior", "length_bp", "length_cM"]:
        if c not in hap.columns:
            hap[c] = pd.NA
    front = ["sample", "haplotype", "tree_node_id", "pop", "super_pop", "anc", "chrom", "start", "end", "mean_posterior", "length_bp", "length_cM"]
    hap = hap[front + [c for c in hap.columns if c not in front]]

    final = root / "final"
    report = root / "report"
    final.mkdir(parents=True, exist_ok=True)
    report.mkdir(parents=True, exist_ok=True)
    hap.to_csv(final / "trace_haplotype_segments.tsv.gz", sep="\t", index=False, compression="gzip")
    hap.to_csv(final / "trace_person_segments.tsv.gz", sep="\t", index=False, compression="gzip")

    if sample_file and Path(sample_file).exists():
        panel = read_panel(sample_file)
        sample_den = panel[["sample", "pop", "super_pop"]].copy()
    else:
        sample_den = sample_map[["sample", "pop", "super_pop"]].drop_duplicates().copy()
    sample_den = sample_den.dropna(subset=["sample"]).drop_duplicates("sample")
    chroms = sorted(set(hap["chrom"].dropna().astype(str))) if not hap.empty else [str(i) for i in range(1, 23)] + ["X"]
    denom = sum(chr_lengths_grch37(chroms).values())
    if denom <= 0:
        denom = 3031042417
    diploid = 2 * denom

    if not hap.empty:
        sample_len = hap.groupby(["sample", "anc"], dropna=False).agg(
            n_segments=("length_bp", "size"),
            total_bp=("length_bp", "sum"),
            mean_posterior=("mean_posterior", "mean") if "mean_posterior" in hap.columns else ("length_bp", "size"),
        ).reset_index()
    else:
        sample_len = pd.DataFrame(columns=["sample", "anc", "n_segments", "total_bp", "mean_posterior"])
    if sample_len.empty:
        sample_len = sample_den[["sample"]].copy()
        sample_len["anc"] = "TRACE_Archaic"
        sample_len["n_segments"] = 0
        sample_len["total_bp"] = 0.0
        sample_len["mean_posterior"] = pd.NA
    sample_len = sample_den.merge(sample_len, on="sample", how="left")
    sample_len["anc"] = sample_len["anc"].fillna("TRACE_Archaic")
    sample_len["n_segments"] = sample_len["n_segments"].fillna(0).astype(int)
    sample_len["total_bp"] = sample_len["total_bp"].fillna(0.0)
    sample_len["total_mb"] = sample_len["total_bp"] / 1e6
    sample_len["segment_percent_haploid"] = 100 * sample_len["total_bp"] / denom
    sample_len["segment_percent_diploid"] = 100 * sample_len["total_bp"] / diploid
    sample_len.to_csv(report / "trace_all_length.csv", index=False)

    group_cols = [x.strip() for x in str(stats_by_group or "").split(",") if x.strip()]
    group_cols = [g for g in group_cols if g in sample_len.columns]
    group_tables = {}
    for g in group_cols:
        gp = sample_len.dropna(subset=[g]).groupby([g, "anc"], dropna=False).agg(
            n_samples_total=("sample", "nunique"),
            n_samples_with_segment=("n_segments", lambda z: int((z > 0).sum())),
            mean_total_bp_per_sample=("total_bp", "mean"),
            median_total_bp_per_sample=("total_bp", "median"),
            mean_segments_per_sample=("n_segments", "mean"),
            segment_percent_haploid=("segment_percent_haploid", "mean"),
            segment_percent_diploid=("segment_percent_diploid", "mean"),
        ).reset_index()
        gp["sample_prevalence_percent"] = 100 * gp["n_samples_with_segment"] / gp["n_samples_total"].clip(lower=1)
        group_tables[g] = gp
        gp.to_csv(report / f"trace_all_percent_by_{g}.csv", index=False)
    if "super_pop" in group_tables:
        group_percent = group_tables["super_pop"]
    elif group_tables:
        group_percent = next(iter(group_tables.values()))
    else:
        group_percent = pd.DataFrame()
    if not group_percent.empty:
        group_percent.to_csv(report / "trace_all_percent.csv", index=False)
    else:
        pd.DataFrame().to_csv(report / "trace_all_percent.csv", index=False)
    plot_percent(report, group_percent)

    win_burden = pd.DataFrame()
    win_group = pd.DataFrame()
    if not hap.empty and window_size:
        w = int(window_size)
        wins = build_windows(hap["chrom"].dropna().unique().tolist(), w)
        if not wins.empty:
            x = hap[["sample", "anc", "chrom", "start", "end"]].copy()
            x["start_win"] = ((x["start"].astype(int) - 1) // w) + 1
            x["end_win"] = ((x["end"].astype(int) - 1) // w) + 1
            expanded = []
            for _, row in x.iterrows():
                for win_id in range(int(row["start_win"]), int(row["end_win"]) + 1):
                    expanded.append((row["sample"], row["anc"], row["chrom"], win_id))
            ov = pd.DataFrame(expanded, columns=["sample", "anc", "chrom", "win_id"]).drop_duplicates()
            ov = ov.merge(wins, on=["chrom", "win_id"], how="left").dropna(subset=["window"])
            den_all = sample_den["sample"].nunique()
            win_burden = ov.groupby(["anc", "chrom", "win_id", "window", "win_start", "win_end", "chr_rank", "genome_pos"]).agg(
                n_samples_with_segment=("sample", "nunique")
            ).reset_index()
            win_burden["n_samples_total"] = den_all
            win_burden["burden_fraction"] = win_burden["n_samples_with_segment"] / max(1, den_all)
            win_burden["burden_percent"] = 100 * win_burden["burden_fraction"]
            win_burden = win_burden.sort_values(["anc", "chr_rank", "win_start"])
            win_burden.to_csv(report / "trace_all_gw.csv", index=False)
            if group_cols:
                g = group_cols[0]
                ovg = ov.merge(sample_den[["sample", g]].drop_duplicates(), on="sample", how="left").dropna(subset=[g])
                deng = sample_den.groupby(g).agg(n_samples_total=("sample", "nunique")).reset_index()
                win_group = ovg.groupby([g, "anc", "chrom", "win_id", "window", "win_start", "win_end", "chr_rank", "genome_pos"]).agg(
                    n_samples_with_segment=("sample", "nunique")
                ).reset_index().merge(deng, on=g, how="left")
                win_group["burden_fraction"] = win_group["n_samples_with_segment"] / win_group["n_samples_total"].clip(lower=1)
                win_group["burden_percent"] = 100 * win_group["burden_fraction"]
                win_group.to_csv(report / "trace_all_gw_by_group.csv", index=False)

    stats = [
        ["haplotype_segment_rows", len(hap)],
        ["samples_with_segments", hap["sample"].dropna().nunique() if "sample" in hap else 0],
        ["tree_nodes_with_segments", hap["tree_node_id"].dropna().nunique() if "tree_node_id" in hap else 0],
        ["n_samples_denominator", sample_den["sample"].nunique()],
        ["haploid_denominator_gb", round(denom / 1e9, 4)],
        ["diploid_denominator_gb", round(diploid / 1e9, 4)],
        ["window_size", int(window_size) if window_size else "NA"],
        ["genome_build", genome_build],
    ]
    summary = pd.DataFrame(stats, columns=["metric", "value"])
    summary.to_csv(report / "trace_summary.tsv", sep="\t", index=False)

    sheets = {"trace_summary": summary, "trace_all_length": sample_len}
    if not group_percent.empty:
        sheets["trace_all_percent"] = group_percent
    if not win_burden.empty:
        sheets["trace_all_gw"] = win_burden.head(100000)
    if not win_group.empty:
        sheets["trace_all_gw_by_group"] = win_group.head(100000)
    write_excel_sheets(report / "trace_report.xlsx", sheets)


def plot_windows(report_dir, output_dir=None):
    """Plot window burden and carrier counts from trace_all_gw.csv."""
    import matplotlib.pyplot as plt
    import pandas as pd

    report = Path(report_dir)
    output = Path(output_dir) if output_dir else report / "figures"
    gw_file = report / "trace_all_gw.csv"
    if not gw_file.is_file():
        raise SystemExit(f"trace_all_gw.csv not found: {gw_file}")
    df = pd.read_csv(gw_file)
    if df.empty:
        raise SystemExit(f"no TRACE windows in: {gw_file}")
    output.mkdir(parents=True, exist_ok=True)

    # A region-focused run normally contains one chromosome. For multiple
    # chromosomes, use genome_pos to avoid joining chromosome ends together.
    one_chrom = df["chrom"].astype(str).nunique() == 1
    x = df["win_start"] / 1e6 if one_chrom else df["genome_pos"] / 1e6
    xlabel = (f'Position on {df["chrom"].iloc[0]} (Mb)' if one_chrom
              else "Genome position (Mb)")
    peak = df.loc[df["burden_percent"].idxmax()]

    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(x, df["burden_percent"], color="#2E86C1", linewidth=1.8)
    ax.scatter([x.loc[peak.name]], [peak["burden_percent"]], color="darkred", zorder=3)
    ax.annotate(f'Peak {peak["burden_percent"]:.1f}%\n{peak["window"]}',
                xy=(x.loc[peak.name], peak["burden_percent"]), xytext=(12, -35),
                textcoords="offset points", arrowprops={"arrowstyle": "->", "color": "darkred"})
    ax.set(xlabel=xlabel, ylabel="Sample burden (%)", title="TRACE archaic segment burden")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    burden_file = output / "burden_line.png"
    fig.savefig(burden_file, dpi=300)
    plt.close(fig)

    window_bp = int((df["win_end"] - df["win_start"] + 1).median())
    window_label = (f"{window_bp / 1000:g} kb" if window_bp < 1_000_000
                    else f"{window_bp / 1_000_000:g} Mb")
    fig, ax = plt.subplots(figsize=(10, 4))
    carriers = df["n_samples_with_segment"]
    ax.hist(carriers, bins=25, color="#28B463", edgecolor="black", alpha=0.7)
    ax.axvline(carriers.median(), color="blue", linestyle="--",
               label=f"Median = {carriers.median():.0f}")
    ax.set(xlabel=f"Samples carrying a segment per {window_label} window",
           ylabel="Number of windows", title="TRACE segment-carrier distribution")
    ax.legend()
    fig.tight_layout()
    carrier_file = output / "carrier_hist.png"
    fig.savefig(carrier_file, dpi=300)
    plt.close(fig)
    print(f"TRACE window plots written to {output}")


def gene_overlap(report_dir, gene_annot, target_genes, output=None):
    """Overlap TRACE burden windows with selected hg19 UCSC refGene genes."""
    import gzip
    import shutil
    import urllib.request
    import pandas as pd

    report = Path(report_dir)
    annotation = Path(gene_annot)
    if not annotation.exists():
        annotation.parent.mkdir(parents=True, exist_ok=True)
        download = annotation.with_suffix(annotation.suffix + ".gz")
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
        print(f"Downloading hg19 refGene annotation to {annotation}")
        try:
            urllib.request.urlretrieve(url, download)
            with gzip.open(download, "rb") as src, annotation.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            download.unlink()
        except Exception as exc:
            raise SystemExit(f"failed to download gene annotation: {exc}") from exc

    raw = pd.read_csv(annotation, sep="\t", header=None, usecols=[2, 4, 5, 12])
    raw.columns = ["chrom", "start", "end", "gene"]
    raw["chrom"] = raw["chrom"].astype(str).str.replace("^chr", "", regex=True)
    genes = raw.groupby(["chrom", "gene"], as_index=False).agg(start=("start", "min"), end=("end", "max"))
    wanted = [x.strip() for x in target_genes.split(",") if x.strip()]
    genes = genes[genes["gene"].isin(wanted)]
    if genes.empty:
        raise SystemExit(f"none of the requested genes were found: {', '.join(wanted)}")

    gw_file = report / "trace_all_gw.csv"
    if not gw_file.is_file():
        raise SystemExit(f"trace_all_gw.csv not found: {gw_file}")
    windows = pd.read_csv(gw_file)
    windows["chrom"] = windows["chrom"].astype(str).str.replace("^chr", "", regex=True)
    rows = []
    for _, window in windows.iterrows():
        matched = genes[(genes["chrom"] == window["chrom"]) &
                        (genes["start"] < window["win_end"]) &
                        (genes["end"] > window["win_start"] - 1)]
        if not matched.empty:
            row = window.to_dict()
            row["genes"] = ",".join(sorted(matched["gene"].unique()))
            row["gene_count"] = matched["gene"].nunique()
            rows.append(row)
    result = pd.DataFrame(rows)
    output_path = Path(output) if output else report / "trace_gene_overlap.csv"
    result.to_csv(output_path, index=False)
    print(f"TRACE gene overlap written to {output_path} ({len(result)} windows)")


def main():
    parser = argparse.ArgumentParser(description="Helper commands for integrated ibdmix.sh TRACE secondary analysis")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("normalize-sample-map")
    p.add_argument("--sample-map", required=True)
    p.add_argument("--sample-file", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("auto-sample-map")
    p.add_argument("--tree-file", required=True)
    p.add_argument("--sample-file", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("combine")
    p.add_argument("--root", required=True)
    p.add_argument("--sample-file")
    p.add_argument("--window-size", type=int, default=1000000)
    p.add_argument("--stats-by-group", default="super_pop")
    p.add_argument("--genome-build", default="GRCh37")
    p = sub.add_parser("plot", help="plot TRACE window burden and carrier counts")
    p.add_argument("--report-dir", required=True)
    p.add_argument("--output-dir")
    p = sub.add_parser("gene-overlap", help="overlap TRACE windows with selected hg19 refGene genes")
    p.add_argument("--report-dir", required=True)
    p.add_argument("--target-genes", default="OSM,LIF,GAL3ST1,NF2")
    p.add_argument("--gene-annot", default="refGene.txt")
    p.add_argument("--output")
    args = parser.parse_args()
    if args.cmd == "normalize-sample-map":
        normalize_sample_map(args.sample_map, args.sample_file, args.out)
    elif args.cmd == "auto-sample-map":
        auto_sample_map(args.tree_file, args.sample_file, args.out)
    elif args.cmd == "combine":
        combine(args.root, sample_file=args.sample_file, window_size=args.window_size,
                stats_by_group=args.stats_by_group, genome_build=args.genome_build)
    elif args.cmd == "plot":
        plot_windows(args.report_dir, output_dir=args.output_dir)
    elif args.cmd == "gene-overlap":
        gene_overlap(args.report_dir, args.gene_annot, args.target_genes, output=args.output)


if __name__ == "__main__":
    main()

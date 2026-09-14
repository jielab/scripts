#!/usr/bin/env python3
import csv
import gzip
import math
import os
import re
import shutil
import subprocess
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


def norm_chrom(x):
    chrom = re.sub(r"^chr", "", str(x), flags=re.IGNORECASE).upper()
    chrom = re.sub(r"\.0$", "", chrom)
    chrom = {"X": "23", "Y": "24", "M": "25", "MT": "25"}.get(chrom, chrom)
    return chrom.lstrip("0") or "0"


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


def sumstat_reader(path):
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
    return {
        "delimiter": delim,
        "usecols": set(usecols),
        "snp": c_snp,
        "chr": c_chr,
        "pos": c_pos,
        "ea": c_ea,
        "nea": c_nea,
        "beta": c_beta,
        "se": c_se,
        "p": c_p,
    }


def normalize_sumstat_chunk(ch, spec):
    ch = ch.rename(columns={spec["chr"]: "CHR", spec["pos"]: "POS", spec["ea"]: "EA",
                            spec["nea"]: "NEA", spec["beta"]: "BETA", spec["se"]: "SE"})
    if spec["snp"]:
        ch = ch.rename(columns={spec["snp"]: "SNP"})
    else:
        ch["SNP"] = ""
    if spec["p"]:
        ch = ch.rename(columns={spec["p"]: "P"})
    else:
        ch["P"] = pd.NA
    ch["CHR"] = ch["CHR"].map(norm_chrom)
    ch["POS"] = pd.to_numeric(ch["POS"], errors="coerce")
    ch["BETA"] = pd.to_numeric(ch["BETA"], errors="coerce")
    ch["SE"] = pd.to_numeric(ch["SE"], errors="coerce")
    ch["P"] = pd.to_numeric(ch["P"], errors="coerce")
    ch = ch.dropna(subset=["CHR", "POS", "BETA", "SE"])
    ch = ch[ch["SE"] > 0]
    ch["EA"] = ch["EA"].astype(str).str.upper()
    ch["NEA"] = ch["NEA"].astype(str).str.upper()
    return ch[(ch["EA"] != "") & (ch["NEA"] != "")]


def iter_sumstat_chunks(path, chunksize=500000):
    spec = sumstat_reader(path)
    for ch in pd.read_csv(path, sep=spec["delimiter"], usecols=lambda c: c in spec["usecols"],
                          chunksize=chunksize, compression="infer"):
        ch = normalize_sumstat_chunk(ch, spec)
        if not ch.empty:
            yield ch


def tabix_index(path):
    for suffix in (".tbi", ".csi"):
        index = str(path) + suffix
        if os.path.isfile(index) and os.path.getsize(index) > 0 and os.path.getmtime(index) >= os.path.getmtime(path):
            return index
    return None


def load_regions_tabix(path, interval_map, label=None, chunksize=500000):
    tabix = shutil.which("tabix")
    if not tabix or not tabix_index(path):
        return None
    contigs = subprocess.check_output([tabix, "-l", str(path)], text=True).splitlines()
    contig_map = {norm_chrom(x): x for x in contigs}
    region_args = []
    for chrom, (starts, ends) in interval_map.items():
        source_chrom = contig_map.get(norm_chrom(chrom))
        if source_chrom is None:
            continue
        region_args.extend(f"{source_chrom}:{int(start)}-{int(end)}" for start, end in zip(starts, ends))
    if label:
        print(f"{label}: tabix query {path} for {len(region_args)} merged intervals", flush=True)
    if not region_args:
        return empty_sumstats()
    spec = sumstat_reader(path)
    cols, _ = read_header(path)
    proc = subprocess.Popen([tabix, str(path), *region_args], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    chunks = []
    try:
        for ch in pd.read_csv(proc.stdout, sep=spec["delimiter"], names=cols, header=None,
                              usecols=lambda c: c in spec["usecols"], chunksize=chunksize):
            ch = normalize_sumstat_chunk(ch, spec)
            if not ch.empty:
                chunks.append(ch)
    except pd.errors.EmptyDataError:
        pass
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
    rc = proc.wait()
    if rc:
        raise RuntimeError(f"tabix exited with status {rc}: {path}")
    result = pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()
    if label:
        print(f"{label}: tabix retained {len(result):,} usable rows", flush=True)
    return result


def empty_sumstats():
    return pd.DataFrame(columns=["SNP", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P"])


def merge_regions(regions):
    """Return non-overlapping intervals grouped by canonical chromosome."""
    grouped = {}
    for chrom, start, end in regions:
        chrom = norm_chrom(chrom)
        start, end = int(float(start)), int(float(end))
        if end < start:
            start, end = end, start
        grouped.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, intervals in grouped.items():
        out = []
        for start, end in sorted(intervals):
            if out and start <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = (
            np.asarray([x[0] for x in out], dtype=np.int64),
            np.asarray([x[1] for x in out], dtype=np.int64),
        )
    return merged


def load_regions(path, regions, label=None):
    """Scan a summary-statistics file once and retain the union of regions.

    The input may be ordinary gzip, so this is the bounded-memory fallback for
    files that cannot be queried with tabix. Merged intervals and binary
    searches avoid testing every input row against every requested locus.
    """
    interval_map = merge_regions(regions)
    if not interval_map:
        return empty_sumstats()
    try:
        indexed = load_regions_tabix(path, interval_map, label=label)
        if indexed is not None:
            return indexed
    except Exception as exc:
        print(f"WARNING: tabix query failed; falling back to one full scan: {exc}", file=sys.stderr, flush=True)
    if label:
        n_intervals = sum(len(x[0]) for x in interval_map.values())
        print(f"{label}: scanning {path} once for {n_intervals} merged intervals", flush=True)
    chunks = []
    input_rows = 0
    for ch in iter_sumstat_chunks(path):
        input_rows += len(ch)
        chrom_values = ch["CHR"].to_numpy(dtype=str)
        positions = ch["POS"].to_numpy(dtype=np.float64)
        keep = np.zeros(len(ch), dtype=bool)
        for chrom in np.unique(chrom_values):
            interval = interval_map.get(chrom)
            if interval is None:
                continue
            row_index = np.flatnonzero(chrom_values == chrom)
            pos = positions[row_index]
            starts, ends = interval
            interval_index = np.searchsorted(starts, pos, side="right") - 1
            valid = interval_index >= 0
            selected = np.zeros(len(row_index), dtype=bool)
            selected[valid] = pos[valid] <= ends[interval_index[valid]]
            keep[row_index[selected]] = True
        if keep.any():
            chunks.append(ch.loc[keep].copy())
    result = pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()
    if label:
        print(f"{label}: retained {len(result):,} of {input_rows:,} usable rows", flush=True)
    return result


def select_region(d, chrom, start, end):
    if d.empty:
        return d.copy()
    chrom = norm_chrom(chrom)
    return d[(d["CHR"].astype(str) == chrom) & (d["POS"] >= start) & (d["POS"] <= end)].copy()


def load_region(path, chrom=None, start=None, end=None):
    if chrom is not None:
        return select_region(load_regions(path, [(chrom, start, end)]), chrom, start, end)
    chunks = list(iter_sumstat_chunks(path))
    return pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()


def find_lead_and_region(path, region, window_kb):
    if region and re.match(r"^[A-Za-z0-9]+:\d+-\d+$", region):
        chrom, rest = region.split(":", 1)
        start, end = rest.split("-", 1)
        return re.sub(r"^chr", "", chrom, flags=re.IGNORECASE), int(float(start)), int(float(end))
    d = load_region(path)
    if d.empty:
        raise RuntimeError("No usable QTL rows in " + path)
    if d["P"].notna().any():
        lead = d.loc[d["P"].idxmin()]
    else:
        z = (d["BETA"] / d["SE"]).abs()
        lead = d.loc[z.idxmax()]
    chrom = re.sub(r"^chr", "", str(lead["CHR"]), flags=re.IGNORECASE)
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

    qtl_summary, cad_summary, expected_pairs, status_rows = [], [], [], []
    records = []
    qtl_groups = {}
    print(f"GPU-coloc preparation: {len(rows)} manifest rows", flush=True)
    for row_index, row in enumerate(rows):
        omics = row.get("omics", "protein") or "protein"
        trait = row.get("trait", "")
        qtl_file = row.get("file", row.get("exposure", ""))
        typ = row.get("type", "quant") or "quant"
        region = row.get("region", "") or ""
        record = {"row_index": row_index, "omics": omics, "trait": trait, "qtl_file": qtl_file,
                  "type": typ, "input_region": region}
        records.append(record)
        if not trait or not qtl_file or not os.path.exists(qtl_file):
            record["status"] = "not_run"
            record["message"] = "missing qtl file"
            continue
        try:
            chrom, start, end = find_lead_and_region(qtl_file, region, window_kb)
            chrom = norm_chrom(chrom)
            record.update({"chrom": chrom, "start": int(start), "end": int(end),
                           "region_key": f"chr{chrom}:{int(start)}-{int(end)}"})
            qtl_groups.setdefault(qtl_file, []).append(record)
        except Exception as e:
            record["status"] = "failed"
            record["message"] = str(e)

    # Each QTL file commonly supplies up to three loci. Read its requested
    # union once instead of decompressing the same file once per locus.
    for group_index, (qtl_file, group) in enumerate(qtl_groups.items(), start=1):
        try:
            regions = [(x["chrom"], x["start"], x["end"]) for x in group]
            qtl_union = load_regions(qtl_file, regions)
            for record in group:
                qtl = select_region(qtl_union, record["chrom"], record["start"], record["end"])
                qtl_signal = safe_name(
                    f"QTL__{record['omics']}__{record['trait']}__chr{record['chrom']}_{record['start']}_{record['end']}"
                )
                qtl_s, qtl_lead, qtl_strength = to_signal(qtl, qtl_signal, record["type"])
                feather.write_feather(qtl_s, qtl_dir / f"{qtl_signal}.feather")
                qtl_summary.append({"signal": qtl_signal, "chromosome": record["chrom"],
                                    "location_min": record["start"], "location_max": record["end"],
                                    "signal_strength": qtl_strength, "lead_variant": qtl_lead,
                                    "omics": record["omics"], "trait": record["trait"],
                                    "source_file": qtl_file})
                record["qtl_signal"] = qtl_signal
                record["qtl_n"] = int(qtl_s.shape[1])
        except Exception as e:
            for record in group:
                if "qtl_signal" not in record:
                    record["status"] = "failed"
                    record["message"] = f"QTL preparation failed: {e}"
        if group_index % 25 == 0 or group_index == len(qtl_groups):
            print(f"GPU-coloc preparation: QTL files {group_index}/{len(qtl_groups)}", flush=True)

    # The previous implementation called load_region(cad_gwas, ...) for every
    # unique locus. With hundreds of loci that decompressed the same outcome
    # file hundreds of times. Retain the union in one sequential scan.
    unique_regions = {}
    for record in records:
        if "qtl_signal" in record:
            unique_regions.setdefault(record["region_key"],
                                      (record["chrom"], record["start"], record["end"]))
    cad_cache, cad_errors = {}, {}
    if unique_regions:
        try:
            cad_union = load_regions(cad_gwas, list(unique_regions.values()), label="Outcome GWAS")
            for region_index, (region_key, (chrom, start, end)) in enumerate(unique_regions.items(), start=1):
                try:
                    cad = select_region(cad_union, chrom, start, end)
                    # One outcome signal per exact region is sufficient. The old
                    # trait-specific copies created a large cross-product of
                    # identical CAD signals whenever several omics traits shared a locus.
                    cad_signal = safe_name(f"CAD__chr{chrom}_{start}_{end}")
                    cad_s, cad_lead, cad_strength = to_signal(cad, cad_signal, outcome_type)
                    feather.write_feather(cad_s, cad_dir / f"{cad_signal}.feather")
                    cad_summary.append({"signal": cad_signal, "chromosome": chrom,
                                        "location_min": start, "location_max": end,
                                        "signal_strength": cad_strength, "lead_variant": cad_lead,
                                        "source_file": cad_gwas})
                    cad_cache[region_key] = {"signal": cad_signal, "n": int(cad_s.shape[1])}
                except Exception as e:
                    cad_errors[region_key] = str(e)
                if region_index % 100 == 0 or region_index == len(unique_regions):
                    print(f"GPU-coloc preparation: outcome signals {region_index}/{len(unique_regions)}", flush=True)
        except Exception as e:
            for region_key in unique_regions:
                cad_errors[region_key] = f"Outcome GWAS preparation failed: {e}"

    for record in records:
        base_status = {"omics": record["omics"], "trait": record["trait"],
                       "region": record.get("region_key", record["input_region"])}
        if "qtl_signal" not in record:
            status_rows.append({**base_status, "status": record.get("status", "failed"),
                                "message": record.get("message", "QTL preparation failed")})
            continue
        cad_info = cad_cache.get(record["region_key"])
        if cad_info is None:
            status_rows.append({**base_status, "status": "failed",
                                "message": cad_errors.get(record["region_key"], "outcome signal missing")})
            continue
        expected_pairs.append({"omics": record["omics"], "trait": record["trait"],
                               "region": record["region_key"], "qtl_signal": record["qtl_signal"],
                               "cad_signal": cad_info["signal"]})
        status_rows.append({**base_status, "status": "ok",
                            "message": f"qtl_n={record['qtl_n']}; cad_n={cad_info['n']}"})

    pd.DataFrame(qtl_summary).to_csv(outdir / "qtl_summary.tsv", sep="\t", index=False)
    pd.DataFrame(cad_summary).to_csv(outdir / "cad_summary.tsv", sep="\t", index=False)
    pd.DataFrame(expected_pairs).to_csv(outdir / "expected_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(status_rows).to_csv(outdir / "signal_preparation_status.tsv", sep="\t", index=False)
    if len(qtl_summary) == 0 or len(cad_summary) == 0:
        raise SystemExit("No GPU-coloc signal files were created. See signal_preparation_status.tsv")
    print(f"GPU-coloc preparation complete: {len(qtl_summary)} QTL signals, "
          f"{len(cad_summary)} outcome signals", flush=True)


def filter_results(raw_results, expected_pairs, output_file, h4_threshold=0.70):
    """Keep only the QTL/outcome pair requested by each C3 manifest row.

    gpu-coloc intentionally tests every overlapping signal pair. C3 creates a
    disease signal for each requested region, so overlapping regions can also
    create scientifically unintended cross-pairs. Retain the exact manifest
    pairs and explicitly record pairs for which gpu-coloc returned no result.
    """
    try:
        h4_threshold = float(h4_threshold)
    except (TypeError, ValueError):
        raise SystemExit("GPU-coloc H4 threshold must be numeric")
    if not 0 <= h4_threshold <= 1:
        raise SystemExit("GPU-coloc H4 threshold must be between 0 and 1")
    expected = pd.read_csv(expected_pairs, sep="\t", dtype=str)
    if expected.empty:
        raise SystemExit("No expected GPU-coloc pairs were recorded: " + expected_pairs)
    raw = pd.DataFrame()
    if os.path.exists(raw_results) and os.path.getsize(raw_results) > 0:
        try:
            raw = pd.read_csv(raw_results, sep="\t")
        except pd.errors.EmptyDataError:
            raw = pd.DataFrame()
    if not raw.empty and {"signal1", "signal2"}.issubset(raw.columns):
        pp_col = "PP.H4" if "PP.H4" in raw.columns else None
        if pp_col:
            raw[pp_col] = pd.to_numeric(raw[pp_col], errors="coerce")
            raw = raw.sort_values(pp_col, ascending=False, na_position="last")
        raw = raw.drop_duplicates(["signal1", "signal2"], keep="first")
        out = expected.merge(raw, how="left", left_on=["qtl_signal", "cad_signal"],
                             right_on=["signal1", "signal2"], indicator="_gpu_merge")
    else:
        out = expected.copy()
        out["signal1"] = out["qtl_signal"]
        out["signal2"] = out["cad_signal"]
        out["PP.H4"] = np.nan
        out["_gpu_merge"] = "left_only"
    if "signal1" in out:
        out["signal1"] = out["signal1"].fillna(out["qtl_signal"])
    if "signal2" in out:
        out["signal2"] = out["signal2"].fillna(out["cad_signal"])
    if "PP.H4" not in out:
        out["PP.H4"] = np.nan
    pp = pd.to_numeric(out["PP.H4"], errors="coerce")
    returned = out["_gpu_merge"].astype(str) == "both"
    finite = np.isfinite(pp)
    out["GPU_status"] = np.select(
        [~returned, returned & finite],
        ["no_gpu_result", "ok"],
        default="gpu_result_na",
    )
    threshold_label = f"{h4_threshold:g}"
    out["GPU_H4_class"] = np.select(
        [~returned, returned & ~finite, pp >= h4_threshold],
        ["no_gpu_result", "no_finite_pp_h4", f"H4>={threshold_label}"],
        default=f"H4<{threshold_label}",
    )
    out["GPU_H4_threshold"] = h4_threshold
    out = out.drop(columns=["_gpu_merge"])
    out.to_csv(output_file, sep="\t", index=False)


def main():
    if len(sys.argv) in (5, 6) and sys.argv[1] == "filter-results":
        threshold = sys.argv[5] if len(sys.argv) == 6 else 0.70
        filter_results(sys.argv[2], sys.argv[3], sys.argv[4], threshold)
        return
    if len(sys.argv) not in (5, 6):
        raise SystemExit("Usage: c3_coloc_GPU.py <manifest.tsv> <outcome_gwas> <outdir> <window_kb> [cc|quant]\n"
                         "   or: c3_coloc_GPU.py filter-results <raw.tsv> <expected_pairs.tsv> <output.tsv> [H4_threshold]")
    outcome_type = sys.argv[5].lower() if len(sys.argv) == 6 else "cc"
    if outcome_type not in {"cc", "quant"}:
        raise SystemExit("outcome_type must be cc or quant")
    prepare_signals(sys.argv[1], sys.argv[2], sys.argv[3], int(float(sys.argv[4])), outcome_type)


if __name__ == "__main__":
    main()

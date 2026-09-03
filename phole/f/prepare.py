#!/usr/bin/env python3
"""Build a gene-by-protein PHOLE score matrix from MAGMA gene results."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    atomic_dataframe,
    atomic_json,
    chromosome_order,
    cis_rows,
    log,
    normalize_chr,
    robust_center_scale,
    str2bool,
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwas-root", required=True)
    parser.add_argument("--category", default="common")
    parser.add_argument("--protein-bed", required=True)
    parser.add_argument("--protein-map", default="")
    parser.add_argument("--gwas", default="")
    parser.add_argument("--out", required=True)
    parser.add_argument("--score", choices=["zstat", "logp"], default="zstat")
    parser.add_argument("--scale", choices=["robust", "none"], default="robust")
    parser.add_argument("--score-cap", type=float, default=12.0)
    parser.add_argument("--cis-flank", type=int, default=100_000)
    parser.add_argument(
        "--protein-filter",
        choices=["none", "any-bonf", "cis-bonf", "coding-bonf", "target-bonf"],
        default="target-bonf",
        help=(
            "protein-level MAGMA gate: none, any genome-wide Bonferroni-significant gene, "
            "a significant gene in the protein interval plus --cis-flank, or a significant "
            "gene overlapping the protein coding interval; target-bonf specifically requires "
            "the measured protein's mapped coding gene"
        ),
    )
    parser.add_argument("--protein-alpha", type=float, default=0.05)
    parser.add_argument("--protein-build", choices=["auto", "37", "38"], default="auto")
    parser.add_argument("--require-magma-done", default="TRUE")
    parser.add_argument("--replace", default="FALSE")
    return parser.parse_args()


def discover_magma(
    root: Path,
    category: str,
    selected: set[str] | None,
    require_done: bool,
) -> list[tuple[str, Path]]:
    # Current gwas_post.sh layout:
    #   <project>/<category>/<GWAS>/magma/<GWAS>.genes.out
    # Retain compatibility with the earlier layout:
    #   <project>/<GWAS>/<category>/magma/<GWAS>.genes.out
    current = root.glob(f"{category}/*/magma/*.genes.out")
    legacy = root.glob(f"*/{category}/magma/*.genes.out")
    files = sorted(set(current) | set(legacy))
    pairs: list[tuple[str, Path]] = []
    skipped_incomplete = 0
    for path in files:
        protein = path.name.removesuffix(".genes.out")
        if selected is not None and protein not in selected:
            continue
        if require_done and not (path.parent / "magma.done").is_file():
            skipped_incomplete += 1
            continue
        if path.parents[2].name == category:
            expected = path.parents[1].name
        elif path.parents[1].name == category:
            expected = path.parents[2].name
        else:
            expected = protein
        if expected != protein:
            log(f"WARNING: trait directory '{expected}' differs from MAGMA prefix '{protein}': {path}")
        pairs.append((protein, path))
    if skipped_incomplete:
        log(f"Skip {skipped_incomplete} MAGMA file(s) without magma.done")
    if not pairs:
        raise SystemExit(
            "ERROR: no MAGMA files found under either "
            f"{root}/{category}/<GWAS>/magma or {root}/<GWAS>/{category}/magma"
        )
    return pairs


def read_bed(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep=r"\s+", comment="#", header=None, dtype=str)
    if raw.shape[1] < 4:
        raise SystemExit(f"ERROR: protein BED needs at least CHR START END PROTEIN: {path}")
    raw = raw.iloc[:, : min(raw.shape[1], 5)].copy()
    header_like = (
        raw.iloc[0, 0].upper().lstrip("#") in {"CHR", "CHROM", "CHROMOSOME"}
        or raw.iloc[0, 1].upper() == "START"
    )
    header = [str(x).upper().lstrip("#") for x in raw.iloc[0].tolist()] if header_like else []
    if header_like:
        raw = raw.iloc[1:].copy()
    gene = None
    if raw.shape[1] >= 5:
        fifth = raw.iloc[:, 4].copy()
        if header_like and header[4] in {"GENE", "GENE_SYMBOL", "HGNC", "HGNC.SYMBOL"}:
            gene = fifth
        elif not header_like:
            # Headerless five-column files are ambiguous. A genuine gene column
            # normally has high cardinality; the local protein BED instead stores a
            # small repeated Olink panel label in column five.
            unique_fifth = fifth.dropna().nunique()
            if unique_fifth > max(25, int(0.10 * len(fifth))):
                gene = fifth
            else:
                log(
                    f"WARNING: ignore low-cardinality BED column 5 ({unique_fifth} values); "
                    "treat it as annotation rather than GENE"
                )
    raw = raw.iloc[:, :4].copy()
    raw.columns = ["CHR", "START", "END", "PROTEIN"]
    raw["CHR"] = raw["CHR"].map(normalize_chr)
    raw["START"] = pd.to_numeric(raw["START"], errors="coerce")
    raw["END"] = pd.to_numeric(raw["END"], errors="coerce")
    raw["POS"] = (raw["START"] + raw["END"]) / 2.0
    raw["PROTEIN"] = raw["PROTEIN"].astype(str)
    raw["GENE"] = raw["PROTEIN"] if gene is None else gene.to_numpy()
    raw["GENE"] = raw["GENE"].fillna(raw["PROTEIN"]).astype(str)
    raw = raw.drop_duplicates("PROTEIN", keep="first")
    return raw


def read_protein_map(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as handle:
        first_data = next((line for line in handle if line.strip() and not line.startswith("#")), "")
    separator = "\t" if "\t" in first_data else r"\s+"
    mapping = pd.read_csv(path, sep=separator, comment="#", dtype=str)
    mapping.columns = [str(x).upper().lstrip("#") for x in mapping.columns]
    if {"ASSAY", "HGNC.SYMBOL"}.issubset(mapping.columns):
        mapping = mapping.rename(columns={"ASSAY": "PROTEIN", "HGNC.SYMBOL": "GENE"})
    if not {"PROTEIN", "GENE"}.issubset(mapping.columns):
        mapping = pd.read_csv(path, sep=separator, comment="#", header=None, dtype=str)
        if mapping.shape[1] < 2:
            raise SystemExit(f"ERROR: protein map needs PROTEIN and GENE columns: {path}")
        mapping = mapping.iloc[:, :2]
        mapping.columns = ["PROTEIN", "GENE"]
    return mapping[["PROTEIN", "GENE"]].drop_duplicates("PROTEIN", keep="first")


def magma_table(path: Path) -> pd.DataFrame:
    dat = pd.read_csv(path, sep=r"\s+", comment="#")
    dat.columns = [str(x).upper().lstrip("#") for x in dat.columns]
    required = {"GENE", "CHR", "START", "STOP", "P"}
    missing = sorted(required - set(dat.columns))
    if missing:
        raise SystemExit(f"ERROR: MAGMA file misses {','.join(missing)}: {path}")
    if dat["GENE"].duplicated().any():
        raise SystemExit(f"ERROR: duplicated GENE identifiers in MAGMA file: {path}")
    return dat


def magma_source_table(path: Path) -> pd.DataFrame:
    wanted = {"GENE", "CHR", "START", "STOP", "P"}
    dat = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        usecols=lambda column: str(column).upper().lstrip("#") in wanted,
    )
    dat.columns = [str(x).upper().lstrip("#") for x in dat.columns]
    missing = sorted(wanted - set(dat.columns))
    if missing:
        raise SystemExit(f"ERROR: MAGMA file misses {','.join(missing)}: {path}")
    dat["GENE"] = dat["GENE"].astype(str)
    if dat["GENE"].duplicated().any():
        raise SystemExit(f"ERROR: duplicated GENE identifiers in MAGMA file: {path}")
    return dat


def read_magma_meta(path: Path) -> dict[str, str]:
    meta = path.parent / "magma.meta.tsv"
    if not meta.exists():
        return {}
    values: dict[str, str] = {}
    try:
        with meta.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t", maxsplit=1)
                if len(fields) == 2 and fields[0].strip().lower() != "key":
                    values[fields[0].strip().lower()] = fields[1].strip()
    except Exception:
        return {}
    return values


def read_magma_build(path: Path) -> str | None:
    return read_magma_meta(path).get("grch")


def read_gene_loc_symbols(paths: set[Path]) -> dict[str, set[str]]:
    symbol_to_ids: dict[str, set[str]] = {}
    for path in sorted(paths):
        if not path.is_file():
            log(f"WARNING: MAGMA gene-loc is unavailable for target-gene ID mapping: {path}")
            continue
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) >= 6:
                    symbol_to_ids.setdefault(fields[5].upper(), set()).add(fields[0])
    return symbol_to_ids


def infer_bed_build(path: Path) -> str | None:
    name = path.name.lower()
    hit = re.search(r"(?:^|[._-])b?(37|38)(?:[._-]|$)", name)
    return hit.group(1) if hit else None


def main() -> None:
    args = arguments()
    out = Path(args.out)
    replace = str2bool(args.replace)
    done = out / "prepare.done"
    if done.exists() and (out / "score.npy").exists() and not replace:
        log(f"PHOLE matrix exists: {out / 'score.npy'}")
        return
    if args.cis_flank < 0:
        raise SystemExit("ERROR: --cis-flank must be non-negative")
    if not 0 < args.protein_alpha < 1:
        raise SystemExit("ERROR: --protein-alpha must be in (0,1)")

    selected = {x.strip() for x in args.gwas.split(",") if x.strip()} or None
    require_magma_done = str2bool(args.require_magma_done)
    pairs = discover_magma(Path(args.gwas_root), args.category, selected, require_magma_done)
    bed = read_bed(Path(args.protein_bed))
    if args.protein_map:
        mapping = read_protein_map(Path(args.protein_map))
        bed = bed.drop(columns=["GENE"]).merge(mapping, on="PROTEIN", how="left")
        bed["GENE"] = bed["GENE"].fillna(bed["PROTEIN"])

    bed_exact = bed.set_index("PROTEIN", drop=False)
    bed_upper = bed.assign(PROTEIN_UPPER=bed["PROTEIN"].str.upper()).drop_duplicates("PROTEIN_UPPER").set_index("PROTEIN_UPPER")
    target_rows = []
    retained = []
    for protein, path in pairs:
        if protein in bed_exact.index:
            row = bed_exact.loc[protein]
        elif protein.upper() in bed_upper.index:
            row = bed_upper.loc[protein.upper()]
            log(f"WARNING: case-insensitive protein BED match: {protein} -> {row['PROTEIN']}")
        else:
            log(f"WARNING: no coding-gene position; skip protein: {protein}")
            continue
        target_rows.append({
            "PROTEIN": protein,
            "GENE": str(row["GENE"]),
            "CHR": normalize_chr(row["CHR"]),
            "START": float(row["START"]),
            "END": float(row["END"]),
            "POS": float(row["POS"]),
            "MAGMA_FILE": str(path),
        })
        retained.append((protein, path))
    if len(retained) < 1:
        raise SystemExit("ERROR: no MAGMA traits matched the protein BED")
    pairs = retained
    target_candidates = pd.DataFrame(target_rows)

    builds = {x for _, path in pairs if (x := read_magma_build(path)) is not None}
    if len(builds) > 1:
        raise SystemExit(f"ERROR: mixed MAGMA genome builds: {sorted(builds)}")
    magma_build = next(iter(builds), None)
    bed_build = infer_bed_build(Path(args.protein_bed)) if args.protein_build == "auto" else args.protein_build
    if magma_build and bed_build and magma_build != bed_build:
        raise SystemExit(
            f"ERROR: MAGMA is GRCh{magma_build}, but protein BED appears GRCh{bed_build}: {args.protein_bed}"
        )
    if magma_build and not bed_build:
        log(f"WARNING: protein BED build not encoded in filename; assuming GRCh{magma_build}")
        bed_build = magma_build
    gene_loc_paths = {
        Path(value)
        for _, path in pairs
        if (value := read_magma_meta(path).get("gene_loc")) is not None
    }
    gene_symbol_to_ids = read_gene_loc_symbols(gene_loc_paths)

    # Audit every protein before it enters the graph. Significance is defined
    # within each MAGMA file using alpha / number of finite gene P values.
    # cis-bonf uses interval overlap against the measured protein's coding
    # interval expanded by --cis-flank. The audit includes excluded proteins.
    out.mkdir(parents=True, exist_ok=True)
    source_parts = []
    seen_source: set[str] = set()
    protein_qc_rows = []
    filtered_pairs = []
    filtered_target_rows = []
    for index, ((protein, path), target_row) in enumerate(zip(pairs, target_rows)):
        dat_source = magma_source_table(path)
        p_value = pd.to_numeric(dat_source["P"], errors="coerce").to_numpy(float)
        valid_p = np.isfinite(p_value) & (p_value >= 0) & (p_value <= 1)
        n_gene = int(valid_p.sum())
        bonf_threshold = args.protein_alpha / n_gene if n_gene else np.nan
        significant = valid_p & (p_value <= bonf_threshold) if n_gene else np.zeros(len(p_value), dtype=bool)

        source_chr_qc = dat_source["CHR"].map(normalize_chr).to_numpy(object)
        source_start_qc = pd.to_numeric(dat_source["START"], errors="coerce").to_numpy(float)
        source_stop_qc = pd.to_numeric(dat_source["STOP"], errors="coerce").to_numpy(float)
        target_chr_qc = normalize_chr(target_row["CHR"])
        target_start_qc = float(target_row["START"]) - args.cis_flank
        target_end_qc = float(target_row["END"]) + args.cis_flank
        cis = (
            (source_chr_qc == target_chr_qc)
            & np.isfinite(source_start_qc)
            & np.isfinite(source_stop_qc)
            & (source_stop_qc >= target_start_qc)
            & (source_start_qc <= target_end_qc)
        )
        coding = (
            (source_chr_qc == target_chr_qc)
            & np.isfinite(source_start_qc)
            & np.isfinite(source_stop_qc)
            & (source_stop_qc >= float(target_row["START"]))
            & (source_start_qc <= float(target_row["END"]))
        )
        cis_valid = cis & valid_p
        cis_significant = cis & significant
        coding_valid = coding & valid_p
        coding_significant = coding & significant
        target_gene = str(target_row["GENE"])
        target_magma_ids = {target_gene}
        target_magma_ids.update(gene_symbol_to_ids.get(target_gene.upper(), set()))
        target_gene_mask = dat_source["GENE"].isin(target_magma_ids).to_numpy(bool)
        target_gene_valid = target_gene_mask & valid_p
        target_gene_significant = target_gene_mask & significant

        min_index = int(np.nanargmin(np.where(valid_p, p_value, np.nan))) if n_gene else None
        cis_n_gene = int(cis_valid.sum())
        cis_min_index = (
            int(np.nanargmin(np.where(cis_valid, p_value, np.nan))) if cis_n_gene else None
        )
        coding_n_gene = int(coding_valid.sum())
        coding_min_index = (
            int(np.nanargmin(np.where(coding_valid, p_value, np.nan))) if coding_n_gene else None
        )
        target_gene_n = int(target_gene_valid.sum())
        target_gene_min_index = (
            int(np.nanargmin(np.where(target_gene_valid, p_value, np.nan))) if target_gene_n else None
        )
        n_significant = int(significant.sum())
        n_cis_significant = int(cis_significant.sum())
        n_coding_significant = int(coding_significant.sum())
        n_target_gene_significant = int(target_gene_significant.sum())
        if args.protein_filter == "none":
            keep = n_gene > 0
            reason = "pass" if keep else "no_valid_gene_p"
        elif args.protein_filter == "any-bonf":
            keep = n_significant > 0
            reason = "pass" if keep else "no_bonferroni_significant_gene"
        elif args.protein_filter == "cis-bonf":
            keep = n_cis_significant > 0
            reason = "pass" if keep else "no_cis_bonferroni_significant_gene"
        elif args.protein_filter == "coding-bonf":
            keep = n_coding_significant > 0
            if keep:
                reason = "pass"
            elif coding_n_gene == 0:
                reason = "coding_gene_absent_from_magma_output"
            else:
                reason = "coding_gene_not_bonferroni_significant"
        else:
            keep = n_target_gene_significant > 0
            if keep:
                reason = "pass"
            elif target_gene_n == 0:
                reason = "target_gene_absent_from_magma_output"
            else:
                reason = "target_gene_not_bonferroni_significant"

        protein_qc_rows.append({
            "PROTEIN": protein,
            "GENE": target_row["GENE"],
            "CHR": target_chr_qc,
            "START": target_row["START"],
            "END": target_row["END"],
            "N_GENE_TESTED": n_gene,
            "BONF_THRESHOLD": bonf_threshold,
            "N_SIGNIFICANT_GENE": n_significant,
            "MIN_P": p_value[min_index] if min_index is not None else np.nan,
            "MIN_P_SOURCE": dat_source.iloc[min_index]["GENE"] if min_index is not None else "",
            "N_CIS_GENE": cis_n_gene,
            "N_CIS_SIGNIFICANT_GENE": n_cis_significant,
            "CIS_MIN_P": p_value[cis_min_index] if cis_min_index is not None else np.nan,
            "CIS_MIN_P_SOURCE": dat_source.iloc[cis_min_index]["GENE"] if cis_min_index is not None else "",
            "N_CODING_GENE": coding_n_gene,
            "N_CODING_SIGNIFICANT_GENE": n_coding_significant,
            "CODING_MIN_P": p_value[coding_min_index] if coding_min_index is not None else np.nan,
            "CODING_MIN_P_SOURCE": dat_source.iloc[coding_min_index]["GENE"] if coding_min_index is not None else "",
            "TARGET_MAGMA_IDS": ",".join(sorted(target_magma_ids)),
            "N_TARGET_GENE": target_gene_n,
            "N_TARGET_GENE_SIGNIFICANT": n_target_gene_significant,
            "TARGET_GENE_MIN_P": p_value[target_gene_min_index] if target_gene_min_index is not None else np.nan,
            "TARGET_GENE_MIN_P_SOURCE": dat_source.iloc[target_gene_min_index]["GENE"] if target_gene_min_index is not None else "",
            "PROTEIN_FILTER": args.protein_filter,
            "KEEP": bool(keep),
            "EXCLUSION_REASON": reason,
            "MAGMA_FILE": str(path),
        })

        if keep:
            filtered_pairs.append((protein, path))
            filtered_target_rows.append(target_row)
            novel = ~dat_source["GENE"].isin(seen_source)
            if novel.any():
                new_source = dat_source.loc[novel, ["GENE", "CHR", "START", "STOP"]].copy()
                source_parts.append(new_source)
                seen_source.update(new_source["GENE"].tolist())
        if (index + 1) % 250 == 0 or index + 1 == len(pairs):
            log(
                f"MAGMA protein QC: {index + 1}/{len(pairs)} "
                f"kept={len(filtered_pairs):,} union={len(seen_source):,}"
            )
    protein_qc = pd.DataFrame(protein_qc_rows)
    atomic_dataframe(protein_qc, out / "protein.qc.tsv.gz")
    if not filtered_pairs:
        raise SystemExit(
            f"ERROR: --protein-filter {args.protein_filter} excluded all {len(pairs)} matched proteins; "
            f"inspect {out / 'protein.qc.tsv.gz'}"
        )
    pairs = filtered_pairs
    target = pd.DataFrame(filtered_target_rows)
    target.insert(0, "TARGET_INDEX", np.arange(len(target), dtype=int))
    log(
        f"Protein filter {args.protein_filter}: kept {len(target):,}/{len(target_candidates):,}; "
        f"excluded {len(target_candidates) - len(target):,}"
    )

    # Local MAGMA outputs are almost, but not perfectly, identical gene grids.
    # The union is built only from proteins that passed the declared gate.
    source = pd.concat(source_parts, ignore_index=True)
    source.columns = ["SOURCE", "CHR", "START", "STOP"]
    source["SOURCE"] = source["SOURCE"].astype(str)
    source["CHR"] = source["CHR"].map(normalize_chr)
    source["START"] = pd.to_numeric(source["START"], errors="coerce")
    source["STOP"] = pd.to_numeric(source["STOP"], errors="coerce")
    source["POS"] = (source["START"] + source["STOP"]) / 2.0
    source["_CHR_ORDER"] = source["CHR"].map(chromosome_order)
    source = source.sort_values(["_CHR_ORDER", "START", "STOP", "SOURCE"], kind="mergesort")
    source = source.drop(columns="_CHR_ORDER").reset_index(drop=True)
    source.insert(0, "SOURCE_INDEX", np.arange(len(source), dtype=int))
    source_ids = source["SOURCE"].to_numpy(str)
    source_chr = source["CHR"].to_numpy(object)
    source_pos = source["POS"].to_numpy(float)

    tmp_matrix = out / f".score.{os.getpid()}.npy"
    matrix = np.lib.format.open_memmap(
        tmp_matrix, mode="w+", dtype="float32", shape=(len(source), len(target))
    )
    matrix[:] = np.nan
    qc_rows = []
    log(f"Build MAGMA graph matrix: {len(source):,} source genes x {len(target):,} protein children")
    for index, (protein, path) in enumerate(pairs):
        dat = magma_table(path)
        if args.score == "zstat" and "ZSTAT" not in dat.columns:
            raise SystemExit(f"ERROR: --score zstat requested but ZSTAT is absent: {path}")
        dat["GENE"] = dat["GENE"].astype(str)
        dat = dat.set_index("GENE", drop=False).reindex(source_ids)
        missing = int(dat["GENE"].isna().sum())
        if missing > max(100, int(0.05 * len(source))):
            raise SystemExit(f"ERROR: too many missing MAGMA genes ({missing}) in {path}")
        if args.score == "zstat":
            values = pd.to_numeric(dat["ZSTAT"], errors="coerce").to_numpy(float)
        else:
            p = pd.to_numeric(dat["P"], errors="coerce").to_numpy(float)
            values = -np.log10(np.clip(p, np.finfo(float).tiny, 1.0))
        if args.score_cap > 0:
            if args.score == "zstat":
                values = np.clip(values, -args.score_cap, args.score_cap)
            else:
                values = np.clip(values, 0.0, args.score_cap)
        target_row = target.iloc[index]
        cis = cis_rows(source_chr, source_pos, target_row["CHR"], target_row["POS"], args.cis_flank)
        center, scale = 0.0, 1.0
        if args.scale == "robust":
            center, scale = robust_center_scale(values, mask=~cis)
            values = (values - center) / scale
        matrix[:, index] = values.astype("float32")
        qc_rows.append({
            "PROTEIN": protein,
            "N_SOURCE": len(source),
            "N_MISSING": missing,
            "N_CIS": int(cis.sum()),
            "CENTER": center,
            "SCALE": scale,
            "MAX_SCORE": float(np.nanmax(values)),
        })
        if (index + 1) % 50 == 0 or index + 1 == len(pairs):
            log(f"MAGMA files: {index + 1}/{len(pairs)}")
    matrix.flush()
    del matrix
    os.replace(tmp_matrix, out / "score.npy")

    atomic_dataframe(source, out / "source.tsv.gz")
    atomic_dataframe(target, out / "target.tsv.gz")
    atomic_dataframe(pd.DataFrame(qc_rows), out / "prepare.qc.tsv.gz")
    atomic_json(out / "manifest.json", {
        "format": "phole-matrix-v1",
        "input_mode": "magma",
        "score": args.score,
        "score_semantics": (
            "MAGMA gene-level association Z statistic; it is not a signed protein-abundance effect"
            if args.score == "zstat"
            else "standardized -log10(MAGMA gene-level P); it is not a signed protein-abundance effect"
        ),
        "directional": False,
        "scaled": args.scale,
        "score_cap": args.score_cap,
        "cis_flank": args.cis_flank,
        "protein_filter": args.protein_filter,
        "protein_alpha": args.protein_alpha,
        "magma_build": magma_build,
        "protein_bed_build": bed_build,
        "n_source": len(source),
        "n_target_candidates": len(target_candidates),
        "n_target_excluded": len(target_candidates) - len(target),
        "n_target_gene_missing": int((protein_qc["N_TARGET_GENE"] == 0).sum()),
        "n_target_gene_bonferroni_significant": int(
            (protein_qc["N_TARGET_GENE_SIGNIFICANT"] > 0).sum()
        ),
        "n_target": len(target),
        "gwas_root": str(Path(args.gwas_root).resolve()),
        "category": args.category,
        "protein_bed": str(Path(args.protein_bed).resolve()),
        "protein_map": str(Path(args.protein_map).resolve()) if args.protein_map else None,
        "magma_gene_loc_files": sorted(str(path) for path in gene_loc_paths),
        "require_magma_done": require_magma_done,
        "source_grid": "union_of_selected_completed_magma_genes",
    })
    done.write_text("done\n", encoding="utf-8")
    log(f"PHOLE matrix done: {out / 'score.npy'}")


if __name__ == "__main__":
    main()

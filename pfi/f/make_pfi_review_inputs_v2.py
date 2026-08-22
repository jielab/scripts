#!/usr/bin/env python3
"""Robustly create blinded PFI review inputs.

This version is deliberately schema-tolerant. It does not assume a fixed text
column such as `model_text`; it detects the available text columns and uses only
columns that actually exist in the parquet/CSV inputs. This avoids DuckDB binder
errors when older/newer pipeline outputs use `blinded_text` rather than
`model_text`.
"""
from __future__ import annotations

import argparse
import json
import random
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception as e:  # pragma: no cover
    raise SystemExit(f"ERROR: pyarrow is required for this step: {e}")


def norm_col(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(x).lower())).strip("_")


def first_existing(cols: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {norm_col(c): c for c in cols}
    for c in candidates:
        if norm_col(c) in lookup:
            return lookup[norm_col(c)]
    return None


def dataset_for(path: str | Path) -> ds.Dataset:
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"ERROR: input does not exist: {p}")
    if p.is_dir():
        return ds.dataset(str(p), format="parquet", partitioning="hive")
    if p.suffix.lower() == ".parquet":
        return ds.dataset(str(p), format="parquet")
    raise SystemExit(f"ERROR: expected parquet path or directory: {p}")


def safe_read_small(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    if p.is_dir() or p.suffix.lower() == ".parquet":
        d = dataset_for(p)
        avail = list(d.schema.names)
        use = [c for c in (columns or avail) if c in avail]
        if not use:
            return pd.DataFrame()
        return d.to_table(columns=use).to_pandas()
    try:
        return pd.read_csv(p, dtype=str, low_memory=False)
    except Exception:
        return pd.DataFrame()


def make_article_id(row: dict) -> str:
    for k in ["article_id", "pmid", "pmcid", "doi"]:
        v = row.get(k)
        if v is not None and str(v).strip() not in {"", "nan", "None"}:
            return str(v).strip()
    return ""


def infer_domain(mesh: str, journal: str = "") -> str:
    x = f"{mesh} {journal}".lower()
    rules = [
        ("public_health", r"public health|epidemiol|population|environment|health policy|occupational|social medicine|global health"),
        ("clinical_science", r"clinical|trial|patient|therapy|surgery|oncology|cardiology|neurology|pediatr|psychiatry|radiology|emergency"),
        ("basic_medical_research", r"cell|molecular|gene|genom|protein|immunol|microbiol|biochem|physiol|animal|mouse|mice"),
        ("methods_data_software", r"method|statistics|bioinformatics|database|software|algorithm|machine learning|artificial intelligence"),
    ]
    for label, pat in rules:
        if re.search(pat, x):
            return label
    return "other_biomedical"


def publication_type_proxy(text: str) -> str:
    x = text.lower()
    if re.search(r"systematic review|meta[- ]analysis", x):
        return "review_meta_analysis"
    if re.search(r"randomi[sz]ed|clinical trial|trial registration|nct\d{8}", x):
        return "clinical_trial"
    if re.search(r"cohort|case[- ]control|cross[- ]sectional|registry|biobank", x):
        return "observational_study"
    if re.search(r"software|database|web server|algorithm|pipeline|package", x):
        return "methods_resource"
    return "research_article"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--blinded", required=True)
    ap.add_argument("--features", default="")
    ap.add_argument("--journal_rank", default="")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_xlsx", required=True)
    ap.add_argument("--key_out", required=True)
    ap.add_argument("--summary_out", required=True)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--min_text_chars", type=int, default=500)
    ap.add_argument("--max_scan_rows", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--high_anchor_csv", default="")
    ap.add_argument("--low_anchor_csv", default="")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    blinded_ds = dataset_for(args.blinded)
    bcols = list(blinded_ds.schema.names)

    id_cols = [c for c in [
        first_existing(bcols, ["article_id", "articleid", "id"]),
        first_existing(bcols, ["pmid", "pubmed_id"]),
        first_existing(bcols, ["pmcid", "pmc_id"]),
        first_existing(bcols, ["doi"]),
        first_existing(bcols, ["year", "publication_year", "pub_year"]),
        first_existing(bcols, ["journal_key"]),
        first_existing(bcols, ["journal", "journal_title", "source_title"]),
        first_existing(bcols, ["issn"]),
        first_existing(bcols, ["eissn"]),
        first_existing(bcols, ["mesh_terms", "mesh", "keywords"]),
        first_existing(bcols, ["has_pmc_full_text", "has_full_text", "full_text_available"]),
    ] if c]

    # Use only available text-bearing columns. Do not include nonexistent names in SQL/coalesce.
    text_candidates = [
        "model_text", "blinded_text", "full_text", "text", "article_text", "body_text",
        "title", "abstract", "intro_text", "methods_text", "results_text", "discussion_text",
        "data_availability_text", "software_availability_text", "ethics_text", "trial_registration_text",
    ]
    text_cols = []
    for c in text_candidates:
        hit = first_existing(bcols, [c])
        if hit and hit not in text_cols:
            text_cols.append(hit)
    if not text_cols:
        text_cols = [c for c in bcols if re.search(r"text|abstract|title|method|result|discussion", c, re.I)]
    if not text_cols:
        raise SystemExit(f"ERROR: no usable text column found in blinded dataset. Available columns: {bcols}")

    scan_cols = list(dict.fromkeys(id_cols + text_cols))
    scanner = blinded_ds.scanner(columns=scan_cols, batch_size=4096, use_threads=True)

    reservoir: list[dict] = []
    n_seen = 0
    n_scanned = 0
    full_text_col = first_existing(scan_cols, ["has_pmc_full_text", "has_full_text", "full_text_available"])
    next_report = 100000

    print(
        "INFO: sampling blinded review rows with "
        f"text_cols={text_cols}; full_text_col={full_text_col or 'none'}; "
        f"target_n={args.n}; max_scan_rows={args.max_scan_rows or 'unlimited'}",
        flush=True,
    )

    for batch in scanner.to_batches():
        df = batch.to_pandas()
        n_scanned += len(df)
        # Normalize column names inside selected records, but keep original while building text.
        text = df[text_cols].fillna("").astype(str).agg("\n".join, axis=1)
        ok = text.str.len() >= args.min_text_chars
        if full_text_col and full_text_col in df.columns:
            ft = df[full_text_col].astype(str).str.lower().isin(["1", "true", "yes", "y"])
            ok = ok & ft
        if not ok.any():
            continue
        sub = df.loc[ok].copy()
        sub_text = text.loc[ok]
        for idx, row in sub.iterrows():
            rec_raw = {norm_col(k): (None if pd.isna(v) else v) for k, v in row.to_dict().items()}
            tx = str(sub_text.loc[idx])
            rec = dict(rec_raw)
            rec["blinded_text"] = tx[:120000]
            rec["text_chars"] = len(tx)
            aid = make_article_id(rec)
            if not aid:
                continue
            rec["article_id"] = aid
            mesh = str(rec.get("mesh_terms", "") or "")
            journal = str(rec.get("journal", rec.get("journal_title", "")) or "")
            rec.setdefault("domain", infer_domain(mesh, journal))
            rec.setdefault("publication_type_proxy", publication_type_proxy(tx[:20000]))
            n_seen += 1
            if len(reservoir) < args.n:
                reservoir.append(rec)
            else:
                j = rng.randrange(n_seen)
                if j < args.n:
                    reservoir[j] = rec
        if n_scanned >= next_report:
            print(
                f"INFO: scanned {n_scanned:,} rows; eligible {n_seen:,}; selected {len(reservoir):,}",
                flush=True,
            )
            while next_report <= n_scanned:
                next_report *= 2
        if args.max_scan_rows > 0 and n_scanned >= args.max_scan_rows:
            print(
                f"INFO: reached review sampling scan cap ({args.max_scan_rows:,} rows); "
                f"eligible {n_seen:,}; selected {len(reservoir):,}",
                flush=True,
            )
            break

    if not reservoir:
        raise SystemExit(
            "ERROR: Could not build candidate pool. No rows passed the full-text/text-length filter. "
            f"Available columns={bcols}; text_cols={text_cols}; min_text_chars={args.min_text_chars}"
        )

    review = pd.DataFrame(reservoir)
    # Stable order for reproducible downstream files.
    review = review.sample(frac=1, random_state=args.seed).reset_index(drop=True)
    review.insert(0, "review_id", [f"PFI{args.seed:03d}_{i+1:05d}" for i in range(len(review))])

    # Ensure core columns exist.
    for c in ["article_id", "pmid", "pmcid", "doi", "year", "domain", "publication_type_proxy", "mesh_terms", "blinded_text", "text_chars"]:
        if c not in review.columns:
            review[c] = ""

    # Key file includes unblinded metadata if present. If not present, identifiers still allow later joins.
    key_cols = [c for c in [
        "review_id", "article_id", "pmid", "pmcid", "doi", "year", "journal_key", "journal", "journal_title",
        "issn", "eissn", "domain", "publication_type_proxy", "mesh_terms", "text_chars"
    ] if c in review.columns]
    key = review[key_cols].copy()

    # Blinded review file: keep identifiers needed for machine merge but no author/institution/journal if available.
    score_cols = [
        "new_results", "new_methods", "new_data", "new_software", "clinical_public_health_value",
        "evidence_strength", "reproducibility", "claim_evidence_gap", "paper_score", "review_notes"
    ]
    blinded_cols = ["review_id", "article_id", "pmid", "pmcid", "year", "domain", "publication_type_proxy", "mesh_terms", "text_chars", "blinded_text"]
    out = review[[c for c in blinded_cols if c in review.columns]].copy()
    for c in score_cols:
        out[c] = ""

    out_csv = Path(args.out_csv); out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_xlsx = Path(args.out_xlsx); out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    key_out = Path(args.key_out); key_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out = Path(args.summary_out); summary_out.parent.mkdir(parents=True, exist_ok=True)

    out.to_csv(out_csv, index=False)
    try:
        xlsx_out = out.copy()
        if "blinded_text" in xlsx_out.columns:
            xlsx_out["blinded_text"] = xlsx_out["blinded_text"].fillna("").astype(str).str.slice(0, 32000)
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
            xlsx_out.to_excel(writer, index=False, sheet_name="blinded_review")
    except Exception as e:
        print(f"WARNING: could not write XLSX ({e}); CSV was written successfully: {out_csv}")

    pq.write_table(pa.Table.from_pandas(key, preserve_index=False), key_out)

    summary = pd.DataFrame([{
        "n_scanned_rows": n_scanned,
        "n_eligible_rows": n_seen,
        "n_selected_rows": len(out),
        "seed": args.seed,
        "min_text_chars": args.min_text_chars,
        "available_blinded_columns": ";".join(bcols),
        "text_columns_used": ";".join(text_cols),
        "full_text_column_used": full_text_col or "",
        "note": "Schema-tolerant no-expert review input builder; journal/author metadata excluded from blinded CSV when possible.",
    }])
    summary.to_csv(summary_out, index=False)

    print(f"OK: wrote blinded review CSV: {out_csv} ({len(out):,} rows)")
    print(f"OK: wrote blinded review XLSX: {out_xlsx}")
    print(f"OK: wrote review key parquet: {key_out}")
    print(f"OK: wrote sampling summary: {summary_out}")


if __name__ == "__main__":
    main()

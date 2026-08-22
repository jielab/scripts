#!/usr/bin/env python3
"""Run a fine-tuned PubMedBERT PFI model over all blinded PubMed/PMC text."""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

def sql_string(value: object) -> str:
    return str(value).replace("'", "''")

def quote_ident(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'

def parquet_source(path: str) -> str:
    p = Path(path)
    if p.is_dir() or p.suffix != ".parquet":
        return sql_string((p / "**" / "*.parquet").as_posix())
    return sql_string(p.as_posix())

def open_parquet_dataset(path: str) -> ds.Dataset:
    p = Path(path)
    if p.is_dir():
        return ds.dataset(p, format="parquet", partitioning="hive")
    return ds.dataset(path, format="parquet")

PREDICTION_SCHEMA = pa.schema(
    [
        ("pmid", pa.string()),
        ("pmcid", pa.string()),
        ("doi", pa.string()),
        ("year", pa.int64()),
        ("has_pmc_full_text", pa.int64()),
        ("publication_types", pa.string()),
        ("mesh_terms", pa.string()),
        ("keywords", pa.string()),
        ("expected_journal_score_0_100", pa.float32()),
        ("predicted_paper_score", pa.float32()),
    ]
)

def prediction_table(rows: dict[str, object]) -> pa.Table:
    return pa.Table.from_pydict(rows, schema=PREDICTION_SCHEMA)

def parquet_row_count(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        return pq.ParquetFile(path).metadata.num_rows
    except Exception:
        return None

def duckdb_row_count(path: Path, where_sql: str | None = None) -> int | None:
    try:
        import duckdb
    except Exception:
        return None
    if not path.exists():
        return None
    q = f"SELECT count(*) FROM read_parquet('{sql_string(path.as_posix())}', union_by_name=true)"
    if where_sql:
        q += f" WHERE {where_sql}"
    con = duckdb.connect()
    try:
        return int(con.execute(q).fetchone()[0])
    except Exception:
        return None
    finally:
        con.close()

def build_sample_parquet(
    *,
    source_path: str,
    sample_path: Path,
    columns: list[str],
    sample_rows: int,
    seed: int,
    reuse: bool,
    full_text_only: bool,
) -> None:
    if sample_rows <= 0:
        return
    try:
        import duckdb
    except Exception as e:
        raise RuntimeError("DuckDB is required for random prediction sampling.") from e

    if reuse:
        existing_rows = parquet_row_count(sample_path)
        if existing_rows == sample_rows:
            if full_text_only:
                existing_full_text_rows = duckdb_row_count(sample_path, "has_pmc_full_text = 1")
                if existing_full_text_rows != existing_rows:
                    print(
                        f"Existing sample is not full-text-only ({existing_full_text_rows or 0:,}/{existing_rows:,}); "
                        f"rebuilding: {sample_path}",
                        flush=True,
                    )
                else:
                    print(f"Reusing existing full-text random sample: {sample_path} ({existing_rows:,} rows)", flush=True)
                    return
            else:
                print(f"Reusing existing random sample: {sample_path} ({existing_rows:,} rows)", flush=True)
                return
        if existing_rows is not None:
            print(
                f"Existing sample has {existing_rows:,} rows, expected {sample_rows:,}; rebuilding: {sample_path}",
                flush=True,
            )

    sample_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = Path(str(sample_path) + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    cols_sql = ", ".join(quote_ident(c) for c in columns)
    source_sql = parquet_source(source_path)
    tmp_sql = sql_string(tmp_path.as_posix())
    where_sql = "WHERE has_pmc_full_text = 1" if full_text_only else ""
    q = f"""
    COPY (
      WITH candidates AS (
        SELECT {cols_sql}
        FROM read_parquet('{source_sql}', union_by_name=true, hive_partitioning=1)
        {where_sql}
      )
      SELECT {cols_sql}
      FROM candidates
      USING SAMPLE reservoir({sample_rows} ROWS) REPEATABLE ({seed})
    )
    TO '{tmp_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """
    sample_label = "full-text random prediction sample" if full_text_only else "random prediction sample"
    print(f"Building {sample_label}: {sample_rows:,} rows -> {sample_path}", flush=True)
    con = duckdb.connect()
    try:
        con.execute("PRAGMA threads=8")
        con.execute(q)
    finally:
        con.close()
    tmp_path.replace(sample_path)
    actual_rows = parquet_row_count(sample_path)
    print(f"Wrote random prediction sample: {actual_rows or 0:,} rows -> {sample_path}", flush=True)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blinded", required=True)
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--max_length", type=int, default=512)
    ap.add_argument("--num_workers", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0, help="Debug limit; 0 means all rows")
    ap.add_argument("--sample_rows", type=int, default=0, help="Randomly sample this many input rows before prediction; 0 disables sampling")
    ap.add_argument("--sample_seed", type=int, default=1)
    ap.add_argument("--sample_path", default=None)
    ap.add_argument("--fresh_sample", action="store_true")
    ap.add_argument("--full_text_only", action="store_true", help="Restrict prediction candidates to rows with PMC full text")
    ap.add_argument("--progress_every", type=int, default=10000, help="Print progress every N predicted rows; 0 disables progress")
    ap.add_argument("--cpu_full_warning_rows", type=int, default=1_000_000)
    ap.add_argument("--allow_cpu_full", action="store_true")
    ap.add_argument("--fp16", action="store_true")
    args = ap.parse_args()

    try:
        import torch
        from transformers import AutoModelForSequenceClassification, AutoTokenizer
    except Exception as e:
        raise RuntimeError("Missing torch/transformers. Install transformer dependencies first.") from e

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    tmp_out = str(args.out) + ".tmp"
    if Path(tmp_out).exists():
        Path(tmp_out).unlink()

    dataset = open_parquet_dataset(args.blinded)
    columns = ["pmid", "pmcid", "doi", "year", "has_pmc_full_text", "publication_types", "mesh_terms", "keywords", "blinded_text"]
    available = [c for c in columns if c in dataset.schema.names]
    total_rows = dataset.count_rows()
    scan_filter = None
    candidate_rows = total_rows
    if args.full_text_only:
        if "has_pmc_full_text" not in dataset.schema.names:
            raise SystemExit("ERROR: --full_text_only requires a has_pmc_full_text column in the blinded text parquet.")
        scan_filter = ds.field("has_pmc_full_text") == 1
        candidate_rows = dataset.count_rows(filter=scan_filter)
        if candidate_rows <= 0:
            raise SystemExit("ERROR: no full-text rows found in blinded text parquet.")
    effective_sample_rows = args.sample_rows
    if args.limit and args.limit > 0 and effective_sample_rows > args.limit:
        effective_sample_rows = args.limit
    if effective_sample_rows > candidate_rows:
        effective_sample_rows = candidate_rows
    target_rows = effective_sample_rows if effective_sample_rows > 0 else candidate_rows
    if args.limit and args.limit > 0:
        target_rows = min(target_rows, args.limit)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Blinded table rows: {total_rows:,}", flush=True)
    if args.full_text_only:
        print(f"Full-text eligible rows: {candidate_rows:,}", flush=True)
    else:
        print(f"Prediction candidate rows: {candidate_rows:,}", flush=True)
    if effective_sample_rows > 0:
        print(f"Random sample rows: {effective_sample_rows:,}; seed: {args.sample_seed}", flush=True)
    if args.limit and args.limit > 0:
        print(f"Prediction limit: {args.limit:,}; target rows this run: {target_rows:,}", flush=True)
    print(f"Device: {device}", flush=True)
    if device.type == "cpu" and target_rows > args.cpu_full_warning_rows and not args.allow_cpu_full:
        raise SystemExit(
            "ERROR: refusing this PubMedBERT prediction on CPU because it is expected to be very slow. "
            f"Target rows: {target_rows:,}. Use a CUDA environment, set PFI_PRED_LIMIT for a test run, "
            "lower PFI_PRED_SAMPLE_ROWS, or set PFI_ALLOW_CPU_FULL_PRED=1 to run this CPU job intentionally."
        )

    prediction_input = args.blinded
    if effective_sample_rows > 0:
        sample_path = Path(args.sample_path) if args.sample_path else Path(str(args.out) + f".sample_{effective_sample_rows}_seed{args.sample_seed}.parquet")
        build_sample_parquet(
            source_path=args.blinded,
            sample_path=sample_path,
            columns=available,
            sample_rows=effective_sample_rows,
            seed=args.sample_seed,
            reuse=not args.fresh_sample,
            full_text_only=args.full_text_only,
        )
        prediction_input = str(sample_path)
        dataset = open_parquet_dataset(prediction_input)
        available = [c for c in columns if c in dataset.schema.names]
        scan_filter = None

    print(f"Loading model from {args.model_dir}", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(args.model_dir)
    model = AutoModelForSequenceClassification.from_pretrained(args.model_dir)
    model.to(device)
    model.eval()

    writer = None
    total = 0
    next_progress = args.progress_every if args.progress_every and args.progress_every > 0 else None

    def write_table(table: pa.Table):
        nonlocal writer
        if writer is None:
            writer = pq.ParquetWriter(tmp_out, table.schema, compression="zstd")
        writer.write_table(table)

    with torch.no_grad():
        for rb in dataset.to_batches(columns=available, batch_size=args.batch_size, filter=scan_filter):
            d = rb.to_pydict()
            texts = [(x or "") for x in d.get("blinded_text", [])]
            if not texts:
                continue
            enc = tokenizer(
                texts,
                truncation=True,
                padding=True,
                max_length=args.max_length,
                return_tensors="pt",
            )
            enc = {k: v.to(device) for k, v in enc.items()}
            autocast_enabled = bool(args.fp16 and device.type == "cuda")
            if hasattr(torch, "amp") and hasattr(torch.amp, "autocast"):
                autocast_context = torch.amp.autocast(device_type="cuda", enabled=autocast_enabled)
            else:
                autocast_context = torch.cuda.amp.autocast(enabled=autocast_enabled)
            with autocast_context:
                logits = model(**enc).logits
            pred = torch.clamp(logits, 0.0, 1.0).detach().cpu().numpy() * 100.0

            n = len(texts)
            rows = {
                "pmid": d.get("pmid", [None] * n),
                "pmcid": d.get("pmcid", [None] * n),
                "doi": d.get("doi", [None] * n),
                "year": d.get("year", [None] * n),
                "has_pmc_full_text": d.get("has_pmc_full_text", [None] * n),
                "publication_types": d.get("publication_types", [None] * n),
                "mesh_terms": d.get("mesh_terms", [None] * n),
                "keywords": d.get("keywords", [None] * n),
                "expected_journal_score_0_100": pred[:, 0].astype(np.float32),
                "predicted_paper_score": pred[:, 1].astype(np.float32),
            }
            write_table(prediction_table(rows))
            total += n
            if next_progress is not None and total >= next_progress:
                print(f"predicted {total:,} rows", flush=True)
                while next_progress <= total:
                    next_progress += args.progress_every
            if args.limit and total >= args.limit:
                break

    if writer is not None:
        writer.close()
        Path(tmp_out).replace(args.out)
    else:
        pq.write_table(prediction_table({name: [] for name in PREDICTION_SCHEMA.names}), args.out, compression="zstd")
    print(f"Wrote predictions for {total:,} rows -> {args.out}", flush=True)

if __name__ == "__main__":
    main()

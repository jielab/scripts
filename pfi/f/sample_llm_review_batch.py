#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import duckdb


# 🚩 Project utilities
import importlib.util as _module_util
import sys as _module_sys
from pathlib import Path as _ModulePath

_common_name = "pfi_common"
if _common_name not in _module_sys.modules:
	_spec = _module_util.spec_from_file_location(_common_name, _ModulePath(__file__).with_name("0.common.py"))
	_common = _module_util.module_from_spec(_spec)
	_module_sys.modules[_common_name] = _common
	_spec.loader.exec_module(_common)
else:
	_common = _module_sys.modules[_common_name]
parquet_expr = _common.parquet_expr


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--blinded", required=True)
	ap.add_argument("--features", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--n", type=int, default=1000)
	ap.add_argument("--seed", type=int, default=1)
	args = ap.parse_args()
	con = duckdb.connect()
	b = parquet_expr(args.blinded)
	df = con.execute(
		f"SELECT coalesce(CAST(article_id AS VARCHAR), CAST(pmid AS VARCHAR), CAST(pmcid AS VARCHAR)) AS article_id, year, coalesce(domain,'unknown') AS domain, coalesce(blinded_text, model_text, full_text, abstract, title, '') AS blinded_text FROM read_parquet('{b}', hive_partitioning=1) WHERE length(coalesce(blinded_text, model_text, full_text, abstract, title, ''))>500 USING SAMPLE reservoir({args.n} ROWS) REPEATABLE ({args.seed})"
	).df()
	con.close()
	Path(args.out).parent.mkdir(parents=True, exist_ok=True)
	with open(args.out, "w", encoding="utf-8") as f:
		for _, r in df.iterrows():
			f.write(json.dumps(r.to_dict(), ensure_ascii=False) + "\n")
	print(f"Wrote {len(df)} records to {args.out}")


if __name__ == "__main__":
	main()

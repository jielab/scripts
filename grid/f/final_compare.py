#!/usr/bin/env python3
"""FINAL comparison of verified current ABM, PRSformer and PRSformer2 runs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import pandas as pd
from final_outputs import aggregate, read_aggregate, current_abm_run


def workbook(d, path):
    from openpyxl import Workbook
    from openpyxl.utils.dataframe import dataframe_to_rows
    wb = Workbook()
    ws = wb.active
    ws.title = "performance"
    for row in dataframe_to_rows(d, index=False, header=True):
        ws.append(row)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    wb.save(path)


LABELS = {"PRSformer": "PRSformer", "PRSformer2": "PRSformer2 exploratory"}


def benchmark_current(run, root):
    from kinship import active, compatible
    path = run / "evaluation.json"
    return path.is_file() and compatible(json.loads(path.read_text()), active(root))


def current_tables(root, traits):
    from abm_core.model import file_sha256
    from abm_core.benchmark_results import performance
    from abm_core.prsformer_adapter import verify_complete
    from benchmark_paths import prsformer_directory, prsformer2_current
    root = Path(root)
    tables, status = [], []
    for trait in traits:
        found = False
        for pointer in [root / "abm/FINAL" / trait / "blocknet.CURRENT.json"]:
            run = current_abm_run(pointer)
            if run is None:
                continue
            d = read_aggregate(run / "main_performance.tsv")
            d["analysis_version"], d["run_id"], d["source"] = "FINAL", run.name, "ABM"
            d["ld_layout"] = json.loads((run / "RUN_MANIFEST.json").read_text())["contract"]["representation"]
            tables.append(d)
            q = pd.read_csv(run / "comparator_status.tsv", sep="\t")
            q["trait"] = trait
            status.extend(q.to_dict("records"))
            found = True
        if not found:
            status.append(dict(trait=trait, method="ABM", status="not_completed"))
    cache = root / "prsformer/cache"
    for trait in traits:
        run = prsformer_directory(root, trait) / "run"
        if benchmark_current(run, root):
            d = performance(run, cache, [trait], "PRSformer", exact_traits=True)
            d["source"], d["analysis_version"], d["run_id"] = "PRSformer", "FINAL", str(run.relative_to(root / "prsformer"))
            d["training_mode"] = "independent_trait"
            tables.append(d)
            status.append(dict(trait=trait, method="PRSformer", status="completed"))
        else:
            status.append(dict(trait=trait, method="PRSformer", status="not_completed"))
        run = prsformer2_current(root, trait)
        if run is not None and (run / "COMPLETE.json").is_file() and benchmark_current(run / "_LOCAL_ONLY", root):
            d = verify_complete(run, cache, [trait])
            d["source"], d["analysis_version"], d["run_id"] = "PRSformer2", "FINAL", str(run.relative_to(root / "prsformer2"))
            d["training_mode"] = "independent_trait"
            tables.append(d)
            status.append(dict(trait=trait, method="PRSformer2", status="completed"))
        else:
            status.append(dict(trait=trait, method="PRSformer2", status="not_completed"))
    from disco_publication import read_performance, METHOD
    for table in tables:
        table["method"] = table.method.replace(LABELS)
    for item in status:
        item["method"] = LABELS.get(item["method"], item["method"])
    for trait in traits:
        table, reason = read_performance(root, trait)
        if table is not None:
            table["analysis_version"], table["run_id"] = "FINAL", "disco_development_adapter"
            tables.append(table)
        status.append(dict(trait=trait, method=METHOD, status="completed" if table is not None else "unavailable", reason=reason))
    return tables, status


def run(root, traits):
    from abm_core.model import file_sha256
    tables, status = current_tables(root, traits)
    d = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=["trait", "method", "stage", "source"])
    aggregate(d)
    out = Path(root) / "final"
    out.mkdir(parents=True, exist_ok=True)
    sensitivity = []
    for trait in traits:
        for pointer in (Path(root) / "abm/FINAL" / trait).glob("*.CURRENT.json"):
            if pointer.name == "blocknet.CURRENT.json":
                continue
            directory = current_abm_run(pointer)
            if directory:
                table = read_aggregate(directory / "main_performance.tsv")
                table["run_id"], table["source"] = directory.name, "ABM sensitivity"
                sensitivity.append(table)
    (pd.concat(sensitivity, ignore_index=True) if sensitivity else pd.DataFrame(columns=["trait", "method", "source"])).to_csv(out / "sensitivity_performance.tsv", sep="\t", index=False)
    d.to_csv(out / "performance.tsv", sep="\t", index=False)
    workbook(d, out / "performance.xlsx")
    pd.DataFrame(status).to_csv(out / "performance_status.tsv", sep="\t", index=False)
    # Diagnostics retain their own targets and run identities, outside method performance.
    contrasts = []
    for trait in traits:
        for pointer in [Path(root) / "abm/FINAL" / trait / "blocknet.CURRENT.json"]:
            directory = current_abm_run(pointer)
            if directory is None:
                continue
            path = directory / "count_incremental_contrasts.tsv"
            if path.is_file():
                table = read_aggregate(path)
                table["trait"], table["run_id"] = trait, directory.name
                table["report_scope"] = "completed_run"
                contrasts.append(table)
    for marker in (Path(root) / "abm_reports").glob("*/COMPLETE.json"):
        folder = marker.parent
        receipt = json.loads(marker.read_text())
        if receipt.get("smoke"):
            continue
        if (receipt.get("status") != "completed" or
            file_sha256(folder / "REPORT_MANIFEST.json") != receipt["manifest_sha256"] or
            file_sha256(folder / "table_manifest.json") != receipt["table_manifest_sha256"]):
            raise ValueError("Derived report identity changed")
        provenance = json.loads((folder / "REPORT_MANIFEST.json").read_text())
        if provenance["trait"] not in traits:
            continue
        for item in json.loads((folder / "table_manifest.json").read_text()):
            if item["name"] == "count_incremental_contrasts":
                path = folder / item["path"]
                if file_sha256(path) != item["sha256"]:
                    raise ValueError("Derived count contrast changed")
                table = read_aggregate(path)
                table["report_scope"] = "derived_report"
                table["parent_run"] = provenance["parent"]
                table["parent_manifest_sha256"] = provenance["parent_manifest_sha256"]
                contrasts.append(table)
    (pd.concat(contrasts, ignore_index=True) if contrasts else
     pd.DataFrame(columns=["trait", "run_id", "objective", "contrast", "status"])).to_csv(
         out / "count_incremental_contrasts.tsv", sep="\t", index=False)
    history = []
    for marker in (Path(root) / "historical").glob("*/COMPLETE.json"):
        folder = marker.parent
        receipt = json.loads(marker.read_text())
        if receipt.get("status") != "completed" or file_sha256(folder / "IMPORT_MANIFEST.json") != receipt["manifest_sha256"]:
            raise ValueError("Historical import identity changed")
        for item in json.loads((folder / "table_manifest.json").read_text()):
            if item["name"] == "historical_performance":
                path = folder / item["path"]
                if file_sha256(path) != item["sha256"]:
                    raise ValueError("Historical aggregate changed")
                history.append(read_aggregate(path))
    (pd.concat(history, ignore_index=True) if history else pd.DataFrame(columns=["method", "status"])).to_csv(
        out / "historical_performance.tsv", sep="\t", index=False)
    print(json.dumps(dict(output=str(out / "performance.tsv"), rows=len(d), status=status)))


def main(argv=None):
    p = argparse.ArgumentParser(description="Build the FINAL comparison from completed current runs.")
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--traits", default="height,ldl,t2dm")
    args = p.parse_args(argv)
    traits = args.traits.split(",")
    if len(set(traits)) != len(traits) or not set(traits) <= {"height", "ldl", "t2dm"}:
        raise ValueError("Use unique traits from height,ldl,t2dm")
    run(args.analysis_root, traits)


if __name__ == "__main__":
    main()

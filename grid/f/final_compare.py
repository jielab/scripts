#!/usr/bin/env python3
"""FINAL comparison of verified current ABM, PRSformer and PRSformer2 runs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import pandas as pd


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


def current_tables(root, traits):
    from abm_core.model import file_sha256
    from abm_core.benchmark_results import performance
    from abm_core.prsformer_adapter import verify_complete
    root = Path(root)
    tables, status = [], []
    for trait in traits:
        found = False
        for pointer in sorted((root / "abm/FINAL" / trait).glob("*.CURRENT.json")):
            run = pointer.parent / json.loads(pointer.read_text())["run"]
            receipt = json.loads((run / "COMPLETE.json").read_text())
            if receipt.get("status") != "completed" or receipt.get("smoke"):
                raise ValueError("ABM current run is not a completed production run")
            if file_sha256(run / "RUN_MANIFEST.json") != receipt["manifest_sha256"]:
                raise ValueError("ABM run manifest changed")
            phases = json.loads((run / "PHASES.json").read_text())["phases"]
            for phase in ("evaluate", "report"):
                for name, expected in phases[phase].items():
                    if file_sha256(run / name) != expected:
                        raise ValueError("ABM report changed: " + name)
            d = pd.read_csv(run / "main_performance.tsv", sep="\t")
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
    run = root / "prsformer/run"
    if (run / "evaluation.json").is_file():
        d = performance(run, cache, traits, "PRSformer")
        d["source"], d["analysis_version"], d["run_id"] = "PRSformer", "FINAL", run.name
        tables.append(d)
        status.extend(dict(trait=t, method="PRSformer", status="completed") for t in traits)
    else:
        status.extend(dict(trait=t, method="PRSformer", status="not_completed") for t in traits)
    pointer = root / "prsformer2/CURRENT.json"
    if pointer.is_file():
        run = pointer.parent / json.loads(pointer.read_text())["run"]
        d = verify_complete(run, cache, traits)
        d["source"], d["analysis_version"], d["run_id"] = "PRSformer2", "FINAL", run.name
        tables.append(d)
        status.extend(dict(trait=t, method="PRSformer2", status="completed") for t in traits)
    else:
        status.extend(dict(trait=t, method="PRSformer2", status="not_completed") for t in traits)
    return tables, status


def run(root, traits):
    from abm_core.model import file_sha256
    tables, status = current_tables(root, traits)
    d = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=["trait", "method", "stage", "source"])
    out = Path(root) / "final"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "performance.tsv", sep="\t", index=False)
    workbook(d, out / "performance.xlsx")
    pd.DataFrame(status).to_csv(out / "performance_status.tsv", sep="\t", index=False)
    # Diagnostics retain their own targets and run identities, outside method performance.
    contrasts = []
    for trait in traits:
        for pointer in sorted((Path(root) / "abm/FINAL" / trait).glob("*.CURRENT.json")):
            directory = pointer.parent / json.loads(pointer.read_text())["run"]
            path = directory / "count_incremental_contrasts.tsv"
            if path.is_file():
                table = pd.read_csv(path, sep="\t")
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
                table = pd.read_csv(path, sep="\t")
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
                history.append(pd.read_csv(path, sep="\t"))
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

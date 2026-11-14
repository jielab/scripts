#!/usr/bin/env python3
"""Publish aggregate module results separately from retained runtime artifacts.

Callers validate their scientific completion receipts before publishing. Private
rosters, predictions, donor edges and training checkpoints stay under tmp/.
"""
from __future__ import annotations

from pathlib import Path
import json
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_layout import publish_result


def initialize(root):
    root = Path(root)
    for name in ("pca", "csx", "disco", "prsformer", "prepare", "twins",
                 "evaluate", "architecture", "final", "logs", "jobs", "tmp", "qc"):
        (root / name).mkdir(parents=True, exist_ok=True)


def publish_files(root, module, prefix, source, names):
    source, root = Path(source), Path(root)
    for name in names:
        publish_result(source / name, root / module / (prefix + name))


def publish_count_stage(root, trait, phase, folder):
    import pandas as pd
    root, folder = Path(root), Path(folder)
    initialize(root)
    stages = {"geometry": ("geometry",), "radius": ("prepare",),
              "radius_twins": ("prepare", "twins"), "prepare": ("geometry", "prepare"),
              "twins": ("twins",), "evaluate": ("evaluate",)}[phase]
    for stage in stages:
        receipt = folder / (stage.upper() + "_COMPLETE.json")
        publish_result(receipt, root / "qc/abm" / trait / receipt.name)
        if stage in {"geometry", "prepare"}:
            publish_files(root, "prepare", trait + ".", folder, ("features.tsv",))
        if stage == "prepare":
            model = json.loads((folder / "model.json").read_text())
            rows = [dict(trait=trait, role=role, n=value["n"], groups=value["groups"],
                         radius=model["radius"], task=model["task"])
                    for role, value in model["lineage"].items()]
            pd.DataFrame(rows).to_csv(folder / "cohort_summary.tsv", sep="\t", index=False)
            publish_files(root, "prepare", trait + ".", folder, ("cohort_summary.tsv",))
        if stage == "twins":
            rows = []
            for role in ("calibration", "test"):
                frame = pd.read_csv(folder / "_LOCAL_ONLY" / (role + ".predictions.tsv.gz"),
                                    sep="\t", usecols=["supported"])
                rows.append(dict(trait=trait, role=role, n=len(frame),
                                 supported_n=int(frame.supported.sum()),
                                 actual_coverage=float(frame.supported.mean())))
            pd.DataFrame(rows).to_csv(folder / "support_summary.tsv", sep="\t", index=False)
            publish_files(root, "twins", trait + ".", folder, ("support_summary.tsv",))
        if stage == "evaluate":
            publish_files(root, "evaluate", trait + ".", folder,
                          ("performance.tsv", "risk_coverage.tsv", "ancestry_audit.tsv",
                           "risk_coverage.png", "report.html"))


def publish_count_report(root, stage, folder, trait=None):
    folder, root = Path(folder), Path(root)
    prefix = (trait or "all") + "."
    receipt = folder / (stage.upper() + "_COMPLETE.json")
    record = json.loads(receipt.read_text())
    # The receipt explicitly enumerates aggregate report artifacts.
    status_names = {"status.tsv", "benchmark_comparison_status.tsv"}
    for name in record["artifacts"]:
        destination = root / ("qc" if name in status_names else stage)
        publish_result(folder / name, destination / (prefix + name))
    publish_result(receipt, root / "qc/abm" / (trait or "all") / receipt.name)


def publish_prsformer(root, trait, method, run, benchmark):
    root, run, benchmark = Path(root), Path(run), Path(benchmark)
    initialize(root)
    prefix = trait + (".prsformer2." if method == "prsformer2" else ".")
    publish_files(root, "prsformer", prefix, benchmark,
                  ("performance.tsv", "paper_style_metrics.tsv"))
    if method == "prsformer":
        publish_files(root, "prsformer", prefix, run.parent,
                      ("performance.png", "performance.xlsx", "training.png", "training.xlsx"))
    # The selected trained model is a result; optimizer checkpoints remain tmp.
    publish_result(run / "model.pt", root / "prsformer" / (prefix + "model.pt"), immutable=False)

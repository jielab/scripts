#!/usr/bin/env python3
"""Run the complementary checks for an existing split, collecting every outcome.

No Kraken library is built, no original assignment is rewritten, and failures
in one check do not discard evidence from the remaining checks. Threshold
diagnostics alone query the large index, serially with demand paging.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from common import (CODE_ROOT, atomic_json, code_identity, fingerprint, load_config, locked,
                    require_compute, utcnow)
from qc_paths import configured_split_root
from result_output import (managed_outputs, prepare as prepare_result,
                           finish as finish_result, replace_bool)
from storage_layout import metadata_dir


def sample_parts(value):
    parts = value.split("/")
    if len(parts) != 2 or any(not p or p in (".", "..") or "\\" in p for p in parts):
        raise ValueError("--sample must be DATASET/SAMPLE")
    return parts


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return default


def intact_result(folder):
    """Cheap completed-result check; no multi-terabyte input hashing."""
    folder = Path(folder)
    record = read_json(folder / "complete.json", {})
    files = record.get("files", {})
    if record.get("status") != "complete" or not files:
        return False
    for name, size in files.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            return False
        file = folder / relative
        if not file.is_file() or file.is_symlink() or file.stat().st_size != size:
            return False
    # Competitive validation records only the references it actually used.
    # Inspect metadata, not FASTA/index contents, before reusing that evidence.
    for item in record.get("input_states", []):
        try:
            state = Path(item["path"]).stat()
        except (OSError, KeyError):
            return False
        for key, value in (("size", state.st_size), ("mtime_ns", state.st_mtime_ns), ("ctime_ns", state.st_ctime_ns)):
            if key in item and item[key] != value:
                return False
    return True


def request_input_states(cfg, meta):
    """Small source receipts and cheap stats of relevant external inputs."""
    paths = [meta / name for name in ("complete.json", "summary.json", "provenance.json", "counts.tsv")]
    paths += [Path(cfg["metadata_dir"]) / "reference_taxids.tsv",
              Path(cfg["library_dir"]) / "plan.json",
              Path(cfg["db_dir"]) / "database.json"]
    paths += list(Path(cfg["taxonomy_dir"]).glob("*.dmp"))
    paths += list(Path(cfg["metadata_dir"]).joinpath("summaries").glob("*assembly_summary*"))
    human = Path(cfg["human_reference_dir"])
    paths += list(human.glob("GCA_000001405.15_GRCh38_no_alt_analysis_set.fa*"))
    paths += [human / "BUILD.COMPLETE"]
    sample = (read_json(meta / "provenance.json", {}).get("context") or {}).get("sample", {})
    paths += [Path(sample[k]) for k in ("r1", "r2") if sample.get(k)]
    paths += [Path(cfg["tools_dir"]) / "minimap2/2.28/minimap2",
              Path(cfg["tools_dir"]) / "bwa-mem2/2.2.1/bwa-mem2-launch"]
    result = {}
    for path in sorted(set(paths)):
        try:
            stat = path.stat()
            result[str(path)] = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns)
        except FileNotFoundError:
            result[str(path)] = {"exists": False}
    return result


def execute_step(step, *, logs, reuse=False, runner=subprocess.run):
    """Keep an explicit failed/blocked result, then allow the caller to continue."""
    result = {"step": step["name"], "status": "pending", "exit_code": None,
              "seconds": 0, "output": str(step["out"]),
              "log": str(logs / (step["name"] + ".log")), "note": "",
              "command": step["command"]}
    if step.get("disabled"):
        result.update(status="not_requested", note=step["disabled"])
        return result
    if step.get("blocked"):
        result.update(status="blocked", note=step["blocked"])
        return result
    if reuse and intact_result(step["out"]):
        result.update(status="reused", exit_code=0,
                      note="Same recorded suite request; complete files reused")
        return result
    started = time.monotonic()
    logs.mkdir(parents=True, exist_ok=True)
    try:
        with Path(result["log"]).open("w") as handle:
            completed = runner(step["command"], stdout=handle, stderr=subprocess.STDOUT)
        result["exit_code"] = completed.returncode
        result["status"] = "complete" if completed.returncode == 0 else "failed"
        if completed.returncode:
            lines = Path(result["log"]).read_text(errors="replace").splitlines()
            result["note"] = " | ".join(lines[-4:])[-1600:]
    except OSError as exc:
        result.update(status="failed", note=f"{type(exc).__name__}: {exc}")
    result["seconds"] = round(time.monotonic() - started, 3)
    return result


def _link(base, destination, label):
    relative = os.path.relpath(destination, base)
    return '<a href="' + html.escape(relative, quote=True) + '">' + html.escape(label) + '</a>'


def render_report(out, report):
    """Small standalone, readable report; tables remain exact TSV/JSON data."""
    sample = report["sample_id"]
    summary = report.get("split_summary") or {}
    total = int(summary.get("fragments", 0))
    statuses = summary.get("status", summary.get("statuses", {}))
    if not isinstance(statuses, dict):
        statuses = {}
    rows = []
    for name, value in (report.get("status_fragments") or statuses).items():
        count = int(value)
        rows.append(f"<tr><td>{html.escape(name)}</td><td>{count:,}</td><td>{count / total:.6%}</td></tr>" if total else "")
    stage_rows = []
    for item in report["steps"]:
        path = Path(item["output"])
        preferred = next((path / f for f in ("human-validation.tsv", "human-competitive-validation.json",
                         "plant-validation.json", "threshold-summary.tsv", "sample-qc.tsv", "report.html",
                         "reference-audit.json") if (path / f).is_file()), path)
        stage_rows.append("<tr>" + "".join([
            "<td>" + html.escape(item["step"]) + "</td>",
            '<td class="' + html.escape(item["status"]) + '">' + html.escape(item["status"]) + "</td>",
            "<td>" + _link(out, preferred, "results") + "</td>",
            "<td>" + (_link(out, item["log"], "log") if Path(item["log"]).is_file() else "—") + "</td>",
            "<td>" + html.escape(item["note"]) + "</td>",
        ]) + "</tr>")
    text = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>Metag validation — {html.escape(sample)}</title>
<style>body{{font:16px/1.6 system-ui,sans-serif;max-width:1250px;margin:38px auto;padding:0 24px;color:#193043;background:#f8fafc}}h1,h2{{line-height:1.25}}table{{border-collapse:collapse;width:100%;background:white;margin:18px 0 30px}}td,th{{text-align:left;padding:10px 12px;border-bottom:1px solid #dce5ec;vertical-align:top}}th{{background:#e8eff5}}td:last-child{{overflow-wrap:anywhere}}.failed{{color:#aa2632;font-weight:700}}.complete,.reused{{color:#1f6953}}.partially_tested{{color:#996000;font-weight:700}}code{{background:#e8eff5;padding:2px 5px}}a{{color:#075d8c}}.note{{padding:16px;background:#fff3d4;border-left:4px solid #d6971e}}</style>
<h1>Metag validation: {html.escape(sample)}</h1>
<p>Status: <strong>{html.escape(report['status'])}</strong> · {total:,} paired fragments in the saved split.</p>
<p>One fragment is one read pair for paired-end data. File size is not a classification fraction.</p>
<h2>Saved classification</h2><table><tr><th>Status</th><th>Fragments</th><th>Fraction</th></tr>{''.join(rows)}</table>
<h2>Independent checks and diagnostics</h2><table><tr><th>Check</th><th>Execution / coverage</th><th>Output</th><th>Log</th><th>Explanation</th></tr>{''.join(stage_rows)}</table>
<div class="note">Alignment support is conditional on the reference panel. Human-only GRCh38 matches mean human-compatible. Competitive panels strengthen the evidence but do not prove species identity against omitted references. Plant reads do not by themselves prove food intake. Threshold diagnostics measure label sensitivity, not accuracy.</div>
<h2>How to read the outputs</h2><ul>
<li><b>human_candidates:</b> support among the human pairs already extracted by Kraken.</li>
<li><b>human_full:</b> search all original pairs; the crosstab shows compatible reads missed by the human bin.</li>
<li><b>human_competitive / plants:</b> paired support, competing matches, ambiguous/unmapped reads, loci and omitted references. MAPQ is within the declared panel.</li>
<li><b>reads:</b> Q30, N content, complexity and adapter signals by saved classification status. The QC scope is recorded.</li>
<li><b>thresholds:</b> the same seeded random sample across all input ordinals, compared at c=0.1/g=2, c=0/g=2 and c=0/g=1.</li>
</ul><p>Machine-readable status: <a href="suite-report.json">suite-report.json</a> and <a href="steps.tsv">steps.tsv</a>.</p>
</html>"""
    (out / "report.html").write_text(text, encoding="utf-8")


def write_outputs(out, report):
    atomic_json(out / "suite-report.json", report)
    fields = ["step", "status", "exit_code", "seconds", "output", "log", "note"]
    with (out / "steps.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report["steps"])
    render_report(out, report)


def annotate_coverage(result):
    paths = [Path(result["output"]) / "plant-validation.json",
             Path(result["output"]) / "human-competitive-validation.json"]
    record = next((read_json(path) for path in paths if path.is_file()), None)
    if record is not None and result["status"] in ("complete", "reused"):
        missing = int(record.get("untested_fragments", 0))
        if missing or record.get("status") == "partially_tested":
            result.update(status="partially_tested",
                          note=f"Execution finished; {missing:,} candidate fragments not tested. See panel/reference omissions.")


@managed_outputs
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--root")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out")
    parser.add_argument("--replace", type=replace_bool, default=False,
                        help="TRUE recomputes checks; FALSE resumes the same completed checks [FALSE]")
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--qc-fragments", type=int, default=0,
                        help="0 checks every original pair [0]; a positive limit is explicitly a prefix")
    parser.add_argument("--threshold-fragments", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--max-reference-gib", type=float, default=24)
    parser.add_argument("--max-species", type=int, default=60)
    parser.add_argument("--candidates-only", action="store_true",
                        help="Only saved candidate checks; omit all raw-input scans and threshold queries")
    parser.add_argument("--without-thresholds", action="store_true",
                        help="Run all independent checks without opening the large Kraken index")
    parser.add_argument("--allow-local", action="store_true")
    args = parser.parse_args(argv)
    require_compute(args.allow_local)
    if args.threads < 1 or args.qc_fragments < 0 or args.threshold_fragments < 1 or args.max_species < 1 or args.max_reference_gib <= 0:
        parser.error("Invalid thread, fragment or reference limits")
    parts = sample_parts(args.sample)
    cfg = load_config(args.config, args.root)
    folder = configured_split_root(cfg).joinpath(*parts)
    meta = metadata_dir(folder)
    summary = read_json(meta / "summary.json")
    if not (meta / "complete.json").is_file() or not isinstance(summary, dict):
        parser.error("--all needs an existing committed split (qc/complete.json and qc/summary.json)")
    if summary.get("sample_id") not in (None, args.sample):
        parser.error("Requested sample differs from saved summary")
    qc = Path(cfg["work_dir"]) / "qc"
    out = Path(args.out).expanduser().absolute() if args.out else qc.joinpath("validate", *parts)
    old = read_json(out / "suite-report.json", {})
    request = {"version": "20261114", "sample": args.sample,
               "split_identity": read_json(meta / "complete.json", {}).get("identity"),
               "split_input_digest": summary.get("input_digest"),
               "parameters": {k: v for k, v in vars(args).items() if k not in ("config", "root", "out", "replace", "allow_local")},
               "effective_config": {k: v for k, v in cfg.items() if k not in ("config_path", "code_root")},
               "source_identity": code_identity(),
               "input_states": request_input_states(cfg, meta),
               "paths": {k: cfg[k] for k in ("root", "reference_raw", "metadata_dir", "library_dir", "output_dir", "tools_dir", "human_reference_dir")}}
    request_id = fingerprint(request)  # Small JSON identity only, no FASTQ/index read.
    # Recreate only this aggregate report. The individual stages live in their
    # existing canonical directories and can be resumed below.
    prepare_result(out, replace=True, cfg=cfg, protected=[folder])
    snapshot = out / "effective-config.json"
    atomic_json(snapshot, cfg)
    base = ["--config", str(snapshot), "--sample", args.sample]
    local = ["--allow-local"] if args.allow_local else []
    stages = []

    def add(name, module, destination, options=(), disabled="", has_sample=True):
        command = [sys.executable, str(CODE_ROOT / "f" / (module + ".py"))]
        command += base if has_sample else ["--config", str(snapshot)]
        command += ["--out", str(destination), *map(str, options)]
        if module not in ("split_summary", "reference_audit"):
            command += local
        if module != "reference_audit":
            command += ["--replace", "TRUE"]
        stages.append(dict(name=name, command=command, out=destination, disabled=disabled))

    add("report", "split_summary", qc.joinpath("split", *parts), ["--require-complete"])
    add("references", "reference_audit", qc / "reference-audit", has_sample=False)
    no_raw = "--candidates-only: raw-input analysis not requested" if args.candidates_only else ""
    add("human_candidates", "human_validate", qc.joinpath("human", *parts, "candidates"),
        ["--scope", "candidates", "--threads", args.threads])
    panel_args = ["--threads", args.threads, "--max-reference-gib", args.max_reference_gib,
                  "--max-species", args.max_species]
    add("human_competitive", "plant_validate", qc.joinpath("human-competitive", *parts),
        ["--group", "human", *panel_args])
    add("plants", "plant_validate", qc.joinpath("plants", *parts),
        ["--include-broad", *panel_args])
    add("reads", "sample_diagnose", qc.joinpath("diagnose", *parts),
        ["--limit", args.qc_fragments], no_raw)
    add("human_full", "human_validate", qc.joinpath("human", *parts, "full"),
        ["--scope", "full", "--threads", args.threads], no_raw)
    add("thresholds", "threshold_diagnose", qc.joinpath("threshold", *parts),
        ["--fragments", args.threshold_fragments, "--seed", args.seed, "--threads", args.threads],
        no_raw or ("--without-thresholds: no global index query requested" if args.without_thresholds else ""))
    status_counts = {}
    with (meta / "counts.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            status_counts[row["status"]] = status_counts.get(row["status"], 0) + int(row["fragments"])
    if sum(status_counts.values()) != int(summary["fragments"]):
        raise RuntimeError("Saved counts do not conserve split fragments")
    report = {"schema": "metag-validation-suite-v1", "status": "running",
              "sample_id": args.sample, "started_utc": utcnow(), "request": request,
              "request_identity": request_id, "split_summary": summary,
              "status_fragments": status_counts, "steps": [],
              "meaning": "Execution status and reference coverage are reported separately; candidate evidence is not biological ground truth."}
    old_steps = {s["step"]: s for s in old.get("steps", [])}
    same_request = old.get("request_identity") == request_id
    for step in stages:
        print(f"VALIDATE {args.sample}: {step['name']}", flush=True)
        if step["name"] in ("plants", "human_competitive") and any(
                s["step"] == "references" and s["status"] == "failed" for s in report["steps"]):
            step["blocked"] = "Reference metadata audit failed; stale competitive panels were not used. See references log."
        previous = old_steps.get(step["name"], {})
        reuse = (not args.replace and same_request and previous.get("status") in ("complete", "reused")
                 and step["name"] not in ("report", "references"))
        if step["name"] == "references":
            # Multiple sample jobs may share this small metadata report.
            with locked(qc / ".reference-audit.lock"):
                result = execute_step(step, logs=out / "logs", reuse=False)
        else:
            result = execute_step(step, logs=out / "logs", reuse=reuse)
        annotate_coverage(result)
        report["steps"].append(result)
        write_outputs(out, report)
    failed = [s["step"] for s in report["steps"] if s["status"] in ("failed", "blocked")]
    partial = [s["step"] for s in report["steps"] if s["status"] == "partially_tested"]
    report.update(status="incomplete" if failed else "partially_tested" if partial else "complete",
                  failed_steps=failed, partially_tested_steps=partial, completed_utc=utcnow())
    write_outputs(out, report)
    if not failed:
        finish_result(out)
    print(f"Validation {report['status']}: {out / 'report.html'}", flush=True)
    if failed:
        print("Failed checks (other checks were still attempted): " + ", ".join(failed), file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()

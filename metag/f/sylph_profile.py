#!/usr/bin/env python3
"""Optional genome-level species profiling with sylph, using metag representatives.

This module is intentionally separate from Kraken read assignment.  Sylph is used
as a genome-level containment/ANI profiler to discover plausible species; its
output can then seed ``species_refine.py`` for competitive read-level mapping.
It never directly writes species FASTQs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess

from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, load_config, require_compute


def _exe(value):
    resolved = shutil.which(value) if os.sep not in value else str(Path(value).resolve(strict=True))
    if not resolved:
        raise RuntimeError(f"Cannot find sylph executable: {value}")
    return resolved


def _reps(path, reference_raw):
    rows = []
    taxids, paths = set(), set()
    with Path(path).open(newline="", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        required = {"species_taxid", "group", "path", "name"}
        if not r.fieldnames or not required <= set(r.fieldnames):
            raise RuntimeError(f"Representative table missing {required}")
        for row in r:
            p = Path(row["path"])
            p = p if p.is_absolute() else Path(reference_raw) / p
            p = p.resolve(strict=True)
            tid = int(row["species_taxid"])
            if tid <= 0 or tid in taxids or str(p) in paths or not p.is_file():
                raise RuntimeError(f"Invalid/duplicate species representative: taxid={tid}, path={p}")
            taxids.add(tid); paths.add(str(p))
            rows.append({**row, "species_taxid": tid, "resolved_path": str(p)})
    return rows


def plan_db(cfg, reps_path, out, threads):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    rows = _reps(reps_path, cfg["reference_raw"])
    groups = {"prok": [], "virus": [], "euk": []}
    for row in rows:
        group = row["group"]
        bucket = "prok" if group in {"bacteria", "archaea"} else ("virus" if group == "virus" else "euk")
        groups[bucket].append(row)
    manifest = []
    for bucket, values in groups.items():
        path = out / f"{bucket}.genomes.txt"
        path.write_text("".join(row["resolved_path"] + "\n" for row in values), encoding="utf-8")
        manifest.append({"bucket": bucket, "genomes": len(values), "list": str(path)})
    script = out / "build-sylph-db.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", 'SYLPH="${SYLPH:-sylph}"', f'THREADS="${{THREADS:-{threads}}}"',
             f'cd -- {shlex.quote(str(out))}']
    if groups["prok"]:
        lines += [f'$SYLPH sketch -l {shlex.quote(str(out / "prok.genomes.txt"))} -t "$THREADS" -c 200 -o prok',
                  '$SYLPH convert-db-two-screen prok.syldb -o prok_v2']
    if groups["virus"]:
        lines += [f'$SYLPH sketch -l {shlex.quote(str(out / "virus.genomes.txt"))} -t "$THREADS" -c 100 -o virus']
    if groups["euk"]:
        lines += [f'$SYLPH sketch -l {shlex.quote(str(out / "euk.genomes.txt"))} -t "$THREADS" -c 200 -o euk']
    script.write_text("\n".join(lines) + "\n", encoding="utf-8"); script.chmod(0o755)
    report = {"schema": "metag-sylph-plan-v1", "representatives": len(rows), "groups": manifest,
              "build_script": str(script),
              "note": "prok uses syl2db two-stage conversion; viruses remain c100 syldb; euk remains c200 syldb"}
    atomic_json(out / "sylph-plan.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))


def profile(cfg, reps_path, dbs, r1, r2, out, threads, sylph, min_ani):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    exe = _exe(sylph)
    raw = out / "sylph-profile.tsv"
    command = [exe, "profile", "-d", *[str(Path(db).resolve(strict=True)) for db in dbs]]
    command += ["-1", str(Path(r1).resolve(strict=True)), "-2", str(Path(r2).resolve(strict=True)),
                "-t", str(threads), "--minimum-ani", str(min_ani), "-o", str(raw)]
    log = out / "sylph.log"
    with log.open("w", encoding="utf-8") as h:
        h.write("COMMAND " + json.dumps(command) + "\n"); h.flush()
        run = subprocess.run(command, stdout=h, stderr=subprocess.STDOUT)
    if run.returncode:
        raise RuntimeError(f"sylph exited {run.returncode}; inspect {log}")
    reps = _reps(reps_path, cfg["reference_raw"])
    by_path = {str(Path(row["resolved_path"]).resolve()): row for row in reps}
    by_basename = {}
    for row in reps:
        by_basename.setdefault(Path(row["resolved_path"]).name, []).append(row)
    candidates = []
    filtered_rows = invalid_ani_rows = 0
    with raw.open(newline="", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        if not reader.fieldnames or "Genome_file" not in reader.fieldnames or "Adjusted_ANI" not in reader.fieldnames:
            raise RuntimeError(f"Unexpected sylph output columns: {reader.fieldnames}")
        for row in reader:
            try:
                ani = float(row["Adjusted_ANI"])
            except (TypeError, ValueError):
                invalid_ani_rows += 1
                continue
            if not math.isfinite(ani) or not 0 <= ani <= 100:
                invalid_ani_rows += 1
                continue
            if ani < min_ani:
                filtered_rows += 1
                continue
            genome = str(Path(row["Genome_file"]).resolve())
            ref = by_path.get(genome)
            if ref is None:
                # sylph may preserve the literal list path rather than normalized path.
                matches = by_basename.get(Path(row["Genome_file"]).name, [])
                if len(matches) > 1:
                    raise RuntimeError(f"Ambiguous sylph genome basename: {row['Genome_file']}")
                ref = matches[0] if matches else None
            if ref is None:
                raise RuntimeError(f"Qualified sylph genome has no representative mapping: {row['Genome_file']}; "
                                   "supply the representative table used to build this database")
            candidates.append({"species_taxid": ref["species_taxid"], "name": ref["name"], "group": ref["group"],
                               "adjusted_ani": ani,
                               "taxonomic_abundance": row.get("Taxonomic_abundance", ""),
                               "sequence_abundance": row.get("Sequence_abundance", ""),
                               "effective_coverage": row.get("Eff_cov", ""),
                               "genome_file": row["Genome_file"]})
    candidates.sort(key=lambda r: (-float(r["adjusted_ani"]), r["species_taxid"]))
    cand = out / "sylph-candidates.tsv"
    fields = ["species_taxid", "name", "group", "adjusted_ani", "taxonomic_abundance", "sequence_abundance", "effective_coverage", "genome_file"]
    with cand.open("w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(candidates)
    report = {"schema": "metag-sylph-profile-v1", "command": command, "min_ani": min_ani,
              "candidates": len(candidates), "candidate_table": str(cand), "raw_profile": str(raw),
              "rows_below_minimum_ani": filtered_rows, "rows_with_invalid_ani": invalid_ani_rows,
              "interpretation": "Genome-level candidate discovery; ANI and abundance do not assign individual reads or validate refinement accuracy"}
    atomic_json(out / "sylph-summary.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))


@managed_outputs
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config"); p.add_argument("--root")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("plan-db")
    a.add_argument("--config", default=argparse.SUPPRESS)
    a.add_argument("--root", default=argparse.SUPPRESS)
    a.add_argument("--representatives", required=True); a.add_argument("--out", required=True); a.add_argument("--threads", type=int, default=50)
    a.add_argument("--allow-local", action="store_true")
    b = sub.add_parser("profile")
    b.add_argument("--config", default=argparse.SUPPRESS)
    b.add_argument("--root", default=argparse.SUPPRESS)
    b.add_argument("--representatives", required=True); b.add_argument("--db", action="append", required=True)
    b.add_argument("--r1", required=True); b.add_argument("--r2", required=True); b.add_argument("--out", required=True)
    b.add_argument("--threads", type=int, default=32); b.add_argument("--sylph"); b.add_argument("--min-ani", type=float, default=95.0)
    b.add_argument("--replace", type=replace_bool, default=True)
    b.add_argument("--allow-local", action="store_true")
    args = p.parse_args(argv)
    cfg = load_config(args.config, args.root)
    if args.threads < 1:
        p.error("--threads must be positive")
    if args.command == "plan-db":
        # Planning reads metadata and creates lists only; it is allowed on the submit host.
        plan_db(cfg, args.representatives, args.out, args.threads)
    else:
        require_compute(args.allow_local)
        if not math.isfinite(args.min_ani) or not 0 <= args.min_ani <= 100:
            p.error("--min-ani must be finite and in 0..100")
        args.sylph = args.sylph or str(Path(cfg["tools_dir"]) / "sylph/1.0.0/sylph")
        if not prepare_result(args.out, replace=args.replace, cfg=cfg, protected=[args.representatives, args.r1, args.r2, *args.db]):
            return
        profile(cfg, args.representatives, args.db, args.r1, args.r2, args.out, args.threads, args.sylph, args.min_ani)
        finish_result(args.out)


if __name__ == "__main__":
    main()

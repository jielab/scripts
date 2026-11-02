#!/usr/bin/env python3
"""Plan a balanced HERV census subset and audit completed individual censuses.

The existing ``1.herv.py cohort burden`` remains the scalable all-sample
genotype analysis. This module never runs RepeatMasker or changes a completed
census. Percentages describe reference-conditioned sequence models, not
independent assemblies or the prevalence of polymorphic HERV insertions.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import re
import shlex
import statistics
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.0.0"
CATEGORIES = ("ERV_without_MaLR", "ERV_with_MaLR", "all_LTR")
COUNT_FIELDS = ("union_bp", "genome_bp", "union_non_N_bp", "genome_non_N_bp")
PERCENT_FIELDS = ("percent_genome", "percent_non_N")
LABELS = ("diploid", "reference_callable_diploid")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if r.get(k) is None else r.get(k, "") for k in fields} for r in rows)


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def metadata(path, group_column="auto", sex_column="auto"):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        lines = [line for line in f if line.strip() and not line.startswith("##")]
    if not lines:
        raise ValueError("Empty sample metadata")
    values = list(csv.reader(lines, delimiter="\t")) if "\t" in lines[0] else [line.split() for line in lines]
    header = [h.lower().lstrip("#").strip() for h in values[0]]
    if len(set(header)) != len(header):
        raise ValueError("Duplicate metadata header")
    def column(request, candidates):
        names = candidates if request == "auto" else (request.lower(),)
        return next((name for name in names if name in header), None)
    sample_column = column("auto", ("sample_id", "sample", "iid", "id"))
    group_column = column(group_column, ("super_pop", "superpop", "super_population", "population", "pop"))
    sex_column = column(sex_column, ("sex", "gender"))
    if not sample_column or not group_column or not sex_column:
        raise ValueError("Metadata requires sample/sample_id/IID, a population column, and sex/gender; auto prefers super_pop")
    if any(len(row) != len(header) for row in values[1:]):
        raise ValueError("Metadata row does not match its header")
    rows = [dict(zip(header, row)) for row in values[1:]]
    seen, result = set(), []
    for row in rows:
        sample = row[sample_column].strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"Invalid sample_id: {sample!r}")
        if sample in seen:
            raise ValueError(f"Duplicate sample_id: {sample}")
        seen.add(sample)
        group, sex = row[group_column].strip(), row[sex_column].strip()
        if not group or not sex or group in {"NA", "."} or sex in {"NA", "."}:
            raise ValueError(f"Missing population/sex metadata for {sample}")
        include = row.get("include", "1").strip()
        unrelated = row.get("unrelated", "unknown").strip()
        if include not in {"0", "1"} or unrelated not in {"0", "1", "unknown"}:
            raise ValueError("include must be 0/1; unrelated must be 0/1/unknown")
        result.append(dict(sample_id=sample, population=group, sex=sex,
                           include=include, unrelated=unrelated))
    if not result:
        raise ValueError("Empty sample metadata")
    return result


def balanced_subset(rows, total=20, per_stratum=2, seed=20261102, always=()):
    if total < 1 or per_stratum < 1:
        raise ValueError("total and per-stratum must be positive")
    eligible = {r["sample_id"]: r for r in rows if r["include"] == "1" and r["unrelated"] != "0"}
    if len(set(always)) != len(always) or set(always) - eligible.keys():
        raise ValueError("Explicit included samples must be unique and eligible in metadata")
    if len(always) > total:
        raise ValueError("Explicit included samples exceed total")
    # Sorting first makes the result independent of input row order.
    grouped = defaultdict(list)
    for sample in sorted(eligible):
        row = eligible[sample]
        grouped[row["population"], row["sex"]].append(row)
    rng = random.Random(seed)
    for group in sorted(grouped):
        rng.shuffle(grouped[group])
    chosen = [eligible[s] for s in always]
    selected = set(always)
    counts = defaultdict(int)
    for row in chosen:
        counts[row["population"], row["sex"]] += 1
    # Round robin gives each population/sex stratum a turn before a second pick.
    for _ in range(per_stratum):
        for key in sorted(grouped):
            if len(chosen) >= total:
                break
            candidates = [r for r in grouped[key] if r["sample_id"] not in selected]
            if candidates and counts[key] < per_stratum:
                row = candidates[0]
                chosen.append(row)
                selected.add(row["sample_id"])
                counts[key] += 1
    if not chosen:
        raise ValueError("No eligible samples")
    support = [dict(population=k[0], sex=k[1], eligible=len(v), selected=counts[k])
               for k, v in sorted(grouped.items())]
    return chosen, support


def validate_fraction(row):
    parsed = dict(row)
    for field in COUNT_FIELDS:
        raw = row[field]
        if isinstance(raw, bool) or isinstance(raw, float) and (not math.isfinite(raw) or not raw.is_integer()):
            raise ValueError(f"Noninteger {field}")
        val = int(raw)
        if val < 0:
            raise ValueError(f"Negative {field}")
        parsed[field] = val
    if (parsed["union_bp"] > parsed["genome_bp"] or
            parsed["union_non_N_bp"] > parsed["genome_non_N_bp"] or
            parsed["genome_non_N_bp"] > parsed["genome_bp"] or
            parsed["union_non_N_bp"] > parsed["union_bp"]):
        raise ValueError("HERV union or non-N counts exceed their denominator")
    for field, numerator, denominator in (
            ("percent_genome", "union_bp", "genome_bp"),
            ("percent_non_N", "union_non_N_bp", "genome_non_N_bp")):
        expected = 100 * parsed[numerator] / parsed[denominator] if parsed[denominator] else None
        raw = row.get(field)
        actual = None if raw in (None, "", "NA") else float(raw)
        if expected is None:
            if actual is not None:
                raise ValueError("A zero denominator requires an unavailable percentage")
        elif actual is None or not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(f"{field} disagrees with its counts")
        parsed[field] = actual
    return parsed


def inspect_census(root, sample):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    done = json.loads((root / "analysis.done.json").read_text())
    result = json.loads((root / "genome_fraction.json").read_text())
    settings = run.get("settings", {})
    if run.get("status") != "complete":
        raise ValueError("Census run is not complete")
    if settings.get("sample") != sample:
        raise ValueError("Census sample differs from metadata")
    if done.get("settings") != settings or result.get("settings") != settings:
        raise ValueError("Summary/completion settings differ from the census run")
    required = {"method", "fasta_sha256", "library_sha256", "repeatmasker_sha256",
                "rmblast_version", "chunk_bp", "halo", "uncertain_flank", "parameters", "code_sha256"}
    if required - settings.keys():
        raise ValueError("Census lacks comparison provenance: " + ", ".join(sorted(required - settings.keys())))
    scope = [str(c).removeprefix("chr") for c in run.get("scope", [])]
    if len(scope) != len(set(scope)) or not scope or any(c not in {str(i) for i in range(1, 23)} for c in scope):
        raise ValueError("This census implementation requires an explicit unique autosomal scope")
    table = root / "genome_fraction.tsv"
    if done.get("results_sha256") != sha256(table):
        raise ValueError("Census table SHA256 differs from the completion receipt")
    with table.open() as f:
        rows = [validate_fraction(r) for r in csv.DictReader(f, delimiter="\t")]
    indexed = {(r["label"], r["category"]): r for r in rows}
    expected_rows = [validate_fraction(r) for r in result.get("categories", [])]
    expected = {(r["label"], r["category"]): r for r in expected_rows}
    if len(indexed) != len(rows) or len(expected) != len(expected_rows) or set(indexed) != set(expected):
        raise ValueError("Duplicate or inconsistent fraction categories")
    for key, row in indexed.items():
        if any(row[f] != expected[key][f] for f in COUNT_FIELDS + PERCENT_FIELDS):
            raise ValueError("Census JSON and TSV disagree")
    if any((label, cat) not in indexed for label in LABELS for cat in CATEGORIES):
        raise ValueError("Missing diploid/matched-reference categories")
    # Code versions remain separate unless proven comparable by a later audit.
    comparison = {k: settings[k] for k in sorted(settings) if k not in {"sample", "library", "repeatmasker"}}
    comparison["chromosomes"] = sorted(scope, key=int)
    comparison_id = hashlib.sha256(json.dumps(comparison, sort_keys=True).encode()).hexdigest()
    output = []
    for cat in CATEGORIES:
        model, ref = (indexed[label, cat] for label in LABELS)
        row = dict(sample_id=sample, category=cat, comparison_id=comparison_id,
                   scope="autosomes:" + ",".join(comparison["chromosomes"]),
                   model_genome_bp=model["genome_bp"], model_non_N_bp=model["genome_non_N_bp"],
                   callable_model_fraction=model["genome_non_N_bp"] / model["genome_bp"] if model["genome_bp"] else None,
                   model_percent=model["percent_genome"], model_percent_non_N=model["percent_non_N"],
                   matched_reference_percent=ref["percent_genome"], matched_reference_percent_non_N=ref["percent_non_N"],
                   delta_percentage_points=(model["percent_non_N"] - ref["percent_non_N"]
                                            if model["percent_non_N"] is not None and ref["percent_non_N"] is not None else None),
                   source=str(root), results_sha256=done["results_sha256"])
        output.append(row)
    return output, comparison


def plan(args):
    roster = metadata(args.metadata, args.group_column, args.sex_column)
    chosen, support = balanced_subset(roster, args.total, args.per_stratum, args.seed,
                                     tuple(s for s in args.include_samples.split(",") if s))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    planned, commands = [], ["#!/usr/bin/env bash", "set -euo pipefail", "# Sequential, resumable individual HERV validation; this file is a plan only."]
    for row in chosen:
        target = Path(args.census_root).resolve() / row["sample_id"]
        action, reason = "run_or_resume", ""
        if (target / "analysis.done.json").exists():
            try:
                inspect_census(target, row["sample_id"])
                action = "reuse_complete"
            except (ValueError, OSError, KeyError, TypeError) as exc:
                action, reason = "blocked_invalid_completed_result", str(exc)
        command = ""
        if action == "run_or_resume":
            command = shlex.join(["bash", str(Path(args.entry).resolve()), "--dat", args.dataset,
                "--build", args.build, "--sample", row["sample_id"], "--out", str(target),
                "--replace", "FALSE", "--no-reuse", "--jobs", str(args.jobs), "--cores", str(args.cores)])
            commands.append(command)
        planned.append(dict(row, run_dir=str(target), action=action, reason=reason, command=command))
    write_tsv(out / "sample_plan.tsv", planned, ["sample_id", "population", "sex", "unrelated", "run_dir", "action", "reason", "command"])
    write_tsv(out / "stratum_support.tsv", support, ["population", "sex", "eligible", "selected"])
    (out / "run_selected.sh").write_text("\n".join(commands) + "\n")
    save_json(out / "plan.json", dict(status="planned", version=VERSION, metadata_sha256=sha256(args.metadata),
        selected=len(planned), seed=args.seed, requested_maximum=args.total, per_stratum=args.per_stratum,
        note="No analyses launched. Zero-selected strata are disclosed; this small validation subset is not a population prevalence sample.",
        all_sample_route="1.herv.py cohort burden: shared reference atlas and cohort PGEN; no per-person RepeatMasker required."))
    return {"status": "planned", "selected": len(planned), "out": str(out)}


def summarize(args):
    roster = metadata(args.metadata, args.group_column, args.sex_column)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows, audit, comparisons = [], [], {}
    for meta in roster:
        root = Path(args.census_root) / meta["sample_id"]
        if meta["include"] != "1":
            audit.append(dict(sample_id=meta["sample_id"], status="excluded_metadata", reason="include=0"))
            continue
        try:
            found, comparison = inspect_census(root, meta["sample_id"])
            for row in found:
                row.update(population=meta["population"], sex=meta["sex"], unrelated=meta["unrelated"])
                rows.append(row)
                comparisons[row["comparison_id"]] = comparison
            audit.append(dict(sample_id=meta["sample_id"], status="verified_summary", reason=""))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            audit.append(dict(sample_id=meta["sample_id"], status="unavailable", reason=str(exc)))
    fields = ["sample_id", "population", "sex", "unrelated", "category", "comparison_id", "scope",
              "model_genome_bp", "model_non_N_bp", "callable_model_fraction", "model_percent", "model_percent_non_N",
              "matched_reference_percent", "matched_reference_percent_non_N", "delta_percentage_points", "source", "results_sha256"]
    write_tsv(out / "individual_census.tsv", rows, fields)
    write_tsv(out / "sample_audit.tsv", audit, ["sample_id", "status", "reason"])
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["comparison_id"], row["population"], row["sex"], row["category"]].append(row)
    summary = []
    for key, rr in sorted(grouped.items()):
        for metric in ("model_percent_non_N", "matched_reference_percent_non_N", "delta_percentage_points", "callable_model_fraction"):
            vals = [r[metric] for r in rr if r[metric] is not None]
            summary.append(dict(comparison_id=key[0], population=key[1], sex=key[2], category=key[3], metric=metric,
                n=len(vals), mean=statistics.mean(vals) if vals else None,
                median=statistics.median(vals) if vals else None,
                sd=statistics.stdev(vals) if len(vals) > 1 else None,
                minimum=min(vals) if vals else None, maximum=max(vals) if vals else None,
                inference="descriptive_validation_subset"))
    write_tsv(out / "group_summary.tsv", summary, ["comparison_id", "population", "sex", "category", "metric", "n", "mean", "median", "sd", "minimum", "maximum", "inference"])
    state = "complete" if rows else "unavailable"
    save_json(out / "manifest.json", dict(status=state, version=VERSION, metadata_sha256=sha256(args.metadata),
        verified_samples=sum(r["status"] == "verified_summary" for r in audit),
        unavailable_samples=sum(r["status"] == "unavailable" for r in audit), comparisons=comparisons,
        verification="Checks saved summary SHA256, count/percentage consistency, sample, completion and scientific settings; does not rehash raw genome/chunk files.",
        scope_note="The current individual census is autosomal. X/Y are retained in the separate allele3 chromosome analyses; no autosomal estimate is labeled whole-genome.",
        interpretation="Reference-conditioned sequence-model fractions and matched callable-reference contrasts. Not independent assemblies, insertion counts, or proof that HERV insertion created a third nucleotide allele."))
    if not rows:
        raise ValueError("No verified complete individual census; see sample_audit.tsv")
    return {"status": state, "verified_samples": len(rows) // len(CATEGORIES), "out": str(out)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("plan", "summarize"):
        q = sub.add_parser(name)
        q.add_argument("--metadata", required=True, help="Existing 1KG panel or TSV with sample/sample_id, population/super_pop, sex/gender")
        q.add_argument("--group-column", default="auto", help="Exact column; auto prefers super_pop then population/pop")
        q.add_argument("--sex-column", default="auto", help="Exact column; auto accepts sex or gender")
        q.add_argument("--census-root", required=True, help="Contains one census directory per sample ID")
        q.add_argument("--out", required=True)
        if name == "plan":
            q.add_argument("--entry", default=str(Path(__file__).resolve().parents[1] / "2.herv.sh"))
            q.add_argument("--dataset", default="1kg")
            q.add_argument("--build", choices=("37", "38"), default="38")
            q.add_argument("--total", type=int, default=20)
            q.add_argument("--per-stratum", type=int, default=2)
            q.add_argument("--include-samples", default="", help="Comma-separated eligible IDs counted within --total")
            q.add_argument("--seed", type=int, default=20261102)
            q.add_argument("--jobs", type=int, default=2)
            q.add_argument("--cores", type=int, default=8)
    args = p.parse_args(argv)
    if args.command == "plan" and (not 1 <= args.jobs <= 8 or args.cores < 1):
        p.error("jobs must be 1..8 and cores must be positive")
    try:
        print(json.dumps(plan(args) if args.command == "plan" else summarize(args), indent=2))
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"herv-cohort: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

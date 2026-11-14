"""Summarize committed split results without reading FASTQs or loading indexes."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import html
import json
from pathlib import Path
from storage_layout import metadata_dir
from qc_paths import configured_split_root

from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, load_config, utcnow


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def table(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def load_sample(folder, sample_id):
    """Check completion, artifact sizes and metadata count conservation cheaply."""
    folder = Path(folder)
    metadata = metadata_dir(folder)
    receipt = read_json(metadata / "complete.json")
    if receipt.get("status") != "complete":
        raise ValueError("Completion receipt is not complete")
    outputs = receipt.get("outputs", {})
    prefix = "qc/" if metadata != folder else ""
    required = {prefix + name for name in ("summary.json", "counts.tsv", "bin-counts.tsv", "provenance.json")}
    if not required <= outputs.keys():
        raise ValueError("Completion receipt lacks required metadata")
    for name, saved in outputs.items():
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Unsafe receipt path")
        path = folder / rel
        if path.is_symlink() or not path.is_file() or path.stat().st_size != saved["state"]["size"]:
            raise ValueError(f"Missing/changed committed artifact: {name}")
    s = read_json(metadata / "summary.json")
    if s["sample_id"] != sample_id:
        raise ValueError(f"Sample ID mismatch: {s['sample_id']} != {sample_id}")
    if not s.get("count_conservation") or s.get("requested_limit") != 0 or s.get("validated_scope") != "entire_input":
        raise ValueError("Result is not an entire-input, count-conserving split")
    n = int(s["fragments"])
    if n <= 0 or int(s["reads"]) != n * (2 if s["paired"] else 1):
        raise ValueError("Invalid fragment/read total")
    with (metadata / "counts.tsv").open(newline="", encoding="utf-8") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    counts, statuses, groups, bins = Counter(), Counter(), Counter(), Counter()
    for row in rows:
        amount = int(row["fragments"])
        if amount < 0 or int(row["reads"]) != amount * (2 if s["paired"] else 1):
            raise ValueError("Invalid counts row")
        counts["fragments"] += amount; counts["reads"] += int(row["reads"])
        statuses[row["status"]] += amount; groups[row["group"]] += amount; bins[row["bin"]] += amount
    if counts["fragments"] != n or counts["reads"] != int(s["reads"]):
        raise ValueError("counts.tsv totals differ from summary")
    if dict(statuses) != s["statuses"] or dict(groups) != s["groups"]:
        raise ValueError("counts.tsv status/group totals differ from summary")
    with (metadata / "bin-counts.tsv").open(newline="", encoding="utf-8") as h:
        bin_rows = list(csv.DictReader(h, delimiter="\t"))
    saved_bins = {r["bin"]: int(r["fragments"]) for r in bin_rows}
    if len(saved_bins) != len(bin_rows) or saved_bins != dict(bins):
        raise ValueError("Physical bin counts differ from taxon counts")
    kraken = s.get("kraken_statuses", s["statuses"])
    if sum(kraken.values()) != n:
        raise ValueError("Kraken status totals do not conserve fragments")
    provenance = read_json(metadata / "provenance.json")
    return s, rows, provenance["context"]


SAMPLE_FIELDS = ["sample_id", "fragments", "reads", "unclassified_fragments", "unclassified_fraction",
    "higher_rank_fragments", "higher_rank_fraction", "species_fragments", "species_fraction",
    "classified_fragments", "classified_fraction", "kraken_human_fragments", "final_human_fragments",
    "human_screen_enabled", "host_positive_fragments", "host_rescued_fragments", "species_taxa_total",
    "species_taxa_materialized", "species_other_fragments", "physical_bins", "fastq_files",
    "classification_seconds", "split_seconds", "archive_seconds", "artifact_bytes",
    "confidence", "minimum_hit_groups", "database_identity", "input_digest", "output_dir"]

# Append fields; keep the existing table prefix and stored status codes unchanged.
SAMPLE_FIELDS += ["unresolved_high_fragments", "unresolved_high_fraction",
                  "known_group_high_fragments", "known_group_high_fraction",
                  "kraken_human_fraction", "final_human_fraction"]


def sample_row(s, context, folder, counts=None):
    row = {k: s.get(k, "") for k in SAMPLE_FIELDS}
    row.update(output_dir=str(folder), confidence=context.get("confidence"),
               minimum_hit_groups=context.get("minimum_hit_groups"),
               database_identity=context.get("database_identity"),
               human_screen_enabled=s.get("split_policy", {}).get("human_host", {}).get("enabled", False))
    status = s.get("kraken_statuses", s["statuses"])
    for name, key in (("unclassified", "unclassified"), ("higher_rank", "higher_rank"),
                      ("species", "assigned_species")):
        row[name + "_fragments"] = status.get(key, 0)
        row[name + "_fraction"] = status.get(key, 0) / s["fragments"]
    row["classified_fragments"] = s["fragments"] - row["unclassified_fragments"]
    row["classified_fraction"] = row["classified_fragments"] / s["fragments"]
    # The unresolved group includes taxid=0 as well as cross-group ancestors;
    # only the latter belongs to root taxon_unresolved.
    group_status = s.get("kraken_group_status")
    if group_status is not None or (counts is not None and not s.get("host_rescued_fragments", 0)):
        if group_status is not None:
            unresolved = int(group_status.get("unresolved", {}).get("higher_rank", 0))
        else:
            unresolved = sum(int(r["fragments"]) for r in counts
                             if r["group"] == "unresolved" and r["status"] == "higher_rank")
        row["unresolved_high_fragments"] = unresolved
        row["known_group_high_fragments"] = row["higher_rank_fragments"] - unresolved
        for key in ("unresolved_high", "known_group_high"):
            row[key + "_fraction"] = row[key + "_fragments"] / s["fragments"]
    for key in ("kraken_human", "final_human"):
        value = s.get(key + "_fragments")
        row[key + "_fraction"] = value / s["fragments"] if isinstance(value, int) else ""
    return row


def write_readable_report(out, report, rows, states):
    """A single overview; the TSVs retain exact values for downstream analysis."""
    columns = [("sample_id", "Sample"), ("fragments", "Pairs / fragments"),
               ("unclassified_fraction", "Unclassified"),
               ("higher_rank_fraction", "Above species"),
               ("species_fraction", "Species assigned"),
               ("unresolved_high_fraction", "Cross-group ancestor"),
               ("kraken_human_fragments", "Kraken human candidates"),
               ("kraken_human_fraction", "Human candidate fraction")]
    def cell(key, value):
        if value is None or value == "":
            return "—"
        if key.endswith("_fraction"):
            return f"{float(value):.6%}"
        if isinstance(value, int):
            return f"{value:,}"
        return html.escape(str(value))
    header = "".join(f"<th>{html.escape(title)}</th>" for _, title in columns)
    body = "".join("<tr>" + "".join(f"<td>{cell(key, row.get(key))}</td>" for key, _ in columns) + "</tr>"
                   for row in rows)
    problems = [state for state in states if state["state"] != "complete"]
    issues = "".join(f"<li><b>{html.escape(s['sample_id'])}</b>: {html.escape(s['state'])} — "
                     f"{html.escape(s['detail'])}</li>" for s in problems)
    text = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>Metagenomic split report</title><style>
body{{font:16px/1.6 system-ui,sans-serif;max-width:1450px;margin:32px auto;padding:0 24px;color:#172b3a}}
h1,h2{{line-height:1.25}}table{{border-collapse:collapse;font-size:14px;width:100%}}
th,td{{padding:10px;border-bottom:1px solid #dae1e7;text-align:right}}th:first-child,td:first-child{{text-align:left}}
th{{background:#eef4f7}}.note{{background:#f4f7f9;padding:16px;border-left:4px solid #3782a3}}
.scroll{{overflow:auto}}code{{background:#eef2f4;padding:2px 5px}}</style>
<h1>Metagenomic split report</h1>
<p>{report['complete_samples']} / {report['expected_samples']} requested samples complete.
State: <b>{html.escape(report['state'])}</b>. Generated {html.escape(report['generated_utc'])}.</p>
<div class="note">Every fraction uses all input fragments (one fragment is one pair for paired input).
<b>Species assigned, above species, and unclassified are mutually exclusive.</b>
<code>taxon_unresolved</code> contains higher-rank calls; a root file contains calls without a unique group.
<code>TAXID_name</code> files contain separate species; <code>species_other</code> contains species below the
file-materialization threshold in that group. Species counts include both. These are storage/taxonomic categories, not confidence levels.
Every exact taxid remains in counts.tsv, and bin-counts.tsv lists the taxids sharing each physical file stem.
Human candidates require independent alignment validation. A disabled human screen is not a negative test.</div>
<h2>Sample counts</h2><div class="scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>
<h2>Completeness</h2><ul>{issues or '<li>All requested samples passed completion metadata, artifact sizes and count conservation.</li>'}</ul>
<p>Verification does not reread or rehash large FASTQs. Read assignment fractions are not cell or organism abundances.</p>
<p>Exact tables: <a href="samples.tsv">samples.tsv</a>, <a href="statuses.tsv">statuses.tsv</a>,
<a href="groups.tsv">groups.tsv</a>, <a href="species-by-sample.tsv">species-by-sample.tsv</a>,
<a href="sample-status.tsv">sample-status.tsv</a>.</p></html>
"""
    (Path(out) / "report.html").write_text(text, encoding="utf-8")


def summarize(root, sample_ids, out, compare_root=None):
    root, out = Path(root), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sample_rows, states, status_rows, group_rows, comparisons = [], [], [], [], []
    species_totals = Counter(); species_meta = {}; cohort = Counter()
    fields = ["sample_id", "species_taxid", "name", "group", "status", "fragments", "fraction_of_sample"]
    with (out / "species-by-sample.tsv").open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sid in sample_ids:
            folder = root / sid
            if not (metadata_dir(folder) / "complete.json").is_file():
                states.append({"sample_id": sid, "state": "missing_or_incomplete", "detail": str(folder)})
                continue
            try:
                s, rows, ctx = load_sample(folder, sid)
                row = sample_row(s, ctx, folder, rows)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                states.append({"sample_id": sid, "state": "invalid", "detail": str(exc)})
                continue
            states.append({"sample_id": sid, "state": "complete", "detail": "completion metadata + artifact sizes + count conservation"})
            sample_rows.append(row)
            for key in ("fragments", "reads", "unclassified_fragments", "higher_rank_fragments", "species_fragments"):
                cohort[key] += row[key]
            for layer, values in (("final", s["statuses"]), ("kraken", s.get("kraken_statuses", s["statuses"]))):
                for key, n in sorted(values.items()):
                    status_rows.append(dict(sample_id=sid, layer=layer, status=key, fragments=n, fraction=n/s["fragments"]))
            for layer, values in (("final", s["groups"]), ("kraken", s.get("kraken_groups", s["groups"]))):
                for key, n in sorted(values.items()):
                    group_rows.append(dict(sample_id=sid, layer=layer, group=key, fragments=n, fraction=n/s["fragments"]))
            local = Counter()
            for r in rows:
                tid = int(r["species_taxid"] or 0)
                if tid:
                    key = (tid, r["status"])
                    local[key] += int(r["fragments"]); species_meta[key] = (r["name"], r["group"])
            for key, n in sorted(local.items()):
                tid, status = key; name, group = species_meta[key]
                writer.writerow(dict(sample_id=sid, species_taxid=tid, name=name, group=group,
                    status=status, fragments=n, fraction_of_sample=n/s["fragments"]))
                species_totals[key] += n
            if compare_root:
                comparison = {"sample_id": sid, "state": "missing_or_invalid"}
                try:
                    other, other_rows, other_ctx = load_sample(Path(compare_root) / sid, sid)
                    if (s["input_digest"] != other["input_digest"] or s["fragments"] != other["fragments"] or
                            ctx["database_identity"] != other_ctx["database_identity"]):
                        raise ValueError("Comparison requires identical input digest, fragment count and database")
                    b = sample_row(other, other_ctx, Path(compare_root) / sid, other_rows)
                    comparison.update(state="comparable", baseline_confidence=row["confidence"],
                        comparison_confidence=b["confidence"], baseline_hit_groups=row["minimum_hit_groups"],
                        comparison_hit_groups=b["minimum_hit_groups"])
                    for k in ("unclassified", "higher_rank", "species", "classified"):
                        comparison[k + "_delta_fragments"] = b[k + "_fragments"] - row[k + "_fragments"]
                        comparison[k + "_delta_fraction"] = b[k + "_fraction"] - row[k + "_fraction"]
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    comparison["detail"] = str(exc)
                comparisons.append(comparison)
    table(out / "samples.tsv", SAMPLE_FIELDS, sample_rows)
    table(out / "sample-status.tsv", ["sample_id", "state", "detail"], states)
    table(out / "statuses.tsv", ["sample_id", "layer", "status", "fragments", "fraction"], status_rows)
    table(out / "groups.tsv", ["sample_id", "layer", "group", "fragments", "fraction"], group_rows)
    table(out / "species-total.tsv", ["species_taxid", "name", "group", "status", "fragments"],
        [dict(species_taxid=k[0], status=k[1], name=species_meta[k][0], group=species_meta[k][1], fragments=n)
         for k, n in species_totals.most_common()])
    if compare_root:
        table(out / "sensitivity.tsv", ["sample_id", "state", "detail", "baseline_confidence", "comparison_confidence",
              "baseline_hit_groups", "comparison_hit_groups"] +
              [k + suffix for k in ("unclassified", "higher_rank", "species", "classified")
               for suffix in ("_delta_fragments", "_delta_fraction")], comparisons)
    complete = len(sample_rows)
    report = dict(schema="metag-cohort-summary-v1", generated_utc=utcnow(), input_root=str(root),
        expected_samples=len(sample_ids), sample_ids=list(sample_ids), complete_samples=complete,
        incomplete_or_invalid_samples=len(sample_ids)-complete, totals=dict(cohort),
        state="complete" if complete == len(sample_ids) and complete else "partial",
        comparison_failures=sum(x["state"] != "comparable" for x in comparisons),
        validation="Metadata hashes, artifact existence/sizes and count conservation; large FASTQs/tars not rehashed.",
        interpretation="Fractions use all fragments. Species totals are read-assignment counts, not organism abundances. Disabled human screens are not negative tests. Sensitivity differences are not accuracy estimates.")
    atomic_json(out / "split-summary.json", report)
    write_readable_report(out, report, sample_rows, states)
    return report


@managed_outputs
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replace", type=replace_bool, default=True)
    p.add_argument("--config"); p.add_argument("--root")
    p.add_argument("--manifest", help="Frozen split jobs.json defines expected sample set and output root")
    p.add_argument("--input-root", help="Directory directly containing DATASET/SAMPLE results")
    p.add_argument("--compare-root", help="Sensitivity result root; input digest and database must match")
    p.add_argument("--sample", action="append", help="Filter exact DATASET/SAMPLE; may repeat")
    p.add_argument("--out", help="Default: work_dir/qc/RUN_LABEL/[DATASET/SAMPLE/]summary")
    p.add_argument("--require-complete", action="store_true", help="Write reports, then exit nonzero if samples or comparisons are missing/invalid")
    args = p.parse_args(argv)
    cfg = load_config(args.config, args.root)
    root = args.input_root
    if args.manifest:
        plan = read_json(args.manifest)
        context = plan["context"].get("split", plan["context"])
        samples = [s["sample_id"] for s in context["samples"]]
        frozen = read_json(plan["config"])
        root = root or frozen["output_dir"]
        if not args.config and not args.root:
            cfg = frozen
    else:
        root = root or str(configured_split_root(cfg))
        # Include interrupted/uncommitted output directories. Looking only for
        # complete.json would silently drop precisely the failures we must show.
        samples = sorted(str(path.relative_to(root)) for path in Path(root).glob("*/*") if path.is_dir())
    if not Path(root).is_dir() and not args.sample and not args.manifest:
        p.error(f"Split root does not exist: {root}")
    if args.sample:
        if args.manifest and not set(args.sample) <= set(samples):
            p.error("Requested sample is absent from manifest")
        samples = args.sample
    if len(samples) != len(set(samples)):
        p.error("Duplicate sample IDs")
    if any(Path(s).is_absolute() or ".." in Path(s).parts or len(Path(s).parts) != 2 for s in samples):
        p.error("Sample IDs must have form DATASET/SAMPLE")
    report_root = Path(cfg["work_dir"]) / "qc" / cfg.get("run_label", "split")
    if args.sample:
        if len(samples) == 1:
            report_root /= samples[0]
        else:
            report_root /= "selected"
    out = args.out or str(report_root / "summary")
    if not prepare_result(out, replace=args.replace, cfg=cfg, protected=[root, args.compare_root, args.manifest]):
        saved = read_json(Path(out) / "split-summary.json")
        print(json.dumps(saved, indent=2))
        if args.require_complete and (saved["state"] != "complete" or saved["comparison_failures"]):
            raise SystemExit(2)
        return
    report = summarize(Path(root).resolve(), sorted(samples), out, args.compare_root)
    if report["state"] == "complete" and not report["comparison_failures"]:
        finish_result(out)
    print(json.dumps(report, indent=2))
    if args.require_complete and (report["state"] != "complete" or report["comparison_failures"]):
        raise SystemExit(2)


if __name__ == "__main__":
    main()

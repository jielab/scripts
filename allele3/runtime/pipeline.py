#!/usr/bin/env python3
"""Resume modern results, then run the existing ancient/PCA/optional Atlas workflow.

Keep orchestration outside f/: adding an entry point must not invalidate the
modern scanner's code signature or rewrite a completed run.json used by joins.
"""
from __future__ import annotations

import fcntl
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
import cli
import evolve_cli
from core import VERSION, fingerprint, stable_hash
from gnomad_io import discover_inputs


SUMMARY_FILES = (
    "REPORT.md", "qc_summary.json", "frequency_counts.tsv", "population_counts.tsv",
    "population_summary.tsv", "population_classes.tsv", "relationship_counts.tsv",
    "annotation_counts.tsv", "site_ac_spectrum.tsv", "singleton_sensitivity.tsv",
    "site_ac_spectrum_availability.json",
)


def validate_completed_inputs(a, options, run):
    """Reuse recorded scientific results, checking inputs/settings and shard identity.

    Like evolve's load_run(), completed scans retain their recorded code version;
    a newer entry point alone is not a reason to repeat a whole-genome scan.
    Incomplete scans still use run_scan's strict original code signature check.
    """
    identity = run["identity"]
    if stable_hash(identity) != run["signature"]:
        raise ValueError("Modern manifest identity/signature mismatch; restore the original run.json")
    if options != run["options"] or options != identity["options"]:
        raise ValueError("Modern scan settings differ from completed results; use the original settings or a new --dir-modern")
    paths = discover_inputs(a.vcf, None if a.vcf else a.vcf_dir,
                            None if a.vcf else (a.chromosomes or "autosomes"))
    current = dict(inputs=[fingerprint(p) for p in paths], fasta=fingerprint(options["fasta"]),
                   fasta_index=fingerprint(options["fasta"] + ".fai"),
                   callable_bed=fingerprint(options["callable_bed"]) if options["callable_bed"] else None)
    if Path(options["fasta"] + ".gzi").exists():
        current["fasta_gzi"] = fingerprint(options["fasta"] + ".gzi")
    if (any(current.get(key) != identity.get(key) for key in (*current, "fasta_gzi"))
            or len(run["shards"]) != len(paths)
            or len(set(run["shards"])) != len(paths)):
        raise ValueError("Modern inputs differ from completed results; restore the inputs or use a new --dir-modern")


def summary_complete(out, run, a):
    if run["status"] != "complete" or run.get("summary_parameters") != dict(
            min_training_sites=a.min_training_sites, make_plots=not a.no_plots):
        return False
    directory = out / "summary"
    if run['options'].get('with_architecture') and not (directory/'architecture_manifest.json').is_file():
        return False
    manifest = directory / "summarize_manifest.json"
    if not manifest.is_file():
        return False
    summary = json.loads(manifest.read_text())
    expected_inputs = [str(Path(shard) / "strata.tsv.gz") for shard in run["shards"]]
    return (summary.get("status") == "complete" and summary.get("inputs") == expected_inputs
            and bool(summary.get("outputs"))
            and all(Path(path).is_file() for path in summary["outputs"])
            and all((directory / name).is_file() for name in SUMMARY_FILES))


def ensure_modern(a):
    options = cli.validated_options(a)
    out = Path(a.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".allele3.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(f"Another allele3 process is using {out}") from exc
        manifest = out / "run.json"
        previous = json.loads(manifest.read_text()) if manifest.exists() else None
        if previous and previous.get("status") in ("scanned", "complete"):
            validate_completed_inputs(a, options, previous)
            run = cli.load_run(out)  # Check completion markers and all shard output fingerprints.
            print(f"[resume] modern scan: {len(run['shards'])} completed shards verified; skipped", flush=True)
        else:
            print("[pipeline] modern scan: start/resume unfinished chromosomes", flush=True)
            run = cli.run_scan(a, options)
        if summary_complete(out, run, a):
            print(f"[resume] modern summary: {out / 'summary'}; skipped", flush=True)
        else:
            print("[pipeline] modern summary", flush=True)
            cli.summarize_run(out, a.min_training_sites, not a.no_plots)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["--version"]:
        print(VERSION)
        return 0
    if "--help" in argv or not argv:
        print("Default pipeline: resume modern scan/summary -> ancient analysis -> PCA -> report.\n"
              "Use allele3.sh --with-atlas to add AlphaGenome; --no-pca to skip PCA.\n"
              "--dir-modern selects the modern output; --dir-out selects the evolution output.\n")
        cli.parser().print_help()
        evolve_cli.parser().print_help()
        return 0
    try:
        if argv[0] != "pipeline" or "--" not in argv:
            raise ValueError("Use allele3.sh to supply modern and evolution arguments")
        boundary = argv.index("--")
        modern = cli.parser().parse_args(argv[1:boundary])
        evolution_args = argv[boundary + 1:]
        evolution = evolve_cli.parser().parse_args(evolution_args)
        if modern.command != "run" or evolution.command != "evolve":
            raise ValueError("Pipeline requires modern run followed by evolve")
        if Path(modern.out).resolve() != Path(evolution.modern_out).resolve():
            raise ValueError("Pipeline's modern output and evolution input must match")
        if Path(modern.out).resolve() == Path(evolution.out).resolve():
            raise ValueError("Modern and evolution outputs must be different directories")
        ensure_modern(modern)
        print(f"[pipeline] ancient analysis; PCA={'on' if evolution.with_pca else 'off'}; "
              f"AlphaGenome={'on' if evolution.with_atlas else 'off (optional; no API key read)'}", flush=True)
        return evolve_cli.main(evolution_args)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f"allele3 pipeline: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

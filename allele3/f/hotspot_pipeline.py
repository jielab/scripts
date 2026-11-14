#!/usr/bin/env python3
"""Genome-first FINAL workflow: discover, freeze controls, annotate, report.

Existing scan stores are read-only. No HERV label or Atlas score is used to
select hotspots. Whole-window inference stays at the cache's real resolution.
Optional local Atlas/ancient stages report availability explicitly. Completed
stages survive a later failure and are reused only with matching receipts.
"""
from __future__ import annotations

# Hotspot Pipeline
# ======================================================================
import argparse
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import html
import importlib
import json
import math
import os
from pathlib import Path
import sys

from analysis_integrity import InputPath, atomic_json, check_output_paths, identity, inventory, validate_manifest

ROOT = Path(__file__).resolve().parents[1]
PATH_ROOT = next(p for p in ROOT.parents if (p / "0f/workflow_paths.py").is_file())
sys.path.insert(0, str(PATH_ROOT / "0f"))
from workflow_paths import configured_path

VERSION = "allele3-FINAL-7.0.0"
LIMITATIONS = [
    "Discovery counts sites with at least two count-supported nonreference alleles among sites with at least one. Joint REF presence and a common four-allele AN are not established by these scan counters.",
    "Working-model NB2/BY nominations require independent technical review and replication. AC5/AC10 are sensitivity tiers in the same cohort, not independent replication.",
    "Source resolution is preserved. Coarse counts cannot establish finer hotspots or be proportionally corrected for excluded bases.",
    "HERV, gene and other BED overlaps are post-discovery descriptions. Whole-window exclusion masks veto nominations without modifying source counts.",
    "A bounded sample of multi-ALT sites describes that target panel, not all possible mutations in the entire window. Window-wide tolerance needs an all-eligible or independent prespecified sampling frame and stronger callable population evidence.",
    "AVI is a frequency-trained impact ranking; its 18 exported feature contributions are SHAP attributions, not measured molecular functions. Neither constitutes a clinical pathogenicity classification or proof of pleiotropy/balancing selection.",
    "Published ancient selection contrasts are reused only after exact reciprocal build/allele matching. A published biallelic effect is not assigned to a third nucleotide; imputed-panel absence is not extinction.",
]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def code_fingerprint():
    paths = [*ROOT.glob("*.sh"), ROOT / "allele3.env", *sorted((ROOT / "f").glob("*")),
             *sorted((ROOT / "config").glob("*.tsv"))]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
            if p.is_file() and p.suffix in {".py", ".sh", ".env", ".html", ".js", ".tsv"}}


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def count_rows(path, maximum=None):
    n = 0
    for _ in read_table(path):
        n += 1
        if maximum is not None and n > maximum:
            raise ValueError(f"Target table exceeds the declared --max-targets={maximum}: {path}")
    return n


def source_paths(args, metadata=True):
    """Bound large source identities without reading genomes or the Atlas."""
    paths, protected = [], []
    if args.summary:
        paths.append(Path(args.summary))
    else:
        store = Path(args.scan_dir)
        protected.append(InputPath(store, "directory"))
        paths.extend(store / n for n in ("run.json", "relocation.json", "adoption_receipt.json", "input_locations.json")
                     if n == "run.json" or (store / n).is_file())
        if metadata and (store / "run.json").is_file():
            from scan_store import load_run_metadata
            run = load_run_metadata(store)
            for entry in run.get("shards", []):
                shard = Path(entry)
                if not shard.is_absolute():
                    shard = store / shard
                protected.append(InputPath(shard, "directory"))
                paths.extend(shard / n for n in ("COMPLETE.json", "scan.json", "strata.tsv.gz", "scan_scope.tsv",
                                                  "multiallelic_sites.tsv.gz") if (shard / n).is_file())
    paths.extend(Path(p) for p in args.exclude_bed + args.site_table)
    paths.extend(Path(s.split("=", 1)[1]) for s in args.annotation_bed if "=" in s)
    for key in ("covariates", "atlas_resources", "atlas_variants", "atlas_metrics", "atlas_targets",
                "selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38"):
        if getattr(args, key, None):
            paths.append(Path(getattr(args, key)))
    for key in ("fasta38", "fasta37"):
        if getattr(args, key, None):
            paths.extend(Path(getattr(args, key) + s) for s in (".fai", ".gzi") if Path(getattr(args, key) + s).is_file())
    # Expand resource paths using the documented schema, retaining missing paths
    # as explicit identities. Do not inspect BGZF contents during this step.
    from atlas_scores import default_resources
    from atlas_scores import manifest_source_paths
    if not args.atlas_resources:
        args.atlas_resources = str(default_resources(args.dir_data, all_resources=True))
    config = Path(args.atlas_resources)
    paths.append(config)
    if config.is_file() and args.with_atlas != "FALSE":
        paths.extend(manifest_source_paths(config, args.dir_data))
    paths = sorted({p.expanduser().resolve() for p in paths})
    protected.extend(InputPath(p) for p in paths)
    return paths, protected


def auto_site_tables(args):
    if args.site_table or args.summary:
        return list(args.site_table)
    from scan_store import load_run_metadata
    run = load_run_metadata(args.scan_dir)
    result = []
    for entry in run.get("shards", []):
        shard = Path(entry)
        if not shard.is_absolute():
            shard = Path(args.scan_dir) / shard
        p = shard / "multiallelic_sites.tsv.gz"
        if p.is_file():
            result.append(str(p.resolve()))
    if result and args.site_sampling != "selected_multialt_candidates":
        raise ValueError("Auto-detected multiallelic caches are selected_multialt_candidates; provide an explicit independently justified site table for another sampling frame")
    # Completed stores can preserve lexical chromosome/shard order. Read only
    # the first data row of each cache and sort disjoint source shards by their
    # actual coordinates; the exporter still rejects overlaps/duplicates.
    def first_coordinate(path):
        row = next(read_table(path), None)
        if row is None:
            return (999, 0, path)
        c = row["chrom"].lower().removeprefix("chr").upper()
        order = int(c) if c.isdigit() else {"X":23,"Y":24,"M":25,"MT":25}.get(c,998)
        return (order, int(row["pos"]), path)
    return sorted(result, key=first_coordinate)


class StageRunner:
    def __init__(self, out, signature):
        self.out, self.signature = Path(out), signature
        self.results = {}

    def run(self, name, module, argv):
        target = self.out / "runtime" / name
        expected = digest(dict(pipeline_request=self.signature, module=module, argv=argv))
        manifest = target / "manifest.json"
        if manifest.is_file():
            saved = json.loads(manifest.read_text())
            if saved.get("pipeline_stage_signature") == expected and saved.get("status") in {"complete", "complete_with_exclusions"}:
                saved = validate_manifest(target)
                self.results[name] = dict(status=saved["status"], reused=True, manifest=str(manifest))
                print(f"[hotspots] reuse {name}", flush=True)
                return saved
        print(f"[hotspots] start {name}", flush=True)
        if module == "hotspot_characterize":
            parse, run_stage = characterize_parser, characterize_run
        elif module == "window_selection":
            parse, run_stage = selection_parser, selection_run
        else:
            worker = importlib.import_module(module)
            parse, run_stage = worker.parser, worker.run
        result = run_stage(parse().parse_args([*argv, "--out", str(target)]))
        saved = validate_manifest(target)
        saved["pipeline_stage_signature"] = expected
        atomic_json(manifest, saved)
        self.results[name] = dict(status=saved["status"], reused=False, manifest=str(manifest))
        print(f"[hotspots] complete {name}", flush=True)
        return result

    def unavailable(self, name, reason, status="not_available"):
        self.results[name] = dict(status=status, reason=reason)


def characterize_args(args, windows, sites):
    result = ["--windows", str(windows), "--build", "GRCh38", "--threshold", str(args.primary_threshold),
              "--controls-per-case", str(args.controls_per_case), "--match-columns", args.match_columns,
              "--max-match-distance", str(args.max_match_distance), "--control-gap-bp", str(args.control_gap_bp),
              "--sites-per-window", str(args.sites_per_window), "--max-targets", str(args.max_targets),
              "--seed", str(args.seed)]
    if args.characterize_window_bp:
        result += ["--window-bp", str(args.characterize_window_bp)]
    if args.covariates:
        result += ["--covariates", args.covariates]
    for item in args.annotation_bed:
        result += ["--annotation-bed", item]
    for item in sites:
        result += ["--site-table", item]
    if sites:
        result += ["--site-sampling", args.site_sampling]
    return result


def finite(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def report(out, args, stages, window_manifest, atlas_panel=None):
    """Self-contained review of actual output; no remote assets or placeholder findings."""
    all_rows = list(read_table(out / "runtime/windows/all_windows.tsv.gz"))
    primary = [r for r in all_rows if r.get("is_primary") == "1"]
    nominations = [r for r in all_rows if r.get("nominated") == "1"]
    plotted = [r for r in primary if finite(r.get("oe")) is not None]
    lengths = {str(c): max([int(r["end0"]) for r in primary if r["chrom"].removeprefix("chr") == str(c)], default=0) for c in range(1,23)}
    offsets, total = {}, 0
    for c in lengths:
        offsets[c], total = total, total + lengths[c]
    top = max([float(r["oe"]) for r in plotted], default=2.)
    ymax = max(3., math.ceil(top))
    elements = []
    for y in (0., 1., 2., ymax):
        py = 310 - 270 * y / ymax
        elements.append(f'<line x1="56" x2="1085" y1="{py}" y2="{py}" stroke="{"#ec9e4b" if y == 2 else "#d9e1eb"}"/><text x="45" y="{py+4}" text-anchor="end">{y:g}</text>')
    for c in lengths:
        if not lengths[c] or not total:
            continue
        x = 56 + 1029 * (offsets[c] + lengths[c] / 2) / total
        elements.append(f'<text x="{x}" y="332" text-anchor="middle">{c}</text>')
    for r in plotted:
        c = r["chrom"].removeprefix("chr")
        x = 56 + 1029 * (offsets[c] + (int(r["start0"]) + int(r["end0"])) / 2) / max(1,total)
        y = 310 - 270 * float(r["oe"]) / ymax
        label = f'{r["window_id"]} | O/E={float(r["oe"]):.3f} | q={r.get("q_value")} | {r.get("nomination_reason")}'
        color = "#c93e50" if r.get("nominated") == "1" else ("#247b91" if int(c)%2 else "#6d98a4")
        elements.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{4 if r.get("nominated")=="1" else 2}" fill="{color}" opacity=".8"><title>{html.escape(label)}</title></circle>')
    svg = '<svg viewBox="0 0 1120 365" role="img" aria-label="Primary threshold genome-wide observed over expected"><text x="16" y="20">O/E</text>' + ''.join(elements) + '<text x="560" y="359" text-anchor="middle">Autosome · cumulative scanned coordinate</text></svg>'
    columns = ["window_id", "threshold", "observed_multi", "expected_multi", "oe", "ci95_low", "q_value", "mask_status"]
    row_html = []
    for row in sorted(nominations, key=lambda r: float(r["oe"]), reverse=True):
        cells = []
        for k in columns:
            value = row.get(k, "NA")
            number = finite(value)
            if k in {"expected_multi", "oe", "ci95_low"} and number is not None:
                value = f"{number:.3f}"
            elif k == "q_value" and number is not None:
                value = f"{number:.3g}"
            cells.append('<td>' + html.escape(str(value)) + '</td>')
        row_html.append('<tr>' + ''.join(cells) + '</tr>')
    table = '<table id="hits"><thead><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in columns)+'</tr></thead><tbody>'+''.join(row_html)+'</tbody></table>'
    stage_html = ''.join(f'<tr><td>{html.escape(name)}</td><td>{html.escape(s["status"])}</td><td>{html.escape(s.get("reason", "receipt verified; reused" if s.get("reused") else "completed and receipted"))}</td></tr>' for name,s in stages.items())
    limitations = ''.join('<li>'+html.escape(x)+'</li>' for x in LIMITATIONS)
    parameters = html.escape(json.dumps({k:getattr(args,k) for k in ("window_sizes","step_bp","thresholds","primary_threshold","min_sites","min_observed","min_expected","min_support","min_oe","fdr","site_sampling","seed")}, indent=2))
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Multi-ALT hotspot screening</title><style>
    :root{font:16px/1.6 system-ui,sans-serif;color:#192c3e;background:#f3f6fa}body{margin:0}main{max-width:1250px;margin:auto;padding:40px 24px}h1{font-size:38px;line-height:1.15}h2{margin-top:0}header p{max-width:950px;color:#435a6b}.eyebrow{font-weight:750;letter-spacing:.12em;color:#247b91}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}.card,section{background:white;border:1px solid #dfe7ef;border-radius:12px;padding:24px;margin:22px 0}.card strong{display:block;font-size:32px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:10px;border-bottom:1px solid #e4e9ee;text-align:left}th{background:#f1f5f8}svg{width:100%;font:12px system-ui}input{padding:10px;border:1px solid #bccbd5;border-radius:6px;width:300px;max-width:90%;margin:0 0 16px}pre{overflow:auto;font-size:13px}a{color:#14657e}.scroll{overflow:auto}.note{border-left:4px solid #ec9e4b;padding:8px 16px;background:#fff7eb}@media(max-width:700px){.cards{grid-template-columns:1fr}.card{margin:0}h1{font-size:29px}main{padding:20px 12px}}
    </style><main><header><div class="eyebrow">ALLELE3 / FINAL 7.0</div><h1>Multi-ALT hotspot screening</h1><p>Frozen genome-first screening, post-discovery matched characterization and exact optional molecular/ancient evidence. This report shows computed outputs from this run.</p></header>'''
    document += f'<div class="cards"><div class="card"><strong>{len(all_rows):,}</strong>window × AC-tier rows</div><div class="card"><strong>{len(nominations):,}</strong>working-model nominations</div><div class="card"><strong>{max([float(r["oe"]) for r in nominations],default=0):.2f}×</strong>highest nominated O/E</div></div>'
    document += '<section><h2>Discovery landscape</h2><p>Primary AC tier; points show source-resolution windows with a finite O/E. Inference availability is recorded in the tables. Hover for identity. Red points meet the nomination gates. The line at 2× is an effect-size reference.</p>'+svg+'</section>'
    document += '<section><h2>Nominated windows</h2><p class="note">These are working-model candidates. A missing technical mask means not assessed. Confidence intervals are conditional on fitted expected values and residual dispersion.</p><input id="filter" placeholder="Filter chromosome, window or tier" aria-label="Filter nominated windows"><div class="scroll">'+table+'</div>' + ('<p>No windows met every gate. Empty nominations are a valid result.</p>' if not nominations else '') + '</section>'
    document += '<section><h2>Stage availability</h2><table><thead><tr><th>Stage</th><th>Status</th><th>Evidence</th></tr></thead><tbody>'+stage_html+'</tbody></table><p>Atlas target panel: '+html.escape(json.dumps(atlas_panel or {"status":"not_available"},ensure_ascii=False))+'</p></section>'
    document += '<section><h2>Fixed analysis settings</h2><pre>'+parameters+'</pre><p>Source mode: '+html.escape(window_manifest.get('mode','unknown'))+'. The complete testing family, support audits and model assumptions are recorded in the window manifest.</p></section>'
    document += '<section><h2>Interpretation limits</h2><ul>'+limitations+'</ul></section><section><h2>Reviewable outputs</h2><p><a href="../results/windows.nominated_windows.tsv">Nominated windows</a> · <a href="../results/windows.all_windows.tsv.gz">All windows</a> · <a href="../results/windows.calibration.tsv">Calibration</a> · <a href="../results/characterize.window_matches.tsv">Matched controls</a> · <a href="../results/characterize.site_sampling_audit.tsv">Target sampling</a> · <a href="../qc/manifest.json">Run manifest</a></p></section></main><script>document.getElementById("filter").addEventListener("input",function(){const v=this.value.toLowerCase();document.querySelectorAll("#hits tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(v));});</script></html>'
    (out / "reports").mkdir(parents=True, exist_ok=True)
    (out / "reports/hotspots.html").write_text(document, encoding="utf-8")


def execute(args, out, signature):
    runner = StageRunner(out, signature)
    windows = ["--build", "GRCh38", "--window-sizes", args.window_sizes, "--step-bp", str(args.step_bp),
               "--thresholds", args.thresholds, "--primary-threshold", str(args.primary_threshold)]
    windows += ["--summary", args.summary] if args.summary else ["--store", args.scan_dir]
    for key in ("min_sites", "min_observed", "min_expected", "min_support", "min_training_sites", "min_calibration_windows", "min_oe", "fdr"):
        windows += ["--" + key.replace("_", "-"), str(getattr(args,key))]
    for path in args.exclude_bed:
        windows += ["--exclude-bed", path]
    window_manifest = runner.run("windows", "genome_windows", windows)
    sites = auto_site_tables(args)
    char_args = characterize_args(args, out / "runtime/windows/all_windows.tsv.gz", sites)
    runner.run("characterize", "hotspot_characterize", char_args)
    targets = out / "runtime/characterize/targets.tsv"
    n_targets = count_rows(targets)
    atlas_panel = dict(status="not_available", internally_sampled_target_sites=n_targets)
    atlas_variants, atlas_metrics, atlas_targets = args.atlas_variants, args.atlas_metrics, args.atlas_targets
    if args.with_atlas == "FALSE":
        runner.unavailable("atlas", "Disabled by request", "not_run")
    elif atlas_variants:
        if not atlas_targets:
            raise ValueError("External --atlas-variants requires its exact --atlas-targets table")
        n_external = count_rows(atlas_targets, args.max_targets)
        atlas_panel = dict(status="provided", source="external_supplied_panel_ascertainment_unverified",
                           target_sites=n_external, internally_sampled_target_sites=n_targets,
                           target_table=str(Path(atlas_targets).resolve()))
        runner.unavailable("atlas", f"Using supplied bounded derivatives for {n_external} external target sites; internal frozen panel has {n_targets} sites", "not_run")
    elif not n_targets:
        if args.with_atlas == "TRUE":
            raise ValueError("Atlas requested but the frozen target table is empty; provide an eligible site cache")
        runner.unavailable("atlas", "No target sites available; window discovery remains complete")
    else:
        import atlas_scores as atlas_scores
        atlas_argv = ["--targets", str(targets), "--build", "GRCh38", "--resources", args.atlas_resources,
                      "--dir-data", args.dir_data, "--max-targets", str(args.max_targets), "--require-all-resources"]
        if args.atlas_source_version:
            atlas_argv += ["--source-version", args.atlas_source_version]
        parsed = atlas_scores.parser().parse_args([*atlas_argv, "--out", str(out / "runtime/atlas")])
        try:
            atlas_scores.load_resources(parsed)
        except ValueError as exc:
            # Only missing local files/indices may defer auto mode. Malformed
            # schemas, assembly and semantic errors always stop the workflow.
            if args.with_atlas != "auto" or not str(exc).startswith(("Required Atlas file missing", "Matching Tabix")):
                raise
            runner.unavailable("atlas", str(exc))
        else:
            runner.run("atlas", "atlas_scores", atlas_argv)
            atlas_variants, atlas_metrics, atlas_targets = str(out / "runtime/atlas/atlas_variants.tsv.gz"), str(out / "runtime/atlas/atlas_metrics.tsv.gz"), str(targets)
            atlas_panel = dict(status="computed", source="internally_frozen_panel", source_sampling=args.site_sampling,
                               target_sites=n_targets, internally_sampled_target_sites=n_targets, target_table=str(targets))
    if atlas_variants and args.with_atlas != "FALSE":
        # The exact frozen target table supplies counts. Do not rescan every
        # source-site cache a second time just to attach bounded Atlas scores.
        impact_argv = characterize_args(args, out / "runtime/windows/all_windows.tsv.gz", []) + ["--atlas-variants", atlas_variants, "--atlas-targets", atlas_targets]
        if atlas_metrics:
            impact_argv += ["--atlas-metrics", atlas_metrics]
        runner.run("impact", "hotspot_characterize", impact_argv)
        # Matching and target selection precede all functional outcomes.
        for filename in ("selected_windows.tsv", "window_matches.tsv"):
            if (out / "runtime/characterize" / filename).read_bytes() != (out / "runtime/impact" / filename).read_bytes():
                raise ValueError("Functional annotation changed the frozen window/target design")
    else:
        runner.unavailable("impact", "No bounded Atlas derivatives available", "not_run")
    selection_keys = ("selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38")
    missing = [str(getattr(args,k)) for k in selection_keys if not Path(getattr(args,k)).is_file()]
    if args.with_selection == "FALSE":
        runner.unavailable("selection", "Disabled by request", "not_run")
    elif not n_targets or missing:
        reason = "No target sites" if not n_targets else "Missing local ancient/reference resources: " + ", ".join(missing)
        if args.with_selection == "TRUE":
            raise ValueError(reason)
        runner.unavailable("selection", reason)
    else:
        selection_argv = ["--catalog", str(out / "runtime/characterize/selection_catalog.tsv"),
                          "--window-membership", str(out / "runtime/characterize/window_membership.tsv"),
                          "--effect-allele", args.selection_effect_allele]
        for key in selection_keys:
            selection_argv += ["--" + key.replace("_", "-"), str(getattr(args,key))]
        runner.run("selection", "window_selection", selection_argv)
    report(out, args, runner.results, window_manifest, atlas_panel)
    from analysis_integrity import publish_hotspots
    publish_hotspots(out)
    return dict(status="complete", version=VERSION, stages=runner.results, window_counts=window_manifest.get("counts", {}),
                n_targets=n_targets, atlas_target_panel=atlas_panel, report=str(out / "reports/hotspots.html"), limitations=LIMITATIONS)


def run(args):
    if not args.summary and not args.scan_dir:
        args.scan_dir = os.environ.get("ALLELE_SCAN_STORE", configured_path("ANALYSIS_ROOT", "allele3/gnomad_scan"))
    if args.scan_dir:
        from analysis_integrity import resolve_result_input
        args.scan_dir = str(resolve_result_input(args.scan_dir))
    paths, protected = source_paths(args, metadata=not args.preflight)
    out = check_output_paths(args.out, protected)
    required = [Path(args.summary)] if args.summary else [Path(args.scan_dir) / "run.json"]
    required += [Path(p) for p in args.site_table + args.exclude_bed]
    required += [Path(s.split("=",1)[1]) for s in args.annotation_bed if "=" in s]
    missing = [str(p) for p in required if not p.is_file()]
    if args.preflight:
        return dict(status="planned", mode="metadata_only_no_genomic_reads", missing_required=missing,
                    inputs=[identity(p) for p in paths], out=str(out), limitations=LIMITATIONS)
    if missing:
        raise ValueError("Required inputs missing: " + ", ".join(missing))
    inputs = [identity(p) for p in paths]
    request = dict(version=VERSION, parameters={k:v for k,v in vars(args).items() if k not in {"preflight", "out"}},
                   inputs=inputs, code=code_fingerprint())
    signature = digest(request)
    out.mkdir(parents=True, exist_ok=True)
    (out / "qc/locks").mkdir(parents=True, exist_ok=True)
    with (out / "qc/locks/pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Pipeline output is in use: {out}") from exc
        manifest = out / "qc/manifest.json"
        if manifest.is_file():
            old = json.loads(manifest.read_text())
            if old.get("request_signature") != signature:
                raise ValueError("Input, settings or code changed; use a new --out directory to preserve the frozen discovery design")
            if old.get("status") == "complete":
                return validate_manifest(out, "qc/manifest.json")
        elif out.exists() and any(p.name not in {"jobs", "logs", "qc"} for p in out.iterdir()):
            raise ValueError("Unreceipted files exist in pipeline output; use a new --out directory")
        out.mkdir(parents=True, exist_ok=True)
        started = dict(status="running", request_signature=signature, request=request,
                       started_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        atomic_json(manifest, started)
        try:
            result = execute(args, out, signature)
            if [identity(p) for p in paths] != inputs or code_fingerprint() != request["code"]:
                raise ValueError("Scientific source files or code changed during the pipeline")
            # The shared scheduler validates this top-level input dictionary.
            # Include external child dependencies as well as optional missing
            # resources so source changes cannot reuse a stale completed job.
            scientific_inputs = {item["path"]:item for item in inputs}
            for stage in result.get("stages", {}).values():
                if stage.get("manifest"):
                    child = json.loads(Path(stage["manifest"]).read_text())
                    for path, saved in child.get("inputs", {}).items():
                        if not Path(path).resolve().is_relative_to(out):
                            scientific_inputs[path] = saved
            result.update(inventory_base="module", request_signature=signature, request=request, inputs=scientific_inputs, started_utc=started["started_utc"],
                          completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(), output_inventory=inventory(out,"qc/manifest.json"))
            atomic_json(manifest, result)
            return result
        except BaseException as exc:
            atomic_json(manifest, dict(started, status="failed", error=str(exc),
                                      recovery="Completed child stages are preserved. Rerun the same command to resume; changed inputs/settings/code require a new output."))
            raise


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group()
    source.add_argument("--scan-dir", "--store", dest="scan_dir", help="Completed immutable scan store; default ALLELE_SCAN_STORE")
    source.add_argument("--summary", help="Legacy window_summary.tsv; preserves its existing expected model and resolution")
    p.add_argument("--out", required=True)
    p.add_argument("--build", choices=("GRCh38", "38", "hg38"), default="GRCh38")
    p.add_argument("--window-sizes", default="1000000")
    p.add_argument("--step-bp", type=int, default=1000000)
    p.add_argument("--thresholds", default="5,10")
    p.add_argument("--primary-threshold", type=int, choices=(5,10), default=5)
    for key, kind, default in (("min-sites",int,500),("min-observed",int,25),("min-expected",float,10.),
                                ("min-support",float,.99),("min-training-sites",int,100),
                                ("min-calibration-windows",int,50),("min-oe",float,2.),("fdr",float,.05)):
        p.add_argument("--"+key, type=kind, default=default)
    p.add_argument("--exclude-bed", action="append", default=[], help="BED0 whole-window nomination veto; repeatable")
    p.add_argument("--annotation-bed", action="append", default=[], metavar="NAME=PATH", help="Post-discovery sorted GRCh38 BED4; repeatable")
    p.add_argument("--characterize-window-bp", type=int, help="One scale for matched characterization; required for multiscale discovery")
    p.add_argument("--covariates", help="Prespecified window_id-keyed matching covariates")
    p.add_argument("--match-columns", default="n_eligible")
    p.add_argument("--controls-per-case", type=int, default=2)
    p.add_argument("--max-match-distance", type=float, default=1.)
    p.add_argument("--control-gap-bp", type=int, default=0)
    p.add_argument("--site-table", action="append", default=[], help="Sorted full ALT AC/AN cache; auto uses scan multiallelic caches if available")
    p.add_argument("--site-sampling", choices=("selected_multialt_candidates", "all_eligible_source_sites", "prespecified_sites"), default="selected_multialt_candidates")
    p.add_argument("--sites-per-window", type=int, default=25)
    p.add_argument("--max-targets", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20261106)
    p.add_argument("--with-atlas", choices=("auto", "TRUE", "FALSE"), default="auto")
    p.add_argument("--atlas-resources", help="Auto prefers indexed per-chromosome VCF; legacy all-three TSV manifest remains supported")
    p.add_argument("--atlas-source-version", help="Explicit downloaded release identity; never inferred from a filename")
    p.add_argument("--dir-data", default=os.environ.get("DIR_DATA", configured_path("DIR_BIG")))
    p.add_argument("--atlas-variants", help="Optional already-computed bounded derivative")
    p.add_argument("--atlas-metrics")
    p.add_argument("--atlas-targets", help="Exact counts-bearing target table for supplied derivatives")
    p.add_argument("--with-selection", choices=("auto", "TRUE", "FALSE"), default="auto")
    ancient = os.environ.get("ALLELE_ANCIENT_DIR", configured_path("DIR_GEN", "ancient"))
    p.add_argument("--selection-stats", default=str(Path(ancient) / "Selection_Summary_Statistics_01OCT2025.tsv.gz"))
    p.add_argument("--selection-effect-allele", choices=("raw_only", "ALT"), default="raw_only")
    for key, env, default in (("fasta38","ALLELE_FASTA",configured_path("DIR_GEN","1hgp/38/GRCH38.fasta.gz")),
                              ("fasta37","ALLELE_FASTA37",configured_path("DIR_GEN","1hgp/37/GRCH37.fasta.gz")),
                              ("chain38to37","ALLELE_CHAIN38TO37",configured_path("DIR_FILES","liftOver/hg38ToHg19.over.chain.gz")),
                              ("chain37to38","ALLELE_CHAIN37TO38",configured_path("DIR_FILES","liftOver/hg19ToHg38.over.chain.gz"))):
        p.add_argument("--"+key, default=os.environ.get(env,default))
    p.add_argument("--preflight", action="store_true", help="Read-only metadata inventory; no genomic reads or scientific computation")
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit("hotspots: " + str(exc))
    print(json.dumps({k:result[k] for k in ("status","version","window_counts","n_targets","report","missing_required") if k in result}, sort_keys=True))
    return 0


# Hotspot Characterize
# ======================================================================
_description_hotspot_characterize = 'Characterize frozen genome-window nominations without using HERV to find them.\n\nThe input is the COMPLETE all_windows.tsv.gz from genome_windows, including\nnominated=0 rows. Nominations are never changed here. Overlapping nominated\nwindows are reduced to a deterministic coordinate-first representative for\nmatching, while all nominations remain in the annotation output. Controls are\nnon-nominated, tested, nonoverlapping windows from the same chromosome, window\nsize, actual span and AC threshold. Matching uses only declared pre-discovery\ncovariates; annotation, molecular scores, observed multi-ALT counts and scan\nstatistics never enter the distance or tie-breaker.\n\nAnnotation input is sorted BED4 (0-based half-open) supplied as NAME=PATH. The\nfourth field is a category, e.g. ERV1, ERVK, exon, promoter. Categories can\noverlap: each category and __ANY__ are union-counted independently, and category\nfractions must not be added. This is post-discovery descriptive characterization,\nnot an independent replication or a causal enrichment/selection test.\n\nOptional --site-table inputs are streamed and sampled by a fixed hash without\nlooking at allele counts, Atlas scores, HERV labels or ancient outcomes. Supply\nthe sampling frame explicitly. Native multiallelic_sites.tsv.gz caches are\naccepted, but their restricted sampling cannot measure whole-window tolerance.\ntargets.tsv is directly accepted by base_occupancy atlas. selection_catalog.tsv is\naccepted by the published-selection overlay. ancient_target_requests.tsv is a\nrequest list: a WINDOW match is not a valid one-to-one ancient SITE pair.\n\nOptional bounded Atlas derivatives summarize complete three-ALT vectors only;\nmissing values remain missing. AVI and its signed SHAP attributions are trained\nusing gnomAD-frequency proxy labels. These profiles are not pathogenicity\nprobabilities, molecular pleiotropy counts, tolerance or positive selection.\n'

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
import heapq
import re
import statistics

from analysis_integrity import transactional

CHARACTERIZE_VERSION = "hotspot-characterization-1.0"
CHARACTERIZE_BASES = set("ACGT")
CHARACTERIZE_NA = {"", ".", "NA", "NAN", "NONE", "NULL"}
CHARACTERIZE_WINDOW_REQUIRED = {"window_id", "chrom", "start0", "end0", "build", "window_bp",
                   "threshold", "n_eligible", "observed_multi", "expected_multi", "oe",
                   "p_value", "q_value", "nominated", "inference_status"}
CHARACTERIZE_ALLOWED_MATCH_COLUMNS = {
    "n_eligible", "gc_fraction", "cpg_fraction", "mappability_mean", "callable_fraction",
    "mean_an", "median_an", "recombination_rate", "context_entropy",
}
CHARACTERIZE_ATTRIBUTIONS = {
    "MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",
    "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",
    "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",
    "IS_INSERTION", "IS_DELETION",
}
CHARACTERIZE_LIMITS = [
    "Nominations and their multiple-testing family are inherited from the full window scan; this module does not optimize thresholds or rerun discovery using annotations.",
    "Matched controls are an exposure-defined, non-nominated comparison in the discovery dataset, not independent replication. Matching covariates, calipers, chromosome constraints and unmatched cases must be reported.",
    "Default matching uses log1p(n_eligible) only. It does not establish balance for sequence context, recombination, callability or other unprovided covariates; none are invented.",
    "Controls respect the same whole-window exclusion-mask veto as nominations. A not_assessed mask status remains unassessed; zero overlap is not a certification of technical reliability.",
    "One coordinate-first representative per overlapping nomination component is used for matching. Remaining nominated windows are annotated but do not become additional independent matched cases.",
    "BED annotations are declared GRCh38 by the caller; coordinate syntax and interval arithmetic are checked, not the underlying annotation's biological accuracy or source completeness.",
    "Annotation fractions use union overlap within each category and within each source's __ANY__ category. Different categories may overlap and must not be summed.",
    "Annotation contrasts and Atlas profiles are descriptive after genome-wide discovery. No annotation P value, causal tolerance estimate or selection coefficient is produced here.",
    "Optional target sampling is hash-based within selected windows and records its frame and coverage. A multiallelic-only input remains a multiallelic-only panel; no absent source ALT is imputed to AC=0, and REF presence is not inferred.",
    "Every selected case/control window receives the same preallocated maximum target quota; a global-cap remainder is deliberately unused. Windows with fewer provided source sites retain all of them and expose the coverage difference. Source-ALT AF order is not a joint-sample four-base rank or mutation order.",
    "AVI and SHAP use gnomAD-frequency-derived training labels. A lower AVI or fewer high-rank predictions in selected windows is not independent evidence of tolerance, phenotypic pleiotropy or balancing selection.",
    "Complete three-ALT score vectors provide a common per-site theoretical substitution denominator only for represented target sites, not for every genomic base in a window. Sparse splicing and unavailable scores are never set to zero.",
    "SHAP sign denotes a contribution to AVI, not expression up/down. Nonzero attribution columns do not count biological functions.",
    "The ancient request table is not a paired-site catalog. The selection catalog supports exact published nucleotide-contrast lookup only; a published biallelic signal must not be assigned to another third base.",
]


def characterize_canonical_chrom(value):
    value = re.sub(r"^chr", "", str(value).strip(), flags=re.I).upper()
    if value == "M": value = "MT"
    if value not in {*(str(i) for i in range(1, 23)), "X", "Y", "MT"}:
        raise ValueError(f"Noncanonical chromosome: {value!r}")
    return "chrM" if value == "MT" else "chr" + value


def characterize_chrom_order(value):
    c = characterize_canonical_chrom(value)[3:]
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "M": 25}[c]


def characterize_numeric(value, label, *, integer=False, optional=False):
    text = str(value).strip()
    if text.upper() in CHARACTERIZE_NA:
        if optional: return None
        raise ValueError(f"Missing {label}")
    try: result = float(text)
    except ValueError as exc: raise ValueError(f"Non-numeric {label}: {text!r}") from exc
    if not math.isfinite(result) or integer and (result < 0 or result > 2**53 or result != int(result)):
        raise ValueError(f"Invalid {label}: {text!r}")
    return int(result) if integer else result


def characterize_build38(value):
    if str(value).strip().lower() not in {"grch38", "38", "hg38"}:
        raise ValueError("All inputs must explicitly use GRCh38/hg38; no liftover or relabelling")
    return "GRCh38"


def characterize_read_tsv(path, required=()):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        def lines():
            for line in handle:
                if not line.strip() or line.startswith("##"): continue
                yield line
        reader = csv.DictReader(lines(), delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)) or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"{path}: need unique TSV columns {','.join(sorted(required))}")
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path}:{line}: malformed TSV row")
            yield {k: v.strip() for k, v in row.items()}


def characterize_write_tsv(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "NA") for k in fields})


def characterize_annotation_specs(values):
    specs = []
    for text in values or []:
        if "=" not in text: raise ValueError("--annotation-bed requires NAME=PATH")
        name, path = text.split("=", 1)
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", name) or not path:
            raise ValueError("Annotation names must be simple distinct labels with a nonempty local path")
        if name in {s[0] for s in specs}: raise ValueError(f"Duplicate annotation source label: {name}")
        specs.append((name, Path(path).expanduser().resolve()))
    return specs


def characterize_protected_sources(raw):
    result = [InputPath(p) for _, p in characterize_annotation_specs(raw.get("annotation_bed"))]
    result += [InputPath(p) for p in raw.get("site_table", []) or []]
    return result


def characterize_matching_columns(value):
    names = [x.strip() for x in value.split(",")]
    if not names or len(set(names)) != len(names) or any(x not in CHARACTERIZE_ALLOWED_MATCH_COLUMNS for x in names):
        raise ValueError("Match columns must be distinct pre-discovery covariates from: " + ",".join(sorted(CHARACTERIZE_ALLOWED_MATCH_COLUMNS)))
    return names


def characterize_load_windows(path, threshold=5, window_bp=None, maximum=500000, covariates=None):
    rows, ids, locations, primary_sizes, observed_sizes = [], set(), set(), set(), set()
    for row in characterize_read_tsv(path, CHARACTERIZE_WINDOW_REQUIRED):
        characterize_build38(row["build"])
        sid = row["window_id"]
        if sid.upper() in CHARACTERIZE_NA or sid in ids: raise ValueError("Window IDs must be unique and nonmissing in the full scan table")
        ids.add(sid)
        if len(ids) > maximum: raise ValueError(f"--max-windows={maximum} exceeded")
        row["chrom"] = characterize_canonical_chrom(row["chrom"])
        for k in ("start0", "end0", "window_bp", "threshold", "n_eligible", "nominated"):
            row[k] = characterize_numeric(row[k], k, integer=True)
        if row["end0"] <= row["start0"] or row["window_bp"] < row["end0"] - row["start0"]:
            raise ValueError("Windows must be 0-based half-open with a positive span <= window_bp")
        if row["threshold"] not in (5, 10) or row["nominated"] not in (0, 1):
            raise ValueError("Require AC5/AC10 and binary frozen nominations")
        key = row["chrom"], row["start0"], row["end0"], row["window_bp"], row["threshold"]
        if key in locations: raise ValueError("Duplicate window coordinates/scale/threshold")
        locations.add(key)
        row["span_bp"] = row["end0"] - row["start0"]
        row["mask_overlap_bp"] = characterize_numeric(row.get("mask_overlap_bp", "0"), "mask_overlap_bp", integer=True)
        row["mask_status"] = row.get("mask_status", "not_assessed")
        if row["mask_status"] not in {"not_assessed", "no_overlap_with_supplied_masks", "overlapping_exclusion"}:
            raise ValueError("Unknown window exclusion-mask status")
        if row["mask_overlap_bp"] > row["span_bp"] or (row["mask_overlap_bp"] > 0) != (row["mask_status"] == "overlapping_exclusion"):
            raise ValueError("Window mask status and overlap are inconsistent")
        if row["nominated"] and row["mask_overlap_bp"]:
            raise ValueError("A frozen nomination cannot overlap an exclusion-mask veto")
        for k in ("observed_multi", "expected_multi", "oe", "p_value", "q_value"):
            row[k] = characterize_numeric(row[k], k, optional=True, integer=k == "observed_multi")
        if row["observed_multi"] is not None and row["observed_multi"] > row["n_eligible"]:
            raise ValueError("Observed multi-ALT sites cannot exceed eligible sites")
        if any(row[k] is not None and not 0 <= row[k] <= 1 for k in ("p_value", "q_value")):
            raise ValueError("Window P/Q must lie within [0,1]")
        if row["expected_multi"] is not None and not 0 <= row["expected_multi"] <= row["n_eligible"] + 1e-8:
            raise ValueError("Expected multi-ALT sites must lie within the eligible denominator")
        if row["observed_multi"] is not None and row["expected_multi"] and row["oe"] is not None and not math.isclose(row["oe"], row["observed_multi"]/row["expected_multi"], rel_tol=1e-4, abs_tol=1e-7):
            raise ValueError("Window O/E disagrees with its observed and expected counts")
        if row["nominated"] and (row["oe"] is None or row["oe"] < 2 or not row["expected_multi"] or row["observed_multi"] is None or row["p_value"] is None or row["q_value"] is None):
            raise ValueError("A frozen >=2x nomination needs a finite >=2 O/E and tested P/Q")
        row["tested"] = row["p_value"] is not None and row["q_value"] is not None
        if row.get("is_primary", "0") not in {"0", "1"}: raise ValueError("is_primary must be binary when supplied")
        if row.get("is_primary") == "1": primary_sizes.add(row["window_bp"])
        if row["threshold"] == threshold:
            observed_sizes.add(row["window_bp"])
            rows.append(row)
    if window_bp is None:
        choices = primary_sizes & observed_sizes if primary_sizes else observed_sizes
        if len(choices) != 1: raise ValueError("Supply --window-bp to select one scan scale; multiple scales must not be pooled")
        window_bp = next(iter(choices))
    rows = [r for r in rows if r["window_bp"] == window_bp]
    if not rows: raise ValueError("No windows match the requested threshold and scale")
    if covariates:
        lookup, seen = {r["window_id"]: r for r in rows}, set()
        for c in characterize_read_tsv(covariates, {"window_id"}):
            wid = c["window_id"]
            if wid in seen: raise ValueError("Covariate table contains duplicate window IDs")
            seen.add(wid)
            for name in set(c) - {"window_id"}:
                if name not in CHARACTERIZE_ALLOWED_MATCH_COLUMNS:
                    raise ValueError(f"Unapproved matching covariate column: {name}")
            if wid in lookup:
                for name, value in c.items():
                    if name != "window_id":
                        if name == "n_eligible" and characterize_numeric(value, name, integer=True) != lookup[wid][name]:
                            raise ValueError("Covariate table cannot replace the scan's eligible denominator")
                        lookup[wid][name] = value
    rows.sort(key=lambda r: (characterize_chrom_order(r["chrom"]), r["start0"], r["end0"], r["window_id"]))
    return rows, dict(input_window_rows=len(ids), selected_threshold=threshold, selected_window_bp=window_bp,
                      selected_window_rows=len(rows), source_primary_sizes=sorted(primary_sizes))


class characterize_WindowIndex:
    def __init__(self, windows):
        groups = defaultdict(list)
        for row in windows: groups[row["chrom"]].append(row)
        self.groups = {}
        for chrom, rows in groups.items():
            rows.sort(key=lambda r: (r["start0"], r["end0"], r["window_id"]))
            ends, last = [], -1
            for row in rows:
                last = max(last, row["end0"]); ends.append(last)
            self.groups[chrom] = rows, [r["start0"] for r in rows], ends

    def overlaps(self, chrom, start, end):
        rows, starts, ends = self.groups.get(characterize_canonical_chrom(chrom), ([], [], []))
        lo, hi = bisect_right(ends, start), bisect_left(starts, end)
        return [r for r in rows[lo:hi] if r["end0"] > start]


def characterize_nomination_representatives(windows):
    nominated = [r for r in windows if r["nominated"]]
    representatives, leads = [], {}
    last_chrom, last_end, lead = None, -1, None
    for row in nominated:
        if row["chrom"] != last_chrom or row["start0"] >= last_end:
            representatives.append(row); lead = row["window_id"]
            last_end = row["end0"]
        else:
            last_end = max(last_end, row["end0"])
        last_chrom = row["chrom"]
        leads[row["window_id"]] = lead
    return representatives, leads


def characterize_match_windows(windows, columns, controls_per_case=2, max_distance=1., gap_bp=0, seed=20261106, max_cases=2000):
    """Greedy nearest controls; all distances/tie-breaks are outcome-blind."""
    import numpy as np
    from scipy.spatial import cKDTree

    representatives, leads = characterize_nomination_representatives(windows)
    if len(representatives) > max_cases: raise ValueError(f"--max-cases={max_cases} exceeded; explicitly set a larger bounded follow-up")
    nominated_index = characterize_WindowIndex([r for r in windows if r["nominated"]])
    values, all_values = {}, []
    for row in windows:
        vector = [characterize_numeric(row.get(k, "NA"), k, optional=True) for k in columns]
        if any(v is None for v in vector): continue
        if any(v < 0 for k, v in zip(columns, vector) if k in {"n_eligible", "mean_an", "median_an", "recombination_rate", "context_entropy"}):
            raise ValueError("Matching counts, rates and entropy must be nonnegative")
        if any(not 0 <= v <= 1 for k, v in zip(columns, vector) if k.endswith("_fraction") or k == "mappability_mean"):
            raise ValueError("Matching fractions/mappability_mean must be within [0,1]")
        vector = [math.log1p(v) if k == "n_eligible" else v for k, v in zip(columns, vector)]
        values[row["window_id"]] = vector
        if row["tested"] and not row.get("mask_overlap_bp", 0): all_values.append(vector)
    with np.errstate(over="ignore", invalid="ignore"):
        sd = np.std(np.asarray(all_values, dtype=float), axis=0) if all_values else np.ones(len(columns))
    if not np.isfinite(sd).all(): raise ValueError("Matching covariate scale is non-finite; verify units and values")
    sd = np.where(sd > 1e-12, sd, 1.)
    pools = defaultdict(list)
    def key(row): return row["chrom"], row["window_bp"], row["span_bp"], row["threshold"]
    for row in windows:
        if row["nominated"] or not row["tested"] or row.get("mask_overlap_bp", 0) or row["window_id"] not in values: continue
        if nominated_index.overlaps(row["chrom"], max(0, row["start0"]-gap_bp), row["end0"]+gap_bp): continue
        pools[key(row)].append(row)
    trees = {k: cKDTree(np.asarray([values[r["window_id"]] for r in rows]) / sd) for k, rows in pools.items()}
    used, selected, matches, balance = set(), [], [], []
    used_by_chrom = defaultdict(list)
    for case in representatives:
        cid = case["window_id"]
        match_id = "GW_" + hashlib.sha256((str(seed)+"|"+cid).encode()).hexdigest()[:16]
        selected.append(dict(case, group="case", window_match_id=match_id, match_role="representative_nomination"))
        choices = []
        if cid in values and key(case) in trees:
            cv = np.asarray(values[cid]) / sd
            indices = trees[key(case)].query_ball_point(cv, max_distance * math.sqrt(len(columns)) + 1e-12)
            for i in indices:
                control = pools[key(case)][i]
                wid = control["window_id"]
                if wid in used: continue
                distance = float(np.sqrt(np.mean(((np.asarray(values[wid])-values[cid])/sd)**2)))
                tie = hashlib.sha256((str(seed)+"|"+cid+"|"+wid).encode()).hexdigest()
                choices.append((distance, tie, control))
        chosen = []
        for distance, _, control in sorted(choices, key=lambda x: (x[0], x[1])):
            if any(control["start0"] < end+gap_bp and control["end0"]+gap_bp > start for start, end in used_by_chrom[control["chrom"]]): continue
            used.add(control["window_id"])
            used_by_chrom[control["chrom"]].append((control["start0"], control["end0"]))
            chosen.append((control, distance))
            selected.append(dict(control, group="control", window_match_id=match_id, match_role="non_nominated_matched_control"))
            for j, name in enumerate(columns):
                a, b = values[cid][j], values[control["window_id"]][j]
                balance.append(dict(window_match_id=match_id, case_window_id=cid, control_window_id=control["window_id"],
                                    covariate=name, transformation="log1p" if name == "n_eligible" else "identity",
                                    case_value=a, control_value=b, pooled_tested_sd=float(sd[j]),
                                    standardized_difference=(a-b)/float(sd[j])))
            if len(chosen) == controls_per_case: break
        status = ("missing_matching_covariates" if cid not in values else "unmatched" if not chosen else
                  "matched" if len(chosen) == controls_per_case else "partial_match")
        if chosen:
            for rank, (control, distance) in enumerate(chosen, 1):
                matches.append(dict(window_match_id=match_id, case_window_id=cid, control_window_id=control["window_id"],
                                    control_rank=rank, matching_distance=distance, n_controls=len(chosen), status=status,
                                    representative_window_id=cid))
        else:
            matches.append(dict(window_match_id=match_id, case_window_id=cid, control_window_id=None, control_rank=None,
                                matching_distance=None, n_controls=0, status=status, representative_window_id=cid))
    for row in windows:
        if row["nominated"] and leads[row["window_id"]] != row["window_id"]:
            matches.append(dict(window_match_id=None, case_window_id=row["window_id"], control_window_id=None,
                                control_rank=None, matching_distance=None, n_controls=0,
                                status="overlapping_nomination_not_independent_matched_case", representative_window_id=leads[row["window_id"]]))
    info = dict(nominated_windows=len(leads), representative_cases=len(representatives), matched_controls=len(used),
                fully_matched_cases=sum(r["status"] == "matched" and r["control_rank"] == 1 for r in matches),
                tested_mask_veto_windows=sum(r["tested"] and bool(r.get("mask_overlap_bp", 0)) for r in windows),
                algorithm="greedy_coordinate_order_nearest_caliper_without_replacement_or_control_overlap",
                columns=columns, scales={name: float(v) for name, v in zip(columns, sd)},
                scaling_population="all tested, mask-eligible windows with complete declared covariates before nomination/control filtering",
                scope="descriptive_discovery_dataset_comparison_not_independent_replication")
    return selected, matches, balance, leads, info


@dataclass
class characterize_UnionCount:
    total: int = 0
    last_end: int = -1
    records: int = 0
    def add(self, start, end):
        self.total += max(0, end-max(start, self.last_end))
        self.last_end = max(self.last_end, end)
        self.records += 1


def characterize_annotate_windows(windows, selected, leads, specs):
    selected_map = {r["window_id"]: r for r in selected}
    annotated = [dict(r, group="case", match_role="secondary_overlapping_nomination", window_match_id=None)
                 for r in windows if r["nominated"] and r["window_id"] not in selected_map]
    annotated += selected
    index = characterize_WindowIndex(annotated)
    output, audits = [], []
    for name, path in specs:
        accum, categories, audit = {}, set(), Counter()
        last_chrom, last_start, seen_chrom = None, -1, set()
        if annotated:
            opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
            with opener(path, "rt", encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith(("#", "track ", "browser ")): continue
                    cells = line.rstrip("\r\n").split("\t")
                    if len(cells) < 4: raise ValueError(f"{path}:{line_number}: need sorted BED4 with category")
                    chrom = characterize_canonical_chrom(cells[0]); start = characterize_numeric(cells[1], "BED start", integer=True); end = characterize_numeric(cells[2], "BED end", integer=True)
                    category = cells[3].strip()
                    if end <= start or category.upper() in CHARACTERIZE_NA or category == "__ANY__": raise ValueError("BED categories must be explicit and intervals positive; __ANY__ is reserved")
                    if chrom != last_chrom:
                        if chrom in seen_chrom: raise ValueError("BED chromosomes must be contiguous")
                        seen_chrom.add(chrom); last_start = -1
                    if start < last_start: raise ValueError("BED records must be sorted by start within each chromosome; no implicit whole-file sorting")
                    last_chrom, last_start = chrom, start
                    categories.add(category); audit["source_records"] += 1
                    for w in index.overlaps(chrom, start, end):
                        lo, hi = max(start, w["start0"]), min(end, w["end0"])
                        for cat in (category, "__ANY__"):
                            accum.setdefault((w["window_id"], cat), characterize_UnionCount()).add(lo, hi)
                        audit["record_window_overlaps"] += 1
        # Sparse positive category rows plus one explicit __ANY__ row/window
        # avoid an O(windows * genes) zero table for gene-level BED annotations.
        # For a supplied source, an omitted category/window combination has 0 bp.
        categories_by_window = defaultdict(set)
        for wid, cat in accum: categories_by_window[wid].add(cat)
        for w in annotated:
            observed_categories = categories_by_window[w["window_id"]]
            for category in sorted(observed_categories | {"__ANY__"}):
                state = accum.get((w["window_id"], category), characterize_UnionCount())
                assert 0 <= state.total <= w["span_bp"]
                output.append(dict(window_id=w["window_id"], chrom=w["chrom"], start0=w["start0"], end0=w["end0"],
                                   group=w["group"], window_match_id=w["window_match_id"], match_role=w["match_role"],
                                   representative_window_id=leads.get(w["window_id"], w["window_id"]),
                                   annotation_source=name, category=category, union_overlap_bp=state.total,
                                   window_bp=w["span_bp"], overlap_fraction=state.total/w["span_bp"], overlapping_records=state.records,
                                   status="post_discovery_descriptive_union_overlap"))
        audits.append(dict(annotation_source=name, source=identity(path, full=path.stat().st_size <= 64*1024*1024),
                           categories=len(categories), counts=dict(audit),
                           output_encoding="positive category overlaps plus explicit __ANY__ for every annotated window; absent category/window = zero overlap",
                           status="complete" if annotated else "not_read_no_selected_windows"))
    return output, audits


def characterize_annotation_contrasts(annotations, matches):
    lookup = {(r["window_id"], r["annotation_source"], r["category"]): r for r in annotations}
    labels_by_window = defaultdict(set)
    for r in annotations: labels_by_window[r["window_id"]].add((r["annotation_source"], r["category"]))
    pairs = defaultdict(list)
    for m in matches:
        if m["control_window_id"]: pairs[(m["window_match_id"], m["case_window_id"])].append(m["control_window_id"])
    result = []
    for (mid, case), controls in pairs.items():
        labels = sorted(set().union(*(labels_by_window[w] for w in [case, *controls])))
        for source, category in labels:
            cv = lookup.get((case, source, category), {}).get("overlap_fraction", 0.)
            values = [lookup.get((c, source, category), {}).get("overlap_fraction", 0.) for c in controls]
            mean = statistics.fmean(values)
            result.append(dict(window_match_id=mid, case_window_id=case, annotation_source=source, category=category,
                               n_controls=len(values), case_fraction=cv, mean_control_fraction=mean, difference=cv-mean,
                               ratio=cv/mean if mean > 0 else None, status="descriptive_no_annotation_hypothesis_test"))
    return result


def characterize_parse_site(row):
    for label in ("build", "genome_build"):
        if row.get(label, "").upper() not in CHARACTERIZE_NA: characterize_build38(row[label])
    chrom, pos, ref = characterize_canonical_chrom(row["chrom"]), characterize_numeric(row["pos"], "site POS", integer=True), row["ref"].upper()
    if pos < 1 or ref not in CHARACTERIZE_BASES: raise ValueError("Sites need a positive 1-based POS and canonical REF")
    if {"observed_alts_all", "alt_ac_all", "alt_an_all"} <= set(row):
        alt_text, ac_text, an_text = row["observed_alts_all"], row["alt_ac_all"], row["alt_an_all"]
    elif {"alts_by_af", "alt_ac", "alt_an"} <= set(row):
        alt_text, ac_text, an_text = row["alts_by_af"], row["alt_ac"], row["alt_an"]
    elif {"alts", "all_alt_ac", "all_alt_an"} <= set(row):
        alt_text, ac_text, an_text = row["alts"], row["all_alt_ac"], row["all_alt_an"]
    else:
        raise ValueError("Site tables need aligned reported ALT/AC/AN arrays; coordinate-only all_sites summaries do not establish observed allele counts")
    alts = [x.strip().upper() for x in alt_text.split(",")]
    ac = [characterize_numeric(x, "ALT AC", integer=True) for x in ac_text.split(",")]
    an = [characterize_numeric(x, "ALT AN", integer=True) for x in an_text.split(",")]
    if not 1 <= len(alts) <= 3 or len(set(alts)) != len(alts) or any(a not in CHARACTERIZE_BASES or a == ref for a in alts) or not len(alts) == len(ac) == len(an):
        raise ValueError("Invalid or misaligned reported ALT arrays")
    if any(n < 1 or c > n for c, n in zip(ac, an)): raise ValueError("ALT counts must satisfy 0<=AC<=positive AN")
    sid = row.get("site_id", "")
    if sid.upper() in CHARACTERIZE_NA: sid = f"{chrom}:{pos}:{ref}"
    ranked = sorted(zip(alts, ac, an), key=lambda x: (-x[1]/x[2], x[0]))
    strict = row.get("strict_filter_complete", "0")
    if strict not in {"0", "1"}: raise ValueError("strict_filter_complete must be binary when supplied")
    return dict(chrom=chrom, pos=pos, ref=ref, site_id=sid, alts=",".join(alts), build="GRCh38",
                observed_alts_all=",".join(alts), alt_ac_all=",".join(map(str, ac)), alt_an_all=",".join(map(str, an)),
                all_alt_ac=",".join(map(str, ac)), all_alt_an=",".join(map(str, an)),
                n_alt_ac5=sum(c >= 5 for c in ac), n_alt_ac10=sum(c >= 10 for c in ac),
                primary_alt=ranked[0][0], secondary_alt=ranked[1][0] if len(ranked)>1 else None,
                strict_filter_complete=int(strict), ref_presence_status="not_inferable_from_reported_ALT_counts",
                shared_AN=int(len(set(an)) == 1), _counts=dict(zip(alts, zip(ac, an))))


CHARACTERIZE_TARGET_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "alts", "family", "group", "match_id", "window_id", "window_match_id",
                 "partition", "comparison_partition", "observed_alts_all", "alt_ac_all", "alt_an_all", "all_alt_ac", "all_alt_an",
                 "primary_alt", "secondary_alt", "n_alt_ac5", "n_alt_ac10", "strict_filter_complete", "ref_presence_status", "shared_AN",
                 "source_sampling", "sampling_fraction", "selection_basis", "catalog_scope"]
CHARACTERIZE_REQUEST_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "alts", "window_id", "window_match_id", "group",
                  "observed_alts_all", "alt_ac_all", "alt_an_all", "source_sampling", "request_scope", "site_pairing_status"]


def characterize_export_targets(paths, sampling, selected, out, per_window=25, maximum=5000, seed=20261106):
    if paths and not sampling: raise ValueError("--site-table requires an explicit --site-sampling declaration")
    if not paths: return [], [], dict(status="not_requested", source_sampling=None)
    rows = sorted(selected, key=lambda r: (characterize_chrom_order(r["chrom"]), r["start0"], r["window_id"]))
    if len(rows) > maximum: raise ValueError("--max-targets cannot provide even one possible target per selected window; increase it explicitly")
    base = min(per_window, maximum//len(rows)) if rows else 0
    quotas = {r["window_id"]: base for r in rows}
    index = characterize_WindowIndex(rows); heaps = defaultdict(list); seen = Counter(); counts = Counter()
    last_key = None
    for path in paths if rows else []:
        for raw in characterize_read_tsv(path, {"chrom", "pos", "ref"}):
            chrom, pos = characterize_canonical_chrom(raw["chrom"]), characterize_numeric(raw["pos"], "POS", integer=True)
            key = characterize_chrom_order(chrom), pos
            if last_key is not None and key <= last_key: raise ValueError("Site tables must be globally chromosome/POS sorted with unique positions across files")
            last_key = key; counts["source_rows_streamed"] += 1
            matches = index.overlaps(chrom, pos-1, pos)
            if not matches: continue
            site = characterize_parse_site(raw); counts["rows_in_selected_windows"] += 1
            if len(matches) != 1: raise ValueError("Target export requires nonoverlapping representative/control windows")
            w = matches[0]; wid = w["window_id"]; seen[wid] += 1
            coordinate_identity = f"{site['chrom']}:{site['pos']}:{site['ref']}"
            digest = int(hashlib.sha256((str(seed)+"|"+coordinate_identity).encode()).hexdigest(), 16)
            heap = heaps[wid]; item = (-digest, site["pos"], site)
            if len(heap) < quotas[wid]: heapq.heappush(heap, item)
            elif quotas[wid] and item > heap[0]: heapq.heapreplace(heap, item)
    targets, audit = [], []
    selected_ids = set()
    for w in rows:
        wid = w["window_id"]
        sites = sorted((item[2] for item in heaps[wid]), key=lambda x: x["pos"])
        for site in sites:
            if site["site_id"] in selected_ids: raise ValueError("A site_id identifies duplicate sampled coordinates")
            selected_ids.add(site["site_id"])
            targets.append(dict(site, family="genome_window", group=w["group"], match_id=w["window_match_id"],
                                window_id=wid, window_match_id=w["window_match_id"], partition="genome_window", comparison_partition="genome_window",
                                source_sampling=sampling, sampling_fraction=len(sites)/seen[wid] if seen[wid] else None,
                                selection_basis="fixed_seed_coordinate_identity_hash_not_allele_count_or_score",
                                catalog_scope="sampled_reported_modern_ALT_states_in_frozen_genome_windows"))
        audit.append(dict(window_id=wid, group=w["group"], window_match_id=w["window_match_id"], source_sampling=sampling,
                          source_rows_in_window=seen[wid], requested_quota=quotas[wid], selected_sites=len(sites),
                          sampling_fraction=len(sites)/seen[wid] if seen[wid] else None,
                          scan_supported_eligible_sites=w["n_eligible"], scan_raw_eligible_sites=w.get("raw_n_eligible", "NA"),
                          status="no_source_sites" if not seen[wid] else "sampled" if len(sites)<seen[wid] else "all_provided_source_sites"))
    targets.sort(key=lambda r: (characterize_chrom_order(r["chrom"]), r["pos"]))
    requests = [dict(t, request_scope="exact_nucleotide_request_not_selection_evidence",
                     site_pairing_status="not_assigned_window_match_is_not_an_ancient_site_pair") for t in targets]
    characterize_write_tsv(out/"targets.tsv", targets, CHARACTERIZE_TARGET_FIELDS)
    characterize_write_tsv(out/"selection_catalog.tsv", targets, CHARACTERIZE_TARGET_FIELDS)
    characterize_write_tsv(out/"ancient_target_requests.tsv", requests, CHARACTERIZE_REQUEST_FIELDS)
    characterize_write_tsv(out/"window_membership.tsv", targets, ["site_id", "window_id", "window_match_id", "group", "family", "chrom", "pos", "ref", "source_sampling"])
    return targets, audit, dict(status="complete", counts=dict(counts), selected_sites=len(targets), source_sampling=sampling,
                               completeness="sampling-frame declaration; not proof of whole-window or all-callable-site coverage",
                               equal_preallocated_quota_per_case_and_control_window=base,
                               unused_global_capacity_before_missing_source_sites=maximum-base*len(rows),
                               sampling_algorithm="smallest SHA256(seed|canonical_chrom:pos:ref) within each window; bounded fixed allocation across windows")


def characterize_grouped_rows(path, required):
    previous, rows, seen = None, [], set()
    for row in characterize_read_tsv(path, {*required, "site_id"}):
        sid = row["site_id"]
        if sid.upper() in CHARACTERIZE_NA: raise ValueError("Atlas site_id is missing")
        if sid != previous:
            if previous is not None: yield previous, rows
            if sid in seen: raise ValueError("Atlas derivative rows must be contiguous by site_id")
            seen.add(sid); previous, rows = sid, []
        rows.append(row)
        if len(rows) > 1000: raise ValueError("An Atlas site has >1000 rows; supply bounded scalar derivative summaries")
    if previous is not None: yield previous, rows


def characterize_score_identity(row):
    characterize_build38(row["build"])
    chrom, pos, ref, allele = characterize_canonical_chrom(row["chrom"]), characterize_numeric(row["pos"], "Atlas POS", integer=True), row["ref"].upper(), row["allele"].upper()
    if pos < 1 or ref not in CHARACTERIZE_BASES or allele not in CHARACTERIZE_BASES: raise ValueError("Invalid Atlas nucleotide identity")
    return chrom, pos, ref, allele


def characterize_atlas_profiles(variants, metrics, target_path, internal_targets, selected, threshold, maximum=50000):
    if not variants:
        if metrics or target_path: raise ValueError("--atlas-metrics/--atlas-targets require --atlas-variants")
        return [], dict(status="not_requested")
    index = characterize_WindowIndex(selected); selected_map = {w["window_id"]: w for w in selected}
    target_rows = internal_targets
    if target_path:
        target_rows = []
        for row in characterize_read_tsv(target_path, {"site_id", "chrom", "pos", "ref"}):
            if len(target_rows) >= maximum: raise ValueError("--max-atlas-sites exceeded by target counts")
            target_rows.append(characterize_parse_site(row))
    targets, positions = {}, {}
    for t in target_rows:
        key = t["chrom"], t["pos"]
        if t["site_id"] in targets or key in positions: raise ValueError("Atlas target-count table requires unique site IDs and positions")
        if len(targets) >= maximum: raise ValueError("--max-atlas-sites exceeded by target counts")
        targets[t["site_id"]] = t; positions[key] = t["site_id"]
    expected = defaultdict(set)
    for t in target_rows:
        for w in index.overlaps(t["chrom"], t["pos"]-1, t["pos"]): expected[w["window_id"]].add(t["site_id"])
    has_target_denominator = bool(target_path or internal_targets)
    available_sites, source_sites, source_versions, audit = defaultdict(set), {}, {}, Counter()
    values, complete, any_scored = defaultdict(list), defaultdict(set), defaultdict(set)
    metric_specs = {
        "AVI_raw": ("gnomad_frequency_trained_variant_impact", "AVI score; not fitness, pathogenicity probability or tolerance"),
        "AVI_PHRED": ("gnomad_frequency_trained_avi_phred_rank", "Genome-wide AVI rank; PHRED20 means top1%, not probability"),
    }
    avi_sources = set()
    def register_source(resource, version):
        if resource in source_versions and source_versions[resource] != version:
            raise ValueError(f"Mixed source versions for {resource}; no cross-version window profile")
        source_versions[resource] = version
    def add_scores(sid, site, wid, metric, scores):
        if scores: any_scored[wid, metric].add(sid)
        if set(scores) != CHARACTERIZE_BASES-{site[2]}: return
        complete[wid, metric].add(sid)
        for allele, value in scores.items():
            values[wid, metric, "all_three_possible"].append(value)
            counts = targets.get(sid, {}).get("_counts")
            if counts is None:
                values[wid, metric, "source_counts_unavailable"].append(value)
            elif allele not in counts:
                values[wid, metric, "not_in_source_rows"].append(value)
            else:
                ac = counts[allele][0]
                if ac > 0: values[wid, metric, "reported_positive_AC"].append(value)
                label = "reported_zero" if ac == 0 else f"AC{threshold}_supported" if ac >= threshold else f"reported_below_AC{threshold}"
                values[wid, metric, label].append(value)
    for sid, rows in characterize_grouped_rows(variants, {"chrom", "pos", "build", "ref", "allele", "avi_raw", "avi_phred", "status", "resource_id", "source_version"}):
        if len(source_sites) >= maximum: raise ValueError("--max-atlas-sites exceeded; do not read the 8.8B-row source as a derivative")
        first = characterize_score_identity(rows[0])[:3]
        if has_target_denominator and sid not in targets: raise ValueError("Atlas contains a site absent from the supplied target-count table")
        if sid in targets and (targets[sid]["chrom"], targets[sid]["pos"], targets[sid]["ref"]) != first:
            raise ValueError("Atlas REF/coordinate differs from target counts")
        source_sites[sid] = first
        matches = index.overlaps(first[0], first[1]-1, first[1])
        seen, raw_scores, phred_scores = set(), {}, {}
        for row in rows:
            chrom, pos, ref, allele = characterize_score_identity(row)
            if (chrom, pos, ref) != first or allele in seen: raise ValueError("Contradictory/duplicate Atlas allele rows")
            seen.add(allele); register_source(row["resource_id"], row["source_version"])
            avi_sources.add((row["resource_id"], row["source_version"]))
            if len(avi_sources) > 1: raise ValueError("Atlas variant profiles require one AVI resource/version across sites")
            raw = characterize_numeric(row["avi_raw"], "AVI raw", optional=True); phred = characterize_numeric(row["avi_phred"], "AVI PHRED", optional=True)
            if phred is not None and phred < 0: raise ValueError("AVI PHRED must be nonnegative")
            if allele == ref:
                if raw is not None or phred is not None: raise ValueError("REF has no REF-to-REF AVI substitution score")
                continue
            if row["status"] not in {"available", "available_partial"}:
                if raw is not None or phred is not None: raise ValueError("Unavailable Atlas row carries a score")
                continue
            if raw is not None: raw_scores[allele] = raw
            if phred is not None: phred_scores[allele] = phred
        audit["variant_sites"] += 1
        for w in matches:
            wid = w["window_id"]; available_sites[wid].add(sid)
            if not has_target_denominator: expected[wid].add(sid)
            add_scores(sid, first, wid, "AVI_raw", raw_scores)
            add_scores(sid, first, wid, "AVI_PHRED", phred_scores)
    coordinate_ids = {}
    for sid, loc in source_sites.items():
        if loc[:2] in coordinate_ids: raise ValueError("Atlas site IDs duplicate a genomic position")
        coordinate_ids[loc[:2]] = sid
    if metrics:
        for sid, rows in characterize_grouped_rows(metrics, {"chrom", "pos", "build", "ref", "allele", "resource_id", "source_version", "metric", "value", "status"}):
            if sid not in source_sites: raise ValueError("Atlas metrics contain a site absent from Atlas variants")
            first = source_sites[sid]; scores = defaultdict(dict); keys = set()
            for row in rows:
                chrom, pos, ref, allele = characterize_score_identity(row)
                if (chrom, pos, ref) != first: raise ValueError("Atlas metric identity/REF conflicts with variants")
                resource, metric = row["resource_id"], row["metric"]
                register_source(resource, row["source_version"])
                key = resource, metric, allele
                if key in keys: raise ValueError("Duplicate Atlas metric/base row")
                keys.add(key)
                if metric in CHARACTERIZE_ATTRIBUTIONS:
                    name = resource+":"+metric
                    spec = ("signed_shap_contribution_to_avi", "Signed AVI attribution; not raw molecular effect, direction, function count or tolerance")
                elif metric == "alphagenome_splicing":
                    name = resource+":"+metric
                    spec = ("reference_relative_splicing_prediction", "Sparse predicted splicing disruption; missing is not zero")
                else:
                    audit["ignored_unrecognized_or_duplicate_metric_rows"] += 1
                    continue
                declared = row.get("metric_semantics", "")
                if declared and declared != spec[0]: raise ValueError("Known Atlas metric is labelled with conflicting semantics")
                metric_specs[name] = spec
                value = characterize_numeric(row["value"], "Atlas metric", optional=True)
                if allele == ref:
                    if value is not None: raise ValueError("Structural REF metric must remain unavailable")
                    continue
                if row["status"] != "available":
                    if value is not None: raise ValueError("Unavailable Atlas metric carries a value")
                    continue
                if value is not None:
                    if metric == "alphagenome_splicing" and value < 0: raise ValueError("Splicing summary cannot be negative")
                    scores[name][allele] = value
            for w in index.overlaps(first[0], first[1]-1, first[1]):
                for name, v in scores.items(): add_scores(sid, first, w["window_id"], name, v)
    groups = ["all_three_possible", "reported_positive_AC", f"AC{threshold}_supported", f"reported_below_AC{threshold}", "reported_zero", "not_in_source_rows", "source_counts_unavailable"]
    result = []
    for wid, w in selected_map.items():
        for metric, (semantics, limit) in sorted(metric_specs.items()):
            for group in groups:
                data = values[wid, metric, group]
                result.append(dict(window_id=wid, window_match_id=w["window_match_id"], group=w["group"], metric=metric,
                                   allele_scope=group, target_sites=len(expected[wid]), sites_in_atlas_input=len(available_sites[wid]),
                                   sites_with_any_metric=len(any_scored[wid, metric]), sites_with_complete_three_metric=len(complete[wid, metric]),
                                   theoretical_substitutions=3*len(expected[wid]), analyzed_substitutions=len(data),
                                   mean=statistics.fmean(data) if data else None, median=statistics.median(data) if data else None,
                                   mean_absolute=statistics.fmean(abs(v) for v in data) if data else None,
                                   fraction_avi_phred_ge20=sum(v>=20 for v in data)/len(data) if data and metric=="AVI_PHRED" else None,
                                   metric_semantics=semantics, interpretation_limit=limit,
                                   denominator_scope="provided_target_table" if has_target_denominator else "represented_Atlas_sites_only",
                                   status="descriptive_complete_three_ALT_profile" if data else "unavailable_no_complete_three_ALT_values"))
    return result, dict(status="complete", counts=dict(audit), source_versions=source_versions,
                        target_denominator="provided_target_table" if has_target_denominator else "represented_Atlas_sites_only",
                        data_scope="bounded_target_derivatives_not_all_genomic_bases", target_count_sites=len(targets))


CHARACTERIZE_SELECTED_FIELDS = ["window_id", "chrom", "start0", "end0", "build", "window_bp", "span_bp", "threshold", "n_eligible",
                   "observed_multi", "expected_multi", "oe", "p_value", "q_value", "nominated", "group", "window_match_id", "match_role",
                   "mask_overlap_bp", "mask_status"]
CHARACTERIZE_MATCH_FIELDS = ["window_match_id", "case_window_id", "control_window_id", "control_rank", "matching_distance", "n_controls", "status", "representative_window_id"]
CHARACTERIZE_BALANCE_FIELDS = ["window_match_id", "case_window_id", "control_window_id", "covariate", "transformation", "case_value", "control_value", "pooled_tested_sd", "standardized_difference"]
CHARACTERIZE_ANNOTATION_FIELDS = ["window_id", "chrom", "start0", "end0", "group", "window_match_id", "match_role", "representative_window_id",
                     "annotation_source", "category", "union_overlap_bp", "window_bp", "overlap_fraction", "overlapping_records", "status"]
CHARACTERIZE_CONTRAST_FIELDS = ["window_match_id", "case_window_id", "annotation_source", "category", "n_controls", "case_fraction", "mean_control_fraction", "difference", "ratio", "status"]
CHARACTERIZE_SAMPLING_FIELDS = ["window_id", "group", "window_match_id", "source_sampling", "source_rows_in_window", "requested_quota", "selected_sites", "sampling_fraction", "scan_supported_eligible_sites", "scan_raw_eligible_sites", "status"]
CHARACTERIZE_ATLAS_FIELDS = ["window_id", "window_match_id", "group", "metric", "allele_scope", "target_sites", "sites_in_atlas_input", "sites_with_any_metric",
                "sites_with_complete_three_metric", "theoretical_substitutions", "analyzed_substitutions", "mean", "median", "mean_absolute",
                "fraction_avi_phred_ge20", "metric_semantics", "interpretation_limit", "denominator_scope", "status"]


@transactional("manifest.json", input_files=("windows", "covariates", "atlas_variants", "atlas_metrics", "atlas_targets"), input_prefixes=(), resolve_inputs=characterize_protected_sources)
def characterize_run(args):
    characterize_build38(args.build)
    if args.threshold not in (5,10) or args.controls_per_case < 1 or args.controls_per_case > 20 or not 0 < args.max_match_distance <= 10 or args.control_gap_bp < 0:
        raise ValueError("Require AC5/AC10, 1..20 controls, positive match distance<=10 and nonnegative gap")
    if min(args.max_windows, args.max_cases, args.sites_per_window, args.max_targets, args.max_atlas_sites) < 1:
        raise ValueError("Analysis and sampling caps must be positive")
    columns = characterize_matching_columns(args.match_columns)
    specs = characterize_annotation_specs(args.annotation_bed)
    windows, scope = characterize_load_windows(args.windows, args.threshold, args.window_bp, args.max_windows, args.covariates)
    selected, matches, balance, leads, matching = characterize_match_windows(windows, columns, args.controls_per_case, args.max_match_distance,
                                                               args.control_gap_bp, args.seed, args.max_cases)
    annotations, annotation_audit = characterize_annotate_windows(windows, selected, leads, specs)
    out = Path(args.out)
    characterize_write_tsv(out/"selected_windows.tsv", selected, CHARACTERIZE_SELECTED_FIELDS)
    with (out/"selected_windows.bed").open("w") as handle:
        for row in sorted(selected, key=lambda r: (characterize_chrom_order(r["chrom"]), r["start0"], r["window_id"])):
            handle.write(f"{row['chrom']}\t{row['start0']}\t{row['end0']}\t{row['window_id']}\t0\t.\n")
    characterize_write_tsv(out/"window_matches.tsv", matches, CHARACTERIZE_MATCH_FIELDS)
    characterize_write_tsv(out/"window_match_balance.tsv", balance, CHARACTERIZE_BALANCE_FIELDS)
    characterize_write_tsv(out/"window_annotations.tsv", annotations, CHARACTERIZE_ANNOTATION_FIELDS)
    characterize_write_tsv(out/"annotation_contrasts.tsv", characterize_annotation_contrasts(annotations, matches), CHARACTERIZE_CONTRAST_FIELDS)
    targets, sampling_audit, sampling = characterize_export_targets(args.site_table, args.site_sampling, selected, out,
                                                       args.sites_per_window, args.max_targets, args.seed)
    if not args.site_table:
        for filename, fields in [("targets.tsv", CHARACTERIZE_TARGET_FIELDS), ("selection_catalog.tsv", CHARACTERIZE_TARGET_FIELDS),
                                 ("ancient_target_requests.tsv", CHARACTERIZE_REQUEST_FIELDS),
                                 ("window_membership.tsv", ["site_id", "window_id", "window_match_id", "group", "family", "chrom", "pos", "ref", "source_sampling"])]:
            characterize_write_tsv(out/filename, [], fields)
    characterize_write_tsv(out/"site_sampling_audit.tsv", sampling_audit, CHARACTERIZE_SAMPLING_FIELDS)
    atlas, atlas_audit = characterize_atlas_profiles(args.atlas_variants, args.atlas_metrics, args.atlas_targets, targets, selected, args.threshold, args.max_atlas_sites)
    characterize_write_tsv(out/"window_atlas_summary.tsv", atlas, CHARACTERIZE_ATLAS_FIELDS)
    return dict(status="complete", schema=CHARACTERIZE_VERSION, build="GRCh38", coordinate_system="windows_BED0_half_open; sites_1_based",
                nominations_status="present" if matching["nominated_windows"] else "no_candidates",
                scope=scope, matching=matching, annotations=annotation_audit, target_sampling=sampling, atlas=atlas_audit,
                source=identity(args.windows, full=True),
                parameters={k: getattr(args,k) for k in ("threshold", "window_bp", "controls_per_case", "match_columns", "max_match_distance", "control_gap_bp", "seed", "sites_per_window", "max_targets")},
                limitations=CHARACTERIZE_LIMITS)


def characterize_parser():
    p = argparse.ArgumentParser(description=_description_hotspot_characterize)
    p.add_argument("--windows", required=True, help="Complete genome scan all_windows.tsv.gz, including non-nominated tested windows")
    p.add_argument("--build", required=True, choices=("GRCh38", "38", "hg38"))
    p.add_argument("--out", required=True)
    p.add_argument("--threshold", type=int, choices=(5,10), default=5)
    p.add_argument("--window-bp", type=int, help="Select one scale; defaults to the scan's unique primary size")
    p.add_argument("--annotation-bed", action="append", default=[], metavar="NAME=PATH", help="Repeatable GRCh38 sorted BED4 category source; annotations never select cases")
    p.add_argument("--covariates", help="Optional window_id keyed TSV of approved pre-discovery matching covariates")
    p.add_argument("--match-columns", default="n_eligible")
    p.add_argument("--controls-per-case", type=int, default=2)
    p.add_argument("--max-match-distance", type=float, default=1., help="RMS standardized covariate distance; n_eligible is log1p transformed")
    p.add_argument("--control-gap-bp", type=int, default=0, help="Minimum extra gap from every nomination and already-used control")
    p.add_argument("--max-windows", type=int, default=500000)
    p.add_argument("--max-cases", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261106)
    p.add_argument("--site-table", action="append", default=[], help="Repeatable sorted site count caches; full ALT AC/AN arrays required")
    p.add_argument("--site-sampling", choices=("selected_multialt_candidates", "all_eligible_source_sites", "prespecified_sites"))
    p.add_argument("--sites-per-window", type=int, default=25)
    p.add_argument("--max-targets", type=int, default=5000)
    p.add_argument("--atlas-variants", help="Bounded atlas_variants.tsv.gz derivative, never the original 8.8B row source")
    p.add_argument("--atlas-metrics", help="Optional bounded atlas_metrics.tsv.gz derivative")
    p.add_argument("--atlas-targets", help="Exact target table with actual reported ALT AC/AN arrays for observation categories")
    p.add_argument("--max-atlas-sites", type=int, default=50000)
    return p


def characterize_main(argv=None):
    try: result = characterize_run(characterize_parser().parse_args(argv))
    except (ValueError, OSError, RuntimeError) as exc: raise SystemExit("hotspot-characterize: " + str(exc))
    print(json.dumps(dict(status=result["status"], nominations_status=result["nominations_status"],
                          representative_cases=result["matching"]["representative_cases"], matched_controls=result["matching"]["matched_controls"])))
    return 0


# Window Selection
# ======================================================================
_description_window_selection = 'Descriptive published ancient-selection evidence in a frozen window panel.\n\nExact biallelic contrasts are mapped with the existing reciprocal FASTA-checked\nmatcher. REF and ALT views of the same contrast count once in each window.\nThis command does not infer selection on an unmeasured third nucleotide.\n'


import chronology as hs

SELECTION_VERSION = "window-published-selection-1.0"
SELECTION_LIMITATIONS = [
    "Descriptive overlap with published Akbari et al. 2026 selection contrasts; no new selection model or window-level selection P value is fitted.",
    "The frozen target panel is sampled from nominated/control windows; availability and ascertainment limit comparisons between windows.",
    "One published REF/ALT contrast is counted once per window, even when both nucleotide views are present.",
    "A published contrast is never transferred to another ALT at the same position. Missing ancient alleles are unavailable, not extinct or AF=0.",
    "Published FDR values remain the source study's FDR. A positive S direction is interpreted only under an explicit caller-verified ALT effect contract.",
    "Higher multi-ALT occupancy does not by itself establish balancing selection, directional positive selection, stability or mutational birth time.",
]
SELECTION_FIELDS = ["window_id", "group", "match_id", "catalog_sites", "catalog_alt_labels",
          "sites_with_published_contrast", "published_exact_contrasts",
          "valid_accepted_filter_contrasts", "published_FDR_le_0_05_contrasts",
          "published_FDR_le_0_05_fraction_of_valid_contrasts",
          "published_ALT_positive_S_FDR_le_0_05_contrasts",
          "exact_secondary_label_contrasts", "catalog_site_matching_fraction",
          "effect_allele_contract", "analysis_status"]
SELECTION_DETAIL = ["window_id", "group", "match_id", "site_id", "published_contrast",
          "statistics_status", "published_effect_contract", "source_doi", *hs.SELECTION_PUBLISHED]


def selection_rows(path, required=()):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: missing required columns {sorted(set(required) - set(reader.fieldnames or []))}")
        yield from reader


def selection_write(path, fields, values):
    with open(path, "w", newline="") as handle:
        output = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        output.writeheader()
        for row in values:
            output.writerow({k: "NA" if row.get(k) is None else row.get(k, "NA") for k in fields})


def selection_finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def selection_summarize(catalog, membership, matched, out, effect_allele="raw_only"):
    """Bounded selected-panel aggregation; reject inconsistent duplicate views."""
    sites = {}
    for row in hs.selection_catalog_rows(catalog):
        key = row["site_id"]
        if key in sites:
            raise ValueError(f"Duplicate catalog site_id: {key}")
        sites[key] = row
    windows, site_windows = {}, defaultdict(list)
    seen = set()
    for row in selection_rows(membership, ("site_id", "window_id", "group")):
        row["match_id"] = row.get("match_id", row.get("window_match_id", ""))
        key, wid = row["site_id"], row["window_id"]
        if key not in sites:
            raise ValueError(f"Window membership has no catalog site: {key}")
        if (key, wid) in seen:
            raise ValueError(f"Duplicate site/window membership: {key}, {wid}")
        seen.add((key, wid))
        if row["group"] not in {"case", "control"} or not wid or not row["match_id"]:
            raise ValueError("Window membership needs window_id, case/control group and match_id")
        group, match = row["group"], row["match_id"]
        if (group, match) != (sites[key]["group"], sites[key]["match_id"]):
            raise ValueError(f"Window membership disagrees with catalog group/match_id: {key}")
        state = windows.setdefault(wid, dict(group=group, match_id=match, sites=set(), alts=0, contrasts={}))
        if (state["group"], state["match_id"]) != (group, match):
            raise ValueError(f"Inconsistent group or match_id for window: {wid}")
        state["sites"].add(key)
        state["alts"] += len(sites[key]["alts"].split(","))
        site_windows[key].append(wid)
    if set(sites) != set(site_windows):
        raise ValueError("Every selection catalog site must have a window membership")
    n_rows = 0
    for row in selection_rows(matched, ("site_id", "published_contrast", "statistics_status", "published_effect_contract", *hs.SELECTION_PUBLISHED)):
        n_rows += 1
        site, contrast = row["site_id"], row["published_contrast"]
        if site not in sites:
            raise ValueError(f"Published match has no catalog site: {site}")
        if not contrast:
            raise ValueError("Missing published contrast identity")
        expected_contract = "caller_declared_ALT" if effect_allele == "ALT" else "raw_only_unverified"
        if row["published_effect_contract"] != expected_contract:
            raise ValueError("Published effect contract disagrees with the aggregation request")
        signature = tuple(row.get(k) for k in ("site_id", "statistics_status", "published_effect_contract", *hs.SELECTION_PUBLISHED))
        for wid in site_windows[site]:
            contrasts = windows[wid]["contrasts"]
            previous = contrasts.get(contrast)
            if previous is None:
                contrasts[contrast] = dict(signature=signature, row=row, secondary=False, views=set())
                previous = contrasts[contrast]
            elif previous["signature"] != signature:
                raise ValueError(f"Conflicting duplicate published contrast: {contrast}")
            view = (row.get("counted_base38"), row.get("counted_base_relation"))
            if view in previous["views"]:
                raise ValueError(f"Duplicate nucleotide view for published contrast: {contrast}")
            previous["views"].add(view)
            previous["secondary"] |= row.get("role") == "secondary"
    summaries, details = [], []
    for wid, state in sorted(windows.items()):
        valid, significant, positive, secondary = 0, 0, 0, 0
        matched_sites = set()
        for contrast, record in sorted(state["contrasts"].items()):
            row = record["row"]
            good = row["statistics_status"] == "valid" and hs.selection_statistics_valid(row) and row["FILTER"] in {"PASS", "."}
            sig = good and selection_finite(row["FDR"]) <= 0.05
            valid += good
            significant += sig
            positive += bool(sig and selection_finite(row["S"]) > 0)
            secondary += record["secondary"]
            matched_sites.add(row["site_id"])
            details.append(dict(row, window_id=wid, group=state["group"], match_id=state["match_id"]))
        n = len(state["sites"])
        summaries.append(dict(window_id=wid, group=state["group"], match_id=state["match_id"],
            catalog_sites=n, catalog_alt_labels=state["alts"], sites_with_published_contrast=len(matched_sites),
            published_exact_contrasts=len(state["contrasts"]), valid_accepted_filter_contrasts=valid,
            published_FDR_le_0_05_contrasts=significant,
            published_FDR_le_0_05_fraction_of_valid_contrasts=significant / valid if valid else None,
            published_ALT_positive_S_FDR_le_0_05_contrasts=positive if effect_allele == "ALT" else None,
            exact_secondary_label_contrasts=secondary,
            catalog_site_matching_fraction=len(matched_sites) / n if n else None,
            effect_allele_contract=effect_allele,
            analysis_status="descriptive_published_evidence" if valid else "no_valid_published_contrast"))
    out = Path(out)
    selection_write(out / "window_selection_summary.tsv", SELECTION_FIELDS, summaries)
    selection_write(out / "window_selection_contrasts.tsv", SELECTION_DETAIL, details)
    return dict(catalog_sites=len(sites), windows=len(windows), nucleotide_view_rows=n_rows,
                unique_window_contrasts=sum(len(s["contrasts"]) for s in windows.values()),
                window_summary_rows=len(summaries))


@transactional("manifest.json")
def selection_run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    matcher_args = hs.selection_parser().parse_args([
        "--catalog", args.catalog, "--selection-stats", args.selection_stats,
        "--fasta38", args.fasta38, "--fasta37", args.fasta37,
        "--chain38to37", args.chain38to37, "--chain37to38", args.chain37to38,
        "--effect-allele", args.effect_allele, "--out", str(out / "exact_matches")])
    match_result = hs.selection_run(matcher_args)
    qc = selection_summarize(args.catalog, args.window_membership,
        out / "exact_matches/published_selection_alleles.tsv.gz", out, args.effect_allele)
    return dict(status="complete", version=SELECTION_VERSION, qc=qc, source_doi=hs.SELECTION_SOURCE,
                dataset_doi=hs.SELECTION_DATASET, matcher_status=match_result["status"],
                scope="descriptive exact published contrasts in a frozen bounded window panel",
                effect_allele_contract=args.effect_allele, limitations=SELECTION_LIMITATIONS)


def selection_parser():
    p = argparse.ArgumentParser(description=_description_window_selection)
    for name in ("catalog", "window-membership", "selection-stats", "fasta38", "fasta37", "chain38to37", "chain37to38", "out"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--effect-allele", choices=("raw_only", "ALT"), default="raw_only",
                   help="ALT requires independent verification of the published S/X orientation")
    return p


def selection_main(argv=None):
    result = selection_run(selection_parser().parse_args(argv))
    print(json.dumps({"status": result["status"], "qc": result["qc"]}, sort_keys=True))
    return 0


def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'characterize': characterize_main, 'selection': selection_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == "__main__":
    raise SystemExit(cli())


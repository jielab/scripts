#!/usr/bin/env python3
"""Actual GRCh38 PCA reruns, invariant geometry comparisons and driver analysis.

Time is metadata, never an input column of the genotype PCA. All scenarios
project the same eligible people. A sensitivity changes training people,
markers or PLINK parameters before fitting; it never edits plotted coordinates.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import contextlib
import csv
import fcntl
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np
from scipy.linalg import subspace_angles
from scipy.spatial import cKDTree
from scipy.stats import chi2, spearmanr

from allele_measurement import (stats_technology, stats_number,
    stats_bootstrap_mean, stats_fixed_composition, stats_adjust_fdr)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "0f"))
from workflow_paths import configured_path

VERSION = "pca-study-1"
REGIONS = {"C": "Central", "E": "East", "N": "North", "SE": "Southeast", "SW": "Southwest"}
AUTOSOMES = {str(x) for x in range(1, 23)}
BASE = dict(maf=.05, geno=.10, ld_r2=.10, ld_kb=200, basis="joint_balanced",
            cap=200, seed_offset=0, transversions=False, long_ld=True,
            min_coverage=0., training_technology="all", omit_region="",
            omit_chrom="", top_windows=0)
SCENARIOS = {
    "baseline": {},
    "maf01": {"maf": .01}, "maf10": {"maf": .10},
    "missing02": {"geno": .02},
    "ld05": {"ld_r2": .05}, "ld20": {"ld_r2": .20},
    "seed2": {"seed_offset": 1}, "seed3": {"seed_offset": 2},
    "cap100": {"cap": 100},
    "ancient": {"basis": "ancient_balanced"},
    "modern": {"basis": "modern_projection"},
    "transversions": {"transversions": True},
    "coverage1": {"min_coverage": 1.},
    "shotgun": {"training_technology": "shotgun_label"},
    "long_ld_included": {"long_ld": False},
    "no_top_windows": {"top_windows": 5},
    **{"omit_" + r: {"omit_region": r} for r in REGIONS},
    **{"omit_chr" + str(c): {"omit_chrom": str(c)} for c in range(1, 23)},
}
CORE = ["baseline", "maf01", "maf10", "missing02", "ld05", "ld20", "seed2",
        "cap100", "ancient", "modern", "transversions", "coverage1", "shotgun",
        "long_ld_included", "no_top_windows"]
LIMITATIONS = [
    "PC1/PC2 are genetic coordinates; age BP is metadata, not PC3 or a predicted date.",
    "Rotations, signs and one global scale are aligned using genetic PCs only; age is excluded from every geometry metric.",
    "Stability is conditional on the supplied imputed panel. It is not independent ancient-read validation.",
    "Regional changes can reflect migration, mixture, sampling and technical differences; no selection coefficient is inferred.",
    "Squared allele-loading shares locate influential markers; they are not causal contributions or total variance explained.",
    "Technical and age association tests are descriptive conditional associations; collinear blocks may be unidentifiable.",
    "Bootstrap CIs resample people conditional on fitted coordinates, without refitting PCA or inventing dating uncertainty.",
    "Modern samples remain separate. Sparse regions/periods are missing, never zero or interpolated.",
]


@contextlib.contextmanager
def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as h: yield h
    elif path.suffix == ".zst":
        import io, zstandard
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as stream:
            with io.TextIOWrapper(stream) as h: yield h
    else:
        with path.open() as h: yield h


def rows(path):
    with open_text(path) as h:
        source = (x for x in h if x.strip() and not x.startswith("##"))
        try: header = next(source).strip().lstrip("#").split()
        except StopIteration: raise ValueError("Empty table: " + str(path))
        # Canonical metadata has spaces in field values; all generated tables use tabs.
        for line in source:
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) != len(header): cells = line.split()
            if len(cells) != len(header): raise ValueError("Malformed row in " + str(path))
            yield dict(zip(header, cells))


def table(path, records, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    records = list(records)
    fields = fields or list(dict.fromkeys(k for r in records for k in r))
    fields = fields or ["status"]
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(partial, "wt") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            clean = {}
            for k in fields:
                v = r.get(k)
                if isinstance(v, np.generic): v = v.item()
                if isinstance(v, float) and not math.isfinite(v): v = None
                clean[k] = "." if v is None else json.dumps(v, separators=(",", ":")) if isinstance(v, (list, dict)) else v
            w.writerow(clean)
    partial.replace(path)


def save_json(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(path.name + ".partial")
    pending.write_text(json.dumps(obj, indent=2, allow_nan=False, default=str) + "\n")
    pending.replace(path)


def state(path):
    path = Path(path).absolute()
    s = path.stat()
    if not s.st_size: raise ValueError("Empty input: " + str(path))
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns)


def complete(paths):
    return all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in paths)


def read_scores(path):
    rr = list(rows(path))
    if not rr: raise ValueError("No PCA coordinates: " + str(path))
    names = sorted((k for k in rr[0] if re.fullmatch(r"PC\d+_AVG", k)), key=lambda x: int(x[2:-4]))
    if len(names) < 2: raise ValueError("Projected PC*_AVG scores required")
    ids = [r["IID"] for r in rr]
    if len(set(ids)) != len(ids): raise ValueError("Duplicate PCA sample IDs")
    values = np.array([[float(r[k]) for k in names] for r in rr])
    if not np.isfinite(values).all(): raise ValueError("Nonfinite PCA coordinates")
    return ids, values


def load_metadata(path):
    records = list(rows(path))
    if not records or not {"sample_id", "cohort", "age_bp", "include", "unrelated", "region"} <= records[0].keys():
        raise ValueError("Use canonical sample_metadata.tsv from 1.ancient_pca.sh")
    result = {}
    for r in records:
        sid = r["sample_id"]
        if sid in result: raise ValueError("Duplicate metadata ID: " + sid)
        if r["cohort"] not in {"ancient","modern"} or any(r[k] not in {"0","1"} for k in ("include","unrelated")):
            raise ValueError("Invalid canonical cohort/include/unrelated metadata: " + sid)
        if r["include"]=="1" and r["cohort"]=="ancient" and r["region"] not in REGIONS:
            raise ValueError("Included ancient sample needs an official broad region: " + sid)
        r = dict(r); r["technology"] = stats_technology(r)
        r["age"] = stats_number(r["age_bp"])
        r["coverage"] = stats_number(r.get("target_coverage"))
        if r["include"] == "1" and (r["age"] is None or r["age"] < 0):
            raise ValueError("Included sample has invalid age: " + sid)
        if r["coverage"] is not None and r["coverage"]<0:
            raise ValueError("Negative target coverage: " + sid)
        result[sid] = r
    return result


def scenario_settings(args):
    definitions = {name: dict(BASE, **options) for name, options in SCENARIOS.items()}
    if args.scenario_file:
        supplied = json.loads(Path(args.scenario_file).read_text())
        if not isinstance(supplied, dict): raise ValueError("Scenario JSON must map names to settings")
        for name, values in supplied.items():
            if not re.fullmatch(r"[A-Za-z0-9_]+", name) or not isinstance(values, dict):
                raise ValueError("Invalid scenario name/settings")
            if set(values) - BASE.keys(): raise ValueError("Unknown scenario settings: " + ",".join(set(values) - BASE.keys()))
            definitions[name] = dict(BASE, **values)
    wanted = CORE if args.scenarios == "core" else list(definitions) if args.scenarios == "full" else args.scenarios.split(",")
    wanted = list(dict.fromkeys(["baseline", *wanted]))
    unknown = set(wanted) - definitions.keys()
    if unknown: raise ValueError("Unknown scenarios: " + ",".join(sorted(unknown)))
    result = {}
    for name in wanted:
        s = definitions[name]
        if not 0 < s["maf"] < .5 or not 0 <= s["geno"] < 1 or not 0 < s["ld_r2"] <= 1:
            raise ValueError("Invalid MAF, missingness or LD threshold: " + name)
        if s["basis"] not in {"joint_balanced", "ancient_balanced", "modern_projection"}:
            raise ValueError("Unknown PCA basis: " + name)
        if any(not isinstance(s[k], int) or s[k] < 1 for k in ("cap", "ld_kb")):
            raise ValueError("cap and ld_kb must be positive integers")
        if any(not isinstance(s[k], int) or s[k] < 0 for k in ("seed_offset", "top_windows")) or s["min_coverage"] < 0:
            raise ValueError("Seed offset, top windows and minimum coverage must be nonnegative")
        if any(not isinstance(s[k], bool) for k in ("transversions", "long_ld")):
            raise ValueError("transversions and long_ld must be JSON booleans")
        if not 0 <= args.seed+s["seed_offset"] <= 4294967295:
            raise ValueError("PLINK seed must fit an unsigned 32-bit integer")
        if s["training_technology"] not in {"all", "shotgun_label", "capture_label"}:
            raise ValueError("Unknown training_technology")
        if s["omit_region"] not in {"", *REGIONS} or s["omit_chrom"] not in {"", *map(str, range(1, 23))}:
            raise ValueError("Invalid omitted region/chromosome")
        s = dict(s, seed=args.seed + s["seed_offset"], components=args.components, time_bin_years=2000)
        result[name] = s
    if result["baseline"].get("top_windows"): raise ValueError("Baseline cannot depend on its own loadings")
    return result


def source_inputs(args, out, runtime):
    prefix = Path(args.pgen_prefix or Path(args.ancient_dir) / "38/ancient").absolute()
    from ancient_panel import require_prepared_panel
    require_prepared_panel(prefix)
    pvar = Path(str(prefix) + ".pvar")
    for suffix in (".pgen", ".pvar", ".psam"):
        state(str(prefix) + suffix)
    from ancient_panel import panel_build, hq_in_panel_coordinates
    if panel_build(pvar) != "GRCh38":
        raise ValueError("Study requires the current ancient/38/ancient GRCh38 panel")
    meta_path = Path(args.metadata) if args.metadata else Path(args.ancient_dir) / "38/pca3d/tables/sample_metadata.tsv"
    metadata_sources = {}
    if (not args.metadata and args.supplement_xlsx) or not meta_path.is_file():
        if args.metadata: raise ValueError("Explicit metadata does not exist: " + str(meta_path))
        supplement = Path(args.supplement_xlsx or Path(args.ancient_dir) / "files/Nature.xlsx")
        metadata_sources = dict(supplement_xlsx=state(supplement), psam=state(str(prefix)+".psam"))
        metadata_dir = runtime / ("plan_metadata" if args.command=="study-plan" else "metadata")
        meta_path = metadata_dir / "sample_metadata.tsv"
        receipt = metadata_dir / "input_state.json"
        prior = json.loads(receipt.read_text()) if receipt.is_file() else {}
        if prior != metadata_sources or not complete([meta_path, metadata_dir/"metadata_qc.json"]):
            from evolve import prepare_metadata
            created = prepare_metadata(str(prefix) + ".psam", supplement, metadata_dir)
            meta_path = Path(created["metadata_path"])
            save_json(receipt, metadata_sources)
    metadata = load_metadata(meta_path)
    sample_rows = list(rows(str(prefix) + ".psam"))
    ids = [r["IID"] for r in sample_rows]
    if len(set(ids)) != len(ids) or set(ids) != set(metadata):
        raise ValueError("Metadata must cover the original PSAM IDs exactly")
    eligible = [sid for sid in ids if metadata[sid]["include"] == "1" and metadata[sid]["unrelated"] == "1"]
    if len(eligible) <= args.components + 1: raise ValueError("Insufficient eligible samples")
    fids = {r["IID"]: r.get("FID", "0") for r in sample_rows}
    hq = Path(args.hq_variants or Path(args.ancient_dir) / "files/HighQualityVariants.tsv")
    lifted = hq_in_panel_coordinates(pvar, hq, Path(os.environ.get("PCA3D_RUNTIME_ROOT", configured_path("ANALYSIS_ROOT","allele3/1.ancient_pca"))) / "qc")
    snapshots = {k: state(p) for k, p in dict(pgen=str(prefix)+".pgen", pvar=pvar,
        psam=str(prefix)+".psam", metadata=meta_path, high_quality=hq).items()}
    if metadata_sources: snapshots["supplement_xlsx"] = metadata_sources["supplement_xlsx"]
    if args.exclude_bed: snapshots["exclude_bed"] = state(args.exclude_bed)
    snapshots["long_ld_bed"] = state(ROOT / "config/long_ld.GRCh38.bed")
    return dict(prefix=prefix, metadata=metadata, metadata_path=meta_path.absolute(),
                ids=ids, eligible=eligible, fids=fids, hq=Path(lifted), snapshots=snapshots)


def keep_file(path, ids, fids):
    Path(path).write_text("#FID\tIID\n" + "".join(f"{fids[s]}\t{s}\n" for s in ids))


def choose_training(source, settings):
    # Reuse the same cohort/region/time-group selection as the standard workflow.
    from evolve import select_pca_basis
    md = {}; eligible = set(source["eligible"])
    for sid, r in source["metadata"].items():
        row = dict(r, stratum=r.get("stratum", r["region"]))
        if sid not in eligible: row["include"] = "0"
        md[sid] = row
    _, basis, groups, exclusions = select_pca_basis(md, source["ids"], basis=settings["basis"],
        max_per_group=settings["cap"], time_bin_years=settings["time_bin_years"], seed=settings["seed"])
    retained = []
    for sid in basis:
        r = md[sid]
        if settings["omit_region"] and r["region"] == settings["omit_region"]: continue
        if settings["min_coverage"] and r["cohort"] == "ancient" and (r["coverage"] is None or r["coverage"] < settings["min_coverage"]): continue
        if settings["training_technology"] != "all" and r["cohort"] == "ancient" and r["technology"] != settings["training_technology"]: continue
        retained.append(sid)
    return retained


def marker_ids(hq, path, settings, exclude_bed, top_intervals):
    from evolve import _masks
    # Reduced long-LD exclusion still excludes MHC, as in the standard workflow.
    masks = _masks(exclude_bed, None, 0, settings["long_ld"], "GRCh38")
    for c, start, end in top_intervals: masks.setdefault(c, []).append((start, end))
    for c in masks:
        merged = []
        for start, end in sorted(masks[c]):
            if merged and start <= merged[-1][1]: merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else: merged.append((start, end))
        masks[c] = merged
    starts = {c: [a for a, b in v] for c, v in masks.items()}
    counts = Counter()
    partial = Path(str(path) + ".partial")
    with partial.open("w") as dest:
        for r in rows(hq):
            counts["hq_rows"] += 1
            c = r["CHROM"].removeprefix("chr")
            if c not in AUTOSOMES: continue
            c, pos = int(c), int(r["POS"]); ref, alt = r["REF"], r["ALT"]
            if pos <= 0 or ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt: continue
            if settings["omit_chrom"] == str(c): counts["omitted_chromosome"] += 1; continue
            if settings["transversions"] and {ref, alt} in ({"A", "G"}, {"C", "T"}):
                counts["transitions_removed"] += 1; continue
            j = bisect_right(starts.get(c, []), pos-1)-1
            if j >= 0 and pos-1 < masks[c][j][1]: counts["masked"] += 1; continue
            if r["ID"] in {"", "."} or any(x.isspace() for x in r["ID"]): raise ValueError("Invalid HQ marker ID")
            dest.write(r["ID"] + "\n"); counts["selected"] += 1
    partial.replace(path)
    return dict(counts), masks


def plink_commands(args, source, folder, settings, runtime):
    common, pruned, ref = runtime/"basis_common", folder/"all_pruned", folder/"ref"
    shared = [str(args.plink2), "--threads", str(args.threads), "--memory", str(args.memory_mb)]
    def cmd(*parts): return shared + [str(x) for x in parts]
    # Both genotype writes sort the lifted coordinates; variant IDs are never parsed as positions.
    return [
        ("filter", cmd("--pfile", source["prefix"], "--keep", folder/"basis.keep", "--extract", runtime/"background.ids",
            "--chr", "1-22", "--snps-only", "just-acgt", "--min-alleles", 2, "--max-alleles", 2,
            "--maf", settings["maf"], "--geno", settings["geno"], "--sort-vars", "--make-pgen", "--out", common),
            [str(common)+s for s in (".pgen", ".pvar", ".psam")]),
        ("prune", cmd("--pfile", common, "--indep-pairwise", str(settings["ld_kb"])+"kb", 1,
            settings["ld_r2"], "--out", folder/"ld"), [str(folder/"ld.prune.in")]),
        ("panel", cmd("--pfile", source["prefix"], "--keep", folder/"all.keep", "--extract", folder/"ld.prune.in",
            "--sort-vars", "--make-pgen", "--out", pruned), [str(pruned)+s for s in (".pgen", ".pvar", ".psam")]),
        ("fit", cmd("--pfile", pruned, "--keep", folder/"basis.keep", "--freq", "counts", "--pca", args.components,
            *(["approx"] if settings["n_basis"] > 5000 else ["meanimpute"]), "allele-wts", "vcols=chrom,ref,alt",
            "--seed", settings["seed"], "--out", ref), [str(ref)+s for s in (".eigenvec", ".eigenvec.allele", ".eigenval", ".acount")]),
        ("project", cmd("--pfile", pruned, "--read-freq", str(ref)+".acount", "--score", str(ref)+".eigenvec.allele",
            2, 5, "header-read", "variance-standardize", "cols=+scoresums",
            "--score-col-nums", "6-"+str(5+args.components), "--out", folder/"all"), [str(folder/"all.sscore")]),
        ("missingness", cmd("--pfile", pruned, "--missing", "sample-only", "--out", runtime/"measured"),
            [str(runtime/"measured.smiss")]),
    ]


def cached_markers(source, path, settings, exclude_bed, top_intervals, runtime):
    """Share training-independent candidate sets across MAF/basis/seed reruns."""
    identity=dict(inputs={k:v for k,v in source["snapshots"].items() if k in
                         {"pvar","high_quality","exclude_bed","long_ld_bed"}},
                  transversions=settings["transversions"],long_ld=settings["long_ld"],
                  omit_chrom=settings["omit_chrom"],top_intervals=[list(v) for v in top_intervals])
    # The digest names a small settings object; no genotype/repository hashing.
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:20]
    cache=runtime/"marker_sets";cache.mkdir(parents=True,exist_ok=True)
    ids=cache/(key+".ids");receipt=cache/(key+".json")
    old=json.loads(receipt.read_text()) if receipt.is_file() else {}
    if old.get("identity")==identity and ids.is_file() and ids.stat().st_size==old.get("bytes"):
        counts=old["counts"];masks={int(c):v for c,v in old["masks"].items()}
    else:
        counts,masks=marker_ids(source["hq"],ids,settings,exclude_bed,top_intervals)
        save_json(receipt,dict(identity=identity,counts=counts,masks=masks,bytes=ids.stat().st_size))
    path.unlink(missing_ok=True);shutil.copy2(ids, path)
    return counts,masks


def execute_scenario(args, name, settings, source, out, runtime, top_intervals=()):
    # A plan must never overwrite keeps/markers belonging to an existing fit.
    folder = out / ("plan/fits" if args.command=="study-plan" else "fits") / name
    work = runtime / ("plan/"+name if args.command=="study-plan" else name)
    folder.mkdir(parents=True, exist_ok=True); work.mkdir(parents=True, exist_ok=True)
    state_path = work / "state.json"
    identity = dict(version=VERSION, settings=settings, inputs=source["snapshots"],
                    top_intervals=[list(interval) for interval in top_intervals],
                    components=args.components, output_dir=str(folder.absolute()),
                    plink2=str(args.plink2), plink_version=getattr(args,"plink_version","unknown"))
    old = json.loads(state_path.read_text()) if state_path.is_file() else {}
    final_outputs = [folder/"all.sscore", folder/"ref.eigenvec.allele", folder/"ref.eigenval",
                     folder/"ref.eigenvec", folder/"ref.acount", folder/"all_pruned.pvar", folder/"basis.keep",
                     folder/"basis_composition.tsv", folder/"ld.prune.in", work/"measured.smiss"]
    if args.command != "study-plan" and args.replace == "FALSE" and old.get("identity") == identity and old.get("status") == "complete" and complete(final_outputs):
        ids, values = read_scores(folder/"all.sscore")
        if ids == source["eligible"] and values.shape[1] == args.components:
            print("[pca-study] reuse " + name, flush=True)
            return old
    basis = choose_training(source, settings)
    if len(basis) < max(20, args.components+2):
        result = dict(status="not_estimable", name=name, identity=identity, n_basis=len(basis),
                      reason="insufficient_training_people", folder=str(folder))
        if args.command != "study-plan": save_json(state_path, result)
        return result
    keep_file(folder/"basis.keep", basis, source["fids"])
    keep_file(folder/"all.keep", source["eligible"], source["fids"])
    composition = Counter((source["metadata"][s]["cohort"], source["metadata"][s]["region"],
                           int(source["metadata"][s]["age"]//2000)*2000,
                           source["metadata"][s]["technology"], source["metadata"][s].get("data_source","unknown")) for s in basis)
    table(folder/"basis_composition.tsv", (dict(cohort=g[0],region=g[1],start_bp=g[2],end_bp=g[2]+2000,
          technology=g[3],data_source=g[4],n=n) for g,n in sorted(composition.items())))
    counts, masks = cached_markers(source, work/"background.ids", settings, args.exclude_bed, top_intervals, runtime)
    if counts.get("selected", 0) <= args.components+1:
        result = dict(status="not_estimable", name=name, identity=identity, n_basis=len(basis),
                      marker_counts=counts, reason="insufficient_candidate_markers", folder=str(folder))
        if args.command != "study-plan": save_json(state_path, result)
        return result
    table(work/"excluded_regions.tsv", (dict(chrom=c, start0=a, end=b) for c in sorted(masks) for a, b in masks[c]))
    stages = plink_commands(args, source, folder, dict(settings, n_basis=len(basis)), work)
    plan = dict(name=name, status="planned", build="GRCh38", identity=identity, n_basis=len(basis),
        n_projected=len(source["eligible"]), marker_counts=counts, folder=str(folder),
        command_plan=[dict(stage=n, command=c, outputs=p) for n, c, p in stages])
    save_json(work/"commands.json", plan)
    (work/"commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" +
        "\n".join(shlex.join(c) for n, c, p in stages) + "\n")
    if args.command == "study-plan": return plan
    save_json(state_path, dict(plan, status="running"))
    reusable = args.replace == "FALSE" and old.get("identity") == identity
    finished = list(old.get("finished_stages", [])) if reusable else []
    try:
        for stage, command, products in stages:
            if reusable and stage in finished and complete(products): continue
            # Once one stage changes, all dependent stages are recomputed.
            reusable = False
            print("[pca-study] " + name + " / " + stage, flush=True)
            logfile = work / (stage + ".log")
            with logfile.open("w") as log:
                log.write(shlex.join(command) + "\n"); log.flush()
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            native_log = Path(command[command.index("--out")+1] + ".log")
            if native_log.is_file() and native_log != logfile:
                native_log.replace(work / (stage + ".plink.log"))
            if not complete(products): raise RuntimeError("PLINK did not produce complete " + stage + " outputs")
            if stage not in finished: finished.append(stage)
            save_json(state_path, dict(plan, status="running", finished_stages=finished))
        ids, scores = read_scores(folder/"all.sscore")
        if ids != source["eligible"] or scores.shape[1] != args.components:
            raise ValueError("Projected sample order/count or component count disagrees with plan")
        n_markers = sum(1 for _ in rows(folder/"all_pruned.pvar"))
        result = dict(plan, status="complete", finished_stages=finished, n_markers=n_markers,
                      score_path=str(folder/"all.sscore"), missingness_path=str(work/"measured.smiss"))
        save_json(state_path, result)
        # Re-creatable, large intermediates stay out of persistent result folders.
        if not args.keep_intermediates:
            for prefix in (work/"basis_common", folder/"all_pruned"):
                for ext in (".pgen", ".pgen.pgi"):
                    Path(str(prefix)+ext).unlink(missing_ok=True)
            # all_pruned PVAR/PSAM and every final scientific output are retained.
        return result
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        result = dict(plan, status="failed", finished_stages=finished, reason=str(exc), logs=str(work))
        save_json(state_path, result)
        if name == "baseline": raise RuntimeError("Baseline PCA failed; inspect " + str(work)) from exc
        return result


def geometry(reference, candidate, neighbors=15, seed=20261115, max_neighbors=5000):
    """Compare genetic coordinates only; invariant to sign/rotation/global scale."""
    x, y = np.asarray(reference, float), np.asarray(candidate, float)
    if x.shape != y.shape or x.ndim != 2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Geometry comparison requires equal finite n-by-k matrices")
    n, k = x.shape; xc, yc = x-x.mean(0), y-y.mean(0)
    if n <= k+1 or min(np.linalg.matrix_rank(xc), np.linalg.matrix_rank(yc)) < k:
        raise ValueError("Degenerate genetic coordinate subspace")
    u, singular, vt = np.linalg.svd(yc.T @ xc)
    rotation = u @ vt
    scale = float(singular.sum()/np.sum(yc*yc))
    aligned = yc @ rotation * scale + x.mean(0)
    radius = math.sqrt(float(np.mean(np.sum(xc*xc, axis=1))))
    distance = np.linalg.norm(aligned-x, axis=1)/radius
    angles = np.degrees(subspace_angles(xc, yc))
    rng = np.random.default_rng(seed)
    ix = np.sort(rng.choice(n, min(n, max_neighbors), replace=False))
    nn = min(neighbors, len(ix)-1)
    def nearest(a):
        _, idx = cKDTree(a[ix]).query(a[ix], k=nn+1)
        return [set([int(j) for j in row if j != i][:nn]) for i, row in enumerate(idx)]
    a, b = nearest(x), nearest(aligned)
    overlap = [len(aa & bb)/nn for aa, bb in zip(a, b)]
    jaccard = [len(aa & bb)/len(aa | bb) if aa | bb else 1. for aa, bb in zip(a, b)]
    pairs = min(50000, n*20)
    ii = rng.integers(n, size=pairs); jj = (ii+rng.integers(1, n, size=pairs)) % n
    dx = np.linalg.norm(x[ii]-x[jj], axis=1); dy = np.linalg.norm(y[ii]-y[jj], axis=1)
    corr = float(np.corrcoef(dx, dy)[0, 1]) if min(np.std(dx), np.std(dy)) > 1e-15 else None
    rankcorr = float(spearmanr(dx, dy).statistic) if corr is not None else None
    result = dict(n_shared=n, dimensions=k, procrustes_rms=float(np.sqrt(np.mean(distance**2))),
        normalized_displacement_median=float(np.median(distance)), normalized_displacement_p95=float(np.quantile(distance,.95)),
        principal_angle_max_degrees=float(max(angles)), principal_angle_mean_degrees=float(np.mean(angles)),
        mean_squared_subspace_cosine=float(np.mean(np.cos(np.radians(angles))**2)),
        neighbor_overlap=float(np.mean(overlap)), neighbor_jaccard=float(np.mean(jaccard)),
        neighbor_n=len(ix), neighbors=nn, distance_pearson=corr, distance_spearman=rankcorr,
        scale=scale, comparison="centered_orthogonal_rotation_and_single_scale; age_not_used")
    return result, aligned, distance, dict(rotation=rotation.tolist(), scale=scale,
        reference_mean=x.mean(0).tolist(), candidate_mean=y.mean(0).tolist())


def loading_tables(folder, destination, window_bp=1000000):
    positions = {}
    for r in rows(folder/"all_pruned.pvar"):
        if r["ID"] in positions: raise ValueError("Duplicate retained PVAR ID: " + r["ID"])
        positions[r["ID"]] = (r["CHROM"].removeprefix("chr"), int(r["POS"]), r["ALT"])
    seen = set(); marker_rows = []
    for r in rows(folder/"ref.eigenvec.allele"):
        vid = r["ID"]
        if vid not in positions: raise ValueError("Allele weights lack matching PVAR: " + vid)
        chrom, pos, alt = positions[vid]
        if r.get("A1", r.get("ALLELE")) != alt: continue
        if vid in seen: raise ValueError("Duplicate ALT loading: " + vid)
        seen.add(vid)
        pc = sorted((k for k in r if re.fullmatch(r"PC\d+", k)), key=lambda x: int(x[2:]))
        values = [float(r[k]) for k in pc]
        if not np.isfinite(values).all(): raise ValueError("Nonfinite allele loading")
        marker_rows.append(dict(id=vid, chrom=chrom, pos=pos, window_start0=(pos-1)//window_bp*window_bp,
                                **{k: v for k, v in zip(pc, values)}))
    if not marker_rows or set(positions) != seen: raise ValueError("Need one ALT weight per retained biallelic SNP")
    pc = sorted((k for k in marker_rows[0] if re.fullmatch(r"PC\d+", k)), key=lambda x: int(x[2:]))
    totals = {k: sum(r[k]**2 for r in marker_rows) for k in pc}
    windows, chromosomes = {}, {}
    for r in marker_rows:
        key = r["chrom"], r["window_start0"]
        w = windows.setdefault(key, dict(chrom=key[0], start0=key[1], end=key[1]+window_bp, n_markers=0, **{k:0. for k in pc}))
        c = chromosomes.setdefault(r["chrom"], dict(chrom=r["chrom"], n_markers=0, **{k:0. for k in pc}))
        for target in (w, c):
            target["n_markers"] += 1
            for k in pc: target[k] += r[k]**2/totals[k] if totals[k] else 0.
    window_rows = []
    for r in windows.values():
        r["PC1_PC2_share"] = (r["PC1"]+r["PC2"])/2
        window_rows.append(r)
    window_rows.sort(key=lambda r: (-r["PC1_PC2_share"], int(r["chrom"]), r["start0"]))
    for rank, r in enumerate(window_rows, 1): r["rank"] = rank
    destination.mkdir(parents=True, exist_ok=True)
    table(destination/"window_loadings.tsv", window_rows)
    table(destination/"chromosome_loadings.tsv", sorted(chromosomes.values(), key=lambda r: int(r["chrom"])))
    top = sorted(marker_rows, key=lambda r: -(r["PC1"]**2/totals["PC1"] + r["PC2"]**2/totals["PC2"]))[:1000]
    table(destination/"top_marker_loadings.tsv", top)
    return window_rows


def design_blocks(metadata, ids, missingness):
    selected = [s for s in ids if metadata[s]["cohort"] == "ancient" and metadata[s]["include"] == "1"
                and metadata[s]["unrelated"] == "1" and 0 <= metadata[s]["age"] <= 10000]
    coverage_available = any(metadata[s]["coverage"] is not None for s in selected)
    if coverage_available: selected = [s for s in selected if metadata[s]["coverage"] is not None]
    if missingness: selected = [s for s in selected if stats_number(missingness.get(s)) is not None]
    if len(selected) < 20: return selected, {}, dict(reason="insufficient_complete_cases")
    blocks = {}; age = np.array([metadata[s]["age"]/1000 for s in selected])
    age = (age-age.mean())/(age.std() or 1.)
    blocks["age"] = np.column_stack([age, age**2, age**3])
    for label, field in (("geography", "region"), ("technology", "technology"), ("data_source", "data_source")):
        values = [metadata[s].get(field, "unknown") for s in selected]
        levels = sorted(set(values))
        blocks[label] = np.column_stack([[float(x == level) for x in values] for level in levels[1:]]) if len(levels)>1 else np.empty((len(selected),0))
    if coverage_available: blocks["coverage"] = np.log1p([metadata[s]["coverage"] for s in selected])[:,None]
    if missingness: blocks["genotype_missingness"] = np.array([float(missingness[s]) for s in selected])[:,None]
    for key, values in blocks.items():
        if values.shape[1]:
            values = values-values.mean(0); sd = values.std(0)
            blocks[key] = values[:,sd>1e-12]/sd[sd>1e-12]
    return selected, blocks, dict(age_terms="cubic_polynomial_of_centered_age_in_ka; no_date_prediction",
        coverage="log1p_reported_target_coverage", modern="excluded", genotype_missingness="actual PLINK smiss",
        unavailable_fields=[k for k in ("coverage", "genotype_missingness") if k not in blocks])


def association_blocks(values, blocks):
    """Conditional block partial R2 and HC3 Wald test, allowing nuisance collinearity."""
    y = np.asarray(values, float); n = len(y)
    if y.ndim == 1: y = y[:,None]
    full = np.column_stack([np.ones(n), *blocks.values()])
    u, singular, _ = np.linalg.svd(full, full_matrices=False)
    rank = int(sum(singular > max(full.shape)*np.finfo(float).eps*singular[0]))
    qfull = u[:,:rank]; residual = y-qfull@(qfull.T@y)
    leverage = np.sum(qfull*qfull, axis=1); sse = np.sum(residual**2, axis=0)
    output = []
    for block, x in blocks.items():
        reduced = np.column_stack([np.ones(n), *(v for k,v in blocks.items() if k != block)])
        ur, sr, _ = np.linalg.svd(reduced, full_matrices=False)
        rr = int(sum(sr > max(reduced.shape)*np.finfo(float).eps*sr[0])); qr = ur[:,:rr]
        delta = rank-rr
        reduced_residual = y-qr@(qr.T@y); reduced_sse = np.sum(reduced_residual**2, axis=0)
        q = np.empty((n,0))
        if delta:
            ux, _, _ = np.linalg.svd(x-qr@(qr.T@x), full_matrices=False); q = ux[:,:delta]
        for i in range(y.shape[1]):
            base = dict(block=block, pc="PC"+str(i+1), n=n, added_df=delta, requested_df=x.shape[1],
                full_rank=rank, partial_r2=float(max(0., 1-sse[i]/reduced_sse[i])) if reduced_sse[i]>0 else None,
                status="not_estimable", p=None, reason="", model="conditional_HC3_block_Wald")
            if not delta: base["reason"] = "no_unique_information_after_other_blocks"
            elif n <= rank+5 or max(leverage) >= .999: base["reason"] = "insufficient_residual_support"
            else:
                beta = q.T@y[:,i]; weighted = q*(residual[:,i]/(1-leverage))[:,None]
                cov = weighted.T@weighted
                if np.linalg.matrix_rank(cov) < delta: base["reason"] = "singular_HC3_covariance"
                else:
                    statistic = float(beta @ np.linalg.solve(cov,beta))
                    base.update(status="estimated", wald_chi2=statistic, p=float(chi2.sf(statistic,delta)))
            output.append(base)
    return stats_adjust_fdr(output)


def metadata_associations(ids, scores, metadata, missing_path, dest):
    missingness = {}
    if missing_path and Path(missing_path).is_file():
        for r in rows(missing_path):
            value = r.get("F_MISS", r.get("F_MISS_DOSAGE"))
            if stats_number(value) is not None: missingness[r["IID"]] = value
    selected, blocks, settings = design_blocks(metadata, ids, missingness)
    if blocks:
        index = {s:i for i,s in enumerate(ids)}
        results = association_blocks(scores[[index[s] for s in selected]], blocks)
    else: results = [dict(status="not_estimable", reason="insufficient_complete_cases", n=len(selected))]
    table(dest/"metadata_associations.tsv", results)
    save_json(dest/"association_settings.json", dict(settings, n=len(selected),
        hypothesis_family="all_available_PC_by_prespecified_metadata_block", interpretation="conditional association, not causal explanation"))
    return results


def regional_changes(ids, coordinates, metadata, scenario, bootstrap, seed, min_samples=20):
    from allele_measurement import stats_seed_for
    rr = []; index = {s:i for i,s in enumerate(ids)}
    for region in REGIONS:
        regional = [s for s in ids if metadata[s]["cohort"] == "ancient" and metadata[s]["include"]=="1"
                    and metadata[s]["unrelated"]=="1" and metadata[s]["region"]==region and 0<=metadata[s]["age"]<=10000]
        for old_lo, young_lo in ((6000,4000),(4000,2000),(2000,0)):
            group = [[s for s in regional if lo <= metadata[s]["age"] < lo+2000] for lo in (old_lo,young_lo)]
            for subset in ("all", "capture_label", "shotgun_label"):
                a,b = [[s for s in g if subset=="all" or metadata[s]["technology"]==subset] for g in group]
                base = dict(scenario=scenario, region=region, region_name=REGIONS[region], subset=subset,
                    older_start_bp=old_lo, older_end_bp=old_lo+2000, younger_start_bp=young_lo, younger_end_bp=young_lo+2000,
                    primary_contrast=int(old_lo==6000), n_older=len(a), n_younger=len(b), direction="younger_minus_older",
                    uncertainty="individual_bootstrap_conditional_on_fixed_genetic_coordinates",
                    status="estimated" if min(len(a),len(b))>=min_samples else "not_estimable")
                delta=lo=hi=None
                if base["status"]=="estimated":
                    x=coordinates[[index[s] for s in a],:2]; y=coordinates[[index[s] for s in b],:2]
                    ax=stats_bootstrap_mean(x,bootstrap,stats_seed_for(seed,(region,old_lo,subset,"old")))[3]
                    by=stats_bootstrap_mean(y,bootstrap,stats_seed_for(seed,(region,old_lo,subset,"young")))[3]
                    delta=y.mean(0)-x.mean(0); lo,hi=np.quantile(by-ax,[.025,.975],axis=0)
                for i in range(2):
                    rr.append(dict(base,pc="PC"+str(i+1),estimate=None if delta is None else float(delta[i]),
                                   ci_low=None if lo is None else float(lo[i]),ci_high=None if hi is None else float(hi[i])))
    return rr


def fixed_composition(ids, scores, metadata, bootstrap, seed, min_samples):
    primary = [i for i,s in enumerate(ids) if metadata[s]["cohort"]=="ancient" and metadata[s]["include"]=="1"
               and metadata[s]["unrelated"]=="1" and 0<=metadata[s]["age"]<=10000]
    result=[]
    if not primary: return [dict(status="not_estimable", reason="no_eligible_ancient_samples")]
    for label,field in (("fixed_geography","region"),("fixed_technology","technology")):
        counts=Counter(metadata[ids[i]][field] for i in primary)
        weights={g:n/len(primary) for g,n in counts.items()}
        values=stats_fixed_composition(scores[primary,:2],
            [min(4,int(metadata[ids[i]]["age"]//2000)) for i in primary],
            [metadata[ids[i]][field] for i in primary], weights, range(5), min_samples, bootstrap, seed)
        for row in values:
            for i in range(2):
                rr={k:v for k,v in row.items() if k not in {"estimate","ci_low","ci_high","se"}}
                rr.update(measure=label,pc="PC"+str(i+1),start_bp=row["time_bin"]*2000,end_bp=(row["time_bin"]+1)*2000,
                    reference_weights=weights, **{k:None if row[k] is None else float(row[k][i]) for k in ("estimate","ci_low","ci_high","se")})
                result.append(rr)
    return result


def run_drivers(args, out, records, metadata):
    baseline=records.get("baseline",{})
    if baseline.get("status")!="complete": raise ValueError("A completed baseline study is required")
    folder=Path(baseline["folder"]); dest=out/"drivers"; dest.mkdir(exist_ok=True)
    windows=loading_tables(folder,dest,args.window_bp)
    ids,scores=read_scores(folder/"all.sscore")
    metadata_associations(ids,scores,metadata,baseline.get("missingness_path"),dest)
    table(dest/"regional_changes.tsv",regional_changes(ids,scores,metadata,"baseline",args.bootstrap,args.seed,args.min_samples))
    table(dest/"fixed_composition.tsv",fixed_composition(ids,scores,metadata,args.bootstrap,args.seed,args.min_samples))
    return windows


def id_set(path, field=None):
    if field: return {r[field] for r in rows(path)}
    with Path(path).open() as h: return {line.strip() for line in h if line.strip()}


def set_overlap(a,b):
    return len(a & b)/len(a | b) if a | b else None


def eigenvalue_gap(folder, k):
    """Relative separation at the plotted subspace boundary, not variance explained."""
    values=[float(x) for x in (folder/"ref.eigenval").read_text().split()]
    if len(values)<=k or values[k-1]<=0:return None
    return (values[k-1]-values[k])/values[k-1]


def run_stability(args, out, records, metadata):
    base=records.get("baseline",{})
    if base.get("status")!="complete": raise ValueError("Baseline has not completed")
    ref_ids,x=read_scores(Path(base["folder"])/"all.sscore")
    ref_index={s:i for i,s in enumerate(ref_ids)}
    ref_markers=id_set(Path(base["folder"])/"ld.prune.in")
    ref_training=id_set(Path(base["folder"])/"basis.keep","IID")
    summary=[]; changed=[]; group_rows=[]; contrasts=[]
    for name, record in records.items():
        if record.get("status")!="complete":
            summary.append(dict(scenario=name,status=record.get("status","missing"),reason=record.get("reason","")));continue
        ids,y=read_scores(Path(record["folder"])/"all.sscore")
        lookup={s:i for i,s in enumerate(ids)}; common=[s for s in ref_ids if s in lookup]
        ix=[ref_index[s] for s in common]; iy=[lookup[s] for s in common]
        markers=id_set(Path(record["folder"])/"ld.prune.in")
        training=id_set(Path(record["folder"])/"basis.keep","IID")
        folder=out/"aligned"/name;folder.mkdir(parents=True,exist_ok=True)
        for k in sorted({2,*[v for v in (5,10) if v<=min(x.shape[1],y.shape[1])]}):
            try:
                metrics,aligned,displacement,transform=geometry(x[ix,:k],y[iy,:k],args.neighbors,args.seed,args.max_neighbor_samples)
            except ValueError as exc:
                summary.append(dict(scenario=name,dimensions=k,status="not_estimable",reason=str(exc)));continue
            table(folder/("coordinates_k"+str(k)+".tsv.gz"),
                (dict(sample_id=s,**{"PC"+str(i+1):float(v) for i,v in enumerate(aligned[j])}) for j,s in enumerate(common)))
            save_json(folder/("transform_k"+str(k)+".json"),transform)
            difference={key:v for key,v in record.get("identity",{}).get("settings",{}).items()
                        if v!=base.get("identity",{}).get("settings",{}).get(key)}
            summary.append(dict(scenario=name,status="estimated",**metrics,
                n_markers=len(markers),n_training=len(training),marker_jaccard=set_overlap(ref_markers,markers),
                training_jaccard=set_overlap(ref_training,training),changed_settings=difference,
                reference_eigenvalue_boundary_gap=eigenvalue_gap(Path(base["folder"]),k),
                candidate_eigenvalue_boundary_gap=eigenvalue_gap(Path(record["folder"]),k)))
            if k==2:
                table(folder/"displacement.tsv.gz",(dict(sample_id=s,normalized_displacement=float(displacement[j]),
                    cohort=metadata[s]["cohort"],region=metadata[s]["region"],age_bp=metadata[s]["age"],
                    technology=metadata[s]["technology"],source=metadata[s].get("data_source","unknown")) for j,s in enumerate(common)))
                for field in ("region","technology","data_source"):
                    for level in sorted({metadata[s].get(field,"unknown") for s in common}):
                        indices=[j for j,s in enumerate(common) if metadata[s].get(field,"unknown")==level]
                        group_rows.append(dict(scenario=name,field=field,group=level,n=len(indices),
                            median_displacement=float(np.median(displacement[indices])),p95_displacement=float(np.quantile(displacement[indices],.95))))
                contrasts.extend(regional_changes(common,aligned,metadata,name,args.bootstrap,args.seed,args.min_samples))
                changed.append(dict(scenario=name,procrustes_rms=metrics["procrustes_rms"],
                    marker_jaccard=set_overlap(ref_markers,markers),training_jaccard=set_overlap(ref_training,training),
                    tested_change=difference,interpretation="sensitivity to the declared intervention; not a causal attribution",
                    reference_PC2_PC3_relative_gap=eigenvalue_gap(Path(base["folder"]),2),
                    candidate_PC2_PC3_relative_gap=eigenvalue_gap(Path(record["folder"]),2),
                    review_geometry=int(metrics["procrustes_rms"]>.15 or metrics["neighbor_overlap"]<.7),
                    review_rule="descriptive only: RMS>0.15 or neighbor overlap<0.70; not a significance test"))
    table(out/"stability.tsv",summary);table(out/"shape_changes.tsv",changed)
    table(out/"group_displacement.tsv",group_rows);table(out/"regional_changes.tsv",contrasts)
    return summary


def time_stability_payload(out, records, metadata):
    """All aligned genetic points plus unchanged supplied dates, stored once."""
    completed={name:r for name,r in records.items() if r.get("status")=="complete"}
    if "baseline" not in completed: raise ValueError("3D stability requires a completed baseline")
    base_path=out/"aligned/baseline/coordinates_k2.tsv.gz"
    if not base_path.is_file():
        raise ValueError("Run stability before rendering 3D; aligned baseline coordinates are missing")
    base=list(rows(base_path));ids=[r["sample_id"] for r in base]
    if not ids or len(set(ids))!=len(ids) or set(ids)-set(metadata):
        raise ValueError("3D baseline IDs must be unique and covered by metadata")
    age=[metadata[s]["age"]/1000 for s in ids]
    if not np.isfinite(age).all() or min(age)<0: raise ValueError("Invalid supplied dates for 3D stability")
    groups={key:[i for i,s in enumerate(ids) if metadata[s]["cohort"]=="ancient" and metadata[s]["region"]==key] for key in REGIONS}
    groups["modern"]=[i for i,s in enumerate(ids) if metadata[s]["cohort"]=="modern"]
    if sum(map(len,groups.values()))!=len(ids): raise ValueError("3D cohort/region groups do not cover all samples")
    hover=[[s,metadata[s]["cohort"],metadata[s]["region"],metadata[s].get("data_source","unknown"),
            metadata[s]["technology"],metadata[s].get("coverage") if metadata[s].get("coverage") is not None else "NA"] for s in ids]
    scenarios=[];limits=[]
    for name,r in completed.items():
        path=out/"aligned"/name/"coordinates_k2.tsv.gz"
        if not path.is_file(): raise ValueError("Missing aligned 3D coordinates for completed scenario: "+name)
        points=list(rows(path));lookup={p["sample_id"]:p for p in points}
        if len(lookup)!=len(points) or set(lookup)!=set(ids):
            raise ValueError("3D scenarios must contain exactly the same sample IDs: "+name)
        xy=np.array([[float(lookup[s]["PC1"]),float(lookup[s]["PC2"])] for s in ids])
        if not np.isfinite(xy).all():raise ValueError("Nonfinite aligned 3D coordinates: "+name)
        limits.append([xy.min(0),xy.max(0)])
        centroids={}
        for region in ("N","C"):
            values=[]
            for lo,hi in ((6,8),(4,6)):
                ix=[i for i in groups[region] if lo<=age[i]<hi]
                if ix:
                    values.append(dict(period=f"{lo}–{hi} ka",n=len(ix),pc1=float(xy[ix,0].mean()),pc2=float(xy[ix,1].mean()),
                                       age_ka=float(np.mean([age[i] for i in ix]))))
            centroids[region]=values
        changed={key:v for key,v in r.get("identity",{}).get("settings",{}).items()
                 if v!=completed["baseline"].get("identity",{}).get("settings",{}).get(key)}
        scenarios.append(dict(name=name,pc1=xy[:,0].tolist(),pc2=xy[:,1].tolist(),centroids=centroids,
            n_markers=r.get("n_markers"),n_training=r.get("n_basis"),changed_settings=changed))
    extrema=np.array(limits);low=extrema[:,0,:].min(0);high=extrema[:,1,:].max(0)
    margin=np.maximum((high-low)*.04,1e-5)
    return dict(schema="pca-time-stability-1",n_samples=len(ids),sample_ids=ids,age_ka=age,hover=hover,groups=groups,
        scenarios=scenarios,region_names={**REGIONS,"modern":"Modern"},
        colors=dict(zip([*REGIONS,"modern"],["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#777777"])),
        axis_ranges=dict(pc1=[float(low[0]-margin[0]),float(high[0]+margin[0])],
                         pc2=[float(low[1]-margin[1]),float(high[1]+margin[1])],age_ka=[0.,max(1.,max(age)*1.03)]),
        camera=dict(eye=dict(x=1.55,y=1.55,z=.85)),
        age_source="Supplied metadata age_bp / 1000; identical in every scenario, never inferred from PCA",
        centroids="Arithmetic means of actual PC1, PC2 and supplied age in [6,8) and [4,6) ka; labels show actual n",
        alignment="Saved k=2 genetic alignment: rotation/reflection, translation and one global scale; age excluded")


def write_time_stability(out, records, metadata, make_static=True):
    try:
        from plotly.offline import get_plotlyjs
    except ImportError as exc:
        raise RuntimeError("Interactive PCA stability requires plotly (included in env.yml). Install that environment, or use --no-plots for tables only.") from exc
    payload=time_stability_payload(out,records,metadata)
    encoded=json.dumps(payload,ensure_ascii=False,separators=(",",":"),allow_nan=False).replace("<","\\u003c")
    template=r'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Ancient PCA: 3D shape stability</title>
<style>body{font:15px system-ui;margin:24px;color:#18354a}h1{margin-bottom:8px}.controls{display:flex;gap:22px;flex-wrap:wrap;align-items:center;padding:14px;background:#eef3f7;border-radius:8px}select,button{font:inherit;padding:6px}#plot{width:100%;height:76vh;min-height:580px}#details{font-size:13px;white-space:pre-wrap}.note{max-width:1200px;line-height:1.5}</style>
<h1>Ancient PCA: 3D shape stability</h1>
<p class="note">Baseline and candidate share every sample, the supplied dates, axis ranges and camera. Only the aligned genetic PC1/PC2 coordinates change. This is a time-indexed genetic visualization; the vertical axis is not a third genetic PC or a predicted date.</p>
<div class="controls"><label>Candidate <select id="scenario"></select></label>
<label>6–8 → 4–6 ka centroids <select id="centroids"><option value="both">North + Central</option><option value="N">North</option><option value="C">Central</option><option value="none">Off</option></select></label>
<button id="reset" type="button">Reset shared camera</button><span id="counts"></span></div>
<div id="plot"></div><p id="details"></p>
<p class="note">Drag either panel to rotate both. The scenario menu preserves the camera. Click a region in the shared legend to hide/show it in both panels. Modern controls are a separate gray group. Centroid heights are the actual mean supplied dates, not arbitrary period midpoints; lines connect descriptive period means and do not trace individual ancestry or estimate selection. Hover for sample ID, date, region, technology, source and coverage. All data and Plotly are embedded; no network is needed.</p>
<p><a href="report.html">Study tables and figures</a> · <a href="figures/pca_time_stability.png">Static 3D PNG</a></p>
<script id="pca-study-data" type="application/json">__PAYLOAD__</script>
<script>__PLOTLY_BUNDLE__</script>
<script>
const data=JSON.parse(document.getElementById('pca-study-data').textContent);
const target=document.getElementById('plot'), menu=document.getElementById('scenario'), centroidMenu=document.getElementById('centroids');
const baseline=data.scenarios.find(s=>s.name==='baseline');
let camera=JSON.parse(JSON.stringify(data.camera)), syncing=false, initialized=false;
data.scenarios.forEach(s=>{const option=document.createElement('option');option.value=s.name;option.textContent=s.name;menu.appendChild(option)});
menu.value=(data.scenarios.find(s=>s.name!=='baseline')||baseline).name;
function scene(domain){return {domain:{x:domain,y:[0,1]},camera:camera,aspectmode:'manual',aspectratio:{x:1,y:1,z:1.2},
 xaxis:{title:{text:'Aligned genetic PC1'},range:data.axis_ranges.pc1,autorange:false},
 yaxis:{title:{text:'Aligned genetic PC2'},range:data.axis_ranges.pc2,autorange:false},
 zaxis:{title:{text:'Supplied age (ka BP)'},range:data.axis_ranges.age_ka,autorange:false}}}
function traces(s,which,legend){
 const result=Object.entries(data.groups).map(([group,ix])=>({type:'scatter3d',mode:'markers',scene:which,name:data.region_names[group],legendgroup:group,showlegend:legend,
  x:ix.map(i=>s.pc1[i]),y:ix.map(i=>s.pc2[i]),z:ix.map(i=>data.age_ka[i]),customdata:ix.map(i=>data.hover[i]),
  marker:{size:group==='modern'?2.5:2.1,color:data.colors[group],opacity:group==='modern'?.65:.45},
  hovertemplate:'%{customdata[0]}<br>PC1 %{x:.5f} · PC2 %{y:.5f}<br>Supplied age %{z:.3f} ka BP<br>%{customdata[1]} · region %{customdata[2]}<br>Source: %{customdata[3]}<br>Technology: %{customdata[4]}<br>Target coverage: %{customdata[5]}<extra></extra>'}));
 for(const region of ['N','C']){
  const cc=s.centroids[region]||[],choice=centroidMenu.value;
  result.push({type:'scatter3d',scene:which,mode:'lines+markers+text',name:data.region_names[region]+' period means',showlegend:false,
   visible:choice==='both'||choice===region,x:cc.map(c=>c.pc1),y:cc.map(c=>c.pc2),z:cc.map(c=>c.age_ka),
   text:cc.map(c=>region+' '+c.period+' (n='+c.n+')'),textposition:'top center',textfont:{size:11,color:data.colors[region]},
   customdata:cc.map(c=>[c.period,c.n]),line:{color:data.colors[region],width:7},marker:{color:data.colors[region],size:6,symbol:'diamond'},
   hovertemplate:data.region_names[region]+' %{customdata[0]}<br>n=%{customdata[1]}<br>Mean PC1 %{x:.5f} · PC2 %{y:.5f}<br>Mean supplied age %{z:.3f} ka BP<extra></extra>'});
 }
 return result;
}
async function render(){
 const candidate=data.scenarios.find(s=>s.name===menu.value),settings=Object.entries(candidate.changed_settings).map(([k,v])=>k+'='+JSON.stringify(v)).join(', ');
 const layout={scene:scene([0,.47]),scene2:scene([.53,1]),margin:{l:0,r:0,t:40,b:0},uirevision:'shared-genetic-time-view',
  legend:{orientation:'h',x:0,y:1.08,groupclick:'togglegroup'},paper_bgcolor:'#ffffff',
  annotations:[{text:'Baseline',xref:'paper',yref:'paper',x:.23,y:1.02,showarrow:false,font:{size:17}},
               {text:candidate.name,xref:'paper',yref:'paper',x:.77,y:1.02,showarrow:false,font:{size:17}}]};
 await Plotly.react(target,[...traces(baseline,'scene',true),...traces(candidate,'scene2',false)],layout,
  {responsive:true,scrollZoom:true,displaylogo:false,toImageButtonOptions:{format:'png',filename:'pca_time_stability',width:1600,height:900,scale:2}});
 document.getElementById('counts').textContent=data.n_samples.toLocaleString()+' people in each panel · '+data.scenarios.length+' fitted scenarios';
 document.getElementById('details').textContent='Candidate: '+candidate.name+' | retained SNPs: '+candidate.n_markers+' | training people: '+candidate.n_training+'\nChanged settings: '+(settings||'none (baseline)');
 target.dataset.ready='true';target.dataset.scenario=candidate.name;
 if(!initialized){initialized=true;target.on('plotly_relayout',async event=>{
   const moved=event['scene.camera']||event['scene2.camera'];if(!moved||syncing)return;
   syncing=true;camera=JSON.parse(JSON.stringify(moved));try{await Plotly.relayout(target,{'scene.camera':camera,'scene2.camera':camera})}finally{syncing=false}
 });}
}
menu.addEventListener('change',render);centroidMenu.addEventListener('change',render);
document.getElementById('reset').addEventListener('click',async()=>{camera=JSON.parse(JSON.stringify(data.camera));await Plotly.relayout(target,{'scene.camera':camera,'scene2.camera':camera})});
render().catch(error=>{target.dataset.error=error.message;target.textContent='3D rendering failed: '+error.message;console.error(error)});
</script></html>'''
    page=template.replace("__PAYLOAD__",encoded,1).replace("__PLOTLY_BUNDLE__",get_plotlyjs(),1)
    path=out/"pca_time_stability.html";partial=path.with_suffix(".html.partial")
    partial.write_text(page);partial.replace(path)
    if make_static:
        import matplotlib.pyplot as plt
        summary=list(rows(out/"stability.tsv")) if (out/"stability.tsv").is_file() else []
        ranked=sorted((r for r in summary if r.get("dimensions")=="2" and r.get("status")=="estimated" and r["scenario"]!="baseline"),key=lambda r:-float(r["procrustes_rms"]))
        selected=ranked[0]["scenario"] if ranked else payload["scenarios"][-1]["name"]
        base=next(s for s in payload["scenarios"] if s["name"]=="baseline");candidate=next(s for s in payload["scenarios"] if s["name"]==selected)
        fig=plt.figure(figsize=(14,7),layout="constrained")
        for number,s in enumerate((base,candidate),1):
            ax=fig.add_subplot(1,2,number,projection="3d")
            for group,ix in payload["groups"].items():
                ax.scatter([s["pc1"][i] for i in ix],[s["pc2"][i] for i in ix],[payload["age_ka"][i] for i in ix],
                    s=2,alpha=.35,color=payload["colors"][group],label=payload["region_names"][group],rasterized=True)
            for region in ("N","C"):
                cc=s["centroids"][region]
                ax.plot([c["pc1"] for c in cc],[c["pc2"] for c in cc],[c["age_ka"] for c in cc],"o-",color=payload["colors"][region],lw=2,ms=4)
            ax.set(xlim=payload["axis_ranges"]["pc1"],ylim=payload["axis_ranges"]["pc2"],zlim=payload["axis_ranges"]["age_ka"],
                xlabel="Aligned genetic PC1",ylabel="Aligned genetic PC2",zlabel="Supplied age (ka BP)",title=s["name"])
            ax.view_init(elev=20,azim=42);ax.set_box_aspect((1,1,1.2))
            if number==1:ax.legend(loc="upper left",fontsize=8,markerscale=3,frameon=False)
        fig.suptitle("3D shape comparison: fixed dates, axis ranges and view\nNorth/Central lines connect actual 6–8 ka and 4–6 ka period means",fontsize=14)
        (out/"figures").mkdir(exist_ok=True)
        fig.savefig(out/"figures"/("pca_time_stability.png"),dpi=160)
        plt.close(fig)
    return dict(path=str(path),n_samples=payload["n_samples"],scenarios=len(payload["scenarios"]))


def figures_and_report(out, records, metadata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figs=out/"figures";figs.mkdir(exist_ok=True);images=[]
    three_d=write_time_stability(out,records,metadata)
    images.append("figures/pca_time_stability.png")
    summary=list(rows(out/"stability.tsv")) if (out/"stability.tsv").is_file() else []
    good=[r for r in summary if r.get("status")=="estimated" and r.get("dimensions")=="2"]
    if good:
        fig,axes=plt.subplots(1,2,figsize=(13, max(4,.28*len(good)+1)),constrained_layout=True)
        labels=[r["scenario"] for r in good];yy=np.arange(len(labels))
        axes[0].barh(yy,[float(r["procrustes_rms"]) for r in good]);axes[0].set(yticks=yy,yticklabels=labels,xlabel="Normalized Procrustes RMS (lower = more similar)")
        axes[1].barh(yy,[float(r["neighbor_overlap"]) for r in good]);axes[1].set(yticks=yy,yticklabels=labels,xlim=(0,1),xlabel="Nearest-neighbor overlap (higher = more similar)")
        fig.suptitle("Genetic geometry only: time is not part of the comparison")
        fig.savefig(figs/("stability.png"),dpi=160)
        plt.close(fig);images.append("figures/stability.png")
        chosen=sorted((r for r in good if r["scenario"]!="baseline"),key=lambda r:-float(r["procrustes_rms"]))[:5]
        chosen=[{"scenario":"baseline"},*chosen]
        fig,axes=plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
        colors=dict(zip(REGIONS,["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00"]))
        limits=[]
        for ax,item in zip(axes.flat,chosen):
            name=item["scenario"];points=list(rows(out/"aligned"/name/"coordinates_k2.tsv.gz"))
            for region in REGIONS:
                pp=[p for p in points if metadata[p["sample_id"]]["region"]==region and metadata[p["sample_id"]]["cohort"]=="ancient"]
                ax.scatter([float(p["PC1"]) for p in pp],[float(p["PC2"]) for p in pp],s=2,alpha=.35,color=colors[region],label=REGIONS[region],rasterized=True)
            pp=[p for p in points if metadata[p["sample_id"]]["cohort"]=="modern"]
            ax.scatter([float(p["PC1"]) for p in pp],[float(p["PC2"]) for p in pp],s=3,alpha=.5,color="#666666",label="Modern",rasterized=True)
            limits.extend((float(p["PC1"]),float(p["PC2"])) for p in points)
            ax.set(title=name,xlabel="Aligned genetic PC1",ylabel="Aligned genetic PC2")
        if limits:
            arr=np.array(limits);lo=arr.min(0);hi=arr.max(0);margin=np.maximum((hi-lo)*.04,1e-5)
            for ax in axes.flat:
                ax.set_xlim(lo[0]-margin[0],hi[0]+margin[0]);ax.set_ylim(lo[1]-margin[1],hi[1]+margin[1]);ax.set_aspect("equal")
        for ax in list(axes.flat)[len(chosen):]:ax.axis("off")
        axes.flat[0].legend(frameon=False,fontsize=8,markerscale=3,ncol=2)
        fig.suptitle("Baseline and largest observed geometry differences, after alignment")
        fig.savefig(figs/("aligned_comparison.png"),dpi=160)
        plt.close(fig);images.append("figures/aligned_comparison.png")
    regional_path=out/"drivers/regional_changes.tsv"
    if regional_path.is_file():
        rr=[r for r in rows(regional_path) if r.get("primary_contrast")=="1"]
        if any(r.get("status")=="estimated" for r in rr):
            fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
            groups=list(REGIONS)
            for ax,pc in zip(axes,("PC1","PC2")):
                for subset,offset,color,label in (("all",-.20,"#263745","All"),
                        ("capture_label",0,"#0072B2","Capture"),("shotgun_label",.20,"#D55E00","Shotgun")):
                    first=True
                    for r in rr:
                        if r["status"]!="estimated" or r["pc"]!=pc or r["subset"]!=subset:continue
                        y=groups.index(r["region"])+offset
                        ax.hlines(y,float(r["ci_low"]),float(r["ci_high"]),color=color,lw=1.6)
                        ax.plot(float(r["estimate"]),y,"o",ms=5,color=color,label=label if first else None);first=False
                ax.axvline(0,color="#aaaaaa",lw=1)
                counts={r["region"]:r for r in rr if r["subset"]=="all" and r["pc"]==pc}
                labels=[REGIONS[g]+" ("+counts[g]["n_older"]+"→"+counts[g]["n_younger"]+")" if g in counts else REGIONS[g] for g in groups]
                ax.set(yticks=range(len(groups)),yticklabels=labels,ylabel="Region (older → younger n)",xlabel=pc+": younger − older",ylim=(-.6,len(groups)-.4))
                ax.invert_yaxis();ax.grid(axis="x",alpha=.2);ax.legend(frameon=False,fontsize=9)
            fig.suptitle("6–8 ka to 4–6 ka: regional changes\n95% individual-bootstrap CIs conditional on the fitted genetic coordinates")
            fig.savefig(figs/("regional_6to4ka.png"),dpi=160)
            plt.close(fig);images.append("figures/regional_6to4ka.png")
    association_path=out/"drivers/metadata_associations.tsv"
    if association_path.is_file():
        rr=[r for r in rows(association_path) if r.get("status")=="estimated"]
        if rr:
            blocks=list(dict.fromkeys(r["block"] for r in rr));pcs=sorted({r["pc"] for r in rr},key=lambda k:int(k[2:]))
            arr=np.full((len(pcs),len(blocks)),np.nan)
            for r in rr:arr[pcs.index(r["pc"]),blocks.index(r["block"])]=float(r["partial_r2"])
            fig,ax=plt.subplots(figsize=(max(9,len(blocks)*1.6),max(3,len(pcs)*.48)),constrained_layout=True)
            cmap=plt.get_cmap("Blues").copy();cmap.set_bad("#eeeeee")
            im=ax.imshow(arr,aspect="auto",vmin=0,vmax=max(.05,float(np.nanmax(arr))),cmap=cmap)
            for i in range(len(pcs)):
                for j in range(len(blocks)):
                    ax.text(j,i,"NA" if not np.isfinite(arr[i,j]) else f"{100*arr[i,j]:.1f}%",ha="center",va="center",fontsize=9,
                            color="white" if np.isfinite(arr[i,j]) and arr[i,j]>np.nanmax(arr)*.6 else "#162735")
            ax.set(xticks=range(len(blocks)),xticklabels=[s.replace("_"," ") for s in blocks],yticks=range(len(pcs)),yticklabels=pcs,
                   title="Conditional partial R² by metadata block (ancient samples, 0–10 ka)\nAssociation after other blocks; no causal interpretation")
            plt.setp(ax.get_xticklabels(),rotation=15,ha="right");fig.colorbar(im,ax=ax,label="Partial R²")
            fig.savefig(figs/("metadata_associations.png"),dpi=160)
            plt.close(fig);images.append("figures/metadata_associations.png")
    composition_path=out/"drivers/fixed_composition.tsv"
    if composition_path.is_file():
        rr=[r for r in rows(composition_path) if r.get("status")=="estimated"]
        if rr:
            fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
            for ax,pc in zip(axes,("PC1","PC2")):
                for measure,color in (("fixed_geography","#0072B2"),("fixed_technology","#D55E00")):
                    pp=sorted((r for r in rr if r["measure"]==measure and r["pc"]==pc),key=lambda r:float(r["start_bp"]))
                    if not pp:continue
                    xx=[(float(r["start_bp"])+float(r["end_bp"]))/2000 for r in pp]
                    yy=[float(r["estimate"]) for r in pp]
                    # Connect adjacent supported bins only; absent strata are never bridged.
                    for i in range(1,len(pp)):
                        if int(pp[i]["start_bp"])-int(pp[i-1]["start_bp"])==2000:
                            ax.plot(xx[i-1:i+1],yy[i-1:i+1],color=color,lw=1)
                    for x,r in zip(xx,pp):ax.vlines(x,float(r["ci_low"]),float(r["ci_high"]),color=color,lw=1.6)
                    ax.scatter(xx,yy,color=color,s=22,label=measure.replace("_"," "))
                ax.set(xlabel="Age (ka BP; older at left)",ylabel=pc,title=pc+" at fixed composition",xlim=(10,0));ax.legend(frameon=False);ax.grid(alpha=.2)
            fig.suptitle("Whole-window reference weights; unsupported periods remain missing")
            fig.savefig(figs/("fixed_composition.png"),dpi=160)
            plt.close(fig);images.append("figures/fixed_composition.png")
    window_path=out/"drivers/window_loadings.tsv"
    if window_path.is_file():
        ww=list(rows(window_path))[:20]
        fig,ax=plt.subplots(figsize=(11,6),constrained_layout=True)
        ax.barh(np.arange(len(ww)),[100*float(r["PC1_PC2_share"]) for r in ww])
        ax.set(yticks=np.arange(len(ww)),yticklabels=[r["chrom"]+":"+str(int(r["start0"])+1)+"–"+r["end"] for r in ww],
               xlabel="Mean PC1/PC2 squared-loading share (%)",title="Most influential marker windows (not causal loci)")
        ax.invert_yaxis();fig.savefig(figs/"loading_windows.png",dpi=160);plt.close(fig)
        images.append("figures/loading_windows.png")
    status="".join("<tr><td>"+html.escape(n)+"</td><td>"+html.escape(r.get("status","missing"))+"</td><td>"+html.escape(r.get("reason",""))+"</td></tr>" for n,r in records.items())
    links=["stability.tsv","shape_changes.tsv","group_displacement.tsv","regional_changes.tsv",
           "drivers/window_loadings.tsv","drivers/chromosome_loadings.tsv","drivers/metadata_associations.tsv","drivers/fixed_composition.tsv"]
    page='<!doctype html><html lang="en"><meta charset="utf-8"><title>Ancient PCA study</title><style>body{font:16px system-ui;max-width:1300px;margin:36px auto;padding:0 24px;color:#18354a}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd}img{max-width:100%}li{margin:9px 0}</style><h1>Ancient PCA: stability and drivers</h1><p>Real PLINK reruns of the GRCh38 panel. The reported geometry comparisons exclude age. Alignment removes signs, rotations and one global scale, while retaining changes in shape.</p>'
    page+="<table><tr><th>Scenario</th><th>Execution</th><th>Reason</th></tr>"+status+"</table>"
    page+='<p><strong><a href="pca_time_stability.html">Open interactive 3D shape stability</a></strong> — baseline/candidate panels, scenario switching, synchronized camera and fixed axis ranges; all '+str(three_d["n_samples"])+" people retained.</p>"
    page+="".join('<p><img src="'+p+'" alt="'+p+'"></p>' for p in images)
    page+="<h2>Tables</h2><ul>"+"".join('<li><a href="'+p+'">'+p+'</a></li>' for p in links if (out/p).is_file())+"</ul>"
    page+="<h2>Interpretation</h2><ul>"+"".join("<li>"+html.escape(t)+"</li>" for t in LIMITATIONS)+"</ul></html>"
    (out/"report.html").write_text(page)


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    p.add_argument("command",choices=["study","study-plan","stability","drivers"])
    p.add_argument("--out",type=Path,default=Path(configured_path("DIR_GEN","ancient/38/pca3d/study")))
    p.add_argument("--study-dir",type=Path,help="Existing study for stability/drivers (also selects its output)")
    p.add_argument("--runtime-root",type=Path,default=Path(configured_path("ANALYSIS_ROOT","allele3/1.ancient_pca/tmp/study")))
    p.add_argument("--ancient-dir",type=Path,default=Path(configured_path("DIR_GEN","ancient")))
    p.add_argument("--pgen-prefix",type=Path);p.add_argument("--metadata",type=Path);p.add_argument("--hq-variants",type=Path)
    p.add_argument("--supplement-xlsx",type=Path);p.add_argument("--exclude-bed","--pca-exclude-bed",dest="exclude_bed",type=Path)
    p.add_argument("--scenarios",default="core",help="core, full, or comma-separated scenario names; baseline is always first")
    p.add_argument("--scenario-file",type=Path,help="Optional JSON: scenario_name -> overrides of the documented settings")
    p.add_argument("--plink2",default="plink2")
    p.add_argument("--threads","--pca-threads",type=int,default=4)
    p.add_argument("--memory-mb","--pca-memory-mb",type=int,default=16384)
    p.add_argument("--components","--pca-components",type=int,default=10)
    p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--bootstrap",type=int,default=500);p.add_argument("--min-samples",type=int,default=20)
    p.add_argument("--neighbors",type=int,default=15);p.add_argument("--max-neighbor-samples",type=int,default=5000)
    p.add_argument("--window-bp",type=int,default=1000000)
    p.add_argument("--replace",choices=["TRUE","FALSE"],default="FALSE",help="FALSE resumes matching completed stages; TRUE reruns selected scenarios")
    p.add_argument("--keep-intermediates",action="store_true")
    p.add_argument("--no-plots",action="store_true")
    return p


def main(argv=None):
    args=parser().parse_args(argv)
    if min(args.threads,args.components,args.neighbors,args.window_bp)<1 or args.min_samples<2 or args.components<2 or args.bootstrap<2 or args.memory_mb<640 or args.max_neighbor_samples<3:
        raise ValueError("Require >=2 PCs/bootstrap/minimum samples, positive thresholds and >=640 MiB PLINK memory")
    out=(args.study_dir or args.out).absolute();runtime=args.runtime_root.absolute()
    out.mkdir(parents=True,exist_ok=True);runtime.mkdir(parents=True,exist_ok=True)
    with (runtime/"study.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        manifest_path=out/("study_plan.json" if args.command=="study-plan" else "study.json")
        if args.command in {"stability","drivers"}:
            if not manifest_path.is_file(): raise ValueError("Run study first; no study.json at "+str(out))
            saved=json.loads(manifest_path.read_text());records=saved["scenarios"]
            metadata=load_metadata(args.metadata or saved["metadata"])
            if args.command=="drivers":run_drivers(args,out,records,metadata)
            else:run_stability(args,out,records,metadata)
        else:
            settings=scenario_settings(args);source=source_inputs(args,out,runtime)
            if args.command!="study-plan":
                version=subprocess.run([str(args.plink2),"--version"],capture_output=True,text=True,check=True).stdout.strip()
            else:version="not_executed_plan_only"
            args.plink_version=version
            records={};windows=[]
            save_json(manifest_path,dict(version=VERSION,status="running" if args.command=="study" else "planned",
                metadata=str(source["metadata_path"]),parameters=vars(args),plink_version=version,scenarios=records))
            for name,s in settings.items():
                top=[(int(w["chrom"]),int(w["start0"]),int(w["end"])) for w in windows[:s["top_windows"]]]
                if s["top_windows"] and not top:
                    records[name]=dict(name=name,status="awaiting_baseline_loadings",reason="Adaptive exclusion requires actual baseline fit")
                    continue
                records[name]=execute_scenario(args,name,s,source,out,runtime,top)
                if name=="baseline" and records[name]["status"]=="complete":
                    windows=loading_tables(Path(records[name]["folder"]),out/"drivers",args.window_bp)
                save_json(manifest_path,dict(version=VERSION,status="running" if args.command=="study" else "planned",
                    metadata=str(source["metadata_path"]),parameters=vars(args),plink_version=version,scenarios=records))
            if args.command=="study":
                run_drivers(args,out,records,source["metadata"])
                run_stability(args,out,records,source["metadata"])
            metadata=source["metadata"]
            saved=dict(version=VERSION,metadata=str(source["metadata_path"]),parameters=vars(args),
                       inputs=source["snapshots"],plink_version=version,scenarios=records,limitations=LIMITATIONS)
        if args.command!="study-plan" and not args.no_plots:figures_and_report(out,records,metadata)
        failed=any(r.get("status")=="failed" for r in records.values())
        saved.update(status="planned" if args.command=="study-plan" else "complete_with_failed_scenarios" if failed else "complete",
                     last_action=args.command,last_analysis_parameters=vars(args),scenarios=records)
        save_json(manifest_path,saved)
        print(json.dumps(dict(status=saved["status"],out=str(out),scenarios={k:v.get("status") for k,v in records.items()})))
        return 2 if failed else 0


if __name__=="__main__":
    try: raise SystemExit(main())
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print("pca-study: "+str(exc),file=sys.stderr);raise SystemExit(2)

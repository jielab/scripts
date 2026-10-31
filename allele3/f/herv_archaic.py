#!/usr/bin/env python3
"""New GU/IBDmix cohorts against the completed, read-only gnomAD background.

Each cohort retains its own analyzed-person denominator. Query the union of its
lifted tracts with the original SNV QC and ArchitectureAccumulator; never infer
the eligible denominator from the cached multi-ALT-only table.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import fnmatch
import gzip
import itertools
import json
import math
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from types import SimpleNamespace

from analysis_integrity import identity, transactional, validate_manifest

SCHEMA = "allele3-archaic-2"
AUTOSOMES = tuple(str(i) for i in range(1, 23))
CHROMOSOMES = (*AUTOSOMES, "X", "Y")
CODE = ("herv_archaic.py", "herv_modern.py", "scan.py", "scan_store.py", "0.common.py", "analysis_integrity.py")


def sample_ids(path):
    result = set()
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            value = line.split()[0].lstrip("#")
            if value.lower() in {"sample", "sample_id", "iid", "id", "fid"}:
                continue
            if value in result:
                raise ValueError(f"Duplicate analyzed sample in {path}: {value}")
            result.add(value)
    if not result:
        raise ValueError(f"Empty analyzed sample panel: {path}")
    return result


def metadata(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        parts = line.replace("\\t", "\t").split("\t", 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1]
    return result


def archive_complete(path):
    """Older GU publications keep .complete in the leading archive metadata.

    Inspect bounded headers only: never decompress a large tract member during
    discovery. GU's publisher puts completion/provenance before final tables.
    """
    with tarfile.open(path, 'r|gz') as archive:
        consumed = 0
        for member in archive:
            if member.name.removeprefix('./') == '.complete' and member.isfile():
                return True
            consumed += member.size + 512
            if consumed > 1024*1024:
                break
    return False


def inspect_dataset(path, module="auto", chromosomes=AUTOSOMES):
    """Metadata and file inventory only; do not open raw tracts at preflight."""
    path = Path(path).resolve()
    modules = [p for p in (path/"ibdmix", path/"ibdmix-GPU") if p.is_dir()]
    if module != "auto":
        modules = [path/module] if (path/module).is_dir() else []
    if len(modules) != 1:
        raise ValueError(f"{path.name}: expected one IBDmix module; found {[p.name for p in modules]}; select --archaic-module")
    source = modules[0]
    inputs, files, panels = [], [], []
    cohort, contract = None, None
    for chrom in chromosomes:
        directory = source/f"chr{chrom}"
        native = sorted((directory/"final").glob("all_archaic_refs*.segments.tsv.gz"))
        if len(native) > 1:
            raise ValueError(f"Ambiguous finalized tract versions: {directory}")
        tract = native[0] if native else directory/"ibdmix.raw.tar.gz"
        panel, meta = directory/"samples/ALL.txt", directory/"run.meta.tsv"
        required = [tract, panel, meta]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise ValueError("Missing finalized chromosome inputs: " + ", ".join(missing))
        # Published GU runs retain a publication marker after unpacked .complete
        # moves into their archive. Live workspaces instead retain .complete.
        markers = [p for p in (directory/".published-content", directory/".complete") if p.is_file()]
        if not markers and not (tract.name == 'ibdmix.raw.tar.gz' and archive_complete(tract)):
            raise ValueError(f"No GU completion/publication marker: {directory}")
        data = metadata(meta)
        build = data.get("genome_build", "").removeprefix("GRCh").removeprefix("b")
        if build not in {"37", "38"} or data.get("coordinate_system") != "0-based-half-open":
            raise ValueError(f"Unknown genome build or coordinate convention: {meta}")
        refs = tuple(sorted(data.get("refs", "").split()))
        if not refs:
            raise ValueError(f"Missing archaic reference list: {meta}")
        current = (build, refs, float(data["lod_cut"]), int(data["len_cut"]))
        if contract is not None and current != contract:
            raise ValueError(f"Build/reference/filter contract differs between chromosomes: {source}")
        contract = current
        people = sample_ids(panel)
        if cohort is not None and people != cohort:
            raise ValueError(f"Analyzed sample panels differ between chromosomes: {source}")
        cohort = people
        inputs.extend(required + markers)
        if (directory/"cache.meta.tsv").is_file():
            inputs.append(directory/"cache.meta.tsv")
        files.append(str(tract)); panels.append(str(panel))
    return dict(dataset=path.name, source=str(source), status="ready", build=contract[0],
                references=list(contract[1]), native_min_lod=contract[2], native_min_length=contract[3],
                n_samples=len(cohort), sample_panel=panels[0], chromosomes=list(chromosomes),
                files=files, inputs=[str(p) for p in inputs])


def chromosome_groups(selection="auto"):
    """Keep independent sex-chromosome panels and backgrounds out of autosomal LOCO."""
    tokens = [x.strip() for x in selection.split(",")]
    if selection == "auto":
        return [("autosomes", AUTOSOMES), ("X", ("X",)), ("Y", ("Y",))]
    if not tokens or any(x not in {"autosomes", "X", "Y"} for x in tokens) or len(set(tokens)) != len(tokens):
        raise ValueError("Chromosomes must be auto or a comma-separated selection of autosomes,X,Y")
    return [(x, AUTOSOMES if x == "autosomes" else (x,)) for x in tokens]


def discover(root, datasets="auto", module="auto", chromosomes="auto"):
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"GU source directory is missing: {root}")
    names = sorted(p.name for p in root.iterdir() if p.is_dir() and
                   (p.name == "1kg" or re.fullmatch(r"ukb\d+", p.name)))
    auto = datasets == "auto"
    if not auto:
        requested = []
        for pattern in datasets.split(","):
            matches = [n for n in names if fnmatch.fnmatchcase(n, pattern.strip())]
            if not matches:
                raise ValueError(f"No GU dataset matches {pattern!r} under {root}")
            requested.extend(matches)
        names = sorted(set(requested))
    records = []
    for name in names:
        for scope, selected in chromosome_groups(chromosomes):
            identifier = name if scope == "autosomes" else f"{name}.chr{scope}"
            try:
                record = inspect_dataset(root/name, module, selected)
                record.update(dataset=identifier, cohort=name, chromosome_scope=scope)
                records.append(record)
            except (OSError, ValueError, KeyError, tarfile.TarError, EOFError) as exc:
                if not auto and not (chromosomes == "auto" and scope in {"X", "Y"}):
                    raise ValueError(f"Requested dataset {identifier} is not ready: {exc}") from exc
                records.append(dict(dataset=identifier, chromosome_scope=scope, status="not_available", reason=str(exc), inputs=[]))
    if not any(r["status"] == "ready" for r in records):
        raise ValueError(f"No complete selected GU datasets under {root}; " +
                         "; ".join(f"{r['dataset']}: {r['reason']}" for r in records))
    return records


def modern_inputs(modern):
    # Reuse the existing regional-input audit, including relocated FASTA,
    # annotation, shard completion records, raw VCF and its index.
    from herv_pipeline import modern_files
    return modern_files(modern, "windows")


def file_identities(paths):
    return [identity(p) for p in dict.fromkeys(str(Path(p).resolve()) for p in paths)]


def source_identity(dataset, args, common):
    if args.min_lod < dataset["native_min_lod"] or args.min_length < dataset["native_min_length"]:
        raise ValueError(f"{dataset['dataset']}: cannot relax filters below the finalized GU calls")
    extra = [args.chain] if dataset["build"] != "38" else []
    executables = [shutil.which(args.bcftools), shutil.which(args.liftover)]
    return dict(schema=SCHEMA, dataset=dataset,
                settings={k: getattr(args, k) for k in ("min_lod", "min_length", "carrier_thresholds", "bootstrap", "seed", "min_training_sites")},
                inputs=file_identities([*dataset["inputs"], *common, *extra,
                                        *[p for p in executables if p],
                                        *[Path(__file__).with_name(p) for p in CODE]]))


def scan_tract_chrom(vcf, chrom, annotation, track, dataset, out, options, bcftools):
    """Original per-locus QC and exposure hook, restricted to the tract union."""
    import scan
    from herv_modern import read_indexed_header, query_regions
    from allele_common import FastaContext, IntervalMask, canonical_chrom
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    header = read_indexed_header(vcf, bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    tracks = scan.read_label_bed(track, chromosome=chrom)
    intervals = scan.merge_intervals([(x.start, x.end) for x in tracks])
    raw = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
    if len(raw) != 1:
        raise ValueError(f"No unique VCF contig for chr{chrom}")
    bed = out/f"chr{chrom}.regions.bed"
    with bed.open("w") as handle:
        for lo, hi in intervals:
            handle.write(f"{raw[0]}\t{lo}\t{hi}\n")
    masks = {name: IntervalMask(options[key], [chrom]) for name, key in
             (("lcr", "architecture_lcr_bed"), ("segdup", "architecture_segdup_bed"),
              ("low_mappability", "architecture_low_mappability_bed")) if options.get(key)}
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    blocks = []
    if chrom in AUTOSOMES:
        blocks, _ = scan.read_blocks(Path(options["ld_block_dir"])/f"EUR.{options['architecture_build']}.bed",
                                    "EUR", options["architecture_build"])
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    prefix = out/f"chr{chrom}"
    counter = scan.ArchitectureAccumulator({"EUR": blocks}, elements,
        thresholds=options["thresholds"], out_prefix=prefix, masks=masks,
        flag_annotations={"lcr", "segdup"} & set(header["info"]),
        tracks={dataset: tracks}, track_only=True)
    qc = Counter()
    try:
        with FastaContext(options["fasta"]) as fa:
            if any(hi > fa.lengths[chrom] for _, hi in intervals):
                raise ValueError(f"Lifted tract outside reference chromosome {chrom}")
            records = query_regions(vcf, fields, bed, bcftools) if intervals else iter(())
            for (c, pos), group in itertools.groupby(records, key=lambda r: (canonical_chrom(r["chrom"]), r["pos"])):
                if c != chrom:
                    raise ValueError("Unexpected chromosome in targeted archaic query")
                ref, alts, _, flags = scan.group_alleles(group, header, [], options, qc, None, [])
                if ref is None or not alts or (callable_mask and not callable_mask.contains(c, pos)):
                    continue
                observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                                  key=lambda a: (-a["af"], -a["ac"], a["alt"]))
                if not observed:
                    continue
                ctx = fa.context(c, pos).upper()
                if len(ctx) != 3 or ctx[1] != ref:
                    raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
                context, _ = scan.canonical_context(ctx, ref)
                primary = observed[0]
                counter.observe(c, pos, len(observed), context, primary["af"], primary["an"],
                    threshold_counts={t: sum(a["ac"] >= t for a in observed) for t in options["thresholds"]},
                    region_flags=flags)
                qc["eligible_sites"] += 1
            counter.write(prefix)
    finally:
        counter.close()
    audit = dict(chrom=chrom, queried_bp=sum(b-a for a, b in intervals), intervals=len(intervals), qc=dict(qc))
    (out/f"chr{chrom}.audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    return Path(str(prefix)+".exposure_strata.tsv.gz"), audit


@transactional("manifest.json")
def analyze_dataset(args):
    import scan
    from herv_modern import architecture_paths
    out = Path(args.out)
    dataset = args.dataset
    track_args = SimpleNamespace(dir=dataset["source"], format="ibdmix", out=str(out/"tracks"),
        sample_panel=dataset["sample_panel"], n_samples=None, build=dataset["build"], target_build="38",
        chain=args.chain, liftover=args.liftover, min_lod=args.min_lod, min_length=args.min_length,
        carrier_thresholds=args.carrier_thresholds, coordinates="bed0",
        files=dataset.get("files"), chromosomes=dataset.get("chromosomes"))
    track_result = scan.build_tracks(track_args)
    files, audits = [], []
    for chrom, vcf in args.tasks:
        print(f"[archaic {dataset['dataset']}] indexed chr{chrom}", flush=True)
        path, audit = scan_tract_chrom(vcf, chrom, Path(args.annotation)/f"chr{chrom}.tsv.gz",
            track_result["output"], dataset["dataset"], out/"counts", args.options, args.bcftools)
        files.append(path); audits.append(audit)
    scope = dataset.get("chromosome_scope", "autosomes")
    # Keep autosomal LOCO stable, and use a separate within-chromosome
    # descriptive background for X/Y. Permanent counters remain read-only.
    with tempfile.TemporaryDirectory(prefix="archaic-autosomes-", dir=out) as temporary:
        background = Path(temporary)/"global.tsv.gz"
        with gzip.open(background, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, scan.HEADERS["global"], delimiter="\t")
            writer.writeheader()
            for source in architecture_paths(args.modern_dir, "global"):
                with gzip.open(source, "rt") as rows:
                    for row in csv.DictReader(rows, delimiter="\t"):
                        if scan.canonical_chrom(row["chrom"]) in (AUTOSOMES if scope == "autosomes" else (scope,)):
                            writer.writerow(row)
        if scope == "autosomes":
            result = scan.summarize_exposure(background, files, out,
                reps=args.bootstrap, seed=args.seed, min_training_sites=args.min_training_sites)
        else:
            result = summarize_sex_exposure(background, files, out, scope, args.min_training_sites)
    if file_identities([r["path"] for r in args.source_identity["inputs"]]) != args.source_identity["inputs"]:
        raise ValueError("Archaic inputs changed during analysis")
    result.update(schema=SCHEMA, source_identity=args.source_identity, dataset=dataset["dataset"],
        n_samples=dataset["n_samples"], track=track_result, chromosome_audits=audits, chromosome_scope=scope,
        background_signature=args.background_signature, background_read_only=True,
        denominator="Original eligible SNV positions within lifted tracts; includes eligible single-ALT loci",
        interpretation="Unphased tract carrier prevalence, not archaic allele frequency or individual ancestry")
    return result


def summarize_sex_exposure(background, files, out, chrom, min_training_sites=1):
    """Descriptive X/Y O/E against matched loci outside each label on that chromosome.

    AC/AN come from the original scanner, never from an assumed diploid 2*N.
    Global counters cannot resolve PAR or sex strata, so no inferential P/CI is
    produced. Include outside_partition: autosomal LD coverage is not required.
    """
    dims = ("context", "af_bin", "an_bin", "sensitivity", "population_mode", "threshold")
    total = defaultdict(lambda: [0, 0])
    exposed = defaultdict(lambda: [0, 0])
    def read(paths):
        for path in paths:
            with gzip.open(path, "rt") as handle:
                yield from csv.DictReader(handle, delimiter="\t")
    def counts(row):
        n, m = int(row["n_sites"]), int(row["n_multi"])
        if not 0 <= m <= n:
            raise ValueError("Invalid eligible/multi-ALT counts")
        return n, m
    for row in read([background]):
        if row["threshold"] not in {"5", "10"} or row["population_mode"] != "all":
            continue
        if row["chrom"].removeprefix("chr") != chrom:
            raise ValueError("Sex-chromosome background contains another chromosome")
        n, m = counts(row)
        key = tuple(row[x] for x in dims)
        total[key][0] += n; total[key][1] += m
    if not total:
        raise ValueError(f"No original chr{chrom} global background counts")
    for row in read(files):
        if row["threshold"] not in {"5", "10"} or row["population_mode"] != "all":
            continue
        if row["chrom"].removeprefix("chr") != chrom:
            raise ValueError("Exposure chromosome differs from selected background")
        n, m = counts(row)
        key = (row["label"], tuple(row[x] for x in dims))
        exposed[key][0] += n; exposed[key][1] += m
    summaries = defaultdict(lambda: [0, 0, 0., 0])
    out = Path(out)
    with gzip.open(out/"sex_chromosome_strata.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "label", *dims, "exposed_sites", "observed_multi", "control_sites", "control_multi", "expected_multi"])
        for (label, st), (n, m) in sorted(exposed.items()):
            tn, tm = total[st][0] - n, total[st][1] - m
            if not 0 <= tm <= tn:
                raise ValueError(f"chr{chrom} exposure is not a subset of its original background: {label}")
            expected = n * (tm + .5)/(tn + 1) if tn >= min_training_sites else None
            writer.writerow([chrom, label, *st, n, m, tn, tm, expected if expected is not None else "NA"])
            z = summaries[label, st[-1], st[3]]
            z[0] += n; z[1] += m
            if expected is None: z[3] += n
            else: z[2] += expected
    fields = ['label','threshold','sensitivity','observed_multi','expected_multi','oe','ci_low','ci_high','bootstrap_p','fdr','block_count','n_eligible','n_without_loco','outside_partition_sites','status']
    with (out/"enrichment_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader()
        for (label, threshold, sensitivity), (n, m, e, missing) in sorted(summaries.items()):
            writer.writerow(dict(label=label, threshold=threshold, sensitivity=sensitivity,
                observed_multi=m, expected_multi="NA" if missing else e,
                oe=m/e if e and not missing else "NA", ci_low="NA", ci_high="NA",
                bootstrap_p="NA", fdr="NA", block_count=0, n_eligible=n,
                n_without_loco="NA", outside_partition_sites="NA",
                status="insufficient_same_chromosome_control" if missing else "descriptive_same_chromosome"))
    manifest = dict(status="complete", n_tests=0, n_descriptive_rows=len(summaries),
        chromosome_scope=chrom, primary_thresholds=[5, 10],
        inference="Descriptive only; no P, CI or FDR; no autosomal LD-block assumption",
        background="Same chromosome outside each label, matched context/primary-AF/AN/sensitivity; (multi+0.5)/(sites+1)",
        empty_exposure=not bool(summaries),
        limitations=["Original AC/AN retained; no assumption AN=2*N",
            "PAR/non-PAR and sex composition are not resolved by existing counters",
            "Carrier denominator is the chromosome-specific GU analyzed-person panel, including zero-tract people; not allele frequency",
            "X/Y estimates must not be interpreted as autosomal LOCO tests"])
    (out/"enrichment_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def reusable(directory, expected):
    try:
        result = validate_manifest(directory, require_evidence=True)
        return result if result.get("status") == "complete" and result.get("source_identity") == expected else None
    except (OSError, ValueError, KeyError):
        return None


def combine_summaries(sources, out):
    """Recompute multiplicity over the selected cohorts, keeping cohort O/E."""
    import scan
    import scan
    data, fields = [], None
    for dataset, path in sources:
        with Path(path).open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fields is not None and fields != reader.fieldnames:
                raise ValueError("Archaic cohort summary schemas differ")
            fields = reader.fieldnames
            for row in reader:
                row["dataset"] = dataset
                row["fdr"] = "NA"
                data.append(row)
    groups = {}
    for row in data:
        try:
            p = float(row["bootstrap_p"])
        except (ValueError, TypeError):
            continue
        if row["status"] == "complete" and math.isfinite(p):
            groups.setdefault((row["threshold"], row["sensitivity"]), []).append(row)
    for group in groups.values():
        for row, q in zip(group, scan.bh_fdr([float(r["bootstrap_p"]) for r in group])):
            row["fdr"] = q
    with (Path(out)/"archaic_enrichment.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, ["dataset", *(fields or [])], delimiter="\t")
        writer.writeheader(); writer.writerows(data)
    return data


def plot_results(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    selected = [r for r in rows if r["status"] == "complete" and r["threshold"] == "10"
                and r["sensitivity"] == "exclude_combined" and r["label"].endswith(":any")
                and r["ci_low"] not in {"", "NA", "None"}]
    if not selected:
        return
    fig, ax = plt.subplots(figsize=(10, max(4, .24*len(selected))))
    for i, row in enumerate(selected):
        lo, value, hi = (float(row[k]) for k in ("ci_low", "oe", "ci_high"))
        ax.plot([lo, hi], [i, i], color="#176B87"); ax.plot(value, i, "o", color="#176B87", ms=3)
    ax.set_yticks(range(len(selected)), [r["label"] for r in selected], fontsize=7)
    ax.axvline(1, color="black", lw=1, ls="--")
    ax.set_xlabel("Observed / LOCO expected; AC ≥ 10, strict masks; 95% LD-block bootstrap CI")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"archaic_enrichment.{ext}", dpi=160)
    plt.close(fig)


@transactional("manifest.json")
def run(args):
    from herv_modern import targeted_inputs
    selection = getattr(args, "chromosomes", "auto")
    records = discover(args.gu_root, args.datasets, args.module, selection)
    common = modern_inputs(args.modern_dir)
    # This validates the relocated source and all original QC inputs. New
    # exposure queries use the original options, not fresh scanner defaults.
    query_args = SimpleNamespace(out=args.out, modern_dir=args.modern_dir, fasta=None,
        erv_joined=None, segdup_bed=None, lcr_bed=None, low_mappability_bed=None, chromosomes=None)
    original, annotation, tasks, chromosomes = targeted_inputs(query_args, allowed_chromosomes=CHROMOSOMES)
    options = original["options"]
    if not {5, 10} <= set(options.get("thresholds", [])):
        raise ValueError("Archaic analysis requires original AC5/10 counters")
    ready = []
    for record in records:
        if record["status"] != "ready":
            continue
        missing = set(record["chromosomes"]) - chromosomes
        if missing:
            reason = f"Original modern scan lacks chromosomes: {sorted(missing)}"
            if selection != "auto" or record["chromosome_scope"] == "autosomes":
                raise ValueError(reason)
            record.update(status="not_available", reason=reason)
        else:
            ready.append(record)
    if not ready:
        raise ValueError("No selected GU/modern chromosome pairs are ready")
    if any(r["chromosome_scope"] == "autosomes" for r in ready) and not options.get("ld_block_dir"):
        raise ValueError("Autosomal LOCO requires EUR LD blocks")
    summaries, manifests, snapshots = [], {}, []
    for dataset in ready:
        expected = source_identity(dataset, args, common)
        snapshots.append(expected['inputs'])
        directory = Path(args.cache_dir)/dataset["dataset"]
        result = reusable(directory, expected)
        if result is None:
            # cache_dir is this child's output container, not a scientific input.
            # Keep actual sources in the transaction's protected input set.
            child = SimpleNamespace(**{k: v for k, v in vars(args).items() if k != 'cache_dir'})
            child.out = str(directory); child.dataset = dataset; child.source_identity = expected
            child.tasks = [(c, path) for c, path in tasks if c in dataset["chromosomes"]]
            child.annotation = str(annotation); child.options = options
            child.background_signature = original["signature"]
            analyze_dataset(child)
            result = validate_manifest(directory, require_evidence=True)
            dataset["operation"] = "computed"
        else:
            print(f"[reuse] archaic {dataset['dataset']}: inputs and outputs verified", flush=True)
            dataset["operation"] = "reused"
        manifests[dataset["dataset"]] = str(directory/"manifest.json")
        summaries.append((dataset["dataset"], directory/"enrichment_summary.tsv"))
    for saved in snapshots:
        if file_identities([r['path'] for r in saved]) != saved:
            raise ValueError('Archaic inputs changed while preparing the combined result')
    rows = combine_summaries(summaries, args.out)
    with (Path(args.out)/"datasets.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, ["dataset", "chromosome_scope", "status", "operation", "source", "n_samples", "build", "reason"],
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(records)
    plot_results(rows, args.out)
    return dict(status="complete", schema=SCHEMA, datasets=records, dataset_manifests=manifests,
        n_tests=sum(r["status"] == "complete" for r in rows), n_summary_rows=len(rows), background_signature=original["signature"], background_read_only=True,
        fdr_scope="All selected archaic datasets within each AC threshold and sensitivity",
        method="Original scanner QC; autosomes: LOCO/EUR LD-block bootstrap; X/Y: separate same-chromosome outside-tract descriptive O/E",
        interpretation="Cohorts remain separate; overlapping people are not independent replication; zero-tract people remain in carrier denominator")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gu-root", default="/mnt/f/analysis/gu")
    p.add_argument("--datasets", default="auto")
    p.add_argument("--module", choices=("auto", "ibdmix", "ibdmix-GPU"), default="auto")
    p.add_argument("--chromosomes", default="auto", help="auto: autosomes plus available X/Y; or autosomes,X,Y / X,Y / X / Y")
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--chain", default="/mnt/f/gen/chains/hg19ToHg38.over.chain.gz")
    p.add_argument("--liftover", default="/mnt/d/software/bin/liftOver")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--min-lod", type=float, default=4.)
    p.add_argument("--min-length", type=int, default=50000)
    p.add_argument("--carrier-thresholds", default="0.01,0.05,0.10,0.20")
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261024)
    p.add_argument("--min-training-sites", type=int, default=1)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.bootstrap < 2000 or args.min_length < 1 or args.min_lod < 0 or args.min_training_sites < 1:
        raise ValueError("Require ≥2000 bootstrap draws, positive length/training minimum and nonnegative LOD")
    thresholds = [float(x) for x in args.carrier_thresholds.split(",")]
    if not thresholds or any(not 0 < x <= 1 for x in thresholds) or len(set(thresholds)) != len(thresholds):
        raise ValueError("Carrier thresholds must be distinct values in (0,1]")
    protected = [Path(args.gu_root).resolve(), Path(args.modern_dir).resolve()]
    outputs = [Path(args.out).resolve(), Path(args.cache_dir).resolve()]
    for out in outputs:
        for source in protected:
            if out == source or out in source.parents or source in out.parents:
                raise ValueError(f"Generated output overlaps permanent input: {out}")
    if outputs[0] == outputs[1] or outputs[0] in outputs[1].parents or outputs[1] in outputs[0].parents:
        raise ValueError("Summary output and dataset cache must be separate, non-nested directories")
    # Different report directories may share this persistent cohort cache.
    # Keep a stable view while validating, computing and combining its files.
    outputs[1].mkdir(parents=True, exist_ok=True)
    with (outputs[1]/'.archaic.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Archaic dataset cache is already in use: {outputs[1]}') from exc
        run(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"archaic: {exc}")

#!/usr/bin/env python3
"""Portable local entry point. No AlphaGenome API, key or network dependency."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import sys
import time

from core import VERSION, code_hash, dump_json, fingerprint, stable_hash, tsv_writer
from gnomad_io import (FastaContext, _backend, _chromosomes, _file_chrom, discover_inputs,
                       read_header, validate_count_header, parse_region)
from intervals import make_callable
from scan import discover_populations, scan_file
from summarize import summarize

DEFAULT_VCF = "/mnt/f/annot/gnomad/sites"
DEFAULT_FASTA = "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz"
DEFAULT_OUT = "/mnt/d/analysis/allele3/gnomad_local"


def parser():
    p = argparse.ArgumentParser(
        description="allele3 local: gnomAD SNV ALT occupancy, distribution and ancestry summaries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("command", nargs="?", default="run",
                   choices=["run", "scan", "summarize", "check", "annotate", "make-callable"])
    source = p.add_mutually_exclusive_group()
    source.add_argument("--vcf", help="One sites VCF/BCF; --region is optional")
    source.add_argument("--dir-vcf", "--vcf-dir", dest="vcf_dir",
                        help="Directory of disjoint chromosome VCF files")
    p.add_argument("--chr", dest="chromosomes", help="Directory selection: autosomes, all, or 1,2,X")
    p.add_argument("--region", help="1-based inclusive CHR:START-END; requires --vcf")
    p.add_argument("--fasta", default=os.environ.get("ALLELE_FASTA", DEFAULT_FASTA))
    p.add_argument("--dir-out", "--out", dest="out",
                   default=os.environ.get("ALLELE_LOCAL_OUT", os.environ.get("ALLELE_OUT", DEFAULT_OUT)))
    p.add_argument("--dir-cache", "--cache-dir", dest="cache_dir",
                   default=os.environ.get("ALLELE_CACHE_DIR", ""), help="Default: OUT/shards")
    p.add_argument("--jobs", type=int, default=4, help="Concurrent chromosome scans; I/O often limits speed")
    p.add_argument("--backend", choices=["auto", "bcftools", "python"], default="auto")
    p.add_argument("--window-size", type=int, default=1000000)
    p.add_argument("--thresholds", default="1,2,5,10", help="ALT count thresholds; 1 is always included")
    p.add_argument("--min-an", type=int, default=100, help="Per-ALT AN floor; not a universal callability guarantee")
    p.add_argument("--min-af", type=float, default=0.0, help="Per-ALT AF floor (AC/AN); default retains singletons")
    p.add_argument("--populations", default="auto", help="auto, none, or ancestry suffixes e.g afr,eas,nfe")
    p.add_argument("--min-pop-an", type=int, default=100)
    p.add_argument("--min-pop-af", type=float, default=0.0,
                   help="Ancestry-specific AF sensitivity; unequal detection power still remains")
    p.add_argument("--include-unfiltered", action="store_true", help="Explicitly also accept FILTER=.; default PASS only")
    p.add_argument("--exclude-flags", default="", help="Optional INFO region flags: lcr,segdup")
    p.add_argument("--vep", action="store_true", help="Also read existing INFO/vep or CSQ; increases VCF transfer/CPU")
    p.add_argument("--write-all-sites", action="store_true", help="Optional all-site compact TSV; default only multi-ALT details")
    p.add_argument("--candidate-storage", choices=["deferred", "eager"], default="deferred",
                   help="Defer the redundant three-candidate-per-site table until explicit full offline annotate; sparse Atlas does not need it")
    p.add_argument("--callable-bed", help="Optional same-build BED mask; applied to both numerator and denominator")
    p.add_argument("--build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--source-label", default="gnomAD-v4.1-genomes")
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--scores", help="Local precomputed SNV scores for annotate only")
    p.add_argument("--score-column", default="phred")
    p.add_argument("--score-format", choices=["tsv", "avi"], default="tsv",
                   help="avi explicitly selects the official six-column headerless format")
    p.add_argument("--score-label", default="local")
    p.add_argument("--score-version", default="unspecified")
    p.add_argument("--score-build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--score-backend", choices=["auto", "stream", "tabix"], default="auto")
    p.add_argument("--all-sites-an", help="Local official overall AN TSV for make-callable")
    p.add_argument("--callable-min-an", type=int, default=100000)
    p.add_argument("--output-bed", help="New plain BED output for make-callable")
    p.add_argument("--with-architecture", action="store_true")
    p.add_argument("--ld-block-dir", default="/mnt/f/LD/block")
    p.add_argument("--architecture-build", choices=["37", "38"], default="38")
    p.add_argument("--erv-joined", help="Joined all-repeat table; absent annotation is NA")
    p.add_argument("--architecture-prefix", default="architecture", help="Shard-local filename prefix")
    p.add_argument("--architecture-min-sites", type=int, default=10000)
    p.add_argument("--architecture-top-fraction", type=float, default=.01)
    p.add_argument("--architecture-annotation-cache")
    p.add_argument("--architecture-segdup-bed")
    p.add_argument("--architecture-lcr-bed")
    p.add_argument("--architecture-low-mappability-bed")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def validated_options(a):
    if a.jobs < 1 or a.window_size < 1 or a.min_an < 1 or a.min_pop_an < 1 or a.min_training_sites < 1:
        raise ValueError("jobs, window size, AN floors, and minimum training sites must be positive")
    if not 0 <= a.min_af <= 1 or not 0 <= a.min_pop_af <= 1:
        raise ValueError("AF thresholds must be between 0 and 1")
    try:
        thresholds = sorted({1, *(int(t) for t in a.thresholds.split(","))})
    except ValueError as e:
        raise ValueError("--thresholds requires positive comma-separated integers") from e
    if thresholds[0] < 1:
        raise ValueError("AC thresholds must be positive")
    flags = [x for x in a.exclude_flags.split(",") if x]
    if set(flags) - {"lcr", "segdup"}:
        raise ValueError("--exclude-flags accepts lcr and segdup")
    if a.region:
        if not a.vcf:
            raise ValueError("--region requires one --vcf")
        parse_region(a.region)
    if a.vcf and a.chromosomes:
        raise ValueError("Use --region with one --vcf; --chr selects directory inputs")
    options = dict(fasta=str(Path(a.fasta).expanduser().resolve()), backend=a.backend,
                window_size=a.window_size, thresholds=thresholds, min_an=a.min_an,
                min_af=a.min_af, populations=a.populations, min_pop_an=a.min_pop_an,
                min_pop_af=a.min_pop_af, include_unfiltered=a.include_unfiltered,
                exclude_flags=flags, vep=a.vep, write_all_sites=a.write_all_sites,
                candidate_storage=a.candidate_storage,
                callable_bed=str(Path(a.callable_bed).resolve()) if a.callable_bed else None,
                build=a.build, region=a.region, source_label=a.source_label)
    if a.with_architecture:
        if a.build != "GRCh"+a.architecture_build:
            raise ValueError("Architecture annotation build must match the VCF/reference build")
        if a.erv_joined and a.architecture_build != "38":
            raise ValueError("Current joined RepeatMasker input is GRCh38; use an explicitly prepared matching-build track")
        if not a.architecture_prefix or Path(a.architecture_prefix).name!=a.architecture_prefix or a.architecture_prefix in {".",".."}:
            raise ValueError("--architecture-prefix must be a shard-local filename, not a path")
        if a.architecture_min_sites<1 or not 0<a.architecture_top_fraction<=1:
            raise ValueError("Invalid architecture ranking thresholds")
        options.update(with_architecture=True, architecture_schema=2,
            architecture_build=a.architecture_build, ld_block_dir=str(Path(a.ld_block_dir).resolve()),
            erv_joined=str(Path(a.erv_joined).resolve()) if a.erv_joined else None,
            architecture_prefix=a.architecture_prefix, architecture_min_sites=a.architecture_min_sites,
            architecture_top_fraction=a.architecture_top_fraction)
        for key in ('architecture_annotation_cache','architecture_segdup_bed','architecture_lcr_bed','architecture_low_mappability_bed'):
            value=getattr(a,key)
            options[key]=str(Path(value).resolve()) if value else None
        from architecture import annotation_identity
        options['architecture_annotations']=annotation_identity(options)
    return options


def load_run(out):
    path = Path(out) / "run.json"
    if not path.is_file():
        raise ValueError(f"No completed/scanned run manifest: {path}; run scan first")
    data = json.loads(path.read_text())
    if data.get("status") not in ("scanned", "complete"):
        raise ValueError(f"Scan is not complete at {out}; resume the original run command first")
    for shard in data["shards"]:
        verify_shard(shard, data["signature"])
    return data


def verify_shard(shard, signature):
    path = Path(shard)
    marker = path / "COMPLETE.json"
    if not marker.is_file():
        raise ValueError(f"Missing shard completion marker: {marker}")
    complete = json.loads(marker.read_text())
    if complete["signature"] != signature:
        raise ValueError(f"Cache identity differs at {path}; use a separate output/cache directory")
    for name, metadata in complete["outputs"].items():
        f = path / name
        if not f.is_file() or f.stat().st_size != metadata["size"] or fingerprint(f)["edge_sha256"] != metadata["edge_sha256"]:
            raise ValueError(f"Completed cache file is missing or changed: {f}; restore it or use a new cache")
    return complete


def _scan_task(task):
    """Serialize identical shared-cache shards across distinct output folders."""
    target = Path(task[1])
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.with_name(target.name + ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _scan_task_locked(task)


def _scan_task_locked(task):
    source, target, options, expected, signature = task
    target = Path(target)
    if (target / "COMPLETE.json").exists():
        verify_shard(target, signature)
        print(f"[resume] {Path(source).name}: completed shard reused", flush=True)
        return str(target)
    if target.exists():
        raise ValueError(f"Unrecognized incomplete directory: {target}; use a new --dir-cache")
    staging = target.with_name(target.name + ".partial")
    # This precisely named staging directory is owned by this workflow.
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    try:
        result = scan_file(source, staging, options, expected)
        outputs = {p.name: fingerprint(p) for p in staging.iterdir() if p.is_file()}
        dump_json(staging / "COMPLETE.json", dict(signature=signature, outputs=outputs,
                                                 qc=result["qc"], elapsed_seconds=result["elapsed_seconds"]))
        os.replace(staging, target)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    print(f"[done] {Path(source).name}: {result['qc']['eligible_sites']:,} eligible loci; "
          f"{result['qc']['multi_alt_sites']:,} multi-ALT", flush=True)
    return str(target)


def merge_small_reports(shards, out):
    for name in ("frequency_counts.tsv", "population_counts.tsv", "population_classes.tsv",
                 "annotation_counts.tsv", "relationship_counts.tsv"):
        counts, fields = Counter(), None
        for shard in shards:
            with (Path(shard) / name).open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if fields is None:
                    fields = reader.fieldnames
                elif fields != reader.fieldnames:
                    raise ValueError(f"Inconsistent {name} schemas")
                for row in reader:
                    counts[tuple(row[x] for x in fields[:-1])] += int(row["count"])
        handle, writer = tsv_writer(Path(out) / name, fields)
        with handle:
            for key, value in sorted(counts.items()):
                writer.writerow(dict(zip(fields, (*key, value))))
    pop = defaultdict(Counter)
    with (Path(out) / "population_counts.tsv").open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pop[(r["chrom"], r["threshold"], r["population"])][r["metric"]] += int(r["count"])
    cols = ["chrom", "threshold", "population", "n_evaluable", "n_sites", "n_multi", "n_three", "multi_alt_fraction"]
    handle, writer = tsv_writer(Path(out) / "population_summary.tsv", cols)
    with handle:
        for key, values in sorted(pop.items()):
            row = dict(zip(cols[:3], key))
            row.update({m: values[m] for m in cols[3:7]})
            row["multi_alt_fraction"] = values["n_multi"]/values["n_sites"] if values["n_sites"] else "."
            writer.writerow(row)
    from spectrum import merge_site_spectra
    return merge_site_spectra(shards, out)


def summarize_run(out, min_training_sites=100, make_plots=True):
    out = Path(out)
    run = load_run(out)
    options = run["options"]
    mask = options.get("callable_bed")
    if mask and fingerprint(mask) != run["identity"]["callable_bed"]:
        raise ValueError("Callable mask changed since scanning; use the original file or rescan to a new output")
    stage = out / "summary.partial"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    try:
        result = summarize(run["shards"], stage, window_size=options["window_size"],
                           min_training_sites=min_training_sites, make_plots=make_plots,
                           callable_bed=mask, region=options.get("region"), thresholds=options["thresholds"])
        spectrum = merge_small_reports(run["shards"], stage)
        if options.get('with_architecture'):
            from architecture import summarize_architecture, annotation_identity
            if annotation_identity(options)!=options.get('architecture_annotations'):
                raise ValueError('Architecture annotations changed since scanning')
            summarize_architecture(run['shards'],stage,options)
        qc = Counter()
        for shard in run["shards"]:
            record = json.loads((Path(shard) / "scan.json").read_text())
            qc.update(record["qc"])
        dump_json(stage / "qc_summary.json", dict(qc))
        report = f"""# allele3 local analysis

Version: {VERSION}

Source: {options['source_label']}; genome build: {options['build']}.

- Qualifying SNV positions: {qc['eligible_sites']:,}
- Positions with at least two qualifying ALT nucleotides: {qc['multi_alt_sites']:,}
- Positions with three qualifying ALT nucleotides: {qc['three_alt_loci']:,}
- Multi-ALT positions with different ALT AN values: {qc['multi_alt_loci_discordant_an']:,}
- AC thresholds: {options['thresholds']}
- Per-ALT AN floor: {options['min_an']}
- Callable mask: {mask or 'none; proportions are conditional on observed qualifying SNV positions'}
- Site AC spectrum coverage: {spectrum['status']} ({len(spectrum['measured_shards'])}/{spectrum['total_shards']} shards)
- Theoretical candidate storage: {options.get('candidate_storage', 'eager (legacy)')}

The primary quantity is observed non-reference nucleotide occupancy. REF is not
assumed ancestral, major, present, or reconstructable from split AN values.
ALT order is not a mutation timeline. Failed/missing ALT records are not used as
unlisted counterfactuals. Frequencies are each ALT's AC/AN.

The context background uses other scanned autosomes and reports support and
fallback levels. It estimates a conditional observation frequency, not mutation
rate or evolutionary retention. X/Y/MT and regional/single-chromosome scans are
descriptive. All ancestry comparisons remain sensitive to sample size, coverage,
ancestry composition and genotype filtering.

No AlphaGenome network requests were made by this pipeline. Optional offline AVI
annotation is exploratory: AVI used gnomAD frequency proxy labels in training.
This workflow cannot infer haplotype ancestry, allele age, pleiotropy, or balancing
selection from sites-only counts.

Detailed multi-ALT sites and original per-ALT counts are kept in the per-chromosome
shards. With deferred candidate storage, the redundant three-substitution table
is reconstructed losslessly only if full offline annotate is requested. Sparse
Atlas uses its separate bounded plan. Unlisted candidate AC/AN/AF remain missing.

site_ac_spectrum.tsv and singleton_sensitivity.tsv distinguish multi-ALT sites
containing AC=1 ALTs from sites still having two ALTs after AC>=2. The denominator
is the original eligible multi-ALT cohort. This is a count sensitivity, not loss
over evolutionary time. Old caches lacking this spectrum are explicitly unavailable;
mixed cache summaries are partial and cannot be treated as genome-wide zero values.
"""
        (stage / "REPORT.md").write_text(report, encoding="utf-8")
        target, backup = out / "summary", out / "summary.previous"
        # The summarizer writes absolute paths before its directory is committed.
        # Rebase metadata while it is still private, then publish the whole set.
        def rebase(value):
            if isinstance(value, str) and value.startswith(str(stage.resolve()) + os.sep):
                return str(target.resolve()) + value[len(str(stage.resolve())):]
            if isinstance(value, list):
                return [rebase(x) for x in value]
            if isinstance(value, dict):
                return {k: rebase(v) for k, v in value.items()}
            return value
        result = rebase(result)
        for metadata in stage.glob("*.json"):
            dump_json(metadata, rebase(json.loads(metadata.read_text())))
        if backup.exists():
            shutil.rmtree(backup)
        if target.exists():
            os.replace(target, backup)
        os.replace(stage, target)
        if backup.exists():
            shutil.rmtree(backup)
        run["status"] = "complete"
        run["summary"] = str(target.resolve())
        run["summary_parameters"] = dict(min_training_sites=min_training_sites, make_plots=make_plots)
        dump_json(out / "run.json", run)
    except BaseException:
        if stage.exists():
            shutil.rmtree(stage)
        raise
    print(f"[complete] {out / 'summary'}", flush=True)
    return result


def run_scan(a, options):
    out = Path(a.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    directory = None if a.vcf else (a.vcf_dir or os.environ.get("ALLELE_VCF_DIR", os.environ.get("ALLELE_GNOMAD", DEFAULT_VCF)))
    chromosomes = None if a.vcf else (a.chromosomes or "autosomes")
    paths = discover_inputs(a.vcf, directory, chromosomes)
    with FastaContext(options["fasta"]):
        pass
    if options["build"] != "GRCh38" and "gnomad-v4" in options["source_label"].lower():
        raise ValueError("gnomAD v4 uses GRCh38; supply a correct build/source label")
    identity = dict(version=VERSION, code_sha256=code_hash(), options=options,
                    inputs=[fingerprint(p) for p in paths], fasta=fingerprint(options["fasta"]),
                    fasta_index=fingerprint(options["fasta"] + ".fai"),
                    callable_bed=fingerprint(options["callable_bed"]) if options["callable_bed"] else None)
    if Path(options["fasta"] + ".gzi").exists():
        identity["fasta_gzi"] = fingerprint(options["fasta"] + ".gzi")
    signature = stable_hash(identity)
    manifest = out / "run.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if old["signature"] != signature:
            raise ValueError("Inputs, scan options or code differ from this output. Use a new --dir-out; old results are preserved.")
    elif any(out.iterdir()):
        unexpected = [p.name for p in out.iterdir() if p.name != ".allele3.lock"]
        if unexpected:
            raise ValueError(f"Output is not an empty local-workflow directory: {out}; choose a new --dir-out")
    cache = Path(a.cache_dir).expanduser().resolve() if a.cache_dir else out / "shards"
    cache.mkdir(parents=True, exist_ok=True)
    tasks, shards = [], []
    for i, path in enumerate(paths):
        expected = _file_chrom(path)
        if options.get("region"):
            expected = parse_region(options["region"])[0]
        label = f"chr{expected}" if expected else f"input{i+1}"
        target = cache / (label + "." + signature[:16])
        shards.append(str(target))
        tasks.append((str(path), str(target), options, expected, signature))
    run = dict(status="running", signature=signature, identity=identity, options=options,
               shards=shards, version=VERSION, cache=str(cache), inputs=[str(p.resolve()) for p in paths])
    dump_json(manifest, run)
    if options.get('with_architecture') and options.get('architecture_annotation_cache') and options.get('erv_joined'):
        from erv import prepare_joined_cache
        prepare_joined_cache(options['erv_joined'],options['architecture_annotation_cache'])
    started = time.monotonic()
    print(f"[scan] {len(paths)} input(s); {min(a.jobs, len(paths))} concurrent scan(s); "
          f"backend={_backend(paths[0], a.backend)}; no API calls", flush=True)
    if a.jobs == 1 or len(tasks) == 1:
        for task in tasks:
            _scan_task(task)
    else:
        with ProcessPoolExecutor(max_workers=min(a.jobs, len(tasks))) as pool:
            futures = [pool.submit(_scan_task, task) for task in tasks]
            for future in as_completed(futures):
                future.result()
    run["status"] = "scanned"
    run["scan_elapsed_seconds_this_invocation"] = time.monotonic() - started
    dump_json(manifest, run)
    return run


def main(argv=None):
    a = parser().parse_args(argv)
    if a.command == "make-callable":
        if not a.all_sites_an or not a.output_bed:
            raise ValueError("make-callable requires --all-sites-an and --output-bed")
        if str(a.output_bed).endswith((".gz", ".bgz")):
            raise ValueError("--output-bed must be a plain .bed path")
        data = make_callable(a.all_sites_an, a.output_bed, a.callable_min_an, _chromosomes(a.chromosomes))
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 0
    if a.command == "check":
        options = validated_options(a)
        directory = None if a.vcf else (a.vcf_dir or os.environ.get("ALLELE_VCF_DIR", DEFAULT_VCF))
        paths = discover_inputs(a.vcf, directory, None if a.vcf else (a.chromosomes or "autosomes"))
        with FastaContext(options["fasta"]) as fa:
            print(f"FASTA: {len(fa.lengths)} indexed contigs; build={a.build}")
        for path in paths:
            h = read_header(path, a.backend)
            validate_count_header(h)
            print(f"{path.name}: backend={_backend(path,a.backend)}; "
                  f"INFO AC/AN valid; samples_in_header={h['samples']}; "
                  f"ancestries={','.join(discover_populations(h,a.populations)) or 'none'}")
        print("Ready. No API key is read. Header sample count 0 is expected for sites-only files.")
        return 0
    out = Path(a.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".allele3.lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise ValueError(f"Another allele3 process is using {out}") from e
        if a.command in ("run", "scan"):
            run_scan(a, validated_options(a))
            if a.command == "run":
                summarize_run(out, a.min_training_sites, not a.no_plots)
        elif a.command == "summarize":
            summarize_run(out, a.min_training_sites, not a.no_plots)
        elif a.command == "annotate":
            if not a.scores:
                raise ValueError("annotate requires --scores pointing to an existing local file")
            run = load_run(out)
            if a.score_build != run["options"]["build"]:
                raise ValueError("Score genome build must match the scanned VCF/FASTA; no automatic liftover")
            from local_scores import annotate_scores, preflight_scores
            from candidate_cache import candidate_files_for_run
            preflight_scores(a.scores, a.score_column, a.score_build, a.score_backend, a.score_format)
            candidate_cache = candidate_files_for_run(run, out / "candidate_cache")
            result = annotate_scores(
                candidate_cache["files"],
                a.scores, out / "scores", a.score_column, a.score_build,
                a.score_label, a.score_version, a.score_backend, 1000000, a.score_format)
            print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, FileExistsError, OSError) as error:
        print(f"allele3: ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)

#!/usr/bin/env python3
"""HERV common-support summaries, body decomposition and held-out edge scans.

Multi-ALT occupancy is conditional on the original eligible SNV universe.
Reference-allele observation and longitudinal persistence are not inferred.
"""
from __future__ import annotations
from analysis_integrity import transactional
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
import scan
from allele_common import (FastaContext, IntervalMask, canonical_chrom,
                           dump_json, open_text, _parse_header, _record, af_bin, an_bin)

SCHEMA = "herv-focus-modern-1.0"
AUTOSOMES = {str(i) for i in range(1, 23)}


def rows(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_rows(path, data, fields=None):
    data = list(data)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    columns = fields or (list(data[0]) if data else ["status"])
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in data:
            writer.writerow({k: "NA" if v is None else v for k, v in row.items()})


def file_state(path):
    p = Path(path).resolve()
    s = p.stat()
    return {"path": str(p), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns}


def architecture_paths(modern, kind):
    combined = Path(modern) / "summary" / f"architecture.{kind}_strata.tsv.gz"
    if combined.is_file():
        return [combined]
    run = scan.load_run(modern)
    prefix = run["options"].get("architecture_prefix", "architecture")
    files = [Path(s) / f"{prefix}.{kind}_strata.tsv.gz" for s in run["shards"]]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise ValueError(f"Missing existing {kind} counters: {missing[:3]}; original caches are read-only")
    return files


class Background:
    """Fixed base denominator; validation uses only discovery chromosomes."""
    def __init__(self, paths, min_training_sites=100, discovery=None, include_all_chromosomes=False):
        self.counts = defaultdict(lambda: np.zeros(2))
        self.total = defaultdict(lambda: np.zeros(2))
        self.train = defaultdict(lambda: np.zeros(2))
        self.discovery = {canonical_chrom(c) for c in (discovery or [str(i) for i in range(1, 23, 2)])}
        if not self.discovery <= AUTOSOMES:
            raise ValueError("Discovery chromosomes must be autosomes")
        self.minimum = min_training_sites
        for path in paths:
            for r in rows(path):
                if r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}:
                    continue
                chrom = canonical_chrom(r["chrom"])
                if not include_all_chromosomes and chrom not in AUTOSOMES:
                    continue
                key = self.key(r)
                value = np.array([int(r["n_sites"]), int(r["n_multi"])])
                self.counts[(chrom, key)] += value
                self.total[key] += value
                if chrom in self.discovery:
                    self.train[key] += value

    @staticmethod
    def key(r):
        sensitivity = r["sensitivity"].replace("_no_cpg", "")
        if sensitivity == "strict":
            sensitivity = "exclude_combined"
        return (r["context"], r["af_bin"], r["an_bin"], sensitivity, str(r["threshold"]))

    def rate(self, chrom, key, split="all"):
        if split == "all":
            n, m = self.total[key] - self.counts[(chrom, key)]
        elif split == "discovery":
            n, m = self.train[key] - self.counts[(chrom, key)]
        else:
            n, m = self.train[key]
        return None if n < self.minimum else (m + .5) / (n + 1.)


class PrimaryBackground(Background):
    """Exact point-estimate definition in scan.summarize_exposure.

    Include every recorded global chromosome, train with at least one eligible
    site, and never convert an unsupported full-scope total into common support.
    The exposure audit is collected during the existing legacy-counter read.
    """
    def __init__(self, paths, families):
        super().__init__(paths, min_training_sites=1, include_all_chromosomes=True)
        self.families = set(families)
        self.labels = {"family:"+f for f in self.families}
        self.exposure = defaultdict(lambda: np.zeros(5))
        self.exposure_chromosomes = defaultdict(set)

    def observe_exposure(self, r):
        if (r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}
                or r["sensitivity"] not in {"all", "exclude_combined"}
                or r["label"] not in self.labels
                or r["block_id"] == "outside_partition"):
            return
        c = canonical_chrom(r["chrom"])
        sensitivity = "strict" if r["sensitivity"] == "exclude_combined" else "all"
        key = r["label"][7:], sensitivity, int(r["threshold"])
        n, m = int(r["n_sites"]), int(r["n_multi"])
        value = self.exposure[key]
        # Raw observed, partial expected, raw denominator, unsupported M/N.
        value[0] += m
        value[2] += n
        if n:
            self.exposure_chromosomes[key].add(c)
        rate = self.rate(c, self.key(r))
        if rate is None:
            value[3:] += [m, n]
        else:
            value[1] += n*rate


@transactional("summary_manifest.json")
def summarize_common_support(modern, out, minimum=100, coverage=.99, bootstrap=2000, seed=20261027):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    background = Background(architecture_paths(modern, "global"), minimum)
    # Supported observed/expected/denominator, then unsupported obs/denominator.
    blocks = defaultdict(lambda: defaultdict(lambda: np.zeros(5)))
    outside = Counter()
    for path in architecture_paths(modern, "exposure"):
        for r in rows(path):
            if r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}:
                continue
            c = canonical_chrom(r["chrom"])
            if c not in AUTOSOMES:
                continue
            key = r["label"], int(r["threshold"]), r["sensitivity"]
            n, m = int(r["n_sites"]), int(r["n_multi"])
            if r["block_id"] == "outside_partition":
                outside[key] += n
                continue
            v = blocks[key][c + ":" + r["block_id"]]
            rate = background.rate(c, background.key(r))
            if rate is None:
                v[3:] += [m, n]
            else:
                v[:3] += [m, n * rate, n]
    results = []
    for key, vals in sorted(blocks.items()):
        o, e, n, removed_o, removed_n = np.sum(list(vals.values()), axis=0)
        supported = n / (n + removed_n) if n + removed_n else 0.
        usable = {b: v[:2] for b, v in vals.items() if v[2] > 0}
        low, high, p = scan.bootstrap_statistics(usable, bootstrap, seed) if e > 0 else (None, None, None)
        status = "complete" if supported >= coverage and len(usable) >= 20 and e > 0 else "low_support_or_few_blocks"
        results.append(dict(label=key[0], threshold=key[1], sensitivity=key[2],
                            observed_multi=o, expected_multi=e, oe=o/e if e else None,
                            ci_low=low, ci_high=high, bootstrap_p=p, fdr=None,
                            supported_sites=int(n), excluded_sites=int(removed_n), excluded_multi=int(removed_o),
                            supported_fraction=supported, block_count=len(usable),
                            outside_partition_sites=outside[key], status=status))
    tested = [r for r in results if r["status"] == "complete" and r["bootstrap_p"] is not None]
    for r, q in zip(tested, scan.bh_fdr([r["bootstrap_p"] for r in tested])):
        r["fdr"] = q
    write_rows(out / "common_support_enrichment.tsv", results)
    manifest = dict(schema=SCHEMA, status="complete", source=str(Path(modern).resolve()),
                    minimum_training_sites=minimum, minimum_support=coverage, bootstrap=bootstrap,
                    denominator="original eligible SNVs; unsupported strata removed from BOTH numerator and denominator",
                    adjustment="LOCO canonical trinucleotide, primary ALT AF and AN; base denominator fixed across AC tiers",
                    inference="EUR LD-block bootstrap; background offsets treated as fixed; BH over all complete rows",
                    interpretation="cross-sectional multi-ALT occupancy, not selection or temporal persistence")
    dump_json(out / "summary_manifest.json", manifest)
    plot_families(results, out)
    return manifest


def plot_families(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    families = ["ERV1", "ERVK", "ERVL", "ERVL-MaLR"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for t, color, offset in [(5, "#176B87", -.12), (10, "#D16B38", .12)]:
        labelled = False
        for i, family in enumerate(families):
            rr = [r for r in results if r["label"] == "family:"+family and int(r["threshold"]) == t and r["sensitivity"] == "exclude_combined"]
            if rr and rr[0]["oe"] not in (None, "NA"):
                r = rr[0]
                ax.plot(float(r["oe"]), i+offset, "o", color=color, label=f"AC ≥ {t}" if not labelled else None)
                labelled = True
                if r["ci_low"] not in (None, "NA"):
                    ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i+offset]*2, color=color)
    ax.axvline(1, color=".6", ls="--")
    ax.set(yticks=range(4), yticklabels=families, xlabel="Multi-ALT observed / expected", title="HERV families · common support · strict masks")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"family_enrichment.{ext}", dpi=180)
    plt.close(fig)


class Edges:
    """Nearest edge per family/end; base centres avoid strand off-by-one bias."""
    def __init__(self, elements, families):
        self.points = defaultdict(list)
        self.coords = {}
        self.ambiguous_edges = 0
        for e in elements:
            if e.rep_class == "LTR" and e.family in families:
                for end, point in (("5p", e.five_prime), ("3p", e.three_prime)):
                    self.points[e.family, end].append((point, e.strand, e.name, e.source_id))
        for key, v in self.points.items():
            unique = []
            for point, group in itertools.groupby(sorted(v), key=lambda z:z[0]):
                coincident = list(group)
                strands = {z[1] for z in coincident}
                if len(strands) > 1:
                    # The left/right sides must not silently choose different
                    # strands at the same coordinate and create a false peak.
                    self.ambiguous_edges += 1
                    unique.append((point, "?", "ambiguous_orientation", ""))
                else:
                    unique.append(coincident[0])
            self.points[key] = unique
            self.coords[key] = [z[0] for z in unique]

    def locate(self, pos1, flank):
        p = pos1 - .5
        found = []
        for key, vals in self.points.items():
            i = bisect_left(self.coords[key], p)
            candidates = vals[max(0, i-1):min(len(vals), i+1)]
            if not candidates:
                continue
            point, strand, name, source = min(candidates, key=lambda x:(abs(p-x[0]), x))
            if strand == "?":
                continue
            distance = (p-point) * (1 if strand == "+" else -1)
            if -flank <= distance < flank:
                found.append((*key, distance))
        return found

    def intervals(self, flank, length):
        return scan.merge_intervals([(max(0, p-flank), min(length, p+flank))
                                     for vals in self.points.values() for p, *_ in vals
                                     if max(0, p-flank) < min(length, p+flank)])


def query_regions(vcf, fields, bed, bcftools="bcftools"):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER" + "".join("\t%INFO/"+f for f in fields) + "\n"
    cmd = [bcftools, "query", "--allow-undef-tags", "--regions-overlap", "0", "-R", str(bed), "-f", fmt, str(vcf)]
    with tempfile.TemporaryFile(mode="w+t") as errors:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            previous = None
            for line in proc.stdout:
                record = _record(line.rstrip("\r\n").split("\t"), fields)
                key = (canonical_chrom(record["chrom"]), record["pos"])
                if previous and key[0] == previous[0] and key[1] < previous[1]:
                    raise ValueError("Indexed query was not coordinate sorted")
                previous = key
                yield record
            code = proc.wait()
            if code:
                errors.seek(0)
                raise RuntimeError(f"bcftools query failed ({code}): {errors.read()[-2000:]}")
        finally:
            proc.stdout.close()
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()


def strict_site_excluded(chrom, pos, flags, repeats, masks):
    """Exactly match ArchitectureAccumulator's BED-first annotation policy."""
    flags = set(flags)
    excluded = {name: masks[name].contains(chrom, pos) if name in masks else name in flags
                for name in ("lcr", "segdup")}
    excluded["low_mappability"] = (masks["low_mappability"].contains(chrom, pos)
                                  if "low_mappability" in masks else False)
    excluded["simple_repeat"] = any(e.rep_class in {"Simple_repeat", "Low_complexity"} for e in repeats)
    return any(excluded.values())


def read_indexed_header(vcf, bcftools):
    result = subprocess.run([bcftools, "view", "--no-version", "-h", str(vcf)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise ValueError(f"Cannot read VCF header: {result.stderr[-2000:]}")
    return _parse_header(result.stdout)


def scan_boundary_chrom(vcf, chrom, annotation, out, options, args):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"chr{chrom}.bins.tsv.gz"
    marker = out / f"chr{chrom}.json"
    fasta = args.fasta or options["fasta"]
    mask_paths = {n: options.get(k) for n, k in [("segdup", "architecture_segdup_bed"),
                  ("lcr", "architecture_lcr_bed"), ("low_mappability", "architecture_low_mappability_bed")]}
    for name in mask_paths:
        value = getattr(args, name+"_bed", None)
        if value:
            mask_paths[name] = value
    inputs = [vcf, annotation, fasta, fasta+".fai", *[p for p in mask_paths.values() if p]]
    if options.get("callable_bed"):
        inputs.append(options["callable_bed"])
    for suffix in (".tbi", ".csi"):
        if Path(str(vcf)+suffix).is_file():
            inputs.append(str(vcf)+suffix)
            break
    else:
        raise ValueError(f"Indexed VCF required: {vcf}")
    identity = dict(schema=SCHEMA, inputs=[file_state(p) for p in inputs], source_options=options,
                    code=[file_state(p) for p in (Path(__file__), Path(scan.__file__), Path(__file__).with_name("0.common.py"))],
                    families=args.families, flank=args.flank, step=args.step, block_bp=args.block_bp)
    if marker.exists():
        old = json.loads(marker.read_text())
        if old["identity"] != identity:
            raise ValueError(f"Boundary inputs/settings/code changed; use a fresh output: {out}")
        if path.exists() and old.get("output") == file_state(path):
            print(f"[resume] boundary chr{chrom}", flush=True)
            return path
        raise ValueError(f"Completed boundary cache changed: {path}")
    header = read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    edges = Edges(elements, args.families.split(","))
    repeat_sweep = scan.ERVOverlapSweep(elements)
    masks = {n: IntervalMask(p, [chrom]) for n, p in mask_paths.items() if p}
    strict_complete = all(n in masks or n in header["info"] for n in ("lcr", "segdup")) and "low_mappability" in masks
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    qc = Counter()
    partial = Path(str(path)+".partial")
    columns = ["chrom", "block", "family", "boundary", "bin_start", "sensitivity", "threshold",
               "context", "af_bin", "an_bin", "n_sites", "n_multi"]
    bins = defaultdict(lambda: [0, 0])
    last_block = None
    with FastaContext(fasta) as fa, gzip.open(partial, "wt", compresslevel=2, newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(columns)
        def flush():
            writer.writerows([*k, *v] for k, v in sorted(bins.items()))
            bins.clear()
        raw_chroms = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
        if len(raw_chroms) != 1:
            raise ValueError(f"No unique VCF contig for chr{chrom}")
        intervals = edges.intervals(args.flank, fa.lengths[chrom])
        bed = out / f"chr{chrom}.regions.bed"
        with bed.open("w") as h:
            for lo, hi in intervals:
                h.write(f"{raw_chroms[0]}\t{lo}\t{hi}\n")
        records = query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
        clock = time.monotonic()
        for (c, pos), group in itertools.groupby(records, key=lambda r:(canonical_chrom(r["chrom"]), r["pos"])):
            if c != chrom:
                raise ValueError("Unexpected chromosome from targeted query")
            if callable_mask and not callable_mask.contains(c, pos):
                continue
            near = edges.locate(pos, args.flank)
            if not near:
                continue
            ref, alts, non_snv, flags = scan.group_alleles(group, header, [], options, qc, None, [])
            observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                              key=lambda a:(-a["af"], -a["ac"], a["alt"]))
            if not observed:
                continue
            ctx = fa.context(c, pos)
            if ctx[1] != ref:
                raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
            context, _ = scan.canonical_context(ctx, ref)
            if any(b not in "ACGT" for b in context):
                # The original scanner retains these eligible loci. Removing
                # them here would silently change the architecture denominator.
                qc["ambiguous_context_retained"] += 1
            primary = observed[0]
            strata = context, af_bin(primary["af"]), an_bin(primary["an"])
            repeats = repeat_sweep.locate(c, pos)
            excluded = strict_site_excluded(c, pos, flags, repeats, masks)
            sensitivities = ["all"]
            if strict_complete and not excluded:
                sensitivities.append("strict")
            cpg = (ref == "C" and ctx[2] == "G") or (ref == "G" and ctx[0] == "C")
            if not cpg:
                sensitivities += [s+"_no_cpg" for s in list(sensitivities)]
            block = f"{c}:{(pos-1)//args.block_bp}"
            if last_block is not None and last_block != block:
                flush()
            last_block = block
            for family, end, distance in near:
                bin_start = math.floor(distance/args.step)*args.step
                for sensitivity in sensitivities:
                    for threshold in (5, 10):
                        key = (c, block, family, end, bin_start, sensitivity, threshold, *strata)
                        v = bins[key]
                        v[0] += 1
                        v[1] += sum(a["ac"] >= threshold for a in observed) >= 2
            qc["eligible_sites"] += 1
            if time.monotonic()-clock > 60:
                print(f"[boundaries chr{chrom}] {qc['eligible_sites']:,} eligible sites", flush=True)
                clock = time.monotonic()
        flush()
    os.replace(partial, path)
    dump_json(marker, dict(identity=identity, output=file_state(path), qc=dict(qc),
                          strict_masks_complete=strict_complete, ambiguous_oriented_edges=edges.ambiguous_edges,
                          queried_bp=sum(b-a for a,b in intervals)))
    return path


def cluster_estimate(matrix, local=False):
    """Log-ratio cluster sandwich; repeated windows share the same block scores."""
    sums = matrix.sum(axis=0)
    if (sums[:2] <= 0).any() or local and (sums[2:4] <= 0).any():
        return None
    theta = sums[0]/sums[1]
    influence = matrix[:, 0]/sums[0] - matrix[:, 1]/sums[1]
    if local:
        theta /= sums[2]/sums[3]
        influence -= matrix[:, 2]/sums[2] - matrix[:, 3]/sums[3]
    used = (matrix[:, :4 if local else 2].sum(axis=1) > 0).sum()
    correction = math.sqrt(used/max(used-1, 1))
    influence *= correction
    se = float(np.sqrt(np.sum(influence**2)))
    if not se or used < 2:
        return None
    from scipy.stats import norm
    z = math.log(theta)/se
    return dict(estimate=float(theta), ci_low=float(np.exp(math.log(theta)-1.96*se)),
                ci_high=float(np.exp(math.log(theta)+1.96*se)), pointwise_p=float(2*norm.sf(abs(z))),
                z=z, influence=influence/se, n_blocks=int(used))


def profile_counts(paths, background, args):
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: np.zeros(4))))
    for path in paths:
        for r in rows(path):
            c = canonical_chrom(r["chrom"])
            split = "discovery" if c in background.discovery else "validation"
            for scope in ("all", split):
                key = (r["family"], r["boundary"], r["sensitivity"], int(r["threshold"]), scope)
                v = data[key][r["block"]][int(r["bin_start"])]
                n, m = int(r["n_sites"]), int(r["n_multi"])
                rate = background.rate(c, background.key(r), scope)
                v[2] += n
                if rate is None:
                    v[3] += n
                else:
                    v[0] += m
                    v[1] += n*rate
    return data


def infer_profiles(data, args):
    all_blocks = sorted({b for profile in data.values() for b in profile})
    index = {b:i for i,b in enumerate(all_blocks)}
    results, inferential = [], []
    for key, blocks in sorted(data.items()):
        family, end, sensitivity, threshold, split = key
        flank = np.zeros((len(all_blocks), 2))
        flank_n = flank_missing = 0.
        for b, bins in blocks.items():
            for pos, v in bins.items():
                if abs(pos+args.step/2) >= args.flank*.75:
                    flank[index[b]] += v[:2]
                    flank_n += v[2]
                    flank_missing += v[3]
        flank_support = 1-flank_missing/flank_n if flank_n else 0.
        for start in range(-args.flank, args.flank-args.window+1, args.step):
            matrix = np.zeros((len(all_blocks), 4))
            matrix[:, 2:] = flank
            n = missing = 0.
            for b, bins in blocks.items():
                for pos in range(start, start+args.window, args.step):
                    v = bins.get(pos)
                    if v is not None:
                        matrix[index[b], :2] += v[:2]
                        n += v[2]
                        missing += v[3]
            support = 1-missing/n if n else 0.
            # Centre contrasts do not overlap the distal control bins.
            eligible_local = start >= -args.flank*.5 and start+args.window <= args.flank*.5
            global_fit = cluster_estimate(matrix)
            local_fit = cluster_estimate(matrix, True) if eligible_local else None
            enough = support >= args.min_support and matrix[:, 0].sum() >= 20
            local_enough = (flank_support >= args.min_support and matrix[:, 2].sum() >= 20)
            record = dict(family=family, boundary=end, sensitivity=sensitivity, threshold=threshold, split=split,
                          start_bp=start, end_bp=start+args.window, observed_multi=matrix[:, 0].sum(),
                          expected_multi=matrix[:, 1].sum(), eligible_sites=int(n), excluded_sites=int(missing),
                          supported_sites=int(n-missing),
                          supported_fraction=support, flank_supported_fraction=flank_support,
                          oe=None, ci_low=None, ci_high=None, local_ratio=None,
                          local_ci_low=None, local_ci_high=None, local_p=None,
                          discovery_maxT_p=None, n_blocks=0, status="insufficient_support")
            if global_fit:
                record.update(oe=global_fit["estimate"], ci_low=global_fit["ci_low"],
                              ci_high=global_fit["ci_high"], n_blocks=global_fit["n_blocks"])
            if local_fit and local_enough:
                record.update(local_ratio=local_fit["estimate"], local_ci_low=local_fit["ci_low"],
                              local_ci_high=local_fit["ci_high"], local_p=local_fit["pointwise_p"])
            if enough and global_fit and global_fit["n_blocks"] >= args.min_blocks:
                record["status"] = "complete"
            results.append(record)
            if (record["status"] == "complete" and local_fit and local_enough
                    and local_fit["n_blocks"] >= args.min_blocks and sensitivity == "strict"
                    and threshold == 5 and split == "discovery"):
                inferential.append((record, local_fit))
    if inferential:
        influences = np.column_stack([x[1]["influence"] for x in inferential])
        rng = np.random.default_rng(args.seed)
        maximum = []
        for begin in range(0, args.bootstrap, 100):
            multipliers = rng.choice([-1., 1.], size=(min(100, args.bootstrap-begin), len(all_blocks)))
            maximum.extend(np.max(np.abs(multipliers @ influences), axis=1))
        maximum = np.array(maximum)
        for record, fit in inferential:
            record["discovery_maxT_p"] = (1+int(np.sum(maximum >= abs(fit["z"]))))/(args.bootstrap+1)
    peaks = []
    for family in args.families.split(","):
        for end in ("5p", "3p"):
            choices = [r for r,_ in inferential if r["family"] == family and r["boundary"] == end]
            if not choices:
                continue
            # Freeze location on discovery chromosomes before looking at validation.
            best = max(choices, key=lambda r:r["local_ratio"])
            for threshold in (5, 10):
                matches = [r for r in results if r["family"] == family and r["boundary"] == end
                           and r["sensitivity"] == "strict" and r["threshold"] == threshold
                           and r["split"] == "validation" and r["start_bp"] == best["start_bp"]]
                val = matches[0] if matches else None
                peaks.append(dict(family=family, boundary=end, start_bp=best["start_bp"], end_bp=best["end_bp"],
                                  discovery_local_ratio=best["local_ratio"], discovery_maxT_p=best["discovery_maxT_p"],
                                  validation_threshold=threshold, validation_local_ratio=val["local_ratio"] if val else None,
                                  validation_ci_low=val["local_ci_low"] if val else None,
                                  validation_ci_high=val["local_ci_high"] if val else None,
                                  validation_p=val["local_p"] if val and val["status"] == "complete" else None,
                                  validation_fdr=None, status="exploratory_or_unconfirmed"))
    valid = [r for r in peaks if r["validation_threshold"] == 5 and r["validation_p"] is not None]
    for r, q in zip(valid, scan.bh_fdr([r["validation_p"] for r in valid])):
        r["validation_fdr"] = q
        if (r["discovery_maxT_p"] < .05 and q < .05
                and r["discovery_local_ratio"] > 1 and r["validation_local_ratio"] > 1):
            r["status"] = "held_out_chromosome_confirmation"
    return results, peaks


def plot_profiles(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    families = sorted({r["family"] for r in results})
    if not families:
        return
    fig, axes = plt.subplots(len(families), 2, figsize=(11, 3.4*len(families)), squeeze=False)
    for i, family in enumerate(families):
        for j, end in enumerate(("5p", "3p")):
            ax = axes[i, j]
            for split, color in [("discovery", "#176B87"), ("validation", "#D16B38")]:
                rr = sorted([r for r in results if r["family"] == family and r["boundary"] == end
                             and r["sensitivity"] == "strict" and r["threshold"] == 5
                             and r["split"] == split and r["oe"] is not None], key=lambda r:r["start_bp"])
                if rr:
                    x = [(r["start_bp"]+r["end_bp"])/2 for r in rr]
                    ax.plot(x, [r["oe"] for r in rr], color=color, label=split)
                    ax.fill_between(x, [r["ci_low"] for r in rr], [r["ci_high"] for r in rr], color=color, alpha=.13)
            ax.axhline(1, ls="--", color=".6")
            ax.axvline(0, ls=":", color=".5")
            ax.set(title=f"{family} · annotated {end} edge", xlabel="Strand-oriented distance (bp)", ylabel="Multi-ALT O/E")
            if ax.get_legend_handles_labels()[0]:
                ax.legend(frameon=False)
    fig.suptitle("AC ≥ 5 · strict masks · fixed discovery / validation chromosomes", y=1.01)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"boundary_profiles.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def assert_annotation_inputs(options, args):
    """Cached global denominators cannot be paired with a newly filtered scan."""
    mappings = [("fasta", "fasta"), ("erv_joined", "erv_joined"),
                ("segdup_bed", "architecture_segdup_bed"),
                ("lcr_bed", "architecture_lcr_bed"),
                ("low_mappability_bed", "architecture_low_mappability_bed")]
    for argument, option in mappings:
        override = getattr(args, argument, None)
        original = options.get(option)
        if override and (not original or Path(override).expanduser().resolve() != Path(original).expanduser().resolve()):
            raise ValueError(f"--{argument.replace('_', '-')} differs from the original modern scan. "
                             "Use the original input with its cached global denominator, or make a separate modern scan.")


def verify_recorded_file(metadata):
    """Use the old run's metadata receipt; do not introduce full-file hashing."""
    if not isinstance(metadata, dict) or not {"path", "size", "mtime_ns"} <= set(metadata):
        return
    current = file_state(metadata["path"])
    if current["bytes"] != metadata["size"] or current["mtime_ns"] != metadata["mtime_ns"]:
        raise ValueError(f"Original modern input changed since its cached global counts: {metadata['path']}")


def targeted_inputs(args, allowed_chromosomes=None):
    """Resolve an adopted store without changing its scientific identity."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    run = scan.load_run(args.modern_dir)
    options = run["options"]
    for argument in ("fasta", "erv_joined", "segdup_bed", "lcr_bed", "low_mappability_bed"):
        if getattr(args, argument, None):
            setattr(args, argument, str(scan.resolve_run_path(run, getattr(args, argument))))
    assert_annotation_inputs(options, args)
    receipt = run.get("identity", {})
    for name in ("fasta", "fasta_index", "fasta_gzi", "callable_bed"):
        if receipt.get(name):
            scan.verify_recorded_input(run, receipt[name])
    annotation_receipts = (options.get("architecture_annotations") or {}).get("files", {})
    for name in ("erv_joined", "architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed", "EUR"):
        if annotation_receipts.get(name):
            scan.verify_recorded_input(run, annotation_receipts[name])
    if str(options.get("build", "38")) not in {"38", "GRCh38", "hg38"}:
        raise ValueError("HERV targeted analysis requires GRCh38 modern input")
    if options.get("region"):
        raise ValueError("Use completed chromosome/genome shards, not a region-restricted modern scan")
    annotation = options.get("architecture_annotation_cache")
    if not annotation:
        raise ValueError("Completed scan lacks its annotation cache; restore the cache rather than rebuild frozen inputs")
    annotation = Path(scan.resolve_run_path(run, annotation))
    scan.validate_existing_joined_cache(args.erv_joined or options["erv_joined"], annotation)
    allowed = set(allowed_chromosomes) if allowed_chromosomes is not None else AUTOSOMES
    selected = {canonical_chrom(c.strip()) for c in args.chromosomes.split(",")} if args.chromosomes else allowed
    if not selected <= allowed:
        raise ValueError(f"Unsupported targeted chromosomes: {sorted(selected - allowed)}")
    receipts = {str(Path(scan.resolve_run_path(run, v["path"])).resolve()): v
                for v in receipt.get("inputs", []) if isinstance(v, dict) and "path" in v}
    tasks, used, original = [], set(), set()
    for shard in run["shards"]:
        info = json.loads((Path(shard)/"scan.json").read_text())
        for c in info["scope"]:
            c = canonical_chrom(c)
            if c not in allowed:
                continue
            original.add(c)
            if c not in selected:
                continue
            if c in used:
                raise ValueError("Multiple original shards cover one chromosome; use non-overlapping completed chromosome shards")
            used.add(c)
            source = scan.resolve_run_path(run, info["input"])
            metadata = receipts.get(str(Path(source).resolve()))
            if metadata:
                scan.verify_recorded_input(run, metadata)
            tasks.append((c, source))
    if not tasks:
        raise ValueError("No matching completed chromosome shards")
    if args.chromosomes and selected-used:
        raise ValueError(f"Requested chromosomes absent from completed original shards: {sorted(selected-used)}")
    return run, annotation, tasks, original


@transactional("windows_manifest.json")
def run_windows(args):
    if args.step <= 0 or args.window < args.step or args.window % args.step or args.flank % args.step:
        raise ValueError("Positive window/flank must be integer multiples of step")
    if args.window > args.flank or args.block_bp < 100000 or args.bootstrap < 200:
        raise ValueError("window <= flank, block >=100kb and >=200 multiplier draws required")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    run, annotation, tasks, _ = targeted_inputs(args)
    options = run["options"]
    files, used = [], {c for c, _ in tasks}
    for chrom, source in tasks:
        print(f"[boundaries] indexed chr{chrom} query", flush=True)
        files.append(scan_boundary_chrom(source, chrom, annotation/f"chr{chrom}.tsv.gz", out/"cache", options, args))
    background = Background(architecture_paths(args.modern_dir, "global"), args.min_training_sites,
                            args.discovery_chromosomes.split(","))
    data = profile_counts(files, background, args)
    results, peaks = infer_profiles(data, args)
    write_rows(out/"boundary_profiles.tsv", results)
    write_rows(out/"frozen_peak_validation.tsv", peaks)
    plot_profiles(results, out)
    manifest = dict(schema=SCHEMA, status="complete", chromosomes=sorted(used, key=int),
                    families=args.families, window=args.window, step=args.step, flank=args.flank,
                    discovery_chromosomes=sorted(background.discovery, key=int),
                    validation_background="frozen discovery chromosomes only",
                    inference=f"{args.block_bp:g}bp physical-block sandwich CI; Rademacher max-|Z| over strict AC5 discovery family/end/central windows; held-out AC5 BH confirmation",
                    local_control=f"outer 25% of +/-{args.flank} bp; local_ratio = window O/E divided by distal-flank O/E",
                    min_support=args.min_support, minimum_training_sites=args.min_training_sites,
                    caveats=["Observed joined alignment edges are not TSS or intact provirus boundaries",
                             "Nearest edge per family/end; no repeated site in one window",
                             "Coincident same-family/end edges with conflicting strands are excluded from oriented assignment",
                             "Original SNV denominator, not callable-bp mutation rate",
                             "Background offsets treated as fixed; physical blocks approximate dependence",
                             "Local GC, repeat age and recombination not included in this O/E",
                             "Chromosome split is internal validation, not an independent cohort"])
    dump_json(out/"windows_manifest.json", manifest)
    return manifest


SPATIAL_SCHEMA = "herv-spatial-1.0"
PARTITIONS = ("head", "body", "tail", "short_element", "ambiguous_overlap")


class AlignedBodies:
    """One assignment per family/site; no assumption of an intact provirus.

    Distance follows the strand along the concatenated *aligned* bases of a
    joined RepeatMasker record. Gaps do not acquire body positions. Conflicting
    overlapping annotations retain their counts in an explicit residual.
    """
    def __init__(self, elements, families, end_bp=200, bins=20):
        if end_bp < 1 or bins < 2:
            raise ValueError("Positive end-bp and at least two body bins required")
        self.families = set(families)
        self.end_bp, self.bins = end_bp, bins
        self.elements = [e for e in elements if e.rep_class == "LTR" and e.family in self.families]
        self.geometry = {}
        for e in self.elements:
            segments = tuple(scan.merge_intervals(list(e.segments or ((e.start, e.end),))))
            self.geometry[e] = (segments, sum(hi-lo for lo, hi in segments))

    def intervals(self):
        return scan.merge_intervals([segment for segments, _ in self.geometry.values() for segment in segments])

    def coordinate(self, e, pos1):
        segments, length = self.geometry[e]
        prefix, p = 0, int(pos1)-1
        for lo, hi in segments:
            if lo <= p < hi:
                # Twice the base-centre distance is integer and exactly mirrors
                # on the reverse strand, including at normalized bin boundaries.
                distance2 = 2*prefix + 2*(p-lo) + 1
                if e.strand == "-":
                    distance2 = 2*length-distance2
                if length <= 2*self.end_bp:
                    partition = "short_element"
                elif distance2 < 2*self.end_bp:
                    partition = "head"
                elif distance2 >= 2*(length-self.end_bp):
                    partition = "tail"
                else:
                    partition = "body"
                return partition, distance2*self.bins//(2*length)
            prefix += hi-lo
        raise ValueError("A body coordinate was requested for an unaligned base")

    def assign(self, overlaps, pos1):
        by_family = defaultdict(dict)
        for e in overlaps:
            if e.rep_class != "LTR" or e.family not in self.families:
                continue
            segments, _ = self.geometry[e]
            # Name/source aliases at identical geometry do not create replicates.
            geometry = e.strand, segments
            by_family[e.family][geometry] = self.coordinate(e, pos1)
        result = {}
        for family, geometries in by_family.items():
            strands = {g[0] for g in geometries}
            parts = {v[0] for v in geometries.values()}
            bins = {v[1] for v in geometries.values()}
            part = next(iter(parts)) if len(strands) == len(parts) == 1 else "ambiguous_overlap"
            body_bin = next(iter(bins)) if len(strands) == len(bins) == 1 else -1
            overlap = ("single_geometry" if len(geometries) == 1 else
                       "compatible_overlap" if part != "ambiguous_overlap" and body_bin >= 0 else "ambiguous_overlap")
            result[family] = part, body_bin, overlap
        return result


def spatial_annotation_audit(elements, families):
    """Report exact labels, gaps and overlap, without merging biological aliases."""
    output = []
    groups = defaultdict(list)
    for e in elements:
        if e.rep_class == "LTR":
            groups[e.chrom, e.family].append(e)
    for (chrom, family), ee in sorted(groups.items()):
        geometries = [(e.strand, tuple(scan.merge_intervals(list(e.segments or ((e.start, e.end),))))) for e in ee]
        lengths = [sum(t-s for s, t in segments) for _, segments in geometries]
        aligned = sum(t-s for s, t in scan.merge_intervals([p for _, segments in geometries for p in segments]))
        spans = sum(t-s for s, t in scan.merge_intervals([(e.start, e.end) for e in ee]))
        output.append(dict(chrom=chrom, rep_class="LTR", raw_family=family,
                           counter_label="family:"+(family or "unspecified"), selected_exact=family in families,
                           records=len(ee), identical_geometry_extra_records=len(ee)-len(set(geometries)),
                           aligned_record_sum_bp=sum(lengths), aligned_union_bp=aligned,
                           envelope_union_bp=spans, envelope_only_gap_bp=spans-aligned,
                           overlap_multiplicity_extra_bp=sum(lengths)-aligned,
                           min_aligned_bp=min(lengths), median_aligned_bp=float(np.median(lengths)),
                           max_aligned_bp=max(lengths)))
    return output


def scan_spatial_chrom(vcf, chrom, annotation, out, options, args):
    """One indexed query over the family-aligned union, never a genome rescan."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path, marker = out/f"chr{chrom}.spatial.tsv.gz", out/f"chr{chrom}.spatial.json"
    fasta = str(args.fasta or options["fasta"])
    mask_paths = {n: options.get(k) for n, k in (("segdup", "architecture_segdup_bed"),
                  ("lcr", "architecture_lcr_bed"), ("low_mappability", "architecture_low_mappability_bed"))}
    for name in mask_paths:
        mask_paths[name] = getattr(args, name+"_bed", None) or mask_paths[name]
    ld_path = (Path(options["ld_block_dir"])/f"EUR.{options.get('architecture_build', '38')}.bed"
               if options.get("ld_block_dir") else None)
    inputs = [vcf, annotation, fasta, fasta+".fai", *[p for p in mask_paths.values() if p]]
    if ld_path:
        inputs.append(ld_path)
    if options.get("callable_bed"):
        inputs.append(options["callable_bed"])
    for suffix in (".tbi", ".csi"):
        if Path(str(vcf)+suffix).is_file():
            inputs.append(str(vcf)+suffix)
            break
    else:
        raise ValueError(f"Indexed VCF required: {vcf}")
    identity = dict(schema=SPATIAL_SCHEMA, inputs=[file_state(p) for p in inputs], source_options=options,
                    code=[file_state(p) for p in (Path(__file__), Path(scan.__file__), Path(__file__).with_name("0.common.py"))],
                    families=args.families, end_bp=args.end_bp, body_bins=args.body_bins, block_bp=args.block_bp)
    if marker.exists():
        old = json.loads(marker.read_text())
        if old["identity"] != identity:
            raise ValueError(f"Spatial inputs/settings/code changed; use a fresh output: {out}")
        expected = {Path(p).name:state for p, state in old.get("outputs", {}).items()}
        required = [path, out/f"chr{chrom}.annotation_audit.tsv", out/f"chr{chrom}.elements.tsv.gz"]
        if all(p.is_file() and expected.get(p.name) == file_state(p) for p in required):
            print(f"[resume] spatial chr{chrom}", flush=True)
            return path
        raise ValueError(f"Completed spatial cache changed: {path}")
    header = read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"]+[f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    families = args.families.split(",")
    body = AlignedBodies(elements, families, args.end_bp, args.body_bins)
    repeat_sweep = scan.ERVOverlapSweep(elements)
    masks = {n: IntervalMask(p, [chrom]) for n, p in mask_paths.items() if p}
    strict_complete = all(n in masks or n in header["info"] for n in ("lcr", "segdup")) and "low_mappability" in masks
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    eur = scan.BlockSweep(scan.read_blocks(ld_path, "EUR", options.get("architecture_build", "38"))[0]) if ld_path else None
    intervals = body.intervals()
    raw_chroms = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
    if len(raw_chroms) != 1:
        raise ValueError(f"No unique VCF contig for chr{chrom}")
    bed = out/f"chr{chrom}.spatial.regions.bed"
    with bed.open("w") as handle:
        for lo, hi in intervals:
            handle.write(f"{raw_chroms[0]}\t{lo}\t{hi}\n")
    columns = ["chrom", "block", "eur_partition", "family", "partition", "body_bin", "overlap_assignment",
               "sensitivity", "threshold", "context", "af_bin", "an_bin", "n_sites", "n_multi"]
    counts, qc, last_block = defaultdict(lambda: [0, 0]), Counter(), None
    partial = Path(str(path)+".partial")
    with FastaContext(fasta) as fa, gzip.open(partial, "wt", compresslevel=2, newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(columns)
        def flush():
            writer.writerows([*k, *v] for k, v in sorted(counts.items()))
            counts.clear()
        records = query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
        clock = time.monotonic()
        for (c, pos), group in itertools.groupby(records, key=lambda r:(canonical_chrom(r["chrom"]), r["pos"])):
            if c != chrom:
                raise ValueError("Unexpected chromosome from targeted spatial query")
            if callable_mask and not callable_mask.contains(c, pos):
                qc["outside_callable_mask"] += 1
                continue
            overlaps = repeat_sweep.locate(c, pos)
            assignment = body.assign(overlaps, pos)
            if not assignment:
                raise ValueError("Indexed body query returned a position outside the selected aligned family union")
            ref, alts, _, flags = scan.group_alleles(group, header, [], options, qc, None, [])
            observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                              key=lambda a:(-a["af"], -a["ac"], a["alt"]))
            if not observed:
                continue
            ctx = fa.context(c, pos).upper()
            if len(ctx) != 3 or ctx[1] != ref:
                raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
            context, _ = scan.canonical_context(ctx, ref)
            if any(b not in "ACGT" for b in context):
                qc["ambiguous_context_retained"] += 1
            primary = observed[0]
            strata = context, af_bin(primary["af"]), an_bin(primary["an"])
            sensitivities = ["all"]
            if strict_complete and not strict_site_excluded(c, pos, flags, overlaps, masks):
                sensitivities.append("strict")
            if not ((ref == "C" and ctx[2] == "G") or (ref == "G" and ctx[0] == "C")):
                sensitivities += [s+"_no_cpg" for s in list(sensitivities)]
            covered = "inside" if eur and eur.locate(c, pos) else "outside"
            block = f"{c}:{(pos-1)//args.block_bp}"
            if last_block is not None and last_block != block:
                flush()
            last_block = block
            multi = {t: int(sum(a["ac"] >= t for a in observed) >= 2) for t in (5, 10)}
            for family, (part, body_bin, overlap) in assignment.items():
                qc[f"{family}:{overlap}"] += 1
                for sensitivity in sensitivities:
                    for threshold in (5, 10):
                        key = (c, block, covered, family, part, body_bin, overlap, sensitivity, threshold, *strata)
                        counts[key][0] += 1
                        counts[key][1] += multi[threshold]
            qc["eligible_sites_in_union"] += 1
            if time.monotonic()-clock > 60:
                print(f"[spatial chr{chrom}] {qc['eligible_sites_in_union']:,} eligible sites", flush=True)
                clock = time.monotonic()
        flush()
    os.replace(partial, path)
    audit_path, element_path = out/f"chr{chrom}.annotation_audit.tsv", out/f"chr{chrom}.elements.tsv.gz"
    write_rows(audit_path, spatial_annotation_audit(elements, families))
    write_rows(element_path, [dict(chrom=e.chrom, start=e.start, end=e.end, rep_name=e.name,
                                   rep_class=e.rep_class, rep_family=e.family, strand=e.strand, source_id=e.source_id,
                                   aligned_bp=body.geometry[e][1], envelope_bp=e.end-e.start,
                                   unaligned_gap_bp=e.end-e.start-body.geometry[e][1],
                                   name_ends_int=e.name.endswith("-int"),
                                   short_element=body.geometry[e][1] <= 2*args.end_bp,
                                   aligned_segments=",".join(f"{lo}-{hi}" for lo, hi in body.geometry[e][0])) for e in body.elements],
               fields=["chrom", "start", "end", "rep_name", "rep_class", "rep_family", "strand", "source_id",
                       "aligned_bp", "envelope_bp", "unaligned_gap_bp", "name_ends_int", "short_element", "aligned_segments"])
    dump_json(marker, dict(identity=identity, outputs={str(p):file_state(p) for p in (path, audit_path, element_path)},
                          qc=dict(qc), strict_masks_complete=strict_complete, eur_ld_partition_available=eur is not None,
                          queried_aligned_bp=sum(hi-lo for lo, hi in intervals)))
    return path


def spatial_counts(paths, background, with_eur=True, collect_raw=True):
    """Apply the same LOCO rate to every disjoint part of each stratum.

    Values are supported observed/expected/sites and unsupported observed/sites.
    Unsupported loci are never left in the observed numerator. Totals are
    accumulated independently of either partitioning, for accounting checks.
    """
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: np.zeros(5))))
    raw = defaultdict(lambda: np.zeros(2, dtype=np.int64))
    for path in paths:
        for r in rows(path):
            c, n, m = canonical_chrom(r["chrom"]), int(r["n_sites"]), int(r["n_multi"])
            rate = background.rate(c, background.key(r))
            value = np.array([m, n*rate, n, 0, 0]) if rate is not None else np.array([0, 0, 0, m, n])
            scopes = ["all_family_sites"]
            if with_eur and r["eur_partition"] == "inside":
                scopes.append("eur_ld_partition")
            elif with_eur:
                # Keep an explicitly empty legacy-covered profile rather than
                # silently substituting an all-sites result for an absent map.
                prefix = r["family"], r["sensitivity"], int(r["threshold"]), "eur_ld_partition"
                data[(*prefix, "total")]["total"]
            for scope in scopes:
                prefix = r["family"], r["sensitivity"], int(r["threshold"]), scope
                for view, category in (("partition", r["partition"]), ("normalized", int(r["body_bin"])), ("total", "total")):
                    data[(*prefix, view)][category][r["block"]] += value
                if collect_raw and r["sensitivity"] in {"all", "strict"}:
                    raw[(c, *prefix, r["context"], r["af_bin"], r["an_bin"])] += [n, m]
    return data, raw


def _spatial_record(blocks, args):
    vals = np.array(list(blocks.values())) if blocks else np.zeros((0, 5))
    o, e, n, removed_o, removed_n = vals.sum(axis=0)
    eligible = n+removed_n
    fit = cluster_estimate(vals[:, :2]) if e > 0 and o > 0 else None
    count = int(np.sum(vals[:, 2] > 0))
    support = n/eligible if eligible else 0.
    enough = bool(support >= args.min_support and count >= args.min_blocks and o >= 20 and fit)
    return dict(observed_multi=float(o), expected_multi=float(e), oe=float(o/e) if e > 0 else None,
                ci_low=fit["ci_low"] if enough else None, ci_high=fit["ci_high"] if enough else None,
                supported_sites=int(n), excluded_sites=int(removed_n), excluded_multi=int(removed_o),
                eligible_sites=int(eligible), raw_multi=int(o+removed_o), supported_fraction=float(support),
                n_blocks=count, status="complete_descriptive" if enough else "no_sites" if eligible == 0 else "low_support_or_few_events_blocks")


def summarize_spatial(data, args):
    partitions, profiles, identities, contrasts = [], [], [], []
    prefixes = sorted({key[:4] for key in data})
    for prefix in prefixes:
        family, sensitivity, threshold, scope = prefix
        dimensions = dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope)
        total_blocks = data.get((*prefix, "total"), {}).get("total", {})
        total = _spatial_record(total_blocks, args)
        total_e = total["expected_multi"]
        for view, categories, target in (("partition", PARTITIONS, partitions),
                                          ("normalized", tuple(range(args.body_bins))+(-1,), profiles)):
            records = []
            for category in categories:
                blocks = data.get((*prefix, view), {}).get(category, {})
                record = dimensions | _spatial_record(blocks, args)
                if view == "partition":
                    record["partition"] = category
                else:
                    record.update(body_bin=category,
                                  start_fraction=category/args.body_bins if category >= 0 else None,
                                  end_fraction=(category+1)/args.body_bins if category >= 0 else None,
                                  position_label="ambiguous_overlap" if category < 0 else "aligned_body")
                record.update(expected_weight=record["expected_multi"]/total_e if total_e else None,
                              oe_contribution=record["observed_multi"]/total_e if total_e else None,
                              excess_contribution=(record["observed_multi"]-record["expected_multi"])/total_e if total_e else None)
                records.append(record)
            target.extend(records)
            # Verify both supported and unsupported accounting, not just ratios.
            fields = ("observed_multi", "expected_multi", "supported_sites", "excluded_sites", "excluded_multi")
            deltas = {field:sum(r[field] for r in records)-total[field] for field in fields}
            if any(not math.isclose(delta, 0., rel_tol=0., abs_tol=max(1e-8, abs(total[field])*1e-10))
                   for field, delta in deltas.items()):
                raise ValueError(f"Spatial accounting identity failed for {prefix}/{view}: {deltas}")
            weighted = sum(r["oe_contribution"] for r in records) if total_e else None
            ratios = [r["oe"] for r in records if r["oe"] is not None]
            identities.append(dimensions | dict(view=view, overall_oe=total["oe"],
                              expected_weighted_oe=weighted,
                              unweighted_mean_of_available_part_oe=float(np.mean(ratios)) if ratios else None,
                              observed_delta=deltas["observed_multi"], expected_delta=deltas["expected_multi"],
                              supported_sites_delta=deltas["supported_sites"], excluded_sites_delta=deltas["excluded_sites"],
                              excluded_multi_delta=deltas["excluded_multi"],
                              summed_excess_contribution=sum(r["excess_contribution"] for r in records) if total_e else None,
                              status="exact_accounting"))
            if view == "partition":
                target.append(dimensions | total | dict(partition="total", expected_weight=1. if total_e else None,
                              oe_contribution=total["oe"], excess_contribution=total["oe"]-1 if total_e else None))
        # Prespecified comparisons among the unambiguous long-record parts.
        part_data = data.get((*prefix, "partition"), {})
        block_names = sorted(total_blocks)
        for part in ("head", "tail"):
            matrix = np.zeros((len(block_names), 4))
            for i, block in enumerate(block_names):
                matrix[i, :2] = part_data.get(part, {}).get(block, np.zeros(5))[:2]
                matrix[i, 2:] = part_data.get("body", {}).get(block, np.zeros(5))[:2]
            fit = cluster_estimate(matrix, local=True)
            a = _spatial_record(part_data.get(part, {}), args)
            b = _spatial_record(part_data.get("body", {}), args)
            valid = fit and a["status"] == b["status"] == "complete_descriptive"
            contrasts.append(dimensions | dict(contrast=part+"_vs_body", ratio_of_oe=fit["estimate"] if fit else None,
                             ci_low=fit["ci_low"] if valid else None, ci_high=fit["ci_high"] if valid else None,
                             pointwise_p=fit["pointwise_p"] if valid else None, primary_fdr=None,
                             n_blocks=fit["n_blocks"] if fit else 0,
                             status="prespecified_association" if valid else "insufficient_support"))
    primary_scope = "eur_ld_partition" if any(r["scope"] == "eur_ld_partition" for r in contrasts) else "all_family_sites"
    tested = [r for r in contrasts if r["scope"] == primary_scope and r["sensitivity"] == "strict"
              and r["threshold"] == 5 and r["pointwise_p"] is not None]
    for r, q in zip(tested, scan.bh_fdr([r["pointwise_p"] for r in tested])):
        r["primary_fdr"] = q
    return partitions, profiles, identities, contrasts


def _primary_spatial_record(blocks, args, total_without_loco):
    record = _spatial_record(blocks, args)
    supported_m = record["observed_multi"]
    missing_n, missing_m = record["excluded_sites"], record["excluded_multi"]
    partial_e = record["expected_multi"]
    record.update(model="primary_full_scope", observed_multi=float(record["raw_multi"]),
                  supported_observed_multi=supported_m, partial_expected_multi=partial_e,
                  excluded_sites=0, excluded_multi=0, n_without_loco=missing_n,
                  multi_without_loco=missing_m, total_without_loco=int(total_without_loco))
    if total_without_loco:
        # Even a supported head cannot be presented as a weighted share of an
        # unavailable full-scope total. Raw counts and partial E remain auditable.
        record.update(expected_multi=None, oe=None, ci_low=None, ci_high=None,
                      status="unavailable_primary_loco")
    return record


def summarize_primary_spatial(data, args, scoped_exposure=None):
    """Decompose the primary full-scope estimator without substituting common support."""
    partitions, profiles, identities = [], [], []
    prefixes = sorted({key[:4] for key in data if key[1] in {"all", "strict"}
                       and key[3] == "eur_ld_partition"})
    for prefix in prefixes:
        family, sensitivity, threshold, scope = prefix
        dimensions = dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope)
        total_blocks = data.get((*prefix, "total"), {}).get("total", {})
        missing = sum(v[4] for v in total_blocks.values())
        if scoped_exposure is not None:
            # A chromosome subset cannot manufacture a decomposition of an
            # unavailable full-scope total by omitting its unsupported strata.
            missing = max(missing, scoped_exposure.get(prefix[:3], np.zeros(5))[4])
        total = _primary_spatial_record(total_blocks, args, missing)
        total_e = total["expected_multi"]
        for view, categories, target in (("partition", PARTITIONS, partitions),
                                          ("normalized", tuple(range(args.body_bins))+(-1,), profiles)):
            records = []
            for category in categories:
                blocks = data.get((*prefix, view), {}).get(category, {})
                record = dimensions | _primary_spatial_record(blocks, args, missing)
                if view == "partition":
                    record["partition"] = category
                else:
                    record.update(body_bin=category, start_fraction=category/args.body_bins if category >= 0 else None,
                                  end_fraction=(category+1)/args.body_bins if category >= 0 else None,
                                  position_label="ambiguous_overlap" if category < 0 else "aligned_body")
                record.update(expected_weight=record["expected_multi"]/total_e if total_e else None,
                              oe_contribution=record["observed_multi"]/total_e if total_e else None,
                              excess_contribution=(record["observed_multi"]-record["expected_multi"])/total_e if total_e else None)
                records.append(record)
            target.extend(records)
            fields = ("observed_multi", "eligible_sites", "n_without_loco", "multi_without_loco", "partial_expected_multi")
            deltas = {field:sum(r[field] for r in records)-total[field] for field in fields}
            if any(not math.isclose(delta, 0., rel_tol=0., abs_tol=max(1e-8, abs(total[field])*1e-10))
                   for field, delta in deltas.items()):
                raise ValueError(f"Primary spatial accounting failed for {prefix}/{view}: {deltas}")
            identities.append(dimensions | dict(model="primary_full_scope", view=view,
                              overall_oe=total["oe"], expected_weighted_oe=sum(r["oe_contribution"] for r in records) if total_e else None,
                              observed_delta=deltas["observed_multi"], eligible_sites_delta=deltas["eligible_sites"],
                              expected_delta=sum(r["expected_multi"] for r in records)-total_e if total_e is not None else None,
                              partial_expected_delta=deltas["partial_expected_multi"],
                              n_without_loco_delta=deltas["n_without_loco"], total_without_loco=int(missing),
                              summed_excess_contribution=sum(r["excess_contribution"] for r in records) if total_e else None,
                              status="raw_accounting_only_primary_loco_unavailable" if missing else "exact_accounting"))
            if view == "partition":
                target.append(dimensions | total | dict(partition="total", expected_weight=1. if total_e else None,
                              oe_contribution=total["oe"], excess_contribution=total["oe"]-1 if total_e else None))
    return partitions, profiles, identities


def reconcile_spatial_counters(modern, raw, background, chromosomes, families, with_eur=True,
                             comparison_scope="complete_original_autosomes", primary_background=None):
    """Compare raw N/M at every old family stratum before interpreting O/E.

    The old family counters include all aligned body sites, while the old formal
    exposure summary excludes sites outside the EUR partition. Both are checked.
    Expected counts are recomputed under the new common-support settings; this
    does not relabel the primary full-scope O/E as a newly reproduced estimate.
    """
    legacy = defaultdict(lambda: np.zeros(2, dtype=np.int64))
    labels = Counter()
    selected_labels = {"family:"+family for family in families}
    sources = [("erv", "all_family_sites")]
    if with_eur:
        sources.append(("exposure", "eur_ld_partition"))
    for kind, scope in sources:
        for path in architecture_paths(modern, kind):
            for r in rows(path):
                if kind == "exposure" and primary_background is not None:
                    primary_background.observe_exposure(r)
                c = canonical_chrom(r["chrom"])
                if c not in chromosomes or r["population_mode"] != "all":
                    continue
                if kind == "erv" and r["label"].startswith("family:") and r["threshold"] == "5" and r["sensitivity"] == "all":
                    labels[r["label"]] += int(r["n_sites"])
                if r["label"] not in selected_labels or r["threshold"] not in {"5", "10"}:
                    continue
                if r["sensitivity"] not in {"all", "exclude_combined"}:
                    continue
                if kind == "exposure" and r["block_id"] == "outside_partition":
                    continue
                sensitivity = "strict" if r["sensitivity"] == "exclude_combined" else "all"
                key = (c, r["label"][7:], sensitivity, int(r["threshold"]), scope, r["context"], r["af_bin"], r["an_bin"])
                legacy[key] += [int(r["n_sites"]), int(r["n_multi"])]
    sums = defaultdict(lambda: np.zeros((2, 7)))
    mismatched = Counter()
    strata_counts = Counter()
    mismatch_rows = []
    for key in sorted(set(raw) | set(legacy)):
        c, family, sensitivity, threshold, scope, context, af, an = key
        group = family, sensitivity, threshold, scope
        old, new = legacy.get(key, np.zeros(2)), raw.get(key, np.zeros(2))
        strata_counts[group] += 1
        if not np.array_equal(old, new):
            mismatched[group] += 1
            mismatch_rows.append(dict(chrom=c, family=family, sensitivity=sensitivity, threshold=threshold, scope=scope,
                                      context=context, af_bin=af, an_bin=an, counter_sites=int(old[0]), target_sites=int(new[0]),
                                      counter_multi=int(old[1]), target_multi=int(new[1])))
        rate = background.rate(c, background.key(dict(context=context, af_bin=af, an_bin=an,
                                                      sensitivity=sensitivity, threshold=threshold)))
        for i, (n, m) in enumerate((old, new)):
            sums[group][i] += ([n, m, n, m, n*rate, 0, 0] if rate is not None else [n, m, 0, 0, 0, n, m])
    result = []
    for (family, sensitivity, threshold, scope), values in sorted(sums.items()):
        group = family, sensitivity, threshold, scope
        old, new = values
        result.append(dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope,
                           comparison_scope=comparison_scope, compared_chromosomes=",".join(sorted(chromosomes, key=int)),
                           counter_raw_sites=int(old[0]), targeted_raw_sites=int(new[0]), sites_delta=int(new[0]-old[0]),
                           counter_raw_multi=int(old[1]), targeted_raw_multi=int(new[1]), multi_delta=int(new[1]-old[1]),
                           counter_supported_sites=int(old[2]), targeted_supported_sites=int(new[2]),
                           counter_supported_multi=int(old[3]), targeted_supported_multi=int(new[3]),
                           counter_recomputed_expected=float(old[4]), targeted_expected=float(new[4]),
                           expected_delta=float(new[4]-old[4]), excluded_sites=int(new[5]), excluded_multi=int(new[6]),
                           compared_strata=strata_counts[group], mismatched_strata=mismatched[group],
                           status="matched" if not mismatched[group] else "mismatch"))
    label_rows = [dict(counter_label=label, exact_family=label[7:], selected_exact=label[7:] in families,
                       eligible_sites_ac5_all=count, alias_policy="exact family strings only; no alias collapsing")
                  for label, count in sorted(labels.items())]
    return result, mismatch_rows, label_rows


def _optional_number(value):
    if value is None or str(value).strip().lower() in {"", ".", "na", "none", "nan"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite original enrichment value: {value}")
    return result


def plot_spatial(partitions, profiles, out, end_bp, prefix="", model_label="Common-support model"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scope = "eur_ld_partition" if any(r["scope"] == "eur_ld_partition" for r in partitions) else "all_family_sites"
    primary = [r for r in partitions if r["scope"] == scope and r["sensitivity"] == "strict" and r["threshold"] == 5]
    families = sorted({r["family"] for r in primary})
    if not families:
        return
    fig, axes = plt.subplots(len(families), 2, figsize=(12, 3.8*len(families)), squeeze=False)
    categories = (*PARTITIONS, "total")
    names = [f"Head\n{end_bp} bp", "Body", f"Tail\n{end_bp} bp", "Short", "Ambiguous", "Total"]
    for i, family in enumerate(families):
        rr = {r["partition"]:r for r in primary if r["family"] == family}
        for j, part in enumerate(categories):
            r = rr.get(part)
            if not r or r["oe"] is None:
                continue
            color = "#176B87" if part != "total" else ".2"
            axes[i, 0].plot(j, r["oe"], "o", color=color,
                            markerfacecolor=color if r["status"] == "complete_descriptive" else "white")
            if r["ci_low"] is not None:
                axes[i, 0].plot([j, j], [r["ci_low"], r["ci_high"]], color="#176B87")
            axes[i, 1].bar(j, r["excess_contribution"], color="#D16B38" if part != "total" else ".3")
        axes[i, 0].axhline(1, ls="--", color=".6")
        axes[i, 1].axhline(0, color=".6", lw=.8)
        axes[i, 0].set(title=f"{family}: disjoint aligned-body parts", ylabel="Multi-ALT observed / expected")
        axes[i, 1].set(title=f"{family}: weighted contribution to total excess",
                       ylabel="Contribution to overall O/E − 1\n(Observed part − expected part) / expected total")
        if rr.get("total", {}).get("status") == "unavailable_primary_loco":
            for ax in axes[i]:
                ax.text(.5, .5, "Original total unavailable:\nexposed strata lack LOCO training", transform=ax.transAxes, ha="center")
        for ax in axes[i]:
            ax.set(xticks=range(len(categories)), xticklabels=names)
    fig.suptitle(f"{model_label}\nAC ≥ 5 · strict masks · {scope}\nShort ≤ {2*end_bp} aligned bp; physical-block CI; open markers: low support / few events or blocks")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"{prefix}spatial_partitions.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(len(families), 1, figsize=(9, 3.5*len(families)), squeeze=False)
    for i, family in enumerate(families):
        ax = axes[i, 0]
        for threshold, color in ((5, "#176B87"), (10, "#D16B38")):
            rr = sorted([r for r in profiles if r["family"] == family and r["scope"] == scope
                         and r["sensitivity"] == "strict" and r["threshold"] == threshold and r["body_bin"] >= 0],
                        key=lambda r:r["body_bin"])
            if not rr:
                continue
            x = [(r["start_fraction"]+r["end_fraction"])/2 for r in rr]
            y = [r["oe"] if r["oe"] is not None else np.nan for r in rr]
            low = [r["ci_low"] if r["ci_low"] is not None else np.nan for r in rr]
            high = [r["ci_high"] if r["ci_high"] is not None else np.nan for r in rr]
            ax.plot(x, y, "-", color=color, label=f"AC ≥ {threshold}")
            for xx, yy, r in zip(x, y, rr):
                ax.plot(xx, yy, "o", ms=3, color=color,
                        markerfacecolor=color if r["status"] == "complete_descriptive" else "white")
            ax.fill_between(x, low, high, color=color, alpha=.12)
        ax.axhline(1, ls="--", color=".6")
        ax.set(xlim=(0, 1), title=family, xlabel="Strand-oriented position along concatenated aligned bases (5′ → 3′)",
               ylabel="Multi-ALT observed / expected")
        ax.legend(frameon=False)
        if all(r["oe"] is None for r in profiles if r["family"] == family and r["scope"] == scope and r["sensitivity"] == "strict"):
            ax.text(.5, .5, "Original O/E unavailable or no estimable bins", transform=ax.transAxes, ha="center")
    fig.suptitle(f"{model_label}\nWhole joined-repeat profile · short fragments included · gaps excluded\nAmbiguous assignments retained in table; pointwise physical-block CI\nOpen markers: low support / few events or blocks")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"{prefix}whole_element_profiles.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def merge_spatial_tables(paths, destination):
    """Stream large element tables instead of collecting all records in RAM."""
    writer = None
    with open_text(destination, "wt") as target:
        for path in paths:
            with open_text(path) as source:
                reader = csv.DictReader(source, delimiter="\t")
                if reader.fieldnames == ["status"]:
                    continue
                if writer is None:
                    writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t")
                    writer.writeheader()
                elif reader.fieldnames != writer.fieldnames:
                    raise ValueError(f"Inconsistent spatial table fields: {path}")
                writer.writerows(reader)
        if writer is None:
            target.write("status\n")


@transactional("spatial_manifest.json")
def run_spatial(args):
    if args.end_bp < 1 or args.body_bins < 2 or args.block_bp < 100000:
        raise ValueError("Require end-bp > 0, body-bins >= 2 and physical blocks >= 100 kb")
    families = [family.strip() for family in args.families.split(",")]
    if not all(families) or len(set(families)) != len(families):
        raise ValueError("Families must be distinct nonempty exact RepeatMasker family strings")
    args.families = ",".join(families)
    run, annotation, tasks, source_chromosomes = targeted_inputs(args)
    options, out = run["options"], Path(args.out)
    if not {5, 10} <= {int(t) for t in options.get("thresholds", [])}:
        raise ValueError("The source architecture scan must contain AC5 and AC10 counters for spatial reconciliation")
    if not options.get("ld_block_dir"):
        raise ValueError("Primary spatial decomposition requires the EUR LD partition recorded in the source scan")
    files, used = [], {c for c, _ in tasks}
    for chrom, source in tasks:
        print(f"[spatial] indexed chr{chrom} aligned-body query", flush=True)
        files.append(scan_spatial_chrom(source, chrom, annotation/f"chr{chrom}.tsv.gz", out/"cache", options, args))
    with_eur = bool(options.get("ld_block_dir"))
    global_paths = architecture_paths(args.modern_dir, "global")
    background = Background(global_paths, args.min_training_sites)
    primary_background = PrimaryBackground(global_paths, families)
    data, raw = spatial_counts(files, background, with_eur)
    partitions, profiles, identities, contrasts = summarize_spatial(data, args)
    primary_data, _ = spatial_counts(files, primary_background, with_eur, collect_raw=False)
    comparison_scope = "all_22_autosomes" if used == AUTOSOMES else "selected_chromosomes_subset"
    reconciliation, mismatches, labels = reconcile_spatial_counters(args.modern_dir, raw, background, used,
                                                                families, with_eur, comparison_scope, primary_background)
    from herv_primary import scoped_counter_totals, reconcile_spatial_to_counters
    # All recorded global chromosomes train the primary rates; target totals
    # are restricted to the exact chromosomes actually queried by this stage.
    counter_totals = scoped_counter_totals(architecture_paths(args.modern_dir, "exposure"),
                                          primary_background, families, used)
    primary_parts, primary_profiles, primary_identities = summarize_primary_spatial(primary_data, args, counter_totals)
    primary_reconciliation, primary_failures = reconcile_spatial_to_counters(counter_totals, primary_parts)
    primary_checks = {(r["family"], r["sensitivity"], r["threshold"]): r for r in primary_reconciliation}
    for record in [*primary_parts, *primary_profiles, *primary_identities]:
        check = primary_checks.get((record["family"], record["sensitivity"], record["threshold"]), {})
        record.update(comparison_scope=comparison_scope, counter_reconciliation_status=check.get("status"))
    write_rows(out/"spatial_partitions.tsv", partitions)
    write_rows(out/"whole_element_profiles.tsv", profiles)
    write_rows(out/"spatial_identity.tsv", identities)
    write_rows(out/"spatial_contrasts.tsv", contrasts)
    write_rows(out/"primary_spatial_partitions.tsv", primary_parts)
    write_rows(out/"primary_whole_element_profiles.tsv", primary_profiles)
    write_rows(out/"primary_spatial_identity.tsv", primary_identities)
    write_rows(out/"primary_counter_reconciliation.tsv", primary_reconciliation)
    write_rows(out/"family_counter_reconciliation.tsv", reconciliation)
    write_rows(out/"family_counter_mismatches.tsv.gz", mismatches,
               fields=["chrom", "family", "sensitivity", "threshold", "scope", "context", "af_bin", "an_bin",
                       "counter_sites", "target_sites", "counter_multi", "target_multi"])
    write_rows(out/"family_counter_labels.tsv", labels)
    merge_spatial_tables([out/"cache"/f"chr{c}.annotation_audit.tsv" for c, _ in tasks], out/"spatial_annotation_audit.tsv")
    merge_spatial_tables([out/"cache"/f"chr{c}.elements.tsv.gz" for c, _ in tasks], out/"spatial_elements.tsv.gz")
    if not mismatches and not primary_failures and not getattr(args, "no_plots", False):
        label = "Primary full-scope estimator · LOCO min 1"
        if used != source_chromosomes:
            label += "\nSelected target chromosomes"
        plot_spatial(primary_parts, primary_profiles, out, args.end_bp, prefix="primary_", model_label=label)
        plot_spatial(partitions, profiles, out, args.end_bp,
                     model_label=f"Common-support sensitivity · autosomes · LOCO min {args.min_training_sites}")
    manifest = dict(schema=SPATIAL_SCHEMA, status="counter_count_mismatch" if mismatches else "primary_counter_mismatch" if primary_failures else "complete",
                    modern_dir=str(Path(args.modern_dir).resolve()),
                    families=families, exact_family_label_policy="LTR class and exact raw family strings; no alias collapsing",
                    chromosomes=sorted(used, key=int), available_autosomal_chromosomes=sorted(source_chromosomes, key=int),
                    comparison_scope=comparison_scope, full_available_autosomal_scope=used == source_chromosomes,
                    raw_counts_matched=not mismatches, mismatched_strata=len(mismatches),
                    primary_model="primary_full_scope", primary_counter_matches=not primary_failures,
                    primary_counter_mismatches=len(primary_failures), primary_minimum_training_sites=1,
                    primary_background_chromosomes=sorted({c for c, _ in primary_background.counts}, key=scan.chrom_sort_key),
                    primary_missing_rule="Retain raw observed and eligible counts; any unsupported target stratum makes the full-scope total and its weighted decomposition unavailable",
                    models_share_one_count_cache=True, additional_vcf_queries_for_primary_model=0,
                    primary_scope="eur_ld_partition" if with_eur else "all_family_sites",
                    end_bp=args.end_bp, body_bins=args.body_bins, physical_cluster_bp=args.block_bp,
                    denominator="Base-eligible SNV sites, fixed across AC5/AC10",
                    background="LOCO all recorded global chromosomes, minimum 1; canonical trinucleotide / primary ALT AF / AN; (multi+.5)/(sites+1)",
                    common_support_denominator="Unsupported observed and expected counts removed together; autosomal training background",
                    minimum_training_sites=args.min_training_sites, minimum_training_sites_applies_to="common_support_sensitivity",
                    minimum_support=args.min_support,
                    partition_definition={"head":f"first {args.end_bp} aligned bases in strand order",
                                          "tail":f"last {args.end_bp} aligned bases in strand order",
                                          "body":"remaining aligned bases of a long joined record",
                                          "short_element":f"entire record when aligned length <= {2*args.end_bp} bp",
                                          "ambiguous_overlap":"conflicting strand or part assignments; retained in total"},
                    normalized_position="base centres along concatenated merged aligned segments, reversed on minus strand; conflicting bins retained as -1",
                    accounting="overall O/E=sum(E_part/E_total*O_part/E_part); sum excess contributions=overall O/E-1; residual classes retained",
                    counter_reconciliation="N/M checked per chromosome/context/AF/AN/threshold/sensitivity; E recomputed from the same counts and training rule",
                    inference=f"{args.block_bp:g}bp cluster log-ratio sandwich pointwise CI; head/body and tail/body ratios; BH across strict AC5 primary-scope contrasts",
                    caveats=["Joined repeat records are not necessarily complete proviruses, true LTRs, insertion boundaries or TSSs",
                             "The -int suffix is a name audit, not a reconstructed structural assignment",
                             "Sites are counted once within each family; overlapping different families cannot be added",
                             "Normalized bins mix record lengths and are exploratory; short fragments remain separate in the fixed-bp partition",
                             "This is an SNV-denominator statistic, not a mutation rate per callable base",
                             "Background offsets are treated as fixed; local GC/recombination/repeat-age effects are not fitted",
                             "Aligned-family queries exclude external flanks; use the windows stage for flanks",
                             "Spatial physical-block CIs and primary LD-block bootstrap CIs use different uncertainty estimators"],
                    cache_files=[str(p.resolve()) for p in files])
    dump_json(out/"spatial_manifest.json", manifest)
    if mismatches:
        raise ValueError(f"Targeted spatial counts disagree with {len(mismatches)} source counter strata; see {out/'family_counter_mismatches.tsv.gz'}")
    if primary_failures:
        raise ValueError(f"Spatial enrichment disagrees with {len(primary_failures)} same-scope counter totals; see {out/'primary_counter_reconciliation.tsv'}")
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["summary", "spatial", "windows"])
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--erv-joined")
    p.add_argument("--fasta")
    p.add_argument("--families", default="ERV1,ERVK")
    p.add_argument("--end-bp", type=int, default=200, help="Spatial head/tail length in aligned bases; records <= twice this length form the short-element part")
    p.add_argument("--body-bins", type=int, default=20, help="Exploratory strand-oriented normalized aligned-body bins")
    p.add_argument("--flank", type=int, default=2000)
    p.add_argument("--window", type=int, default=200)
    p.add_argument("--step", type=int, default=50)
    p.add_argument("--block-bp", type=int, default=5000000)
    p.add_argument("--min-blocks", type=int, default=20)
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--min-support", type=float, default=.99)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--discovery-chromosomes", default="1,3,5,7,9,11,13,15,17,19,21")
    p.add_argument("--chromosomes")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--no-plots", action="store_true", help="Skip spatial plots")
    for name in ("segdup", "lcr", "low-mappability"):
        p.add_argument("--"+name+"-bed")
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    if not 0 < a.min_support <= 1 or a.min_training_sites < 1:
        raise ValueError("Require 0 < min-support <=1 and positive training count")
    if a.command == "summary":
        return summarize_common_support(a.modern_dir, a.out, a.min_training_sites, a.min_support, a.bootstrap, a.seed)
    if not a.erv_joined:
        raise ValueError("--erv-joined required for targeted spatial/boundary scan")
    if a.command == "spatial":
        return run_spatial(a)
    return run_windows(a)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit(f"herv-modern: {exc}")

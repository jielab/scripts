#!/usr/bin/env python3
"""Map measured nucleotide frequencies onto frozen PC1/PC2/age coordinates.

This is an allele occupancy atlas, not a new PCA or an allele-age estimator.
Use a completed allele3 ``time`` directory with its original PGEN, or supply
group-sorted long observations with actual allele counts and ploidy. Unassayed
bases and missing genotypes are never filled with zero. X/Y are supported by
the long-observation route, which must supply their actual per-person ploidy.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import sys

import numpy as np
from analysis_integrity import atomic_json, identity, transactional, validate_manifest
from pca_contract import historical_identity, read_scores

VERSION = "allele-space-1.0"
NA = {"", ".", "NA", "nan", "None"}
CELL_FIELDS = ["site_id", "chrom", "pos", "genome_build", "allele", "role", "evidence_kind",
    "counting_scope", "cohort", "time_bin", "start_bp", "end_bp", "pc1_bin", "pc2_bin",
    "frame_PC1", "frame_PC2", "frame_age_bp", "n_frame_samples", "n_declared_assayed",
    "n_called", "AC", "AN", "AF", "n_hardcall_carriers", "call_rate_among_declared",
    "frame_measurement_fraction", "usable_cell", "status"]
SUMMARY_FIELDS = ["site_id", "chrom", "pos", "genome_build", "allele", "role", "evidence_kind",
    "counting_scope", "cohort", "n_frame_samples", "n_declared_assayed", "n_called", "AC", "AN", "AF",
    "n_hardcall_carriers", "copy_weighted_PC1", "copy_weighted_PC2", "copy_weighted_age_bp",
    "oldest_observed_carrier_age_bp", "oldest_observed_carrier_ids", "oldest_observed_carrier_regions",
    "older_called_AN", "earliest_supported_time_bin", "earliest_supported_start_bp",
    "earliest_supported_end_bp", "earliest_supported_AC", "earliest_supported_AN",
    "earliest_supported_n_called", "earliest_supported_n_carriers", "earliest_supported_regions",
    "peak_AF_cell", "peak_AF", "peak_AF_AN", "n_usable_cells", "origin_inference"]
LIMITATIONS = [
    "PC1/PC2 are frozen genotype coordinates. Age BP is supplied metadata, not a third principal component.",
    "Agreement with a supplied date table validates joins and transcription, not independent genetic dating.",
    "This atlas describes the supplied, callable, ascertained panel. Its frequencies are not unbiased estimates for all ancient people.",
    "An unrepresented base is unavailable, never zero. Each nucleotide retains its own measured AN; sums of split-ALT counts do not establish tri-allelism.",
    "Oldest observed or encoded carriers do not date the origin of a mutation; imputed hardcalls are not independent ancient read evidence.",
    "Carrier locations and copy-weighted centroids depend on sampling. Cell AF corrects sample-count denominators but not ancestry, assay or archaeological ascertainment.",
    "No selection coefficient, balancing-selection test, or causal effect of introgression is estimated.",
    "X/Y long-table counts use supplied ploidy. No independence model, autosomal LD bootstrap, or significance claim is imposed on them.",
]


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate TSV headers: {path}")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed TSV row in {path}")
            yield row


def write_tsv(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "NA") for k in fields})


def finite(value):
    if value is None or str(value) in NA: return None
    result = float(value)
    return result if math.isfinite(result) else None


def canonical_chrom(value):
    value = str(value).removeprefix("chr")
    if value not in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        raise ValueError("Expected chromosome 1..22, X or Y")
    return value


def load_metadata(path, coordinates, maximum):
    rows = list(read_tsv(path)); seen = set()
    individuals = Counter(r.get("individual_id") for r in rows if r.get("individual_id") not in NA | {None})
    for row in rows:
        sid = row.get("sample_id")
        if not sid or sid in seen: raise ValueError("Missing/duplicate metadata sample_id")
        seen.add(sid)
        age = finite(row.get("age_bp")); cohort = row.get("cohort")
        reasons = []
        if row.get("include") != "1": reasons.append("include_not_1")
        if row.get("unrelated") != "1": reasons.append("unrelated_not_1")
        if sid not in coordinates: reasons.append("coordinates_missing")
        if cohort not in {"ancient", "modern"}: reasons.append("cohort_unknown")
        if cohort == "ancient" and (age is None or not 0 <= age <= maximum): reasons.append("outside_ancient_time_scope")
        if individuals.get(row.get("individual_id"), 0) > 1: reasons.append("duplicate_biological_individual")
        row.update(age_bp=age, eligible=not reasons, exclusion=";".join(reasons),
                   PC1=coordinates.get(sid, [None, None])[0], PC2=coordinates.get(sid, [None, None])[1])
    if not any(r["eligible"] for r in rows): raise ValueError("No eligible samples in the frozen coordinate frame")
    return rows


def date_audit(samples, reference, tolerance):
    if not reference:
        return [], {"status": "not_provided", "interpretation": "No external metadata concordance check was performed"}
    refs = {}
    for row in read_tsv(reference):
        sid = row.get("sample_id")
        if not sid or sid in refs: raise ValueError("Missing/duplicate reference metadata sample_id")
        refs[sid] = row
    result = []
    for sample in samples:
        if sample.get("cohort") != "ancient": continue
        sid = sample["sample_id"]; ref = refs.get(sid)
        observed = sample["age_bp"]; expected = finite(ref.get("age_bp")) if ref else None
        status = ("reference_missing" if ref is None else "date_missing" if observed is None or expected is None
                  else "concordant" if abs(observed - expected) <= tolerance else "date_mismatch")
        fields = []
        if ref:
            for key in ("cohort", "region"):
                if sample.get(key) not in NA | {None} and ref.get(key) not in NA | {None} and sample[key] != ref[key]:
                    fields.append(key)
        if fields: status = "metadata_mismatch"
        result.append(dict(sample_id=sid, age_bp=observed, reference_age_bp=expected,
            difference_years=observed-expected if observed is not None and expected is not None else None,
            mismatched_fields=";".join(fields), status=status))
    compared = sum(r["difference_years"] is not None for r in result)
    return result, {"status": "metadata_concordance_only" if compared else "no_comparable_dates",
        "n_ancient_metadata": len(result), "n_compared": compared,
        "counts": dict(Counter(r["status"] for r in result)),
        "reference": identity(reference, full=True), "tolerance_years": tolerance,
        "interpretation": "Checks ID/date/region transcription; not independent genetic dating or validation by plotting age on Z"}


def build_cells(samples, bins, width, maximum):
    included = [i for i, s in enumerate(samples) if s["eligible"]]
    matrix = np.asarray([[samples[i]["PC1"], samples[i]["PC2"]] for i in included], dtype=float)
    lower, upper = matrix.min(axis=0), matrix.max(axis=0)
    groups = defaultdict(list)
    for i in included:
        s = samples[i]
        xy = [0 if upper[j] == lower[j] else min(bins-1, max(0, int((s[f"PC{j+1}"]-lower[j])/(upper[j]-lower[j])*bins))) for j in (0, 1)]
        epoch = "modern" if s["cohort"] == "modern" else min(int(s["age_bp"]//width), math.ceil(maximum/width)-1)
        groups[(s["cohort"], epoch, *xy)].append(i)
    cells = []
    for (cohort, epoch, x, y), ids in sorted(groups.items(), key=lambda kv: str(kv[0])):
        cells.append(dict(cohort=cohort, time_bin=epoch, pc1_bin=x, pc2_bin=y,
            start_bp=None if cohort == "modern" else epoch*width,
            end_bp=None if cohort == "modern" else min(maximum, (epoch+1)*width),
            frame_PC1=float(np.mean([samples[i]["PC1"] for i in ids])),
            frame_PC2=float(np.mean([samples[i]["PC2"] for i in ids])),
            frame_age_bp=None if cohort == "modern" else float(np.mean([samples[i]["age_bp"] for i in ids])),
            n_frame_samples=len(ids), members=np.asarray(ids, dtype=int)))
    return cells, dict(pc_bins=bins, bin_years=width, max_age_bp=maximum,
        PC1_range=[float(lower[0]), float(upper[0])], PC2_range=[float(lower[1]), float(upper[1])],
        definition="Equal-width PC bins fixed across alleles; occupied metadata cells only; oldest time endpoint inclusive")


def observation_groups(path, samples, build, kind):
    lookup = {s["sample_id"]: i for i, s in enumerate(samples)}
    previous = None; counts = ploidies = declared = info = None
    for row in read_tsv(path):
        required = {"site_id", "chrom", "pos", "allele", "sample_id", "count", "ploidy"}
        if not required <= row.keys(): raise ValueError("Observation TSV requires " + ",".join(sorted(required)))
        key = (row["site_id"], row["allele"])
        if key[0] in NA or key[1] not in "ACGT" or len(key[1]) != 1: raise ValueError("Invalid site_id/nucleotide")
        if previous is not None and key < previous: raise ValueError("Observations must be sorted by site_id, allele; a group cannot recur")
        chrom, pos = canonical_chrom(row["chrom"]), int(row["pos"])
        if pos < 1: raise ValueError("Position must be 1-based and positive")
        if previous and key[0] == previous[0] and (chrom, pos) != (info["chrom"], info["pos"]):
            raise ValueError("One site_id maps to multiple genomic positions across alleles")
        if key != previous:
            if previous is not None: yield info, counts, ploidies, declared
            previous = key
            counts = np.full(len(samples), np.nan); ploidies = np.full(len(samples), np.nan); declared = np.zeros(len(samples), dtype=bool)
            info = dict(site_id=key[0], chrom=chrom, pos=pos, genome_build=build, allele=key[1],
                role=row.get("role", "unspecified"), evidence_kind=kind, counting_scope="explicit_count_of_one_nucleotide")
        elif chrom != info["chrom"] or pos != info["pos"]:
            raise ValueError("One site_id/allele maps to multiple genomic positions")
        sid = row["sample_id"]
        if sid not in lookup: raise ValueError("Observation sample absent from metadata: " + sid)
        i = lookup[sid]
        if declared[i]: raise ValueError("Duplicate sample within an allele observation group: " + sid)
        count, ploidy = finite(row["count"]), finite(row["ploidy"])
        if ploidy not in (1., 2.): raise ValueError("Every declared observation needs actual ploidy 1 or 2")
        if count is not None and (not 0 <= count <= ploidy or (kind != "posterior_dosage" and count != int(count))):
            raise ValueError("Count is outside ploidy or noninteger for a hardcall")
        if str(row["count"]) not in NA and count is None: raise ValueError("Nonfinite allele count")
        declared[i] = True; ploidies[i] = ploidy; counts[i] = np.nan if count is None else count
    if previous is not None: yield info, counts, ploidies, declared


def pgen_groups(directory, prefix, samples, chunk_size):
    from allele_measurement import iter_pgen_chunks, panel_pvar
    directory = Path(directory)
    manifest = validate_manifest(directory, require_evidence=True)
    recorded = manifest.get("inputs", {})
    paths = [Path(str(prefix)+suffix).resolve() for suffix in (".pgen", ".psam")]
    paths.append(panel_pvar(prefix).resolve())
    # The saved plans index the exact original PVAR and sample order. A new panel
    # with similar filenames or dimensions is not an acceptable replacement.
    for path in paths:
        saved = recorded.get(str(path))
        if not saved: raise ValueError(f"Time-stage manifest does not bind original input {path}")
        current = identity(path, full="sha256" in saved)
        keys = ("size", "sha256") if "sha256" in saved else ("size", "mtime_ns")
        if any(current.get(k) != saved.get(k) for k in keys): raise ValueError(f"Original panel changed: {path}")
    with Path(str(prefix)+".psam").open() as handle:
        lines = (x for x in handle if x.strip() and not x.startswith("##"))
        header = next(lines).lstrip("#").split()
        if "IID" not in header: raise ValueError("PSAM lacks IID")
        ix = header.index("IID"); ids = [line.split()[ix] for line in lines]
    if len(ids) != len(set(ids)): raise ValueError("Duplicate PSAM IID")
    lookup = {s["sample_id"]: i for i, s in enumerate(samples)}
    missing = set(ids) - lookup.keys()
    if missing: raise ValueError(f"{len(missing)} PSAM samples lack metadata; preserve the complete sample audit")
    target_indices = np.asarray([lookup[sid] for sid in ids], dtype=int)
    plans = json.loads((directory/"measurement_plans.json").read_text())
    if plans.get("schema") != "nucleotide-count-1": raise ValueError("Unsupported measurement plan schema")
    entries = plans.get("ALT_measurements", [])
    if len({e["allele_id"] for e in entries}) != len(entries): raise ValueError("Duplicate allele measurement plan")
    for start, values, _ in iter_pgen_chunks(str(prefix)+".pgen", len(ids), entries, chunk_size, pvar_path=panel_pvar(prefix)):
        for j, entry in enumerate(entries[start:start+len(values)]):
            if canonical_chrom(entry["chrom"]) in {"X", "Y"}:
                raise ValueError("Sex-chromosome PGEN plans lack actual ploidy; use --observations with explicit ploidy")
            counts = np.full(len(samples), np.nan); ploidies = np.full(len(samples), np.nan); declared = np.zeros(len(samples), dtype=bool)
            counts[target_indices] = np.where(values[j] >= 0, values[j], np.nan)
            ploidies[target_indices] = 2.; declared[target_indices] = True
            info = dict(site_id=entry["site_id"], chrom=entry["chrom"], pos=entry["pos"], genome_build="GRCh38",
                allele=entry["allele"], role=entry.get("role", "unspecified"), evidence_kind="imputed_hardcall",
                counting_scope=entry["counting_scope"])
            yield info, counts, ploidies, declared


def summarize(info, counts, ploidies, declared, samples, cells, minimum, minimum_carriers):
    hardcall = info["evidence_kind"] != "posterior_dosage"
    cell_rows = []
    for cell in cells:
        ids = cell["members"]; called = np.isfinite(counts[ids])
        ac = float(np.sum(counts[ids][called])); an = float(np.sum(ploidies[ids][called]))
        n_called = int(called.sum()); n_declared = int(declared[ids].sum())
        usable = n_called >= minimum
        row = dict(info, **{k: v for k, v in cell.items() if k != "members"}, n_declared_assayed=n_declared,
            n_called=n_called, AC=ac if an else None, AN=an, AF=ac/an if an else None,
            n_hardcall_carriers=int(np.sum(counts[ids][called] > 0)) if hardcall else None,
            call_rate_among_declared=n_called/n_declared if n_declared else None,
            frame_measurement_fraction=n_called/len(ids), usable_cell=int(usable),
            status="descriptive_frequency" if usable else "low_support" if n_called else "unmeasured")
        cell_rows.append(row)
    summaries = []
    for cohort in ("ancient", "modern"):
        ids = np.asarray([i for i, s in enumerate(samples) if s["eligible"] and s["cohort"] == cohort], dtype=int)
        if not len(ids): continue
        called_ids = ids[np.isfinite(counts[ids])]; carrier_ids = called_ids[counts[called_ids] > 0] if hardcall else np.array([], dtype=int)
        ac = float(np.sum(counts[called_ids])); an = float(np.sum(ploidies[called_ids]))
        row = dict(info, cohort=cohort, n_frame_samples=len(ids), n_declared_assayed=int(declared[ids].sum()), n_called=len(called_ids),
            AC=ac if an else None, AN=an, AF=ac/an if an else None,
            n_hardcall_carriers=len(carrier_ids) if hardcall else None, origin_inference="not_estimated")
        for field in ("PC1", "PC2", "age_bp"):
            row["copy_weighted_"+field] = (float(np.average([samples[i][field] for i in called_ids], weights=counts[called_ids]))
                if ac > 0 and (field != "age_bp" or cohort == "ancient") else None)
        if cohort == "ancient" and len(carrier_ids):
            oldest = max(samples[i]["age_bp"] for i in carrier_ids)
            oldest_ids = [i for i in carrier_ids if samples[i]["age_bp"] == oldest]
            row.update(oldest_observed_carrier_age_bp=oldest,
                oldest_observed_carrier_ids=";".join(samples[i]["sample_id"] for i in oldest_ids),
                oldest_observed_carrier_regions=";".join(sorted({samples[i].get("region", "unknown") for i in oldest_ids})),
                older_called_AN=float(sum(ploidies[i] for i in called_ids if samples[i]["age_bp"] > oldest)))
            # A supported first-observed interval is pooled only within the SAME
            # time bin; raw oldest individuals remain separately available.
            grouped = defaultdict(list)
            for c in cell_rows:
                if c["cohort"] == "ancient": grouped[c["time_bin"]].append(c)
            supported = [(e, cc) for e, cc in grouped.items()
                if sum(c["n_called"] for c in cc) >= minimum
                and sum(c["n_hardcall_carriers"] for c in cc) >= minimum_carriers]
            if supported:
                epoch, selected = max(supported, key=lambda x: x[0]); lo, hi = selected[0]["start_bp"], selected[0]["end_bp"]
                in_epoch = {int(i) for c in cells if c["cohort"] == "ancient" and c["time_bin"] == epoch for i in c["members"]}
                members = [i for i in carrier_ids if i in in_epoch]
                row.update(earliest_supported_time_bin=epoch, earliest_supported_start_bp=lo, earliest_supported_end_bp=hi,
                    earliest_supported_AC=sum(c["AC"] or 0 for c in selected), earliest_supported_AN=sum(c["AN"] for c in selected),
                    earliest_supported_n_called=sum(c["n_called"] for c in selected),
                    earliest_supported_n_carriers=sum(c["n_hardcall_carriers"] for c in selected),
                    earliest_supported_regions=";".join(sorted({samples[i].get("region", "unknown") for i in members})))
        valid = [c for c in cell_rows if c["cohort"] == cohort and c["usable_cell"] and c["AC"] > 0]
        row["n_usable_cells"] = sum(c["cohort"] == cohort and c["usable_cell"] for c in cell_rows)
        if valid:
            peak = max(valid, key=lambda c: (c["AF"], c["AN"], str(c["time_bin"]), -c["pc1_bin"], -c["pc2_bin"]))
            row.update(peak_AF_cell=f"{peak['time_bin']}:{peak['pc1_bin']}:{peak['pc2_bin']}", peak_AF=peak["AF"], peak_AF_AN=peak["AN"])
        summaries.append(row)
    return cell_rows, summaries


def plot_examples(examples, out, grid):
    if not examples: return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    columns = min(3, len(examples)); rows = math.ceil(len(examples)/columns)
    fig = plt.figure(figsize=(5*columns, 4.3*rows))
    for i, (info, cells) in enumerate(examples):
        ax = fig.add_subplot(rows, columns, i+1, projection="3d")
        measured = [c for c in cells if c["cohort"] == "ancient" and c["usable_cell"]]
        if measured:
            ax.scatter([c["frame_PC1"] for c in measured], [c["frame_PC2"] for c in measured],
                [c["frame_age_bp"]/1000 for c in measured], c=[c["AF"] for c in measured],
                s=40, vmin=0, vmax=1, cmap="viridis", alpha=.85)
        evidence = {"direct_gt": "direct genotype", "imputed_hardcall": "imputed hardcall", "posterior_dosage": "posterior dosage"}[info["evidence_kind"]]
        ax.set(xlabel="Frozen PC1", ylabel="Frozen PC2", zlabel="Age (ka BP)",
            title=f"{info['chrom']}:{info['pos']} · {info['allele']}\n{evidence}")
        for axis, values in ((ax.set_xlim, grid["PC1_range"]), (ax.set_ylim, grid["PC2_range"])):
            low, high = values; padding = .05*(high-low) if high > low else .5
            axis(low-padding, high+padding)
        ax.set_zlim(grid["max_age_bp"]/1000, 0)
    fig.suptitle("Measured allele frequency in fixed PC/time cells\nMissing cells are omitted; no mutation origin inferred", y=.99)
    fig.subplots_adjust(top=.82 if rows == 1 else .9, right=.89, wspace=.22, hspace=.3)
    scalar = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, 1)); scalar.set_array([])
    color_axis = fig.add_axes([.93, .25, .012, .4])
    fig.colorbar(scalar, cax=color_axis, label="AC / actual callable AN")
    for suffix in ("png", "pdf"): fig.savefig(out/f"allele_space_examples.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


@transactional("manifest.json", input_dirs=("alleles_dir", "basis_dir"))
def run(args):
    if args.pc_bins < 1 or args.bin_years < 1 or args.max_age_bp < args.bin_years or args.min_samples < 1:
        raise ValueError("Invalid cell/sample parameters")
    if args.min_earliest_carriers < 1 or args.chunk_size < 1 or args.plot_alleles < 0 or args.date_tolerance_years < 0:
        raise ValueError("Invalid earliest-observation, chunk or plotting parameters")
    if args.observations and (not args.genome_build or not args.evidence_kind):
        raise ValueError("--observations requires --genome-build and --evidence-kind")
    if args.alleles_dir and not args.pgen_prefix: raise ValueError("--time-dir requires the original --pgen-prefix")
    if args.alleles_dir and (args.genome_build not in (None, "GRCh38") or args.evidence_kind not in (None, "imputed_hardcall")):
        raise ValueError("Saved time-stage plans use GRCh38 nucleotide labels and imputed hardcalls; do not override this evidence identity")
    coords, columns, contract = read_scores(args.coordinates, 2, basis_dir=args.basis_dir,
        receipt=args.coordinate_receipt, require_verified=not args.allow_unverified_coordinates)
    provenance = historical_identity(args.metadata, args.provenance_manifest, contract.get("metadata_recorded"))
    verified = contract["status"] == "verified" and provenance["historical_metadata_verification"] == "verified"
    if not verified and not args.allow_unverified_coordinates:
        raise ValueError("Metadata/coordinates lack verified historical identity; only explicit --allow-unverified-coordinates permits descriptive output")
    samples = load_metadata(args.metadata, coords, args.max_age_bp)
    out = Path(args.out)
    write_tsv(out/"sample_audit.tsv.gz", samples, ["sample_id", "individual_id", "cohort", "age_bp", "region", "include", "unrelated", "PC1", "PC2", "eligible", "exclusion"])
    date_rows, dates = date_audit(samples, args.reference_metadata, args.date_tolerance_years)
    write_tsv(out/"date_metadata_audit.tsv", date_rows, ["sample_id", "age_bp", "reference_age_bp", "difference_years", "mismatched_fields", "status"])
    if args.date_mismatch_policy == "fail" and any(r["status"] in {"date_mismatch", "metadata_mismatch"} for r in date_rows):
        raise ValueError("Date/reference metadata mismatch; see staged date_metadata_audit.tsv diagnostics")
    if args.date_mismatch_policy == "fail" and args.reference_metadata and dates["n_compared"] == 0:
        raise ValueError("Reference metadata contains no comparable ancient dates")
    cells, grid = build_cells(samples, args.pc_bins, args.bin_years, args.max_age_bp)
    write_tsv(out/"coordinate_cells.tsv", cells, ["cohort", "time_bin", "start_bp", "end_bp", "pc1_bin", "pc2_bin", "frame_PC1", "frame_PC2", "frame_age_bp", "n_frame_samples"])
    stream = (observation_groups(args.observations, samples, args.genome_build, args.evidence_kind) if args.observations
              else pgen_groups(args.alleles_dir, args.pgen_prefix, samples, args.chunk_size))
    examples = []; n_alleles = 0
    with gzip.open(out/"allele_cells.tsv.gz", "wt", newline="") as handle, \
            gzip.open(out/"allele_space_summary.tsv.gz", "wt", newline="") as summary_handle:
        writer = csv.DictWriter(handle, CELL_FIELDS, delimiter="\t", extrasaction="ignore", lineterminator="\n"); writer.writeheader()
        summary_writer = csv.DictWriter(summary_handle, SUMMARY_FIELDS, delimiter="\t", extrasaction="ignore", lineterminator="\n"); summary_writer.writeheader()
        for info, counts, ploidies, declared in stream:
            cell_rows, summary = summarize(info, counts, ploidies, declared, samples, cells, args.min_samples, args.min_earliest_carriers)
            for row in cell_rows: writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "NA") for k in CELL_FIELDS})
            for row in summary: summary_writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "NA") for k in SUMMARY_FIELDS})
            n_alleles += 1
            if len(examples) < args.plot_alleles and any(c["cohort"] == "ancient" and c["usable_cell"] and c["AC"] > 0 for c in cell_rows):
                examples.append((info, cell_rows))
    if not n_alleles: raise ValueError("No represented allele observations; refusing a falsely complete empty atlas")
    if not args.no_plots: plot_examples(examples, out, grid)
    result = dict(status="complete" if verified else "complete_descriptive_unverified", version=VERSION,
        n_alleles=n_alleles, n_metadata=len(samples), n_eligible=sum(s["eligible"] for s in samples),
        coordinate_contract=contract, metadata_provenance=provenance, grid=grid, date_audit=dates,
        genotype_buffer="One allele by metadata sample; PGEN decoding uses the selected chunk size",
        inference="Descriptive occupancy and measured frequencies only", limitations=LIMITATIONS, parameters=vars(args))
    atomic_json(out/"manifest.json", result)
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--time-dir", "--alleles-dir", dest="alleles_dir",
        help="Completed allele3 time stage containing measurement_plans.json and manifest.json (not the alleles catalog stage)")
    source.add_argument("--observations", help="TSV sorted by site_id,allele; columns site_id chrom pos allele sample_id count ploidy")
    for name in ("metadata", "coordinates", "out"): p.add_argument("--"+name, required=True)
    for name in ("pgen-prefix", "basis-dir", "coordinate-receipt", "provenance-manifest", "reference-metadata"):
        p.add_argument("--"+name)
    p.add_argument("--genome-build", choices=("GRCh37", "GRCh38"))
    p.add_argument("--evidence-kind", choices=("direct_gt", "imputed_hardcall", "posterior_dosage"))
    p.add_argument("--allow-unverified-coordinates", action="store_true", help="Explicit descriptive-only run without trusted coordinate/metadata receipts")
    for name, default in (("pc-bins", 6), ("bin-years", 1000), ("max-age-bp", 10000), ("min-samples", 20),
            ("min-earliest-carriers", 2), ("chunk-size", 32), ("plot-alleles", 6)):
        p.add_argument("--"+name, type=int, default=default)
    p.add_argument("--date-tolerance-years", type=float, default=1.)
    p.add_argument("--date-mismatch-policy", choices=("fail", "report"), default="fail")
    p.add_argument("--no-plots", action="store_true")
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
        print(json.dumps({k: result[k] for k in ("status", "n_alleles", "n_eligible")}, sort_keys=True))
        return 0
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"allele-space: {exc}")


if __name__ == "__main__":
    raise SystemExit(main())

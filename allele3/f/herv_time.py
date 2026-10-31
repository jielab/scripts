#!/usr/bin/env python3
"""HERV candidate trajectories in a fixed ancient panel; direct-read audit is separate.

The default input is the Reich imputed, biallelic GRCh37 PGEN. Its absence of an
extra base is NOT an observed zero. This module never infers HERV insertions,
first appearances, or selection coefficients from this panel.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import contextlib
import csv
import gzip
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import sys

if "allele_common" not in sys.modules:
    spec = importlib.util.spec_from_file_location("allele_common", Path(__file__).with_name("0.common.py"))
    common = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = common
    spec.loader.exec_module(common)
from allele_common import GenomeMapper, MappingError, canonical_chrom, orient, dump_json
from allele_measurement import (index_target_records, catalog_measurements, panel_pvar,
                                select_primary_measurements, iter_primary_counts)

VERSION = "herv-time-2.0-nucleotide-counts"
BASES = set("ACGT")
NA = {"", ".", "NA", "nan", "None"}
LIMITATIONS = [
    "Targets are ascertained in present-day survivors; trajectories do not estimate ancient multi-allelic discovery or extinction rates.",
    "The imputed biallelic reference panel cannot independently demonstrate ancient persistence of an unrepresented third base.",
    "Hardcall heterozygosity denotes exactly one copy of the counted nucleotide (target versus all other bases), not full-locus heterozygosity, HERV insertion burden or genome-wide diversity.",
    "The prespecified modern primary nucleotide is counted; unavailable primary nucleotides have no alternate-role fallback. External targets without primary_alt use an explicitly labelled deterministic fallback.",
    "Encoded ALT counts remain usable after build reference changes. Native multiallelic PGEN uses its actual allele index, including REF index 0. A biallelic/split row cannot prove its REF nucleotide count; no single-ALT complement is relabelled REF.",
    "PCA, geographic, coverage and data-source adjustment is descriptive; time/ancestry confounding and imputation-reference bias remain.",
    "A time slope is tested two-sided. No increasing direction, peak, or significance threshold is used to choose targets.",
    "Bootstrap intervals resample unrelated people conditional on this fixed matched locus panel; they do not generalize over new HERV loci.",
    "Dates are point estimates in years BP before 1950; dating uncertainty is not propagated.",
]
RUN_OUTPUTS = (
    "mapped_targets.tsv", "allele_availability.tsv", "selected_variants.tsv", "primary_contrast_audit.tsv", "site_availability.tsv", "measurement_plans.json",
    "contrast_trajectories.tsv.gz", "fixed_common_panel.tsv", "sample_bin_composition.tsv", "sample_contrast_burden.tsv",
    "matched_panel_audit.tsv", "fixed_panel_time_summary.tsv", "time_trend.tsv",
    "ancient_fixed_panel_time.png", "ancient_fixed_panel_time.pdf", "ancient_shared_pca_time.png", "ancient_shared_pca_time.pdf",
)


def reset_run_outputs(args, out):
    """Invalidate only this command's generated products, never input files."""
    generated = {out / name for name in (*RUN_OUTPUTS, "manifest.json")}
    inputs = [getattr(args, key, None) for key in ("targets", "metadata", "coordinates", "hq_variants", "fasta38", "fasta37", "chain38to37", "chain37to38")]
    inputs += [str(args.pgen_prefix) + suffix for suffix in (".pgen", ".pvar", ".pvar.zst", ".psam")]
    collisions = {Path(p).resolve() for p in inputs if p} & {p.resolve() for p in generated}
    if collisions:
        raise ValueError("Inputs cannot also be generated time outputs; use a separate input path: " + ", ".join(map(str, sorted(collisions))))
    dump_json(out / "manifest.json", dict(version=VERSION, status="running", note="Previous generated tables and figures invalidated for this run"))
    for name in RUN_OUTPUTS:
        (out / name).unlink(missing_ok=True)


@contextlib.contextmanager
def text_input(path):
    path = Path(path)
    if str(path).endswith(".zst"):
        import zstandard
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as stream:
            with io.TextIOWrapper(stream) as handle:
                yield handle
    elif str(path).endswith((".gz", ".bgz")):
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open() as handle:
            yield handle


def rows(path, required=(), whitespace=False):
    with text_input(path) as handle:
        names = None
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"):
                continue
            cells = line.strip().split() if whitespace else line.rstrip("\n\r").split("\t")
            if names is None:
                names = [x.lstrip("#") for x in cells]
                if len(set(names)) != len(names) or not set(required) <= set(names):
                    raise ValueError(f"{path}: unique header needs {required}")
                continue
            if len(cells) != len(names):
                raise ValueError(f"{path}:{line_no}: malformed row")
            yield dict(zip(names, cells))
        if names is None:
            raise ValueError(f"Empty table: {path}")


def write_rows(path, records, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(temporary, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({k: "NA" if record.get(k) is None else record.get(k, ".") for k in fields})
    os.replace(temporary, path)


def finite(value):
    try:
        value = float(value)
    except (ValueError, TypeError):
        return None
    return value if math.isfinite(value) else None


def load_targets(path, build="38"):
    output, seen = [], set()
    for raw in rows(path, ("chrom", "pos", "ref", "alts")):
        row = dict(raw)
        if build == "37" and "chrom37" in raw:
            if raw.get("mapping_status", "mapped") != "mapped":
                continue
            row.update(chrom=raw["chrom37"], pos=raw["pos37"], ref=raw["ref37"], alts=raw["alts37"])
            if raw.get("primary_alt37", ".") not in NA:
                row["anchor_alt"] = raw["primary_alt37"]
        elif raw.get("primary_alt38", ".") not in NA and raw.get("reference_changed", "0") == "0":
            row["anchor_alt"] = raw["primary_alt38"]
        row["chrom"] = canonical_chrom(row["chrom"])
        row["pos"] = int(row["pos"])
        ref, alts = row["ref"].upper(), row["alts"].upper().split(",")
        if row["chrom"] not in {str(x) for x in range(1, 23)} or row["pos"] < 1:
            raise ValueError("Targets must be autosomal, positive, one-based positions")
        if ref not in BASES or not alts or any(x not in BASES or x == ref for x in alts) or len(set(alts)) != len(alts):
            raise ValueError(f"Invalid canonical SNV alleles: {row}")
        key = row["chrom"], row["pos"]
        if key in seen:
            raise ValueError(f"Duplicate target position: {key}; use one row per position")
        seen.add(key)
        row.update(ref=ref, alts=",".join(sorted(alts)))
        row.setdefault("site_id", f"{row['chrom']}:{row['pos']}:{ref}")
        row.setdefault("family", "unclassified")
        row.setdefault("group", "case")
        row.setdefault("match_id", row.get("match_set", row["site_id"]))
        if row["site_id"] in NA or row["match_id"] in NA:
            raise ValueError("site_id and match_id must be nonmissing")
        output.append(row)
    if len({r["site_id"] for r in output}) != len(output):
        raise ValueError("Duplicate site_id")
    return output


def map_targets(targets, mapper):
    output = []
    for target in targets:
        row = dict(target)
        try:
            mapped = mapper.map(row["chrom"], row["pos"], row["ref"], "GRCh38")
            observed_labels37 = {orient(base, mapped["strand37to38"]) for base in [row["ref"], *row["alts"].split(",")]}
            # A reference change does not prevent counting an encoded nucleotide.
            # Availability is assessed independently for each oriented modern base.
            row.update(mapped)
            row["site_id"] = target["site_id"]
            row.update(mapping_status="mapped", alts37=",".join(sorted(observed_labels37 - {mapped["ref37"]})))
        except MappingError as exc:
            row["mapping_status"] = getattr(exc, "reason", str(exc))
        output.append(row)
    return output


def index_panel(targets, pvar, high_quality):
    """Preserve native allele indices and every same-position companion row."""
    return index_target_records(targets, pvar, high_quality, rows=rows, canonical_chrom=canonical_chrom)


def target_measurements(targets, availability):
    def role(target, base):
        return "primary" if base == target.get("primary_alt") else "other_ALT"
    selected, audit, _ = catalog_measurements(targets, availability, orient=orient, alt_role=role)
    # Original output aliases are retained; explicit panel/count fields disambiguate
    # native multiallelic rows and counted bases that became the panel reference.
    for entry in audit:
        entry.update(alt=entry.get("panel_alt37"), base38=entry["allele"], status=entry["panel_status"])
    return selected, audit


def select_primary(availability, targets):
    """Select the exact modern primary nucleotide using a safe measurement plan."""
    selected, audit = target_measurements(targets, availability)
    return select_primary_measurements(targets, selected, audit)


def choose_primary(availability, targets):
    return select_primary(availability, targets)[0]


def metadata_in_order(path, sample_ids, bin_years):
    table = {}
    for raw in rows(path, ("sample_id", "age_bp", "cohort", "include")):
        sid = raw["sample_id"]
        if sid in table:
            raise ValueError(f"Duplicate metadata sample_id: {sid}")
        table[sid] = raw
    if any(sid not in table for sid in sample_ids):
        raise ValueError("Metadata must identify every original PSAM sample; no positional matching")
    output = []
    for sid in sample_ids:
        row = dict(table[sid])
        age = finite(row.get("age_bp"))
        cohort = row.get("cohort")
        if cohort not in {"ancient", "modern"}:
            raise ValueError(f"{sid}: explicit ancient/modern cohort is required")
        include = row.get("include", "1") == "1" and row.get("unrelated", "1") == "1"
        include = include and age is not None and age >= 0
        if cohort == "ancient" and row.get("region", row.get("stratum", ".")) in NA:
            include = False
        row.update(age_bp=age, include=bool(include), cohort=cohort,
                   region=row.get("region", row.get("stratum", ".")),
                   time_bin=(-1 if cohort == "modern" else int(age // bin_years)) if include else None)
        output.append(row)
    return output


def read_coordinates(path, max_pcs=10):
    from pca_contract import read_scores
    data, columns, audit = read_scores(path, max_pcs)
    return data, columns, audit['coordinate_kind']


def bin_groups(metadata, include_region=True):
    output = defaultdict(list)
    for index, sample in enumerate(metadata):
        if not sample["include"]:
            continue
        output[(sample["cohort"], sample["time_bin"], "all")].append(index)
        if include_region and sample["region"] not in NA:
            output[(sample["cohort"], sample["time_bin"], sample["region"])].append(index)
    return output


def fixed_panel_mask(matrix, metadata, min_samples=20, min_call_rate=.9):
    """Same loci in every eligible pooled ancient epoch and modern cohort."""
    import numpy as np
    groups = {key: values for key, values in bin_groups(metadata, False).items() if len(values) >= min_samples}
    ancient = [key for key in groups if key[0] == "ancient"]
    keep = np.ones(matrix.shape[0], dtype=bool)
    if len(ancient) < 2:
        return np.zeros(matrix.shape[0], dtype=bool), groups, "fewer_than_two_eligible_ancient_bins"
    for members in groups.values():
        keep &= np.isfinite(matrix[:, members]).mean(axis=1) >= min_call_rate
    return keep, groups, "ok"


def bootstrap_mean(values, rng, reps):
    import numpy as np
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return None, None, None
    point = float(values.mean())
    if len(values) < 2:
        return point, None, None
    # Resample people, never read fragments or locus-by-person rows.
    estimates = [float(values[rng.integers(0, len(values), len(values))].mean()) for _ in range(reps)]
    lo, hi = np.quantile(estimates, [.025, .975])
    return point, float(lo), float(hi)


def paired_sample_summaries(primary, targets, matrix, keep, metadata, min_call_rate, eligible_epochs=None):
    import numpy as np
    by_site = {t["site_id"]: t for t in targets}
    matches = defaultdict(lambda: {"case": [], "control": []})
    for index, entry in enumerate(primary):
        if not keep[index]:
            continue
        target = by_site[entry["site_id"]]
        group = target["group"].lower()
        if group not in {"case", "control"}:
            raise ValueError("Matched targets group must be case or control")
        matches[target["match_id"]][group].append(index)
    sums = {}
    diagnostics = []
    for match_id, members in matches.items():
        if len(members["case"]) != 1 or not members["control"]:
            diagnostics.append(dict(match_id=match_id, status="missing_case_or_controls_in_common_panel"))
            continue
        case_index = members["case"][0]
        family = by_site[primary[case_index]["site_id"]]["family"]
        case = np.asarray(matrix[case_index])
        control = np.asarray(matrix[members["control"]])
        valid = np.isfinite(case) & np.isfinite(control).all(axis=0)
        for measure in ("biallelic_hardcall_heterozygosity", "contrast_allele_fraction"):
            c = (case == 1).astype(float) if measure.startswith("biallelic") else case / 2
            t = (control == 1).mean(axis=0) if measure.startswith("biallelic") else control.mean(axis=0) / 2
            record = sums.setdefault((family, measure), [np.zeros(len(metadata)), np.zeros(len(metadata)), np.zeros(len(metadata), dtype=int), 0])
            record[0][valid] += c[valid]
            record[1][valid] += t[valid]
            record[2] += valid
            record[3] += 1
        diagnostics.append(dict(match_id=match_id, status="retained", family=family, n_controls=len(members["control"])))
    output = []
    for (family, measure), (case, control, called, total) in sorted(sums.items()):
        for index, sample in enumerate(metadata):
            n = int(called[index])
            usable = sample["include"] and n > 0 and n / total >= min_call_rate
            output.append(dict(sample_id=sample["sample_id"], family=family, measure=measure,
                cohort=sample["cohort"], age_bp=sample["age_bp"], region=sample["region"], time_bin=sample["time_bin"],
                n_matched_sets=total, n_complete_sets=n, call_rate=n / total,
                case_value=float(case[index] / n) if usable else None,
                control_value=float(control[index] / n) if usable else None,
                difference=float((case[index] - control[index]) / n) if usable else None,
                primary_trend_eligible=int(eligible_epochs is None or (sample["cohort"], sample["time_bin"], "all") in eligible_epochs),
                status="ok" if usable else "excluded_or_insufficient_matched_call_rate"))
    return output, diagnostics


def summarize_time(samples, bin_years, min_samples, reps, seed):
    import numpy as np
    groups = defaultdict(list)
    for row in samples:
        if row["status"] != "ok":
            continue
        for region in ("all", row["region"]):
            if region in NA:
                continue
            groups[row["family"], row["measure"], row["cohort"], row["time_bin"], region].append(row)
    rng = np.random.default_rng(seed)
    output = []
    for (family, measure, cohort, time_bin, region), records in sorted(groups.items()):
        n = len(records)
        point, lo, hi = bootstrap_mean([r["difference"] for r in records], rng, reps)
        output.append(dict(family=family, measure=measure, cohort=cohort, time_bin=time_bin, region=region,
            start_bp=0 if cohort == "modern" else time_bin * bin_years,
            end_bp=0 if cohort == "modern" else (time_bin + 1) * bin_years,
            n_samples=n, case_mean=float(np.mean([r["case_value"] for r in records])),
            control_mean=float(np.mean([r["control_value"] for r in records])),
            difference=point, ci_low=lo, ci_high=hi,
            status="ok" if n >= min_samples else "small_bin_descriptive_only"))
    return output


def adjusted_slope(records, metadata, coordinates, pc_columns):
    """Ancient people only; HC3 OLS on paired case-minus-control index."""
    import numpy as np
    from scipy.stats import t as student_t
    meta = {r["sample_id"]: r for r in metadata}
    eligible = [r for r in records if r["status"] == "ok" and r["cohort"] == "ancient" and r.get("primary_trend_eligible", 1) == 1]
    if pc_columns:
        eligible = [r for r in eligible if r["sample_id"] in coordinates]
    use_coverage = bool(eligible) and all(finite(meta[r["sample_id"]].get("target_coverage")) is not None for r in eligible)
    regions = sorted({r["region"] for r in eligible})
    sources = sorted({meta[r["sample_id"]].get("data_source", ".") for r in eligible})
    names = ["intercept", "toward_present_per_1000y"] + pc_columns
    names += ["region:" + x for x in regions[1:]] + ["source:" + x for x in sources[1:]]
    if use_coverage:
        names.append("log1p_target_coverage")
    values, y = [], []
    for row in eligible:
        sample = meta[row["sample_id"]]
        values.append([1., -float(row["age_bp"]) / 1000] + coordinates.get(row["sample_id"], [])[:len(pc_columns)] +
            [float(row["region"] == x) for x in regions[1:]] +
            [float(sample.get("data_source", ".") == x) for x in sources[1:]] +
            ([math.log1p(float(sample["target_coverage"]))] if use_coverage else []))
        y.append(row["difference"])
    base = dict(n=len(y), model="HC3_OLS_ancient_only", direction="positive = case-control difference increases toward present",
                covariates=";".join(names), pcs=len(pc_columns), coverage_adjusted=int(use_coverage))
    if len(y) < len(names) + 20 or len({r["time_bin"] for r in eligible}) < 3:
        return dict(base, status="insufficient_samples_or_time_bins")
    X = np.asarray(values, dtype=float)
    # Drop only constant nuisance columns; retain time and report collinearity.
    active = [0, 1] + [i for i in range(2, X.shape[1]) if np.std(X[:, i]) > 1e-12]
    X = X[:, active]
    base["covariates"] = ";".join(names[i] for i in active)
    for i in range(1, X.shape[1]):
        X[:, i] -= X[:, i].mean()
        if i > 1:
            X[:, i] /= X[:, i].std()
    condition = float(np.linalg.cond(X))
    base["condition_number"] = condition if math.isfinite(condition) else None
    if np.linalg.matrix_rank(X) < X.shape[1] or condition > 1e8:
        return dict(base, status="time_covariate_collinearity_not_estimable")
    y = np.asarray(y)
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ y
    residual = y - X @ beta
    leverage = np.einsum("ij,jk,ik->i", X, inv, X)
    if np.max(leverage) >= .999:
        return dict(base, status="high_leverage_not_estimable")
    weighted = X * (residual / (1 - leverage))[:, None]
    variance = inv @ weighted.T @ weighted @ inv
    se = math.sqrt(max(0., float(variance[1, 1])))
    if se <= 0 or not math.isfinite(se):
        return dict(base, status="zero_or_invalid_residual_variance", slope=float(beta[1]))
    df = len(y) - X.shape[1]
    p = float(2 * student_t.sf(abs(beta[1] / se), df))
    crit = float(student_t.ppf(.975, df))
    # Time VIF warns about weaker identifiability below hard singularity.
    other = np.delete(X, 1, axis=1)
    time_res = X[:, 1] - other @ np.linalg.lstsq(other, X[:, 1], rcond=None)[0]
    vif = float(np.sum(X[:, 1] ** 2) / np.sum(time_res ** 2))
    return dict(base, status="exploratory_adjusted_association", slope=float(beta[1]), se=se,
                ci_low=float(beta[1] - crit * se), ci_high=float(beta[1] + crit * se), p=p, time_vif=vif)


def bh(records, field="p"):
    ordered = sorted((float(row[field]), i) for i, row in enumerate(records) if row.get(field) is not None)
    running = 1.
    for rank in range(len(ordered), 0, -1):
        p, index = ordered[rank - 1]
        running = min(running, p * len(ordered) / rank)
        records[index]["q"] = running


def figures(time_rows, sample_rows, coordinates, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usable = [r for r in time_rows if r["region"] == "all" and r["status"] == "ok"]
    if usable:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for ax, measure in zip(axes, ("biallelic_hardcall_heterozygosity", "contrast_allele_fraction")):
            for family in sorted({r["family"] for r in usable}):
                records = sorted([r for r in usable if r["family"] == family and r["measure"] == measure and r["cohort"] == "ancient"], key=lambda r: r["start_bp"])
                if not records:
                    continue
                x = [(r["start_bp"] + r["end_bp"]) / 2000 for r in records]
                y = [r["difference"] for r in records]
                line = ax.plot(x, y, marker="o", label=family)[0]
                lo = [r["ci_low"] for r in records]
                hi = [r["ci_high"] for r in records]
                if all(v is not None for v in lo + hi):
                    ax.fill_between(x, lo, hi, color=line.get_color(), alpha=.17)
            ax.axhline(0, color="grey", lw=.8)
            ax.invert_xaxis()
            ax.set(xlabel="Age (thousand years BP)", ylabel="HERV minus matched control", title=measure.replace("_", " "))
            if ax.lines:
                ax.legend(fontsize=8)
        fig.suptitle("Fixed modern-locus imputed nucleotide panel · descriptive")
        fig.tight_layout()
        fig.savefig(out / "ancient_fixed_panel_time.png", dpi=180)
        fig.savefig(out / "ancient_fixed_panel_time.pdf")
        plt.close(fig)
    else:
        for ext in ("png", "pdf"):
            (out / ("ancient_fixed_panel_time." + ext)).unlink(missing_ok=True)
    projected = [r for r in sample_rows if r["status"] == "ok" and r["measure"] == "biallelic_hardcall_heterozygosity" and len(coordinates.get(r["sample_id"], [])) >= 2]
    if projected:
        families = sorted({r["family"] for r in projected})
        fig = plt.figure(figsize=(6 * len(families), 5))
        for index, family in enumerate(families, 1):
            rs = [r for r in projected if r["family"] == family and r["cohort"] == "ancient"]
            if not rs:
                continue
            ax = fig.add_subplot(1, len(families), index, projection="3d")
            point = ax.scatter([coordinates[r["sample_id"]][0] for r in rs], [coordinates[r["sample_id"]][1] for r in rs],
                [r["age_bp"] for r in rs], c=[r["difference"] for r in rs], cmap="coolwarm", s=3, alpha=.5)
            ax.set(xlabel="Shared PC1", ylabel="Shared PC2", zlabel="Age BP", title=family)
            fig.colorbar(point, ax=ax, shrink=.5, label="Paired target-versus-rest heterozygosity difference")
        fig.savefig(out / "ancient_shared_pca_time.png", dpi=160)
        fig.savefig(out / "ancient_shared_pca_time.pdf")
        plt.close(fig)
    else:
        # A rerun with --pcs 0/1 must not leave an old, apparently current plot.
        for ext in ("png", "pdf"):
            (out / ("ancient_shared_pca_time." + ext)).unlink(missing_ok=True)


def run(args):
    import numpy as np
    from evolve import read_psam_ids
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    reset_run_outputs(args, out)
    targets = load_targets(args.targets)
    prefix = str(args.pgen_prefix)
    pvar = panel_pvar(prefix)
    sample_ids = read_psam_ids(prefix + ".psam")
    metadata = metadata_in_order(args.metadata, sample_ids, args.bin_years)
    coordinates, pc_columns, coord_kind = read_coordinates(args.coordinates, args.pcs)
    with GenomeMapper(args.fasta38, args.fasta37, args.chain38to37, args.chain37to38) as mapper:
        mapped = map_targets(targets, mapper)
    map_fields = ["site_id", "chrom", "pos", "ref", "alts", "primary_alt", "family", "group", "match_id", "mapping_status",
                  "chrom37", "pos37", "ref37", "alts37", "strand37to38", "reference_changed"]
    write_rows(out / "mapped_targets.tsv", mapped, map_fields)
    indexed = index_panel(mapped, pvar, args.hq_variants)
    selected_alts, available = target_measurements(mapped, indexed)
    fields = ["site_id", "variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "base38", "status", "source_id",
              "allele_id", "allele_idx", "allele", "role", "panel_status", "panel_ref37", "panel_alleles37", "counted_base37",
              "measurement_mode", "counting_scope", "reference_count_scope", "comparison_kind"]
    write_rows(out / "allele_availability.tsv", available, fields)
    primary, primary_audit = select_primary_measurements(mapped, selected_alts, available)
    primary_fields = fields + ["primary_rule", "counted_base38", "invert_dosage"]
    write_rows(out / "selected_variants.tsv", primary, primary_fields)
    dump_json(out / "measurement_plans.json", dict(schema="nucleotide-count-1", selected_primary=primary,
        rule="Direct encoded ALT or native multiallelic allele index; unknown split REF is unavailable"))
    write_rows(out / "primary_contrast_audit.tsv", primary_audit,
               ["site_id", "predefined_primary_alt38", "primary_rule", "primary_status", "measurement_status",
                "counted_base38", "counted_base37", "pvar_index", "allele_idx", "measurement_mode", "counting_scope", "invert_dosage"])
    chosen = {r["site_id"] for r in primary}
    chosen_rows = {r["site_id"]: r for r in primary}
    primary_status = {r["site_id"]: r for r in primary_audit}
    for row in mapped:
        row["primary_contrast_available"] = int(row["site_id"] in chosen)
        row["primary_status"] = primary_status[row["site_id"]]["primary_status"]
        row["primary_rule"] = primary_status[row["site_id"]]["primary_rule"]
        if row["site_id"] in chosen_rows:
            selected = chosen_rows[row["site_id"]]
            # Read validation needs two distinct anchor labels. The counted primary
            # nucleotide may itself be the GRCh37 REF after a reference change.
            anchor_alt37 = selected["counted_base37"]
            if anchor_alt37 == row["ref37"]:
                anchor_alt37 = orient(row["ref"], row["strand37to38"])
            row.update(primary_ref37=row["ref37"], primary_alt37=anchor_alt37,
                       primary_alt38=selected["allele"], primary_counted_base37=selected["counted_base37"])
    map_fields += ["primary_contrast_available", "primary_status", "primary_rule", "primary_ref37", "primary_alt37", "primary_alt38", "primary_counted_base37"]
    write_rows(out / "mapped_targets.tsv", mapped, map_fields)
    write_rows(out / "site_availability.tsv", mapped, map_fields)
    base_manifest = dict(version=VERSION, input_targets=len(targets), mapping_counts=dict(Counter(r["mapping_status"] for r in mapped)),
        allele_availability=dict(Counter(r["panel_status"] for r in available)), primary_contrasts=len(primary),
        primary_contrast_status=dict(Counter(r["primary_status"] for r in primary_audit)),
        primary_contrast_rules=dict(Counter(r["primary_rule"] for r in primary)),
        count_kind="imputed_hardcall_nucleotide_count", heterozygosity_definition="counted_nucleotide_dosage_equals_one",
        measurement_schema="nucleotide-count-1", coordinate_kind=coord_kind, limitations=LIMITATIONS,
        measure_definitions={"biallelic_hardcall_heterozygosity":"counted nucleotide dosage equals one; target versus all remaining nucleotides",
                             "contrast_allele_fraction":"counted nucleotide dosage divided by two"},
        parameters={k: v for k, v in vars(args).items() if k != "handler"})
    if not primary:
        result = dict(base_manifest, status="no_usable_exact_high_quality_contrast")
        dump_json(out / "manifest.json", result)
        return result
    temporary = out / "contrast_dosage.partial.npy"
    matrix = np.lib.format.open_memmap(temporary, mode="w+", dtype="float32", shape=(len(primary), len(sample_ids)))
    matrix[:] = np.nan
    try:
        measurement_qc = Counter()
        groups = bin_groups(metadata)
        trajectory_fields = ["site_id", "variant_id", "counted_base38", "counted_base37", "allele_idx", "invert_dosage", "primary_rule",
                             "counting_scope", "cohort", "time_bin", "region", "start_bp", "end_bp", "n_eligible", "n_called", "an", "counted_ac", "counted_af", "hardcall_heterozygosity", "call_rate", "status"]
        def trajectory_records():
            for index, (entry, counted) in enumerate(iter_primary_counts(prefix, len(sample_ids), primary, audit=measurement_qc)):
                matrix[index] = np.where(counted >= 0, counted, np.nan)
                for (cohort, time_bin, region), members in sorted(groups.items()):
                    values = matrix[index, members]
                    values = values[np.isfinite(values)]
                    n = len(values)
                    yield dict(site_id=entry["site_id"], variant_id=entry["variant_id"], cohort=cohort,
                        counted_base38=entry["counted_base38"], counted_base37=entry["counted_base37"],
                        allele_idx=entry["allele_idx"], invert_dosage=0, primary_rule=entry["primary_rule"],
                        counting_scope=entry["counting_scope"],
                        time_bin=time_bin, region=region, start_bp=0 if cohort == "modern" else time_bin * args.bin_years,
                        end_bp=0 if cohort == "modern" else (time_bin + 1) * args.bin_years,
                        n_eligible=len(members), n_called=n, an=2*n, counted_ac=int(values.sum()) if n else None,
                        counted_af=float(values.mean()/2) if n else None, hardcall_heterozygosity=float(np.mean(values == 1)) if n else None,
                        call_rate=n/len(members), status="observed_imputed_nucleotide_count" if n else "no_called_genotypes")
        write_rows(out / "contrast_trajectories.tsv.gz", trajectory_records(), trajectory_fields)
        keep, epochs, panel_status = fixed_panel_mask(matrix, metadata, args.min_samples, args.min_call_rate)
        panel_rows = [dict(entry, in_fixed_panel=int(keep[index])) for index, entry in enumerate(primary)]
        write_rows(out / "fixed_common_panel.tsv", panel_rows, primary_fields + ["in_fixed_panel"])
        write_rows(out / "sample_bin_composition.tsv", [dict(cohort=c, time_bin=t, region=r, n=len(v)) for (c,t,r),v in sorted(groups.items())], ["cohort", "time_bin", "region", "n"])
        sample_rows, match_rows = paired_sample_summaries(primary, targets, matrix, keep, metadata, args.min_call_rate, epochs)
        sample_fields = ["sample_id", "family", "measure", "cohort", "age_bp", "region", "time_bin", "n_matched_sets", "n_complete_sets", "call_rate", "case_value", "control_value", "difference", "primary_trend_eligible", "status"]
        write_rows(out / "sample_contrast_burden.tsv", sample_rows, sample_fields)
        write_rows(out / "matched_panel_audit.tsv", match_rows, ["match_id", "status", "family", "n_controls"])
        time_rows = summarize_time(sample_rows, args.bin_years, args.min_samples, args.bootstrap, args.seed)
        time_fields = ["family", "measure", "cohort", "time_bin", "region", "start_bp", "end_bp", "n_samples", "case_mean", "control_mean", "difference", "ci_low", "ci_high", "status"]
        write_rows(out / "fixed_panel_time_summary.tsv", time_rows, time_fields)
        trend_groups = defaultdict(list)
        for record in sample_rows:
            trend_groups[record["family"], record["measure"]].append(record)
        trend = []
        for (family, measure), records in sorted(trend_groups.items()):
            trend.append(dict(family=family, measure=measure, **adjusted_slope(records, metadata, coordinates, pc_columns)))
        bh(trend)
        trend_fields = ["family", "measure", "status", "n", "model", "direction", "covariates", "pcs", "coverage_adjusted", "condition_number", "time_vif", "slope", "se", "ci_low", "ci_high", "p", "q"]
        write_rows(out / "time_trend.tsv", trend, trend_fields)
        if not args.no_plots:
            figures(time_rows, sample_rows, coordinates, out)
        result = dict(base_manifest, status="complete_descriptive_imputed_panel", fixed_panel_status=panel_status,
            measurement_qc=dict(measurement_qc), fixed_panel_contrasts=int(keep.sum()), fixed_panel_epochs=len(epochs), matched_sets_retained=sum(r["status"] == "retained" for r in match_rows),
            primary_trend_epoch_rule="Only pooled epochs meeting min_samples and used to define the fixed common panel; sparse epochs remain descriptive with primary_trend_eligible=0",
            qualified_epochs=[dict(cohort=c, time_bin=t, n_metadata=len(members)) for (c,t,_),members in sorted(epochs.items())],
            pca_overlay_status="available_at_least_two_pcs" if coordinates and len(pc_columns) >= 2 else "not_rendered_fewer_than_two_pcs_or_missing_coordinates",
            trend_status=dict(Counter(r["status"] for r in trend)))
        dump_json(out / "manifest.json", result)
        return result
    finally:
        del matrix
        temporary.unlink(missing_ok=True)


def read_evidence(target, counts, min_depth=8, min_support=3, min_fraction=.1):
    """A read count is not an allele count; flag >2 bases within one person."""
    if any(not isinstance(v, int) or v < 0 for v in counts.values()):
        raise ValueError("Read counts must be nonnegative integers")
    depth = sum(counts.values())
    anchor = (target["ref"], target.get("anchor_alt", sorted(target["alts"].split(","))[0]))
    if anchor[1] not in target["alts"].split(","):
        raise ValueError("Read-audit anchor ALT must belong to target ALTs")
    extra = sorted(set(target["alts"].split(",")) - set(anchor))
    supported = {base for base in "ACGT" if counts[base] >= min_support and counts[base] / max(1, depth) >= min_fraction}
    damage = {base for base in extra if (base == "T" and "C" in anchor) or (base == "A" and "G" in anchor)}
    transitions = {frozenset("AG"), frozenset("CT")}
    transversions = {base for base in extra if all(frozenset((base, a)) not in transitions for a in anchor)}
    callable_ = depth >= min_depth
    status = "low_depth_uninformative" if not callable_ else "excess_bases_in_one_diploid_sample" if len(supported) > 2 else "read_support_only"
    return dict(depth=depth, anchor_ref=anchor[0], anchor_alt=anchor[1], extra_bases=",".join(extra),
        supported_bases=",".join(sorted(supported)), n_supported_bases=len(supported),
        extra_supported=",".join(sorted(set(extra) & supported)) if callable_ else None,
        extra_damage_susceptible=",".join(sorted(damage)),
        transversion_extra_supported=",".join(sorted(transversions & supported)) if callable_ else None,
        callable=int(callable_), status=status)


def read_audit(args):
    targets = load_targets(args.targets, args.build)
    lookup = {(r["chrom"], r["pos"]): r for r in targets}
    metadata = {}
    if args.metadata:
        for record in rows(args.metadata, ("sample_id", "age_bp", "cohort")):
            if record["sample_id"] in metadata:
                raise ValueError("Duplicate metadata sample_id in read audit")
            metadata[record["sample_id"]] = record
    seen = set()
    summary = Counter()
    output = []
    for record in rows(args.counts, ("sample_id", "chrom", "pos", "A", "C", "G", "T")):
        sid, chrom, pos = record["sample_id"], canonical_chrom(record["chrom"]), int(record["pos"])
        if args.metadata and sid not in metadata:
            raise ValueError(f"Read-count sample absent from supplied metadata: {sid}")
        key = sid, chrom, pos
        if key in seen:
            raise ValueError(f"Duplicate sample/locus read-count row: {key}")
        seen.add(key)
        target = lookup.get((chrom, pos))
        if target is None:
            raise ValueError(f"Read count outside declared target panel/build: {key}")
        counts = {base: int(record[base]) for base in "ACGT"}
        evidence = read_evidence(target, counts, args.min_depth, args.min_support, args.min_fraction)
        if record.get("depth_capped", "0") == "1":
            evidence.update(status="depth_capped_uninformative", callable=0, extra_supported=None, transversion_extra_supported=None)
        age = finite(metadata.get(sid, {}).get("age_bp"))
        cohort = metadata.get(sid, {}).get("cohort", ".")
        summary[evidence["status"]] += 1
        output.append(dict(sample_id=sid, site_id=target["site_id"], chrom=chrom, pos=pos, build=args.build,
            family=target["family"], group=target["group"], age_bp=age, cohort=cohort, **counts, **evidence))
    out = Path(args.out)
    fields = ["sample_id", "site_id", "chrom", "pos", "build", "family", "group", "cohort", "age_bp", "A", "C", "G", "T", "depth", "anchor_ref", "anchor_alt", "extra_bases", "supported_bases", "n_supported_bases", "extra_supported", "extra_damage_susceptible", "transversion_extra_supported", "callable", "status"]
    write_rows(out / "direct_read_evidence.tsv", output, fields)
    result = dict(version=VERSION, status="read_evidence_audit_only", n_sample_loci=len(output), counts=dict(summary),
        rules=dict(min_depth=args.min_depth, min_support=args.min_support, min_fraction=args.min_fraction),
        limitations=["Fragment counts are not chromosomes or genotype likelihoods and are never added to imputed allele counts.",
                    "Zero support is not proof of absence. More than two supported bases within a diploid person flags mapping/contamination/error.",
                    "Anchors use the primary imputed contrast when provided in mapped targets; otherwise reference plus alphabetical first ALT. They are not ancestral/derived or time-ordered alleles.",
                    "External count tables need independently documented mapping, base, fragment and damage filtering.",
                    "These small-panel support records do not estimate population tri-allelic prevalence or a temporal enrichment trend."])
    dump_json(out / "read_audit_manifest.json", result)
    return result


def read_counts(args):
    import pysam
    targets = load_targets(args.targets, args.build)
    if len(targets) > args.max_targets:
        raise ValueError("Read counting is a bounded validation step; reduce targets or explicitly raise --max-targets")
    samples = list(rows(args.bam_manifest, ("sample_id", "bam")))
    if len({r["sample_id"] for r in samples}) != len(samples):
        raise ValueError("Duplicate BAM sample_id")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    def records():
        with pysam.FastaFile(args.fasta) as reference:
            for sample in samples:
                if sample.get("build", args.build) != args.build:
                    raise ValueError("BAM manifest build disagrees with --build")
                mode = "rc" if sample["bam"].endswith(".cram") else "rb"
                with pysam.AlignmentFile(sample["bam"], mode, reference_filename=args.fasta) as bam:
                    if not bam.has_index():
                        raise ValueError(f"Indexed BAM/CRAM required: {sample['bam']}")
                    for target in targets:
                        chrom = target["chrom"] if target["chrom"] in bam.references else "chr" + target["chrom"]
                        fasta_chrom = target["chrom"] if target["chrom"] in reference.references else "chr" + target["chrom"]
                        if chrom not in bam.references or fasta_chrom not in reference.references:
                            raise ValueError("BAM/reference chromosome absent; check build")
                        pos = target["pos"] - 1
                        if reference.fetch(fasta_chrom, pos, pos + 1).upper() != target["ref"]:
                            raise ValueError("Target REF disagrees with declared read-count FASTA")
                        fragments, conflicts, capped = {}, set(), False
                        for column in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="nofilter", min_base_quality=0,
                                                ignore_overlaps=False, ignore_orphans=False, max_depth=args.max_depth):
                            if column.reference_pos != pos:
                                continue
                            capped = column.nsegments >= args.max_depth
                            for pileup in column.pileups:
                                read = pileup.alignment
                                q = pileup.query_position
                                if q is None or pileup.is_del or pileup.is_refskip or read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
                                    continue
                                if read.mapping_quality < args.mapq or not read.query_qualities or read.query_qualities[q] < args.baseq:
                                    continue
                                if q - read.query_alignment_start < args.trim_ends or read.query_alignment_end - q - 1 < args.trim_ends:
                                    continue
                                base = read.query_sequence[q].upper()
                                if base not in BASES:
                                    continue
                                name = (read.get_tag("RG") if read.has_tag("RG") else ".", read.query_name)
                                if name in conflicts:
                                    continue
                                if name in fragments and fragments[name] != base:
                                    fragments.pop(name)
                                    conflicts.add(name)
                                else:
                                    fragments[name] = base
                        count = Counter(fragments.values())
                        yield dict(sample_id=sample["sample_id"], chrom=target["chrom"], pos=target["pos"], **{x:count[x] for x in "ACGT"},
                                   depth_capped=int(capped), discordant_fragments=len(conflicts))
    write_rows(out / "read_counts.tsv.gz", records(), ["sample_id", "chrom", "pos", "A", "C", "G", "T", "depth_capped", "discordant_fragments"])
    result = dict(status="complete_fragment_counts", build=args.build, targets=len(targets), samples=len(samples),
                  mapq=args.mapq, baseq=args.baseq, trim_ends=args.trim_ends,
                  filtering="Primary nonduplicate aligned reads; both aligned ends trimmed; overlapping mates counted once; discordant mates discarded",
                  limitations=["Short-read repeat mapping and contamination may remain despite quality filters.", "Raw evidence, not diploid genotypes or insertion calls."])
    dump_json(out / "read_counts_manifest.json", result)
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sp = p.add_subparsers(dest="command", required=True)
    run_p = sp.add_parser("run", help="Exact ancient-panel audit and fixed-panel matched descriptive trajectories")
    for name in ("targets", "pgen-prefix", "metadata", "hq-variants", "fasta38", "fasta37", "chain38to37", "chain37to38", "out"):
        run_p.add_argument("--" + name, required=True)
    run_p.add_argument("--coordinates")
    run_p.add_argument("--pcs", type=int, default=10)
    run_p.add_argument("--bin-years", type=int, default=2000)
    run_p.add_argument("--min-samples", type=int, default=20)
    run_p.add_argument("--min-call-rate", type=float, default=.9)
    run_p.add_argument("--bootstrap", type=int, default=500)
    run_p.add_argument("--seed", type=int, default=20261027)
    run_p.add_argument("--no-plots", action="store_true")
    run_p.set_defaults(handler=run)
    audit = sp.add_parser("read-audit", help="Separate direct-read support audit; never treat reads as allele counts")
    for name in ("targets", "counts", "out"):
        audit.add_argument("--" + name, required=True)
    audit.add_argument("--metadata")
    audit.add_argument("--build", choices=["37", "38"], required=True)
    audit.add_argument("--min-depth", type=int, default=8)
    audit.add_argument("--min-support", type=int, default=3)
    audit.add_argument("--min-fraction", type=float, default=.1)
    audit.set_defaults(handler=read_audit)
    count = sp.add_parser("read-counts", help="Optional small-panel fragment counting from indexed local BAM/CRAM")
    for name in ("targets", "bam-manifest", "fasta", "out"):
        count.add_argument("--" + name, required=True)
    count.add_argument("--build", choices=["37", "38"], required=True)
    count.add_argument("--mapq", type=int, default=30)
    count.add_argument("--baseq", type=int, default=30)
    count.add_argument("--trim-ends", type=int, default=5)
    count.add_argument("--max-targets", type=int, default=2000)
    count.add_argument("--max-depth", type=int, default=100000)
    count.set_defaults(handler=read_counts)
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.command == "run" and (args.bin_years < 1 or args.min_samples < 2 or args.bootstrap < 100 or not 0 < args.min_call_rate <= 1 or args.pcs < 0):
        p.error("Positive time bins, min-samples>=2, bootstrap>=100, 0<call-rate<=1 and nonnegative PCs required")
    if args.command == "read-audit" and (args.min_depth < 1 or args.min_support < 1 or not 0 < args.min_fraction <= 1):
        p.error("Positive depth/support and 0<min-fraction<=1 required")
    if args.command == "read-counts" and (min(args.mapq, args.baseq, args.trim_ends) < 0 or min(args.max_targets, args.max_depth) < 1):
        p.error("Read thresholds must be nonnegative and bounds positive")
    result = args.handler(args)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"herv-time: {exc}", file=sys.stderr)
        raise SystemExit(2)

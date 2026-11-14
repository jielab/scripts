#!/usr/bin/env python3
"""Consolidated chronology workflow; internal stages share one explicit module."""
from __future__ import annotations

# Chronology
# ======================================================================
# ========================================================================
# herv time
# ========================================================================
_description_herv_time = 'HERV candidate trajectories in a fixed ancient panel; direct-read audit is separate.\n\nThe default input is the Reich imputed, biallelic GRCh37 PGEN. Its absence of an\nextra base is NOT an observed zero. This module never infers HERV insertions,\nfirst appearances, or selection coefficients from this panel.'
#!/usr/bin/env python3
"""HERV candidate trajectories in a fixed ancient panel; direct-read audit is separate.

The default input is the Reich imputed, biallelic GRCh37 PGEN. Its absence of an
extra base is NOT an observed zero. This module never infers HERV insertions,
first appearances, or selection coefficients from this panel.
"""

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
    from pca import frame_read_scores as read_scores
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
        def trajectory_records(matrix=matrix):
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
    p = argparse.ArgumentParser(description=_description_herv_time)
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


def _entry_herv_time():
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"herv-time: {exc}", file=sys.stderr)
        raise SystemExit(2)


# ========================================================================
# herv allele time
# ========================================================================
_description_herv_allele_time = 'Full-catalog, full-PSAM ancient allele trajectories with bounded genotype RAM.\n\nEvery modern nucleotide is audited independently against its actual panel allele\nindex. An unrepresented base is unavailable, never AC=0. Native multiallelic\nPGEN preserves direct REF and ALT counts; biallelic/split REF is unavailable\nwithout a joint-genotype provenance contract. A single ALT complement is never\nreported as REF. The existing analysis remains an ALT-trajectory analysis.'
#!/usr/bin/env python3
"""Full-catalog, full-PSAM ancient allele trajectories with bounded genotype RAM.

Every modern nucleotide is audited independently against its actual panel allele
index. An unrepresented base is unavailable, never AC=0. Native multiallelic
PGEN preserves direct REF and ALT counts; biallelic/split REF is unavailable
without a joint-genotype provenance contract. A single ALT complement is never
reported as REF. The existing analysis remains an ALT-trajectory analysis.
"""
from analysis_integrity import transactional
from collections import OrderedDict
import numpy as np
from scipy.stats import t as student_t
from allele_measurement import iter_pgen_chunks

ALLELE_VERSION = "herv-all-allele-time-2.0-nucleotide-counts"
ALLELE_LIMITATIONS = [
    "The supplied catalog defines the inference scope: HERV cases/controls or sampled targets in frozen genome windows. Neither is all human genomic SNVs or all ancient variation; a window match is not a matched ancient site pair.",
    "Modern survivor ascertainment and the imputation reference restrict interpretation; unrepresented extra bases are unavailable, not zero, and missing hardcalls are not absence.",
    "Every represented ALT is counted through its actual PGEN allele index, including native multiallelic indices 0..3 after build-reference changes. Unknown split REF is unavailable; no single-ALT complement is called REF. Counts describe encoded panel genotypes, not unconditional population frequencies.",
    "Hardcall heterozygosity is the fraction carrying exactly one copy of the counted nucleotide (target versus all other nucleotides), not full-locus multiallelic heterozygosity.",
    "The primary chronology is explicitly unrelated ancient people aged 0 to max_age_bp inclusive. Older, related, relationship-unknown and modern samples remain separately labelled descriptive groups.",
    "Modern people are a separate endpoint, not ancient observations assigned an artificial continuous date. The available PSAM is not necessarily the full publication cohort.",
    "HC3 OLS slopes of imputed allele fraction are two-sided descriptive associations, not selection coefficients, allele ages or evidence of increasing genomic diversity.",
    "Only shared fixed PCs are adjusted; missing requested PCs/region/source/coverage produce complete-case exclusions. Confounding, imputation bias and date uncertainty remain.",
    "Global BH and conservative Benjamini-Yekutieli FDR include all tested ALT contrasts; REF availability is audited without additional trajectory tests. AC5 primary FDR excludes low-count ALT and incomplete strict filters.",
    "modern_ac10 is an individual ALT count threshold. The original AC10 locus subset instead requires locus_ac10_supported=1: at least two modern ALT bases each have AC>=10.",
    "Family summaries resample physical blocks; matched summaries resample entire chromosomes containing both case and controls. Intervals condition on the current cohort, not a new random sample of people.",
    "Family means are unadjusted descriptions of a fixed represented locus panel. The separate per-allele model supplies covariate adjustment; ancestry composition can still change descriptive curves.",
    "Matching is at the locus/role level; group AF differences are not same-person complete-pair genotype estimates. The modern AF-matched candidate design cannot independently test enrichment or modern AF differences.",
]
ALLELE_OWNED = ("mapped_catalog.tsv.gz", "all_allele_audit.tsv.gz", "allele_trajectories.tsv.gz", "measurement_plans.json", "allele_callability.tsv.gz",
         "adjusted_allele_trends.tsv.gz", "sample_coverage.tsv.gz", "sample_time_region_coverage.tsv",
         "panel_availability_counts.tsv", "catalog_coverage.tsv", "fixed_panel_family_trajectories.tsv",
         "matched_locus_role_trajectories.tsv", "matched_locus_audit.tsv", "model_covariates.json",
         "family_allele_time.png", "family_allele_time.pdf", "allele_time_heatmap.png", "allele_time_heatmap.pdf",
         "allele_panel_availability.png", "allele_panel_availability.pdf", "trends.raw.tsv.gz")
ALLELE_LOCUS_FIELDS = ["locus_n_alt_ac5", "locus_n_alt_ac10", "locus_ac10_supported"]
ALLELE_AUDIT_FIELDS = ["site_id", "chrom", "pos", "ref", "allele", "role", "modern_ac", "modern_an", "modern_af",
                "modern_ac5", "modern_ac10", *ALLELE_LOCUS_FIELDS, "family", "group", "partition", "comparison_partition", "match_id",
                "strict_filter_complete", "mapping_status", "panel_status", "contrast_id", "pvar_index",
                "panel_ref37", "panel_alt37", "counted_base37", "invert_dosage", "reference_count_scope",
                "panel_alleles37", "allele_id", "allele_idx", "measurement_mode", "counting_scope", "comparison_kind", "source_id"]
ALLELE_TRAJECTORY_FIELDS = ["site_id", "contrast_id", "family", "group", "partition", "match_id", "role", "counted_base38",
                     "paired_modern_alt", "counting_scope", "modern_ac", "modern_ac5", "modern_ac10", *ALLELE_LOCUS_FIELDS, "strict_filter_complete",
                     "cohort", "age_scope", "time_bin", "region", "measurement_source", "start_bp", "end_bp", "n_eligible", "n_called",
                     "ac", "an", "af", "hardcall_heterozygosity", "call_rate", "fixed_panel", "status",
                     "counted_base37", "allele_id", "pvar_index", "allele_idx", "measurement_mode", "reference_count_scope"]
ALLELE_TREND_FIELDS = ["allele_index", "site_id", "contrast_id", "family", "group", "partition", "role", "counted_base38",
                "modern_ac", "modern_ac5", "modern_ac10", *ALLELE_LOCUS_FIELDS, "strict_filter_complete", "fixed_panel", "status", "n", "n_bins",
                "model_call_rate", "reason", "expected_pc_count", "covariates", "condition_number", "time_vif", "slope_per_1000yr_toward_present", "se",
                "ci_low", "ci_high", "p_two_sided", "q_bh_all_alt", "q_by_all_alt", "q_bh_primary_ac5", "q_by_primary_ac5",
                "allele_id", "allele_idx", "measurement_mode", "counting_scope"]


@contextlib.contextmanager
def allele_table(path, fields):
    path = Path(path)
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(partial, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        def emit(row):
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, ".") for k in fields})
        yield emit
    os.replace(partial, path)


def allele_load_catalog(path):
    result, positions, ids = [], set(), set()
    for raw in rows(path, ("chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id")):
        r = dict(raw)
        r["chrom"], r["pos"], r["ref"] = canonical_chrom(r["chrom"]), int(r["pos"]), r["ref"].upper()
        all_alts = r.get("observed_alts_all", r["alts"])
        if all_alts in NA: all_alts = r["alts"]
        # New full catalogs explicitly align all_alt_ac/an to their alts field.
        if "all_alt_ac" in r: all_alts = r["alts"]
        bases = all_alts.upper().split(",")
        if r["chrom"] not in {str(x) for x in range(1, 23)} or r["pos"] < 1 or r["ref"] not in BASES:
            raise ValueError("Catalog must contain canonical autosomal GRCh38 SNVs")
        if not bases or len(set(bases)) != len(bases) or any(b not in BASES or b == r["ref"] for b in bases):
            raise ValueError(f"Invalid catalog ALTs: {r['site_id']}")
        key = r["chrom"], r["pos"]
        if key in positions or r["site_id"] in ids: raise ValueError("Catalog requires unique positions and site_id values")
        positions.add(key); ids.add(r["site_id"])
        if r["group"] not in {"case", "control"}: raise ValueError("Catalog group must be case or control")
        if any(r[name] in NA for name in ("site_id","family","match_id")):
            raise ValueError("Catalog site_id, family and match_id must be nonmissing")
        counts = {}
        ac_text, an_text = r.get("all_alt_ac", r.get("alt_ac_all", "")), r.get("all_alt_an", r.get("alt_an_all", ""))
        if ac_text not in NA and an_text not in NA:
            ac, an = [int(x) for x in ac_text.split(",")], [int(x) for x in an_text.split(",")]
            if len(ac) != len(bases) or len(an) != len(bases) or any(a < 0 or n < a or n < 1 for a,n in zip(ac, an)):
                raise ValueError(f"Catalog all-ALT AC/AN mismatch: {r['site_id']}")
            counts = {b: (a,n) for b,a,n in zip(bases, ac, an)}
        locus_counts = {}
        for threshold in (5, 10):
            name = f"n_alt_ac{threshold}"
            recorded = int(r[name]) if r.get(name, ".") not in NA else None
            derived = sum(a >= threshold for a,_ in counts.values()) if counts else recorded
            if recorded is not None and (recorded < 0 or recorded > len(bases) or recorded != derived):
                raise ValueError(f"Catalog {name} disagrees with all-ALT counts: {r['site_id']}")
            locus_counts[f"locus_{name}"] = derived
        n10 = locus_counts["locus_n_alt_ac10"]
        supported = int(n10 >= 2) if n10 is not None else None
        if r.get("ac10_supported", ".") not in NA:
            recorded = int(r["ac10_supported"])
            if recorded not in {0, 1} or (supported is not None and recorded != supported):
                raise ValueError(f"Catalog ac10_supported disagrees with locus ALT counts: {r['site_id']}")
            supported = recorded
        locus_counts["locus_ac10_supported"] = supported
        for role in ("primary_alt", "secondary_alt"):
            if r.get(role, ".") not in NA and r[role] not in bases:
                raise ValueError(f"Catalog {role} is not an observed ALT: {r['site_id']}")
        if r.get("primary_alt", ".") not in NA and r.get("primary_alt") == r.get("secondary_alt"):
            raise ValueError("Primary and secondary ALT roles must be distinct")
        r.update(alts=",".join(bases), _counts=counts, **locus_counts,
                 partition=r.get("partition", r.get("case_region", "unclassified")),
                 strict_filter_complete=int(r.get("strict_filter_complete", "0") == "1"))
        result.append(r)
    return result


def allele_alt_role(target, base):
    if base == target.get("primary_alt"): return "primary"
    if base == target.get("secondary_alt"): return "secondary"
    ac = target["_counts"].get(base, (None,None))[0]
    return "additional_lowcount" if ac is not None and ac < 5 else "additional"


def allele_window_matched_catalog(targets):
    return any(t.get("window_id", "") not in NA and t.get("window_id")
               or "genome_window" in t.get("catalog_scope", "") for t in targets)


def allele_catalog_contrasts(mapped, availability):
    """Return ALT measurements and their audit, preserving the two-result API."""
    selected, audit, _ = catalog_measurements(mapped, availability, orient=orient, alt_role=allele_alt_role)
    # Keep the existing modern REF/ALT catalog identifier. This label identifies
    # the catalog nucleotide; counting_scope describes the actual measurement.
    for entry in audit:
        if entry["role"] != "REF":
            entry["contrast_id"] = f"{entry['site_id']}|{entry['ref']}>{entry['allele']}"
    return selected, audit


def allele_load_samples(path, psam, bin_years, max_age):
    ids = [r["IID"] for r in rows(psam, ("IID",), whitespace=True)]
    if len(set(ids)) != len(ids): raise ValueError("PSAM IID values must be unique")
    meta = {}
    for r in rows(path, ("sample_id", "age_bp", "cohort", "include")):
        if r["sample_id"] in meta: raise ValueError("Duplicate metadata sample_id")
        meta[r["sample_id"]] = r
    individual_counts = Counter(r.get('individual_id') for r in meta.values() if r.get('individual_id') not in NA | {None})
    result = []
    for sid in ids:
        raw = dict(meta.get(sid, {}))
        age = finite(raw.get("age_bp")); cohort = raw.get("cohort", "unknown")
        if cohort not in {"ancient", "modern", "unknown"}: raise ValueError("Metadata cohort must be ancient or modern")
        relation = raw.get("unrelated", "unknown")
        allowed = raw.get("include") == "1"
        if not raw: scope = "metadata_missing"
        elif not allowed: scope = "metadata_excluded"
        elif cohort == "unknown": scope = "invalid_or_unknown_cohort"
        elif cohort == "modern": scope = "modern"
        elif age is None or age < 0: scope = "invalid_ancient_date"
        elif relation == "0": scope = "related_ancient_descriptive"
        elif relation != "1": scope = "relationship_unknown_descriptive"
        elif individual_counts.get(raw.get("individual_id"), 0) > 1: scope = "duplicate_individual_descriptive"
        elif age > max_age: scope = "older_ancient_descriptive"
        else: scope = "primary_ancient"
        eligible = scope in {"modern", "primary_ancient", "older_ancient_descriptive", "related_ancient_descriptive", "relationship_unknown_descriptive"}
        time_bin = "modern" if scope == "modern" else (min(int(age // bin_years), math.ceil(max_age/bin_years)-1) if scope == "primary_ancient" else int(age // bin_years)) if eligible else None
        raw.update(sample_id=sid, age_bp=age, cohort=cohort, unrelated=relation, age_scope=scope, included=eligible,
                   time_bin=time_bin, region=raw.get("region", raw.get("stratum", "unknown")))
        if raw["region"] in NA: raw["region"] = "unknown"
        result.append(raw)
    return result, len(set(meta) - set(ids))


def allele_sample_groups(samples, bin_years, max_age, minimum):
    grouped = defaultdict(list)
    for index, s in enumerate(samples):
        if not s["included"]: continue
        key = s["cohort"], s["age_scope"], s["time_bin"]
        grouped[(*key, "all")].append(index)
        if s["region"] != "all": grouped[(*key, s["region"])].append(index)
    groups = []
    for (cohort, scope, epoch, region), members in sorted(grouped.items(), key=lambda x: str(x[0])):
        groups.append(dict(cohort=cohort, age_scope=scope, time_bin=epoch, region=region,
            start_bp=None if cohort == "modern" else epoch*bin_years,
            end_bp=None if cohort == "modern" else min((epoch+1)*bin_years,max_age) if scope == "primary_ancient" else (epoch+1)*bin_years,
            n_eligible=len(members), members=np.asarray(members, dtype=np.int32),
            qualified_epoch=int(region == "all" and scope in {"primary_ancient", "modern"} and len(members) >= minimum)))
    from allele_measurement import stats_technology as technology
    technical = defaultdict(list)
    for i,s in enumerate(samples):
        if s['age_scope'] not in {'primary_ancient','modern'}: continue
        source = technology(s) + '|' + s.get('data_source','unknown')
        technical[s['cohort'],s['age_scope'],s['time_bin'],s['region'],source].append(i)
    for (cohort,scope,epoch,region,source),members in sorted(technical.items(),key=lambda x:str(x[0])):
        groups.append(dict(cohort=cohort,age_scope=scope,time_bin=epoch,region=region,measurement_source=source,
            start_bp=None if cohort=='modern' else epoch*bin_years,
            end_bp=None if cohort=='modern' else min((epoch+1)*bin_years,max_age),
            n_eligible=len(members),members=np.asarray(members,dtype=np.int32),qualified_epoch=0))
    return groups


class allele_TrendModel:
    """Reuse design projections for equal call masks; never retain locus matrices."""
    def __init__(self, samples, groups, coordinates, pc_names, min_call_rate, expected_pc_count=10, unavailable_reason=None):
        self.unavailable_reason = unavailable_reason
        self.expected_pc_count = len(pc_names) if expected_pc_count is None else expected_pc_count
        if self.expected_pc_count and pc_names != [f'PC{i}_AVG' for i in range(1, self.expected_pc_count + 1)] and not unavailable_reason:
            raise ValueError('TrendModel requires all requested consecutive projected AVG PCs')
        self.samples, self.minimum = samples, min_call_rate
        epochs = {g["time_bin"] for g in groups if g["age_scope"] == "primary_ancient" and g["qualified_epoch"]}
        eligible = np.array([s["age_scope"] == "primary_ancient" and s["time_bin"] in epochs for s in samples])
        candidates = [s for s,k in zip(samples,eligible) if k]
        self.use_region = any(s["region"] != "unknown" for s in candidates)
        self.use_source = any(s.get("data_source", ".") not in NA for s in candidates)
        self.use_coverage = any(finite(s.get("target_coverage")) is not None for s in candidates)
        for i,s in enumerate(samples):
            eligible[i] &= (not pc_names or s["sample_id"] in coordinates)
            eligible[i] &= (not self.use_region or s["region"] != "unknown")
            eligible[i] &= (not self.use_source or s.get("data_source", ".") not in NA)
            coverage = finite(s.get("target_coverage"))
            eligible[i] &= (not self.use_coverage or coverage is not None and coverage >= 0)
        regions = sorted({s["region"] for s,k in zip(samples,eligible) if k}) if self.use_region else []
        sources = sorted({s.get("data_source", ".") for s,k in zip(samples,eligible) if k}) if self.use_source else []
        self.names = ["intercept", "toward_present_per_1000yr", *pc_names, *["region:"+v for v in regions[1:]], *["source:"+v for v in sources[1:]]]
        if self.use_coverage: self.names.append("log1p_target_coverage")
        self.X = np.zeros((len(samples), len(self.names)))
        for i,s in enumerate(samples):
            if not eligible[i]: continue
            self.X[i] = [1., -s["age_bp"]/1000, *coordinates.get(s["sample_id"], [])[:len(pc_names)],
                         *[float(s["region"] == v) for v in regions[1:]], *[float(s.get("data_source") == v) for v in sources[1:]],
                         *([math.log1p(float(s["target_coverage"]))] if self.use_coverage else [])]
        self.eligible, self.cache = (np.zeros_like(eligible) if unavailable_reason else eligible), OrderedDict()

    def prepare(self, mask):
        if self.unavailable_reason:
            return dict(status='adjustment_unavailable', reason=self.unavailable_reason, expected_pc_count=self.expected_pc_count, n=0)
        key = np.packbits(mask).tobytes()
        if key in self.cache:
            self.cache.move_to_end(key); return self.cache[key]
        indices = np.flatnonzero(mask)
        bins = {self.samples[i]["time_bin"] for i in indices}
        base = dict(n=len(indices), n_bins=len(bins), model_call_rate=len(indices)/max(1,int(self.eligible.sum())))
        if base["model_call_rate"] < self.minimum:
            return dict(base, status="low_model_call_rate")
        if len(bins) < 3 or len(indices) < len(self.names)+20:
            return dict(base, status="insufficient_samples_or_time_bins")
        X = self.X[indices].copy()
        active = [0,1] + [i for i in range(2,X.shape[1]) if np.std(X[:,i]) > 1e-12]
        X = X[:,active]
        for i in range(1,X.shape[1]):
            X[:,i] -= X[:,i].mean()
            if i > 1: X[:,i] /= X[:,i].std()
        base["covariates"] = ";".join(self.names[i] for i in active)
        condition = float(np.linalg.cond(X)); base["condition_number"] = condition if math.isfinite(condition) else None
        if np.linalg.matrix_rank(X) < X.shape[1] or condition > 1e8:
            return dict(base, status="time_covariate_collinearity_not_estimable")
        # Work on X directly: forming X'X would square the condition number and
        # destabilize the HC3 time estimate when age and ancestry nearly align.
        projection = np.linalg.pinv(X)
        leverage = np.einsum("ij,ji->i", X, projection)
        if max(leverage) >= .999: return dict(base, status="high_leverage_not_estimable")
        other = np.delete(X,1,axis=1)
        residual_time = X[:,1] - other@np.linalg.lstsq(other,X[:,1],rcond=None)[0]
        base["time_vif"] = float(np.sum(X[:,1]**2)/np.sum(residual_time**2))
        value = dict(base, status="ready", indices=indices, X=X, projection=projection,
                     influence=projection[1]/(1-leverage), df=len(indices)-X.shape[1])
        self.cache[key] = value
        while len(self.cache) > 4: self.cache.popitem(last=False)
        return value

    def fit_chunk(self, values, allowed):
        output, batches = [{} for _ in values], defaultdict(list)
        for i,row in enumerate(values):
            if not allowed[i]:
                output[i] = dict(status="outside_fixed_panel_or_strict_filters")
                continue
            mask = self.eligible & (row >= 0)
            batches[np.packbits(mask).tobytes()].append(i)
        for key, members in batches.items():
            mask = np.unpackbits(np.frombuffer(key,dtype=np.uint8),count=len(self.samples)).astype(bool)
            model = self.prepare(mask)
            public = {k:v for k,v in model.items() if k not in {"indices","X","projection","influence","df"}}
            if model["status"] != "ready":
                for i in members: output[i] = dict(public)
                continue
            Y = values[np.ix_(members,model["indices"])].astype(float)/2
            beta = Y@model["projection"].T
            residual = Y-beta@model["X"].T
            ses = np.sqrt(np.sum((residual*model["influence"])**2,axis=1))
            crit = float(student_t.ppf(.975,model["df"]))
            for j,i in enumerate(members):
                slope, se = float(beta[j,1]), float(ses[j])
                output[i] = dict(public, status=("descriptive_adjusted_imputed_association" if self.expected_pc_count else "descriptive_no_genetic_adjustment") if se > 1e-15 else "zero_residual_variance",
                    slope_per_1000yr_toward_present=slope, se=se,
                    ci_low=slope-crit*se if se > 1e-15 else None, ci_high=slope+crit*se if se > 1e-15 else None,
                    p_two_sided=float(2*student_t.sf(abs(slope/se),model["df"])) if se > 1e-15 else None)
        return output


def allele_fdr(pvalues, mask=None):
    pvalues = np.asarray(pvalues)
    indices = np.flatnonzero(np.isfinite(pvalues) & (np.ones(len(pvalues),dtype=bool) if mask is None else mask))
    bh, by = np.full(len(pvalues),np.nan), np.full(len(pvalues),np.nan)
    if not len(indices): return bh,by
    order = indices[np.argsort(pvalues[indices],kind="stable")]
    q = pvalues[order]*len(order)/np.arange(1,len(order)+1)
    q = np.minimum(1,np.minimum.accumulate(q[::-1])[::-1])
    bh[order] = q
    by[order] = np.minimum(1,q*np.sum(1/np.arange(1,len(order)+1)))
    return bh,by


def allele_aggregate_key(entry):
    return entry["family"],entry["group"],entry["partition"],entry["role"]


def allele_block_interval(blocks, bootstrap, seed):
    points = np.asarray(list(blocks.values()),dtype=float)
    if not len(points) or points[:,1].sum() == 0: return None,None,None,0
    point = float(points[:,0].sum()/points[:,1].sum())
    if len(points) < 2 or bootstrap < 20: return point,None,None,len(points)
    rng = np.random.default_rng(seed)
    ratios = []
    for _ in range(bootstrap):
        chosen = points[rng.integers(0,len(points),len(points))].sum(axis=0)
        ratios.append(chosen[0]/chosen[1])
    low,high = np.quantile(ratios,[.025,.975])
    return point,float(low),float(high),len(points)


@transactional("manifest.json")
def allele_run(args):
    if args.bin_years < 1 or args.max_age_bp < args.bin_years or args.min_samples < 2 or not 0 < args.min_call_rate <= 1:
        raise ValueError("Invalid chronology or callability settings")
    if not 1 <= args.chunk_size <= 1024 or args.pcs < 0 or args.block_bp < 100000:
        raise ValueError("Require chunk-size 1..1024, nonnegative PCs, and blocks >=100 kb")
    out = Path(args.out).resolve(); out.mkdir(parents=True,exist_ok=True)
    inputs = [args.catalog,args.metadata,args.hq_variants,args.coordinates,args.fasta38,args.fasta37,args.chain38to37,args.chain37to38]
    inputs += [str(args.pgen_prefix)+s for s in (".pgen",".pvar",".pvar.zst",".psam")]
    owned = {out/p for p in (*ALLELE_OWNED,"manifest.json")}
    if any(Path(p).resolve() in owned for p in inputs if p): raise ValueError("Inputs and generated time outputs must not overlap")
    dump_json(out/"manifest.json",dict(status="running",version=ALLELE_VERSION))
    for name in ALLELE_OWNED: (out/name).unlink(missing_ok=True)
    prefix = str(args.pgen_prefix)
    samples, extra_meta = allele_load_samples(args.metadata,prefix+".psam",args.bin_years,args.max_age_bp)
    groups = allele_sample_groups(samples,args.bin_years,args.max_age_bp,args.min_samples)
    from pca import frame_read_scores as read_scores, frame_historical_identity as historical_identity, frame_coordinate_join as coordinate_join
    unavailable = None
    try:
        coords, pc_names, coord_audit = read_scores(args.coordinates, args.pcs,
            basis_dir=getattr(args, 'basis_dir', None), receipt=getattr(args, 'coordinate_receipt', None))
        metadata_provenance = historical_identity(args.metadata, getattr(args, 'provenance_manifest', None), coord_audit.get('metadata_recorded'))
        if args.pcs and metadata_provenance['historical_metadata_verification'] != 'verified':
            raise ValueError('Metadata lacks trusted full historical SHA; adjusted time model unavailable')
    except (ValueError, OSError) as exc:
        unavailable = str(exc)
        coords, pc_names = {}, []
        coord_audit = dict(status='unavailable', reason=unavailable, expected_pc_count=args.pcs, actual_pc_count=0,
                           coordinate_kind='unavailable', basis_id=None)
        metadata_provenance = {'historical_metadata_verification': 'unverified'}
    coord_kind = coord_audit['coordinate_kind']
    coord_audit['join'] = coordinate_join(samples, coords)
    coord_audit['metadata_provenance'] = metadata_provenance
    dump_json(out / 'coordinate_contract.json', coord_audit)
    model = allele_TrendModel(samples, groups, coords, pc_names, args.min_call_rate, args.pcs, unavailable)
    dump_json(out/"model_covariates.json",dict(model="HC3_OLS_of_imputed_ALT_fraction_ancient_only", direction="two_sided; positive increases toward present",
        coordinate_contract=coord_audit, genetic_adjustment="unavailable" if unavailable else "none_explicit_pcs_0" if args.pcs == 0 else "verified_shared_PCs", candidate_covariates=model.names, n_complete_case_samples=int(model.eligible.sum()), coordinate_kind=coord_kind,
        coverage_adjusted=model.use_coverage, region_adjusted=model.use_region, source_adjusted=model.use_source,
        missing_requested_covariates="complete_case_exclusion", cache_max_design_masks=4))
    targets = allele_load_catalog(args.catalog)
    with GenomeMapper(args.fasta38,args.fasta37,args.chain38to37,args.chain37to38) as mapper:
        targets = map_targets(targets,mapper)
    map_fields = sorted(set().union(*(r.keys() for r in targets)) - {"_counts"}) if targets else ["site_id","mapping_status"]
    write_rows(out/"mapped_catalog.tsv.gz",targets,map_fields)
    pvar = panel_pvar(prefix)
    availability = index_panel(targets,pvar,args.hq_variants)
    selected,audit = allele_catalog_contrasts(targets,availability)
    dump_json(out/"measurement_plans.json", dict(schema="nucleotide-count-1", ALT_measurements=selected,
        unknown_split_reference="unavailable; no complement is relabelled REF"))
    del availability
    write_rows(out/"all_allele_audit.tsv.gz",audit,ALLELE_AUDIT_FIELDS)
    available_counts = Counter((r["family"],r["group"],r["partition"],r["role"],r["strict_filter_complete"],r["panel_status"]) for r in audit)
    write_rows(out/"panel_availability_counts.tsv",(dict(zip(["family","group","partition","role","strict_filter_complete","panel_status"],k),n_alleles=v) for k,v in sorted(available_counts.items())),
                  ["family","group","partition","role","strict_filter_complete","panel_status","n_alleles"])
    by_site = defaultdict(list)
    for entry in selected: by_site[entry["site_id"]].append(entry)
    coverage = Counter()
    for t in targets:
        for level,family,group,partition in [("catalog","all","all","all"),("family",t["family"],t["group"],"all"),("partition",t["family"],t["group"],t["partition"])]:
            key = level,family,group,partition
            coverage[(*key,"catalog_sites")] += 1
            coverage[(*key,"mapped_sites")] += t.get("mapping_status") == "mapped"
            coverage[(*key,"any_HQ_exact_ALT_sites")] += bool(by_site[t["site_id"]])
            coverage[(*key,"modern_ALT_bases")] += len(t["alts"].split(","))
            coverage[(*key,"represented_HQ_ALT_contrasts")] += len(by_site[t["site_id"]])
    write_rows(out/"catalog_coverage.tsv",(dict(zip(["level","family","group","partition","metric"],k),value=v) for k,v in sorted(coverage.items())),
                  ["level","family","group","partition","metric","value"])
    del audit,by_site
    n_samples,n_variants = len(samples),len(selected)
    fixed_flags = np.zeros(n_variants, dtype=bool)
    family_expected = Counter()
    sample_called = np.zeros(n_samples,dtype=np.int64); sample_fixed = np.zeros(n_samples,dtype=np.int64)
    pvalues = np.full(n_variants,np.nan); primary_tests = np.zeros(n_variants,dtype=bool)
    fixed_count, variant_count, max_rows = 0,None,0
    measurement_qc = Counter()
    fixed_epochs = [i for i,g in enumerate(groups) if g["qualified_epoch"]]
    primary_epochs = [i for i in fixed_epochs if groups[i]["age_scope"] == "primary_ancient"]
    pooled = list(fixed_epochs)
    family_group_indices = {i for i,g in enumerate(groups) if g["age_scope"] in {"primary_ancient", "modern"}}
    family_blocks = defaultdict(lambda: defaultdict(lambda: np.zeros(4)))
    pair_records = defaultdict(lambda: {"case":[],"control":[]})
    expected_pairs = defaultdict(Counter)
    for target in targets:
        for base in target["alts"].split(","):
            expected_pairs[(target["match_id"],allele_alt_role(target,base))][target["group"]] += 1
    call_fields = ["site_id","contrast_id","role","fixed_panel","min_qualified_epoch_call_rate","n_called_all_PSAM","n_PSAM","strict_filter_complete"]
    with allele_table(out/"allele_trajectories.tsv.gz",ALLELE_TRAJECTORY_FIELDS) as trajectory, \
         allele_table(out/"allele_callability.tsv.gz",call_fields) as call_audit, \
         allele_table(out/"trends.raw.tsv.gz",ALLELE_TREND_FIELDS) as trend:
        for start,values,variant_count in iter_pgen_chunks(prefix+".pgen",n_samples,selected,args.chunk_size,pvar_path=pvar,audit=measurement_qc) if selected else []:
            max_rows = max(max_rows,len(values))
            called = values >= 0
            sample_called += called.sum(axis=0)
            stats = []
            for group in groups:
                data = values[:,group["members"]]
                valid = data >= 0
                n = valid.sum(axis=1)
                ac = np.where(valid,data,0).sum(axis=1)
                af = np.divide(ac,2*n,out=np.full(len(values),np.nan),where=n>0)
                het = np.divide((data == 1).sum(axis=1),n,out=np.full(len(values),np.nan),where=n>0)
                stats.append((n,ac,af,het,n/group["n_eligible"]))
            min_rates = np.min(np.stack([stats[i][4] for i in fixed_epochs]),axis=0) if fixed_epochs else np.zeros(len(values))
            keep = (min_rates >= args.min_call_rate) & (len(primary_epochs) >= 2)
            sample_fixed += called[keep].sum(axis=0)
            fixed_count += int(keep.sum())
            fixed_flags[start:start+len(values)] = keep
            strict = np.array([bool(e["strict_filter_complete"]) for e in selected[start:start+len(values)]])
            fits = model.fit_chunk(values,keep & strict)
            for local,entry in enumerate(selected[start:start+len(values)]):
                index = start+local
                basic = dict(site_id=entry["site_id"],contrast_id=entry["contrast_id"],family=entry["family"],group=entry["group"],
                    partition=entry["partition"],match_id=entry["match_id"],role=entry["role"],counted_base38=entry["allele"],
                    paired_modern_alt=entry["allele"],modern_ac=entry["modern_ac"],modern_ac5=entry["modern_ac5"],modern_ac10=entry["modern_ac10"],
                    strict_filter_complete=entry["strict_filter_complete"],fixed_panel=int(keep[local]))
                basic.update({k:entry[k] for k in ALLELE_LOCUS_FIELDS})
                basic.update({k:entry.get(k) for k in ("allele_id","counted_base37","pvar_index","allele_idx",
                    "measurement_mode","counting_scope","reference_count_scope")})
                call_audit(dict(basic,min_qualified_epoch_call_rate=float(min_rates[local]),n_called_all_PSAM=int(called[local].sum()),n_PSAM=n_samples))
                for group_index,group in enumerate(groups):
                    n,ac,af,het,rate = [s[local] for s in stats[group_index]]
                    record = dict(basic,**{k:v for k,v in group.items() if k not in {"members","qualified_epoch"}},n_called=int(n),
                        ac=int(ac) if n else None,an=2*int(n),af=float(af) if n else None,hardcall_heterozygosity=float(het) if n else None,
                        call_rate=float(rate),status="imputed_nucleotide_count" if n else "no_called_genotypes")
                    trajectory(record)
                    if keep[local] and strict[local] and group_index in family_group_indices and n:
                        block = entry["chrom"]+":"+str((int(entry["pos"])-1)//args.block_bp)
                        keys = [allele_aggregate_key(entry),(entry["family"],entry["group"],"all_partitions",entry["role"])]
                        for key in keys:
                            family_blocks[(*key,group_index)][block] += [af,1,het,2*n]
                if keep[local] and strict[local]:
                    for key in (allele_aggregate_key(entry), (entry['family'],entry['group'],'all_partitions',entry['role'])):
                        family_expected[key] += 1
                fit = fits[local]
                if fit.get("p_two_sided") is not None: pvalues[index] = fit["p_two_sided"]
                primary_tests[index] = bool(entry["modern_ac5"] and strict[local])
                trend(dict(basic,allele_index=index,**fit))
                if keep[local] and strict[local]:
                    pair_records[(entry["match_id"],entry["role"])][entry["group"]].append((entry,
                        np.asarray([[stats[g][2][local],stats[g][3][local]] for g in pooled],dtype=np.float32)))
            if start == 0 or (start+len(values)) % 10000 < args.chunk_size:
                print(f"[allele-time] {start+len(values):,}/{n_variants:,} measured ALT nucleotides; all {n_samples:,} PSAM samples; genotype buffer <= {args.chunk_size} rows",flush=True)
    qbh,qby = allele_fdr(pvalues); pbh,pby = allele_fdr(pvalues,primary_tests)
    def adjusted_rows():
        for row in rows(out/"trends.raw.tsv.gz"):
            i = int(row["allele_index"])
            for name,values_ in (("q_bh_all_alt",qbh),("q_by_all_alt",qby),("q_bh_primary_ac5",pbh),("q_by_primary_ac5",pby)):
                row[name] = float(values_[i]) if np.isfinite(values_[i]) else None
            yield row
    write_rows(out/"adjusted_allele_trends.tsv.gz",adjusted_rows(),ALLELE_TREND_FIELDS)
    (out/"trends.raw.tsv.gz").unlink()
    family_rows = []
    for (family,group,partition,role,gi),blocks in sorted(family_blocks.items()):
        mean,lo,hi,n_blocks = allele_block_interval({b:v[:2] for b,v in blocks.items()},args.bootstrap,args.seed)
        summed = np.sum(list(blocks.values()),axis=0)
        g = groups[gi]
        family_rows.append(dict(family=family,group=group,partition=partition,role=role,region=g["region"],measurement_source=g.get("measurement_source","all"),expected_fixed_contrasts=family_expected[family,group,partition,role],cohort=g["cohort"],age_scope=g["age_scope"],time_bin=g["time_bin"],
            start_bp=g["start_bp"],end_bp=g["end_bp"],n_metadata=g["n_eligible"],n_fixed_contrasts=int(summed[1]),mean_AF=mean,ci_low=lo,ci_high=hi,
            mean_heterozygosity=float(summed[2]/summed[1]),summed_effective_AN=int(summed[3]),n_blocks=n_blocks,
            status="fixed_panel_descriptive" if n_blocks >= 20 else "few_blocks_descriptive"))
    write_rows(out/"fixed_panel_family_trajectories.tsv",[r for r in family_rows if r['region']=='all' and r['measurement_source']=='all'],list(family_rows[0]) if family_rows else ["status"])
    write_rows(out/"regional_fixed_panel_family_trajectories.tsv",family_rows,list(family_rows[0]) if family_rows else ["status"])
    joint_family_composition(family_rows, samples, args)
    pair_blocks = defaultdict(lambda:defaultdict(lambda:np.zeros(2)))
    match_audit = []
    window_matching = allele_window_matched_catalog(targets)
    for (match_id,role),expected in expected_pairs.items():
        pair = pair_records.get((match_id,role),{"case":[],"control":[]})
        cases,controls = pair["case"],pair["control"]
        valid = len(cases) == expected["case"] == 1 and len(controls) == expected["control"] and bool(controls)
        same_chromosome = bool(cases) and all(c["chrom"] == cases[0][0]["chrom"] for c,_ in controls)
        valid &= same_chromosome and not window_matching
        match_audit.append(dict(match_id=match_id,role=role,n_case_contrasts=len(cases),n_control_contrasts=len(controls),
            expected_case_contrasts=expected["case"],expected_control_contrasts=expected["control"],status="paired_fixed_loci" if valid else "missing_or_nonunique_exact_role"))
        if not valid: continue
        case,case_values = cases[0]
        difference = case_values-np.mean([v for _,v in controls],axis=0)
        # Controls may be far from their matched case. Whole chromosomes retain
        # dependence on both sides; a case-only 5 Mb block would not do so.
        block = case["chrom"]
        for j,gi in enumerate(pooled):
            for m,measure in enumerate(("pair_conditional_AF_difference","biallelic_heterozygosity_difference")):
                if np.isfinite(difference[j,m]):
                    pair_blocks[(case["family"],case["partition"],role,measure,gi)][block] += [difference[j,m],1]
    matched_rows = []
    for (family,partition,role,measure,gi),blocks in sorted(pair_blocks.items()):
        point,lo,hi,n_blocks = allele_block_interval(blocks,args.bootstrap,args.seed)
        g = groups[gi]
        matched_rows.append(dict(family=family,partition=partition,role=role,measure=measure,cohort=g["cohort"],age_scope=g["age_scope"],time_bin=g["time_bin"],
            start_bp=g["start_bp"],end_bp=g["end_bp"],n_matched_loci=int(sum(v[1] for v in blocks.values())),n_blocks=n_blocks,difference=point,ci_low=lo,ci_high=hi,
            denominator="locus/role matched group estimates; not same-person paired complete calls",
            bootstrap_unit="whole_chromosome_case_and_controls_together",status="descriptive_chromosome_bootstrap"))
    write_rows(out/"matched_locus_role_trajectories.tsv",matched_rows,list(matched_rows[0]) if matched_rows else ["status"])
    write_rows(out/"matched_locus_audit.tsv",match_audit,["match_id","role","n_case_contrasts","n_control_contrasts","expected_case_contrasts","expected_control_contrasts","status"])
    sample_fields = ["sample_id","individual_id","cohort","age_bp","age_scope","unrelated","include","region","time_bin","data_source","target_coverage",
                     "n_HQ_exact_contrasts","n_called_HQ_contrasts","call_rate_HQ","n_fixed_contrasts","n_called_fixed_contrasts","model_covariate_complete"]
    write_rows(out/"sample_coverage.tsv.gz",(dict(s,n_HQ_exact_contrasts=n_variants,n_called_HQ_contrasts=int(sample_called[i]),
        call_rate_HQ=float(sample_called[i]/n_variants) if n_variants else None,n_fixed_contrasts=fixed_count,
        n_called_fixed_contrasts=int(sample_fixed[i]),model_covariate_complete=int(model.eligible[i])) for i,s in enumerate(samples)),sample_fields)
    write_rows(out/"sample_time_region_coverage.tsv",(dict({k:v for k,v in g.items() if k != "members"},
        n_model_covariate_complete=int(model.eligible[g["members"]].sum()),summed_called_HQ_contrasts=int(sample_called[g["members"]].sum())) for g in groups),
        ["cohort","age_scope","time_bin","region","start_bp","end_bp","n_eligible","qualified_epoch","n_model_covariate_complete","summed_called_HQ_contrasts"])
    joint = joint_run_joint(args, targets, selected, samples, coord_audit, coords, fixed_flags)
    if not args.no_plots: allele_figures([r for r in family_rows if r['region']=='all' and r['measurement_source']=='all'],available_counts,out)
    catalog_receipt = Path(args.catalog).parent/"candidate_manifest.json"
    if not catalog_receipt.is_file() and allele_window_matched_catalog(targets):
        catalog_receipt = Path(args.catalog).parent/"manifest.json"
    catalog_provenance = json.loads(catalog_receipt.read_text()) if catalog_receipt.is_file() else {}
    audit_status_counts = Counter()
    for key,value in available_counts.items(): audit_status_counts[key[5]] += value
    result = dict(joint_analysis=joint, coordinate_contract=coord_audit, genetic_adjustment="unavailable" if unavailable else "none_explicit_pcs_0" if args.pcs == 0 else "verified_shared_PCs", version=ALLELE_VERSION,status=("complete_adjustment_unavailable" if unavailable else "complete_descriptive_all_allele_panel") if selected else "no_usable_exact_high_quality_contrast",
        n_PSAM_samples=n_samples,sample_cohorts=dict(Counter(s["cohort"] for s in samples)),sample_age_scopes=dict(Counter(s["age_scope"] for s in samples)),
        sample_unrelated_flags=dict(Counter(str(s["unrelated"]) for s in samples)),
        metadata_extra_ids=extra_meta,n_model_covariate_complete=int(model.eligible.sum()),catalog_sites=len(targets),catalog_ALT_bases=sum(len(t["alts"].split(",")) for t in targets),
        exact_HQ_ALT_contrasts=n_variants,n_fixed_contrasts=fixed_count,n_tested_ALT_contrasts=int(np.isfinite(pvalues).sum()),
        n_tested_primary_strict_AC5=int((np.isfinite(pvalues)&primary_tests).sum()),original_PGEN_variant_count=variant_count,
        maximum_genotype_chunk_rows=max_rows,maximum_genotype_buffer_bytes=max_rows*n_samples,
        buffer_accounting="Decoded int8 genotype buffer only; models use additional bounded O(chunk_size * samples) numeric work arrays",
        genotype_storage="int8 chunk by full PSAM; no full allele-by-sample matrix or memmap",sample_subset="none; every PSAM individual audited/read",
        count_kind="imputed_diploid_nucleotide_hardcall",measurement_schema="nucleotide-count-1",
        measurement_qc=dict(measurement_qc),
        REF_rule="REF availability is audited only; no single-ALT complement is emitted as a REF trajectory",
        contrast_id_definition="Modern catalog REF/ALT identity; the actual count is specified by counting_scope and allele_idx",
        heterozygosity_definition="fraction with counted nucleotide dosage exactly one; target versus all other nucleotides, not full-locus heterozygosity; existing biallelic_heterozygosity_difference identifier retained",
        matched_AF_definition="case minus control counted-nucleotide frequency, conditional on represented calls; existing pair_conditional_AF_difference identifier retained",
        main_time_scope=f"unrelated=1, include=1, 0<=age_bp<={args.max_age_bp}; qualified pooled epochs; complete covariates",
        allele_audit_counts=dict(audit_status_counts),coordinate_kind=coord_kind,
        qualified_epochs=[{k:g[k] for k in ("cohort","age_scope","time_bin","n_eligible")} for g in groups if g["qualified_epoch"]],
        fixed_panel_rule="Same contrasts pass call-rate in every qualified unrelated ancient epoch and qualified modern endpoint; sparse epochs are individual-allele descriptions only",
        catalog_source_manifest=str(catalog_receipt) if catalog_provenance else None,
        catalog_source_scope=("sampled_reported_modern_ALT_states_in_frozen_genome_windows; not genome-wide ancient allele frequencies"
            if allele_window_matched_catalog(targets) else catalog_provenance.get("scope",catalog_provenance.get("case_selection","supplied catalog; scope must be checked"))),
        locus_pairing="disabled_window_matches_are_not_site_pairs" if allele_window_matched_catalog(targets) else "explicit_locus_role_matches",
        outputs=[name for name in ALLELE_OWNED if (out/name).is_file()],
        limitations=ALLELE_LIMITATIONS,parameters={k:v for k,v in vars(args).items()})
    dump_json(out/"manifest.json",result)
    return result


def allele_figures(rows,availability,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usable = [r for r in rows if r["partition"] == "all_partitions"]
    if usable:
        measured_epochs = {r["time_bin"] for r in usable if r["cohort"] != "modern"}
        epochs = list(range(max(measured_epochs),min(measured_epochs)-1,-1)) if measured_epochs else []
        ranges = {r["time_bin"]:(r["start_bp"],r["end_bp"]) for r in usable if r["cohort"] != "modern"}
        width = max((hi-lo for lo,hi in ranges.values()),default=1000)
        labels = [f"{ranges.get(e,(e*width,(e+1)*width))[0]/1000:g}–{ranges.get(e,(e*width,(e+1)*width))[1]/1000:g}k" for e in epochs]
        index = {e:i for i,e in enumerate(epochs)}
        has_modern = any(r["cohort"] == "modern" for r in usable)
        if has_modern:
            # A blank column/separator prevents a modern cohort being shown as
            # the next equally spaced ancient date or a continuation of its line.
            index["modern"] = len(labels)+1
            labels += ["", "modern"]
        families = sorted({r["family"] for r in usable})
        fig,axes = plt.subplots(len(families),1,figsize=(11,4*len(families)),squeeze=False)
        for ax,family in zip(axes[:,0],families):
            for group,role in sorted({(r["group"],r["role"]) for r in usable if r["family"] == family}):
                data = sorted([r for r in usable if r["family"] == family and r["group"] == group and r["role"] == role],key=lambda r:index[r["time_bin"]])
                ancient = [r for r in data if r["cohort"] != "modern"]
                current = [r for r in data if r["cohort"] == "modern"]
                line = ax.plot([index[r["time_bin"]] for r in ancient],[r["mean_AF"] for r in ancient],marker="o",ls="-" if group == "case" else "--",label=f"{group}: {role}")[0]
                if current:
                    ax.scatter([index["modern"]]*len(current),[r["mean_AF"] for r in current],color=line.get_color(),marker="D",s=30,zorder=3)
                for row in data:
                    if row["ci_low"] is not None and row["ci_high"] is not None:
                        x = index[row["time_bin"]]
                        ax.vlines(x,row["ci_low"],row["ci_high"],color=line.get_color(),alpha=.65,lw=1)
                        ax.hlines([row["ci_low"],row["ci_high"]],x-.07,x+.07,color=line.get_color(),alpha=.65,lw=1)
            if has_modern: ax.axvline(len(epochs),color="grey",ls=":",lw=1)
            ax.set(xticks=range(len(labels)),xticklabels=labels,ylabel="Mean counted-nucleotide ALT frequency",title=family)
            ax.legend(fontsize=8,ncol=2)
        fig.suptitle("Fixed modern-catalog imputed contrasts · modern is a separate endpoint")
        fig.tight_layout()
        for ext in ("png","pdf"): fig.savefig(out/f"family_allele_time.{ext}",dpi=170)
        plt.close(fig)
        case_rows = [r for r in rows if r["group"] == "case" and r["partition"] != "all_partitions"]
        keys = sorted({(r["family"],r["partition"],r["role"]) for r in case_rows})
        if keys:
            matrix = np.full((len(keys),len(labels)),np.nan); lookup = {k:i for i,k in enumerate(keys)}
            for r in case_rows: matrix[lookup[r["family"],r["partition"],r["role"]],index[r["time_bin"]]] = r["mean_AF"]
            fig,ax = plt.subplots(figsize=(11,max(4,.25*len(keys))))
            im = ax.imshow(matrix,aspect="auto",vmin=0,vmax=1,cmap="viridis")
            ax.set(xticks=range(len(labels)),xticklabels=labels,yticks=range(len(keys)),yticklabels=[" / ".join(k) for k in keys],title="All fixed-panel family/partition/ALT roles; no peak-based selection")
            ax.tick_params(axis="y",labelsize=7); fig.colorbar(im,ax=ax,label="Mean counted-nucleotide AF"); fig.tight_layout()
            for ext in ("png","pdf"): fig.savefig(out/f"allele_time_heatmap.{ext}",dpi=170)
            plt.close(fig)
    counts = Counter()
    for key,n in availability.items():
        if key[1] == "case" and key[3] != "REF": counts[key[3],key[5] == "assayed_pass"] += n
    roles = sorted({k[0] for k in counts})
    if roles:
        fig,ax = plt.subplots(figsize=(9,4))
        measured = [counts[r,True] for r in roles]; absent = [counts[r,False] for r in roles]
        ax.bar(roles,measured,label="HQ encoded nucleotide"); ax.bar(roles,absent,bottom=measured,label="Unavailable / unsupported")
        ax.set(ylabel="Modern ALT identities",title="Full catalog panel availability: unavailable never means AF=0"); ax.legend(); fig.tight_layout()
        for ext in ("png","pdf"): fig.savefig(out/f"allele_panel_availability.{ext}",dpi=170)
        plt.close(fig)


def allele_parser():
    p = argparse.ArgumentParser(description=_description_herv_allele_time)
    p.add_argument("command",choices=["run"],nargs="?",default="run")
    p.add_argument("--catalog",required=True)
    for name in ("pgen-prefix","metadata","hq-variants","fasta38","fasta37","chain38to37","chain37to38","out"):
        p.add_argument("--"+name,required=True)
    p.add_argument("--coordinates")
    for name in ("basis-dir", "coordinate-receipt", "provenance-manifest"):
        p.add_argument("--" + name)
    for name,value in (("bin-years",1000),("max-age-bp",10000),("min-samples",20),("pcs",10),("chunk-size",64),("block-bp",5000000),("bootstrap",200),("seed",20261028)):
        p.add_argument("--"+name,type=int,default=value)
    p.add_argument("--min-call-rate",type=float,default=.9)
    p.add_argument("--no-plots",action="store_true")
    return p


def allele_main(argv=None):
    allele_run(allele_parser().parse_args(argv))
    return 0


def _entry_herv_allele_time():
    try: raise SystemExit(allele_main())
    except (ValueError,OSError,RuntimeError) as allele_exc: raise SystemExit(f"herv-allele-time: {allele_exc}")


# ========================================================================
# herv selection
# ========================================================================
_description_herv_selection = "Attach published selection statistics to exact modern-catalog allele labels.\n\nThis optional, offline command does not fit a selection model. It independently\nmaps the complete supplied GRCh38 catalog, streams the published GRCh37 table,\nand uses a disk index to bound memory. Raw statistics retain their original\nREF/ALT contrast. Signed allele effects are unavailable unless the caller\nexplicitly declares that the source table's S and X refer to ALT."
#!/usr/bin/env python3
"""Attach published selection statistics to exact modern-catalog allele labels.

This optional, offline command does not fit a selection model. It independently
maps the complete supplied GRCh38 catalog, streams the published GRCh37 table,
and uses a disk index to bound memory. Raw statistics retain their original
REF/ALT contrast. Signed allele effects are unavailable unless the caller
explicitly declares that the source table's S and X refer to ALT.
"""
import hashlib
import sqlite3
from allele_common import fingerprint

SELECTION_VERSION = "herv-published-selection-1.0"
SELECTION_AUTOSOMES = frozenset(str(i) for i in range(1, 23))
SELECTION_SOURCE = "https://doi.org/10.1038/s41586-026-10358-1"
SELECTION_DATASET = "https://doi.org/10.7910/DVN/7RVV9N"
SELECTION_PUBLISHED = ("CHROM", "POS", "REF", "ALT", "ANC", "ID", "RSID", "AF", "S", "SE",
             "X", "P_X", "POSTERIOR", "FDR", "CHI2_BE", "FILTER")
SELECTION_BASE_FIELDS = ["site_id", "chrom", "pos", "ref", "counted_base38", "role", "modern_ac", "modern_an",
               "modern_af", "family", "group", "partition", "comparison_partition", "match_id",
               "strict_filter_complete", "mapping_status", "chrom37", "pos37", "ref37", "counted_base37",
               "strand37to38", "reference_changed"]
SELECTION_MATCH_FIELDS = SELECTION_BASE_FIELDS + ["published_contrast", "counted_base_relation", "statistics_status",
                             "published_effect_contract", "oriented_s", "oriented_x", "source_doi"] + list(SELECTION_PUBLISHED)
SELECTION_AVAIL_FIELDS = SELECTION_BASE_FIELDS + ["availability_status", "n_published_contrasts", "n_accepted_filter_contrasts",
                             "n_valid_accepted_filter_contrasts", "n_FDR_le_0_01", "n_FDR_le_0_05", "n_FDR_le_0_10",
                             "n_FDR_le_0_50", "frequency_zero_inferred"]
SELECTION_SUMMARY_FIELDS = ["family", "group", "partition", "comparison_partition", "role", "catalog_bases",
                  "mapped_bases", "matched_bases", "accepted_filter_bases", "valid_accepted_filter_bases", "FDR_le_0_01_bases",
                  "FDR_le_0_05_bases", "FDR_le_0_10_bases", "FDR_le_0_50_bases", "matching_fraction",
                  "FDR_le_0_05_fraction_of_valid_accepted_filter_bases", "unit"]
SELECTION_OWNED = ("catalog_index.sqlite", "published_selection_alleles.tsv.gz", "availability.tsv",
         "published_selection_summary.tsv", "manifest.json")
SELECTION_LIMITATIONS = [
    "Published Akbari 2026 statistics are reused; no new GLMM, selection coefficient or independent replication is inferred.",
    "The supplied modern catalog defines ascertainment. A missing published contrast is unavailable, never evidence of ancient AF=0.",
    "Every statistic belongs to its original exact biallelic REF/ALT contrast. It is not transferred to a third base at the same position.",
    "REF and ALT views of one published contrast are the same evidence, not independent tests. REF labels are not assumed ancestral.",
    "The modern REF is a coordinate label, not automatically an observed gnomAD allele; its modern AC/AN are left unavailable.",
    "Descriptive FDR summaries count catalog base identities with at least one numerically valid contrast whose FILTER is PASS or '.'. A dot is accepted as unflagged, not proof of a separate QC pass. These are not independent loci and no enrichment P value is computed.",
    "Published FDR, posterior and calibrated X are preserved. S/SE is not substituted for calibrated X and P_X is not recalculated.",
    "HERV mapping, frequency, context, imputation coverage, LD and population structure can confound cross-family comparisons.",
    "The default raw_only contract leaves oriented_s and oriented_x unavailable. ALT orientation is only enabled by an explicit caller declaration.",
]


def selection_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_code_identity():
    return {name: selection_sha256(Path(__file__).with_name(name))
            for name in ('chronology.py', 'chronology.py', "0.common.py")}


def selection_identity(args):
    paths = {key: getattr(args, key) for key in
             ("catalog", "selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38")}
    inputs = {key: fingerprint(path) for key, path in paths.items()}
    inputs["catalog"]["sha256"] = selection_sha256(args.catalog)
    for key in ("fasta38", "fasta37"):
        index = Path(str(paths[key]) + ".fai")
        if index.is_file():
            inputs[key + "_fai"] = fingerprint(index)
        compressed_index = Path(str(paths[key]) + ".gzi")
        if str(paths[key]).endswith((".gz", ".bgz")) and compressed_index.is_file():
            inputs[key + "_gzi"] = fingerprint(compressed_index)
    return dict(version=SELECTION_VERSION, code=selection_code_identity(), inputs=inputs, effect_allele=args.effect_allele,
                mapping="independent reciprocal GRCh38-to-GRCh37 with FASTA checks")


@contextlib.contextmanager
def selection_writer(path, fields):
    path = Path(path)
    partial = path.with_name(path.name + ".partial")
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(partial, "wt", newline="") as handle:
        output = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        output.writeheader()
        def emit(row):
            output.writerow({key: "NA" if row.get(key) is None else row.get(key, ".") for key in fields})
        yield emit
    os.replace(partial, path)


def selection_catalog_rows(path):
    """Read one catalog site at a time; SQLite verifies global uniqueness."""
    required = ("chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id")
    for raw in rows(path, required):
        row = dict(raw)
        row["chrom"], row["pos"], row["ref"] = canonical_chrom(row["chrom"]), int(row["pos"]), row["ref"].upper()
        alts_text = row.get("observed_alts_all", row["alts"])
        if alts_text in NA or "all_alt_ac" in row:
            alts_text = row["alts"]
        alts = alts_text.upper().split(",")
        if (row["chrom"] not in SELECTION_AUTOSOMES or row["pos"] < 1 or row["ref"] not in BASES
                or not alts or len(set(alts)) != len(alts) or any(b not in BASES or b == row["ref"] for b in alts)):
            raise ValueError(f"Catalog requires canonical autosomal SNVs: {row['site_id']}")
        if row["site_id"] in NA or row["match_id"] in NA or row["group"] not in {"case", "control"}:
            raise ValueError("Catalog needs nonmissing site_id/match_id and case/control group")
        counts = {}
        ac_text, an_text = row.get("all_alt_ac", row.get("alt_ac_all", ".")), row.get("all_alt_an", row.get("alt_an_all", "."))
        if (ac_text in NA) != (an_text in NA):
            raise ValueError("Catalog AC and AN must either both be present or both be unavailable")
        if ac_text not in NA:
            ac, an = [int(x) for x in ac_text.split(",")], [int(x) for x in an_text.split(",")]
            if len(ac) != len(alts) or len(an) != len(alts) or any(a < 0 or n < 1 or a > n for a, n in zip(ac, an)):
                raise ValueError(f"Catalog ALT AC/AN order or bounds invalid: {row['site_id']}")
            counts = {base: (a, n) for base, a, n in zip(alts, ac, an)}
        for key in ("primary_alt", "secondary_alt"):
            if row.get(key, ".") not in NA and row[key] not in alts:
                raise ValueError(f"{key} is not a catalog ALT: {row['site_id']}")
        row.update(alts=",".join(alts), _counts=counts,
                   partition=row.get("partition", row.get("case_region", "unclassified")),
                   comparison_partition=row.get("comparison_partition", row.get("partition", "unclassified")))
        yield row


def selection_base_info(target, base):
    ac, an = target["_counts"].get(base, (None, None))
    role = ("REF" if base == target["ref"] else "primary" if base == target.get("primary_alt")
            else "secondary" if base == target.get("secondary_alt")
            else "additional_lowcount" if ac is not None and ac < 5 else "additional")
    row = {key: target.get(key) for key in SELECTION_BASE_FIELDS}
    row.update(counted_base38=base, role=role, modern_ac=ac, modern_an=an, modern_af=ac / an if an else None,
               strict_filter_complete=int(str(target.get("strict_filter_complete", "0")) == "1"))
    row["counted_base37"] = orient(base, target["strand37to38"]) if target["mapping_status"] == "mapped" else None
    return row


def selection_build_index(path, catalog, mapper, qc):
    db = sqlite3.connect(path)
    db.execute("PRAGMA journal_mode=DELETE")
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA cache_size=-32768")
    db.executescript("""
      CREATE TABLE sites(site_id TEXT PRIMARY KEY, chrom TEXT, pos INTEGER, chrom37 TEXT, pos37 INTEGER, payload TEXT,
                         UNIQUE(chrom,pos));
      CREATE UNIQUE INDEX mapped_coordinate ON sites(chrom37,pos37);
      CREATE TABLE alleles(site_id TEXT, base38 TEXT, payload TEXT, PRIMARY KEY(site_id,base38));
      CREATE TABLE contrasts(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,PRIMARY KEY(chrom,pos,ref,alt));
      CREATE TABLE matches(site_id TEXT,base38 TEXT,contrast TEXT,passed INTEGER,valid INTEGER,fdr REAL,
                           PRIMARY KEY(site_id,base38,contrast));
    """)
    try:
        for target in selection_catalog_rows(catalog):
            target = map_targets([target], mapper)[0]
            mapped = target["mapping_status"] == "mapped"
            qc["catalog_sites"] += 1
            qc["mapping_" + target["mapping_status"]] += 1
            db.execute("INSERT INTO sites VALUES(?,?,?,?,?,?)", (target["site_id"], target["chrom"], target["pos"],
                       target.get("chrom37") if mapped else None, target.get("pos37") if mapped else None, json.dumps(target)))
            for base in [target["ref"], *target["alts"].split(",")]:
                info = selection_base_info(target, base)
                db.execute("INSERT INTO alleles VALUES(?,?,?)", (target["site_id"], base, json.dumps(info)))
                qc["catalog_bases"] += 1
            if qc["catalog_sites"] % 10000 == 0:
                db.commit()
        if not qc["catalog_sites"]:
            raise ValueError("Empty catalog")
        db.commit()
        return db
    except BaseException:
        db.close()
        raise


def selection_statistics_valid(row):
    # CHI2_BE is retained verbatim; missing batch diagnostics do not erase the
    # published inferential fields, and no new batch threshold is invented.
    values = {key: finite(row[key]) for key in ("AF", "S", "SE", "X", "P_X", "POSTERIOR", "FDR")}
    return (all(v is not None for v in values.values()) and values["SE"] > 0
            and all(0 <= values[k] <= 1 for k in ("AF", "P_X", "POSTERIOR", "FDR")))


def selection_join_statistics(db, path, emit, effect_allele, qc):
    for row in rows(path, SELECTION_PUBLISHED):
        qc["published_rows_streamed"] += 1
        if qc["published_rows_streamed"] % 1000000 == 0:
            print(f"[published-selection] {qc['published_rows_streamed']:,} source rows streamed", flush=True)
        chrom = canonical_chrom(row["CHROM"])
        try:
            pos = int(row["POS"])
        except ValueError as exc:
            raise ValueError("Published POS must be an integer") from exc
        ref, alt = row["REF"], row["ALT"]
        if chrom not in SELECTION_AUTOSOMES or pos < 1 or ref not in BASES or alt not in BASES or ref == alt:
            qc["published_noncanonical_snv_skipped"] += 1
            continue
        found = db.execute("SELECT payload FROM sites WHERE chrom37=? AND pos37=?", (chrom, pos)).fetchone()
        if found is None:
            continue
        target = json.loads(found[0])
        qc["published_rows_at_catalog_coordinate"] += 1
        if ref != target["ref37"]:
            qc["published_reference_mismatch"] += 1
            continue
        if alt not in target["alts37"].split(","):
            qc["published_alt_not_in_modern_labels"] += 1
            continue
        try:
            db.execute("INSERT INTO contrasts VALUES(?,?,?,?)", (chrom, pos, ref, alt))
        except sqlite3.IntegrityError as exc:
            raise ValueError(f"Duplicate published exact contrast at catalog site: {chrom}:{pos}:{ref}:{alt}") from exc
        key = f"{chrom}:{pos}:{ref}:{alt}"
        valid, passed = selection_statistics_valid(row), row["FILTER"] in {"PASS", "."}
        qc["matched_published_contrasts"] += 1
        qc["matched_valid_accepted_filter_contrasts"] += int(valid and passed)
        for base37, relation in ((ref, "REF"), (alt, "ALT")):
            base38 = orient(base37, target["strand37to38"])
            info = selection_base_info(target, base38)
            if base38 not in {target["ref"], *target["alts"].split(",")}:
                raise ValueError("Mapped contrast contains a base absent from the modern catalog")
            info.update(row)
            signed = effect_allele == "ALT" and valid
            sign = 1 if relation == "ALT" else -1
            info.update(published_contrast=key, counted_base_relation=relation,
                        statistics_status="valid" if valid else "invalid_numeric_fields",
                        published_effect_contract="caller_declared_ALT" if effect_allele == "ALT" else "raw_only_unverified",
                        oriented_s=sign * float(row["S"]) if signed else None,
                        oriented_x=sign * float(row["X"]) if signed else None, source_doi=SELECTION_SOURCE)
            emit(info)
            db.execute("INSERT INTO matches VALUES(?,?,?,?,?,?)", (target["site_id"], base38, key, int(passed), int(valid),
                       float(row["FDR"]) if valid else None))
            qc["matched_base_contrast_rows"] += 1
        if qc["matched_published_contrasts"] % 10000 == 0:
            db.commit()
    db.commit()


def selection_export_availability(db, out, qc):
    grouped = defaultdict(Counter)
    query = """SELECT a.payload,COUNT(m.contrast),COALESCE(SUM(m.passed),0),
      COALESCE(SUM(m.passed AND m.valid),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.01),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.05),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.10),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.50),0)
      FROM alleles a LEFT JOIN matches m ON a.site_id=m.site_id AND a.base38=m.base38
      GROUP BY a.site_id,a.base38 ORDER BY a.site_id,a.base38"""
    with selection_writer(out / "availability.tsv", SELECTION_AVAIL_FIELDS) as emit:
        for payload, n, passed, valid, f01, f05, f10, f50 in db.execute(query):
            row = json.loads(payload)
            mapped = row["mapping_status"] == "mapped"
            status = ("mapping_unavailable" if not mapped else "not_in_published_summary" if not n
                      else "available_valid_filter_accepted" if valid else "available_filter_accepted_invalid_statistics" if passed
                      else "available_flagged_only")
            row.update(availability_status=status, n_published_contrasts=n, n_accepted_filter_contrasts=passed,
                       n_valid_accepted_filter_contrasts=valid, n_FDR_le_0_01=f01, n_FDR_le_0_05=f05,
                       n_FDR_le_0_10=f10, n_FDR_le_0_50=f50, frequency_zero_inferred=0)
            emit(row)
            qc["availability_" + status] += 1
            key = tuple(row[k] or "unclassified" for k in ("family", "group", "partition", "comparison_partition", "role"))
            c = grouped[key]
            c.update(catalog_bases=1, mapped_bases=int(mapped), matched_bases=int(n > 0), accepted_filter_bases=int(passed > 0),
                     valid_accepted_filter_bases=int(valid > 0), FDR_le_0_01_bases=int(f01 > 0), FDR_le_0_05_bases=int(f05 > 0),
                     FDR_le_0_10_bases=int(f10 > 0), FDR_le_0_50_bases=int(f50 > 0))
    with selection_writer(out / "published_selection_summary.tsv", SELECTION_SUMMARY_FIELDS) as emit:
        for key, counts in sorted(grouped.items()):
            row = dict(zip(SELECTION_SUMMARY_FIELDS[:5], key))
            row.update(counts)
            row.update(matching_fraction=counts["matched_bases"] / counts["catalog_bases"],
                       FDR_le_0_05_fraction_of_valid_accepted_filter_bases=(counts["FDR_le_0_05_bases"] / counts["valid_accepted_filter_bases"]
                                                               if counts["valid_accepted_filter_bases"] else None),
                       unit="modern_base_identity_not_independent_locus")
            emit(row)


@transactional("manifest.json")
def selection_run(args, mapper=None):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    generated = {out / name for name in SELECTION_OWNED}
    inputs = {Path(getattr(args, key)).resolve() for key in
              ("catalog", "selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38")}
    if inputs & {p.resolve() for p in generated}:
        raise ValueError("Input/output collision; use a separate output directory")
    ident = selection_identity(args)
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if old.get("identity") != ident:
            raise ValueError("Input, code or effect contract changed; choose a new output directory")
        if old.get("status") == "complete":
            expected = set(SELECTION_OWNED) - {"manifest.json"}
            if set(old.get("outputs", {})) != expected:
                raise ValueError("Incomplete output receipt; choose a new output directory")
            for name, saved in old["outputs"].items():
                if not (out / name).is_file() or fingerprint(out / name) != saved:
                    raise ValueError(f"Completed output changed: {out / name}")
            return old
    elif any((out / name).exists() for name in SELECTION_OWNED if name != "manifest.json"):
        raise ValueError("Unreceipted output exists; choose a new output directory")
    dump_json(manifest, dict(status="running", identity=ident))
    for name in SELECTION_OWNED:
        if name != "manifest.json":
            (out / name).unlink(missing_ok=True)
            (out / (name + ".partial")).unlink(missing_ok=True)
    qc, db = Counter(), None
    try:
        context = contextlib.nullcontext(mapper) if mapper is not None else GenomeMapper(
            args.fasta38, args.fasta37, args.chain38to37, args.chain37to38)
        with context as active_mapper:
            db = selection_build_index(out / "catalog_index.sqlite", args.catalog, active_mapper, qc)
        with selection_writer(out / "published_selection_alleles.tsv.gz", SELECTION_MATCH_FIELDS) as emit:
            selection_join_statistics(db, args.selection_stats, emit, args.effect_allele, qc)
        selection_export_availability(db, out, qc)
        db.close()
        db = None
        result = dict(status="complete", identity=ident, source_doi=SELECTION_SOURCE, dataset_doi=SELECTION_DATASET, qc=dict(qc),
                      published_fields=list(SELECTION_PUBLISHED), effect_allele_contract=args.effect_allele,
                      accepted_FILTER_values=["PASS", "."],
                      numeric_validation_fields=["AF", "S", "SE", "X", "P_X", "POSTERIOR", "FDR"],
                      orientation_verification="explicit caller declaration, not automatically inferred from field names",
                      statistical_scope="descriptive reuse of published exact biallelic contrasts", limitations=SELECTION_LIMITATIONS,
                      outputs={name: fingerprint(out / name) for name in SELECTION_OWNED if name != "manifest.json"})
        dump_json(manifest, result)
        return result
    except BaseException as exc:
        if db is not None:
            db.close()
        dump_json(manifest, dict(status="failed", identity=ident, qc=dict(qc), error=str(exc)))
        raise


def selection_parser():
    p = argparse.ArgumentParser(description=_description_herv_selection)
    p.add_argument("--catalog", required=True, help="Complete modern site_catalog.tsv.gz; one unique GRCh38 site per row")
    p.add_argument("--selection-stats", required=True, help="Local published GRCh37 TSV/TSV.GZ; exact 16-field header checked")
    for name in ("fasta38", "fasta37", "chain38to37", "chain37to38"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--effect-allele", choices=("raw_only", "ALT"), default="raw_only",
                   help="ALT explicitly declares that raw published S/X refer to its GRCh37 ALT; use only after verifying source documentation. Default preserves raw fields without signed reinterpretation.")
    return p


def _entry_herv_selection():
    selection_result = selection_run(selection_parser().parse_args())
    print(json.dumps({"status": selection_result["status"], "qc": selection_result["qc"]}, sort_keys=True))


# ========================================================================
# herv joint
# ========================================================================
_description_herv_joint = 'Bounded same-person matched HERV statistics on a frozen measurable panel.\n\nPass 1 freezes pair membership using callability alone. Pass 2 reads one pair\nat a time and accumulates sums/counts. No sample-by-catalog matrix is retained.'
"""Bounded same-person matched HERV statistics on a frozen measurable panel.

Pass 1 freezes pair membership using callability alone. Pass 2 reads one pair
at a time and accumulates sums/counts. No sample-by-catalog matrix is retained.
"""
import base64
from analysis_integrity import atomic_json, identity
from allele_measurement import open_pgen_with_pvar, read_counted_dosage
from allele_measurement import stats_number as number, stats_technology as technology, stats_linear_model as linear_model, stats_adjust_fdr as adjust_fdr, stats_fixed_composition as fixed_composition, stats_write_table as write_table, stats_read_table as read_table, stats_seed_for as seed_for


def joint_pair_candidates(targets, selected, fixed_flags):
    if allele_window_matched_catalog(targets):
        return [], [dict(status='not_estimable', reason='window_matches_are_not_ancient_site_pairs')]
    expected=defaultdict(Counter); by_target={t['site_id']:t for t in targets}; found=defaultdict(lambda:{'case':[],'control':[]})
    for target in targets:
        for base in target['alts'].split(','):expected[target['match_id'],allele_alt_role(target,base)][target['group']]+=1
    for entry,keep in zip(selected,fixed_flags):
        if keep and entry['strict_filter_complete']:found[entry['match_id'],entry['role']][entry['group']].append(entry)
    pairs=[];audit=[]
    for key,count in sorted(expected.items()):
        cases,controls=found[key]['case'],found[key]['control']
        good=len(cases)==count['case']==1 and len(controls)==count['control'] and bool(controls)
        good=good and all(e['chrom']==cases[0]['chrom'] for e in controls)
        audit.append(dict(match_id=key[0],role=key[1],expected_cases=count['case'],expected_controls=count['control'],
                          measured_cases=len(cases),measured_controls=len(controls),status='candidate' if good else 'not_estimable',
                          reason='' if good else 'missing_nonunique_or_cross_chromosome_exact_pair'))
        if not good:continue
        case=cases[0];target=by_target[case['site_id']]
        memberships=[v for v in target.get('candidate_families',case['family']).split(',') if v not in {'','.'}]
        scope='boundary_only_sensitivity' if case['partition']=='boundary_only' else 'nonoverlap_main' if len(memberships)<=1 and str(target.get('spatial_overlap','0')) in {'0','.','False','false','none','','single_geometry'} and case['partition']!='ambiguous_overlap' else 'overlap_sensitivity'
        pairs.append(dict(pair_id=key[0]+'|'+key[1],case=case,controls=controls,scope=scope,
                          group=(case['family'],case['partition'],case['role'],scope),membership=memberships))
    return pairs,audit


def joint_iter_pair_values(prefix, n, pairs, *, audit=None):
    """Use the same nucleotide count as individual and family trajectories."""
    with open_pgen_with_pvar(str(prefix)+'.pgen', panel_pvar(prefix), n) as reader:
        def read(entry):
            counted = read_counted_dosage(reader, entry['measurement_plan'], n, audit)
            value = counted.astype(float)/2
            value[counted < 0] = np.nan
            return value
        for pair in pairs:
            if not pair['controls']:
                raise ValueError('A matched pair requires at least one control')
            case = read(pair['case'])
            total = np.zeros(n)
            valid = np.isfinite(case)
            for control in pair['controls']:
                value = read(control)
                valid &= np.isfinite(value)
                total += np.nan_to_num(value)
            controls = total/len(pair['controls'])
            case[~valid] = np.nan
            controls[~valid] = np.nan
            yield pair, case, controls, valid


def joint_sample_cells(samples,args):
    cells=defaultdict(list);formal=[]
    for i,s in enumerate(samples):
        row=dict(s,formal=s.get('include')=='1' and s.get('unrelated')=='1' and s['region']!='unknown',
                 measurement_source=technology(s)+'|'+s.get('data_source','unknown'))
        formal.append(row)
        if not row['formal']:continue
        if s['age_scope']=='primary_ancient':
            cells['all',s['time_bin']].append(i);cells[s['region'],s['time_bin']].append(i)
        elif s['cohort']=='modern':cells[s['region'],'modern'].append(i)
    return formal,{key:np.array(v,dtype=np.int32) for key,v in cells.items()}


def joint_finalize_panel_cells(blocks, panel_counts, cells, args):
    """Each cell must contain every frozen pair; no per-cell panel renormalization."""
    rows=[]
    for group,n_pairs in sorted(panel_counts.items()):
        for cell,indices in sorted(cells.items(),key=lambda x:str(x[0])):
            data=blocks.get((group,cell),{})
            accepted=sum(v[3] for v in data.values())
            ok=accepted==n_pairs and len(indices)>=args.min_samples
            total=np.sum(list(data.values()),axis=0) if data else np.zeros(4)
            row=dict(family=group[0],partition=group[1],role=group[2],membership_scope=group[3],region=cell[0],time_bin=cell[1],
                cohort='modern' if cell[1]=='modern' else 'ancient',n=len(indices),n_fixed_pairs=n_pairs,n_estimable_pairs=int(accepted),
                estimate=float(total[0]/n_pairs) if ok else None,case_mean=float(total[1]/n_pairs) if ok else None,
                control_mean=float(total[2]/n_pairs) if ok else None,status='estimated' if ok else 'not_estimable',
                reason='' if ok else 'some_frozen_pairs_lack_cell_support',measure='fixed_measurable_panel_same_person_pair_difference',
                weights_id='equal_frozen_pair',support_id='|'.join(group),ci_low=None,ci_high=None,
                uncertainty='whole_chromosome_bootstrap_conditional_on_people',n_genomic_blocks=len(data),seed=args.seed)
            if ok and len(data)>=3 and args.bootstrap>=2:
                a=np.array(list(data.values()));rng=np.random.default_rng(seed_for(args.seed,(group,cell)));draws=[]
                for _ in range(args.bootstrap):
                    x=a[rng.integers(len(a),size=len(a))].sum(axis=0);draws.append(x[0]/x[3])
                row.update(ci_low=float(np.quantile(draws,.025)),ci_high=float(np.quantile(draws,.975)))
            rows.append(row)
    return rows


def joint_marker_overlap(args,targets,out):
    basis=getattr(args,'basis_dir',None) or (str(Path(args.coordinates).parent) if args.coordinates else None)
    marker=Path(basis)/'all_pruned.pvar' if basis else None
    if marker is None or not marker.is_file():return {'status':'not_available','needed':'original basis all_pruned.pvar'}
    positions={(str(t.get('chrom37')),int(t['pos37'])) for t in targets if t.get('mapping_status')=='mapped' and t.get('pos37') is not None}
    counts=Counter();overlap=[]
    for row in rows(marker,('CHROM','POS','ID','REF','ALT'),whitespace=True):
        chrom=row['CHROM'].removeprefix('chr');counts[chrom]+=1
        if (chrom,int(row['POS'])) in positions:overlap.append(row)
    write_table(out/'pca_marker_target_overlap.tsv',overlap)
    return dict(status='audited',marker_source=identity(marker,full=True),n_original_markers=sum(counts.values()),
                markers_by_chromosome=dict(counts),n_exact_coordinate_overlap=len(overlap),
                exclusions_applied=False,flanks_added=False,LOCO='not_run_requires_separate_basis_design')


def joint_run_joint(args,targets,selected,samples,coord_contract,coordinates,fixed_flags):
    out=Path(args.out);n=len(samples)
    people,cells=joint_sample_cells(samples,args)
    pairs,pair_audit=joint_pair_candidates(targets,selected,fixed_flags)
    epochs=sorted({s['time_bin'] for s in samples if s['age_scope']=='primary_ancient' and s.get('include')=='1' and s.get('unrelated')=='1'})
    qualified={ep:cells.get(('all',ep),np.array([],dtype=np.int32)) for ep in epochs}
    qualified={ep:idx for ep,idx in qualified.items() if len(idx)>=args.min_samples}
    frozen=[];freeze_audit=[]
    measurement_qc=Counter()
    # Measurability-only selection. No age slope, effect sign, peak or p value is examined.
    for pair,case,control,valid in joint_iter_pair_values(args.pgen_prefix,n,pairs,audit=measurement_qc) if pairs else []:
        rates={ep:float(valid[idx].mean()) for ep,idx in qualified.items()}
        keep=len(rates)>=2 and all(rate>=args.min_call_rate for rate in rates.values())
        freeze_audit.append(dict(pair_id=pair['pair_id'],family=pair['group'][0],partition=pair['group'][1],role=pair['group'][2],
            membership_scope=pair['scope'],candidate_families=','.join(pair['membership']),
            epoch_call_rates=json.dumps(rates,sort_keys=True),frozen=int(keep),
            case_contrast=pair['case']['contrast_id'],control_contrasts=';'.join(c['contrast_id'] for c in pair['controls']),
            status='fixed' if keep else 'not_estimable',reason='' if keep else 'insufficient_same_person_cross_epoch_callability'))
        if keep:frozen.append(pair)
    panel_counts=Counter(p['group'] for p in frozen)
    for row in freeze_audit:
        key=(row['family'],row['partition'],row['role'],row['membership_scope'])
        row['fixed_pair_weight']=1/panel_counts[key] if row['frozen'] else None
    write_table(out/'joint_pair_candidate_audit.tsv',pair_audit)
    write_table(out/'joint_fixed_pair_panel.tsv',freeze_audit)
    atomic_json(out/'joint_panel_plan.json',dict(status='planned',n_fixed_pairs=len(frozen),
        selection='HQ strict counted-nucleotide measurements; same-person pair callability >= threshold in every qualified pooled ancient epoch; frozen before effect estimates',
        qualified_epochs={str(k):len(v) for k,v in qualified.items()},min_call_rate=args.min_call_rate,
        weights='equal fixed pairs; controls equal within each fixed match',
        formal_estimand='pair-by-region-by-time means, then fixed pair weights; any missing pair makes that cell unavailable',
        personal_estimand='observed_pair_conditional unless every frozen pair is measured; formal personal model uses complete panel only'))
    totals={g:np.zeros((n,3)) for g in panel_counts} # case sum, control sum, paired count
    chromosome_totals=defaultdict(dict)
    blocks=defaultdict(lambda:defaultdict(lambda:np.zeros(4)))
    cell_fields=['pair_id','family','partition','role','membership_scope','region','time_bin','cohort','n','n_eligible','call_rate','case_mean','control_mean','difference','status']
    with allele_table(out/'joint_pair_cell_estimates.tsv.gz',cell_fields) as emit, allele_table(out/'joint_pair_call_masks.tsv.gz',['pair_id','PSAM_order_bitmask_base64','n_called','n_PSAM']) as mask_emit:
        for pair,case,control,valid in joint_iter_pair_values(args.pgen_prefix,n,frozen,audit=measurement_qc) if frozen else []:
            g=pair['group'];chrom=pair['case']['chrom'];data=totals[g]
            data[:,0]+=np.nan_to_num(case);data[:,1]+=np.nan_to_num(control);data[:,2]+=valid
            if chrom not in chromosome_totals[g]:chromosome_totals[g][chrom]=np.zeros((n,3))
            chromosome_totals[g][chrom]+=np.column_stack((np.nan_to_num(case),np.nan_to_num(control),valid))
            mask_emit(dict(pair_id=pair['pair_id'],PSAM_order_bitmask_base64=base64.b64encode(np.packbits(valid).tobytes()).decode(),n_called=int(valid.sum()),n_PSAM=n))
            for cell,idx in sorted(cells.items(),key=lambda x:str(x[0])):
                called=idx[valid[idx]];good=len(called)>=args.min_samples and len(called)/len(idx)>=args.min_call_rate
                c=float(case[called].mean()) if good else None;z=float(control[called].mean()) if good else None
                emit(dict(pair_id=pair['pair_id'],family=g[0],partition=g[1],role=g[2],membership_scope=g[3],region=cell[0],time_bin=cell[1],
                    cohort='modern' if cell[1]=='modern' else 'ancient',n=len(called),n_eligible=len(idx),call_rate=len(called)/len(idx),case_mean=c,control_mean=z,
                    difference=c-z if good else None,status='estimated' if good else 'not_estimable'))
                if good:blocks[g,cell][chrom]+=[c-z,c,z,1]
    cell_rows=joint_finalize_panel_cells(blocks,panel_counts,cells,args)
    write_table(out/'joint_fixed_panel_trajectories.tsv',cell_rows)
    # Fixed regional composition, with identical pair sets in every region/epoch.
    region_counts=Counter(s['region'] for s in people if s['formal'] and s['age_scope']=='primary_ancient')
    region_weights={r:n/sum(region_counts.values()) for r,n in sorted(region_counts.items())} if region_counts else {}
    write_table(out/'joint_region_weights.tsv',[dict(region=r,weight=w,weights_id='formal_ancient_region_composition') for r,w in region_weights.items()])
    standardized=[]
    for g in panel_counts:
        for ep in epochs:
            found={r['region']:r for r in cell_rows if (r['family'],r['partition'],r['role'],r['membership_scope'])==g and r['time_bin']==ep}
            missing=[r for r in region_weights if r not in found or found[r]['estimate'] is None]
            good=bool(region_weights) and not missing
            standardized.append(dict(family=g[0],partition=g[1],role=g[2],membership_scope=g[3],time_bin=ep,
                measure='fixed_region_and_pair_composition',n_fixed_pairs=panel_counts[g],
                estimate=sum(w*found[r]['estimate'] for r,w in region_weights.items()) if good else None,
                status='estimated' if good else 'not_estimable',reason='' if good else 'missing_region_pair_support:'+';'.join(missing),
                weights_id='formal_ancient_region_composition',support_id='|'.join(g)))
    write_table(out/'joint_standardized_trajectories.tsv',standardized)
    fits=[];missingness=[]
    verified=coord_contract.get('status')=='verified' and coord_contract.get('metadata_provenance',{}).get('historical_metadata_verification')=='verified'
    explicit_unadjusted=args.pcs==0
    personal_fields=['sample_id','cohort','age_bp','region','family','partition','role','membership_scope','n_fixed_pairs','n_valid_pairs','case_mean','control_mean','paired_difference','callability','estimand','formal_model_eligible','support_id']
    with allele_table(out/'joint_person_statistics.tsv.gz',personal_fields) as emit:
        for g,data in sorted(totals.items()):
            P=panel_counts[g];counts=data[:,2];complete=counts==P
            difference=np.divide(data[:,0]-data[:,1],counts,out=np.full(n,np.nan),where=counts>0)
            strict=np.where(complete,difference,np.nan)
            for i,s in enumerate(samples):
                emit(dict(sample_id=s['sample_id'],cohort=s['cohort'],age_bp=s['age_bp'],region=s['region'],family=g[0],partition=g[1],role=g[2],membership_scope=g[3],
                    n_fixed_pairs=P,n_valid_pairs=int(counts[i]),case_mean=float(data[i,0]/counts[i]) if counts[i] else None,
                    control_mean=float(data[i,1]/counts[i]) if counts[i] else None,paired_difference=float(difference[i]) if counts[i] else None,callability=float(counts[i]/P),
                    estimand='strict_complete_fixed_panel' if complete[i] else 'observed_pair_conditional',
                    formal_model_eligible=int(complete[i] and people[i]['formal'] and s['age_scope']=='primary_ancient'),support_id='|'.join(g)))
            for region in ['all',*sorted(region_counts)]:
                members=[i for i,s in enumerate(people) if s['age_scope']=='primary_ancient' and (region=='all' or s['region']==region)]
                ss=[people[i] for i in members]
                fit=linear_model(ss,strict[members],pcs=coordinates if args.pcs else None,minimum=args.min_samples) if verified or explicit_unadjusted else dict(status='not_estimable',reason='genetic_adjustment_unavailable',n=0,estimate=None,p=None)
                fit.update(family=g[0],partition=g[1],role=g[2],membership_scope=g[3],region=region,n_fixed_pairs=P,
                    measure='strict_complete_fixed_panel_pair_time_association',estimand='complete_panel_people_only; transport_to_incomplete_people_not_assumed',
                    genetic_adjustment='verified_shared_PCs' if verified and args.pcs else 'none_explicit_pcs_0' if explicit_unadjusted else 'unavailable',
                    inference_scope='person_HC3_conditional_on_fixed_loci; chromosome_sensitivity_reported_separately')
                # Whole chromosomes keep all related ALTs and both matched arms together.
                slopes=[]
                if fit.get('estimate') is not None and len(chromosome_totals[g])>=3:
                    for chromosome,removed in chromosome_totals[g].items():
                        n_removed=sum(p['group']==g and p['case']['chrom']==chromosome for p in frozen)
                        denom=P-n_removed
                        if denom<=0:continue
                        response=np.where(complete,(data[:,0]-data[:,1]-removed[:,0]+removed[:,1])/denom,np.nan)
                        loo=linear_model(ss,response[members],pcs=coordinates if args.pcs else None,minimum=args.min_samples)
                        if loo.get('estimate') is not None:slopes.append(loo['estimate'])
                fit.update(n_chromosome_sensitivity=len(slopes),chromosome_jackknife_se=float(np.std(slopes,ddof=0)*math.sqrt(len(slopes)-1)) if len(slopes)>=3 else None,
                           chromosome_sensitivity_status='estimated' if len(slopes)>=3 else 'insufficient_genomic_blocks')
                fits.append(fit)
            members=[i for i,s in enumerate(people) if s['age_scope']=='primary_ancient'];ss=[people[i] for i in members]
            fit=linear_model(ss,counts[members]/P,pcs=coordinates if args.pcs else None,minimum=args.min_samples) if verified or explicit_unadjusted else dict(status='not_estimable',reason='genetic_adjustment_unavailable')
            missingness.append(dict(fit,family=g[0],partition=g[1],role=g[2],membership_scope=g[3],measure='paired_callability_time_association'))
    adjust_fdr(fits);write_table(out/'joint_person_models.tsv',fits);write_table(out/'joint_missingness_models.tsv',missingness)
    funnel=joint_measurability_funnel(out,selected,fixed_flags,freeze_audit,fits,args)
    overlap=joint_marker_overlap(args,targets,out)
    figure_files=[] if args.no_plots else joint_plot_joint(funnel,cell_rows,standardized,out,args)
    result=dict(status='complete' if frozen else 'no_usable_contrast',n_candidate_pairs=len(pairs),n_fixed_pairs=len(frozen),n_PSAM=n,
        per_person_storage='O(samples * family_partition_role * chromosomes); never samples * loci',
        formal_person_model='strict complete fixed panel only; incomplete personal means descriptive observed_pair_conditional',
        formal_cell_model='all frozen pairs with original fixed weights; any unsupported pair produces NA',
        dependencies='person HC3 conditional on fixed loci; genomic chromosome bootstrap and leave-chromosome-out sensitivity separately',
        extra_ALT='unencoded is NA; actual native allele indices are retained; unknown split REF is NA, not an ALT complement',
        measurement_schema='nucleotide-count-1',measurement_qc=dict(measurement_qc),
        measurement_qc_scope='read operations across freeze and estimation passes, not unique biological observations',
        spatial_membership='nonoverlapping main; ambiguous/multiple memberships separately labelled sensitivity; seeded assignment is not unique biology',
        basis_id=coord_contract.get('basis_id'),coordinate_kind=coord_contract.get('coordinate_kind'),
        pca_background_audit=overlap,figures=figure_files,interpretation='fixed measurable allele panel; not ancient original gnomAD O/E or selection coefficients')
    atomic_json(out/'joint_manifest.json',result)
    return result


def joint_measurability_funnel(out,selected,fixed_flags,panel,fits,args):
    counts=defaultdict(lambda:defaultdict(set))
    for row in read_table(out/'all_allele_audit.tsv.gz'):
        if row['role']=='REF':continue
        g=(row['family'],row['partition'],row['role']);allele=row['site_id']+'|'+row['allele']
        counts[g]['catalog_ALT'].add(allele)
        if row['mapping_status']=='mapped':counts[g]['mapped_GRCh37_ALT'].add(allele)
        if row['panel_status'] in {'assayed_pass','not_in_high_quality','pvar_filter_failed'} or row.get('pvar_index') not in {'NA','.','',None}:
            counts[g]['exact_PVAR_REF_ALT'].add(allele)
        if row['panel_status']=='assayed_pass':counts[g]['HQ_exact_ALT'].add(allele)
    for entry,keep in zip(selected,fixed_flags):
        if keep:counts[entry['family'],entry['partition'],entry['role']]['fixed_epoch_callable_ALT'].add(entry['contrast_id'])
    stages_order=['catalog_ALT','mapped_GRCh37_ALT','exact_PVAR_REF_ALT','HQ_exact_ALT','fixed_epoch_callable_ALT','modeled_ALT_contrasts']
    for row in read_table(out/'adjusted_allele_trends.tsv.gz'):
        if row.get('status') in {'descriptive_adjusted_imputed_association','descriptive_no_genetic_adjustment'}:
            counts[row['family'],row['partition'],row['role']]['modeled_ALT_contrasts'].add(row['site_id']+'|'+row['counted_base38'])
    rows=[dict(family=g[0],partition=g[1],role=g[2],stage=stage,n=len(stages[stage]),unit='ALT_contrasts') for g,stages in sorted(counts.items()) for stage in stages_order]
    epochs=defaultdict(set)
    for row in read_table(out/'allele_trajectories.tsv.gz'):
        if row['role']=='REF' or row['region']!='all' or row.get('measurement_source','.') not in {'.','all','NA',''}: continue
        if (number(row.get('n_called')) or 0) >= args.min_samples and (number(row.get('call_rate')) or 0) >= args.min_call_rate:
            epochs[row['family'],row['partition'],row['role'],row['cohort'],row['time_bin']].add(row['contrast_id'])
    rows += [dict(family=g[0],partition=g[1],role=g[2],cohort=g[3],time_bin=g[4],stage='epoch_callable_ALT',n=len(values),unit='ALT_contrasts') for g,values in sorted(epochs.items())]
    for row in panel:
        if row['frozen']:rows.append(dict(family=row['family'],partition=row['partition'],role=row['role'],stage='fixed_same_person_pair',n=1,unit='matched_pair',membership_scope=row['membership_scope']))
    for row in fits:
        rows.append(dict(family=row['family'],partition=row['partition'],role=row['role'],stage='strict_person_model_complete_cases',n=row.get('n',0),unit='people',region=row['region'],membership_scope=row['membership_scope'],status=row['status']))
    write_table(out/'joint_measurability_funnel.tsv',rows)
    return rows


def joint_plot_joint(funnel,cells,standardized,out,args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(16,9))
    stages=['catalog_ALT','mapped_GRCh37_ALT','exact_PVAR_REF_ALT','HQ_exact_ALT','fixed_epoch_callable_ALT']
    for family in sorted({r['family'] for r in funnel}):
        axes[0].plot(range(len(stages)),[sum(r['n'] for r in funnel if r['family']==family and r['stage']==stage) for stage in stages],marker='o',label=family)
    axes[0].set(xticks=range(len(stages)),xticklabels=stages,ylabel='ALT contrasts',title='Measurability funnel');axes[0].tick_params(axis='x',rotation=25);axes[0].legend()
    for family in sorted({r['family'] for r in standardized}):
        for partition in sorted({r['partition'] for r in standardized if r['family']==family}):
            rr=sorted([r for r in standardized if r['family']==family and r['partition']==partition and r['role']=='primary' and r['membership_scope']=='nonoverlap_main'],key=lambda r:r['time_bin'])
            axes[1].plot([(r['time_bin']+.5)*args.bin_years for r in rr],[r['estimate'] if r['estimate'] is not None else np.nan for r in rr],marker='o',label=family+':'+partition)
    axes[1].set(xlabel='Age BP',ylabel='Same-person case minus matched control allele fraction',title='Fixed regional and pair composition; primary ALT');axes[1].legend(fontsize=8)
    fig.suptitle('HERV conditional measurable-panel statistics; not ancient gnomAD enrichment or selection')
    fig.text(.01,.02,'Machine tables: joint_measurability_funnel.tsv / joint_standardized_trajectories.tsv. N, exclusions and provenance: manifest.json.',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    for ext in ('png','pdf'):fig.savefig(out/('F6_joint_measurability.'+ext),dpi=180)
    plt.close(fig)
    return ['F6_joint_measurability.png','F6_joint_measurability.pdf']


def joint_family_composition(family_rows,samples,args):
    """Additional fixed-composition versions of the same full ALT family panel."""
    primary=[s for s in samples if s['age_scope']=='primary_ancient']
    definitions={
        'fixed_region':Counter((s['region'],'all') for s in primary),
        'fixed_region_technology_source':Counter((s['region'],technology(s)+'|'+s.get('data_source','unknown')) for s in primary),
    }
    output=[];weights=[]
    keys=sorted({(r['family'],r['group'],r['partition'],r['role']) for r in family_rows})
    lookup={(r['family'],r['group'],r['partition'],r['role'],r['region'],r.get('measurement_source','all'),r['time_bin']):r for r in family_rows if r['age_scope']=='primary_ancient'}
    epochs=range(math.ceil(args.max_age_bp/args.bin_years))
    for method,counts in definitions.items():
        if not counts:continue
        total=sum(counts.values());ws={key:n/total for key,n in counts.items()}
        weights.extend(dict(method=method,region=key[0],measurement_source=key[1],weight=w,n_reference=counts[key]) for key,w in ws.items())
        for key in keys:
            for ep in epochs:
                strata={cell:lookup.get((*key,*cell,ep)) for cell in ws}
                missing=[cell for cell,r in strata.items() if not r or r['n_metadata']<args.min_samples or r['n_fixed_contrasts']!=r['expected_fixed_contrasts']]
                valid=not missing
                output.append(dict(family=key[0],group=key[1],partition=key[2],role=key[3],time_bin=ep,cohort='ancient',measure=method,
                    estimate=sum(ws[cell]*r['mean_AF'] for cell,r in strata.items()) if valid else None,
                    n=sum(r['n_metadata'] for r in strata.values() if r),n_eff=1/sum(ws[cell]**2/r['n_metadata'] for cell,r in strata.items()) if valid else None,
                    n_fixed_contrasts=next(iter(strata.values()))['expected_fixed_contrasts'] if valid else None,
                    status='descriptive_fixed_composition' if valid else 'not_estimable',
                    reason='' if valid else 'missing_stratum_or_contrast_support:'+json.dumps(missing),
                    uncertainty='not_PC_adjusted; no_person_or_locus_independence_assumed',weights_id=method))
    write_table(Path(args.out)/'family_composition_weights.tsv',weights)
    write_table(Path(args.out)/'standardized_family_trajectories.tsv',output)
    return output


# Mapped Targets
# ======================================================================
_description_mapped_targets = 'Export GRCh37 target labels from an intact completed time-stage mapping.\n\nNo liftover, reference lookup, or genotype decoding is repeated. Every requested\nGRCh38 nucleotide, including the original REF label, is translated explicitly.\nThe GRCh37 REF is the reference already checked by the time-stage mapper, which\ncan differ from the oriented GRCh38 REF. ALT order has no temporal meaning.\n'


from analysis_integrity import validate_manifest

MAPPED_SCHEMA = "allele3-mapped-targets-1"
MAPPED_BASES = frozenset("ACGT")
MAPPED_AUTOSOMES = frozenset(str(c) for c in range(1, 23))
MAPPED_NA = {"", ".", "NA", "N/A", "None", "null"}
MAPPED_COMPLEMENT = str.maketrans("ACGT", "TGCA")
MAPPED_TARGET_FIELDS = ["site_id", "chrom", "pos", "ref", "alts", "requested_bases", "build", "family", "group", "match_id",
                 "chrom38", "pos38", "ref38", "requested_alts38", "strand37to38", "reference_changed",
                 "grch37_ref_among_requested_bases"]
MAPPED_BASE_FIELDS = ["site_id", "chrom38", "pos38", "base38", "role38", "chrom37", "pos37", "base37",
               "role37", "strand37to38", "complemented", "reference_changed"]
MAPPED_EXCLUSION_FIELDS = ["site_id", "chrom38", "pos38", "ref38", "requested_alts38", "reason", "detail"]


def mapped_text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else Path(path).open()


def mapped_chrom(value):
    value = str(value).removeprefix("chr")
    if value.isdigit():
        value = str(int(value))
    if value not in MAPPED_AUTOSOMES:
        raise ValueError("The completed time-stage mapping contract is autosomal")
    return value


def mapped_position(value):
    value = str(value)
    if not value.isdigit() or int(value) < 1:
        raise ValueError("Coordinates must be positive 1-based integers")
    return int(value)


def mapped_allele_set(value, ref, label):
    parts = str(value).upper().split(",")
    if not parts or len(parts) != len(set(parts)) or not set(parts) <= MAPPED_BASES or ref in parts:
        raise ValueError(f"Invalid distinct SNV ALT labels in {label}")
    return frozenset(parts)


def mapped_fingerprint(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def mapped_read_candidates(path, maximum):
    candidates, seen_coordinates, raw_hashes = {}, {}, {}
    duplicates = 0
    with mapped_text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"chrom", "pos", "ref", "alts"} <= set(reader.fieldnames or []):
            raise ValueError("Candidates require chrom,pos,ref,alts and optional site_id")
        for raw in reader:
            c, p, ref = mapped_chrom(raw["chrom"]), mapped_position(raw["pos"]), raw["ref"].upper()
            if ref not in MAPPED_BASES:
                raise ValueError("Candidate REF must be one canonical nucleotide")
            alts = mapped_allele_set(raw["alts"], ref, "candidate")
            site_id = raw.get("site_id", f"{c}:{p}:{ref}")
            if site_id in MAPPED_NA:
                raise ValueError("A supplied candidate site_id cannot be missing")
            for key in ("build", "genome_build", "assembly", "target_build"):
                if raw.get(key) not in MAPPED_NA | {None} and str(raw[key]).lower() not in {"38", "grch38", "hg38"}:
                    raise ValueError("Candidates must be explicitly GRCh38 when a build column is present")
            normalized = dict(raw, site_id=site_id, chrom=c, pos=p, ref=ref, alts=",".join(sorted(alts)))
            digest = mapped_fingerprint(normalized)
            if site_id in candidates:
                if raw_hashes[site_id] != digest:
                    raise ValueError("Conflicting duplicate candidate site_id: " + site_id)
                duplicates += 1
                continue
            if (c, p) in seen_coordinates:
                raise ValueError("Multiple candidate site IDs describe the same original coordinate")
            candidates[site_id], raw_hashes[site_id] = normalized, digest
            seen_coordinates[c, p] = site_id
            if len(candidates) > maximum:
                raise ValueError(f"More than {maximum} candidate sites; use a bounded panel or explicitly raise --max-targets")
    if not candidates:
        raise ValueError("No candidate targets")
    return candidates, duplicates


def mapped_verified_time(time_dir):
    result = validate_manifest(time_dir, require_evidence=True)
    if not str(result.get("version", "")).startswith("herv-all-allele-time-"):
        raise ValueError("A completed herv_allele_time producer manifest is required")
    if "mapped_catalog.tsv.gz" not in result["output_inventory"]:
        raise ValueError("mapped_catalog.tsv.gz is absent from the verified time output inventory")
    return result


def mapped_read_selected_mapping(path, candidates):
    selected, hashes = {}, {}
    counts = Counter()
    with mapped_text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"site_id", "chrom", "pos", "ref", "alts", "mapping_status"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Time mapped_catalog lacks original GRCh38 identity/mapping fields")
        for raw in reader:
            counts["catalog_rows_streamed"] += 1
            sid = raw["site_id"]
            if sid not in candidates:
                continue
            counts["selected_catalog_rows"] += 1
            digest = mapped_fingerprint(raw)
            if sid in selected:
                if hashes[sid] != digest:
                    raise ValueError("Conflicting duplicate source site_id: " + sid)
                counts["identical_source_duplicates"] += 1
                continue
            selected[sid], hashes[sid] = raw, digest
    return selected, dict(counts)


def mapped_oriented(base, strand):
    return base if strand == "+" else base.translate(MAPPED_COMPLEMENT)


def mapped_validated_mapping(row):
    """Match herv_time.map_targets, including genuine build-reference changes."""
    required = {"chrom37", "pos37", "ref37", "alts37", "strand37to38"}
    if not required <= set(row):
        raise ValueError("Mapped source row lacks GRCh37 coordinate/allele/strand fields")
    c38, p38, ref38 = mapped_chrom(row["chrom"]), mapped_position(row["pos"]), row["ref"].upper()
    c37, p37, ref37, strand = mapped_chrom(row["chrom37"]), mapped_position(row["pos37"]), row["ref37"].upper(), row["strand37to38"]
    if ref38 not in MAPPED_BASES or ref37 not in MAPPED_BASES or strand not in {"+", "-"}:
        raise ValueError("Mapped source row has invalid nucleotide or strand")
    source_alts = mapped_allele_set(row["alts"], ref38, "source GRCh38")
    source_bases37 = {mapped_oriented(base, strand) for base in {ref38, *source_alts}}
    labels37 = mapped_allele_set(row["alts37"], ref37, "source GRCh37")
    if set(labels37) != source_bases37 - {ref37}:
        raise ValueError("Source alts37 disagrees with oriented original REF+ALT labels and the actual ref37")
    changed = int(mapped_oriented(ref37, strand) != ref38)
    checks = {"chrom38":c38, "pos38":p38, "ref38":ref38, "source_ref38":ref38,
              "source_to38":"+", "reference_changed":changed}
    for key, expected in checks.items():
        if key in row and row[key] not in MAPPED_NA:
            observed = mapped_chrom(row[key]) if key == "chrom38" else mapped_position(row[key]) if key == "pos38" else int(row[key]) if key == "reference_changed" else row[key]
            if observed != expected:
                raise ValueError("Source mapping metadata disagrees with its original identity: " + key)
    if row.get("context38") not in MAPPED_NA | {None}:
        context = row["context38"].upper()
        if len(context) != 3 or context[1] != ref38:
            raise ValueError("Source context38 disagrees with ref38")
    return dict(chrom38=c38, pos38=p38, ref38=ref38, chrom37=c37, pos37=p37, ref37=ref37,
                strand37to38=strand, reference_changed=changed, source_bases38={ref38, *source_alts})


def mapped_exclusion(target, reason, detail=""):
    return dict(site_id=target["site_id"], chrom38=target["chrom"], pos38=target["pos"], ref38=target["ref"],
                requested_alts38=target["alts"], reason=reason, detail=detail)


def mapped_translate(candidates, selected):
    proposed, excluded = {}, []
    for sid, target in candidates.items():
        row = selected.get(sid)
        if row is None:
            excluded.append(mapped_exclusion(target, "not_in_time_catalog"))
            continue
        if (mapped_chrom(row["chrom"]), mapped_position(row["pos"])) != (target["chrom"], target["pos"]):
            excluded.append(mapped_exclusion(target, "source38_coordinate_mismatch"))
            continue
        if row["ref"].upper() != target["ref"]:
            excluded.append(mapped_exclusion(target, "source38_reference_mismatch", "Source REF="+row["ref"]))
            continue
        source_alts = mapped_allele_set(row["alts"], row["ref"].upper(), "source GRCh38")
        requested = {target["ref"], *target["alts"].split(",")}
        missing = requested - {row["ref"].upper(), *source_alts}
        if missing:
            excluded.append(mapped_exclusion(target, "requested_bases_missing_from_source", ",".join(sorted(missing))))
            continue
        if row["mapping_status"] != "mapped":
            excluded.append(mapped_exclusion(target, "source_mapping_unavailable", row["mapping_status"]))
            continue
        mapping = mapped_validated_mapping(row)
        proposed[sid] = dict(mapping, requested_bases38=requested)
    by_coordinate = defaultdict(list)
    for sid, mapping in proposed.items():
        by_coordinate[mapping["chrom37"], mapping["pos37"]].append(sid)
    for coordinate, ids in by_coordinate.items():
        if len(ids) < 2:
            continue
        detail = f"GRCh37 {coordinate[0]}:{coordinate[1]}; candidates=" + ",".join(sorted(ids))
        for sid in ids:
            excluded.append(mapped_exclusion(candidates[sid], "mapped_coordinate_collision", detail))
            del proposed[sid]
    targets37, base_rows = [], []
    for sid, mapping in sorted(proposed.items(), key=lambda x: (int(x[1]["chrom37"]), x[1]["pos37"], x[0])):
        target, strand = candidates[sid], mapping["strand37to38"]
        bases37 = {mapped_oriented(base, strand) for base in mapping["requested_bases38"]}
        targets37.append(dict(site_id=sid, chrom=mapping["chrom37"], pos=mapping["pos37"], ref=mapping["ref37"],
            alts=",".join(sorted(bases37-{mapping["ref37"]})), requested_bases=",".join(sorted(bases37)), build="GRCh37",
            family=target.get("family"), group=target.get("group"), match_id=target.get("match_id"),
            chrom38=target["chrom"], pos38=target["pos"], ref38=target["ref"], requested_alts38=target["alts"],
            strand37to38=strand, reference_changed=mapping["reference_changed"],
            grch37_ref_among_requested_bases=int(mapping["ref37"] in bases37)))
        for base in [target["ref"], *sorted(mapping["requested_bases38"]-{target["ref"]})]:
            base37 = mapped_oriented(base, strand)
            base_rows.append(dict(site_id=sid, chrom38=target["chrom"], pos38=target["pos"], base38=base,
                role38="REF_label" if base == target["ref"] else "ALT", chrom37=mapping["chrom37"], pos37=mapping["pos37"],
                base37=base37, role37="REF_label" if base37 == mapping["ref37"] else "ALT",
                strand37to38=strand, complemented=int(strand == "-"), reference_changed=mapping["reference_changed"]))
    return targets37, base_rows, sorted(excluded, key=lambda x: (int(x["chrom38"]), x["pos38"], x["site_id"]))


def mapped_write_table(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:"NA" if row.get(k) is None else row[k] for k in fields})


@transactional("manifest.json", input_dirs=("time_dir",), input_files=("targets",), input_prefixes=())
def mapped_run(args):
    if args.max_targets < 1:
        raise ValueError("--max-targets must be a positive explicit limit")
    time_dir = Path(args.time_dir)
    source_manifest = time_dir/"manifest.json"
    mapped = time_dir/"mapped_catalog.tsv.gz"
    print(f"[dataset] mapped-targets: time-dir={time_dir}; candidates={args.targets}", flush=True)
    receipt = mapped_verified_time(time_dir)
    print("[build-check] GRCh38 -> GRCh37: reusing the completed time-stage mapping and checked reference labels; no new liftOver, FASTA or PGEN read", flush=True)
    source_before = identity(source_manifest, full=True)
    mapped_before = identity(mapped)
    targets_before = identity(args.targets, full=True)
    candidates, candidate_duplicates = mapped_read_candidates(args.targets, args.max_targets)
    selected, scan_counts = mapped_read_selected_mapping(mapped, candidates)
    targets37, base_rows, exclusions = mapped_translate(candidates, selected)
    if len(targets37) + len(exclusions) != len(candidates):
        raise ValueError("Internal candidate accounting failure")
    if (source_before != identity(source_manifest, full=True) or mapped_before != identity(mapped)
            or targets_before != identity(args.targets, full=True)):
        raise ValueError("Mapping receipt or candidate inputs changed during export")
    out = Path(args.out)
    mapped_write_table(out/"targets_grch37.tsv", MAPPED_TARGET_FIELDS, targets37)
    mapped_write_table(out/"base_mapping.tsv", MAPPED_BASE_FIELDS, base_rows)
    mapped_write_table(out/"exclusions.tsv", MAPPED_EXCLUSION_FIELDS, exclusions)
    status = "not_available" if not targets37 else "complete_with_exclusions" if exclusions else "complete"
    print(f"[mapped-targets] requested={len(candidates)} exported={len(targets37)} excluded={len(exclusions)}; {status}", flush=True)
    print(f"[build-check] mapped-targets: exported={len(targets37)}; excluded={len(exclusions)}; unreported/missing mappings remain unavailable", flush=True)
    return dict(status=status, schema=MAPPED_SCHEMA, n_candidates=len(candidates), n_exported=len(targets37), n_excluded=len(exclusions),
                n_base_mappings=len(base_rows), identical_candidate_duplicates=candidate_duplicates, stream_counts=scan_counts,
                exclusion_reasons=dict(Counter(r["reason"] for r in exclusions)), source_time_status=receipt["status"],
                source_time_version=receipt["version"], source_time_manifest=source_before, source_mapped_catalog=mapped_before,
                candidate_identity=targets_before, source_build="GRCh38", target_build="GRCh37",
                source_validation="Completed time output_inventory checked before reading; mapped_catalog must be inventoried; small manifest/candidate SHA and mapped-catalog size/mtime rechecked after streaming",
                selected_record_validation="Original coordinates and REF must match; every requested nucleotide must be in the source; selected duplicate IDs must have identical records",
                mapping_semantics="Reuse time mapper's reciprocal single-base mapping and FASTA37 REF; orient every original REF/ALT nucleotide; no new liftOver or reference query",
                requested_bases_contract="The requested_bases column is the exact oriented original REF+requested ALT set; a structural GRCh37 REF outside it must not become an additional measurement",
                limits=["A mapping/reference label does not establish that an allele is observed or encoded in the UKB panel.",
                        "Original GRCh38 site_id is preserved; read chrom/pos/ref fields for GRCh37 identity.",
                        "ALT output is a sorted set, never mutation order, age, priority, or allele-frequency rank.",
                        "The GRCh37 REF can be outside the requested modern bases; it is an assembly label, not a newly observed allele.",
                        "Unmapped and colliding candidates are excluded, never converted into zero counts.",
                        "Only selected candidate records receive semantic checks; unrelated catalog rows are streamed and discarded."])


def mapped_parser():
    p = argparse.ArgumentParser(description=_description_mapped_targets)
    p.add_argument("--time-dir", required=True, help="Completed time stage containing its verified manifest.json and mapped_catalog.tsv.gz")
    p.add_argument("--targets", required=True, help="Existing GRCh38 candidates/targets.tsv")
    p.add_argument("--max-targets", type=int, default=5000)
    p.add_argument("--out", required=True)
    return p


def mapped_main(argv=None):
    result = mapped_run(mapped_parser().parse_args(argv))
    return 0 if result["n_exported"] else 2


if __name__ == '__main__':
    if sys.argv[1:2] == ['mapped-targets']:
        raise SystemExit(mapped_main(sys.argv[2:]))
    routes = {'time': _entry_herv_time, 'allele': _entry_herv_allele_time, 'selection': _entry_herv_selection}
    arguments = list(sys.argv[1:])
    action = arguments.pop(0) if arguments and arguments[0] in routes else 'time'
    if action is None:
        raise SystemExit('Choose a stage: ' + ', '.join(routes))
    sys.argv[1:] = arguments
    routes[action]()


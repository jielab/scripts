#!/usr/bin/env python3
"""Full-catalog, full-PSAM ancient allele trajectories with bounded genotype RAM.

Every modern nucleotide is audited independently against its actual panel allele
index. An unrepresented base is unavailable, never AC=0. Native multiallelic
PGEN preserves direct REF and ALT counts; biallelic/split REF is unavailable
without a joint-genotype provenance contract. A single ALT complement is never
reported as REF. The existing analysis remains an ALT-trajectory analysis.
"""
from __future__ import annotations
from analysis_integrity import transactional

import argparse
from collections import Counter, OrderedDict, defaultdict
import contextlib
import csv
import gzip
import json
import math
import os
from pathlib import Path
import sys
import numpy as np
from scipy.stats import t as student_t
import herv_time as ht
from allele_common import GenomeMapper, canonical_chrom, orient, dump_json
from allele_measurement import catalog_measurements, panel_pvar, iter_pgen_chunks

VERSION = "herv-all-allele-time-2.0-nucleotide-counts"
LIMITATIONS = [
    "The catalog is a complete supplied modern AC5 multi-ALT HERV/body-or-edge case universe and matched controls, not all human genomic SNVs or all ancient variation.",
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
OWNED = ("mapped_catalog.tsv.gz", "all_allele_audit.tsv.gz", "allele_trajectories.tsv.gz", "measurement_plans.json", "allele_callability.tsv.gz",
         "adjusted_allele_trends.tsv.gz", "sample_coverage.tsv.gz", "sample_time_region_coverage.tsv",
         "panel_availability_counts.tsv", "catalog_coverage.tsv", "fixed_panel_family_trajectories.tsv",
         "matched_locus_role_trajectories.tsv", "matched_locus_audit.tsv", "model_covariates.json",
         "family_allele_time.png", "family_allele_time.pdf", "allele_time_heatmap.png", "allele_time_heatmap.pdf",
         "allele_panel_availability.png", "allele_panel_availability.pdf", "trends.raw.tsv.gz")
LOCUS_FIELDS = ["locus_n_alt_ac5", "locus_n_alt_ac10", "locus_ac10_supported"]
AUDIT_FIELDS = ["site_id", "chrom", "pos", "ref", "allele", "role", "modern_ac", "modern_an", "modern_af",
                "modern_ac5", "modern_ac10", *LOCUS_FIELDS, "family", "group", "partition", "comparison_partition", "match_id",
                "strict_filter_complete", "mapping_status", "panel_status", "contrast_id", "pvar_index",
                "panel_ref37", "panel_alt37", "counted_base37", "invert_dosage", "reference_count_scope",
                "panel_alleles37", "allele_id", "allele_idx", "measurement_mode", "counting_scope", "comparison_kind", "source_id"]
TRAJECTORY_FIELDS = ["site_id", "contrast_id", "family", "group", "partition", "match_id", "role", "counted_base38",
                     "paired_modern_alt", "counting_scope", "modern_ac", "modern_ac5", "modern_ac10", *LOCUS_FIELDS, "strict_filter_complete",
                     "cohort", "age_scope", "time_bin", "region", "measurement_source", "start_bp", "end_bp", "n_eligible", "n_called",
                     "ac", "an", "af", "hardcall_heterozygosity", "call_rate", "fixed_panel", "status",
                     "counted_base37", "allele_id", "pvar_index", "allele_idx", "measurement_mode", "reference_count_scope"]
TREND_FIELDS = ["allele_index", "site_id", "contrast_id", "family", "group", "partition", "role", "counted_base38",
                "modern_ac", "modern_ac5", "modern_ac10", *LOCUS_FIELDS, "strict_filter_complete", "fixed_panel", "status", "n", "n_bins",
                "model_call_rate", "reason", "expected_pc_count", "covariates", "condition_number", "time_vif", "slope_per_1000yr_toward_present", "se",
                "ci_low", "ci_high", "p_two_sided", "q_bh_all_alt", "q_by_all_alt", "q_bh_primary_ac5", "q_by_primary_ac5",
                "allele_id", "allele_idx", "measurement_mode", "counting_scope"]


@contextlib.contextmanager
def table(path, fields):
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


def load_catalog(path):
    result, positions, ids = [], set(), set()
    for raw in ht.rows(path, ("chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id")):
        r = dict(raw)
        r["chrom"], r["pos"], r["ref"] = canonical_chrom(r["chrom"]), int(r["pos"]), r["ref"].upper()
        all_alts = r.get("observed_alts_all", r["alts"])
        if all_alts in ht.NA: all_alts = r["alts"]
        # New full catalogs explicitly align all_alt_ac/an to their alts field.
        if "all_alt_ac" in r: all_alts = r["alts"]
        bases = all_alts.upper().split(",")
        if r["chrom"] not in {str(x) for x in range(1, 23)} or r["pos"] < 1 or r["ref"] not in ht.BASES:
            raise ValueError("Catalog must contain canonical autosomal GRCh38 SNVs")
        if not bases or len(set(bases)) != len(bases) or any(b not in ht.BASES or b == r["ref"] for b in bases):
            raise ValueError(f"Invalid catalog ALTs: {r['site_id']}")
        key = r["chrom"], r["pos"]
        if key in positions or r["site_id"] in ids: raise ValueError("Catalog requires unique positions and site_id values")
        positions.add(key); ids.add(r["site_id"])
        if r["group"] not in {"case", "control"}: raise ValueError("Catalog group must be case or control")
        if any(r[name] in ht.NA for name in ("site_id","family","match_id")):
            raise ValueError("Catalog site_id, family and match_id must be nonmissing")
        counts = {}
        ac_text, an_text = r.get("all_alt_ac", r.get("alt_ac_all", "")), r.get("all_alt_an", r.get("alt_an_all", ""))
        if ac_text not in ht.NA and an_text not in ht.NA:
            ac, an = [int(x) for x in ac_text.split(",")], [int(x) for x in an_text.split(",")]
            if len(ac) != len(bases) or len(an) != len(bases) or any(a < 0 or n < a or n < 1 for a,n in zip(ac, an)):
                raise ValueError(f"Catalog all-ALT AC/AN mismatch: {r['site_id']}")
            counts = {b: (a,n) for b,a,n in zip(bases, ac, an)}
        locus_counts = {}
        for threshold in (5, 10):
            name = f"n_alt_ac{threshold}"
            recorded = int(r[name]) if r.get(name, ".") not in ht.NA else None
            derived = sum(a >= threshold for a,_ in counts.values()) if counts else recorded
            if recorded is not None and (recorded < 0 or recorded > len(bases) or recorded != derived):
                raise ValueError(f"Catalog {name} disagrees with all-ALT counts: {r['site_id']}")
            locus_counts[f"locus_{name}"] = derived
        n10 = locus_counts["locus_n_alt_ac10"]
        supported = int(n10 >= 2) if n10 is not None else None
        if r.get("ac10_supported", ".") not in ht.NA:
            recorded = int(r["ac10_supported"])
            if recorded not in {0, 1} or (supported is not None and recorded != supported):
                raise ValueError(f"Catalog ac10_supported disagrees with locus ALT counts: {r['site_id']}")
            supported = recorded
        locus_counts["locus_ac10_supported"] = supported
        for role in ("primary_alt", "secondary_alt"):
            if r.get(role, ".") not in ht.NA and r[role] not in bases:
                raise ValueError(f"Catalog {role} is not an observed ALT: {r['site_id']}")
        if r.get("primary_alt", ".") not in ht.NA and r.get("primary_alt") == r.get("secondary_alt"):
            raise ValueError("Primary and secondary ALT roles must be distinct")
        r.update(alts=",".join(bases), _counts=counts, **locus_counts,
                 partition=r.get("partition", r.get("case_region", "unclassified")),
                 strict_filter_complete=int(r.get("strict_filter_complete", "0") == "1"))
        result.append(r)
    return result


def alt_role(target, base):
    if base == target.get("primary_alt"): return "primary"
    if base == target.get("secondary_alt"): return "secondary"
    ac = target["_counts"].get(base, (None,None))[0]
    return "additional_lowcount" if ac is not None and ac < 5 else "additional"


def catalog_contrasts(mapped, availability):
    """Return ALT measurements and their audit, preserving the two-result API."""
    selected, audit, _ = catalog_measurements(mapped, availability, orient=orient, alt_role=alt_role)
    # Keep the existing modern REF/ALT catalog identifier. This label identifies
    # the catalog nucleotide; counting_scope describes the actual measurement.
    for entry in audit:
        if entry["role"] != "REF":
            entry["contrast_id"] = f"{entry['site_id']}|{entry['ref']}>{entry['allele']}"
    return selected, audit


def load_samples(path, psam, bin_years, max_age):
    ids = [r["IID"] for r in ht.rows(psam, ("IID",), whitespace=True)]
    if len(set(ids)) != len(ids): raise ValueError("PSAM IID values must be unique")
    meta = {}
    for r in ht.rows(path, ("sample_id", "age_bp", "cohort", "include")):
        if r["sample_id"] in meta: raise ValueError("Duplicate metadata sample_id")
        meta[r["sample_id"]] = r
    individual_counts = Counter(r.get('individual_id') for r in meta.values() if r.get('individual_id') not in ht.NA | {None})
    result = []
    for sid in ids:
        raw = dict(meta.get(sid, {}))
        age = ht.finite(raw.get("age_bp")); cohort = raw.get("cohort", "unknown")
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
        if raw["region"] in ht.NA: raw["region"] = "unknown"
        result.append(raw)
    return result, len(set(meta) - set(ids))


def sample_groups(samples, bin_years, max_age, minimum):
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
    from temporal_stats import technology
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


class TrendModel:
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
        self.use_source = any(s.get("data_source", ".") not in ht.NA for s in candidates)
        self.use_coverage = any(ht.finite(s.get("target_coverage")) is not None for s in candidates)
        for i,s in enumerate(samples):
            eligible[i] &= (not pc_names or s["sample_id"] in coordinates)
            eligible[i] &= (not self.use_region or s["region"] != "unknown")
            eligible[i] &= (not self.use_source or s.get("data_source", ".") not in ht.NA)
            coverage = ht.finite(s.get("target_coverage"))
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


def fdr(pvalues, mask=None):
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


def aggregate_key(entry):
    return entry["family"],entry["group"],entry["partition"],entry["role"]


def block_interval(blocks, bootstrap, seed):
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
def run(args):
    if args.bin_years < 1 or args.max_age_bp < args.bin_years or args.min_samples < 2 or not 0 < args.min_call_rate <= 1:
        raise ValueError("Invalid chronology or callability settings")
    if not 1 <= args.chunk_size <= 1024 or args.pcs < 0 or args.block_bp < 100000:
        raise ValueError("Require chunk-size 1..1024, nonnegative PCs, and blocks >=100 kb")
    out = Path(args.out).resolve(); out.mkdir(parents=True,exist_ok=True)
    inputs = [args.catalog,args.metadata,args.hq_variants,args.coordinates,args.fasta38,args.fasta37,args.chain38to37,args.chain37to38]
    inputs += [str(args.pgen_prefix)+s for s in (".pgen",".pvar",".pvar.zst",".psam")]
    owned = {out/p for p in (*OWNED,"manifest.json")}
    if any(Path(p).resolve() in owned for p in inputs if p): raise ValueError("Inputs and generated time outputs must not overlap")
    dump_json(out/"manifest.json",dict(status="running",version=VERSION))
    for name in OWNED: (out/name).unlink(missing_ok=True)
    prefix = str(args.pgen_prefix)
    samples, extra_meta = load_samples(args.metadata,prefix+".psam",args.bin_years,args.max_age_bp)
    groups = sample_groups(samples,args.bin_years,args.max_age_bp,args.min_samples)
    from pca_contract import read_scores, historical_identity, coordinate_join
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
    model = TrendModel(samples, groups, coords, pc_names, args.min_call_rate, args.pcs, unavailable)
    dump_json(out/"model_covariates.json",dict(model="HC3_OLS_of_imputed_ALT_fraction_ancient_only", direction="two_sided; positive increases toward present",
        coordinate_contract=coord_audit, genetic_adjustment="unavailable" if unavailable else "none_explicit_pcs_0" if args.pcs == 0 else "verified_shared_PCs", candidate_covariates=model.names, n_complete_case_samples=int(model.eligible.sum()), coordinate_kind=coord_kind,
        coverage_adjusted=model.use_coverage, region_adjusted=model.use_region, source_adjusted=model.use_source,
        missing_requested_covariates="complete_case_exclusion", cache_max_design_masks=4))
    targets = load_catalog(args.catalog)
    with GenomeMapper(args.fasta38,args.fasta37,args.chain38to37,args.chain37to38) as mapper:
        targets = ht.map_targets(targets,mapper)
    map_fields = sorted(set().union(*(r.keys() for r in targets)) - {"_counts"}) if targets else ["site_id","mapping_status"]
    ht.write_rows(out/"mapped_catalog.tsv.gz",targets,map_fields)
    pvar = panel_pvar(prefix)
    availability = ht.index_panel(targets,pvar,args.hq_variants)
    selected,audit = catalog_contrasts(targets,availability)
    dump_json(out/"measurement_plans.json", dict(schema="nucleotide-count-1", ALT_measurements=selected,
        unknown_split_reference="unavailable; no complement is relabelled REF"))
    del availability
    ht.write_rows(out/"all_allele_audit.tsv.gz",audit,AUDIT_FIELDS)
    available_counts = Counter((r["family"],r["group"],r["partition"],r["role"],r["strict_filter_complete"],r["panel_status"]) for r in audit)
    ht.write_rows(out/"panel_availability_counts.tsv",(dict(zip(["family","group","partition","role","strict_filter_complete","panel_status"],k),n_alleles=v) for k,v in sorted(available_counts.items())),
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
    ht.write_rows(out/"catalog_coverage.tsv",(dict(zip(["level","family","group","partition","metric"],k),value=v) for k,v in sorted(coverage.items())),
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
            expected_pairs[(target["match_id"],alt_role(target,base))][target["group"]] += 1
    call_fields = ["site_id","contrast_id","role","fixed_panel","min_qualified_epoch_call_rate","n_called_all_PSAM","n_PSAM","strict_filter_complete"]
    with table(out/"allele_trajectories.tsv.gz",TRAJECTORY_FIELDS) as trajectory, \
         table(out/"allele_callability.tsv.gz",call_fields) as call_audit, \
         table(out/"trends.raw.tsv.gz",TREND_FIELDS) as trend:
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
                basic.update({k:entry[k] for k in LOCUS_FIELDS})
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
                        keys = [aggregate_key(entry),(entry["family"],entry["group"],"all_partitions",entry["role"])]
                        for key in keys:
                            family_blocks[(*key,group_index)][block] += [af,1,het,2*n]
                if keep[local] and strict[local]:
                    for key in (aggregate_key(entry), (entry['family'],entry['group'],'all_partitions',entry['role'])):
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
    qbh,qby = fdr(pvalues); pbh,pby = fdr(pvalues,primary_tests)
    def adjusted_rows():
        for row in ht.rows(out/"trends.raw.tsv.gz"):
            i = int(row["allele_index"])
            for name,values_ in (("q_bh_all_alt",qbh),("q_by_all_alt",qby),("q_bh_primary_ac5",pbh),("q_by_primary_ac5",pby)):
                row[name] = float(values_[i]) if np.isfinite(values_[i]) else None
            yield row
    ht.write_rows(out/"adjusted_allele_trends.tsv.gz",adjusted_rows(),TREND_FIELDS)
    (out/"trends.raw.tsv.gz").unlink()
    family_rows = []
    for (family,group,partition,role,gi),blocks in sorted(family_blocks.items()):
        mean,lo,hi,n_blocks = block_interval({b:v[:2] for b,v in blocks.items()},args.bootstrap,args.seed)
        summed = np.sum(list(blocks.values()),axis=0)
        g = groups[gi]
        family_rows.append(dict(family=family,group=group,partition=partition,role=role,region=g["region"],measurement_source=g.get("measurement_source","all"),expected_fixed_contrasts=family_expected[family,group,partition,role],cohort=g["cohort"],age_scope=g["age_scope"],time_bin=g["time_bin"],
            start_bp=g["start_bp"],end_bp=g["end_bp"],n_metadata=g["n_eligible"],n_fixed_contrasts=int(summed[1]),mean_AF=mean,ci_low=lo,ci_high=hi,
            mean_heterozygosity=float(summed[2]/summed[1]),summed_effective_AN=int(summed[3]),n_blocks=n_blocks,
            status="fixed_panel_descriptive" if n_blocks >= 20 else "few_blocks_descriptive"))
    ht.write_rows(out/"fixed_panel_family_trajectories.tsv",[r for r in family_rows if r['region']=='all' and r['measurement_source']=='all'],list(family_rows[0]) if family_rows else ["status"])
    ht.write_rows(out/"regional_fixed_panel_family_trajectories.tsv",family_rows,list(family_rows[0]) if family_rows else ["status"])
    from herv_joint import family_composition
    family_composition(family_rows, samples, args)
    pair_blocks = defaultdict(lambda:defaultdict(lambda:np.zeros(2)))
    match_audit = []
    for (match_id,role),expected in expected_pairs.items():
        pair = pair_records.get((match_id,role),{"case":[],"control":[]})
        cases,controls = pair["case"],pair["control"]
        valid = len(cases) == expected["case"] == 1 and len(controls) == expected["control"] and bool(controls)
        same_chromosome = bool(cases) and all(c["chrom"] == cases[0][0]["chrom"] for c,_ in controls)
        valid &= same_chromosome
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
        point,lo,hi,n_blocks = block_interval(blocks,args.bootstrap,args.seed)
        g = groups[gi]
        matched_rows.append(dict(family=family,partition=partition,role=role,measure=measure,cohort=g["cohort"],age_scope=g["age_scope"],time_bin=g["time_bin"],
            start_bp=g["start_bp"],end_bp=g["end_bp"],n_matched_loci=int(sum(v[1] for v in blocks.values())),n_blocks=n_blocks,difference=point,ci_low=lo,ci_high=hi,
            denominator="locus/role matched group estimates; not same-person paired complete calls",
            bootstrap_unit="whole_chromosome_case_and_controls_together",status="descriptive_chromosome_bootstrap"))
    ht.write_rows(out/"matched_locus_role_trajectories.tsv",matched_rows,list(matched_rows[0]) if matched_rows else ["status"])
    ht.write_rows(out/"matched_locus_audit.tsv",match_audit,["match_id","role","n_case_contrasts","n_control_contrasts","expected_case_contrasts","expected_control_contrasts","status"])
    sample_fields = ["sample_id","individual_id","cohort","age_bp","age_scope","unrelated","include","region","time_bin","data_source","target_coverage",
                     "n_HQ_exact_contrasts","n_called_HQ_contrasts","call_rate_HQ","n_fixed_contrasts","n_called_fixed_contrasts","model_covariate_complete"]
    ht.write_rows(out/"sample_coverage.tsv.gz",(dict(s,n_HQ_exact_contrasts=n_variants,n_called_HQ_contrasts=int(sample_called[i]),
        call_rate_HQ=float(sample_called[i]/n_variants) if n_variants else None,n_fixed_contrasts=fixed_count,
        n_called_fixed_contrasts=int(sample_fixed[i]),model_covariate_complete=int(model.eligible[i])) for i,s in enumerate(samples)),sample_fields)
    ht.write_rows(out/"sample_time_region_coverage.tsv",(dict({k:v for k,v in g.items() if k != "members"},
        n_model_covariate_complete=int(model.eligible[g["members"]].sum()),summed_called_HQ_contrasts=int(sample_called[g["members"]].sum())) for g in groups),
        ["cohort","age_scope","time_bin","region","start_bp","end_bp","n_eligible","qualified_epoch","n_model_covariate_complete","summed_called_HQ_contrasts"])
    from herv_joint import run_joint
    joint = run_joint(args, targets, selected, samples, coord_audit, coords, fixed_flags)
    if not args.no_plots: figures([r for r in family_rows if r['region']=='all' and r['measurement_source']=='all'],available_counts,out)
    catalog_receipt = Path(args.catalog).parent/"candidate_manifest.json"
    catalog_provenance = json.loads(catalog_receipt.read_text()) if catalog_receipt.is_file() else {}
    audit_status_counts = Counter()
    for key,value in available_counts.items(): audit_status_counts[key[5]] += value
    result = dict(joint_analysis=joint, coordinate_contract=coord_audit, genetic_adjustment="unavailable" if unavailable else "none_explicit_pcs_0" if args.pcs == 0 else "verified_shared_PCs", version=VERSION,status=("complete_adjustment_unavailable" if unavailable else "complete_descriptive_all_allele_panel") if selected else "no_usable_exact_high_quality_contrast",
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
        catalog_source_scope=catalog_provenance.get("scope",catalog_provenance.get("case_selection","supplied catalog; scope must be checked")),
        outputs=[name for name in OWNED if (out/name).is_file()],
        limitations=LIMITATIONS,parameters={k:v for k,v in vars(args).items()})
    dump_json(out/"manifest.json",result)
    return result


def figures(rows,availability,out):
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


def parser():
    p = argparse.ArgumentParser(description=__doc__)
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


if __name__ == "__main__":
    try: run(parser().parse_args())
    except (ValueError,OSError,RuntimeError) as exc: raise SystemExit(f"herv-allele-time: {exc}")

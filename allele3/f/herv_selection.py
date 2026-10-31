#!/usr/bin/env python3
"""Attach published selection statistics to exact modern-catalog allele labels.

This optional, offline command does not fit a selection model. It independently
maps the complete supplied GRCh38 catalog, streams the published GRCh37 table,
and uses a disk index to bound memory. Raw statistics retain their original
REF/ALT contrast. Signed allele effects are unavailable unless the caller
explicitly declares that the source table's S and X refer to ALT.
"""
from __future__ import annotations
from analysis_integrity import transactional

import argparse
from collections import Counter, defaultdict
import contextlib
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3

import herv_time as ht
from allele_common import GenomeMapper, canonical_chrom, dump_json, fingerprint, orient

VERSION = "herv-published-selection-1.0"
AUTOSOMES = frozenset(str(i) for i in range(1, 23))
SOURCE = "https://doi.org/10.1038/s41586-026-10358-1"
DATASET = "https://doi.org/10.7910/DVN/7RVV9N"
PUBLISHED = ("CHROM", "POS", "REF", "ALT", "ANC", "ID", "RSID", "AF", "S", "SE",
             "X", "P_X", "POSTERIOR", "FDR", "CHI2_BE", "FILTER")
BASE_FIELDS = ["site_id", "chrom", "pos", "ref", "counted_base38", "role", "modern_ac", "modern_an",
               "modern_af", "family", "group", "partition", "comparison_partition", "match_id",
               "strict_filter_complete", "mapping_status", "chrom37", "pos37", "ref37", "counted_base37",
               "strand37to38", "reference_changed"]
MATCH_FIELDS = BASE_FIELDS + ["published_contrast", "counted_base_relation", "statistics_status",
                             "published_effect_contract", "oriented_s", "oriented_x", "source_doi"] + list(PUBLISHED)
AVAIL_FIELDS = BASE_FIELDS + ["availability_status", "n_published_contrasts", "n_accepted_filter_contrasts",
                             "n_valid_accepted_filter_contrasts", "n_FDR_le_0_01", "n_FDR_le_0_05", "n_FDR_le_0_10",
                             "n_FDR_le_0_50", "frequency_zero_inferred"]
SUMMARY_FIELDS = ["family", "group", "partition", "comparison_partition", "role", "catalog_bases",
                  "mapped_bases", "matched_bases", "accepted_filter_bases", "valid_accepted_filter_bases", "FDR_le_0_01_bases",
                  "FDR_le_0_05_bases", "FDR_le_0_10_bases", "FDR_le_0_50_bases", "matching_fraction",
                  "FDR_le_0_05_fraction_of_valid_accepted_filter_bases", "unit"]
OWNED = ("catalog_index.sqlite", "published_selection_alleles.tsv.gz", "availability.tsv",
         "published_selection_summary.tsv", "manifest.json")
LIMITATIONS = [
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


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def code_identity():
    return {name: sha256(Path(__file__).with_name(name))
            for name in ("herv_selection.py", "herv_time.py", "0.common.py")}


def identity(args):
    paths = {key: getattr(args, key) for key in
             ("catalog", "selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38")}
    inputs = {key: fingerprint(path) for key, path in paths.items()}
    inputs["catalog"]["sha256"] = sha256(args.catalog)
    for key in ("fasta38", "fasta37"):
        index = Path(str(paths[key]) + ".fai")
        if index.is_file():
            inputs[key + "_fai"] = fingerprint(index)
        compressed_index = Path(str(paths[key]) + ".gzi")
        if str(paths[key]).endswith((".gz", ".bgz")) and compressed_index.is_file():
            inputs[key + "_gzi"] = fingerprint(compressed_index)
    return dict(version=VERSION, code=code_identity(), inputs=inputs, effect_allele=args.effect_allele,
                mapping="independent reciprocal GRCh38-to-GRCh37 with FASTA checks")


@contextlib.contextmanager
def writer(path, fields):
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


def catalog_rows(path):
    """Read one catalog site at a time; SQLite verifies global uniqueness."""
    required = ("chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id")
    for raw in ht.rows(path, required):
        row = dict(raw)
        row["chrom"], row["pos"], row["ref"] = canonical_chrom(row["chrom"]), int(row["pos"]), row["ref"].upper()
        alts_text = row.get("observed_alts_all", row["alts"])
        if alts_text in ht.NA or "all_alt_ac" in row:
            alts_text = row["alts"]
        alts = alts_text.upper().split(",")
        if (row["chrom"] not in AUTOSOMES or row["pos"] < 1 or row["ref"] not in ht.BASES
                or not alts or len(set(alts)) != len(alts) or any(b not in ht.BASES or b == row["ref"] for b in alts)):
            raise ValueError(f"Catalog requires canonical autosomal SNVs: {row['site_id']}")
        if row["site_id"] in ht.NA or row["match_id"] in ht.NA or row["group"] not in {"case", "control"}:
            raise ValueError("Catalog needs nonmissing site_id/match_id and case/control group")
        counts = {}
        ac_text, an_text = row.get("all_alt_ac", row.get("alt_ac_all", ".")), row.get("all_alt_an", row.get("alt_an_all", "."))
        if (ac_text in ht.NA) != (an_text in ht.NA):
            raise ValueError("Catalog AC and AN must either both be present or both be unavailable")
        if ac_text not in ht.NA:
            ac, an = [int(x) for x in ac_text.split(",")], [int(x) for x in an_text.split(",")]
            if len(ac) != len(alts) or len(an) != len(alts) or any(a < 0 or n < 1 or a > n for a, n in zip(ac, an)):
                raise ValueError(f"Catalog ALT AC/AN order or bounds invalid: {row['site_id']}")
            counts = {base: (a, n) for base, a, n in zip(alts, ac, an)}
        for key in ("primary_alt", "secondary_alt"):
            if row.get(key, ".") not in ht.NA and row[key] not in alts:
                raise ValueError(f"{key} is not a catalog ALT: {row['site_id']}")
        row.update(alts=",".join(alts), _counts=counts,
                   partition=row.get("partition", row.get("case_region", "unclassified")),
                   comparison_partition=row.get("comparison_partition", row.get("partition", "unclassified")))
        yield row


def base_info(target, base):
    ac, an = target["_counts"].get(base, (None, None))
    role = ("REF" if base == target["ref"] else "primary" if base == target.get("primary_alt")
            else "secondary" if base == target.get("secondary_alt")
            else "additional_lowcount" if ac is not None and ac < 5 else "additional")
    row = {key: target.get(key) for key in BASE_FIELDS}
    row.update(counted_base38=base, role=role, modern_ac=ac, modern_an=an, modern_af=ac / an if an else None,
               strict_filter_complete=int(str(target.get("strict_filter_complete", "0")) == "1"))
    row["counted_base37"] = orient(base, target["strand37to38"]) if target["mapping_status"] == "mapped" else None
    return row


def build_index(path, catalog, mapper, qc):
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
        for target in catalog_rows(catalog):
            target = ht.map_targets([target], mapper)[0]
            mapped = target["mapping_status"] == "mapped"
            qc["catalog_sites"] += 1
            qc["mapping_" + target["mapping_status"]] += 1
            db.execute("INSERT INTO sites VALUES(?,?,?,?,?,?)", (target["site_id"], target["chrom"], target["pos"],
                       target.get("chrom37") if mapped else None, target.get("pos37") if mapped else None, json.dumps(target)))
            for base in [target["ref"], *target["alts"].split(",")]:
                info = base_info(target, base)
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


def statistics_valid(row):
    # CHI2_BE is retained verbatim; missing batch diagnostics do not erase the
    # published inferential fields, and no new batch threshold is invented.
    values = {key: ht.finite(row[key]) for key in ("AF", "S", "SE", "X", "P_X", "POSTERIOR", "FDR")}
    return (all(v is not None for v in values.values()) and values["SE"] > 0
            and all(0 <= values[k] <= 1 for k in ("AF", "P_X", "POSTERIOR", "FDR")))


def join_statistics(db, path, emit, effect_allele, qc):
    for row in ht.rows(path, PUBLISHED):
        qc["published_rows_streamed"] += 1
        if qc["published_rows_streamed"] % 1000000 == 0:
            print(f"[published-selection] {qc['published_rows_streamed']:,} source rows streamed", flush=True)
        chrom = canonical_chrom(row["CHROM"])
        try:
            pos = int(row["POS"])
        except ValueError as exc:
            raise ValueError("Published POS must be an integer") from exc
        ref, alt = row["REF"], row["ALT"]
        if chrom not in AUTOSOMES or pos < 1 or ref not in ht.BASES or alt not in ht.BASES or ref == alt:
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
        valid, passed = statistics_valid(row), row["FILTER"] in {"PASS", "."}
        qc["matched_published_contrasts"] += 1
        qc["matched_valid_accepted_filter_contrasts"] += int(valid and passed)
        for base37, relation in ((ref, "REF"), (alt, "ALT")):
            base38 = orient(base37, target["strand37to38"])
            info = base_info(target, base38)
            if base38 not in {target["ref"], *target["alts"].split(",")}:
                raise ValueError("Mapped contrast contains a base absent from the modern catalog")
            info.update(row)
            signed = effect_allele == "ALT" and valid
            sign = 1 if relation == "ALT" else -1
            info.update(published_contrast=key, counted_base_relation=relation,
                        statistics_status="valid" if valid else "invalid_numeric_fields",
                        published_effect_contract="caller_declared_ALT" if effect_allele == "ALT" else "raw_only_unverified",
                        oriented_s=sign * float(row["S"]) if signed else None,
                        oriented_x=sign * float(row["X"]) if signed else None, source_doi=SOURCE)
            emit(info)
            db.execute("INSERT INTO matches VALUES(?,?,?,?,?,?)", (target["site_id"], base38, key, int(passed), int(valid),
                       float(row["FDR"]) if valid else None))
            qc["matched_base_contrast_rows"] += 1
        if qc["matched_published_contrasts"] % 10000 == 0:
            db.commit()
    db.commit()


def export_availability(db, out, qc):
    grouped = defaultdict(Counter)
    query = """SELECT a.payload,COUNT(m.contrast),COALESCE(SUM(m.passed),0),
      COALESCE(SUM(m.passed AND m.valid),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.01),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.05),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.10),0),
      COALESCE(SUM(m.passed AND m.valid AND m.fdr<=0.50),0)
      FROM alleles a LEFT JOIN matches m ON a.site_id=m.site_id AND a.base38=m.base38
      GROUP BY a.site_id,a.base38 ORDER BY a.site_id,a.base38"""
    with writer(out / "availability.tsv", AVAIL_FIELDS) as emit:
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
    with writer(out / "published_selection_summary.tsv", SUMMARY_FIELDS) as emit:
        for key, counts in sorted(grouped.items()):
            row = dict(zip(SUMMARY_FIELDS[:5], key))
            row.update(counts)
            row.update(matching_fraction=counts["matched_bases"] / counts["catalog_bases"],
                       FDR_le_0_05_fraction_of_valid_accepted_filter_bases=(counts["FDR_le_0_05_bases"] / counts["valid_accepted_filter_bases"]
                                                               if counts["valid_accepted_filter_bases"] else None),
                       unit="modern_base_identity_not_independent_locus")
            emit(row)


@transactional("manifest.json")
def run(args, mapper=None):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    generated = {out / name for name in OWNED}
    inputs = {Path(getattr(args, key)).resolve() for key in
              ("catalog", "selection_stats", "fasta38", "fasta37", "chain38to37", "chain37to38")}
    if inputs & {p.resolve() for p in generated}:
        raise ValueError("Input/output collision; use a separate output directory")
    ident = identity(args)
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if old.get("identity") != ident:
            raise ValueError("Input, code or effect contract changed; choose a new output directory")
        if old.get("status") == "complete":
            expected = set(OWNED) - {"manifest.json"}
            if set(old.get("outputs", {})) != expected:
                raise ValueError("Incomplete output receipt; choose a new output directory")
            for name, saved in old["outputs"].items():
                if not (out / name).is_file() or fingerprint(out / name) != saved:
                    raise ValueError(f"Completed output changed: {out / name}")
            return old
    elif any((out / name).exists() for name in OWNED if name != "manifest.json"):
        raise ValueError("Unreceipted output exists; choose a new output directory")
    dump_json(manifest, dict(status="running", identity=ident))
    for name in OWNED:
        if name != "manifest.json":
            (out / name).unlink(missing_ok=True)
            (out / (name + ".partial")).unlink(missing_ok=True)
    qc, db = Counter(), None
    try:
        context = contextlib.nullcontext(mapper) if mapper is not None else GenomeMapper(
            args.fasta38, args.fasta37, args.chain38to37, args.chain37to38)
        with context as active_mapper:
            db = build_index(out / "catalog_index.sqlite", args.catalog, active_mapper, qc)
        with writer(out / "published_selection_alleles.tsv.gz", MATCH_FIELDS) as emit:
            join_statistics(db, args.selection_stats, emit, args.effect_allele, qc)
        export_availability(db, out, qc)
        db.close()
        db = None
        result = dict(status="complete", identity=ident, source_doi=SOURCE, dataset_doi=DATASET, qc=dict(qc),
                      published_fields=list(PUBLISHED), effect_allele_contract=args.effect_allele,
                      accepted_FILTER_values=["PASS", "."],
                      numeric_validation_fields=["AF", "S", "SE", "X", "P_X", "POSTERIOR", "FDR"],
                      orientation_verification="explicit caller declaration, not automatically inferred from field names",
                      statistical_scope="descriptive reuse of published exact biallelic contrasts", limitations=LIMITATIONS,
                      outputs={name: fingerprint(out / name) for name in OWNED if name != "manifest.json"})
        dump_json(manifest, result)
        return result
    except BaseException as exc:
        if db is not None:
            db.close()
        dump_json(manifest, dict(status="failed", identity=ident, qc=dict(qc), error=str(exc)))
        raise


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", required=True, help="Complete modern site_catalog.tsv.gz; one unique GRCh38 site per row")
    p.add_argument("--selection-stats", required=True, help="Local published GRCh37 TSV/TSV.GZ; exact 16-field header checked")
    for name in ("fasta38", "fasta37", "chain38to37", "chain37to38"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--effect-allele", choices=("raw_only", "ALT"), default="raw_only",
                   help="ALT explicitly declares that raw published S/X refer to its GRCh37 ALT; use only after verifying source documentation. Default preserves raw fields without signed reinterpretation.")
    return p


if __name__ == "__main__":
    result = run(parser().parse_args())
    print(json.dumps({"status": result["status"], "qc": result["qc"]}, sort_keys=True))

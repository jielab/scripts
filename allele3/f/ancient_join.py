"""Streaming selection joins and bounded candidate design for allele3 evolution.

The sampling universe is modern multi-ALT sites. It is NOT an all-site
case/control enrichment analysis. SQLite bounds RAM; no genotype matrix or
AlphaGenome request occurs in this module.
"""
from __future__ import annotations

from array import array
from bisect import bisect_left
from collections import Counter, OrderedDict, defaultdict
import csv
import gzip
import heapq
import itertools
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
import time

from core import dump_json, fingerprint, fmt, open_text, stable_hash, tsv_writer
from gnomad_io import canonical_chrom
from harmonize import AUTOSOMES, BASES, GenomeMapper, MappingError, orient

SCHEMA = "allele3-evolution-1"
SELECTION_FIELDS = ("CHROM", "POS", "REF", "ALT", "ANC", "ID", "RSID", "AF", "S", "SE",
                    "X", "P_X", "POSTERIOR", "FDR", "CHI2_BE", "FILTER")
JOIN_FIELDS = ["site_id", "chrom", "pos", "ref", "chrom37", "pos37", "ref37", "alt37",
               "strand37to38", "reference_changed", "source_ref38", "anchor_ref", "anchor_alt",
               "extra_alt", "extra_ac", "extra_an", "extra_af", "extra_nfe_ac", "extra_nfe_an",
               "extra_nfe_af", "anchor_alt_modern_status", "anchor_alt_modern_af", "selection_id",
               "rsid", "ancestral_base37", "selection_pooled_af", "s", "se", "x", "p_x", "posterior",
               "fdr", "chi2_be", "selection_filter", "context", "canonical_context", "region_flags",
               "has_non_snv", "an_relative_spread", "n_modern_alt", "consequence_terms", "gene_symbols",
               "estimated_increasing_base", "estimated_decreasing_base", "extra_selection_status", "population_class"]
CANDIDATE_FIELDS = JOIN_FIELDS + ["group", "match_set", "match_stratum", "locus_block",
                                "selection_evidence", "comparison_universe"]


def table_rows(path, required=()):
    """Read public TSVs with ## metadata and a #CHROM or ordinary header."""
    with open_text(path) as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("##"))
        reader = csv.DictReader(lines, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty table: {path}")
        reader.fieldnames = [name.lstrip("#") for name in reader.fieldnames]
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {', '.join(sorted(missing))}")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed/truncated TSV row in {path}: line {reader.line_num}")
            yield row


def number(value, default=None):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return default
    return x if math.isfinite(x) else default


def _modern_sites(shards):
    for shard in shards:
        site_file = Path(shard) / "multiallelic_sites.tsv.gz"
        allele_file = Path(shard) / "multiallelic_alleles.tsv.gz"
        allele_groups = itertools.groupby(table_rows(allele_file, ("chrom", "pos", "ref", "alt", "status")),
                                          key=lambda r: (r["chrom"], r["pos"], r["ref"]))
        for site in table_rows(site_file, ("chrom", "pos", "ref", "n_alt")):
            try:
                key, group = next(allele_groups)
            except StopIteration as exc:
                raise ValueError("Modern site/allele files have different site counts") from exc
            if key != (site["chrom"], site["pos"], site["ref"]):
                raise ValueError("Modern site/allele file order or keys disagree")
            alleles = [r for r in group if r["status"] == "observed_pass"]
            if len(alleles) != int(site["n_alt"]):
                raise ValueError("Modern occupancy differs between site and allele tables")
            yield site, alleles
        if next(allele_groups, None) is not None:
            raise ValueError("Unmatched trailing modern allele records")


def _stage_identity(kind, inputs, parameters):
    return dict(schema=SCHEMA, join_schema=2 if kind == "join" else None, kind=kind, inputs={k: fingerprint(v) for k, v in inputs.items()},
                parameters=parameters)


def completed_stage(out, identity):
    path = Path(out) / "manifest.json"
    if not path.exists():
        return None
    result = json.loads(path.read_text())
    if result.get("identity") != identity:
        raise ValueError(f"Inputs/settings changed at {out}; choose a new output directory")
    for name, fp in result.get("outputs", {}).items():
        p = Path(out) / name
        if not p.exists() or fingerprint(p) != fp:
            raise ValueError(f"Completed output changed: {p}; restore it or choose a new output directory")
    return result


def finish_stage(out, identity, result, filenames):
    out = Path(out)
    result.update(identity=identity, outputs={name: fingerprint(out / name) for name in filenames},
                  manifest_path=str((out / "manifest.json").resolve()))
    dump_json(out / "manifest.json", result)
    return result


def _contrast_rows(row, mapped, site, alleles, qc):
    """Shared exact allele/numeric rules for both join strategies. REF is checked by the caller."""
    ref, alt = row["REF"], row["ALT"]
    chrom = canonical_chrom(row["CHROM"])
    strand = mapped["strand37to38"]
    anchor_ref, anchor_alt = orient(ref, strand), orient(alt, strand)
    by_base = {a["alt"]: a for a in alleles}
    anchor_modern = by_base.get(anchor_alt)
    s, posterior, pooled_af = number(row["S"]), number(row["POSTERIOR"]), number(row["AF"])
    x, se = number(row["X"]), number(row["SE"])
    if (s is None or posterior is None or pooled_af is None or x is None or se is None or se <= 0
            or not 0 <= posterior <= 1 or not 0 <= pooled_af <= 1):
        qc["invalid_selection_numeric"] += 1
        return
    for allele in alleles:
        extra = allele["alt"]
        if extra in (anchor_ref, anchor_alt):
            continue
        qc["extra_allele_contrasts"] += 1
        result = dict(site_id=mapped["site_id"], chrom=mapped["chrom38"], pos=mapped["pos38"],
                      ref=mapped["ref38"], chrom37=chrom, pos37=row["POS"], ref37=ref, alt37=alt,
                      strand37to38=strand, reference_changed=mapped["reference_changed"],
                      source_ref38=mapped["source_ref38"], anchor_ref=anchor_ref, anchor_alt=anchor_alt,
                      extra_alt=extra, extra_ac=allele["ac"], extra_an=allele["an"], extra_af=allele["af"],
                      extra_nfe_ac=allele.get("ac_nfe", "."), extra_nfe_an=allele.get("an_nfe", "."),
                      extra_nfe_af=allele.get("af_nfe", "."),
                      anchor_alt_modern_status=("observed_pass" if anchor_modern else
                                               "source_reference_frequency_unknown" if anchor_alt == mapped["source_ref38"]
                                               else "not_observed_pass"),
                      anchor_alt_modern_af=anchor_modern["af"] if anchor_modern else ".",
                      selection_id=row["ID"], rsid=row["RSID"], ancestral_base37=row["ANC"],
                      selection_pooled_af=row["AF"], s=row["S"], se=row["SE"], x=row["X"],
                      p_x=row["P_X"], posterior=row["POSTERIOR"], fdr=row["FDR"], chi2_be=row["CHI2_BE"],
                      selection_filter=row["FILTER"], context=mapped["context38"],
                      canonical_context=mapped["canonical_context"], region_flags=site["region_flags"],
                      has_non_snv=site["has_non_snv"], an_relative_spread=site["an_relative_spread"],
                      n_modern_alt=site["n_alt"], consequence_terms=allele["consequence_terms"],
                      gene_symbols=allele["gene_symbols"],
                      estimated_increasing_base=anchor_alt if s > 0 else anchor_ref if s < 0 else ".",
                      estimated_decreasing_base=anchor_ref if s > 0 else anchor_alt if s < 0 else ".",
                      extra_selection_status="not_estimated_by_anchor_contrast",
                      population_class=site.get("population_class", "NA"))
        yield result


def _join_modern_first(modern_out, selection_path, out_dir, fasta38, fasta37,
                       chain38to37, chain37to38, *, mapper=None):
    """Join all QC published selection contrasts to modern multi-ALT loci.

    Different ancient ALT contrasts at one coordinate remain different rows.
    A contrast's s is NEVER assigned to extra_alt. Missing summary entries are
    not declared absent from PVAR. Only point SNVs and autosomes are accepted.
    """
    from cli import load_run
    modern = load_run(modern_out)
    source_build = modern["options"]["build"]
    inputs = dict(modern_manifest=Path(modern_out) / "run.json", selection=selection_path,
                  fasta38=fasta38, fasta37=fasta37, chain38to37=chain38to37, chain37to38=chain37to38)
    inputs.update({"fai38": str(fasta38) + ".fai", "fai37": str(fasta37) + ".fai"})
    identity = _stage_identity("join", inputs, dict(modern_signature=modern["signature"], source_build=source_build,
                                                                 strategy="modern-first"))
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if old := completed_stage(out, identity):
        return old
    start = time.monotonic()
    db_path = out / "modern_index.sqlite"
    tmp = out / "modern_index.sqlite.partial"
    if tmp.exists():
        tmp.unlink()
    db = sqlite3.connect(tmp)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA cache_size=-65536")
    db.execute("CREATE TABLE modern (chrom37 TEXT, pos37 INTEGER, payload TEXT, matched INTEGER DEFAULT 0, "
               "PRIMARY KEY(chrom37,pos37)) WITHOUT ROWID")
    qc = Counter()
    own_mapper = mapper is None
    mapper = mapper or GenomeMapper(fasta38, fasta37, chain38to37, chain37to38)
    rejects_handle, rejects = tsv_writer(out / "mapping_rejections.tsv.gz", ["chrom", "pos", "ref", "reason"])
    try:
        with rejects_handle:
            for site, alleles in _modern_sites(modern["shards"]):
                qc["modern_multi_alt_sites"] += 1
                try:
                    mapped = mapper.map(site["chrom"], int(site["pos"]), site["ref"], source_build)
                except MappingError as exc:
                    qc["mapping_" + exc.reason] += 1
                    rejects.writerow({**{k: site[k] for k in ("chrom", "pos", "ref")}, "reason": exc.reason})
                    continue
                oriented = []
                for allele in alleles:
                    a = dict(allele)
                    a["alt"] = orient(a["alt"], mapped["source_to38"])
                    oriented.append(a)
                payload = dict(mapped=mapped, site=site, alleles=oriented)
                try:
                    db.execute("INSERT INTO modern(chrom37,pos37,payload) VALUES(?,?,?)",
                               (mapped["chrom37"], mapped["pos37"], json.dumps(payload, separators=(",", ":"))))
                except sqlite3.IntegrityError as exc:
                    raise ValueError("Two modern sites map to the same ancient coordinate; ambiguous source merge") from exc
                qc["mapped_sites"] += 1
                qc["reference_changed_sites"] += mapped["reference_changed"]
                if qc["modern_multi_alt_sites"] % 100000 == 0:
                    db.commit()
                    print(f"[map] {qc['modern_multi_alt_sites']:,} multi-ALT sites", flush=True)
        db.commit()
        # Reuse the exact, bounded coordinate cache to avoid millions of
        # payload SELECT misses on mounted filesystems. Mapping still covers
        # every modern site; no selection statistic or allele is prefiltered.
        db.execute("CREATE TEMP VIEW coordinates AS SELECT chrom37 AS source_chrom, "
                   "pos37 AS source_pos, NULL AS reason FROM modern")
        cache = _CoordinateCache(db)
        # Keep misses inside a transaction too. Before the first matching UPDATE,
        # implicit sqlite transactions otherwise acquire a file lock for every
        # SELECT; millions of misses are especially slow on WSL mounted drives.
        db.execute("BEGIN")
        handle, writer = tsv_writer(out / "selection_intersection.tsv.gz", JOIN_FIELDS)
        with handle:
            for row in table_rows(selection_path, SELECTION_FIELDS):
                qc["selection_rows"] += 1
                if qc["selection_rows"] % 1000000 == 0:
                    db.commit()
                    db.execute("BEGIN")
                    print(f"[selection] {qc['selection_rows']:,} published rows streamed", flush=True)
                ref, alt = row["REF"], row["ALT"]
                if ref not in BASES or alt not in BASES or ref == alt:
                    qc["selection_non_snv"] += 1
                    continue
                chrom = canonical_chrom(row["CHROM"])
                if not cache.contains(chrom, int(row["POS"])):
                    continue
                record = db.execute("SELECT payload FROM modern WHERE chrom37=? AND pos37=?",
                                    (chrom, int(row["POS"]))).fetchone()
                if record is None:
                    continue
                payload = json.loads(record[0])
                mapped, site, alleles = payload["mapped"], payload["site"], payload["alleles"]
                if ref != mapped["ref37"]:
                    qc["selection_ref37_mismatch"] += 1
                    continue
                db.execute("UPDATE modern SET matched=1 WHERE chrom37=? AND pos37=?", (chrom, int(row["POS"])))
                qc["matched_contrasts"] += 1
                writer.writerows(_contrast_rows(row, mapped, site, alleles, qc))
        db.commit()
        qc["mapped_sites_with_summary"] = db.execute("SELECT COUNT(*) FROM modern WHERE matched=1").fetchone()[0]
        qc["mapped_sites_without_summary"] = qc["mapped_sites"] - qc["mapped_sites_with_summary"]
        db.close()
        os.replace(tmp, db_path)
    except BaseException:
        db.close()
        raise
    finally:
        if own_mapper:
            mapper.close()
    return finish_stage(out, identity, dict(qc=dict(qc), elapsed_seconds=time.monotonic() - start,
                         strategy="modern-first", qc_scope="All modern multi-ALT sites receive full mapping; reference implementation",
                         performance=dict(fullmap_calls=qc["modern_multi_alt_sites"], index_bytes=db_path.stat().st_size,
                                          coordinate_cache=cache.metrics()),
                         intersection_path=str(out / "selection_intersection.tsv.gz"),
                         index_path=str(db_path), universe="modern multi-ALT sites; no all-site enrichment denominator; "
                         "extra must be a directly observed passing source ALT, never inferred from the source REF"),
                        ["selection_intersection.tsv.gz", "mapping_rejections.tsv.gz", "modern_index.sqlite"])


def _process_memory():
    """Observed process memory, not a promise about chain-library allocations."""
    result = {"rss_bytes": None, "lifetime_peak_rss_bytes": None}
    try:
        with open("/proc/self/statm", encoding="ascii") as handle:
            result["rss_bytes"] = int(handle.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError, IndexError):
        pass
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["lifetime_peak_rss_bytes"] = int(peak * (1 if sys.platform == "darwin" else 1024))
    except (ImportError, OSError, ValueError):
        pass
    return result


class _CoordinateCache:
    """Exact per-chromosome uint arrays, with a strict array-byte LRU budget.

    No Python set/list of millions of positions is constructed. An individual
    chromosome larger than the budget uses indexed SQLite membership instead,
    preserving every target. Array capacity is allocated exactly before filling.
    The tiny chromosome/count dictionary and LRU bookkeeping are reported apart
    from array bytes. Source shards need not be globally chromosome sorted.
    """
    def __init__(self, db, max_bytes=64 * 1024 * 1024, max_chromosomes=2):
        if max_bytes < 1 or max_chromosomes < 1:
            raise ValueError("Coordinate cache limits must be positive")
        self.db, self.max_bytes, self.max_chromosomes = db, max_bytes, max_chromosomes
        self.arrays, self.array_bytes, self.stats = OrderedDict(), 0, Counter()
        self.info = {}
        for chrom, count, largest in db.execute(
                "SELECT source_chrom,COUNT(*),MAX(source_pos) FROM coordinates "
                "WHERE reason IS NULL GROUP BY source_chrom"):
            code = "I" if largest <= (1 << (8 * array("I").itemsize)) - 1 else "Q"
            self.info[chrom] = (count, code)
        self.stats["all_coordinates_uncompressed_uint_bytes"] = sum(
            count * array(code).itemsize for count, code in self.info.values())

    def _update_peak(self):
        self.stats["peak_array_bytes"] = max(self.stats["peak_array_bytes"], self.array_bytes)
        overhead = sys.getsizeof(self.arrays) + sys.getsizeof(self.info)
        overhead += sum(sys.getsizeof(k) + sys.getsizeof(v) for k, v in self.info.items())
        overhead += sum(sys.getsizeof(k) for k in self.arrays)
        self.stats["peak_bookkeeping_bytes"] = max(self.stats["peak_bookkeeping_bytes"], overhead)

    def contains(self, chrom, pos):
        self.stats["membership_lookups"] += 1
        info = self.info.get(chrom)
        if info is None:
            return False
        values = self.arrays.get(chrom)
        if values is not None:
            self.arrays.move_to_end(chrom)
            self.stats["cached_array_lookups"] += 1
        else:
            count, code = info
            size = sys.getsizeof(array(code)) + count * array(code).itemsize
            if size > self.max_bytes:
                self.stats["disk_membership_queries"] += 1
                return self.db.execute("SELECT 1 FROM coordinates WHERE source_chrom=? AND source_pos=? "
                                       "AND reason IS NULL", (chrom, pos)).fetchone() is not None
            while self.arrays and (len(self.arrays) >= self.max_chromosomes or
                                   self.array_bytes + size > self.max_bytes):
                _, removed = self.arrays.popitem(last=False)
                self.array_bytes -= sys.getsizeof(removed)
                del removed
                self.stats["array_evictions"] += 1
            values = array(code, [0]) * count
            for index, (position,) in enumerate(self.db.execute(
                    "SELECT source_pos FROM coordinates WHERE source_chrom=? AND reason IS NULL ORDER BY source_pos",
                    (chrom,))):
                values[index] = position
            self.arrays[chrom] = values
            self.array_bytes += sys.getsizeof(values)
            if self.array_bytes > self.max_bytes:
                raise RuntimeError("Coordinate array allocator exceeded its declared budget")
            self.stats["array_loads"] += 1
            self._update_peak()
        index = bisect_left(values, pos)
        return index < len(values) and values[index] == pos

    def metrics(self):
        self._update_peak()
        return dict(self.stats, array_budget_bytes=self.max_bytes,
                    max_cached_chromosomes=self.max_chromosomes,
                    retained_array_bytes=self.array_bytes,
                    indexed_coordinates=sum(v[0] for v in self.info.values()),
                    indexed_chromosomes=len(self.info))


def _join_selection_first(modern_out, selection_path, out_dir, fasta38, fasta37,
                          chain38to37, chain37to38, *, mapper=None):
    from cli import load_run
    modern = load_run(modern_out)
    source_build = modern["options"]["build"]
    if source_build not in ("GRCh37", "GRCh38"):
        raise ValueError("source_build must be GRCh37 or GRCh38")
    inputs = dict(modern_manifest=Path(modern_out) / "run.json", selection=selection_path,
                  fasta38=fasta38, fasta37=fasta37, chain38to37=chain38to37, chain37to38=chain37to38)
    inputs.update({"fai38": str(fasta38) + ".fai", "fai37": str(fasta37) + ".fai"})
    identity = _stage_identity("join", inputs, dict(modern_signature=modern["signature"], source_build=source_build,
                                                   strategy="selection-first"))
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if old := completed_stage(out, identity):
        return old
    start = time.monotonic()
    before_memory = _process_memory()
    db_path, tmp = out / "selection_index.sqlite", out / "selection_index.sqlite.partial"
    if tmp.exists():
        tmp.unlink()
    db = sqlite3.connect(tmp)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA cache_size=-65536")
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA mmap_size=0")
    db.execute("CREATE TABLE coordinates (chrom37 TEXT, pos37 INTEGER, source_chrom TEXT, source_pos INTEGER, "
               "strand TEXT, reason TEXT, modern_seen INTEGER DEFAULT 0, "
               "PRIMARY KEY(chrom37,pos37)) WITHOUT ROWID")
    db.execute("CREATE TABLE contrasts (chrom37 TEXT, pos37 INTEGER, row_no INTEGER, payload TEXT, "
               "PRIMARY KEY(chrom37,pos37,row_no)) WITHOUT ROWID")
    qc = Counter()
    own_mapper = mapper is None
    payload_bytes = 0
    try:
        mapper = mapper or GenomeMapper(fasta38, fasta37, chain38to37, chain37to38)
        coord_handle, coord_rejects = tsv_writer(out / "selection_coordinate_rejections.tsv.gz",
                                                ["chrom37", "pos37", "reason"])
        previous_key, previous_state = None, None
        with coord_handle:
            for row_no, row in enumerate(table_rows(selection_path, SELECTION_FIELDS), 1):
                qc["selection_rows"] += 1
                if row_no % 1000000 == 0:
                    db.commit()
                    print(f"[selection index] {row_no:,} published rows streamed", flush=True)
                ref, alt = row["REF"], row["ALT"]
                if ref not in BASES or alt not in BASES or ref == alt:
                    qc["selection_non_snv"] += 1
                    continue
                qc["selection_snv_rows"] += 1
                chrom37, pos37 = canonical_chrom(row["CHROM"]), int(row["POS"])
                key = (chrom37, pos37)
                state = (previous_state if key == previous_key else db.execute(
                    "SELECT source_chrom,source_pos,strand,reason FROM coordinates WHERE chrom37=? AND pos37=?",
                    key).fetchone())
                if state is None:
                    qc["selection_unique_coordinates"] += 1
                    try:
                        if chrom37 not in AUTOSOMES:
                            raise MappingError("non_autosomal_source")
                        if pos37 < 1:
                            raise MappingError("invalid_source_position")
                        if source_build == "GRCh38":
                            qc["selection_coordinate_map_calls"] += 1
                            coords = mapper.map_coordinates(chrom37, pos37, "GRCh37")
                            state = (coords["chrom38"], coords["pos38"], coords["strand37to38"], None)
                        else:
                            # Identity lookup in the input build; do both chains only on a modern hit.
                            qc["selection_coordinate_direct_grch37"] += 1
                            state = (chrom37, pos37, None, None)
                        qc["selection_coordinates_indexed"] += 1
                    except MappingError as exc:
                        state = (None, None, None, exc.reason)
                        qc["selection_coordinate_rejected_" + exc.reason] += 1
                        coord_rejects.writerow(dict(chrom37=chrom37, pos37=pos37, reason=exc.reason))
                    db.execute("INSERT INTO coordinates(chrom37,pos37,source_chrom,source_pos,strand,reason) "
                               "VALUES(?,?,?,?,?,?)", (*key, *state))
                previous_key, previous_state = key, state
                if state[3] is not None:
                    qc["selection_snv_rows_without_usable_coordinate"] += 1
                    continue
                # Fixed-order values avoid copying repeated column names millions of times.
                # Numeric/filter/posterior validity is intentionally NOT a prefilter.
                payload = json.dumps([row[k] for k in SELECTION_FIELDS], separators=(",", ":"))
                payload_bytes += len(payload)  # ensure_ascii=True, hence length is exact UTF-8 byte count.
                db.execute("INSERT INTO contrasts VALUES(?,?,?,?)", (*key, row_no, payload))
                qc["selection_snv_rows_indexed"] += 1
        db.commit()
        db.execute("CREATE UNIQUE INDEX coordinates_source ON coordinates(source_chrom,source_pos) WHERE reason IS NULL")
        db.commit()
        cache = _CoordinateCache(db)
        after_index_memory = _process_memory()
        index_elapsed = time.monotonic() - start
        rejects_handle, rejects = tsv_writer(out / "mapping_rejections.tsv.gz", ["chrom", "pos", "ref", "reason"])
        handle, writer = tsv_writer(out / "selection_intersection.tsv.gz", JOIN_FIELDS)
        with rejects_handle, handle:
            for site, alleles in _modern_sites(modern["shards"]):
                qc["modern_multi_alt_sites"] += 1
                if qc["modern_multi_alt_sites"] % 1000000 == 0:
                    db.commit()
                    print(f"[modern intersect] {qc['modern_multi_alt_sites']:,} sites; "
                          f"{qc['modern_coordinate_hits']:,} coordinate hits", flush=True)
                chrom, pos = canonical_chrom(site["chrom"]), int(site["pos"])
                if chrom not in AUTOSOMES:
                    qc["modern_non_autosomal_sites"] += 1
                    continue
                if not cache.contains(chrom, pos):
                    qc["modern_sites_without_coordinate_hit"] += 1
                    continue
                qc["modern_coordinate_hits"] += 1
                chrom37, pos37, strand, seen = db.execute(
                    "SELECT chrom37,pos37,strand,modern_seen FROM coordinates WHERE source_chrom=? "
                    "AND source_pos=? AND reason IS NULL", (chrom, pos)).fetchone()
                qc["fullmap_calls"] += 1
                try:
                    if source_build == "GRCh38":
                        coords = dict(chrom38=chrom, pos38=pos, chrom37=chrom37, pos37=pos37,
                                      strand37to38=strand)
                        mapped = mapper.map_from_coordinates(coords, site["ref"], source_build)
                    else:
                        mapped = mapper.map(chrom, pos, site["ref"], source_build)
                except MappingError as exc:
                    qc["hit_mapping_" + exc.reason] += 1
                    rejects.writerow({**{k: site[k] for k in ("chrom", "pos", "ref")}, "reason": exc.reason})
                    continue
                if seen:
                    raise ValueError("Two modern sites map to the same ancient coordinate; ambiguous source merge")
                db.execute("UPDATE coordinates SET modern_seen=1 WHERE chrom37=? AND pos37=?", (chrom37, pos37))
                qc["modern_fullmap_passed"] += 1
                qc["hit_reference_changed_sites"] += mapped["reference_changed"]
                oriented = [dict(a, alt=orient(a["alt"], mapped["source_to38"])) for a in alleles]
                matched_ref = False
                for (payload,) in db.execute("SELECT payload FROM contrasts WHERE chrom37=? AND pos37=? ORDER BY row_no",
                                             (chrom37, pos37)):
                    row = dict(zip(SELECTION_FIELDS, json.loads(payload)))
                    if row["REF"] != mapped["ref37"]:
                        qc["selection_ref37_mismatch"] += 1
                        continue
                    matched_ref = True
                    qc["matched_contrasts"] += 1
                    writer.writerows(_contrast_rows(row, mapped, site, oriented, qc))
                qc["modern_hit_sites_with_ref_matched_summary" if matched_ref else
                   "modern_validated_hits_without_ref_matched_summary"] += 1
        db.commit()
        memory = dict(before=before_memory, after_coordinate_index=after_index_memory,
                      after_modern_stream=_process_memory(), coordinate_cache=cache.metrics(),
                      sqlite_page_cache_target_bytes=64 * 1024 * 1024,
                      note="Array budget excludes small Python bookkeeping, FASTA caches and both pyliftover chain indexes. "
                           "RSS is measured for the whole process; lifetime peak can include preceding stages.")
        db.close()
        os.replace(tmp, db_path)
    finally:
        db.close()
        if own_mapper and mapper is not None:
            mapper.close()
    performance = dict(index_bytes=db_path.stat().st_size, selection_payload_bytes=payload_bytes,
                       selection_index_seconds=index_elapsed,
                       coordinate_hits=qc["modern_coordinate_hits"], fullmap_calls=qc["fullmap_calls"],
                       selection_coordinate_map_calls=qc["selection_coordinate_map_calls"], memory=memory)
    return finish_stage(out, identity, dict(qc=dict(qc), elapsed_seconds=time.monotonic() - start,
                         strategy="selection-first", performance=performance,
                         qc_scope="All canonical published SNVs are coordinate-indexed without AF, posterior, S or FILTER preselection. "
                         "All modern multi-ALT rows are streamed. Full mapping/REF QC and mapping_rejections apply ONLY to "
                         "modern coordinate hits. Nonhits are untested, not mapping failures; no global mapped-without-summary count. "
                         "For GRCh37 modern input the coordinate index is direct; reciprocal chains are checked on hits. "
                         "Modern duplicate-coordinate rejection is likewise restricted to validated coordinate hits.",
                         intersection_path=str(out / "selection_intersection.tsv.gz"), index_path=str(db_path),
                         universe="modern multi-ALT sites; no all-site enrichment denominator; "
                         "extra must be a directly observed passing source ALT, never inferred from the source REF"),
                        ["selection_intersection.tsv.gz", "mapping_rejections.tsv.gz",
                         "selection_coordinate_rejections.tsv.gz", "selection_index.sqlite"])


def join_selection(modern_out, selection_path, out_dir, fasta38, fasta37,
                   chain38to37, chain37to38, *, mapper=None, strategy="selection-first"):
    """Intersect complete canonical selection SNVs with complete modern multi-ALT tables.

    The default uses a reciprocal coordinate prefilter and validates FASTA/REF on
    modern hits only. ``modern-first`` preserves the slower reference algorithm.
    Their accepted allele-contrast multisets are identical for valid modern runs;
    row order and mapping-QC scope differ. Both retain the same explicit observed
    source-ALT restriction, numeric validation and selection direction semantics.
    Changing strategy requires a new join directory, not a new modern scan.
    """
    methods = {"selection-first": _join_selection_first, "modern-first": _join_modern_first}
    if strategy not in methods:
        raise ValueError("join strategy must be selection-first or modern-first")
    return methods[strategy](modern_out, selection_path, out_dir, fasta38, fasta37,
                             chain38to37, chain37to38, mapper=mapper)


def _match_stratum(row):
    # Explicit coarse bins, identical between case and control. Matching the
    # extra AF is for FUNCTION comparisons only, not AF-as-outcome analyses.
    af = float(row["selection_pooled_af"])
    maf = min(af, 1 - af)
    maf_bin = min(9, int(maf / .05))
    extra_af = max(float(row["extra_af"]), 1e-12)
    an = max(float(row["extra_an"]), 1.)
    log_af_bin = math.floor(math.log10(extra_af) * 2)
    log_an_bin = math.floor(math.log10(an) * 4)
    coding = ("unknown" if row["consequence_terms"] in ("", ".", "NA") else
              str(int(any(t in row["consequence_terms"] for t in
                          ("missense", "synonymous", "stop_", "splice_", "frameshift")))))
    return f"{row['canonical_context']}|maf{maf_bin}|af{log_af_bin}|an{log_an_bin}|coding{coding}"


def _best_extra_rows(path, seed):
    """The join emits contiguous extra rows for each precise selection contrast.

    Use the same highest-AC extra rule for cases and controls before sampling.
    The expected ordering is part of the join's published file contract.
    """
    rows = table_rows(path, JOIN_FIELDS)
    groups = itertools.groupby(rows, key=lambda r: (r["site_id"], r["selection_id"]))
    for _, group in groups:
        yield max(group, key=lambda r: (float(r["extra_ac"]),
                  stable_hash([seed, r["site_id"], r["selection_id"], r["extra_alt"]])))


def _physical_block(row, size):
    return f"{row['chrom37']}:{(int(row['pos37']) - 1) // size}"


def select_candidates(intersection_path, out_dir, *, max_sites=2000, strong_posterior=.99,
                      weak_posterior=.5, min_extra_ac=5, block_bp=1000000,
                      exclude_mhc=True, allow_filter_warnings=False, seed=20261021):
    """Deterministic, bounded, pre-functional-score case/control selection.

    One strong site per proximity block, one highest-count extra per contrast,
    one matching weak-signal site from a different unused block. Proximity blocks
    are NOT LD-clumped loci; the weak-signal group is NOT proven neutral.
    Unmatched cases are retained with an explicit status and are not paired tests.
    """
    if max_sites < 2 or block_bp < 1 or min_extra_ac < 1 or not 0 <= weak_posterior < strong_posterior <= 1:
        raise ValueError("Invalid selection budget, posterior thresholds, block size or AC floor")
    parameters = dict(max_sites=max_sites, strong_posterior=strong_posterior, weak_posterior=weak_posterior,
                      min_extra_ac=min_extra_ac, block_bp=block_bp, exclude_mhc=exclude_mhc,
                      allow_filter_warnings=allow_filter_warnings, seed=seed)
    identity = _stage_identity("candidate_design", dict(intersection=intersection_path), parameters)
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if old := completed_stage(out, identity):
        return old
    qc, strong_blocks = Counter(), {}
    summary = defaultdict(lambda: [0, 0., 0.])

    def eligible(row):
        qc["intersection_rows"] += 1
        if number(row["extra_ac"], 0) < min_extra_ac:
            qc["below_extra_ac"] += 1
            return False
        if row["selection_filter"] != "PASS" and not allow_filter_warnings:
            qc["selection_filter_excluded"] += 1
            return False
        if row["has_non_snv"] == "1" or row["region_flags"] not in (".", "", "none"):
            qc["complex_or_flagged_site_excluded"] += 1
            return False
        if row["chrom37"] == "6" and 25000000 <= int(row["pos37"]) <= 34000000 and exclude_mhc:
            qc["mhc_excluded"] += 1
            return False
        if float(row["s"]) == 0:
            qc["zero_direction_excluded"] += 1
            return False
        return True

    def rank(row):
        return (float(row["posterior"]), abs(float(row["x"])), float(row["extra_ac"]),
                stable_hash([seed, row["site_id"], row["selection_id"], row["extra_alt"]]))

    for row in _best_extra_rows(intersection_path, seed):
        if not eligible(row):
            continue
        p = float(row["posterior"])
        group = "strong" if p >= strong_posterior else "weak_signal" if p <= weak_posterior else "intermediate"
        key = (group, row["canonical_context"])
        summary[key][0] += 1
        summary[key][1] += float(row["extra_af"])
        summary[key][2] += float(row["extra_ac"])
        qc["eligible_" + group] += 1
        if group == "strong":
            block = _physical_block(row, block_bp)
            if block not in strong_blocks or rank(row) > rank(strong_blocks[block]):
                strong_blocks[block] = row
    cases = sorted(strong_blocks.values(), key=rank, reverse=True)[:max_sites // 2]
    qc["strong_proximity_blocks"] = len(strong_blocks)
    qc["case_sites_selected"] = len(cases)
    # Keep a bounded deterministic pool of controls only for strata actually
    # needed. Reservoirs are capped at max_sites rows per stratum; total <=
    # max_sites * number of case strata (not all genome rows).
    wanted = Counter(_match_stratum(row) for row in cases)
    pools = defaultdict(list)
    case_sites = {row["site_id"] for row in cases}
    for row in _best_extra_rows(intersection_path, seed):
        # Do not double-count QC from the first pass.
        if float(row["posterior"]) > weak_posterior or row["site_id"] in case_sites:
            continue
        if number(row["extra_ac"], 0) < min_extra_ac or float(row["s"]) == 0:
            continue
        if (row["selection_filter"] != "PASS" and not allow_filter_warnings) or row["has_non_snv"] == "1" or row["region_flags"] not in (".", "", "none"):
            continue
        if exclude_mhc and row["chrom37"] == "6" and 25000000 <= int(row["pos37"]) <= 34000000:
            continue
        key = _match_stratum(row)
        if key not in wanted:
            continue
        # Smallest stable hash is an outcome-independent deterministic sample.
        h = int(stable_hash([seed, row["site_id"], row["selection_id"], row["extra_alt"]]), 16)
        item = (-h, row["site_id"] + row["selection_id"] + row["extra_alt"], row)
        pool = pools[key]
        capacity = min(max_sites, max(40, wanted[key] * 20))
        if len(pool) < capacity:
            heapq.heappush(pool, item)
        elif item[:2] > pool[0][:2]:
            heapq.heapreplace(pool, item)
    sorted_pools = {key: [item[2] for item in sorted(values, reverse=True)] for key, values in pools.items()}
    used_blocks = {_physical_block(row, block_bp) for row in cases}
    used_sites = set(case_sites)
    selected, matched = [], 0
    for index, case in enumerate(cases, 1):
        key = _match_stratum(case)
        control = next((row for row in sorted_pools.get(key, [])
                        if row["site_id"] not in used_sites and _physical_block(row, block_bp) not in used_blocks), None)
        match_set = f"pair{index:05d}" if control else "."
        for row, group in ((case, "strong"), (control, "weak_signal")):
            if row is None:
                continue
            block = _physical_block(row, block_bp)
            selected.append({**row, "group": group, "match_set": match_set, "match_stratum": key,
                             "locus_block": block, "selection_evidence": "published_strong" if group == "strong" else "weak_signal_not_neutral",
                             "comparison_universe": "matched_modern_multi_alt_functional_comparison"})
            used_blocks.add(block)
            used_sites.add(row["site_id"])
        matched += int(control is not None)
    qc["matched_pairs"] = matched
    qc["unmatched_cases"] = len(cases) - matched
    qc["selected_sites"] = len(selected)
    handle, writer = tsv_writer(out / "candidates.tsv", CANDIDATE_FIELDS)
    with handle:
        writer.writerows(selected)
    handle, writer = tsv_writer(out / "selection_context_summary.tsv", ["selection_group", "canonical_context", "extra_contrast_rows", "mean_extra_af", "mean_extra_ac"])
    with handle:
        for (group, context), (n, af_sum, ac_sum) in sorted(summary.items()):
            writer.writerow(dict(selection_group=group, canonical_context=context, extra_contrast_rows=n,
                                 mean_extra_af=af_sum / n, mean_extra_ac=ac_sum / n))
    return finish_stage(out, identity, dict(qc=dict(qc), candidate_path=str(out / "candidates.tsv"),
                         caveat="Conditional on modern multi-ALT sites; physical blocks do not establish LD independence. "
                                "Extra AF is matched for functional comparisons, so do not test AF differences in paired candidates."),
                        ["candidates.tsv", "selection_context_summary.tsv"])


def index_target_variants(candidate_path, pvar_path, high_quality_path, out_dir):
    """One sequential PVAR pass; preserve original indices for targeted PGEN reads.

    Targets include each ancient anchor and the extra base versus ancient REF.
    A missing extra is 'not_represented_in_panel', NEVER dosage zero.
    """
    identity = _stage_identity("panel_index", dict(candidates=candidate_path, pvar=pvar_path,
                                                   high_quality=high_quality_path), {})
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if old := completed_stage(out, identity):
        return old
    targets, links = {}, []
    for row in table_rows(candidate_path, ("site_id", "chrom37", "pos37", "ref37", "alt37", "extra_alt", "strand37to38")):
        extra37 = orient(row["extra_alt"], row["strand37to38"])
        for role, alt in (("anchor", row["alt37"]), ("extra", extra37)):
            key = (canonical_chrom(row["chrom37"]), int(row["pos37"]), row["ref37"], alt)
            if alt == key[2]:
                raise ValueError("Extra base equals ancient REF; invalid harmonization")
            targets.setdefault(key, dict(variant_id="_".join(map(str, key)), chrom=key[0], pos=key[1],
                                         ref=key[2], alt=key[3], pvar_index=".", status="not_represented_in_panel"))
            links.append(dict(site_id=row["site_id"], role=role, target_id="_".join(map(str, key)),
                              base38=row["anchor_alt"] if role == "anchor" else row["extra_alt"]))
    hq, seen = set(), set()
    for row in (table_rows(high_quality_path, ("CHROM", "POS", "ID", "REF", "ALT")) if targets else ()):
        key = (canonical_chrom(row["CHROM"]), int(row["POS"]), row["REF"], row["ALT"])
        if key in targets:
            hq.add(key)
    rows = 0
    for index, row in enumerate(table_rows(pvar_path, ("CHROM", "POS", "ID", "REF", "ALT")) if targets else ()):
        rows += 1
        key = (canonical_chrom(row["CHROM"]), int(row["POS"]), row["REF"], row["ALT"])
        if "," in row["ALT"]:
            for alt in row["ALT"].split(","):
                multi_key = (*key[:3], alt)
                if multi_key in targets:
                    # pgenlib allele_idx=1 cannot represent arbitrary multiallelic
                    # rows in this biallelic module. Do not call this absent.
                    if multi_key in seen:
                        raise ValueError(f"Duplicate exact/multiallelic PVAR contrast: {multi_key}")
                    seen.add(multi_key)
                    targets[multi_key].update(variant_id=row["ID"] + ":" + alt, pvar_index=".",
                                              status="represented_multiallelic_unsupported")
        if key in targets:
            if key in seen:
                raise ValueError(f"Duplicate exact PVAR contrast: {key}")
            seen.add(key)
            targets[key].update(variant_id=row["ID"], pvar_index=index,
                                status="assayed_pass" if key in hq else "not_in_high_quality")
        if rows % 5000000 == 0:
            print(f"[PVAR] {rows:,} records streamed; {len(seen):,}/{len(targets):,} requested contrasts located", flush=True)
    fields = ["variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"]
    handle, writer = tsv_writer(out / "selected_variants.tsv", fields)
    with handle:
        writer.writerows(sorted(targets.values(), key=lambda r: (int(r["chrom"]), r["pos"], r["alt"])))
    link_by_id = {"_".join(map(str, k)): v for k, v in targets.items()}
    handle, writer = tsv_writer(out / "site_variant_links.tsv", ["site_id", "role", "target_id", "base38", "variant_id", "status"])
    with handle:
        for link in links:
            target = link_by_id[link["target_id"]]
            writer.writerow({**link, "variant_id": target["variant_id"], "status": target["status"]})
    counts = Counter(row["status"] for row in targets.values())
    return finish_stage(out, identity, dict(pvar_records=rows if targets else None,
                         scan_skipped_no_targets=not bool(targets), target_counts=dict(counts),
                         selected_variants_path=str(out / "selected_variants.tsv"),
                         links_path=str(out / "site_variant_links.tsv")), ["selected_variants.tsv", "site_variant_links.tsv"])

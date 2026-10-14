#!/usr/bin/env python3
"""Frozen NCBI taxonomy, bounded lookups, exclusive routes and index audits."""
from __future__ import annotations
import csv
import fcntl
import hashlib
import json
import sqlite3
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
import re
import subprocess
import uuid

ROOTS = {'human': 9606, 'animal': 33208, 'plant': 33090,
         'fungi': 4751, 'bacteria': 2, 'archaea': 2157, 'virus': 10239}
BASE_GROUPS = tuple(ROOTS) + ('other_eukaryote', 'ambiguous', 'technical', 'other', 'unclassified')
GROUPS = tuple(g for g in ROOTS if g != 'virus') + ('virus_DNA', 'virus_RNA', 'virus_unknown') + BASE_GROUPS[7:]
TAXONOMY_SCHEMA = 'metag-taxonomy-v2'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def taxonomy_identity(db: Path) -> dict:
    out = {}
    for name in ('nodes.dmp', 'names.dmp', 'merged.dmp', 'delnodes.dmp'):
        p = db / 'taxonomy' / name
        if name in ('nodes.dmp', 'names.dmp') and not p.is_file():
            raise RuntimeError(f'Missing same-build taxonomy: {p}')
        out[name] = sha256(p) if p.is_file() else None
    return out


def snapshot(db: Path) -> dict:
    """Index identity uses stat; taxonomy uses content hashes, including optional files."""
    out = {'taxonomy': taxonomy_identity(db)}
    for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d'):
        p = db / name
        if not p.is_file():
            raise RuntimeError(f'Missing {p}')
        s = p.stat()
        out[name] = {'size': s.st_size, 'mtime_ns': s.st_mtime_ns}
    return out


def atomic_json(path, value):
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.partial')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def prepare_sqlite(db: Path) -> Path:
    """Content-addressed, atomically published cache; never hot-update a reader's DB."""
    identity = taxonomy_identity(db)
    key = hashlib.sha256(json.dumps([TAXONOMY_SCHEMA, identity], sort_keys=True).encode()).hexdigest()
    folder = db / '.metag_taxonomy'
    folder.mkdir(exist_ok=True)
    target = folder / f'{key}.sqlite'
    with (folder / '.build.lock').open('a') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if target.exists():
            return target
        tmp = target.with_suffix('.' + uuid.uuid4().hex + '.partial')
        con = sqlite3.connect(tmp)
        try:
            con.executescript('PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-65536; '
                'CREATE TABLE taxa(taxid INTEGER PRIMARY KEY,parent INTEGER,rank TEXT,name TEXT DEFAULT ""); '
                'CREATE TABLE merged(old INTEGER PRIMARY KEY,new INTEGER); '
                'CREATE TABLE deleted(taxid INTEGER PRIMARY KEY); '
                'CREATE TABLE names(name TEXT,taxid INTEGER,kind TEXT); '
                'CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);')
            def rows(name):
                with (db / 'taxonomy' / name).open() as f:
                    for line in f:
                        yield [s.strip() for s in line.split('|')]
            con.executemany('INSERT INTO taxa(taxid,parent,rank) VALUES (?,?,?)',
                            ((int(x[0]), int(x[1]), x[2]) for x in rows('nodes.dmp')))
            # Exact scientific names and official taxonomy synonyms only; no fuzzy/common-name guesses.
            con.executemany('INSERT INTO names VALUES (?,?,?)',
                            ((x[1], int(x[0]), x[3]) for x in rows('names.dmp')
                             if x[3] in ('scientific name', 'synonym', 'equivalent name')))
            con.execute('CREATE INDEX names_lookup ON names(name)')
            con.execute('CREATE INDEX names_taxid ON names(taxid)')
            con.execute('UPDATE taxa SET name=COALESCE((SELECT name FROM names WHERE names.taxid=taxa.taxid '
                        'AND kind="scientific name" LIMIT 1), "")')
            if identity['merged.dmp']:
                con.executemany('INSERT INTO merged VALUES (?,?)',
                                ((int(x[0]), int(x[1])) for x in rows('merged.dmp')))
            if identity['delnodes.dmp']:
                con.executemany('INSERT INTO deleted VALUES (?)', ((int(x[0]),) for x in rows('delnodes.dmp')))
            con.executemany('INSERT INTO metadata VALUES (?,?)',
                            [('schema', TAXONOMY_SCHEMA), ('identity', json.dumps(identity, sort_keys=True))])
            con.commit()
            if taxonomy_identity(db) != identity:
                raise RuntimeError('Taxonomy changed during cache build')
        except BaseException:
            con.close()
            tmp.unlink(missing_ok=True)
            raise
        con.close()
        tmp.replace(target)
    return target


def slug(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('._')[:90] or 'unnamed'


class Taxonomy:
    def __init__(self, sqlite_path: Path, virus_map: Path | None = None):
        self.con = sqlite3.connect(sqlite_path.resolve().as_uri() + '?mode=ro&immutable=1', uri=True)
        self.con.execute('PRAGMA cache_size=-32768')
        self.identity = json.loads(self.con.execute('SELECT value FROM metadata WHERE key="identity"').fetchone()[0])
        self.human_ancestors = set(self.lineage(9606))
        self.group_ancestors = set()
        for tid in ROOTS.values():
            self.group_ancestors.update(self.lineage(tid))
        self.virus = None
        if virus_map is not None:
            from virus_metadata import VirusGenomeResolver
            try:
                self.virus = VirusGenomeResolver(virus_map, self)
            except BaseException:
                self.con.close()
                raise

    @lru_cache(maxsize=100000)
    def canonical_taxid(self, tid: int) -> int:
        if tid == 0:
            return 0
        seen = set()
        while tid not in seen:
            seen.add(tid)
            if self.con.execute('SELECT 1 FROM taxa WHERE taxid=?', (tid,)).fetchone():
                return tid
            m = self.con.execute('SELECT new FROM merged WHERE old=?', (tid,)).fetchone()
            if m:
                tid = m[0]
                continue
            deleted = self.con.execute('SELECT 1 FROM deleted WHERE taxid=?', (tid,)).fetchone()
            raise RuntimeError(f'Taxid {tid} {"deleted" if deleted else "absent"} from same-build taxonomy')
        raise RuntimeError(f'Merged taxid cycle at {tid}')

    @lru_cache(maxsize=100000)
    def node(self, tid: int) -> tuple:
        return self.con.execute('SELECT parent,rank,name FROM taxa WHERE taxid=?',
                                (self.canonical_taxid(tid),)).fetchone()

    @lru_cache(maxsize=50000)
    def lineage(self, tid: int) -> tuple:
        tid = self.canonical_taxid(tid)
        ids, seen = [], set()
        for _ in range(200):
            if tid in seen:
                raise RuntimeError(f'Taxonomy cycle at {tid}')
            seen.add(tid)
            ids.append(tid)
            parent = self.canonical_taxid(self.node(tid)[0])
            if tid == parent:
                return tuple(ids)
            tid = parent
        raise RuntimeError('Taxonomy depth > 200')

    @lru_cache(maxsize=100000)
    def base_route(self, tid: int) -> tuple:
        tid = self.canonical_taxid(tid)
        if tid == 0:
            return ('unclassified', 'taxid_0', 0, '')
        lineage = self.lineage(tid)
        if 9606 in lineage:
            group = 'human'
        elif tid in self.human_ancestors:
            group = 'ambiguous'
        elif 33208 in lineage:
            group = 'animal'
        elif 33090 in lineage:
            group = 'plant'
        elif 4751 in lineage:
            group = 'fungi'
        elif 2 in lineage:
            group = 'bacteria'
        elif 2157 in lineage:
            group = 'archaea'
        elif 10239 in lineage:
            group = 'virus'
        elif tid in self.group_ancestors:
            group = 'ambiguous'
        elif 2759 in lineage:
            group = 'other_eukaryote'
        elif 32630 in lineage:
            group = 'technical'
        else:
            group = 'other'
        species = next((t for t in lineage if self.node(t)[1] == 'species'), 0)
        if species:
            name = self.node(species)[2]
            return (group, f'{species}_{slug(name)}', species, name)
        return (group, f'_unresolved/{tid}_{slug(self.node(tid)[2])}', 0, '')

    @lru_cache(maxsize=100000)
    def annotation(self, tid: int) -> dict:
        base = self.base_route(tid)[0]
        ann = {'canonical_taxid': self.canonical_taxid(tid), 'base_group': base,
               'genome_class': 'not_applicable', 'genome_type_detail': '',
               'resolution_status': 'not_applicable', 'evidence_taxid': ''}
        if base == 'virus':
            if self.virus is None:
                raise RuntimeError('Virus routing requires a frozen --virus-map')
            ann.update(self.virus.resolve(tid))
        return ann

    @lru_cache(maxsize=100000)
    def route(self, tid: int) -> tuple:
        group, bucket, species, name = self.base_route(tid)
        if group == 'virus':
            group = 'virus_' + self.annotation(tid)['genome_class']
        return group, bucket, species, name

    def close(self):
        if self.virus:
            self.virus.close()
        for name in ('node', 'lineage', 'route', 'base_route', 'canonical_taxid', 'annotation'):
            getattr(self, name).cache_clear()
        self.con.close()


def audit_database(db: Path, threads: int = 8, inspect_file: Path | None = None,
                   virus_map: Path | None = None, allow_incomplete_virus: bool = False) -> dict:
    with (db / '.metag.audit.lock').open('a') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _audit_database(db, threads, inspect_file, virus_map, allow_incomplete_virus)


def _audit_database(db, threads, inspect_file, virus_map, allow_incomplete_virus):
    initial = snapshot(db)
    sqlite_file = prepare_sqlite(db)
    inspect = inspect_file or db / 'metag.inspect.tsv'
    tax = Taxonomy(sqlite_file, virus_map)
    try:
        if inspect_file is None:
            tmp = inspect.with_suffix('.' + uuid.uuid4().hex + '.partial')
            cmd = ['kraken2-inspect', '--db', str(db), '--threads', str(threads), '--report-zero-counts']
            with (db / 'metag.inspect.log').open('w') as log, tmp.open('w') as report:
                subprocess.run(cmd, check=True, stderr=log, stdout=report)
            tmp.replace(inspect)
        counts, spp = defaultdict(lambda: [0, 0]), defaultdict(set)
        vc, vs = defaultdict(lambda: [0, 0]), defaultdict(set)
        reasons = defaultdict(int)
        seen = set()
        with inspect.open() as f:
            for line in f:
                if not line.strip():
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) != 6:
                    raise RuntimeError('Expected SIX-column k2 inspect report')
                tid, direct = int(p[4]), int(p[2])
                if direct < 0:
                    raise RuntimeError('Negative direct minimizer count')
                if not tid or direct == 0:
                    continue
                canonical = tax.canonical_taxid(tid)
                if canonical in seen:
                    raise RuntimeError(f'Duplicate inspect taxid {tid}; refusing double counting')
                seen.add(canonical)
                group, _, species, _ = tax.base_route(tid)
                counts[group][0] += 1
                counts[group][1] += direct
                if species:
                    spp[group].add(species)
                if group == 'virus' and tax.virus:
                    ann = tax.annotation(tid)
                    kind = ann['genome_class']
                    vc[kind][0] += 1
                    vc[kind][1] += direct
                    reasons[ann['resolution_status']] += 1
                    if species:
                        vs[kind].add(species)
        coverage = [dict(group=g, taxa_with_direct_minimizers=counts[g][0],
                         species_with_direct_minimizers=len(spp[g]), direct_minimizers=counts[g][1])
                    for g in BASE_GROUPS]
        def save(name, rows):
            with (db / name).open('w', newline='') as f:
                w = csv.DictWriter(f, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
        save('metag.coverage.tsv', coverage)
        virus_coverage = []
        if tax.virus:
            assert counts['virus'][1] == sum(vc[k][1] for k in ('DNA', 'RNA', 'unknown'))
            den_nodes, den_min = counts['virus']
            for kind in ('DNA', 'RNA', 'unknown'):
                virus_coverage.append(dict(genome_class=kind, direct_minimizers=vc[kind][1],
                    taxa_with_direct_minimizers=vc[kind][0], species_with_direct_minimizers=len(vs[kind]),
                    denominator_virus_nodes=den_nodes, denominator_virus_direct_minimizers=den_min,
                    node_fraction=vc[kind][0] / den_nodes if den_nodes else 'NA',
                    minimizer_fraction=vc[kind][1] / den_min if den_min else 'NA',
                    annotated_node_fraction=(den_nodes-vc['unknown'][0])/den_nodes if den_nodes else 'NA',
                    mapping_sha256=tax.virus.identity['mapping_sha256']))
            save('metag.virus_coverage.tsv', virus_coverage)
            tax.virus.assert_unchanged()
        missing = [k for k in ('DNA', 'RNA') if vc[k][1] == 0]
        out = {'schema': 'v2', 'snapshot': initial, 'taxonomy_sqlite': str(sqlite_file.resolve()),
               'taxonomy_sha256': tax.identity, 'coverage': coverage,
               'missing_groups': [g for g in ROOTS if counts[g][1] == 0],
               'virus_mapping': tax.virus.identity if tax.virus else None,
               'virus_coverage': virus_coverage, 'virus_resolution_nodes': dict(reasons),
               'missing_virus_classes': missing, 'virus_types_tested': bool(tax.virus and not missing),
               'virus_coverage_reason': ('mapping_not_supplied' if not tax.virus else
                   ('unannotated_index_nodes_prevent_type_confirmation' if vc['unknown'][0] and missing else
                    'no_index_support_for_missing_types' if missing else 'both_types_have_index_support')),
               'exploratory_virus_coverage': bool(tax.virus and missing and allow_incomplete_virus),
               'inspect_sha256': sha256(inspect),
               'warning': 'Nonzero coverage is NOT complete species coverage; minimizers are not genomes.'}
        if snapshot(db) != initial:
            raise RuntimeError('Database changed during audit')
        atomic_json(db / 'metag.audit.json', out)
        if tax.virus and missing and not allow_incomplete_virus:
            raise RuntimeError(f'Virus DNA/RNA index coverage insufficient: {missing}. See metag.virus_coverage.tsv; '
                               'use --allow-incomplete-virus only for documented exploration.')
        return out
    finally:
        tax.close()

"""Preserve indexed reference hits without inferring unique biological origins."""
from metacache_inputs import original_accession
import csv
import gzip
import itertools
import json
import os
from pathlib import Path
import sqlite3
import subprocess
from contextlib import closing

FIELDS=['query_id','query_header','sequence_accession','reference_file',
        'assembly_accession','species_taxid','feature_hits']


def load_inventory(library):
    with (library/'reference/inventory.tsv').open(newline='') as handle:
        return {row['filename']:row for row in csv.DictReader(handle,delimiter='\t')}


def extract(raw, mapping, evidence, inventory, paired=False):
    """Validate a complete native query and separate reference JSON from merge input."""
    columns=None
    count=0
    summary=None
    ordered=True
    with raw.open() as source, mapping.open('w') as dest, gzip.open(evidence,'wt',compresslevel=1,newline='') as hits:
        writer=csv.DictWriter(hits,fieldnames=FIELDS,delimiter='\t')
        writer.writeheader()
        for line in source:
            if line.startswith('# TABLE_LAYOUT:'):
                columns=[x.strip() for x in line.split(':',1)[1].split('|')]
                required={'query_id','query_header','reference_hits','top_hits','rank','taxid'}
                if not required.issubset(columns):
                    raise RuntimeError('MetaCache lacks -refhits output; rerun ./install.sh')
                kept=[x for x in columns if x!='reference_hits']
                dest.write('# TABLE_LAYOUT: '+'\t|\t'.join(kept)+'\n')
            elif line.startswith('#'):
                dest.write(line)
                if line.startswith('# queries:'):
                    summary=int(line.split(':',1)[1].strip())
            elif line.strip():
                if columns is None:
                    raise RuntimeError('MetaCache output is missing its table header')
                values=[x.strip() for x in line.rstrip('\n').split('|')]
                if len(values)!=len(columns):
                    raise RuntimeError('Malformed MetaCache evidence columns')
                row=dict(zip(columns,values))
                count+=1
                qid=int(row['query_id'])
                if qid<1:
                    raise RuntimeError(f'Invalid MetaCache query ID: {qid}')
                if qid!=count:
                    ordered=False
                dest.write('\t|\t'.join(row[key] for key in kept)+'\n')
                candidates=json.loads(row['reference_hits'])
                if not isinstance(candidates,list):
                    raise RuntimeError('Malformed MetaCache reference evidence')
                for accession,sourcefile,species,hits_count in candidates:
                    filename=Path(sourcefile).name
                    if filename not in inventory:
                        raise RuntimeError(f'Indexed reference absent from active inventory: {filename}')
                    if int(hits_count)<=0 or int(species)<0:
                        raise RuntimeError('Invalid reference feature count/taxid')
                    writer.writerow(dict(query_id=qid,query_header=row['query_header'],
                        sequence_accession=original_accession(accession),reference_file=filename,
                        assembly_accession=inventory[filename]['assembly_accession'],
                        species_taxid=species,feature_hits=hits_count))
    # Upstream's '# queries' summary counts physical reads, but paired table rows
    # and query IDs each represent one fragment.
    expected=count*(2 if paired else 1)
    if columns is None or summary!=expected:
        raise RuntimeError(f'Incomplete MetaCache query output: {count} fragments, expected summary={expected}, got {summary}')
    if not ordered:
        # CPU/GPU workers may finish batches out of order. Sort only the compact
        # mapping, never the much larger reference JSON. GNU sort spills to disk;
        # stable sorting keeps all comment/header lines in their original order.
        print(f'Ordering {count} MetaCache mappings by query ID',flush=True)
        try:
            subprocess.run(['sort','--stable','--numeric-sort','--key=1,1',
                            '--buffer-size=64M','--temporary-directory='+str(mapping.parent),
                            '--output='+str(mapping),str(mapping)],
                           check=True,env={**os.environ,'LC_ALL':'C'})
        except (OSError,subprocess.CalledProcessError) as exc:
            raise RuntimeError(f'Cannot order MetaCache mappings: {exc}') from exc
        with mapping.open() as handle:
            index=0
            for line in handle:
                if line.startswith('#') or not line.strip():continue
                index+=1
                qid=int(line.split('|',1)[0])
                if qid!=index:
                    raise RuntimeError(f'MetaCache query IDs are incomplete or duplicated: '
                                       f'expected {index}, got {qid} after sorting')
        if index!=count:
            raise RuntimeError(f'MetaCache mapping count changed during sorting: {count} -> {index}')
    return count


def build_index(paths, output):
    """Disk-backed join across partitions; memory is bounded independently of read count."""
    temp=output.with_name(output.name+'.partial')
    temp.unlink(missing_ok=True)
    with closing(sqlite3.connect(temp)) as db:
        db.execute('PRAGMA journal_mode=OFF')
        db.execute('PRAGMA synchronous=OFF')
        db.execute('PRAGMA cache_size=-32768')
        db.execute('PRAGMA temp_store=FILE')
        db.executescript('''
            CREATE TABLE refs(id INTEGER PRIMARY KEY, sequence_accession TEXT, reference_file TEXT,
              assembly_accession TEXT, species_taxid INTEGER,
              UNIQUE(sequence_accession,reference_file,assembly_accession,species_taxid));
            CREATE TABLE reads(qid INTEGER PRIMARY KEY, header TEXT NOT NULL);
            CREATE TABLE hits(qid INTEGER, refid INTEGER, features INTEGER);
            CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT);
        ''')
        cache={}
        count=0
        for path in paths:
            with gzip.open(path,'rt',newline='') as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    qid=int(row['query_id'])
                    db.execute('''INSERT INTO reads VALUES(?,?) ON CONFLICT(qid) DO UPDATE SET
                               header=CASE WHEN reads.header=excluded.header THEN reads.header ELSE NULL END''',
                               (qid,row['query_header']))
                    key=(row['sequence_accession'],row['reference_file'],row['assembly_accession'],int(row['species_taxid']))
                    refid=cache.get(key)
                    if refid is None:
                        db.execute('INSERT OR IGNORE INTO refs(sequence_accession,reference_file,assembly_accession,species_taxid) VALUES(?,?,?,?)',key)
                        refid=db.execute('SELECT id FROM refs WHERE sequence_accession=? AND reference_file=? AND assembly_accession=? AND species_taxid=?',key).fetchone()[0]
                        if len(cache)>=100000:cache.clear()
                        cache[key]=refid
                    db.execute('INSERT INTO hits VALUES(?,?,?)',(qid,refid,int(row['feature_hits'])))
                    count+=1
                    if count%50000==0:db.commit()
        db.execute('CREATE INDEX hits_qid ON hits(qid)')
        db.execute('INSERT INTO metadata VALUES(?,?)',('status','complete'))
        db.commit()
    temp.replace(output)
    return count


def grouped(path):
    if not path.is_file():
        raise RuntimeError(f'Per-read reference evidence missing: {path}')
    with closing(sqlite3.connect(path)) as db:
        if db.execute("SELECT value FROM metadata WHERE key='status'").fetchone()!=('complete',):
            raise RuntimeError('Per-read reference evidence is incomplete')
        query='''SELECT h.qid,r.header,f.sequence_accession,f.reference_file,
                        f.assembly_accession,f.species_taxid,h.features
                 FROM hits h JOIN reads r ON r.qid=h.qid JOIN refs f ON f.id=h.refid
                 ORDER BY h.qid,h.refid'''
        for qid,rows in itertools.groupby(db.execute(query),key=lambda row:row[0]):
            yield qid,[dict(zip(FIELDS,row)) for row in rows]


def reference_status(rows):
    assemblies=sorted({row['assembly_accession'] for row in rows})
    if not assemblies:return 'no_retained_reference_hit',assemblies
    if len(assemblies)==1:return 'single_reported_assembly',assemblies
    return 'ambiguous_multiple_assemblies',assemblies

"""Assembly-aware adapters for the prepared ancient GRCh38 pfile.

The published HQ/selection tables remain GRCh37 annotations. Never interpret
positions encoded in a preserved variant ID as GRCh38 positions.
"""
from pathlib import Path
import csv
import fcntl
import hashlib
import json
import os
import re
import sqlite3
import tempfile


def _declared_build(value):
    """Recognize assembly labels without guessing from a containing directory."""
    text = str(value or '').strip().strip('"')
    aliases = {'37': 'GRCh37', 'b37': 'GRCh37', 'hg19': 'GRCh37',
               '38': 'GRCh38', 'b38': 'GRCh38', 'hg38': 'GRCh38'}
    if text.lower() in aliases:
        return aliases[text.lower()]
    match = re.search(r'(?<![A-Za-z0-9])GRCh(37|38)(?![0-9])', text, re.I)
    return 'GRCh' + match.group(1) if match else None


def panel_build(pvar):
    # The coordinate assembly belongs to the PVAR, not its reused PGEN link.
    path = Path(pvar).absolute()
    receipt = path.with_suffix('.prepare.json')
    if receipt.is_file():
        build = _declared_build(json.loads(receipt.read_text()).get('build'))
        if build:
            return build
    if path.is_file():
        with path.open() as handle:
            for line in handle:
                if not line.startswith('##'):
                    break
                if line.startswith(('##reference=', '##assembly=')):
                    build = _declared_build(line.split('=', 1)[1])
                    if build:
                        return build
    raise ValueError('PVAR requires an explicit assembly declaration or preparation receipt: ' + str(path))


def hq_in_panel_coordinates(pvar, hq, cache_root):
    """Join preserved ID/CHROM/REF/ALT; take POS only from the lifted PVAR.

    The disk index bounds RAM. Failed liftover rows (POS=-1) are not eligible.
    Duplicate HQ identities are collapsed; PVAR order/indices are not changed.
    """
    pvar, hq = Path(pvar), Path(hq)
    if panel_build(pvar) != 'GRCh38':
        return hq
    identity = {'schema': 1, 'build': 'GRCh38', 'join': 'ID,CHROM,REF,ALT',
                'inputs': {str(p.absolute()): [p.stat().st_size, p.stat().st_mtime_ns] for p in (pvar, hq)}}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
    folder = Path(cache_root) / ('ancient38-hq-' + key)
    folder.mkdir(parents=True, exist_ok=True)
    output, receipt = folder/'HighQualityVariants.GRCh38.tsv', folder/'complete.json'
    with (folder/'lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if receipt.is_file() and output.is_file():
            prior = json.loads(receipt.read_text())
            if prior['identity'] == identity and prior['bytes'] == output.stat().st_size:
                return output
        with tempfile.TemporaryDirectory(prefix='join-', dir=folder) as tmp:
            db = sqlite3.connect(str(Path(tmp)/'hq.sqlite'))
            try:
                db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
                db.execute('CREATE TABLE hq(id TEXT,chrom TEXT,ref TEXT,alt TEXT,PRIMARY KEY(id,chrom,ref,alt)) WITHOUT ROWID')
                with hq.open() as handle:
                    rows = csv.DictReader(handle, delimiter='\t')
                    if not {'#CHROM','POS','ID','REF','ALT'} <= set(rows.fieldnames or []):
                        raise ValueError('HQ table requires #CHROM POS ID REF ALT')
                    db.executemany('INSERT OR IGNORE INTO hq VALUES(?,?,?,?)',
                        ((r['ID'], r['#CHROM'].removeprefix('chr'),r['REF'],r['ALT']) for r in rows
                         if r.get('FILTER','PASS') in ('PASS','.')))
                db.commit()
                count = 0
                partial = Path(tmp)/'hq.tsv'
                with pvar.open() as src, partial.open('w') as dst:
                    rows = csv.DictReader((line for line in src if not line.startswith('##')), delimiter='\t')
                    if not {'#CHROM','POS','ID','REF','ALT'} <= set(rows.fieldnames or []):
                        raise ValueError('PVAR requires #CHROM POS ID REF ALT')
                    dst.write('#CHROM\tPOS\tID\tREF\tALT\n')
                    for r in rows:
                        if int(r['POS']) <= 0: continue
                        if db.execute('SELECT 1 FROM hq WHERE id=? AND chrom=? AND ref=? AND alt=?',
                           (r['ID'],r['#CHROM'].removeprefix('chr'),r['REF'],r['ALT'])).fetchone():
                            dst.write('\t'.join(r[k] for k in ('#CHROM','POS','ID','REF','ALT'))+'\n'); count += 1
                if not count: raise ValueError('No HQ variants match the lifted ancient PVAR')
                partial.replace(output)
                pending = folder/'complete.partial.json'
                pending.write_text(json.dumps({'identity':identity,'rows':count,'bytes':output.stat().st_size},indent=2)+'\n')
                pending.replace(receipt)
            finally:
                db.close()
    return output


def panel_annotation_rows(pvar, records, targets):
    """Read the 38 PVAR while retaining verified 37 annotation join keys.

    Existing allele-time tables explicitly carry both builds. Their 37 fields
    remain annotation coordinates, not a claim about the genotype file build.
    Only forward, reference-preserving mapped targets can match gen.sh outputs.
    Unmatched/failed rows are yielded unchanged in count with POS=-1, retaining
    the original PGEN row index even when a row cannot contribute a measurement.
    """
    if panel_build(pvar) != 'GRCh38':
        yield from records
        return
    mappings = {}
    for t in targets:
        if t.get('mapping_status', 'mapped') != 'mapped' or t.get('strand37to38') != '+': continue
        c = str(t['chrom']).removeprefix('chr'); p = int(t['pos'])
        old = (str(t['chrom37']).removeprefix('chr'), int(t['pos37']), t['ref37'])
        if c != old[0] or t['ref'] != old[2]: continue
        key = (c,p,t['ref'])
        if key in mappings and mappings[key] != old: raise ValueError('Ambiguous 38-to-37 annotation join')
        mappings[key] = old
    for r in records:
        row = dict(r)
        key = (str(r['CHROM']).removeprefix('chr'),int(r['POS']),r['REF'])
        old = mappings.get(key) if int(r['POS']) > 0 else None
        if old:
            row.update(CHROM=old[0],POS=str(old[1]))
        else:
            row['POS'] = '-1'
        yield row


def require_prepared_panel(prefix):
    """Reject an incomplete rebuild and validate lightweight preparation receipts.

    Current panels require their preparation receipt and
    preserve original rows and use POS=-1 for rejected coordinates.
    """
    prefix = Path(prefix).absolute()
    progress = Path(str(prefix) + '.prepare.inprogress')
    if progress.exists():
        raise ValueError('Ancient GRCh38 preparation is in progress or failed; finish 0.ancient_split.sh before PCA: ' + str(progress))
    receipt = Path(str(prefix) + '.prepare.json')
    paths = [Path(str(prefix)+s) for s in ('.pgen','.pvar','.psam')]
    if not all(p.is_file() for p in paths):
        raise ValueError('Ancient GRCh38 PGEN input is incomplete; run 0.ancient_split.sh and wait for completion: ' + str(prefix))
    if not receipt.is_file():
        raise ValueError('Missing current preparation receipt: ' + str(receipt))
    if receipt.is_file():
        info = json.loads(receipt.read_text())
        actual = [[str(p),p.stat().st_size,p.stat().st_mtime_ns] for p in paths]
        if not info.get('complete') or info.get('build') != 'GRCh38' or info.get('outputs') != actual:
            raise ValueError('Ancient preparation receipt does not match PGEN files; rerun 0.ancient_split.sh: ' + str(prefix))

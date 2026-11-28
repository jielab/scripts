"""Cache keys for published evidence tables. A split receipt binds every step to its inputs."""
from functools import lru_cache
import json
from pathlib import Path
import subprocess

from common import fingerprint, sha256
from metag_refine_tables import REFINE_FORMAT, SCHEMA

HERE = Path(__file__).resolve().parent

@lru_cache(maxsize=None)
def baseline_key():
    return fingerprint([sha256(HERE/'signals.py'), sha256(HERE/'food_sources.tsv'), SCHEMA,
                        sha256(HERE/'run.py'), sha256(HERE/'cache.py')])

def blast_info(prefix):
    return subprocess.check_output(['blastdbcmd', '-db', prefix, '-info'], text=True)

def assay_key(assay, opts):
    if assay == 'amr':
        annotation = opts.get('amr_annotation')
        return fingerprint([sha256(HERE/'resistome.py'), SCHEMA['refine_amr'],
                            sha256(opts['amr_markers']), sha256(annotation) if annotation else None,
                            subprocess.check_output(['diamond','version'],text=True).strip()])
    prefix = Path(opts['blast_db'])
    info = blast_info(str(prefix))
    # BLAST lists physical volumes even for an alias. Include their modification states.
    volumes = [Path(x.strip()) for x in info.partition('Volumes:')[2].splitlines() if x.strip()]
    if not volumes:
        volumes = [prefix]
    files = set()
    for vol in [prefix, *volumes]:
        loc = vol if vol.is_absolute() else prefix.parent / vol
        files.update(p for p in loc.parent.glob(loc.name+'.n*') if p.is_file())
    if not files:
        raise ValueError('Cannot locate local BLAST database volumes for '+str(prefix))
    inventory = prefix.parent/'references.json'
    return fingerprint([sha256(HERE/'blast_audit.py'), SCHEMA['refine_blast'], str(prefix), info,
                        [[str(p),p.stat().st_size,p.stat().st_mtime_ns] for p in sorted(files)],
                        sha256(inventory) if prefix.name=='competition' and inventory.is_file() else None,
                        opts['blast_pairs'], subprocess.check_output(['blastn','-version'],text=True).splitlines()[0]])

ASSAYS = {'blast':'blast_db','amr':'amr_markers'}

def requested(opts):
    return [a for a, key in ASSAYS.items() if opts.get(key)]

def valid_table(dest, old, name):
    p = dest/(name+'.tsv')
    return p.is_file() and old.get('table_hashes',{}).get(name) == sha256(p)

def complete(dest, old, sample, split_sha, opts, condition, keys):
    if (opts.get('replace') or old.get('format') != REFINE_FORMAT or old.get('state') != 'complete' or old.get('sample') != sample
        or old.get('split_sha256') != split_sha or old.get('baseline_key') != baseline_key()
        or (condition is not None and condition != [old.get('sample_type'),old.get('library_type')])):
        return False
    return (all(valid_table(dest, old, name) for name in SCHEMA)
            and all(old.get(a)=='complete' and old.get(a+'_key')==keys[a] for a in requested(opts)))

def receipt(dest):
    p = dest/'receipt.json'
    return json.loads(p.read_text()) if p.is_file() else {}

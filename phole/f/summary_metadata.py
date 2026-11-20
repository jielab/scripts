"""One machine-readable summary manifest per cohort, outside mplot/.
Each row describes one assay; provenance explanations belong in the README.
"""
import csv
from functools import lru_cache
from pathlib import Path

def metadata_path(folder):
    folder=Path(folder)
    return folder/(folder.name+'.sig.meta.tsv')

@lru_cache(maxsize=8)
def _read(path,mtime,size):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not reader.fieldnames or 'trait' not in reader.fieldnames:
            raise ValueError('Cohort summary metadata must have one row per assay and a trait column: '+path)
        key='ASSAY' if 'ASSAY' in reader.fieldnames else 'trait'
        rows={}
        for r in reader:
            if not r[key] or r[key] in rows:
                raise ValueError('Empty or duplicate trait in '+path)
            rows[r[key]]=r
        return rows

def read_metadata(folder):
    path=metadata_path(folder)
    if not path.is_file():return {}
    st=path.stat()
    return _read(str(path.resolve()),st.st_mtime_ns,st.st_size)

def for_summary(path):
    path=Path(path)
    return dict(read_metadata(path.parent.parent).get(path.name.removesuffix('.sig.txt'),{}))

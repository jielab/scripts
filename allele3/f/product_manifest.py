"""Checksummed permanent products; runtime directories are never archived."""
from pathlib import Path
import errno
import hashlib
import json
import os
import shutil
import tempfile


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''):h.update(block)
    return h.hexdigest()


def atomic_replace(source, destination):
    """Publish a file even when analysis and data are on different filesystems."""
    source,destination=Path(source),Path(destination)
    try:os.replace(source,destination)
    except OSError as exc:
        if exc.errno!=errno.EXDEV:raise
        fd,name=tempfile.mkstemp(prefix='.'+destination.name+'.',dir=destination.parent)
        os.close(fd)
        try:
            shutil.copy2(source,name)
            os.replace(name,destination)
            source.unlink()
        finally:
            Path(name).unlink(missing_ok=True)


def write_inventory(output, status='complete'):
    from product_paths import permanent
    output=Path(output).absolute()
    if not permanent(output):return
    items=[]
    for parent,dirs,names in os.walk(output,followlinks=False):
        dirs[:]=sorted(d for d in dirs if not (Path(parent)/d).is_symlink() and d not in {'logs','jobs','tmp','history','previous'} and not (Path(parent)==output and d=='qc'))
        for name in sorted(names):
            path=Path(parent)/name
            if path.is_symlink() or name in {'PRODUCT.json','SHA256SUMS'} or name.endswith(('.partial','.lock','.tmp')):continue
            stat=path.stat()
            items.append(dict(path=path.relative_to(output).as_posix(),bytes=stat.st_size,sha256=digest(path)))
            if (stat.st_size,stat.st_mtime_ns)!=(path.stat().st_size,path.stat().st_mtime_ns):
                raise ValueError('Product changed during inventory: '+str(path))
    record=dict(schema='allele3-permanent-product-1',status=status,genome_build='GRCh38',
        files=items,coordinate_convention='See table headers and scientific metadata; variant positions and BED intervals are not interchangeable.',
        runtime='Runtime logs/jobs/tmp are stored directly under analysis and are not part of this dataset')
    for name,content in [('PRODUCT.json',json.dumps(record,indent=2)+'\n'),
                         ('SHA256SUMS',''.join(r['sha256']+'  '+r['path']+'\n' for r in items))]:
        tmp=output/('.'+name+'.partial');tmp.write_text(content);tmp.replace(output/name)
    return record


def finish(output,result):
    write_inventory(output)
    return result

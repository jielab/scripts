#!/usr/bin/env python3
"""Export a portable current Shiny viewer containing aggregate evidence only."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile
from evidence_index import Indexer

ROOT = Path(__file__).resolve().parents[1]


def share_viewer(root, output, r_bin='Rscript', max_table_mb=128):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.suffix.lower() != '.zip' or output.is_relative_to(root):
        raise ValueError('Use a ZIP path outside the analysis result directory')
    with tempfile.TemporaryDirectory(prefix='le8-share-', dir='/tmp') as scratch:
        bundle = Path(scratch) / 'le8-share'
        shutil.copytree(ROOT / 'shiny', bundle / 'shiny')
        Indexer(root, bundle / 'evidence', max_table_mb).run()
        (bundle / 'README.txt').write_text(
            'LE8 aggregate evidence viewer\n\n'
            'Contains the aggregate database, catalogue and Shiny code. '
            'No participant arrays, model weights or raw input data are included.\n\n'
            "Install R packages: install.packages(c('shiny','bslib','DT','DBI','RSQLite','ggplot2','plotly'))\n"
            'From the extracted folder: Rscript shiny/app.R\n'
            'Open http://127.0.0.1:3839. Keep shiny/ and evidence/ together.\n', encoding='utf-8')
        env = {k: v for k, v in os.environ.items() if k not in {'LE8_EVIDENCE_DB', 'LE8_ANALYSIS_ROOT'}}
        subprocess.run([r_bin, str(bundle / 'shiny/app.R'), '--review'], env=env, check=True)
        output.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix='.le8-share-', suffix='.zip', dir=output.parent)
        os.close(fd)
        archive = Path(name)
        try:
            with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as stream:
                for path in sorted(bundle.rglob('*')):
                    if path.is_file():
                        stream.write(path, path.relative_to(bundle.parent))
            with zipfile.ZipFile(archive) as stream:
                if stream.testzip() is not None:
                    raise ValueError('ZIP verification failed')
            os.replace(archive, output)
        finally:
            archive.unlink(missing_ok=True)
    print(f'[LE8] Share package: {output}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--r-bin', default='Rscript')
    p.add_argument('--max-table-mb', type=int, default=128)
    a = p.parse_args()
    share_viewer(a.analysis_root, a.out, a.r_bin, a.max_table_mb)

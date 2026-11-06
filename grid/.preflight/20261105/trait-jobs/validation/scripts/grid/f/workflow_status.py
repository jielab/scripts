"""Report shared input readiness for the serial FINAL workflow."""
import argparse
from pathlib import Path


def readiness(root, traits, out_root=None, prepared=None, prs_cache=None, show_abm=True):
    cache = Path(prs_cache) if prs_cache else Path(root) / 'prsformer/cache'
    complete = all((cache / name).is_file() for name in
                   ('prepare.json', 'data.tsv.gz', 'variants.tsv.gz', 'genotypes.npy'))
    state = 'READY' if complete else 'WAIT'
    print(f'[GRID prepare] {state}: shared inputs for {traits}: {cache}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--analysis-root', required=True)
    p.add_argument('--traits', required=True)
    p.add_argument('--prs-cache')
    a = p.parse_args()
    readiness(a.analysis_root, a.traits, prs_cache=a.prs_cache)

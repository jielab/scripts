from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess, os

root = Path('/mnt/d/analysis/gu/threshold_review/20260906/legacy_compatible')
code = Path('/mnt/d/scripts/gu/f')
runtime = Path('/home/huangj/anaconda3/envs/gu/bin')
panel = '/mnt/d/data.BIG/refGen/1kg/37/samples.txt'

def run(d):
    env = dict(os.environ, PHYMLCPUS='1', PHYMLMPI='no', OMP_NUM_THREADS='1')
    with (d/'threshold_review.log').open('w') as log:
        def call(args):
            subprocess.run([str(a) for a in args], check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
        call([runtime/'python',code/'phyml_evidence.py','prepare','--out',d,'--sample-file',panel])
        print('PREPARED',d.name,flush=True)
        for phy in (d/'loci').glob('haplotypes.evidence.*.phy'):
            stats = Path(str(phy)+'_phyml_stats.txt')
            if stats.exists() and '3.3.20260528' not in stats.read_text():
                for suffix in ['_phyml_tree.txt', '_phyml_stats.txt', '_phyml_boot_trees.txt',
                               '_phyml_boot_stats.txt', '_phyml_tree.png']:
                    old = Path(str(phy)+suffix)
                    if old.exists():
                        old.rename(Path(str(old)+'.system-version'))
            if not Path(str(phy)+'_phyml_tree.txt').exists():
                print('TREE',d.name,phy.name,flush=True)
                call([runtime/'phyml','-i',phy,'-m','HKY85','-c','4','-a','e','-v','e','-b','100'])
        call([runtime/'python',code/'phyml_evidence.py','summarize','--out',d,'--sample-file',panel])
        call([runtime/'python',code/'phyml_risk.py','--out',d])
        call([runtime/'python',code/'phyml_layered_plot.py','--out',d,'--dpi','160'])
    print('COMPLETE',d.name,flush=True)

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(run,sorted(root.glob('chr*'))))

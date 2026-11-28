"""Post-scan HERV enrichment; pending prerequisites are never a successful job."""
import argparse
import json
from pathlib import Path
from herv_cohort import main, save_json

def run(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--gnomad',type=Path)
    p.add_argument('--allofus',type=Path)
    p.add_argument('--herv-input-dir',type=Path,required=True)
    a=p.parse_args(argv)
    if a.out.resolve()==a.herv_input_dir.resolve() or a.herv_input_dir.resolve() in a.out.resolve().parents:
        raise ValueError('Highlights output must be outside the HERV input directory')
    summary_path=a.herv_input_dir/'qc/summary.json'
    if not summary_path.is_file() or json.loads(summary_path.read_text()).get('status')!='complete':
        raise ValueError('Complete 2.herv_scan.sh before scan highlights; do not rerun successful scans.')
    stores={k:v for k,v in [('gnomad',a.gnomad),('allofus',a.allofus)] if v is not None}
    if not stores:raise ValueError('No summary stores supplied')
    # Run target stores explicitly, then both together for the comparison when
    # they share their parent (the public all-scan layout).
    if len(stores)==2:
        if a.gnomad.parent!=a.allofus.parent:raise ValueError('Summary stores must share their parent')
        args=['--target','gnomad,allofus','--summary-root',str(a.gnomad.parent)]
    else:
        target,store=next(iter(stores.items()));args=['--target',target,'--scan-dir',str(store)]
    result=main(['highlights','--out',str(a.out),'--herv-input-dir',str(a.herv_input_dir),*args])
    if result['status']!='complete' or result.get('comparison',{}).get('status')=='failed':
        raise RuntimeError('Highlights did not complete; see enrichment_status.json. Retry highlights only.')
    return result

if __name__=='__main__':run()

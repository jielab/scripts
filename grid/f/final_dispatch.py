#!/usr/bin/env python3
from __future__ import annotations
import argparse,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main(argv=None):
    a=argparse.ArgumentParser();a.add_argument("--analysis-root",default="/mnt/d/analysis/grid");a.add_argument("--out");a.add_argument("--traits",default="height,ldl,t2dm");a.add_argument("--max-table-mb",type=int,default=128);a.add_argument("--summary-only",action="store_true")
    x=a.parse_args(argv)
    subprocess.run([sys.executable,str(HERE/"final_compare.py"),"--analysis-root",x.analysis_root,"--traits",x.traits],check=True)
    cmd=[sys.executable,str(HERE/"browser_index.py"),"--analysis-root",x.analysis_root,"--traits",x.traits,"--max-table-mb",str(x.max_table_mb)]
    if x.out:cmd+=["--out",x.out]
    if x.summary_only:cmd+=["--summary-only"]
    subprocess.run(cmd,check=True)
    from cache_cleanup import clean_csx
    clean_csx(x.analysis_root, apply=True)
if __name__=="__main__":main()

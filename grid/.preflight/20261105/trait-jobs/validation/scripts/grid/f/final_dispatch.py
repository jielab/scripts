#!/usr/bin/env python3
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main(argv=None):
    a=argparse.ArgumentParser(description="Read-only current-result comparison, article tables and browser; no fitting.");a.add_argument("--analysis-root",default=configured_path('ANALYSIS_ROOT', 'grid'));a.add_argument("--out");a.add_argument("--traits",default="height,ldl,t2dm");a.add_argument("--max-table-mb",type=int,default=128);a.add_argument("--summary-only",action="store_true");a.add_argument("--dry-run",action="store_true")
    x=a.parse_args(argv)
    if x.dry_run:
        print("[GRID final] DRY RUN: final_compare.py -> paper_summary.py -> browser_index.py; read-only inputs, no training or cleanup")
        return
    from kinship import analysis_guard
    with analysis_guard(x.analysis_root):
        subprocess.run([sys.executable,str(HERE/"final_compare.py"),"--analysis-root",x.analysis_root,"--traits",x.traits],check=True)
        subprocess.run([sys.executable,str(HERE/"paper_summary.py"),"--analysis-root",x.analysis_root,"--traits",x.traits],check=True)
        cmd=[sys.executable,str(HERE/"browser_index.py"),"--analysis-root",x.analysis_root,"--traits",x.traits,"--max-table-mb",str(x.max_table_mb)]
        if x.out:cmd+=["--out",x.out]
        if x.summary_only:cmd+=["--summary-only"]
        subprocess.run(cmd,check=True)

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Exercise actual SAIGE GRM construction and FALSE -> TRUE prep resume."""
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
def run(*args): subprocess.run(list(map(str,args)),check=True)

with tempfile.TemporaryDirectory(prefix='saige-prep-test-') as tmp:
    gen=Path(tmp)
    run('/mnt/d/software/bin/plink2','--dummy','500','6000','--seed','823',
        '--memory','1000','--threads','2','--make-pgen','--out',gen/'chr1')
    # SAIGE requires explicit FID and IID in its FAM; PLINK1 conversion provides these.
    cmd=['bash',ROOT/'gwas.sh','prep_gwas','--typed-dir',gen,'--chromosomes','1','--threads','2','--memory','1000']
    run(*cmd,'--sparse-grm','FALSE')
    completed=list((gen/'.prep').glob('*.done.json'))+[gen/('typ.prune'+e) for e in ('.bed','.bim','.fam')]
    stamps={str(p):p.stat().st_mtime_ns for p in completed}
    run(*cmd,'--sparse-grm','TRUE')
    assert all(Path(p).stat().st_mtime_ns==t for p,t in stamps.items()), 'Prune/merge unexpectedly reran'
    assert (gen/'typ.sparseGRM.mtx').stat().st_size>0
    assert (gen/'typ.sparseGRM.mtx.sampleIDs.txt').stat().st_size>0
    grm_stamp=(gen/'.prep/sparseGRM.done.json').stat().st_mtime_ns
    run(*cmd,'--sparse-grm','TRUE')
    assert (gen/'.prep/sparseGRM.done.json').stat().st_mtime_ns==grm_stamp
    print('PASS: actual SAIGE sparse GRM; FALSE -> TRUE resumes prune/merge; repeated TRUE resumes GRM')

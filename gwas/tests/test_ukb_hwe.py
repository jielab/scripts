#!/usr/bin/env python3
"""Read-only large-sample HWE check on 100 typed variants; no GWAS is run."""
from pathlib import Path
import subprocess
import tempfile

src='/mnt/h/ukbGen/37/typ/chr1'
text=subprocess.check_output(['zstd','-dc',src+'.pvar.zst'],text=True)
variants=[line.split()[2] for line in text.splitlines() if line and not line.startswith('#')][:100]
assert len(variants)==100
log=Path(__file__).with_suffix('.log')
with tempfile.TemporaryDirectory(prefix='ukb-hwe-check-') as tmp:
    extract=Path(tmp)/'variants'; extract.write_text('\n'.join(variants)+'\n')
    cmd=['/mnt/d/software/bin/plink2','--pfile',src,'vzs','--extract',str(extract),
         '--memory','1000','--threads','2','--hwe','1e-15','0','--write-snplist','--out',tmp+'/check']
    with log.open('w') as h: subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT,check=True)
    assert Path(tmp+'/check.snplist').stat().st_size>0
print('PASS: --hwe 1e-15 0 accepted with all 488377 UKB typed samples (100-variant read-only check)')

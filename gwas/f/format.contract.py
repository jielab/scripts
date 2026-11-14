"""Per-stage GWAS dependencies and output receipts (no name-only completion)."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile


def state(value):
    path = Path(value).resolve()
    if not path.exists(): return {'path':str(path),'missing':True}
    if path.is_dir():
        return {'path':str(path),'files':[state(str(p)) for p in sorted(path.rglob('*')) if p.is_file()]}
    s=path.stat()
    result=dict(path=str(path),size=s.st_size,mtime=s.st_mtime_ns,ctime=s.st_ctime_ns)
    if s.st_size <= 8*1024*1024:
        result['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def effective_magma_n(env):
    """Resolve the exact worker N argument; no writes or environment mutation."""
    import re
    sidecar = Path(env['GWAS_DIR']) / (env['GWAS'] + '.magma.N')
    def fixed(value, source):
        if not re.fullmatch(r'[1-9][0-9]*', value):
            raise ValueError(f'Invalid MAGMA sample size from {source}: {value!r}')
        return dict(value=value, source=source, mode='N', argument='N=' + value)
    explicit = env.get('MAGMA_N', '').strip()
    if explicit:
        return fixed(explicit, 'MAGMA_N')
    if sidecar.is_file() and sidecar.stat().st_size:
        words = sidecar.read_text().split()
        if not words:
            raise ValueError('Empty MAGMA N sidecar: ' + str(sidecar))
        return fixed(words[0], 'sidecar')
    rows = env.get('GWAS_POST_MAGMA_ROWS', '')
    usable = 0
    if rows:
        with open(rows) as stream:
            for line in stream:
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 3 and re.fullmatch(r'[0-9]+(?:[.][0-9]+)?', fields[2]) and float(fields[2]) >= 50:
                    usable += 1
                    if usable > 1000:
                        return dict(value='', source='per_snp', mode='ncol', argument='ncol=N')
    return fixed(env.get('GWAS_N', '').strip(), 'GWAS_N')


def snapshot(stage, env):
    params={k:env.get(k,'') for k in {
        'format':'GRCH COL_MAP CHR_PREFIX HM3_MODE P_HM3 RSID_MODE RSID_UNMATCHED FILL_N FILL_EAF N_TOTAL DO_LIFTOVER',
        'magma':'GRCH MAGMA_N GWAS_N MAGMA_WINDOW CHRS',
        'lead':'GRCH P_LEAD LEAD_WINDOW CLUMP_KB CHRS REFGEN_POP',
    }[stage].split()}
    keys={
        'format':'RAW PHEF PERF_F DBSNP HM3 HM3_POS REFGEN_CLUMP REFGEN_KEEP CHAIN',
        'magma':'FINAL MAGMA_REF GENE_LOC SYNONYMS DBSNP',
        'lead':'FINAL REFGEN_CLUMP REFGEN_COJO REFGEN_KEEP',
    }[stage].split()
    paths=[]
    for key in keys:
        value=env.get(key,'')
        if not value: continue
        paths.append(value)
        if key in ('MAGMA_REF','REFGEN_CLUMP','REFGEN_COJO') and not Path(value).exists():
            paths.extend(str(p) for p in sorted(Path(value).parent.glob(Path(value).name+'*')) if p.is_file())
    if stage == 'magma':
        params['effective_n'] = effective_magma_n(env)
        paths.append(str(Path(env['GWAS_DIR']) / (env['GWAS'] + '.magma.N')))
    here=Path(__file__).resolve().parent
    paths.extend(str(here / name) for name in ('format.f.sh','format.py','format.contract.py','format.contract.sh'))
    for command in {'format':['bgzip','tabix'], 'magma':['magma'],'lead':['plink2','gcta']}[stage]:
        binary=shutil.which(command)
        if binary: paths.append(binary)
    return dict(schema=1,stage=stage,parameters=params,inputs=[state(x) for x in dict.fromkeys(paths)])


def outputs(stage,e):
    if stage=='format':
        return [e['FINAL'],e['FINAL']+'.tbi',e['FINAL']+'.csi']
    if stage=='magma': return [e['MAGMA_PREFIX']+x for x in ('.genes.out','.genes.raw')] + [e['MAGMA_DIR']+'/magma.done']
    return [e['MERGED']+x for x in ('.awk.snp','.clump.done','.cojo.done','.clumps','.jma.cojo','.ldr.cojo')]


def matches_receipt(saved, current):
    return saved == current


def main():
    if sys.argv[1:] == ['resolve-n']:
        print(effective_magma_n(os.environ)['argument'])
        return 0
    action,stage,receipt=sys.argv[1:]; e=os.environ
    current=snapshot(stage,e)
    if action == 'identity':
        print(hashlib.sha256(json.dumps(current, sort_keys=True).encode()).hexdigest())
        return 0
    value=dict(identity=current,outputs=[state(p) for p in outputs(stage,e)])
    path=Path(receipt)
    if action=='check':
        try: saved=json.loads(path.read_text())
        except (OSError,ValueError): return 1
        return 0 if matches_receipt(saved,value) else 1
    # format must have a data file and an index; empty lead selections still
    # require both independently written completion markers.
    required=outputs(stage,e)
    if stage=='format': required=[e['FINAL']]+[next((p for p in required[1:] if Path(p).is_file()),'')]
    elif stage=='lead': required=required[:3]
    if any(not p or not Path(p).is_file() or Path(p).stat().st_size==0 for p in required):
        raise RuntimeError('Cannot seal incomplete '+stage)
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(dir=path.parent,prefix='.contract-')
    with os.fdopen(fd,'w') as f: json.dump(value,f,sort_keys=True); f.write('\n')
    os.replace(tmp,path)
    return 0

if __name__=='__main__': raise SystemExit(main())

"""Shared inference identity and validation for all PRS-CSx output modes."""
import hashlib,json,math
from pathlib import Path
from pipeline_io import stamp

ROOT=Path(__file__).resolve().parent
POPS=['AFR','EAS','EUR','SAS']

def find_gwas(directory,trait,pop):
    tags=[f'{trait}.{pop}']+(['t2dm.AFA'] if trait=='t2dm' and pop=='AFR' else [])
    for tag in tags:
        for p in (Path(directory)/tag/'gwas'/f'{tag}.gz',Path(directory)/f'{tag}.gz'):
            if p.is_file() and p.stat().st_size:return p.resolve()
    raise FileNotFoundError(f'Missing GWAS: {trait}.{pop} in {directory}')

def validate_mcmc(phi,n,b,t,seed):
    if not(n>b>=0 and t>0 and n//t-b//t>=2 and seed>=0):
        raise ValueError('Require >=2 retained MCMC samples, positive thin, and a nonnegative seed')
    if phi!='auto' and not(math.isfinite(float(phi)) and float(phi)>0):raise ValueError('Invalid phi')

def inference_signature(phi,n,b,t,seed,override,chrs,gwas,snpinfo,bim,ld):
    validate_mcmc(phi,n,b,t,seed)
    code=[ROOT/p for p in ['csx/PRScsx.py','csx/parse_genet.py','csx/mcmc_gtb.py','csx/gigrnd.py','prepare_sumstats.py','split_sumstats.py','normalize_csx_weights.py','csx_config.py']]
    obj={'settings':[str(phi),n,b,t,seed,str(override),list(map(int,chrs))],
         'files':[stamp(p) for p in list(gwas)+[snpinfo,str(bim)+'.bim']+list(ld)],
         'code':[(p.name,hashlib.sha256(p.read_bytes()).hexdigest()) for p in code]}
    return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()

def sample_size(meta,override,pop):
    value=meta['n_gwas_median']
    if override:
        value=dict(x.split('=',1) for x in override.replace(';',',').replace(' ',',').split(',') if x).get(pop,value) if '=' in override else override
    if value is None or not math.isfinite(float(value)) or float(value)<2:raise ValueError(f'Missing/invalid GWAS N: {pop}')
    return round(float(value))

if __name__=='__main__':
    import sys
    phi,n,b,t,seed,override,chrs,snpinfo,bim,*paths=sys.argv[1:]
    print(inference_signature(phi,int(n),int(b),int(t),int(seed),override,chrs.split(),paths[:4],snpinfo,bim,paths[4:]))

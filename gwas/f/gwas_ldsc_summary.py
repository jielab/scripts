#!/usr/bin/env python3
"""Summarize only the requested run's LDSC logs, retaining unavailable traits/pairs."""
import argparse
import csv
import itertools
import math
from pathlib import Path
import re
import subprocess

RG_COLUMNS='p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se'.split()


def read_rg(path):
    text=Path(path).read_text()
    if 'Analysis finished at' not in text: raise ValueError('Incomplete LDSC log: '+str(path))
    marker='Summary of Genetic Correlation Results'
    if marker not in text: raise ValueError('Missing rg summary: '+str(path))
    lines=text.rsplit(marker,1)[1].strip().splitlines()
    header=lines[0].split()
    if header!=RG_COLUMNS: raise ValueError('Unexpected rg columns: '+str(path))
    rows=[]
    for line in lines[1:]:
        if not line.strip(): break
        values=line.split()
        if len(values)!=len(header): raise ValueError('Malformed rg summary: '+line)
        row=dict(zip(header,values))
        for key in ('p1','p2'):
            row[key]=Path(row[key]).name.removesuffix('.sumstats.gz')
        rows.append(row)
    return rows


def read_h2(path):
    text=Path(path).read_text()
    if 'Analysis finished at' not in text: raise ValueError('Incomplete LDSC log: '+str(path))
    fields={}
    for name,pattern in [('h2','Total Observed scale h2'),('intercept','Intercept'),('ratio','Ratio')]:
        match=re.search(re.escape(pattern)+r':\s*([-+\deE.]+)\s*\(([-+\deE.]+)\)',text)
        fields[name],fields[name+'_se']=match.groups() if match else ('NA','NA')
    for name,pattern in [('lambda_gc','Lambda GC'),('mean_chi2',r'Mean Chi\^2')]:
        match=re.search(pattern+r':\s*([-+\deE.]+)',text)
        fields[name]=match.group(1) if match else 'NA'
    if not finite(fields['h2']) or not finite(fields['h2_se']):
        raise ValueError('Missing/non-finite h2 estimate: '+str(path))
    return fields


def finite(value):
    try: return math.isfinite(float(value))
    except (ValueError,TypeError): return False


def write_tsv(path,columns,rows):
    tmp=path.with_suffix('.tmp')
    with tmp.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t',lineterminator='\n',restval='NA')
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def summarize(out,run_h2=True,run_rg=True):
    out=Path(out)
    with (out/'inputs.status.tsv').open() as handle: statuses=list(csv.DictReader(handle,delimiter='\t'))
    traits=[Path(row['FILE']).name.removesuffix('.gz') for row in statuses]
    available=[trait for trait,row in zip(traits,statuses) if row['STATUS']=='SCHEDULED']
    reasons={trait:row['REASON'] for trait,row in zip(traits,statuses) if row['STATUS']!='SCHEDULED'}
    h2=[]
    for trait in traits:
        row=dict(trait=trait,scope='1-22',scale='observed',status='UNAVAILABLE',reason=reasons.get(trait,'h2 disabled'))
        if run_h2 and trait in available:
            row.update(read_h2(out/'h2.log'/(trait+'.h2.log')))
            row.update(status='ESTIMATED',reason='Standard LDSC; X/Y/MT not estimated, see trait qc chromosome audit')
        h2.append(row)
    write_tsv(out/'h2.tsv',['trait','h2','h2_se','intercept','intercept_se','ratio','ratio_se',
                          'lambda_gc','mean_chi2','scope','scale','status','reason'],h2)
    pairs={}
    if run_rg:
        for trait in available[:-1]:
            for row in read_rg(out/'rg.log'/(trait+'.rg.log')):
                key=tuple(sorted((row['p1'],row['p2'])))
                if key in pairs: raise ValueError('Duplicate rg pair: '+str(key))
                pairs[key]=row
    rg=[]
    for a,b in itertools.combinations(traits,2):
        key=tuple(sorted((a,b)))
        row=dict(p1=a,p2=b,status='UNAVAILABLE',reason='; '.join(reasons[t] for t in (a,b) if t in reasons) or 'rg disabled')
        if run_rg and a in available and b in available:
            if key not in pairs: raise ValueError('Missing rg pair: '+str(key))
            row.update(pairs[key])
            valid=finite(row['rg']) and finite(row['p']) and finite(row['se'])
            row.update(status='ESTIMATED' if valid else 'UNAVAILABLE',
                       reason='Standard LDSC, chromosomes 1-22' if valid else 'LDSC returned a non-finite estimate; see rg log')
        rg.append(row)
    write_tsv(out/'rg.tsv',RG_COLUMNS+['status','reason'],rg)
    if run_rg:
        subprocess.run(['Rscript',str(Path(__file__).with_name('gwas_ldsc_plot.R')),str(out)],check=True)
    else:
        (out/'rg.png').unlink(missing_ok=True)
    print(f'Summarized {sum(r["status"]=="ESTIMATED" for r in h2)} h2 estimates and '
          f'{sum(r["status"]=="ESTIMATED" for r in rg)} rg pairs: {out}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--h2',choices=('True','False'),default='True')
    parser.add_argument('--rg',choices=('True','False'),default='True')
    args=parser.parse_args()
    summarize(args.output_dir,args.h2=='True',args.rg=='True')

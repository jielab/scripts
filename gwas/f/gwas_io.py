#!/usr/bin/env python3
"""Streaming, explicit effect-allele normalization; preserve chromosome ordering."""
import gzip
import math
import os
from pathlib import Path

HEADER='SNP CHR POS EA NEA EAF N BETA SE P'.split()

def reader(path):
    return gzip.open(path,'rt') if str(path).endswith(('.gz','.bgz')) else open(path)

def normalize(files,method,out):
    count=0; rejected=0; tmp=str(out)+'.tmp'
    with gzip.open(tmp,'wt') as target:
        target.write('\t'.join(HEADER)+'\n')
        for file in files:
            with reader(file) as source:
                header=source.readline().strip().lstrip('#').split()
                for line in source:
                    vals=line.split()
                    if len(vals)!=len(header): raise ValueError(f'Malformed row in {file}')
                    d=dict(zip(header,vals))
                    if d.get('TEST','ADD')!='ADD': continue
                    if method=='plink2':
                        beta=d.get('BETA')
                        if beta is None and 'OR' in d:
                            try: beta=str(math.log(float(d['OR'])))
                            except (ValueError,OverflowError): beta='NA'
                        row=[d['ID'],d['CHROM'],d['POS'],d['A1'],d.get('OMITTED',d.get('AX','NA')),d['A1_FREQ'],d['OBS_CT'],beta,d.get('SE',d.get('LOG(OR)_SE','NA')),d['P']]
                    elif method=='regenie':
                        # REGENIE ALLELE1 is effect allele; ALLELE0 is reference.
                        try: pv=str(10**(-min(float(d['LOG10P']),300)))
                        except ValueError: pv='NA'
                        row=[d['ID'],d['CHROM'],d['GENPOS'],d['ALLELE1'],d['ALLELE0'],d['A1FREQ'],d['N'],d['BETA'],d['SE'],pv]
                    elif method=='saige':
                        n=d.get('N')
                        if n is None:
                            try: n=str(int(d['N_case'])+int(d['N_ctrl']))
                            except (KeyError,ValueError): n='NA'
                        row=[d.get('MarkerID',d.get('SNPID')),d['CHR'],d['POS'],d['Allele2'],d['Allele1'],d['AF_Allele2'],n,d['BETA'],d['SE'],d['p.value']]
                    else: raise ValueError(method)
                    try:
                        nums=[float(row[j]) for j in (2,5,6,7,8,9)]
                        valid=all(math.isfinite(x) for x in nums) and 0<=nums[1]<=1 and nums[2]>0 and nums[4]>0 and 0<=nums[5]<=1
                        valid=valid and row[3]!=row[4] and all(x not in (None,'NA','.') for x in row[:5])
                    except (ValueError,TypeError): valid=False
                    if not valid: rejected+=1; continue
                    row[9]=str(max(float(row[9]),1e-300))
                    target.write('\t'.join(row)+'\n'); count+=1
    if not count:
        os.unlink(tmp); raise ValueError(f'No valid additive results for {out}')
    os.replace(tmp,out)
    Path(str(out)+'.qc.txt').write_text(f'written\t{count}\nrejected\t{rejected}\n')

if __name__=='__main__':
    import sys
    normalize(sys.argv[3:],sys.argv[1],sys.argv[2])

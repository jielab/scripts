from pathlib import Path
import sys
import numpy as np, pandas as pd, h5py
root=Path(sys.argv[1]);root.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(2026);n=120; samples=80
snps=[f'rs{i+1}' for i in range(n)];pops=['AFR','EAS','EUR','SAS']
ref=root/'ref';ref.mkdir(exist_ok=True)
snp=pd.DataFrame(dict(CHR=[22]*n,SNP=snps,BP=np.arange(n)+100001,A1=['A']*n,A2=['C']*n))
for p in ['AFR','AMR','EAS','EUR','SAS']:snp['FRQ_'+p]=.3
for p in ['AFR','AMR','EAS','EUR','SAS']:snp['FLP_'+p]=1
snp.to_csv(ref/'snpinfo_mult_1kg_hm3',sep='\t',index=False)
for p in pops:
 d=ref/('ldblk_1kg_'+p);d.mkdir(exist_ok=True)
 with h5py.File(d/'ldblk_1kg_chr22.hdf5','w') as h:
  b=h.create_group('blk_1');b['ldblk']=np.eye(n);b['snplist']=np.array(snps,dtype='S')
 g=root/'gwas'/('height.'+p)/'gwas';g.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(dict(SNP=snps,CHR=22,POS=np.arange(n)+100001,EA='A',NEA='C',EAF=.3,N=10000,BETA=rng.normal(0,.02,n),SE=.01,P=.5)).to_csv(g/f'height.{p}.gz',sep='\t',index=False,compression='gzip')
gen=root/'gen';gen.mkdir(exist_ok=True)
with (gen/'input.vcf').open('w') as f:
 f.write('##fileformat=VCFv4.2\n##contig=<ID=22>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
 f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(f's{i:03}' for i in range(samples))+'\n')
 for j in range(n):
  g=rng.binomial(2,.3,samples);gt=np.array(['0/0','0/1','1/1'])[g]
  f.write(f'22\t{100001+j}\t{snps[j]}\tC\tA\t.\tPASS\t.\tGT\t'+'\t'.join(gt)+'\n')
# deliberately scramble IDs to exercise upstream ID/score alignment
order=rng.permutation(samples)
pc=rng.normal(size=(samples,10))
d=pd.DataFrame(pc,columns=[f'PC{i}' for i in range(1,11)]);d.insert(0,'IID',[f's{i:03}' for i in order]);d.to_csv(root/'pca.tsv.gz',sep='\t',index=False)
m=pd.DataFrame(rng.normal(0,3,(4,10)),columns=d.columns[1:]);m.insert(0,'POP',pops);m.to_csv(root/'centers.tsv',sep='\t',index=False)
print(root)


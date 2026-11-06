import json,os,pathlib,importlib.util,numpy as np,pandas as pd
from types import SimpleNamespace
from shared_cohort import preparation_options
import prsformer
root=pathlib.Path('/work/sph-huangj/analysis/grid')
audit=pathlib.Path('/work/sph-huangj/scripts/grid/.preflight/20261105')
a,block=preparation_options(SimpleNamespace(analysis_root=str(root)))
paths={k:getattr(a,k) for k in ['dir_gen','pheno_file','ancestry_file','snp_list','remove','cache_dir']}
paths.update({k:block[k] for k in ['score_root','csx_analysis_root','output_root','cpu_abm_prepared','dosage_supplement_dir','ld_block_dir']})
print('RESOLVED_PATHS',json.dumps(paths,indent=2))
assert not a.remove or a.remove=='none'
assert paths['cache_dir']==str(root/'prsformer/cache')
assert paths['score_root']==os.environ['GRID_SCORE_DIR']
assert paths['dir_gen']==os.environ['GRID_TARGET_DIR']
assert paths['ancestry_file']==str(pathlib.Path(os.environ['GRID_PCA_FILE']).parent/'ukb.ancestry.auto.tsv.gz')
spec=importlib.util.spec_from_file_location('audit_io','f/0.common.py');io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
ids=set(prsformer.sample_ids(pathlib.Path(a.dir_gen)/'chr1.psam','pfile'))
ids={x for x in ids if not str(x).startswith('-')}
pca=pd.read_csv(os.environ['GRID_PCA_FILE'],sep='\t',usecols=['IID'],dtype=str)
pc_ids=set(pca.IID)
imp_ids={x for x in prsformer.sample_ids(pathlib.Path(os.environ['GRID_IMP_DIR'])/'chr1.psam','pfile') if not str(x).startswith('-')}
print('TARGET_VALID_IDS',len(ids),'PCA_IDS',len(pc_ids),'TARGET_WITHOUT_PCA',len(ids-pc_ids),'IMP_IDS',len(imp_ids),'TARGET_WITHOUT_IMP',len(ids-imp_ids))
rows=[]
for trait in ['height','ldl','t2dm']:
 p=pathlib.Path(paths['score_root'])/trait/'1csx.scores.rds'
 d=io.read_result_table(p,dtype={'eid':str,'IID':str})
 idcol='eid' if 'eid' in d else 'IID'
 scoreids=set(d[idcol].astype(str))
 cols=[next((c for c in [f'CSX_{pop}',f'csx.{pop}'] if c in d),None) for pop in ['AFR','EAS','EUR','SAS']]
 assert all(cols),(str(p),d.columns.tolist())
 assert not d[idcol].duplicated().any()
 assert np.isfinite(d[cols].to_numpy(float)).all()
 record=dict(trait=trait,rows=len(d),columns=d.columns.tolist(),missing_target_ids=len(ids-scoreids),extra_score_ids=len(scoreids-ids),scores_without_pca=len(scoreids-pc_ids))
 print('SCORE_CONTRACT',json.dumps(record)); rows.append(record)
 record['roster_matches_new_scores']=ids==scoreids
 record['scores_without_new_projection_input']=len(scoreids-imp_ids)
(audit/'data-contracts.json').write_text(json.dumps(dict(paths=paths,scores=rows),indent=2)+'\n')
print('DATA CONTRACTS INSPECTED')

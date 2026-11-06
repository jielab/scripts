import os,json,pathlib,pandas as pd,numpy as np
from types import SimpleNamespace
from shared_cohort import preparation_options,prepare_frozen_roster
import prsformer
root=pathlib.Path('/work/sph-huangj/scripts/grid/.preflight/20261105/roster-only')
a,block=preparation_options(SimpleNamespace(analysis_root=str(root)))
files=prsformer.genotype_files(a.dir_gen,prsformer.chromosomes(a.chrs))
d=prsformer.make_cohort(a,files)
pca=pd.read_csv(os.environ['GRID_PCA_FILE'],sep='\t',usecols=['IID'],dtype=str)
missing=set(d.eid)-set(pca.IID)
imp_ids=set(prsformer.sample_ids(pathlib.Path(os.environ['GRID_IMP_DIR'])/'chr1.psam','pfile'))
report=dict(cohort=len(d),cohort_without_current_pca=len(missing),cohort_without_imp_genotypes=len(set(d.eid)-imp_ids))
print('COHORT_ALIGNMENT',json.dumps(report),flush=True)
(root.parent/'cohort-alignment.json').write_text(json.dumps(report,indent=2)+'\n')
prepare_frozen_roster(a,root)
print('FROZEN_FAMILY_ROSTER PASS',a.keep,a.split_file,a.split_group_file,flush=True)

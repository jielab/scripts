"""Actual block production and CPU training for both independent extensions."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'f'))
sys.path.insert(0,str(Path(__file__).parent))
import block_inputs as b
from csx_blocks import main
from test_block_producer import fixture


def test_shared_scalar_fusion_is_fixed_external_n_and_allele_aware(tmp_path):
    path,cfg=fixture(tmp_path,'ALL.37.bed',n=40)
    # Flip the effect allele of one population, retaining correct dosage offsets.
    weight=tmp_path/'weights/AFR.tsv'
    w=pd.read_csv(weight,sep='\t')
    w[['A1','A2']]=w[['A2','A1']].to_numpy()
    w.to_csv(weight,sep='\t',index=False)
    prov=tmp_path/'weights/height/provenance.json'
    value=json.loads(prov.read_text());value['models']['csx.AFR']['weights_sha256']=b.sha256(weight)
    prov.write_text(json.dumps(value))
    geno=np.load(tmp_path/'genotype/genotypes.npy')
    import rdata
    scores=rdata.read_rds(tmp_path/'scores/height/1csx.scores.rds')
    scores['csx.AFR']=(2-geno)@np.array([.125,.25])
    from data import write_rds
    write_rds(scores,tmp_path/'scores/height/1csx.scores.rds')
    _,_,receipt=b.prepare_trait(path,'height',device='cpu')
    arr=np.load(receipt['score_file'])
    assert arr.shape==(40,2,1)
    fusion=receipt['shared_weight_fusion']
    assert fusion['target_outcomes_used'] is False and fusion['official_prscsx_meta'] is False
    assert fusion['alpha']==dict(AFR=.1,EAS=.2,EUR=.3,SAS=.4)
    expected=sum(fusion['alpha'][p]*scores['csx.'+p].to_numpy() for p in b.POPS)
    np.testing.assert_allclose(arr.sum((1,2)),expected,rtol=1e-6)
    # Altering N is a scientific feature change and invalidates a completed cache.
    metadata=Path(str(weight)+'.metadata.json')
    metadata.write_text(json.dumps(dict(n_gwas_used=90000)))
    from abm_core.grid_inputs import locate_blocks
    with pytest.raises(ValueError,match='fusion provenance changed'):
        locate_blocks(dict(block_config=str(path),execution={'device':'cpu'}),'height',allow_build=False)


@pytest.mark.parametrize('module',['blockridge','blocknet'])
def test_real_tiny_cpu_block_fit_owns_untouched_test_performance(tmp_path,module):
    path,cfg=fixture(tmp_path,'ALL.37.bed',n=240)
    genotype=np.load(tmp_path/'genotype/genotypes.npy')
    source=tmp_path/'genotype/data.tsv.gz'
    d=pd.read_csv(source,sep='\t',dtype=str)
    rng=np.random.default_rng(114)
    d['split']=np.where(np.arange(len(d))<180,'train','test')
    d['age']=rng.normal(size=len(d));d['sex']=rng.integers(0,2,len(d))
    d['PC1']=rng.normal(size=len(d));d['PC2']=rng.normal(size=len(d))
    # Non-collinear genotype signal is representable in the two shared tokens.
    d['height']=100+4*genotype[:,0]-2*genotype[:,1]+.2*d.age+rng.normal(0,.05,len(d))
    d.to_csv(source,sep='\t',index=False)
    cfg.update(d_model=8,heads=2,latent_tokens=2,layers=1,dropout=0.,epochs=2,batch_size=16,patience=2)
    path.write_text(json.dumps(cfg))
    output=tmp_path/'analysis'
    main([module,'--analysis-root',str(output),'--trait','height','--block-config',str(path),
        '--cohort',str(source),'--validation-n','120','--device','cpu','--threads','1'])
    found=list((output/'csx_blocks').glob('*/height/'+module+'/COMPLETE.json'))
    assert len(found)==1
    receipt=json.loads(found[0].read_text())
    perf=pd.read_csv(found[0].parent/'performance.tsv',sep='\t').iloc[0]
    assert perf.N_total==60 and perf.build_N+perf.tune_N==120
    assert np.isfinite(perf.mse) and perf.test_used_for_fitting==False
    assert perf.official_prscsx==False
    p=pd.read_csv(found[0].parent/'_LOCAL_ONLY/predictions.tsv.gz',sep='\t',dtype={'id':str})
    assert set(p.loc[p.block_role.eq('test'),'id'])==set(d.loc[d.split.eq('test'),'eid'])
    if module=='blockridge':
        assert perf.prediction_r2>.98
    else:
        history=pd.read_csv(found[0].parent/'training.tsv',sep='\t')
        assert len(history)==2 and history.train_loss.notna().all()

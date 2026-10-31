import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import block_inputs
import regional_dimensionality as rd


def bed(n):
    return pd.DataFrame(dict(block_id=[f'b{i}' for i in range(n)],CHR=1,start=np.arange(n)*10,end=np.arange(1,n+1)*10,block_local=np.arange(n)))


def weights(betas, positions=None):
    return pd.DataFrame(dict(BETA=betas,CHR=1,BP=positions if positions is not None else np.arange(len(betas))*10+1))


def test_uniform_and_concentrated():
    m=rd.concentration(np.ones(10))
    assert m['hhi']==pytest.approx(.1)
    assert m['effective_blocks']==pytest.approx(10)
    assert m['effective_blocks_fraction']==pytest.approx(1)
    assert m['gini']==pytest.approx(0)
    m=rd.concentration([1]+[0]*9)
    assert m['hhi']==1 and m['effective_blocks']==1
    assert m['effective_blocks_fraction']==.1
    assert m['gini']==pytest.approx(.9)
    assert m['top10pct_share']==1


@pytest.mark.parametrize('definition,mass',[('beta_squared',8),('abs_beta',4)])
def test_sign_empty_and_unassigned(definition,mass):
    m,top=rd.population_metrics(weights([2,-2,4],[1,11,999]),bed(10),definition)
    assert m['n_blocks']==10
    assert m['total_assigned_mass']==mass
    assert m['assigned_snp_fraction']==pytest.approx(2/3)
    assert m['assigned_mass_fraction']==pytest.approx(1/3 if definition=='beta_squared' else .5)
    assert m['effective_blocks_fraction']==pytest.approx(.2)
    assert len(top)==10 and (top.block_mass==0).sum()==8


def test_maps_and_summary():
    rows=[]
    for pop,n in zip(rd.POPS,[2,4,8,10]):
        m,_=rd.population_metrics(weights([1]),bed(n),'beta_squared')
        assert m['effective_blocks_fraction']==pytest.approx(1/n)
        rows.append(dict(trait='height',population=pop,mass_definition='beta_squared',**m))
    s=rd.summarize(pd.DataFrame(rows)).iloc[0]
    assert s.populations_available==4
    for name,q in [('median',.5),('q25',.25),('q75',.75)]:
        assert s[name+'_effective_blocks_fraction']==pytest.approx(np.quantile([.5,.25,.125,.1],q))


def test_gain():
    rows=[]
    for t in ['height','ldl','t2dm']:
        for definition in rd.MASS_DEFINITIONS:
            m,_=rd.population_metrics(weights([1]),bed(10),definition)
            rows.append(dict(trait=t,population='EUR',mass_definition=definition,**m))
    summary=rd.summarize(pd.DataFrame(rows))
    perf=pd.DataFrame([dict(trait=t,method=method,mse=m,brier=b,stage='test') for t in ['height','ldl','t2dm'] for method,m,b in [('csx-linear',10,.2),('csx-blocknet',8,.1)]])
    g=rd.join_gain(summary,perf).set_index('trait')
    assert g.loc['height','relative_primary_loss_reduction']==pytest.approx(.2)
    assert g.loc['ldl','primary_loss']=='mse'
    assert g.loc['t2dm','relative_primary_loss_reduction']==pytest.approx(.5)
    assert rd.join_gain(summary,pd.DataFrame()).interpretation_status.eq('not_available').all()
    with pytest.raises(ValueError):rd.join_gain(summary,pd.concat([perf,perf]))


@pytest.mark.parametrize('shared',[False,True])
def test_cli_only_weights_blocks_no_hash_no_current(tmp_path,monkeypatch,shared):
    def forbidden(*args,**kwargs):raise AssertionError('Full-file hashing forbidden')
    monkeypatch.setattr(block_inputs,'sha256',forbidden)
    names=['ALL.37.bed'] if shared else [p+'.37.bed' for p in rd.POPS]
    for i,name in enumerate(names):
        (tmp_path/name).write_text(''.join(f'1\t{j*10}\t{(j+1)*10}\n' for j in range(i+2)))
    csx=tmp_path/'csx/height';csx.mkdir(parents=True)
    models={}
    for pop in rd.POPS:
        w=tmp_path/(pop+'.weights.tsv')
        w.write_text('SNP A1 A2 BETA CHR BP\nrs1 A G 2 1 1\nrs2 C T -2 1 11\nrs3 A T 4 1 999\n')
        models['csx.'+pop]=dict(weights_file=str(w),weights_sha256='must_not_scan')
    (csx/'provenance.json').write_text(json.dumps(dict(models=models)))
    config=tmp_path/'config.json'
    config.write_text(json.dumps(dict(genome_build=37,populations=list(rd.POPS),ld_block_dir=str(tmp_path),ld_block_files=','.join(names),csx_analysis_root=str(tmp_path/'csx'))))
    current=tmp_path/'abm/FINAL/height/blocknet.CURRENT.json'
    current.parent.mkdir(parents=True);current.write_text('{"run":"frozen"}')
    before=rd.identity(current)
    rd.main(['--analysis-root',str(tmp_path),'--config',str(config),'--traits','height'])
    assert rd.identity(current)==before
    tables=rd.completed_tables(tmp_path)
    assert len(tables[0])==8 and tables[3].interpretation_status.eq('not_available').all()
    assert len(rd.paper_table(tmp_path,['height']))==1
    assert (tmp_path/'architecture/regional_dimensionality.pdf').is_file()
    w.write_text(w.read_text()+'rs4 A T 1 1 20\n')
    with pytest.raises(ValueError,match='changed'):rd.completed_tables(tmp_path)


def test_zero_mass():
    m,_=rd.population_metrics(weights([0]),bed(10),'beta_squared')
    assert m['status']=='zero_mass' and np.isnan(m['effective_blocks'])

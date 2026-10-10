"""C5 cellulation: annotation, native CIGMA and published evidence ingestion."""
import gzip
import importlib.util
import itertools
import json
import numpy as np
import os
import pandas as pd
import pytest
import subprocess
import sys
import tempfile
from pathlib import Path
from scipy import sparse
from scipy.stats import fisher_exact
from unittest import TestCase, mock

ROOT = Path(__file__).resolve().parents[1]

def load_module(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / 'f' / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

SCRIPT = ROOT / 'f/c5.cellulation.py'
c5 = load_module('le8_validation_c5', 'c5.cellulation.py')
module = load_module('le8_validation_published', 'c5.import_published.py')
shared = load_module('le8_validation_dispatch', '0.common.py')

def csv(path, data):
    pd.DataFrame(data).to_csv(path,index=False)
    return path

@pytest.fixture
def inputs(tmp_path):
    u=csv(tmp_path/'universe.csv',dict(assay=['a','b','c','d','e','nt','bn','unknown'],
        gene=['A','B','C','D','E','NTPROBNP','NPPB',''],match_stratum=['x']*8))
    a=csv(tmp_path/'atlas.csv',dict(gene=['A','B','C','NPPB'],cell_type=['T','T','B','Heart'],
        source=['SYNTHETIC TEST ATLAS']*4,source_version=['test']*4))
    p=csv(tmp_path/'panels.csv',dict(model=['YS','YS','NS','NS'],feature=['a','b','a','d']))
    return u,a,p

def native_inputs(root,n=36,c=2):
    root.mkdir(exist_ok=True)
    rng=np.random.default_rng(12); ds=[f'{i:03d}' for i in range(n)]; cs=[f'CT{i}' for i in range(c)]
    def mat(name,x,cols=cs):
        z=pd.DataFrame(x,index=ds,columns=cols); z.index.name='donor'; z.to_csv(root/name)
    mat('y.csv',rng.normal(size=(n,c)))
    mat('nu.csv',np.full((n,c),.1)); mat('P.csv',np.full((n,c),1/c))
    z=rng.normal(size=(n,70)); K=z@z.T/70
    mat('K.csv',K,ds)
    row=dict(gene='A',tissue='SYNTHETIC',cohort='test',build='38',kinship_scope='cis',
        ctp='y.csv',ctnu='nu.csv',P='P.csv',K='K.csv',ctnu_definition='variance_of_pseudobulk_mean')
    return row, ds, cs

def fake_he(**args):
    """Returns preselected numbers to verify adapter contracts; no scientific fitting."""
    assert args['jk'] is True and args['verbose'] is True
    C=args['Y'].shape[1]
    est=dict(hom_g2=.2,V=np.diag(np.linspace(-.02,.12,C)),hom_e2=.2,W=np.eye(C)*.3)
    p=dict(V=.001,hom_g2=.01,vc=np.full(C,.02),var_hom_g2=.002,var_V=np.eye(C)*.001)
    if 'Kt' in args:
        est.update(hom_g2_b=.1,V_b=np.eye(C)*.1)
        p.update(V_b=.03,hom_g2_b=(.04,),vc_b=(np.full(C,.05),),var_hom_g2_b=.001,var_V_b=np.eye(C)*.002)
    return est,p

def test_bh_missing_family():
    q=c5.bh([.01,np.nan,.04]); np.testing.assert_allclose(q[[0,2]],[.03,.06]); assert np.isnan(q[1])
    assert c5.bh([.01],10)[0]==.1

@pytest.mark.parametrize('p',[-.1,1.1])
def test_bad_p(p):
    with pytest.raises(ValueError): c5.bh([p])

def test_alias_budget_not_collapsed(tmp_path, inputs):
    u,a,_=inputs; p=csv(tmp_path/'p.csv',dict(model=['M','M'],feature=['nt','bn']))
    z=c5.annotate(u,a,p,tmp_path/'out',plots=False)
    assert z['coverage'].iloc[0]['assays']==2
    assert z['coverage'].iloc[0]['unique_genes']==1
    assert z['coverage'].iloc[0]['background_genes']==6

def test_unknown_retained(inputs,tmp_path):
    z=c5.annotate(*inputs,tmp_path/'out',plots=False)
    r=z['annotation']; assert r.loc[r.assay=='unknown','annotation_status'].item()=='gene_unmapped'
    assert r.loc[r.assay=='e','annotation_status'].item()=='no_label_in_this_atlas'

def test_full_fisher_background(inputs,tmp_path):
    z=c5.annotate(*inputs,tmp_path/'out',plots=False)
    r=z['enrichment'].query("model=='YS' and cell_type=='T'").iloc[0]
    assert r['p']==pytest.approx(fisher_exact([[2,0],[0,4]],alternative='greater').pvalue)
    assert r['background_genes']==6

def test_duplicate_assay_rejected(inputs,tmp_path):
    u,a,p=inputs; d=pd.read_csv(u); pd.concat([d,d.iloc[:1]]).to_csv(u,index=False)
    with pytest.raises(ValueError): c5.annotate(u,a,p,tmp_path/'out',plots=False)

def test_mixed_atlas_versions_rejected(inputs,tmp_path):
    u,a,p=inputs; d=pd.read_csv(a); d.loc[0,'source_version']='other'; d.to_csv(a,index=False)
    with pytest.raises(ValueError): c5.annotate(u,a,p,tmp_path/'out',plots=False)

def test_ambiguous_gene_rejected():
    with pytest.raises(ValueError): c5.gene_name('A;B')
    assert c5.gene_name('ENSG000000001.3')=='ENSG000000001'

def test_fold_separation(inputs,tmp_path):
    u,a,_=inputs
    p=csv(tmp_path/'fold.csv',dict(model=['M','M','M'],feature=['a','a','b'],fold=['1','2','2']))
    z=c5.annotate(u,a,p,tmp_path/'out',plots=False)
    assert set(z['coverage'].model)=={'M|fold=1','M|fold=2'}
    stability=pd.read_csv(tmp_path/'out/cell.fold_coverage.csv')
    assert stability.loc[stability.feature=='a','frequency'].item()==1
    assert stability.loc[stability.feature=='b','frequency'].item()==.5

def test_exclusive_gene_randomization_exact():
    A={'a','b','c'}; B={'a','d','e'}; labels={'a','b','c'}
    z=c5.conditional_panel_test(A,B,labels)
    vals=[]; shared=A&B; pool=sorted(A^B)
    for pick in itertools.combinations(pool,len(A-B)):
        X=shared|set(pick); Y=shared|(set(pool)-set(pick))
        vals.append(len(X&labels)/len(X)-len(Y&labels)/len(Y))
    expected=np.mean(np.abs(np.array(vals)-np.mean(vals))>=abs(z['delta']-np.mean(vals))-1e-12)
    assert z['p']==pytest.approx(expected)

def test_identical_panels_not_superiority():
    z=c5.conditional_panel_test({'a'},{'a'},{'a'})
    assert z['delta']==0 and z['p']==1

def test_matched_null_reproducible(inputs,tmp_path):
    c5.annotate(*inputs,tmp_path/'a',plots=False,matched_draws=101,seed=17)
    c5.annotate(*inputs,tmp_path/'b',plots=False,matched_draws=101,seed=17)
    assert (tmp_path/'a/cell.matched_enrichment.csv').read_bytes()==(tmp_path/'b/cell.matched_enrichment.csv').read_bytes()

def test_panel_contrasts_equal_budget(inputs,tmp_path):
    u,a,p=inputs; ct=csv(tmp_path/'ct.csv',dict(model=['YS'],reference=['NS']))
    z=c5.annotate(u,a,p,tmp_path/'o',plots=False,contrasts=ct)
    assert len(z['contrasts'])==4
    d=pd.read_csv(p); d=d.iloc[:-1]; d.to_csv(p,index=False)
    with pytest.raises(ValueError): c5.annotate(u,a,p,tmp_path/'bad',plots=False,contrasts=ct)

def test_donor_identifiers_and_reorder(tmp_path):
    row,ds,cs=native_inputs(tmp_path)
    d=pd.read_csv(tmp_path/'nu.csv',dtype=str).iloc[::-1]; d.to_csv(tmp_path/'nu.csv',index=False)
    args,cells,paths=c5.read_inputs(row,tmp_path)
    assert cells==cs and args['Y'].shape==(36,2)
    assert c5.matrix(tmp_path/'y.csv').index[0]=='000'

def test_duplicate_matrix_header(tmp_path):
    f=tmp_path/'x.csv'; f.write_text('donor,A,A\n001,1,2\n')
    with pytest.raises(ValueError): c5.matrix(f)

@pytest.mark.parametrize('kind',['negative_noise','identity','bad_prop','fixed_intercept','donor_mismatch','bad_noise_definition','same_Kt'])
def test_native_input_rejections(tmp_path,kind):
    row,ds,cs=native_inputs(tmp_path)
    if kind=='negative_noise':
        d=pd.read_csv(tmp_path/'nu.csv',dtype={'donor':str}); d.iloc[0,1]=-.1; d.to_csv(tmp_path/'nu.csv',index=False)
    elif kind=='identity':
        d=pd.DataFrame(np.eye(len(ds)),index=ds,columns=ds); d.index.name='donor'; d.to_csv(tmp_path/'K.csv')
    elif kind=='bad_prop':
        d=pd.read_csv(tmp_path/'P.csv',dtype={'donor':str}); d.iloc[0,1]=1; d.to_csv(tmp_path/'P.csv',index=False)
    elif kind=='fixed_intercept':
        csv(tmp_path/'cov.csv',dict(donor=ds,intercept=np.ones(len(ds)))); row['fixed_covars']='cov.csv'
    elif kind=='donor_mismatch':
        d=pd.read_csv(tmp_path/'nu.csv',dtype=str).iloc[1:]; d.to_csv(tmp_path/'nu.csv',index=False)
    elif kind=='bad_noise_definition': row['ctnu_definition']='SD'
    else: row['Kt']='K.csv'
    with pytest.raises(ValueError): c5.read_inputs(row,tmp_path)

def test_parameter_count_guard(tmp_path):
    row,_,_=native_inputs(tmp_path,n=20,c=6)
    with pytest.raises(ValueError,match='parameter count'): c5.read_inputs(row,tmp_path)

def test_native_contract_keeps_negative_HE_and_uncertainty(tmp_path):
    row,ds,cs=native_inputs(tmp_path); manifest=csv(tmp_path/'manifest.csv',[row])
    out=tmp_path/'out'; rc=c5.run(manifest,out,_fit_function=fake_he)
    assert rc==0
    d=pd.read_csv(out/'cigma.results.csv'); c=pd.read_csv(out/'cigma.cell_types.csv')
    assert pd.isna(d.specificity.iloc[0]) and not d.variance_admissible.iloc[0]
    assert c.specific_variance.iloc[0]==-.02 and c.specific_se.notna().all()
    assert c.FDR_cell_manifest.notna().all()
    assert json.loads((out/'cigma.provenance.json').read_text())['record_type']=='test_double'
    assert len(json.loads((out/'cigma.native_fits.json').read_text()))==1

def test_failed_gene_stays_in_FDR_family(tmp_path):
    row,_,_=native_inputs(tmp_path); bad={**row,'gene':'B','ctnu':'absent.csv'}
    manifest=csv(tmp_path/'manifest.csv',[row,bad]); out=tmp_path/'o'
    assert c5.run(manifest,out,_fit_function=fake_he)==1
    d=pd.read_csv(out/'cigma.results.csv'); assert d.specific_FDR_manifest.iloc[0]==.002
    assert pd.read_csv(out/'cigma.status.csv').status.tolist()==['ok','failed']

def test_native_missing_package_clears_old_outputs(tmp_path,monkeypatch):
    row,_,_=native_inputs(tmp_path); manifest=csv(tmp_path/'manifest.csv',[row]); out=tmp_path/'o'; out.mkdir()
    (out/'cigma.results.csv').write_text('OLD SUCCESS\n')
    def missing(*a): raise ImportError('native package intentionally unavailable in test')
    monkeypatch.setattr(c5,'load_cigma',missing)
    assert c5.run(manifest,out)==1
    assert pd.read_csv(out/'cigma.results.csv').empty

def test_native_trans_contract(tmp_path):
    row,ds,cs=native_inputs(tmp_path)
    rng=np.random.default_rng(23); z=rng.normal(size=(len(ds),60)); k=pd.DataFrame(z@z.T/60,index=ds,columns=ds)
    k.index.name='donor'; k.to_csv(tmp_path/'Kt.csv'); row['Kt']='Kt.csv'
    manifest=csv(tmp_path/'m.csv',[row]); out=tmp_path/'o'
    assert c5.run(manifest,out,_fit_function=fake_he)==0
    assert set(pd.read_csv(out/'cigma.results.csv').component)=={'cis','trans'}

def external_mapping(root,full=False):
    table=csv(root/'external.csv',dict(Gene=['A','B','E'],Pvalue=[.001,.1,.8],Adjusted=[.04,.7,1]))
    mapping=root/'map.json'; mapping.write_text(json.dumps(dict(source='SYNTHETIC published-like example',source_version='test',
        cohort='example',tissue='blood',component='cis',full_test_family=full,
        columns=dict(gene='Gene',specific_p='Pvalue',specific_FDR_manifest='Adjusted'),original_adjustment='Example supplied q')))
    return table,mapping

def test_import_preserves_source_adjustment(tmp_path):
    table,mapping=external_mapping(tmp_path)
    d=c5.import_cigma(table,mapping,tmp_path/'o'); assert d.specific_FDR_manifest.tolist()==[.04,.7,1]

def test_import_partial_no_new_q(tmp_path):
    table,mapping=external_mapping(tmp_path)
    j=json.loads(mapping.read_text()); del j['columns']['specific_FDR_manifest']; mapping.write_text(json.dumps(j))
    d=c5.import_cigma(table,mapping,tmp_path/'o'); assert d.specific_FDR_manifest.isna().all()

def test_import_complete_q_declared_denominator(tmp_path):
    table,mapping=external_mapping(tmp_path,True)
    j=json.loads(mapping.read_text()); del j['columns']['specific_FDR_manifest']; j['n_tests']=10; mapping.write_text(json.dumps(j))
    d=c5.import_cigma(table,mapping,tmp_path/'o'); assert d.specific_FDR_manifest.iloc[0]==.01

def test_precomputed_partial_no_enrichment(inputs,tmp_path):
    table,mapping=external_mapping(tmp_path); c5.import_cigma(table,mapping,tmp_path/'ext')
    out=tmp_path/'o'; out.mkdir(); u,a,p=inputs
    c5.integrate_cigma(tmp_path/'ext/cigma.results.csv',[a,u,p],out)
    assert not (out/'c5.CIGMA_enrichment.csv').exists()
    z=pd.read_csv(out/'c5.CIGMA_annotation.csv'); assert not z.loc[z.assay=='c','CIGMA_available'].item()

def test_precomputed_complete_background_is_tested(inputs,tmp_path):
    table,mapping=external_mapping(tmp_path,True); c5.import_cigma(table,mapping,tmp_path/'ext')
    out=tmp_path/'o'; out.mkdir(); u,a,p=inputs
    c5.integrate_cigma(tmp_path/'ext/cigma.results.csv',[a,u,p],out)
    e=pd.read_csv(out/'c5.CIGMA_enrichment.csv'); assert (e.assayed_tested_background==3).all()

def test_published_gene_calls_are_not_replaced_by_new_bh_threshold(inputs,tmp_path):
    table,mapping=external_mapping(tmp_path,True)
    imported=c5.import_cigma(table,mapping,tmp_path/'ext')
    imported['source_cs_egene']=[False,True,False]
    imported['source_significance_rule']='Explicit source test calls for fixture'
    imported.to_csv(tmp_path/'ext/cigma.results.csv',index=False)
    out=tmp_path/'o';out.mkdir();u,a,p=inputs
    c5.integrate_cigma(tmp_path/'ext/cigma.results.csv',[a,u,p],out)
    e=pd.read_csv(out/'c5.CIGMA_enrichment.csv')
    assert e.significance_rule.str.startswith('Source cs-eGene calls').all()
    assert (e.cs_egene_background==1).all()

def test_test_double_not_native_evidence(inputs,tmp_path):
    row,_,_=native_inputs(tmp_path); manifest=csv(tmp_path/'m.csv',[row]); c5.run(manifest,tmp_path/'ext',_fit_function=fake_he)
    u,a,p=inputs; out=tmp_path/'o'; out.mkdir()
    with pytest.raises(ValueError,match='test double'): c5.integrate_cigma(tmp_path/'ext/cigma.results.csv',[a,u,p],out)

def test_pseudobulk_SEM_squared():
    x=np.array([[1.],[3.],[2.],[6.],[3.],[5.],[4.],[8.]])
    donors=['a']*4+['b']*4; ct=['T','T','B','B']*2
    z=c5.pseudobulk_moments(x,donors,ct,['A'],min_cells=2,min_expression_fraction=0)
    # Cell order B,T; B donor a values 2,6: sample variance=8, variance(mean)=4.
    assert z['ctnu'][0,0,0]==4 and z['ctnu'][0,1,0]==1
    assert z['ctp'][0,0,0]==4
    zz=c5.pseudobulk_moments(sparse.csr_matrix(x),donors,ct,['A'],min_cells=2,min_expression_fraction=0)
    np.testing.assert_equal(z['ctnu'],zz['ctnu'])

def test_pseudobulk_incomplete_donors_not_imputed():
    x=np.arange(1,8).reshape(-1,1)
    z=c5.pseudobulk_moments(x,['a']*4+['b']*3,['T','T','B','B','T','T','B'],['A'],min_cells=2,min_expression_fraction=0)
    assert z['donors'].tolist()==['a'] and z['retained'].tolist()==[True,False]

def test_genotype_kinship_PSD_and_allele_swap():
    rng=np.random.default_rng(17); x=rng.binomial(2,.3,size=(40,30)).astype(float); x[0,0]=np.nan
    K,keep,freq=c5.kinship_from_dosage(x,.05,.1)
    KK,_,_=c5.kinship_from_dosage(2-x,.05,.1)
    np.testing.assert_allclose(K,KK,atol=1e-12)
    assert np.linalg.eigvalsh(K).min()>-1e-10

def test_invalid_genotype_dosage():
    with pytest.raises(ValueError): c5.kinship_from_dosage(np.array([[0,3],[1,1],[2,0]]),.05,.1)

def test_driver_cache_and_met_status(tmp_path,inputs):
    u,a,p=inputs; root=tmp_path/'analysis'
    args=['--Y','cvd_cad','--biom','prot','--analysis-root',str(root),'--atlas',str(a),'--universe',str(u),'--panels',str(p),'--no-plots']
    assert c5.main_driver(args)==0
    out=root/'cvd_cad/prot/c5_cellulation'; t=(out/'c5.cell.enrichment.csv').stat().st_mtime_ns
    assert c5.main_driver(args)==0 and (out/'c5.cell.enrichment.csv').stat().st_mtime_ns==t
    assert c5.main_driver(['--Y','cvd_cad','--biom','met','--analysis-root',str(root),'--no-plots'])==0
    d=pd.read_csv(root/'cvd_cad/met/c5_cellulation/c5.cellulation_status.csv'); assert set(d.status)=={'not_applicable'}

def test_driver_visible_partial_and_complete_status(tmp_path,inputs,capsys):
    u,a,p=inputs; root=tmp_path/'analysis'
    args=['--Y','cad','--biom','prot','--analysis-root',str(root),'--atlas',str(a),'--universe',str(u),'--panels',str(p),'--no-plots']
    for _ in range(2):
        assert c5.main_driver(args)==0
        log=capsys.readouterr().out
        assert '[LE8] START c5_cellulation' in log
        assert '[LE8] DONE C5/cell_expression' in log
        assert '[LE8] SKIP C5/CIGMA' in log
        assert '[LE8] PARTIAL c5_cellulation' in log
    table,mapping=external_mapping(tmp_path)
    c5.import_cigma(table,mapping,tmp_path/'external')
    assert c5.main_driver([*args,'--cigma-results',str(tmp_path/'external/cigma.results.csv')])==0
    log=capsys.readouterr().out
    assert '[LE8] DONE C5/CIGMA' in log and '[LE8] DONE c5_cellulation' in log
    assert '[LE8] PARTIAL' not in log

def test_failed_rerun_no_stale_annotation(tmp_path,inputs):
    u,a,p=inputs; root=tmp_path/'analysis'; args=['--Y','cad','--biom','prot','--analysis-root',str(root),
        '--atlas',str(a),'--universe',str(u),'--panels',str(p),'--no-plots']
    assert c5.main_driver(args)==0
    x=pd.read_csv(p); x.loc[0,'feature']='MISSING'; x.to_csv(p,index=False)
    assert c5.main_driver(args)==1
    out=root/'cad/prot/c5_cellulation'
    assert not (out/'c5.completed.json').exists() and not (out/'c5.cell.enrichment.csv').exists()

def test_cli_annotation(tmp_path,inputs):
    u,a,p=inputs; r=subprocess.run([sys.executable,str(SCRIPT),'annotate','--universe',str(u),
        '--atlas',str(a),'--panels',str(p),'--outdir',str(tmp_path/'cli'),'--no-plots'],capture_output=True,text=True)
    assert r.returncode==0,r.stderr

def test_plot_no_combined_fold_statistics(tmp_path,inputs):
    c5.annotate(*inputs,tmp_path/'plots',prefix='c5.cell',plots=True)
    assert (tmp_path/'plots/c5.cell.enrichment.png').stat().st_size>1000

def test_native_validation_only_not_inference(tmp_path):
    row,_,_=native_inputs(tmp_path); manifest=csv(tmp_path/'m.csv',[row])
    assert c5.run(manifest,tmp_path/'out',validate_only=True)==0
    assert pd.read_csv(tmp_path/'out/cigma.results.csv').empty
    assert json.loads((tmp_path/'out/cigma.provenance.json').read_text())['validate_only'] is True

def test_evidence_keeps_distinct_loci(inputs,tmp_path):
    u,a,p=inputs; root=tmp_path/'disease'; dest=root/'prot/c3_coloc'; dest.mkdir(parents=True)
    csv(dest/'c3.coloc_summary.csv',dict(feature=['a','a'],locus=['locus1','locus2'],
        **{'PP.H4':[.9,.1]},locus_class=['cis','trans']))
    out=tmp_path/'out'; c5.annotate(u,a,p,out,prefix='c5.cell',plots=False)
    c5.evidence_long(root,u,out)
    d=pd.read_csv(out/'c5.evidence_cell_context.csv')
    assert len(d)==2 and set(d.locus)=={'locus1','locus2'}
    assert set(d.cell_type)=={'T'}

def test_pseudobulk_no_phantom_constant_gene():
    x=np.ones((8,1)); donors=['a']*4+['b']*4; ct=['T','T','B','B']*2
    z=c5.pseudobulk_moments(x,donors,ct,['A'],min_cells=2,min_expression_fraction=0)
    assert z['usable'].tolist()==[False]

def test_published_families_mapping_and_missing_cell_tests(tmp_path, monkeypatch):
    annotation=tmp_path/'synthetic.gtf.gz'
    with gzip.open(annotation,'wt') as f:
        for gene,name in [('ENSG1','AAA'),('ENSG2','BBB')]:
            f.write(f'1\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "{gene}"; gene_name "{name}";\n')
    workbook=tmp_path/'synthetic.xlsx'
    rows=pd.DataFrame({'gene':['ENSG1','ENSG2'],'p:V':[.001,.8],'p:sigma_g2':[.1,.5],
                       'sigma_g2':[.2,.1],'se:sigma_g2':[.1,.1],'specificity':[.5,1.4],
                       'V_B':[.2,-.1],'se:V_B':[.1,.1],'p:V_B':[.01,.9]})
    with pd.ExcelWriter(workbook) as w:
        rows.to_excel(w,sheet_name='OneK1K fixture',index=False,startrow=2)
        rows.drop(columns='p:V_B').to_excel(w,sheet_name='CLUES fixture',index=False,startrow=2)
    monkeypatch.setattr(module,'TABLES', [('OneK1K fixture','OneK1K',2,['B']),
                                         ('CLUES fixture','CLUES_ImmVar_meta',2,['B'])])
    out=tmp_path/'outputs'
    module.run(workbook,annotation,out)
    genes=pd.read_csv(out/'cigma.results.csv');cells=pd.read_csv(out/'cigma.cell_types.csv')
    assert len(genes)==4 and set(genes.gene)=={'AAA','BBB'}
    assert genes.source_cs_egene.tolist()==[True,False,True,False]
    assert genes.loc[genes.gene=='BBB','specificity'].isna().all()
    assert cells.loc[cells.gene=='BBB','specific_variance'].eq(-.1).all()
    assert cells.loc[cells.cohort=='CLUES_ImmVar_meta',['p','FDR_cell_manifest']].isna().all().all()
    provenance=json.loads((out/'cigma.provenance.json').read_text())
    assert len(provenance['families'])==2 and provenance['full_test_family']

class C5WorkflowRegression(TestCase):
	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory(prefix='le8-acceptance-', dir='/tmp')
		self.addCleanup(self.tmp.cleanup)
		self.path = Path(self.tmp.name)
		self.clean = mock.patch.dict(os.environ, {k:v for k,v in os.environ.items() if not k.startswith(('LE8_', 'C4_', 'C3_', 'PGS_', 'DATE_FOLLOW_END'))}, clear=True)
		self.clean.start()
		self.addCleanup(self.clean.stop)

	def cell_inputs(self):
		u=pd.DataFrame({'assay':[f'P{i}' for i in range(12)],'gene':[f'GENE{i}' for i in range(12)],'match_stratum':'same'})
		a=pd.DataFrame({'gene':[f'GENE{i}' for i in range(10)],'cell_type':['A']*5+['B']*5,'source':'fixture','source_version':'v1'})
		p=pd.DataFrame({'model':['YS']*3+['NS']*3,'feature':['P0','P1','P2','P5','P6','P7']})
		c=pd.DataFrame({'model':['YS'],'reference':['NS']})
		for name,d in [('universe',u),('atlas',a),('panels',p),('contrasts',c)]: d.to_csv(self.path/f'{name}.csv',index=False)
		return [item for name in ['universe','atlas','panels','contrasts'] for item in ['--'+name,str(self.path/f'{name}.csv')]]

	def test_c5_dispatch_and_direct_signatures_identical_and_seed_changes(self):
		args=self.cell_inputs();out=self.path/'run'
		base=['--Y','cad','--biom','prot','--analysis-root',str(out),'--no-plots','--seed','57','--matched-draws','100',*args]
		self.assertEqual(c5.main_driver(base),0)
		file=out/'cad/prot/c5_cellulation/c5.completed.json'
		first=json.loads(file.read_text())['signature']
		def execute(cmd,env,dry):
			with mock.patch.dict(os.environ,env,clear=True): self.assertEqual(c5.main_driver(list(map(str,cmd[2:]))),0)
		with mock.patch.object(shared,'call',execute): shared.dispatch_main(['c5_cellulation',*base])
		self.assertEqual(first,json.loads(file.read_text())['signature'])
		base[base.index('--seed')+1]='58';self.assertEqual(c5.main_driver(base),0)
		self.assertNotEqual(first,json.loads(file.read_text())['signature'])

	def test_c5_staged_and_published_evidence_have_identical_provenance(self):
		self.cell_inputs()
		published=self.path/'published';work=self.path/'staging'
		d=pd.DataFrame({'feature':['P0'],'locus':['locus1'],'PP.H4':[.8],'locus_class':['cis']})
		for root in [published,work]:
			f=root/'cad/prot/c3_coloc/c3.coloc_summary.csv';f.parent.mkdir(parents=True);d.to_csv(f,index=False)
		with mock.patch.dict(os.environ,LE8_ANALYSIS_ROOT=str(work),LE8_PUBLISHED_ROOT=str(published)):
			a=c5.stamp(work/'cad/prot/c3_coloc/c3.coloc_summary.csv')
			c5.evidence_long(work/'cad',self.path/'universe.csv',self.path/'staged-evidence')
		b=c5.stamp(published/'cad/prot/c3_coloc/c3.coloc_summary.csv')
		c5.evidence_long(published/'cad',self.path/'universe.csv',self.path/'published-evidence')
		self.assertEqual(a,b)
		for file in ['c5.evidence_long.csv','c5.evidence_sources.csv']:
			self.assertEqual((self.path/'staged-evidence'/file).read_bytes(),(self.path/'published-evidence'/file).read_bytes())

	def test_c5_not_requested_is_explicit(self):
		args=self.cell_inputs();c5.annotate(self.path/'universe.csv',self.path/'atlas.csv',self.path/'panels.csv',self.path/'annotation',plots=False)
		d=pd.read_csv(self.path/'annotation/cell.contrast_status.csv')
		self.assertEqual(d.status.iloc[0],'not_requested')

	def test_c5_cross_fold_contrast_rejected(self):
		self.cell_inputs();p=pd.read_csv(self.path/'panels.csv');p['fold']=['a']*3+['b']*3;p.to_csv(self.path/'panels.csv',index=False)
		pd.DataFrame({'model':['YS|fold=a'],'reference':['NS|fold=b']}).to_csv(self.path/'contrasts.csv',index=False)
		with self.assertRaisesRegex(ValueError,'same fold'):
			c5.annotate(self.path/'universe.csv',self.path/'atlas.csv',self.path/'panels.csv',self.path/'annotation',contrasts=self.path/'contrasts.csv',plots=False)

	def test_invalid_c5_settings_fail_at_dispatch(self):
		for args in [['--matched-draws','1'],['--cigma-cells','/tmp/x'],['--universe','/tmp/x']]:
			with self.assertRaises(SystemExit): shared.dispatch_main(['c5_cellulation','--dry-run',*args])


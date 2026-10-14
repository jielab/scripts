from pathlib import Path
import sys,json,sqlite3,subprocess,importlib.util
from dataclasses import replace
import numpy as np,pandas as pd,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
from twin_core import Config,id_hash,ROLE_NAMES
from twin_host import numeric_cohort,load_source,verified_prsformer
from twin_io import write_aggregate,public_table,identified,read_table,sha256,json_write
from twin_cli import benchmark_manifest
from browser_index import Index


def test_string_ids_and_duplicates(tmp_path):
    p=tmp_path/'input.tsv';p.write_text('eid\ty\n0001\t2\n0002\t3\n')
    d=identified(read_table(p));assert d.eid.tolist()==['0001','0002']
    with pytest.raises(ValueError):identified(pd.concat([d,d]))


def test_native_outer_roster_and_missing_test_y_preserved():
    rng=np.random.default_rng(1);n=1000
    d=pd.DataFrame(dict(eid=[str(i) for i in range(n)],split=np.where(np.arange(n)%2,'train','test'),y=rng.normal(size=n),s1=rng.normal(size=n),s2=rng.normal(size=n),age=rng.normal(size=n)))
    d.loc[:9,'y']=np.nan
    raw,q,t,flow,cfg,extra=numeric_cohort(d,{'csx':['s1','s2'],'covariates':['age']},Config())
    assert len(q['ids'])==500 and flow['test_missing_y_N']==5
    assert flow['development_missing_y_N']==5
    assert cfg.grouping_status=='individual_unverified'
    assert not {'y','outcome'}&set(q)
    assert q['ids'].tolist()==d.loc[d.split=='test','eid'].tolist()


def test_supplied_family_cannot_cross_outer_split():
    d=pd.DataFrame(dict(eid=['a','b'],family_id=['F','F'],split=['train','test'],y=[1.,2.],s=[1.,2.]))
    with pytest.raises(ValueError,match='cross'):numeric_cohort(d,{'csx':['s']},Config())


def test_no_outcome_feature_allowed():
    d=pd.DataFrame(dict(eid=['a','b'],split=['train','test'],y=[1.,2.]))
    with pytest.raises(ValueError,match='predictors'):numeric_cohort(d,{'csx':['y']},Config())


def test_public_tables_reject_individual_columns(tmp_path):
    with pytest.raises(ValueError):write_aggregate(tmp_path,{'bad':pd.DataFrame({'eid':['x'],'score':[.5]})},{'trait':'height'})
    with pytest.raises(ValueError):public_table(pd.DataFrame({'donor_id':['x']}))


def test_index_reads_aggregate_only_and_tracks_hashes(tmp_path):
    root=tmp_path/'analysis';out=root/'evidence'
    folder=root/'abm/twin/height/current';folder.mkdir(parents=True)
    json_write(folder.parent/'CURRENT.json',{'run':'current'});json_write(folder/'COMPLETE.json',{'status':'completed'})
    data=pd.DataFrame({'method':['ABM_clean','CSx_reference'],'subset':['selected','selected'],'actual_coverage':[.5,.5],'MSE':[1.,2.]})
    write_aggregate(folder/'test',{'metrics':data},{'trait':'height','module':'Twin test','run_id':'frozen'})
    private=folder/'_LOCAL_ONLY';private.mkdir();(private/'confidence.tsv').write_text('eid\tscore\nx\t1\n')
    old=root/'abm/twin/height/old/test';old.mkdir(parents=True);(old/'metrics.tsv').write_text('method\tMSE\nold\t9\n')
    (root/'final').mkdir();(root/'final/private.tsv').write_text('eid\tvalue\nx\t2\n')
    result=Index(root,out).run(['height']);assert result['tables']==1
    con=sqlite3.connect(out/'evidence.sqlite');cat=pd.read_sql_query('select * from catalogue',con);audit=pd.read_sql_query('select * from audit',con);con.close()
    assert cat.iloc[0]['legacy']==0 and cat.iloc[0]['name']=='metrics'
    assert 'blocked_private' in set(audit.status)


def test_manifest_corruption_is_not_fallback_import(tmp_path):
    # A failed manifest's owned file must NOT be ingested as an unverified plain TSV.
    root=tmp_path/'analysis';run=root/'abm/twin/height/run1';run.mkdir(parents=True)
    json_write(run.parent/'CURRENT.json',{'run':'run1'});json_write(run/'COMPLETE.json',{})
    write_aggregate(run/'test',{'metrics':pd.DataFrame({'method':['M'],'MSE':[1.]})},{'trait':'height'})
    (run/'test/metrics.tsv.gz').write_bytes(b'corrupt')
    result=Index(root,root/'evidence').run()
    con=sqlite3.connect(result['db']);audit=pd.read_sql_query('select * from audit',con);cat=pd.read_sql_query('select * from catalogue',con);con.close()
    assert cat.empty and 'checksum_mismatch' in set(audit.status)


def test_symlink_outside_root_not_imported(tmp_path):
    root=tmp_path/'analysis';root.mkdir();external=tmp_path/'outside.tsv';external.write_text('method\tMSE\nprivate\t3\n')
    (root/'benign.tsv').symlink_to(external)
    assert Index(root,root/'evidence').run()['tables']==0


def test_benchmark_manifest_checks_roster_units_and_hash(tmp_path):
    ids=np.array(['001','002']);d=pd.DataFrame({'eid':ids,'p':[1.,2.],'outcome':[1.,2.]});path=tmp_path/'pred.tsv';d.to_csv(path,sep='\t',index=False)
    entry=dict(name='Original_PRS_CSx',path='pred.tsv',trait='height',endpoint_definition='height_cm',value_type='outcome_scale_prediction',test_roster_sha256=id_hash(ids),predictions_frozen_without_test_labels=True,sha256=sha256(path),prediction_column='p',outcome_column='outcome')
    mf=tmp_path/'methods.json';json_write(mf,{'methods':[entry]});truth=pd.DataFrame({'eid':ids,'y':[1.,2.]})
    p,a=benchmark_manifest(mf,ids,truth,'height','height_cm');np.testing.assert_array_equal(p['Original_PRS_CSx'],[1.,2.])
    entry['value_type']='raw_PRS';json_write(mf,{'methods':[entry]})
    with pytest.raises(ValueError,match='units'):benchmark_manifest(mf,ids,truth,'height','height_cm')


def test_installer_protects_upstream_and_can_reinstall(tmp_path):
    installer=Path(__file__).resolve().parents[2]/'install.py'
    # In the distributable, install.py is one level above grid/.
    if not installer.is_file():pytest.skip('Installer test runs from the distributable, not the installed scripts root')
    original=tmp_path/'scripts/grid';original.mkdir(parents=True)
    paths=['grid.sh','0.pca.sh','1.csx.sh','2.disco.sh','f/main.py','f/data.py','f/prsformer.py','f/prsformer.sh','f/abm.py']
    for name in paths:
        p=original/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('# ORIGINAL '+name+'\n')
    before={name:(original/name).read_bytes() for name in paths}
    subprocess.run([sys.executable,str(installer),'--grid-dir',str(original)],check=True,capture_output=True)
    assert (original/'grid.upstream.sh').read_bytes()==before['grid.sh']
    for name in paths[1:]:assert (original/name).read_bytes()==before[name]
    subprocess.run([sys.executable,str(installer),'--grid-dir',str(original)],check=True,capture_output=True)
    assert (original/'grid.upstream.sh').read_bytes()==before['grid.sh']


def test_legacy_xlsx_cached_values_not_formula_execution(tmp_path):
    import zipfile
    from twin_io import xlsx_values
    p=tmp_path/'legacy.xlsx'
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    payload={
      'xl/workbook.xml':f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="performance" sheetId="1" r:id="rId1"/></sheets></workbook>',
      'xl/_rels/workbook.xml.rels':'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>',
      'xl/sharedStrings.xml':f'<sst xmlns="{ns}"><si><t>method</t></si><si><t>estimate</t></si><si><t>PRS-CSx</t></si></sst>',
      'xl/worksheets/sheet1.xml':f'<worksheet xmlns="{ns}"><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row><row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2"><f>1+1</f><v>2</v></c></row></sheetData></worksheet>'}
    with zipfile.ZipFile(p,'w') as z:
        for k,v in payload.items():z.writestr(k,v)
    d=xlsx_values(p)['performance'];assert d.method.tolist()==['PRS-CSx'] and d.estimate.tolist()==[2.]


def test_shared_PC_kept_in_covariate_baselines_without_duplicate_column():
    rng=np.random.default_rng(10);n=1000
    d=pd.DataFrame(dict(eid=[str(i) for i in range(n)],split=np.where(np.arange(n)%2,'train','test'),y=rng.normal(size=n),s1=rng.normal(size=n),PC1=rng.normal(size=n),age=rng.normal(size=n)))
    raw,_,_,flow,_,_=numeric_cohort(d,{'csx':['s1'],'ancestry':['PC1'],'covariates':['age','PC1']},Config())
    assert raw['features'].count('PC1')==1
    assert raw['features'].index('PC1') in raw['genetic_indices']
    assert raw['features'].index('PC1') in raw['covariate_indices']
    assert flow['clinical_columns']==['age','PC1']


def test_protocol_metadata_does_not_expand_participant_dataframes(tmp_path):
    p=tmp_path/'protocol.json';json_write(p,{'metadata':pd.DataFrame({'eid':['DO_NOT_EMBED_THIS_ID']})})
    text=p.read_text();assert 'DO_NOT_EMBED_THIS_ID' not in text
    assert json.loads(text)['metadata']['contents_omitted']


def test_prepared_cache_cannot_drop_original_test_participant(tmp_path):
    import joblib
    cohort=pd.DataFrame({'eid':['a','b','c'],'split':['train','test','test']})
    table=cohort.iloc[:2].assign(y=[1.,2.],score=[1.,2.])
    path=tmp_path/'prepared.joblib'
    joblib.dump({'tables':{'height':table},'feature_groups':{'height':{'csx':['score']}},'cohort':cohort},path)
    with pytest.raises(ValueError,match='entire frozen original cohort'):load_source(path,'height')


def test_prepared_cache_cannot_change_frozen_family_groups(tmp_path):
    import joblib
    cohort=pd.DataFrame({'eid':['a','b'],'split':['train','test'],'family_id':['F','G']})
    table=cohort.assign(y=[1.,2.],score=[1.,2.]);table.loc[0,'family_id']='replacement'
    path=tmp_path/'prepared.joblib'
    joblib.dump({'tables':{'height':table},'feature_groups':{'height':{'csx':['score']}},'cohort':cohort},path)
    with pytest.raises(ValueError,match='reference groups differ'):load_source(path,'height')


@pytest.mark.parametrize('failure',[None,'missing_endpoint','covariate','missing_test','split'])
def test_prsformer_exact_published_schema(tmp_path,monkeypatch,failure):
    import twin_host
    roster_path=tmp_path/'prsformer/split.rds';roster_path.parent.mkdir();roster_path.touch()
    path=tmp_path/'prsformer/scores/height/predictions.rds';path.parent.mkdir(parents=True);path.touch()
    cohort=pd.DataFrame({'eid':['a','b','c','d'],'split':['train','train','test','test'],'age':[30.,40.,50.,60.]})
    roster=cohort[['eid','split']].copy();roster.loc[1,'split']='validation'
    truth=pd.DataFrame({'eid':['c','d'],'y':[170.,np.nan]})
    pred=pd.DataFrame({'eid':['d','c'],'split':'test','prediction':[172.,171.],'outcome':[np.nan,170.],
                       'covariate.age':[60.,50.],'endpoint_definition':'height','covariate_names':'age'})
    if failure=='missing_endpoint':pred.loc[0,'endpoint_definition']=None
    elif failure=='covariate':pred.loc[0,'covariate.age']=99
    elif failure=='missing_test':pred=pred.iloc[:1]
    elif failure=='split':roster.loc[1,'split']='test'
    monkeypatch.setattr(twin_host,'read_table',lambda p:roster if p==roster_path else pred)
    def read():return verified_prsformer(tmp_path,'height',cohort,truth,Config(),['age'],'height')
    if failure:
        with pytest.raises(ValueError):read()
    else:
        values,status=read();np.testing.assert_array_equal(values,[171.,172.])
        assert status['status']=='same_roster_endpoint_and_covariates_verified'


def test_index_does_not_traverse_cache_or_old_runs(tmp_path,monkeypatch):
    root=tmp_path/'analysis';current=root/'abm/twin/height/now';current.mkdir(parents=True)
    (root/'prsformer/cache').mkdir(parents=True);(current.parent/'old').mkdir()
    json_write(current.parent/'CURRENT.json',{'run':'now'});json_write(current/'COMPLETE.json',{})
    write_aggregate(current/'test',{'metrics':pd.DataFrame({'method':['M'],'MSE':[1.]})},{'trait':'height'})
    visited=[]
    import browser_index
    walk=browser_index.os.walk
    def tracking(*args,**kwargs):
        for item in walk(*args,**kwargs):
            visited.append(Path(item[0]));yield item
    monkeypatch.setattr(browser_index.os,'walk',tracking)
    assert Index(root,root/'evidence').run()['tables']==1
    assert root/'prsformer/cache' not in visited
    assert current.parent/'old' not in visited


def test_broken_twin_manifest_never_falls_back_to_raw_table(tmp_path):
    run=tmp_path/'abm/twin/height/current';run.mkdir(parents=True)
    json_write(run.parent/'CURRENT.json',{'run':'current'});json_write(run/'COMPLETE.json',{})
    (run/'metrics.tsv').write_text('method\tMSE\nM\t1\n')
    (run/'table_manifest.json').write_text('{invalid')
    assert Index(tmp_path,tmp_path/'evidence').run()['tables']==0


def test_failed_reindex_preserves_published_database(tmp_path,monkeypatch):
    (tmp_path/'metrics.tsv').write_text('method\tMSE\nM\t1\n')
    result=Index(tmp_path,tmp_path/'evidence').run()
    database=Path(result['db']);before=database.read_bytes()
    index=Index(tmp_path,tmp_path/'evidence')
    def fail():raise RuntimeError('Interrupted import')
    monkeypatch.setattr(index,'sources',fail)
    with pytest.raises(RuntimeError,match='Interrupted import'):index.run()
    assert database.read_bytes()==before
    assert not index.tmp.exists()

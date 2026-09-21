import csv,gzip,json,sqlite3,tempfile,unittest
from collections import defaultdict,Counter
from pathlib import Path
from ibdmix_final_filter import carrier_mask,filter_rows,load_daf,prepare

class FinalFilterTests(unittest.TestCase):
    def test_cell_summary_is_altai_only_and_excludes_x(self):
        from introgression_density import neanderthal_summary
        con=sqlite3.connect(':memory:')
        con.execute('CREATE TABLE segments(sample_id,chr,start,end,dataset_id,genome_build,method,source)')
        con.executemany('INSERT INTO segments VALUES(?,?,?,?,?,?,?,?)',[
            ('a','1',0,100,'1kg','GRCh37','ibdmix','Altai'),
            ('a','1',50,150,'1kg','GRCh37','ibdmix','Altai'),
            ('a','1',200,400,'1kg','GRCh37','ibdmix','Vindija'),
            ('a','X',0,500,'1kg','GRCh37','ibdmix','Altai')])
        data=neanderthal_summary(con,'1kg','GRCh37',['a','b'],{'1':1000},{'1':{'a','b'}})
        self.assertEqual(data[0]['neanderthal_bp'],150)
        self.assertEqual(data[1]['neanderthal_bp'],0)
        con.close()

    def test_carriers_not_duplicate_calls_and_zero_calls_in_denominator(self):
        calls={'a':[(0,100),(20,80)],'b':[(50,150)],'c':[(70,90)]}
        self.assertEqual(carrier_mask(calls,set('abcdefghij')),[(70,90)])
        self.assertEqual(carrier_mask({'a':[(0,100),(0,100)]},set('abcdefghij')),[])

    def test_boundary_thirty_percent_and_half_open(self):
        self.assertEqual(carrier_mask({'a':[(0,10)],'b':[(10,20)]},set('abc'),.30),[(0,20)])
        self.assertEqual(carrier_mask({'a':[(0,10)]},set('abcdefg'),.30),[])

    def test_masks_split_and_reapply_length_and_keep_parent_score(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'calls.gz'
            fields=['anc','chrom','start','end','length','slod','sites','positive_lods','negative_lods','score_scope']
            with gzip.open(p,'wt') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
                for ref in ['Altai','Denisova','Vindija']:
                    w.writerow(dict(anc=ref,chrom='1',start=0,end=200000,length=200000,slod=8,sites=400,positive_lods=20,negative_lods=5,score_scope='parent_native_call'))
            counts=defaultdict(Counter)
            result=list(filter_rows(p,{'Denisova':[(50000,150001)]},{'1':[(60000,130000)]},'1',counts))
            self.assertEqual([(r['anc'],r['start'],r['end']) for r in result], [('Altai',0,60000),('Altai',130000,200000),('Denisova',0,50000),('Vindija',0,200000)])
            self.assertEqual(result[0]['slod'],'8');self.assertEqual(result[0]['sites'],'')
            self.assertEqual(counts['Denisova']['short_fragment_bp'],49999)

    def test_daf_requires_provenance_and_valid_build(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'mask.bed';p.write_text('chr1\t5\t10\n')
            m=Path(str(p)+'.json');m.write_text(json.dumps(dict(genome_build='GRCh38',reference='Altai',percentile=99.9,source='test fixture')))
            with self.assertRaises(ValueError):load_daf(p)
            m.write_text(json.dumps(dict(genome_build='GRCh37',reference='Altai',percentile=99.9,source='test fixture')))
            self.assertEqual(load_daf(p)[0],{'1':[(5,10)]})
            self.assertEqual(load_daf(None)[1]['status'],'not_applied_missing_author_mask')

    def test_raw_overlap_and_cache_invalidation(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=root/'chr1';(run/'final').mkdir(parents=True);sample=run/'samples/C1';sample.mkdir(parents=True);raw=run/'raw/C1';raw.mkdir(parents=True)
            pops=['ESN','GWD','LWK','MSL','YRI','CEU']
            with (sample/'populations.tsv').open('w') as f:
                f.write('population\tn\n')
                for pop in pops:
                    f.write(f'{pop}\t1\n');(sample/(pop+'.txt')).write_text(pop+'\n')
            for ref in ['Altai','Denisova']:
                for pop in pops:
                    with gzip.open(raw/f'{ref}.{pop}.raw.txt.gz','wt') as f:
                        f.write('ID\tchrom\tstart\tend\tslod\n')
                        if ref=='Altai' and pop=='CEU':f.write(f'{pop}\t1\t50001\t100001\t6\n')
                        if ref=='Denisova' and pop in ['ESN','GWD']:f.write(f'{pop}\t1\t150001\t200001\t6\n')
            source=run/'final/all_archaic_refs.segments.tsv.gz'
            with gzip.open(source,'wt') as f:
                f.write('ID\tchrom\tstart\tend\tlength\tslod\tanc\tscore_scope\nCEU\t1\t0\t250000\t250000\t6\tDenisova\tparent_native_call\n')
            (run/'run.meta.tsv').write_text('fixture')
            meta=dict(genome_build='b37',coordinate_system='0-based-half-open',refs='Altai Denisova')
            p,a=prepare(source,root/'cache',meta)
            with gzip.open(p,'rt') as f:r=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual([(x['start'],x['end']) for x in r],[('0','50000'),('100000','150000'),('200000','250000')])
            self.assertEqual(a['n_african_controls'],5)
            self.assertEqual(prepare(source,root/'cache',meta)[0],p)
            self.assertEqual(source.exists(),True)
            (sample/'CEU.txt').write_text('CEU\n\n')
            self.assertNotEqual(prepare(source,root/'cache',meta)[0],p)

if __name__=='__main__':unittest.main()

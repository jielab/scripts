import csv,sys,sqlite3,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
import phyml_report as report
import dual_lead as dl

class ValidationTests(unittest.TestCase):
    def database(self):
        db=sqlite3.connect(':memory:')
        db.executescript('''CREATE TABLE method_runs(dataset_id,genome_build,chr,method,status,evidence_eligible,availability_note);
        CREATE TABLE segments(dataset_id,genome_build,chr,method,sample_id,source_class,start,end,haplotype);''')
        db.executemany('INSERT INTO method_runs VALUES(?,?,?,?,?,?,?)',[
            ('1kg','GRCh37','6','ibdmix','complete',1,''),('1kg','GRCh37','6','trace','complete',1,'')])
        return db
    def candidate(self):
        return dict(candidate_id='test',locus_key='6',dataset_id='1kg',genome_build='GRCh37',chr='6',lineage='Denisovan',sample_id='A',haplotype='1',candidate_start=100,candidate_end=200)
    def test_same_individual_lineage_and_single_segment(self):
        db=self.database()
        db.executemany('INSERT INTO segments VALUES(?,?,?,?,?,?,?,?,?)',[
            ('1kg','GRCh37','6','ibdmix','B','Denisovan',100,200,None),
            ('1kg','GRCh37','6','ibdmix','A','Neanderthal',100,200,None),
            ('1kg','GRCh37','6','ibdmix','A','Denisovan',100,150,None),
            ('1kg','GRCh37','6','ibdmix','A','Denisovan',150,200,None)])
        rr,_=report.validate([self.candidate()],db)
        r=next(x for x in rr if x['method']=='ibdmix')
        self.assertEqual(r['union_fraction'],1)
        self.assertEqual(r['overlap_fraction'],.5)
        self.assertEqual(r['overlap_pass'],0)
    def test_trace_wrong_haplotype_is_not_phase_support(self):
        db=self.database();db.execute('INSERT INTO segments VALUES(?,?,?,?,?,?,?,?,?)',('1kg','GRCh37','6','trace','A','Ghost/Unknown',100,200,'2'))
        rr,_=report.validate([self.candidate()],db);r=next(x for x in rr if x['method']=='trace')
        self.assertEqual(r['overlap_pass'],1);self.assertEqual(r['phase_overlap_pass'],0)
        self.assertIn('not_lineage_validation',r['lineage_note'])
    def test_missing_method_is_not_zero(self):
        db=self.database();db.execute('DELETE FROM method_runs')
        rr,_=report.validate([self.candidate()],db)
        self.assertTrue(all(r['overlap_pass'] is None for r in rr))
    def test_trace_only_selected_nodes_are_evaluable(self):
        with tempfile.TemporaryDirectory(prefix='gu-report-test-') as tmp:
            root=Path(tmp);(root/'samples').mkdir()
            (root/'run.meta.tsv').write_text('chromosomes\t6\n')
            (root/'samples/trace_sample_map.tsv').write_text('tree_node_id\tsample\thaplotype\n0\tA\t1\n1\tB\t1\n')
            (root/'samples/tree_nodes.txt').write_text('1\n')
            db=self.database();db.execute('ALTER TABLE method_runs ADD COLUMN raw_file')
            db.execute('UPDATE method_runs SET raw_file=? WHERE method="trace"',(str(root/'run.meta.tsv'),))
            rr,_=report.validate([self.candidate()],db)
            r=next(r for r in rr if r['method']=='trace')
            self.assertEqual(r['comparison_available'],0);self.assertIsNone(r['phase_overlap_pass'])
    def test_scoped_run_cannot_imply_whole_chromosome_negative(self):
        with tempfile.TemporaryDirectory(prefix='gu-report-test-') as tmp:
            root=Path(tmp);(root/'run.meta.tsv').write_text('loci_file\\tregion.bed:1:1\n')
            (root/'request.loci.analysis.bed').write_text('6\t500\t1000\tother\n')
            db=self.database();db.execute('ALTER TABLE method_runs ADD COLUMN raw_file')
            db.execute('UPDATE method_runs SET raw_file=? WHERE method="ibdmix"',(str(root/'run.meta.tsv'),))
            rr,_=report.validate([self.candidate()],db)
            r=next(r for r in rr if r['method']=='ibdmix')
            self.assertEqual(r['comparison_available'],0);self.assertIsNone(r['overlap_pass'])
    def test_lineages_do_not_collapse(self):
        rr=[dict(self.candidate(),method='ibdmix',overlap_pass=1,phase_overlap_pass=None)]
        out=[dict(locus_key='6',dataset_id='1kg',genome_build='GRCh37',chr='6',lineage=l) for l in ['Denisovan','Neanderthal']]
        run={('1kg','GRCh37','6','ibdmix'):dict(status='complete',evidence_eligible=1)}
        report.attach_validation(out,rr,run)
        self.assertEqual(out[0]['ibdmix_supported_individuals'],1)
        self.assertIsNone(out[1]['ibdmix_supported_individuals'])

class ActualDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).parent/'final_integration/review'
        with (root/'phyml_locus_report.tsv').open(encoding='utf-8') as f:cls.s=list(csv.DictReader(f,delimiter='\t'))
        with (root/'phyml_haplotype_report.tsv').open(encoding='utf-8') as f:cls.h=list(csv.DictReader(f,delimiter='\t'))
    def test_all_loci_and_all_lineages(self):
        self.assertEqual(len(self.s),16);self.assertEqual(len({r['locus_key'] for r in self.s}),8)
        expected={('1','Neanderthal'):581,('3','Neanderthal'):363,('6','Denisovan'):42,('6','Neanderthal'):100,
                  ('8','Neanderthal'):747,('12','Denisovan'):1016,('X','Neanderthal'):84,('X','Denisovan'):84}
        for key,n in expected.items():
            r=next(r for r in self.s if (r['chr'],r['lineage'])==key)
            self.assertEqual(int(r['n_candidate_copies']),n)
    def test_controls_not_omitted_or_called_positive(self):
        for r in self.s:
            if r['chr'] in ['9','19']:
                self.assertEqual(r['n_candidate_copies'],'0');self.assertNotEqual(r['call'],'supported_candidate')
    def test_x_exact_alias_and_risk(self):
        r=next(r for r in self.s if r['chr']=='X' and r['lineage']=='Neanderthal')
        self.assertEqual(r['input_lead_snp'],'rs6525167');self.assertEqual(r['input_pos_1based'],'66423003')
        self.assertEqual(r['n_risk_copies'],'84')
    def test_candidate_tree_per_lineage(self):
        dd={r['lineage']:r for r in self.s if r['chr']=='6'}
        self.assertEqual(dd['Denisovan']['tree_pass'],'1');self.assertEqual(dd['Neanderthal']['tree_pass'],'0')
    def test_every_candidate_has_dual_marker_fields(self):
        rr=[r for r in self.h if r['diagnostic_candidate_pass']=='1']
        self.assertTrue(rr)
        for r in rr:
            self.assertTrue(r['input_lead_snp']);self.assertTrue(r['best_tag_snp'])
            self.assertGreaterEqual(float(r['best_tag_r2']),0)
    def test_comparators_never_promoted(self):
        for r in self.h:
            if r['diagnostic_candidate_pass']=='0': self.assertEqual(r['call'],'similarity_comparator_only')

if __name__=='__main__':unittest.main()

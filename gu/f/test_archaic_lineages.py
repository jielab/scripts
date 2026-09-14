"""Regression checks for reference-specific calls and two-lineage reporting."""
import csv
import gzip
import json
import io
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from contextlib import ExitStack
import numpy as np
import phyml_gwas as phy
import phyml_report as report
import ibdmix_workflow as ibd
import introgression_density as density
from phyml_gwas_input import write


class LineageTests(unittest.TestCase):
    def test_interrupted_mask_download_resumes_and_validates_gzip(self):
        payload=gzip.compress(b'archaic reference quality mask\n'*100)
        requested=[]
        def response(request,**kwargs):
            lo,hi=map(int,request.get_header('Range').split('=')[1].split('-'))
            hi=min(hi,len(payload)-1);requested.append(lo)
            body=payload[lo:hi+1]
            if len(requested)==1:body=body[:4]  # First connection closes early.
            result=io.BytesIO(body);result.status=206
            result.headers={'Content-Length':str(hi-lo+1),'Content-Range':f'bytes {lo}-{hi}/{len(payload)}','ETag':'"v1"'}
            return result
        with tempfile.TemporaryDirectory() as tmp, patch.object(ibd.urllib.request,'urlopen',response), patch.object(ibd.time,'sleep'):
            target=Path(tmp)/'mask.gz'
            ibd.download('https://example.invalid/mask.gz',target,chunk_bytes=16)
            self.assertEqual(requested[:2],[0,4])
            self.assertEqual(target.read_bytes(),payload)
            self.assertEqual(gzip.decompress(target.read_bytes()),b'archaic reference quality mask\n'*100)
            self.assertFalse(list(Path(tmp).glob('*.part*')))

    def test_five_reference_masks_intersect_common_exclusions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            beds={}
            for i,ref in enumerate(phy.REFS):
                beds[ref]=root/(ref+'.bed');beds[ref].write_text(f'1\t{100+i*20}\t800\n')
            strict=root/'strict.bed';strict.write_text('1\t0\t900\n')
            dups=root/'dups.txt';dups.write_text('0\tchr1\t300\t310\n')
            fasta=root/'reference.fa';fasta.write_text('dummy')
            axt=root/'alignment.axt';axt.write_text('dummy')
            variants=root/'j.cell.2020.01.012-workflow/abridged_variants.gz'
            variants.parent.mkdir();variants.write_text('dummy')
            modern=root/'modern.vcf'
            modern.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n'
                              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                              '1\t401\t.\tAA\tA\t.\tPASS\t.\n')
            def download(url,path):
                if 'strict_mask' in url:return strict
                if 'genomicSuperDups' in url:return dups
                return beds['Altai' if 'AltaiNea' in url else 'Denisova']
            cpg=np.zeros(1000,dtype=bool);cpg[500]=True
            with ExitStack() as stack:
                stack.enter_context(patch.dict(ibd.CHROM_LENGTHS,{'37':{'1':1000}},clear=True))
                stack.enter_context(patch.object(ibd,'download',download))
                stack.enter_context(patch.object(ibd,'published_reference_mask',lambda ref,*args:beds[ref]))
                stack.enter_context(patch.object(ibd,'download_axt',lambda *args:axt))
                stack.enter_context(patch.object(ibd,'cached_cpg_sites',lambda *args:cpg.copy()))
                ibd.prepare_masks(root/'resources','1',modern,fasta,root,root/'out',axt_root=root/'axt',refs=list(phy.REFS))
            for i,ref in enumerate(phy.REFS):
                excluded=np.zeros(1000,dtype=bool)
                for lo,hi in ibd.bed_intervals(root/'out'/f'{ref}.exclude.bed','1',1000):excluded[lo:hi]=True
                self.assertEqual(int((~excluded).sum()),676-i*20)
                self.assertTrue(excluded[305] and excluded[400] and excluded[500])
            self.assertEqual(set(json.loads((root/'out/manifest.json').read_text())['callable_bp']),set(phy.REFS))

    def test_tree_test_does_not_conflate_archaic_lineages(self):
        tree='((H1:1,Denisova:1,Denisova25:1)95:1,(Altai:1,Chagyr:1,Vindija:1):1,H2:1,Ancestral:1);'
        self.assertIsNone(phy.risk_clade(tree,{'H1'},{'H1','H2'}))
        match=phy.risk_clade(tree,{'H1'},{'H1','H2'},phy.LINEAGE_REFS['Denisovan'])
        self.assertEqual(match['bootstrap'],95)
        mixed=tree.replace('Denisova25:1)95','Denisova25:1,Altai:1)95').replace('(Altai:1,Chagyr:1','(Chagyr:1')
        self.assertIsNone(phy.risk_clade(mixed,{'H1'},{'H1','H2'},phy.LINEAGE_REFS['Denisovan']))

    def test_both_reports_and_copy_ids_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            summaries=[dict(locus_id='L',dataset_id='1kg',genome_build='GRCh37',chr='1',
                            core_start=0,core_end=100,core_kb=.1,lineage=lin,
                            n_candidate_haplotypes=1,n_candidate_copies=1) for lin in phy.LINEAGE_REFS]
            write(root/'gwas_loci.tsv',summaries)
            write(root/'gwas_haplotypes.tsv',[dict(locus_id='L',lineage=lin,hap_id='H1',role='risk') for lin in phy.LINEAGE_REFS])
            write(root/'gwas_copies.tsv',[dict(locus_id='L',lineage=lin,hap_id='H1',sample_id='S',haplotype=1) for lin in phy.LINEAGE_REFS])
            ss,dd,cc=report.gwas_rows(root/'gwas_loci.tsv')
            self.assertEqual(len(ss),2); self.assertEqual(len(dd),2); self.assertEqual(len(cc),2)
            self.assertEqual({s['lineage'] for s in ss},set(phy.LINEAGE_REFS))
            self.assertEqual(len({c['candidate_id'] for c in cc}),2)

    def test_control_only_run_is_not_a_negative_denisovan_test(self):
        with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as con:
            root=Path(tmp);(root/'samples/C1').mkdir(parents=True)
            (root/'samples/C1/ALL.txt').write_text('S\n')
            meta=root/'run.meta.tsv';meta.write_text('refs\tAltai Denisova\nbackground_filter\t1\n')
            con.executescript('CREATE TABLE method_runs(dataset_id,genome_build,chr,method,status,evidence_eligible,raw_file,availability_note);'
                              'CREATE TABLE segments(dataset_id,genome_build,method,source_class,sample_id,chr,start,end,haplotype);')
            con.execute('INSERT INTO method_runs VALUES (?,?,?,?,?,?,?,?)',('1kg','GRCh37','1','ibdmix','complete',1,str(meta),'complete'))
            copies=[dict(dataset_id='1kg',genome_build='GRCh37',chr='1',candidate_start=100,candidate_end=200,
                         sample_id='S',haplotype=1,lineage=lin) for lin in phy.LINEAGE_REFS]
            old,_=report.validate(copies,con)
            self.assertEqual([r['comparison_available'] for r in old if r['method']=='ibdmix'],[1,0])
            meta.write_text(meta.read_text()+'export_denisovan\t1\n')
            new,_=report.validate(copies,con)
            self.assertEqual([r['comparison_available'] for r in new if r['method']=='ibdmix'],[1,1])
            self.assertEqual([r['overlap_pass'] for r in new if r['method']=='ibdmix'],[0,0])

    def test_lineage_results_use_separate_bootstrap_and_reference_matches(self):
        nw='((H1:1,Denisova:1,Denisova25:1)95:1,(Altai:1,Chagyr:1,Vindija:1):1,H2:1,Ancestral:1);'
        ident=dict(lineage='Neanderthal',tree_pass=0,tree_bootstrap=None,status='tree_not_supported',call='tree_not_supported')
        haps=[dict(hap_id='H1',role='risk',seq='AA'),dict(hap_id='H2',role='nonrisk',seq='CC')]
        arch={ref:'AA' if ref in ibd.DENISOVANS else 'CC' for ref in phy.REFS}
        details=[dict(hap_id='H1',role='risk',call='risk_haplotype')]
        ss,tt,dd=phy.lineage_results(ident,dict(tree_status='complete',tree_newick=nw),details,haps,arch)
        self.assertEqual([s['tree_pass'] for s in ss],[0,1])
        self.assertEqual([d['n_match'] for d in dd],[0,2])
        self.assertEqual(tt[1]['candidate_clade_archaic_tips'],2)
        self.assertIsNone(ss[1]['ils_probability'])

    def test_background_filter_preserves_independent_denisovan_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); pops=[]
            for pop in ['YRI','CHB']:
                roster=root/(pop+'.txt');roster.write_text(pop+'1\n')
                pops.append(dict(population=pop,super_population='AFR' if pop=='YRI' else 'EAS',sample_file=str(roster)))
                for ref in phy.REFS:
                    with gzip.open(root/f'{ref}.{pop}.raw.txt.gz','wt') as f:
                        f.write('ID\tchrom\tstart\tend\tLOD\n'+f'{pop}1\t1\t1\t101\t5\n')
            write(root/'pop.tsv',pops)
            out=root/'out.gz'
            ibd.finalize(root,root/'pop.tsv',phy.REFS,out,'1','b37','genome',min_bp=10,export_denisovan=True)
            with gzip.open(out,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),4)
            self.assertEqual({r['anc'] for r in rows},{'Denisova','Denisova25'})
            self.assertEqual({r['sample_set'] for r in rows},{'YRI','CHB'})
            # Original control-only profile still exports no Denisovan signal.
            ibd.finalize(root,root/'pop.tsv',['Altai','Denisova'],out,'1','b37','genome',min_bp=10)
            with gzip.open(out,'rt') as f:self.assertEqual(list(csv.DictReader(f,delimiter='\t')),[])

    def test_density_unions_denisovan_references_and_keeps_untested_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=root/'run';(run/'samples').mkdir(parents=True)
            (run/'samples/ALL.txt').write_text('S\n')
            meta=run/'run.meta.tsv'
            meta.write_text('refs\tAltai Denisova Denisova25\nbackground_filter\t1\nexport_denisovan\t1\nmodern_vcf\tinput.vcf\n')
            db=root/'data.sqlite'
            with sqlite3.connect(db) as con:
                con.executescript('CREATE TABLE method_runs(dataset_id,genome_build,chr,method,status,evidence_eligible,raw_file,availability_note);'
                                  'CREATE TABLE sample_populations(dataset_id,sample_id,population,super_population);'
                                  'CREATE TABLE segments(dataset_id,genome_build,method,source_class,sample_id,chr,start,end);')
                con.execute('INSERT INTO method_runs VALUES (?,?,?,?,?,?,?,?)',('1kg','GRCh37','1','ibdmix','complete',1,str(meta),'audited IBDmix test; profile=multi_reference'))
                con.execute("INSERT INTO sample_populations VALUES ('1kg','S','CHB','EAS')")
                for start,end in [(100,300),(200,400)]:
                    con.execute("INSERT INTO segments VALUES ('1kg','GRCh37','ibdmix','Denisovan','S','1',?,?)",(start,end))
            with patch.dict(density.CHROM_LENGTHS,{'37':{'1':1000,'2':1000}},clear=True):
                density.prepare(db,root/'cache',bin_bp=1000)
            pointer=next((root/'cache').glob('*/current.tsv'))
            with pointer.open() as f:version=next(csv.DictReader(f,delimiter='\t'))['directory']
            with gzip.open(pointer.parent/version/'denisovan_matrix.tsv.gz','rt') as f:z=np.loadtxt(f)
            self.assertEqual(z[0],15)  # 300 union bp / 2000 diploid bp, not 400.
            self.assertTrue(np.isnan(z[1]))


if __name__=='__main__':unittest.main()

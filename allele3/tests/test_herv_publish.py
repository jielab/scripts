import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
from herv_publish import publish, hydrate_runtime, PRODUCTS
from herv_export import FIELDS, sha
from herv_cohort import parser, requested_populations
from herv_archive import validate
from herv_layout import product_path, ensure_layout


class PublicationTest(unittest.TestCase):
    def test_highlights_reads_results_and_writes_nested_reports(self):
        from herv_cohort import main
        import research_tests
        from enrichment_regions import RegionHighlighter
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'step3/herv_highlights';ensure_layout(out)
            source=Path(tmp)/'herv.OLD';ensure_layout(source)
            scan=Path(tmp)/'scan';(scan/'qc').mkdir(parents=True)
            (scan/'qc/manifest.json').write_text(json.dumps(dict(status='complete',target='gnomad',build='GRCh38')))
            for name in ('reference_herv.bed','stable_loci.bed'):
                product_path(source,name).write_text('chr1\t10\t20\tERV1|locus1\n')
            def enrichment(command):
                a=research_tests.parser().parse_args(command)
                self.assertEqual(Path(a.herv_bed),source/'results/reference_herv.bed')
                self.assertEqual(Path(a.stable_herv_bed),source/'results/stable_loci.bed')
                _,meta=research_tests.prepare_herv(a,Path(a.out))
                regions=RegionHighlighter(a,Path(a.out),meta)
                self.assertEqual(regions.stable_path,source/'results/stable_loci.bed')
                destination=Path(a.out)
                for name in ('results/highlights.tsv.gz','results/highlights.bed.gz','results/highlights_supported.tsv.gz',
                             'results/regions.tsv.gz','reports/highlights.html','figures/highlights_genome.png','figures/highlights_ranked.png'):
                    p=destination/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture')
                (destination/'qc').mkdir(parents=True,exist_ok=True)
                (destination/'qc/region_highlights.json').write_text(json.dumps(dict(status='complete',counts=dict(highlight_regions=0))))
                return dict(status='complete')
            before={str(p.relative_to(source)):(p.read_bytes(),p.stat().st_mtime_ns) for p in source.rglob('*') if p.is_file()}
            with patch.object(research_tests,'main',side_effect=enrichment):
                result=main(['highlights','--out',str(out),'--herv-input-dir',str(source),'--scan-dir',str(scan),'--target','gnomad','--no-context'])
            self.assertEqual(before,{str(p.relative_to(source)):(p.read_bytes(),p.stat().st_mtime_ns) for p in source.rglob('*') if p.is_file()})
            self.assertEqual(result['status'],'complete')
            self.assertTrue((out/'qc/enrichment_status.json').is_file())
            html=(out/'reports/herv_highlights.html').read_text()
            self.assertIn('../enrichment/gnomad/reports/highlights.html',html)
            self.assertTrue(all(p.is_dir() for p in out.iterdir()))

    def test_roundtrip_missingness_and_failed_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'runtime';out.mkdir();archive=Path(tmp)/'archive'
            ensure_layout(out)
            for name in ('paired.png','paired.pdf','only.pdf'):(out/'figures'/name).write_text('fixture')
            args=parser().parse_args(['run','--out',str(out),'--dir-out-herv_scan',str(archive)])
            for name in PRODUCTS:product_path(out,name).write_text('fixture\n')
            samples=[dict(sample='S1',pop='GBR',super_pop='EUR')]
            identity=dict(build='GRCh38',seed=args.seed,per_pop=args.per_pop,
                populations=requested_populations(args.populations),sample_ids=['S1'],pfile_dir=str(args.pfile_dir.resolve()))
            product_path(out,'plan.json').write_text(json.dumps(dict(samples=samples,identity=identity)))
            product_path(out,'samples.tsv').write_text('sample\tpop\tsuper_pop\nS1\tGBR\tEUR\n')
            settings={k:getattr(args,k) for k in ('min_callable','min_region_coverage','stable_retention',
                'stable_min_sample_fraction','stable_min_population_fraction','bootstrap','seed')}
            sequence=dict(chunk_bp=args.chunk_bp,halo=args.halo)
            product_path(out,'summary.json').write_text(json.dumps(dict(status='complete',missing_samples=[],
                planned_samples=1,completed_samples=1,settings=settings,sequence_annotation_settings=sequence)))
            chroms=['chr'+str(i) for i in range(1,23)];files=[]
            for chrom in chroms:
                for kind in ('segments','uncallable'):
                    p=out/'tmp/export_tables'/kind/(chrom+'.tsv.gz');p.parent.mkdir(parents=True,exist_ok=True)
                    with gzip.open(p,'wt') as f:
                        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(FIELDS)
                        writer.writerow([chrom,10,20,'S1','hap1','ERV1' if kind=='segments' else 'uncallable','GBR','EUR'])
                    files.append(dict(chrom=chrom,kind=kind,path=str(p.relative_to(out)),rows=1,sha256=sha(p)))
            product_path(out,'export_manifest.json').write_text(json.dumps(dict(status='complete',schema='herv-chromosome-export-1',
                build='GRCh38',samples=1,haplotypes=2,chromosomes=chroms,files=files,sequence_annotation_settings=sequence)))
            publish(out,archive)
            self.assertFalse((out/'figures/paired.pdf').exists())
            self.assertTrue((out/'figures/only.pdf').exists())
            self.assertEqual(validate(args)['chromosome_tables'],22)
            with gzip.open(archive/'chr1.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual([r['record_type'] for r in rows],['segments','uncallable'])
            self.assertFalse((archive/'raw').exists());self.assertFalse((archive/'tmp').exists())
            self.assertEqual(len(list(archive.iterdir())),22)
            self.assertFalse((archive/'EXPORT_README.txt').exists())
            product_path(out,'reference_herv.bed').unlink();hydrate_runtime(out,archive)
            self.assertTrue((out/'logs').is_dir());self.assertTrue(product_path(out,'reference_herv.bed').is_file())
            product_path(out,'reference_herv.bed').write_text('stale')
            hydrate_runtime(out,archive,refresh=True)
            self.assertEqual(sha(product_path(out,'reference_herv.bed')),sha(out/'qc/archive/reference_herv.bed'))
            from herv_cohort import write_report, update_highlight_reports
            write_report(args,dict(n_samples=0,n_populations=0),[],[],[],[],[],{})
            update_highlight_reports(args,dict(targets=[]))
            html=(out/'reports/herv_scan.html').read_text()
            self.assertIn('../results/family_stability.tsv',html)
            self.assertIn('../figures/family_stability.png',html)
            self.assertTrue(all(p.is_dir() for p in out.iterdir()))
            original=sha(out/'qc/archive/SHA256SUMS')
            (out/files[0]['path']).write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'checksum'):publish(out,archive)
            self.assertEqual(sha(out/'qc/archive/SHA256SUMS'),original)
            with self.assertRaises(ValueError):publish(archive,archive)


if __name__=='__main__':unittest.main()

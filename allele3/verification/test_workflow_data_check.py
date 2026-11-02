import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'f'))
from data_sources import load_config
from workflow_data_check import actual_sources


class WorkflowDataCheck(unittest.TestCase):
    def test_conflict_prevents_execution_and_matching_declaration_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); prefix=root/'panel'; marker=root/'ran.txt'
            prefix.with_suffix('.pgen').write_bytes(b'not-decoded-by-metadata-audit')
            prefix.with_suffix('.psam').write_text('#IID\nS1\n')
            config=root/'sources.tsv'
            config.write_text('id\ttype\tpath\tbuild\tevidence\trole\n'
                              f'ancient_pgen\tpgen\t{prefix}\t37\timputed_genotypes\tfixture\n')
            entry=root/'entry.py'; entry.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("ran")\n')
            command=[sys.executable,str(ROOT/'f/workflow_data_check.py'),'--config',str(config),
                     '--workflow','pca3d','--',sys.executable,str(entry),'run','--pgen-prefix',str(prefix)]
            for reference, expected in [('GRCh38',2),('GRCh37',0)]:
                with self.subTest(reference=reference):
                    prefix.with_suffix('.pvar').write_text('##reference='+reference+'\n#CHROM\tPOS\tID\tREF\tALT\n1\t10\tx\tA\tG\n')
                    run=subprocess.run(command,text=True,capture_output=True,timeout=30)
                    self.assertEqual(run.returncode,expected,run.stdout+run.stderr)
                    self.assertEqual(marker.exists(),expected==0)
                    self.assertIn('[build-check]',run.stdout+run.stderr)

    def test_actual_pipeline_paths_override_registry(self):
        registry=load_config(ROOT/'config/data_sources.tsv')
        entries, _, _ = actual_sources('genotypes', ['genotypes','--stages','primary,genotypes','--kg-vcf-dir','/different/individual-vcfs',
            '--scan-dir','/different/scan'], registry)
        byid={e['id']:e for e in entries}
        self.assertEqual(byid['1kg_vcf']['path'],'/different/individual-vcfs')
        self.assertEqual(byid['1kg_vcf']['type'],'vcf_directory')
        self.assertEqual(byid['gnomad_scan']['path'],'/different/scan')
        self.assertNotIn('1kg_sites',byid)

    def test_unread_sources_are_not_build_gates(self):
        registry=load_config(ROOT/'config/data_sources.tsv')
        entries, _, _ = actual_sources('selection', ['selection'], registry)
        self.assertEqual(entries, [])
        entries, _, _ = actual_sources('herv', ['--input-format','vcf','--pfile-dir','/unused/pgen'], registry)
        self.assertEqual(entries, [])
        entries, _, _ = actual_sources('pca3d', ['project-plan'], registry)
        self.assertEqual(entries, [])

    def test_gu_module_is_not_silently_changed(self):
        registry=load_config(ROOT/'config/data_sources.tsv')
        entries, _, _ = actual_sources('research-archaic-alleles',
            ['archaic-alleles','--gu-dataset','/actual/ukb012','--module','ibdmix-GPU','--modern-build','37'],registry)
        self.assertEqual(entries[0]['module'],'ibdmix-GPU')
        self.assertEqual(entries[0]['path'],'/actual/ukb012')

    def test_all_build_messages_survive_compact_consoles(self):
        for file, name in [('pipeline_console.py','PipelineScreenLog'),('herv_console.py','HervScreenLog')]:
            with self.subTest(console=file):
                spec=importlib.util.spec_from_file_location('review_console',ROOT/'runtime'/file)
                module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
                log=getattr(module,name)(300); output=io.StringIO()
                with contextlib.redirect_stdout(output):
                    for i in range(3):
                        log.feed(f'[dataset] source{i}\n[build-check] source{i}: GRCh37\n'.encode())
                    log.finish()
                self.assertEqual(output.getvalue().count('[dataset]'),3)
                self.assertEqual(output.getvalue().count('[build-check]'),3)


if __name__=='__main__': unittest.main()

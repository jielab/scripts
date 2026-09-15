"""Completed trees must never be silently replaced during resume."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import phyml_run as runner


class TreeResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.phy = Path(self.tmp.name)/'haplotypes.phy'
        self.text = '2 1\nA A\nB T\n'
        self.phy.write_text(self.text)
        for suffix, text in {'_phyml_tree.txt':'(A,B);', '_phyml_stats.txt':'stats',
                             '_phyml_boot_trees.txt':'(A,B);\n'*100,
                             '_phyml_boot_stats.txt':'stats',
                             '.phyml.log':'100/100\nPrinting the most likely tree\n. Time used 0h0m1s'}.items():
            Path(str(self.phy)+suffix).write_text(text)
        self.binary = Path(self.tmp.name)/'phyml'
        self.binary.write_text('fake executable')
        runner.seal(self.phy, runner.request(self.phy,100,self.binary))

    def snapshot(self):
        return {p.name:(p.read_bytes(),p.stat().st_mtime_ns)
                for p in self.phy.parent.glob('haplotypes*') if not p.name.endswith('.lock')}

    def test_identical_and_whitespace_only_inputs_preserve_all_files(self):
        before = self.snapshot()
        runner.prepare_input(self.phy,self.text)
        runner.prepare_input(self.phy,'2  1\nA     A\nB     T\n')
        self.assertEqual(before,self.snapshot())

    def test_changed_alignment_preserves_completed_files(self):
        before = self.snapshot()
        with self.assertRaisesRegex(RuntimeError,'alignment changed'):
            runner.prepare_input(self.phy,'2 1\nA G\nB T\n')
        self.assertEqual(before,self.snapshot())

    def test_historical_complete_tree_without_receipt_is_protected(self):
        Path(str(self.phy)+'.phyml.complete.json').unlink()
        before = self.snapshot()
        with self.assertRaises(RuntimeError):
            runner.prepare_input(self.phy,'2 1\nA G\nB T\n')
        self.assertEqual(before,self.snapshot())

    def test_missing_completed_product_does_not_authorize_replacement(self):
        Path(str(self.phy)+'_phyml_tree.txt').unlink()
        before = self.snapshot()
        with self.assertRaises(RuntimeError):
            runner.prepare_input(self.phy,'2 1\nA G\nB T\n')
        self.assertEqual(before,self.snapshot())

    def test_unfinished_tree_without_receipt_can_resume(self):
        Path(str(self.phy)+'.phyml.complete.json').unlink()
        Path(str(self.phy)+'.phyml.log').write_text('interrupted at bootstrap 4/100')
        runner.prepare_input(self.phy,'2 1\nA G\nB T\n')
        self.assertEqual(self.phy.read_text(),'2 1\nA G\nB T\n')

    def test_explicit_replace_allows_changed_input(self):
        runner.prepare_input(self.phy,'2 1\nA G\nB T\n',replace=True)
        self.assertEqual(self.phy.read_text(),'2 1\nA G\nB T\n')
        self.assertFalse(Path(str(self.phy)+'.phyml.complete.json').exists())

    def test_new_skip_preserves_existing_locus_tables(self):
        import phyml_gwas as gwas
        from types import SimpleNamespace
        out = self.phy.parent/'unit'
        loc = out/'loci'
        loc.mkdir(parents=True)
        (out/'final').mkdir()
        for path in list(self.phy.parent.glob('haplotypes*')):
            path.rename(loc/path.name)
        (loc/'sites.tsv').write_text('original sites')
        (out/'final/gwas_loci.tsv').write_text('original result')
        before = {str(p):p.read_bytes() for p in out.rglob('*') if p.is_file()}
        row = dict(locus_id='unit',chr='1',lead_pos='10',index_snp='rs1',
                   source_build='37',source_pos='10',p_j='1e-9',beta_j='1',
                   effect_allele='A',search_start='1',search_end='20')
        args = SimpleNamespace(out=out,dataset='1kg',vcf_dir=out)
        with patch.dict(os.environ,PHYML_REPLACE='FALSE'), \
             patch.object(gwas,'vcf_path',side_effect=gwas.SkipLocus('changed inputs','missing_lead')):
            with self.assertRaisesRegex(RuntimeError,'preserved'):
                gwas.run_locus(args,row)
        self.assertEqual(before,{str(p):p.read_bytes() for p in out.rglob('*') if p.is_file()})

    def run_main(self, force=False, verify=False):
        argv=['phyml_run.py','--phy',str(self.phy),'--scope','test']
        if verify: argv.append('--verify-only')
        with patch.object(sys,'argv',argv), patch.dict(os.environ,PHYML_REPLACE='TRUE' if force else 'FALSE'), \
             patch.object(runner.shutil,'which',side_effect=lambda name:str(self.binary) if name=='phyml' else None):
            return runner.main()

    def test_valid_tree_skips_execution(self):
        with patch.object(runner,'run_attempt') as run:
            self.assertEqual(self.run_main(),0)
            run.assert_not_called()

    def test_fresh_run_then_resume_without_recomputing(self):
        products = {suffix:Path(str(self.phy)+suffix).read_text()
                    for suffix in runner.SUFFIXES
                    if Path(str(self.phy)+suffix).is_file()
                    and suffix != '.phyml.complete.json'}
        for path in self.phy.parent.glob('haplotypes*'):
            path.unlink()
        runner.prepare_input(self.phy, self.text)
        def complete(phy, command, deadline):
            for suffix, text in products.items():
                Path(str(phy)+suffix).write_text(text)
            return 0
        with patch.object(runner, 'run_attempt', side_effect=complete) as run:
            self.assertEqual(self.run_main(), 0)
            self.assertIsNone(runner.completion_error(self.phy, 100))
            receipt = Path(str(self.phy)+'.phyml.complete.json').read_bytes()
            runner.prepare_input(self.phy, self.text)
            self.assertEqual(self.run_main(), 0)
            run.assert_called_once()
            self.assertEqual(receipt, Path(str(self.phy)+'.phyml.complete.json').read_bytes())

    def test_changed_binary_does_not_destroy_complete_tree(self):
        self.binary.write_text('changed binary')
        before=self.snapshot()
        with patch.object(runner,'run_attempt') as run:
            with self.assertRaisesRegex(RuntimeError,'preserved'):
                self.run_main()
            run.assert_not_called()
        self.assertEqual(before,self.snapshot())

    def test_force_runs_but_followup_verification_never_replaces(self):
        # Stop before executing any binary, proving force bypasses valid reuse.
        with patch.object(runner,'run_attempt',side_effect=RuntimeError('attempt reached')) as run:
            with self.assertRaisesRegex(RuntimeError,'attempt reached'):
                self.run_main(force=True)
            run.assert_called_once()

    def test_verify_only_ignores_force(self):
        before=self.snapshot()
        with patch.object(runner,'run_attempt') as run:
            self.assertEqual(self.run_main(force=True,verify=True),0)
            run.assert_not_called()
        self.assertEqual(before,self.snapshot())


if __name__=='__main__':
    unittest.main()

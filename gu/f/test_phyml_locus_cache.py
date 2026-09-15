"""Whole-locus resume regressions; no reference tools or real analyses run."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import phyml_locus_cache as cache


class LocusCacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root/'unit'
        (self.out/'final').mkdir(parents=True)
        (self.out/'loci').mkdir()
        self.arch = self.root/'archaic'
        self.arch.mkdir()
        (self.arch/'chr1.vcf.gz').write_text('reference')
        (self.root/'chr1.pgen').write_text('genotypes')
        (self.root/'samples.txt').write_text('sample\tpop\ns1\tEUR\n')
        (self.root/'leads.tsv').write_text('locus_id\tchr\nlead1\t1\n')
        (self.out/'unit.bed').write_text('1\t10\t20\tlead1\n')
        self.cmd = self.out/'unit.cmd'
        self.cmd.write_text(f'export GU_PHYML_LEAD_TABLE={self.root}/leads.tsv\n'
                            f'exec /unused/gu.sh phyml --loci {self.out}/unit.bed --target-dir {self.root}/chr '
                            '--replace-phyml FALSE --plot-phy TRUE --memory-cap 24G\n')
        for name in ('gwas_haplotypes.tsv', 'gwas_copies.tsv', 'loci.tsv', 'evidence_trees.tsv',
                     'haplotypes.tsv', 'haplotype_samples.tsv', 'skipped_loci.tsv'):
            (self.out/'final'/name).write_text('status\n')
        (self.out/'final/gwas_lead.tsv').write_text((self.root/'leads.tsv').read_text())
        (self.out/'final/gwas_loci.tsv').write_text('status\ninsufficient_high_LD_markers\ninsufficient_high_LD_markers\n')
        (self.out/'final/trees.tsv').write_text('tree_status\nnot_run\n')
        (self.out/'final/gwas_parameters.json').write_text(json.dumps({'workflow':'gwas_lead_ld_core_archaic5_v2'}))
        self.log = self.cmd.with_suffix('.log')
        self.log.write_text(f'archaic reference root={self.arch}\nsample metadata={self.root}/samples.txt\n[GU RUN] method=phyml analysis_unit=unit status=complete\n')

    def run_cache(self, mode):
        return cache.process(mode, self.cmd, str(self.arch))

    def test_historical_skip_is_reused(self):
        self.assertTrue(self.run_cache('adopt'))
        self.assertTrue(self.run_cache('check'))

    def test_successful_non_supported_result_is_reused(self):
        # The cache must not restrict reuse to a positive scientific result.
        (self.out/'final/gwas_loci.tsv').write_text('status\ntree_not_supported\ntree_not_supported\n')
        (self.out/'final/trees.tsv').write_text('tree_status\ncomplete\n')
        phy = self.out/'loci/haplotypes.phy'
        phy.write_text('2 1\nA A\nB T\n')
        for suffix, text in {'_phyml_tree.txt':'(A,B);', '_phyml_stats.txt':'stats',
                             '_phyml_boot_trees.txt':'(A,B);\n'*100, '_phyml_boot_stats.txt':'stats',
                             '.phyml.log':'100/100\nPrinting the most likely tree\n. Time used 0h0m1s'}.items():
            Path(str(phy)+suffix).write_text(text)
        for lineage in ('Neanderthal', 'Denisovan'):
            for suffix in ('.png', '.pdf', '.full.png', '.full.pdf'):
                Path(str(phy)+f'_phyml_tree.{lineage}.panelB{suffix}').write_text('plot')
        self.assertTrue(self.run_cache('seal'))
        self.assertTrue(self.run_cache('check'))

    def test_code_only_change_does_not_schedule_completed_analysis(self):
        self.run_cache('seal')
        receipt = self.out/'.phyml.locus.complete.json'
        data = json.loads(receipt.read_text())
        data['request']['code'] = {'phyml_gwas.py':'historical-version'}
        receipt.write_text(json.dumps(data))
        self.assertTrue(self.run_cache('check'))

    def test_check_adopts_unsealed_success_without_overwriting_log(self):
        before = self.log.read_bytes()
        self.assertTrue(self.run_cache('check'))
        self.assertEqual(before, self.log.read_bytes())

    def test_runtime_settings_do_not_invalidate(self):
        self.run_cache('seal')
        self.cmd.write_text(self.cmd.read_text().replace('24G', '16G'))
        with patch.dict(os.environ, PHYML_TREE_CPUS='8', PHYML_TREE_TIMEOUT='123'):
            self.assertTrue(self.run_cache('check'))

    def test_force_or_changed_analysis_invalidates(self):
        self.run_cache('seal')
        original = self.cmd.read_text()
        for text in (original.replace('--replace-phyml FALSE', '--replace-phyml TRUE'),
                     original.replace('--plot-phy TRUE', '--plot-phy FALSE')):
            self.cmd.write_text(text)
            self.assertFalse(self.run_cache('check'))

    def test_changed_source_or_lead_invalidates(self):
        self.run_cache('seal')
        (self.root/'chr1.pgen').write_text('different genotypes')
        self.assertFalse(self.run_cache('check'))
        (self.root/'leads.tsv').write_text('locus_id\tchr\nlead1\t2\n')
        self.assertFalse(self.run_cache('check'))

    def test_missing_or_modified_output_invalidates(self):
        self.run_cache('seal')
        p = self.out/'final/gwas_copies.tsv'
        p.write_text('changed\n')
        self.assertFalse(self.run_cache('check'))
        p.unlink()
        with self.assertRaises(ValueError):
            self.run_cache('check')

    def test_failed_and_interrupted_runs_not_adopted(self):
        self.log.write_text('partial run\n')
        with self.assertRaises(ValueError):
            self.run_cache('adopt')
        self.log.write_text(f'archaic reference root={self.arch}\nanalysis_unit=unit status=complete\n')
        (self.out/'final/gwas_loci.tsv').write_text('status\ntree_failed\ntree_failed\n')
        with self.assertRaises(ValueError):
            self.run_cache('adopt')

    def test_scheduler_skips_without_launch_or_log_overwrite(self):
        self.run_cache('seal')
        before = self.log.read_bytes()
        gu = Path(__file__).resolve().parent.parent/'gu.sh'
        source = gu.read_text()
        functions = source[source.index('gu_completed_analysis_cmd_output(){'):source.index('gu_run_analysis_cmds_local(){')]
        env = dict(os.environ, METHOD='phyml', REPLACE_PHYML_INPUT='FALSE',
                   F=str(gu.parent/'f'), GU_ARCHAIC_ROOT=str(self.arch))
        result = subprocess.run(['bash', '-c', 'set -euo pipefail\n'+functions+'\ngu_run_one_analysis_cmd "$1"',
                                 'test', str(self.cmd)], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('SKIP unit=unit reason=output_complete', result.stderr)
        self.assertEqual(before, self.log.read_bytes())

    def test_partition_finds_completed_locus_after_pending_locus(self):
        self.run_cache('seal')
        pending_cmd = self.root/'unfinished'/'unfinished.cmd'
        pending_cmd.parent.mkdir()
        pending_cmd.write_text('unused')
        listing = self.root/'commands.list'
        listing.write_text(f'{pending_cmd}\n{self.cmd}\n')
        pending = self.root/'pending.list'
        before = self.log.read_bytes()
        result = subprocess.run(['python3', str(Path(cache.__file__)), 'partition', str(listing),
                                 '--pending', str(pending), '--archaic-root', str(self.arch)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(pending.read_text(), f'{pending_cmd}\n')
        self.assertIn('RESUME total=2 reused=1 pending=1', result.stdout)
        self.assertEqual(before, self.log.read_bytes())


if __name__ == '__main__':
    unittest.main()

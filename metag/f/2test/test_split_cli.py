import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_cli as cli
import hpc_query


class SplitCliTests(unittest.TestCase):
    def test_submit_values(self):
        for command in ('split',):
            for suffix, value in (([], False),
                                  (['--submit', 'TRUE'], True), (['--submit', 'FALSE'], False)):
                args = cli.parser().parse_args([command, '--engine', 'kraken', '--ready-only', *suffix, '--queue', 'spec'])
                self.assertIs(args.submit, value)
                self.assertEqual(args.queue, 'spec')
        self.assertFalse(cli.parser().parse_args(['library', '--submit', 'FALSE']).submit)

    def test_old_command_and_bare_submit_are_rejected(self):
        for args in (['classify'], ['split', '--submit']):
            with self.assertRaises(SystemExit):
                cli.parser().parse_args(args)

    def test_queue_applies_after_engine_defaults_to_pilot_and_batches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = root / 'sample.fastq'
            reads.write_text('@read\nA\n+\nI\n')
            cfg = dict(work_dir=str(root / 'work'), output_dir=str(root / 'out'),
                       code_root=str(cli.CODE_ROOT), classifier_engine='kraken',
                       sample_batch_size=1, max_array_size=100, classify_concurrency=2,
                       classify_workers=1, classify_queue='smp', pilot_queue='ser',
                       classify_memory_gib=32, pilot_memory_gib=32,
                       lsf_memory_scope='site_host', lsf_resource_reserve_per_task=False)
            samples = [dict(sample_id='data/sample', dataset='data', sample='sample', r1=str(reads), r2=None)]
            with patch.object(cli, 'engine_resources', return_value=dict(cfg)), \
                 patch.object(cli, 'engine_snapshot', return_value={'identity': 'test'}), \
                 patch.object(cli, '_database_states', return_value={}), \
                 patch.object(cli, 'effective_workers'), \
                 patch.object(hpc_query, 'discover_samples', return_value=samples):
                path = cli.generate_split(cfg, ready_only=True, queue='spec')
            plan = json.loads(path.read_text())
            self.assertEqual(len(plan['commands']), 2)
            for command in plan['commands']:
                text = Path(command['path']).read_text()
                self.assertIn('#BSUB -q spec\n', text)
                self.assertIn('rusage[mem=32GB]', text)

    def test_library_stage_queues_preserve_memory_and_submission_mode(self):
        for conversion_only in (False, True):
            cfg = dict(prepare_queue='ser', prepare_cores=8, build_queue='smp',
                       convert_queue='ser', build_memory_gib=4096, convert_memory_gib=32)
            with patch.object(cli, 'load_config', return_value=cfg), \
                 patch.object(cli, 'generate_library', return_value=Path('/test/jobs.json')) as library, \
                 patch.object(cli, 'generate_conversion', return_value=Path('/test/jobs.json')) as conversion, \
                 patch.object(cli, 'submit_plan') as submit:
                cli.main(['library', '--kunpeng' if conversion_only else '--prepared',
                          '--build-queue', 'smp', '--convert-queue', 'spec', '--submit', 'FALSE'])
                chosen = conversion if conversion_only else library
                chosen.assert_called_once()
                self.assertEqual(cfg['build_queue'], 'smp')
                self.assertEqual(cfg['convert_queue'], 'spec')
                self.assertEqual(cfg['prepare_queue'], 'ser')
                self.assertEqual(cfg['build_memory_gib'], 4096)
                self.assertEqual(cfg['convert_memory_gib'], 32)
                submit.assert_not_called()

    def test_false_does_not_submit(self):
        with patch.object(cli, 'load_config', return_value={}), \
             patch.object(cli, 'generate_split', return_value=Path('/test/jobs.json')) as generate, \
             patch.object(cli, 'submit_plan') as submit:
            cli.main(['split', '--submit', 'FALSE', '--queue', 'spec'])
            generate.assert_called_once()
            submit.assert_not_called()


if __name__ == '__main__':
    unittest.main()

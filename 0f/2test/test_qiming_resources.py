"""Admin's CPU/RAM allocation boundary and frozen-plan compatibility tests."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SHARED = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SHARED))
from qiming_resources import cpu_slots
from workflow_jobs import Job, load_plan, write_plan


class AdminResources(unittest.TestCase):
    def test_queue33_boundaries_use_larger_cpu_or_memory_share(self):
        for cores, memory, expected in [(1, 1, 10), (10, 48, 10), (11, 48, 20),
                (1, 48.001, 20), (20, 96, 20), (21, 96, 30), (1, 96.001, 30),
                (30, 144, 30), (31, 144, 40), (1, 144.001, 40), (24, 160, 40)]:
            with self.subTest(cores=cores, memory=memory):
                self.assertEqual(cpu_slots('33', cores, memory, 24), ('33', expected))

    def test_queue_selection_and_other_node_shapes(self):
        with patch.dict(os.environ, HPC_CPU_QUEUE='33', HPC_MEDIUM_QUEUE='38'):
            self.assertEqual(cpu_slots('auto', 24, 160, 24), ('33', 40))
            self.assertEqual(cpu_slots('auto', 32, 448, 120), ('38', 64))
        self.assertEqual(cpu_slots('38', 8, 64, 120), ('38', 16))
        self.assertEqual(cpu_slots('1t81c', 64, 896, 120), ('1t81c', 120))

    def test_admin_big_memory_exemptions_and_unknown_queues(self):
        for queue in ('v4-3t', 'v3-6t', 'other-site'):
            self.assertEqual(cpu_slots(queue, 24, 160, 24), (queue, 24))

    def test_invalid_or_oversized_requests_fail(self):
        for values in [(0, 32, 24), (1, -1, 24), (1, 32, 0), (1.5, 32, 24),
                       (True, 32, 24), (1, float('nan'), 24), (41, 32, 24), (1, 256, 24)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                cpu_slots('33', *values)

    def test_walltime_keeps_minutes_when_checking_queue_limit(self):
        args = [sys.executable, str(SHARED/'qiming_resources.py'), '--queue', '33',
                '--cores', '24', '--memory-gib', '160', '--format', 'shell', '--walltime']
        p = subprocess.run(args+['0:30'], text=True, capture_output=True)
        self.assertEqual((p.returncode, p.stdout.strip()), (0, '33 40'))
        p = subprocess.run(args+['720:01'], text=True, capture_output=True)
        self.assertNotEqual(p.returncode, 0)

    def test_hpc_rounds_allocation_without_changing_science_or_old_plans(self):
        with tempfile.TemporaryDirectory() as tmp:
            job = Job('cpu', ['tool', '--threads', '24'],
                      resources=dict(queue='33', cores=24, memory_gib=160), env={'THREADS': '24'})
            p = write_plan([job], tmp, mode='hpc')
            _, plan = load_plan(p)
            self.assertEqual(plan['jobs'][0]['resources']['cores'], 40)
            self.assertEqual(plan['jobs'][0]['argv'], job.argv)
            self.assertEqual(plan['jobs'][0]['env']['THREADS'], '24')
            # Legacy frozen plans retain their allocation when workers load them.
            plan['jobs'][0]['resources']['cores'] = 24
            p.write_text(json.dumps(plan))
            self.assertEqual(load_plan(p)[1]['jobs'][0]['resources']['cores'], 24)
            p = write_plan([job], tmp, mode='local')
            self.assertEqual(load_plan(p)[1]['jobs'][0]['resources']['cores'], 24)

    def test_gpu_allocations_are_not_rounded_by_cpu_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = write_plan([Job('gpu', ['true'], resources=dict(queue='4a100-40', cores=8,
                               memory_gib=96, gpus=1))], tmp)
            self.assertEqual(load_plan(p)[1]['jobs'][0]['resources']['cores'], 8)


if __name__ == '__main__':
    unittest.main()

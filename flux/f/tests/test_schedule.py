"""Warmup/cosine endpoints and absolute-step continuation semantics."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from train import scheduled_learning_rate


class ScheduleTests(unittest.TestCase):
    def test_warmup_cosine_endpoints_and_extension(self):
        cfg=dict(learning_rate=6e-4,warmup_steps=100,schedule_steps=2000,final_lr_ratio=.1)
        self.assertAlmostEqual(scheduled_learning_rate(cfg,1),6e-6)
        self.assertAlmostEqual(scheduled_learning_rate(cfg,100),6e-4)
        self.assertAlmostEqual(scheduled_learning_rate(cfg,1050),3.3e-4)
        self.assertAlmostEqual(scheduled_learning_rate(cfg,2000),6e-5)
        self.assertAlmostEqual(scheduled_learning_rate(cfg,2200),6e-5)
        cfg['max_steps']=4000
        self.assertAlmostEqual(scheduled_learning_rate(cfg,1050),3.3e-4)

    def test_constant_and_invalid_schedule(self):
        self.assertEqual(scheduled_learning_rate({'learning_rate':.001},1000),.001)
        with self.assertRaises(ValueError):
            scheduled_learning_rate(dict(learning_rate=.001,warmup_steps=100,schedule_steps=50),1)


if __name__=='__main__': unittest.main()

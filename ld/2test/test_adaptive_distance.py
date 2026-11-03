"""Boundary evidence, not block counts, controls adaptive distance selection."""
import os
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'f'))
import ld_pipeline as ld

class AdaptiveDistanceTests(unittest.TestCase):
    def pipeline(self, directory, signatures, span=160000, auto=True):
        p=object.__new__(ld.Pipeline)
        p.auto_extend_d=auto;p.confirmations=2;p.grid=[10000,20000,40000]
        p.work=Path(directory);p.reference_tolerance=.1
        p.pgen_cache={(pop,1):(None,{'qc':{'first_position':1,'last_position':span*1000+1}}) for pop in ld.POPS}
        p.reuse_selection=lambda c,d:None
        p.reference_preflight=lambda c,d:None
        p.calls=[]
        p.partitions=lambda pop,c,ds:p.calls.append((pop,list(ds)))
        p.summarize=lambda c,ds:{d:{'signature':signatures.get(d,'new'), 'full_span_all_populations':d>=span,
                                  'counts':{'AFR':10},'max_reference_relative_deviation':None} for d in ds}
        p.check_reference_counts=lambda *a,**kw:None
        p.audited_selection_record=lambda c,d,ds,ss:dict(D_kb=d,tested_grid_kb=ds)
        return p

    def test_unstable_40mb_extends_to_full_span(self):
        with TemporaryDirectory() as tmp:
            p=self.pipeline(tmp,{10000:'a',20000:'b',40000:'c',80000:'d',160000:'e'})
            chosen,record=p.evaluate_chromosome(1,p.grid)
            self.assertEqual(chosen,160000)
            self.assertEqual(record['tested_grid_kb'],[10000,20000,40000,80000,160000])
            self.assertEqual(p.grid,[10000,20000,40000])

    def test_requires_two_larger_equal_boundaries_not_equal_counts(self):
        with TemporaryDirectory() as tmp:
            p=self.pipeline(tmp,{10000:'a',20000:'b',40000:'b',80000:'b'},span=250000)
            chosen,record=p.evaluate_chromosome(1,p.grid)
            self.assertEqual(chosen,20000)
            self.assertEqual(record['tested_grid_kb'][-1],80000)

    def test_explicit_fixed_cap_still_stops_without_publication(self):
        with TemporaryDirectory() as tmp:
            p=self.pipeline(tmp,{10000:'a',20000:'b',40000:'c'},auto=False)
            with self.assertRaisesRegex(ld.Stop,'D stability is not yet confirmed'):
                p.evaluate_chromosome(1,p.grid)
            self.assertFalse((p.work/'selected/chr1.json').exists())

    def test_full_span_is_ceil_of_observed_not_hardcoded_build_size(self):
        with TemporaryDirectory() as tmp:
            p=self.pipeline(tmp,{10000:'a',20000:'b',40000:'c'},span=48001)
            chosen,record=p.evaluate_chromosome(1,p.grid)
            self.assertEqual(chosen,48001)

    def test_already_valid_chromosome_does_not_extend(self):
        with TemporaryDirectory() as tmp:
            p=self.pipeline(tmp,{})
            p.reuse_selection=lambda c,d:(10000,{'tested_grid_kb':d})
            self.assertEqual(p.evaluate_chromosome(1,p.grid)[0],10000)
            self.assertEqual(p.calls,[])

if __name__=='__main__': unittest.main()

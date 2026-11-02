"""Real chromosome audits, merged index, and downstream exact-PVAR contract."""
import argparse
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'f'))
import ukb_mfi
import mfi_merge
from analysis_integrity import validate_manifest
from panel_targets import quality_contract,load_quality


def chromosome(root,c,min_info=.8):
    mfi=root/f'chr{c}.mfi';pvar=root/f'chr{c}.pvar';out=root/f'out{c}'
    mfi.write_text(f'{c}:100_A_G rs{c} 100 A G 0.1 G 0.95\n')
    pvar.write_text(f'##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\n{c}\t100\trs{c}\tA\tG\n')
    ukb_mfi.main(['--mfi',str(mfi),'--pvar',str(pvar),'--chromosomes',str(c),'--genome-build','GRCh37','--min-info',str(min_info),'--tabix','--out',str(out)])
    return out,pvar


class MergedMfi(unittest.TestCase):
    def test_real_two_chromosomes_tabix_and_downstream_binding(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,pa=chromosome(root,1);b,pb=chromosome(root,2)
            m=mfi_merge.merge([str(b),str(a)],root/'merged')
            self.assertEqual(m['counts']['quality_pass'],2)
            validate_manifest(root/'merged',require_evidence=True)
            with pysam.TabixFile(str(root/'merged/target_info.tsv.gz')) as tabix:
                self.assertEqual(list(tabix.fetch('2',99,100)),['2\t100\tA\tG\t0.95'])
            args=argparse.Namespace(quality_info=str(root/'merged/target_info.tsv.gz'),panel_build='GRCh37',quality_build=None,evidence_kind='imputed_hardcall')
            proof=quality_contract(args,[dict(pvar=str(pa)),dict(pvar=str(pb))])
            self.assertEqual(proof['status'],'verified_mfi_output_and_original_PVAR')
            self.assertEqual(len(load_quality(args.quality_info,[dict(chrom='1',pos=100),dict(chrom='2',pos=100)])),2)
    def test_modified_original_pvar_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,pa=chromosome(root,1);pa.write_text(pa.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'Original source differs'):
                mfi_merge.merge([str(a)],root/'merged')
    def test_duplicate_chromosomes_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,_=chromosome(root,1)
            with self.assertRaisesRegex(ValueError,'Overlapping chromosome'):
                mfi_merge.merge([str(a),str(a)],root/'merged')
    def test_inconsistent_thresholds_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a,_=chromosome(root,1,.8);b,_=chromosome(root,2,.9)
            with self.assertRaisesRegex(ValueError,'thresholds/layout differ'):
                mfi_merge.merge([str(a),str(b)],root/'merged')


if __name__=='__main__':unittest.main()

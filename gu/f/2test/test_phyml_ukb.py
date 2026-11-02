#!/usr/bin/env python3
"""PhyML/UKB regression tests; optional host CUDA-pruning parity uses PHYML_UPSTREAM."""
import argparse
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "f"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


phyml = load("phyml_audit_test_subject", ROOT / "f/phyml.py")
ukb = load("ukb_audit_test_subject", ROOT / "f/ukb_prepare.py")


class ScientificInputTests(unittest.TestCase):
    def test_ukb_runtime_roots_and_default_references(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, "GU_REF_ROOT is unset"):
                ukb.runtime_directory("GU_REF_ROOT")
            os.environ["GU_REF_ROOT"] = "/configured/references"
            self.assertEqual(ukb.runtime_directory("GU_REF_ROOT"), Path("/configured/references"))
        self.assertEqual(ukb.DEFAULT_ARCHAIC_REFERENCES, ("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija"))

    def test_archaic_duplicate_cannot_erase_heterozygote(self):
        rows = [["1", "100", ".", "A", "G", "0/1"], ["1", "100", ".", "A", "T", "0/0"]]
        audit = {}
        with mock.patch.object(phyml, "query_rows", return_value=(["Altai"], iter(rows))):
            calls = phyml.archaic_calls(Path("unused.vcf"), "1", {"start": 0, "end": 1000}, [{"pos":100,"ref":"A","alt":"G"}], True, audit)
        self.assertEqual(calls[100], "N")
        self.assertEqual(audit[100], "duplicate_archaic_position")

    def test_archaic_heterozygote_with_third_allele_is_unknown(self):
        rows = [["1", "100", ".", "A", "T", "0/1"]]
        with mock.patch.object(phyml, "query_rows", return_value=(["Altai"], iter(rows))):
            calls = phyml.archaic_calls(Path("unused.vcf"), "1", {"start": 0, "end": 1000}, [{"pos":100,"ref":"A","alt":"G"}], True)
        self.assertEqual(calls[100], "N")

    def test_archaic_lookup_does_not_confuse_chr1_and_chr11(self):
        phyml.archaic_vcf.cache_clear()
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "Altai").mkdir()
            (root / "Altai/Altai.chr11.vcf.gz").write_text("unrelated")
            with self.assertRaises(SystemExit):
                phyml.archaic_vcf(root, "Altai", "1")

    def test_ambiguous_or_downloading_archaic_reference_rejected(self):
        phyml.archaic_vcf.cache_clear()
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            directory = root / "Altai"
            directory.mkdir()
            (directory / "old.chr1.vcf.gz").write_text("one")
            (directory / "new.chr1.vcf.gz").write_text("two")
            with self.assertRaises(SystemExit):
                phyml.archaic_vcf(root, "Altai", "1")
            (directory / "old.chr1.vcf.gz.aria2").touch()
            self.assertEqual(phyml.archaic_vcf(root,"Altai","1").name,"new.chr1.vcf.gz")

    def test_effect_allele_stays_correct_when_ref_alt_swapped(self):
        row = dict(lead_pos=100,ref="A",alt="G",index_snp="1:100:A:G",effect_allele="G",beta_j="0.2")
        lead = phyml.exact_lead(row, [["1","100",".","G","A","0|1","1|1"]], 2)
        self.assertEqual((lead["risk_allele"],lead["copies"],lead["allele_order"]),("G",[1,0,0,0],"swapped"))

    def test_haploid_x_uses_one_copy(self):
        row = dict(lead_pos=3000000,ref="A",alt="G",index_snp="X:3000000:A:G",effect_allele="G",beta_j="0.2")
        lead = phyml.exact_lead(row, [["X","3000000",".","A","G","1","0"]], 2, True)
        self.assertEqual(lead["copies"],[1,-1,0,-1])
        with self.assertRaises(SystemExit):
            phyml.called_haploid_base("A","G","1|1")

    def test_phased_ld_uses_requested_copies(self):
        risk = [0,0,1,1,0,1,0,1]
        site = dict(ref="A",alt="G",haps="AAGGGGAA")
        r2, n, allele = phyml.phased_ld(risk,site,[0,1,2,3])
        self.assertEqual((r2,n,allele),(1.0,4,"G"))
        self.assertLess(phyml.phased_ld(risk,site,list(range(8)))[0],1)

    def test_public_resume_options_and_explicit_boolean(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "cmd.sh"
            path.write_text("exec ./gu.sh phyml --dir-target /input/chr --keep-males FALSE --chr 1\n")
            _, opts = phyml.command(path)
            self.assertEqual(opts["--dir-target"],"/input/chr")
            self.assertEqual(opts["--keep-males"],"FALSE")

    def test_duplicate_sex_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "samples.tsv"
            path.write_text("sample\tsex\nperson\t1\nperson\t2\n")
            with self.assertRaises(SystemExit):
                phyml.read_sexes(path)

    def test_pvar_ordinal_selection_and_ambiguous_coordinate_exclusion(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            pvar, out, indexed = folder / "chr1.pvar", folder / "selected", folder / "indexed.pvar"
            pvar.write_text("#CHROM\tPOS\tID\tREF\tALT\n1\t10\trs1\tA\tG\n1\t20\trs1\tA\tT\n1\t20\trs2\tA\tC\n1\t30\trs1\tG\tT\n")
            args = argparse.Namespace(ukb_info_file=None,ukb_qc_pass_variants=None,regions=None,existing_pgen=True,grch=37)
            qc = ukb.select_variants(pvar,"1",args,out,indexed)
            self.assertEqual(out.read_text().splitlines(),["gu_row_1","gu_row_4"])
            self.assertEqual(qc["ambiguous_positions_excluded"],1)
            self.assertEqual(len(indexed.read_text().splitlines()),5)

    def test_large_stderr_does_not_deadlock_site_stream(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            tool = folder / "bcftools"
            tool.write_text("#!/usr/bin/env python3\nimport sys\nsys.stderr.write('warning\\n'*200000)\nsys.stdout.write('1\\t100\\tA\\tG\\n')\n")
            tool.chmod(0o755)
            env = dict(os.environ, PATH=str(folder)+os.pathsep+os.environ["PATH"], PYTHONPATH=str(ROOT / "f"))
            script = "from pathlib import Path; import sys, ukb_prepare; q=ukb_prepare._ukb_scoring_sites(Path('fake.vcf'),'1',37,Path(sys.argv[1])); assert q['modern_sites']==1"
            subprocess.run([sys.executable,"-c",script,str(folder / "sites.tsv")], check=True, timeout=10, capture_output=True, text=True, env=env)


class PruningParityTests(unittest.TestCase):
    def test_cuda_pruning_source_matches_upstream_on_host(self):
        upstream = os.environ.get("PHYML_UPSTREAM")
        if not upstream or not shutil.which("g++"):
            self.skipTest("set PHYML_UPSTREAM and install g++ for the host-only pruning parity test")
        source = (Path(upstream) / "src/lk.c").read_text()
        def function(text, name):
            start = text.index("void " + name + "(")
            opening = text.index("{",start)
            level, end = 1, opening+1
            while level:
                level += (text[end]=="{") - (text[end]=="}")
                end += 1
            return text[start:end]
        cpu = "\n".join(function(source,name) for name in ["Partial_Lk_Exex","Partial_Lk_Exin","Partial_Lk_Inin","Core_Default_Update_Partial_Lk"])
        cuda = (ROOT / "f/phyml_gpu/partials.cu").read_text()
        gpu = cuda[cuda.index("struct Inputs {"):cuda.index('extern "C" void gu_cuda_partial')].replace("__global__","")
        prefix = r'''
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using phydbl=double;
#define YES 1
#define NO 0
#define SMALL DBL_MIN
#define BIG DBL_MAX
#define INV_TWO_TO_THE_LARGE 0x1p-256
#define TWO_TO_THE_LARGE 0x1p256
#define LARGE 256
struct seq { unsigned char *is_ambigu; int *d_state; };
struct t_node { int tax; seq *c_seq; };
struct Dim { int x; } blockIdx, blockDim, threadIdx;
'''
        main = r'''
int main() {
  std::mt19937_64 rng(734578); const int n=17;
  for(int trial=0;trial<128;++trial) for(int k:{1,4}) for(int t1:{0,1}) for(int t2:{0,1}) for(int scaling:{0,1}) {
    std::vector<double> a(n*(t1?1:k)*4),b(n*(t2?1:k)*4),p(k*16),q(k*16),w(n,1),out(n*k*4),expect(n*k*4);
    std::vector<int> s1(n),s2(n),st1(n,-1),st2(n,-1),scale(n),expected_scale(n);
    std::vector<unsigned char> ambiguous1(n,1),ambiguous2(n,1);
    for(auto* v:{&p,&q}) for(int row=0;row<k*4;++row) {
      double total=0;for(int j=0;j<4;++j) total+=((*v)[row*4+j]=1+(rng()%100));
      for(int j=0;j<4;++j) (*v)[row*4+j]/=total;
    }
    for(auto* v:{&a,&b}) for(auto& x:*v) x=(rng()%1000)/1000.0*(trial%3==0?1e-100:1);
    for(int i=0;i<n;++i) {
      s1[i]=t1?0:(rng()%3)*256; s2[i]=t2?0:(rng()%3)*256;
      if(t1) { if(i%3) { st1[i]=rng()%4; ambiguous1[i]=0; } for(int j=0;j<4;++j) a[i*4+j]=st1[i]<0?1:(j==st1[i]); }
      if(t2) { if(i%3) { st2[i]=rng()%4; ambiguous2[i]=0; } for(int j=0;j<4;++j) b[i*4+j]=st2[i]<0?1:(j==st2[i]); }
      if(i%5==0) w[i]=0;
    }
    seq c1{ambiguous1.data(),st1.data()},c2{ambiguous2.data(),st2.data()};t_node n1{t1,&c1},n2{t2,&c2};
    Core_Default_Update_Partial_Lk(&n1,&n2,expect.data(),a.data(),b.data(),p.data(),q.data(),expected_scale.data(),s1.data(),s2.data(),4,k,n,scaling,w.data());
    Inputs input{a.data(),b.data(),p.data(),q.data(),w.data(),s1.data(),s2.data(),st1.data(),st2.data()};
    blockDim.x=128;blockIdx.x=0;
    for(int i=0;i<n;++i) { threadIdx.x=i; prune(n,k,t1,t2,scaling,input,out.data(),scale.data()); }
    if(out!=expect || scale!=expected_scale) {std::cerr<<"pruning mismatch trial="<<trial<<"\n";return 1;}
  }
  std::cout<<"2048 pruning configurations matched the pinned upstream CPU function\n";
}
'''
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            path,binary = folder / "prune.cc",folder / "prune"
            path.write_text(prefix+cpu+gpu+main)
            subprocess.run(["g++","-O2","-std=c++17","-ffp-contract=off",str(path),"-o",str(binary)],check=True,capture_output=True,text=True)
            result = subprocess.run([str(binary)],check=True,capture_output=True,text=True)
            self.assertIn("2048",result.stdout)


if __name__ == "__main__":
    unittest.main()

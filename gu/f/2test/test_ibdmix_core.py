#!/usr/bin/env python3
"""Regression checks for IBDmix integration; no reference datasets required.

Run: python3 -m unittest discover -s tests -p 'test_ibdmix_core.py' -v
Optional upstream differential test: set IBDMIX_UPSTREAM to the official source
checkout. That test executes the actual CUDA Stack source on the CPU and checks
it against upstream C++; it does not constitute CUDA-device validation.
"""
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "f"))
spec = importlib.util.spec_from_file_location("ibdmix_test_subject", ROOT / "f/ibdmix.py")
subject = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subject)


class IntegrationTests(unittest.TestCase):

    def test_axt_explicit_root_accepts_nested_and_flat_layout(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            name = 'chr1.hg19.panTro2.synNet.axt.gz'
            nested = root / 'vsPanTro2'
            nested.mkdir()
            (nested / name).write_bytes(b'fixture')
            with mock.patch.dict(os.environ, {'IBDMIX_AXT_DIR': str(root)}):
                self.assertEqual(subject.axt_reference('1', 'panTro2')[1], nested / name)
                self.assertEqual(subject.axt_reference('1', 'panTro2', root)[1], nested / name)
            (nested / name).unlink()
            (root / name).write_bytes(b'fixture')
            self.assertEqual(subject.axt_reference('1', 'panTro2', root)[1], root / name)

    @staticmethod
    def daf_fixture(folder, references, mask_reference):
        folder = Path(folder)
        run = folder / "run"
        (run / "final").mkdir(parents=True)
        (run / "raw/C1").mkdir(parents=True)
        (run / "samples/C1").mkdir(parents=True)
        meta = dict(genome_build="b37", coordinate_system="0-based-half-open", refs=" ".join(references))
        (run / "run.meta.tsv").write_text("".join(f"{key}\t{value}\n" for key, value in meta.items()))
        populations = ["ESN", "GWD", "LWK", "MSL", "YRI", "EUR"]
        (run / "samples/C1/populations.tsv").write_text("population\tsuper_population\tn\n" + "".join(f"{pop}\t{pop}\t1\n" for pop in populations))
        for pop in populations:
            (run / f"samples/C1/{pop}.txt").write_text(pop + "_person\n")
            for ref in references:
                with gzip.open(run / f"raw/C1/{ref}.{pop}.raw.txt.gz", "wt") as handle:
                    handle.write("ID\tchrom\tstart\tend\tslod\n")
        source = run / "final/calls.tsv.gz"
        with gzip.open(source, "wt") as handle:
            handle.write("sample_id\tchrom\tstart\tend\tlength\tanc\tscore_scope\n")
            for ref in references:
                if ref in subject.DAF_REFERENCES:
                    handle.write(f"EUR_person\t1\t0\t200000\t200000\t{ref}\tnative\n")
        mask = folder / "author-daf.bed"
        mask.write_text("chr1\t75000\t125000\n")
        sidecar = dict(genome_build="GRCh37", reference=mask_reference, percentile=99.9,
                       coordinate_system="0-based-half-open", source="synthetic regression fixture; not an author mask")
        Path(str(mask) + ".json").write_text(json.dumps(sidecar))
        return source, meta, mask, sidecar

    def test_daf_2013_pair_applies_only_explicit_release(self):
        with tempfile.TemporaryDirectory() as folder:
            source, meta, mask, _ = self.daf_fixture(folder, ["Altai.2013", "Denisova.2013"], "Altai.2013")
            result, audit = subject.prepare(source, Path(folder) / "filtered", meta, mask)
            data = list(subject.rows(result))
            self.assertEqual([(int(row["start"]), int(row["end"])) for row in data], [(0, 75000), (125000, 200000)])
            self.assertEqual(audit["daf"]["status"], "applied")
            self.assertEqual(audit["daf"]["applied_references"], ["Altai.2013"])
            self.assertEqual(audit["daf"]["sha256"], subject.checksum(mask))
            self.assertTrue((result.parent / "Altai.2013_DAF.bed").exists())
            self.assertFalse((result.parent / "Altai_DAF.bed").exists())
            self.assertFalse((result.parent / "Altai_overlap.bed").exists())
            self.assertEqual(len(list(subject.rows(source))), 1)  # Native data stays intact.
            self.assertEqual(subject.prepare(source, Path(folder) / "filtered", meta, mask)[0], result)

    def test_five_reference_filter_creates_no_unselected_reference_artifacts(self):
        with tempfile.TemporaryDirectory() as folder:
            source, meta, _, _ = self.daf_fixture(folder, subject.DEFAULT_REFERENCES, "Altai")
            result, audit = subject.prepare(source, Path(folder) / "filtered", meta)
            self.assertEqual(audit["reference_scope"], sorted(subject.DEFAULT_REFERENCES))
            self.assertFalse(any(".2013" in path.name for path in result.parent.iterdir()))
            self.assertEqual(audit["denisovan_overlap_references"], {"Denisova": "Altai", "Denisova25": "Altai"})

    def test_daf_masks_never_cross_reference_releases(self):
        for mask_ref, run_ref in [("Altai", "Altai.2013"), ("Altai.2013", "Altai")]:
            with self.subTest(mask_ref=mask_ref), tempfile.TemporaryDirectory() as folder:
                source, meta, mask, _ = self.daf_fixture(folder, [run_ref], mask_ref)
                result, audit = subject.prepare(source, Path(folder) / "filtered", meta, mask)
                self.assertEqual(list(subject.rows(result)), list(subject.rows(source)))
                self.assertEqual(audit["daf"]["status"], "not_applied_reference_mismatch")
                self.assertEqual(audit["daf"]["applied_references"], [])
                self.assertEqual(audit["daf"]["unfiltered_altai_references"], [run_ref])

    def test_daf_mixed_release_audit_is_partial(self):
        for mask_ref in ["Altai", "Altai.2013"]:
            with self.subTest(mask_ref=mask_ref), tempfile.TemporaryDirectory() as folder:
                source, meta, mask, _ = self.daf_fixture(folder, ["Altai", "Altai.2013"], mask_ref)
                result, audit = subject.prepare(source, Path(folder) / "filtered", meta, mask)
                data = list(subject.rows(result))
                self.assertEqual([int(row["length"]) for row in data if row["anc"] == mask_ref], [75000, 75000])
                self.assertEqual([int(row["length"]) for row in data if row["anc"] != mask_ref], [200000])
                self.assertEqual(audit["daf"]["status"], "partially_applied_reference_mismatch")
                self.assertEqual(audit["daf"]["applied_references"], [mask_ref])

    def test_daf_metadata_requires_exact_reference_and_valid_coordinates(self):
        with tempfile.TemporaryDirectory() as folder:
            _, _, mask, sidecar = self.daf_fixture(folder, ["Altai.2013"], "Altai.2013")
            invalid = [{"reference": "altai.2013"}, {"reference": "Altai 2013"}, {"reference": ["Altai", "Altai.2013"]},
                       {"genome_build": "GRCh38"}, {"percentile": 99}, {"source": " "}, {"coordinate_system": "1-based-inclusive"}]
            for override in invalid:
                with self.subTest(override=override):
                    Path(str(mask) + ".json").write_text(json.dumps(dict(sidecar, **override)))
                    with self.assertRaises(ValueError):
                        subject.load_daf(mask)

    def test_missing_daf_does_not_claim_application(self):
        with tempfile.TemporaryDirectory() as folder:
            source, meta, _, _ = self.daf_fixture(folder, ["Altai.2013", "Denisova.2013"], "Altai.2013")
            result, audit = subject.prepare(source, Path(folder) / "filtered", meta)
            self.assertEqual(list(subject.rows(result)), list(subject.rows(source)))
            self.assertEqual(audit["daf"]["status"], "not_applied_missing_author_mask")
            self.assertEqual(audit["daf"]["unfiltered_altai_references"], ["Altai.2013"])

    def test_daf_filter_requires_explicit_reference(self):
        with self.assertRaisesRegex(ValueError, "explicit exact reference"):
            list(subject.filter_rows("unused", {}, {"1": [(0, 100000)]}, "1", defaultdict(Counter)))

    def test_public_gpu_and_directory_options_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run.sh"
            path.write_text("exec ./gu.sh ibdmix --GPU TRUE --dir-target '/target with space/chr' --archaic-path /refs --keep-males --chr X\n")
            env, opts = subject.resume_command(path)
            self.assertEqual(env["IBDMIX_BACKEND"], "cuda")
            self.assertEqual(opts["--dir-target"], "/target with space/chr")
            self.assertEqual(opts["--archaic-path"], "/refs")
            self.assertEqual(opts["--keep-males"], "TRUE")

    def test_explicit_false_males_and_cpu_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run.sh"
            path.write_text("exec ./gu.sh ibdmix --GPU FALSE --keep-males FALSE --chr 1\n")
            env, opts = subject.resume_command(path)
            self.assertEqual(env["IBDMIX_BACKEND"], "cpu")
            self.assertEqual(opts["--keep-males"], "FALSE")

    def test_unknown_module_and_options_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run.sh"
            for command in ["exec ./gu.sh unknown --GPU FALSE", "exec ./gu.sh ibdmix --unknown-option value",
                            "exec ./gu.sh ibdmix --GPU TRUE --GPU FALSE"]:
                with self.subTest(command=command):
                    path.write_text(command + "\n")
                    with self.assertRaises(ValueError):
                        subject.resume_command(path)

    def test_optional_references_require_explicit_selection(self):
        result = subprocess.run([sys.executable, str(ROOT / "f/ibdmix.py"), "normalize-refs", "altai.2013,denisova.2013"], check=True, capture_output=True, text=True)
        self.assertEqual(result.stdout.strip(), "Altai.2013 Denisova.2013")
        self.assertEqual(subject.DEFAULT_REFERENCES, ("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija"))
        self.assertNotIn("Altai.2013", subject.DEFAULT_REFERENCES)
        self.assertNotIn("Denisova.2013", subject.DEFAULT_REFERENCES)

    def test_mask_reuse_requires_current_implementation_hash(self):
        self.assertTrue(subject.mask_code_matches(subject.sha256_file(subject.__file__)))
        self.assertFalse(subject.mask_code_matches("0" * 64))

    def test_axt_directory_uses_configured_roots(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, "GU_DATA_ROOT is unset"):
                subject.axt_directory()
            os.environ["GU_DATA_ROOT"] = "/configured/data"
            self.assertEqual(subject.axt_directory(), Path("/configured/data/annot/axt/37"))
            os.environ["GU_ANNOT_ROOT"] = "/configured/annotations"
            self.assertEqual(subject.axt_directory(), Path("/configured/annotations/axt/37"))
            os.environ["IBDMIX_AXT_DIR"] = "/configured/specific-axt"
            self.assertEqual(subject.axt_directory(), Path("/configured/specific-axt"))

    def test_provenance_requires_same_reference_selection_and_tabs(self):
        first = "pipeline_version\t2026.1\nrefs\tAltai\nexcluded_mask_sha256\tC1\tAltai\t" + "1" * 64 + "\n"
        second = first.replace("Altai", "Altai.2013")
        self.assertEqual(subject.differences(first, first), "")
        self.assertTrue(subject.differences(first, second))
        with self.assertRaises(ValueError):
            subject.contract(first.replace("\t", "\\t"))

    def test_driver_writes_current_tab_delimited_contract(self):
        driver = (ROOT / "f/ibdmix.sh").read_text()
        function = "run_record() {" + driver.split("run_record() {", 1)[1].split("\ncheck_run_provenance()", 1)[0]
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            (folder / "mask/C1").mkdir(parents=True)
            for name in ["samples.txt", "prepare", "caller", "fasta", "modern.vcf.gz", "modern.vcf.gz.tbi"]:
                (folder / name).write_text("synthetic metadata fixture\n")
            for ref in subject.DEFAULT_REFERENCES:
                (folder / ref).mkdir()
                (folder / ref / "chr1.vcf.gz").write_text("synthetic metadata fixture\n")
                (folder / "mask/C1" / f"{ref}.exclude.bed").write_text("1\t0\t1\n")
            setup = r'''
set -euo pipefail
F=$IBDMIX_TEST_F
work=$IBDMIX_TEST_WORK
pipeline_version=2026-test
profile=multi_reference
background_filter=1
mask_root=$work/mask
custom_masks=
export_denisovan=1
genome_build=b37
IBDMIX_LOCUS_FLANK_BP=0
refs_arg='Altai Chagyr Denisova Denisova25 Vindija'
refs=(Altai Chagyr Denisova Denisova25 Vindija)
analysis_units=(C1)
declare -A unit_chr=([C1]=1) unit_sex=([C1]=all)
lod_cut=4
len_cut=50000
emit_lod_cut=4
minor_allele_count=1
archaic_error=0.01
modern_error_max=0.002
modern_error_proportion=2
sample_file=$work/samples.txt
loci_file=
prepare_bin=$work/prepare
ibdmix_bin=$work/caller
helper=$F/ibdmix.py
ibdmix_mask_root=$work/mask
ibdmix_tmp_root=$work/tmp
ref_fasta=$work/fasta
modern_vcf() { printf '%s/modern.vcf.gz\n' "$work"; }
archaic_vcf() { printf '%s/%s/chr1.vcf.gz\n' "$work" "$1"; }
'''
            script = folder / "contract.sh"
            script.write_text(setup + function + "\nrun_record\n")
            env = dict(os.environ, IBDMIX_TEST_F=str(ROOT / "f"), IBDMIX_TEST_WORK=str(folder), IBDMIX_BACKEND="cpu", GU_BATCH_ID="001", GU_POSTFIX="002")
            result = subprocess.run(["bash", str(script)], env=env, check=True, text=True, capture_output=True)
            self.assertNotIn("\\t", result.stdout)
            self.assertIn("refs\tAltai Chagyr Denisova Denisova25 Vindija", result.stdout)
            self.assertIn("batch_id\t001\n", result.stdout)
            self.assertTrue(subject.contract(result.stdout))
            self.assertEqual(len([row for row in subject.resume_records(result.stdout) if row[0] == "excluded_mask_sha256"]), 5)

    def test_identical_prepare_binary_rebuild_keeps_contract(self):
        before = "\n".join([
            "pipeline_version\tshared-hts-v1", "refs\tAltai",
            "excluded_mask_sha256\tC1\tAltai\t" + "1" * 64,
            "software\t/tools/prepare:123:1", "software\t/tools/ibdmix:456:1",
            "genotype_binary_sha256\t" + "2" * 64,
        ]) + "\n"
        self.assertEqual(subject.contract(before), subject.contract(before.replace("prepare:123:1", "prepare:123:2")))
        self.assertNotEqual(subject.contract(before), subject.contract(before.replace("2" * 64, "3" * 64)))
        self.assertNotEqual(subject.contract(before), subject.contract(before.replace("ibdmix:456:1", "ibdmix:456:2")))

    def test_singleton_is_discarded_without_losing_valid_segment(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            (folder / "samples.txt").write_text("person\n")
            (folder / "populations.tsv").write_text(f"population\tsuper_population\tn\tsample_file\nEUR\tEUR\t1\t{folder / 'samples.txt'}\n")
            with gzip.open(folder / "Altai.EUR.raw.txt.gz", "wt") as handle:
                handle.write("ID\tchrom\tstart\tend\tslod\n")
                handle.write("person\t1\t100\t50100\t6\n")
                handle.write("person\t1\t60000\t60000\t5\n")
            output = folder / "segments.gz"
            subject.finalize(folder, folder / "populations.tsv", ["Altai"], output, "1", "b37", "chr1", background=False)
            with gzip.open(output, "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual((rows[0]["start"], rows[0]["end"], rows[0]["length"]), ("99", "50099", "50000"))

    def test_x_nonpar_components_are_separate(self):
        for build in ("37", "38"):
            with self.subTest(build=build), tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                source, output = folder / "gt.gz", folder / "part"
                spans = subject.x_nonpar_intervals(build)
                with gzip.open(source, "wt") as handle:
                    handle.write("chrom\tpos\tref\talt\tarchaic\tperson\n")
                    for lo, hi in spans:
                        handle.write(f"23\t{lo+1}\tA\tC\t1\t2\n")
                        handle.write(f"23\t{hi}\tA\tC\t1\t0\n")
                subject.split_x_genotypes(source, output, build)
                for i, (lo, hi) in enumerate(spans, 1):
                    with gzip.open(f"{output}.{i}.gz", "rt") as handle:
                        records = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual([int(r["pos"]) for r in records], [lo+1, hi])

    def test_x_par_positions_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "gt.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("chrom\tpos\tref\talt\tarchaic\tperson\n23\t60001\tA\tC\t1\t2\n")
            with self.assertRaises(ValueError):
                subject.split_x_genotypes(source, folder / "part", "37")


class UpstreamDifferentialTests(unittest.TestCase):
    def test_native_shared_prepare_matches_upstream_merger(self):
        upstream = os.environ.get("IBDMIX_UPSTREAM")
        prepare = os.environ.get("IBDMIX_PREPARE_BIN")
        if not upstream or not prepare or not shutil.which("g++"):
            self.skipTest("set IBDMIX_UPSTREAM and IBDMIX_PREPARE_BIN for native merger parity")
        upstream = Path(upstream)
        original = (upstream / "src/generate_gt.cc").read_text()
        # Use the official merger unchanged, replacing only its CLI11 parser.
        source = ('#include <fstream>\n#include <iostream>\n#include "IBDmix/vcf_file.h"\n'
                  'int main(int argc,char** argv) { std::string archaic_file=argv[1], modern_file=argv[2], outfile=argv[3];\n'
                  + original[original.index("  std::ofstream of;"):])
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            cpp, binary = folder / "generate.cc", folder / "generate"
            cpp.write_text(source)
            subprocess.run(["g++", "-std=c++17", "-O2", "-I", str(upstream / "include"), str(cpp), str(upstream / "src/vcf_file.cc"), "-o", str(binary)], check=True, capture_output=True, text=True)
            def vcf(path, names, records):
                path.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                                + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(names) + '\n'
                                + ''.join(f'1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t' + '\t'.join(gt) + '\n' for pos, ref, alt, gt in records))
            modern = folder / "modern.vcf"
            modern_rows = [(90,"A","C"), (100,"A","C"), (110,"A","CT"), (120,"C","G"), (130,"T","G"), (150,"A","G"), (160,"A","C"), (170,"T","."), (180,"A","C"), (190,"G","A")]
            vcf(modern, ["p0","p1","p2","p3"], [(p,r,a,["0/0","0/1","1/1","./."]) for p,r,a in modern_rows])
            archaic_rows = [(95,"A","C","0/1"), (100,"A",".","0/0"), (110,"A","C","0/1"), (115,"A","C","1/1"), (120,"C","T","1/1"), (125,"A",".","0/0"), (130,"A","G","0/1"), (150,"A","G","0/1"), (160,"A","C","./."), (170,"T","C","0/1"), (180,"A",".","0/0"), (185,"A","C","0/1"), (190,"G","A","0/1"), (200,"A","C","0/1")]
            args = [prepare, "--modern", str(modern), "--out", str(folder / "output"), "--unit", "C1", "--chrom", "1", "--threads", "2"]
            for name in ("Altai", "Denisova"):
                path = folder / f"{name}.vcf"
                vcf(path, [name], [(p,r,a,[g]) for p,r,a,g in archaic_rows])
                args.extend(["--archaic", name, str(path), ""])
                subprocess.run([str(binary), str(path), str(modern), str(folder / f"{name}.expected")], check=True, capture_output=True, text=True)
            subprocess.run(args, check=True, capture_output=True, text=True)
            for name in ("Altai", "Denisova"):
                with gzip.open(folder / "output" / f"{name}.C1.gt.txt.gz", "rt") as handle:
                    self.assertEqual(handle.read(), (folder / f"{name}.expected").read_text())

    def test_native_haploid_and_partial_missing_encoding(self):
        prepare = os.environ.get("IBDMIX_PREPARE_BIN")
        if not prepare:
            self.skipTest("set IBDMIX_PREPARE_BIN for native haploid/missing encoding")
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            meta = '##fileformat=VCFv4.2\n##contig=<ID=X,length=155270560>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            modern, archaic = folder / "modern.vcf", folder / "archaic.vcf"
            modern.write_text(meta + columns + 'p0\tp1\tp2\tp3\tp4\tp5\n'
                             + 'X\t3000000\t.\tA\tC\t.\tPASS\t.\tGT\t0\t1\t.\t0/.\t./1\t1|1\n')
            archaic.write_text(meta + columns + 'Altai\nX\t3000000\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n')
            bed = folder / "nonpar.bed"
            bed.write_text(''.join(f'23\t{lo}\t{hi}\n' for lo,hi in subject.x_nonpar_intervals("37")))
            subprocess.run([prepare,"--modern",str(modern),"--out",str(folder),"--unit","X_MALE","--chrom","X","--nonpar-bed",str(bed),"--build","37","--threads","1","--archaic","Altai",str(archaic),""], check=True, capture_output=True, text=True)
            with gzip.open(folder / "Altai.X_MALE.gt.txt.gz", "rt") as handle:
                fields = handle.readlines()[1].strip().split("\t")
            self.assertEqual(fields, ["23","3000000","A","C","1","0","2","9","9","9","2"])

    def test_cuda_frequency_source_matches_upstream_on_host(self):
        upstream = os.environ.get("IBDMIX_UPSTREAM")
        if not upstream or not shutil.which("g++"):
            self.skipTest("set IBDMIX_UPSTREAM and install g++ for the host differential test")
        upstream = Path(upstream)
        cuda = (ROOT / "f/ibdmix.gpu.cu").read_text()
        scoring = cuda[cuda.index("      U64 alt=0, total=0;"):cuda.index("      sites[rows++]=site;")]
        source = r'''
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "IBDmix/Genotype_Reader.h"
using U64 = unsigned long long;
struct Site { U64 pos; double lod[3]; unsigned char bits, archaic, selected; };
int main() {
  std::mt19937_64 rng(47385);
  const char choices[]={'0','1','2','9'};
  for(int mac=0;mac<=3;++mac) for(int subset=0;subset<=1;++subset) {
    std::stringstream input, sample, other_sample;
    std::string header="chrom\tpos\tref\talt\tarchaic";
    for(int i=0;i<37;++i) header+="\ts"+std::to_string(i);
    input<<header<<'\n';
    if(subset) for(int i=35;i>=0;i-=2) sample<<"s"<<i<<'\n';
    other_sample.str(sample.str());
    std::vector<std::string> bodies;
    for(int row=0;row<2048;++row) {
      std::string body;
      for(int col=0;col<38;++col) {
        if(col) body+='\t';
        char gt=choices[rng()%4];
        if(row%17==0 && col) gt='9';
        if(row%19==0 && col) gt='0';
        if(row%23==0 && col) gt='2';
        body+=gt;
      }
      bodies.push_back(body); input<<"1\t"<<row*83+100<<"\tA\tC\t"<<body<<'\n';
    }
    std::stringstream cpu_mask("1\t100\t1500\n1\t98000\t100000\n"), gpu_mask(cpu_mask.str());
    Genotype_Reader reader(&input,&cpu_mask,.01,.002,2,1e-200,mac);
    int n=reader.initialize(sample);
    Sample_Mapper mapping; std::istringstream header_stream(header);
    if(mapping.initialize(header_stream,other_sample,"")!=n) return 2;
    std::vector<int> indexes(n); for(int i=0;i<n;++i) indexes[i]=mapping.getSample(i);
    Mask_Reader masks(&gpu_mask); LodCalculator calculator(.01,.002,2,1e-200);
    std::string chromosome="1";
    std::vector<unsigned char> genotypes(n);
    for(size_t row=0;row<bodies.size();++row) {
      if(!reader.update()) return 3;
      const char* data=bodies[row].data(); int rows=0;
      U64 pos=row*83+100;
''' + scoring + r'''
      if(site.archaic!=reader.getArchaic()-'0' || site.bits!=reader.getLineFilter()) return 4;
      if((total?(double)alt/total:0)!=reader.getAlleleFrequency()) return 5;
      for(int j=0;j<n;++j) {
        unsigned char gt=genotypes[j], bits=site.bits;
        if(!site.selected && site.archaic==0 && gt==2) bits|=16;
        if(!site.selected && site.archaic==2 && gt==0) bits|=8;
        double observed=gt==9?0:site.lod[gt];
        if(observed!=reader.getLodScore(j) || bits!=(reader.getLineFilter()|reader.getRecoverType(j))) return 6;
      }
    }
  }
  std::cout<<"16384 sites: same frequency, mask/MAC filtering, LOD, and recovery bits\n";
}
'''
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            path, binary = folder / "frequency.cc", folder / "frequency"
            path.write_text(source)
            sources = ["Genotype_Reader.cc", "Sample_Mapper.cc", "Mask_Reader.cc", "lod_calculator.cc"]
            subprocess.run(["g++", "-std=c++17", "-O2", "-I", str(upstream / "include"), str(path), *[str(upstream / "src" / p) for p in sources], "-o", str(binary)], check=True, capture_output=True, text=True)
            result = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
            self.assertIn("16384 sites", result.stdout)

    def test_cuda_stack_source_matches_upstream_on_host(self):
        upstream = os.environ.get("IBDMIX_UPSTREAM")
        if not upstream or not shutil.which("g++"):
            self.skipTest("set IBDMIX_UPSTREAM and install g++ for the host differential test")
        upstream = Path(upstream)
        self.assertTrue((upstream / "src/IBD_Segment.cc").is_file())
        cuda = (ROOT / "f/ibdmix.gpu.cu").read_text()
        node = cuda[cuda.index("struct Node {"):cuda.index("// Append fixed-size pages")]
        stack = cuda[cuda.index("struct State {"):cuda.index("__global__ void init")]
        prefix = r'''
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "IBDmix/IBD_Segment.h"
#include "IBDmix/Segment_Recorders.h"
using U64 = unsigned long long;
#define __device__
static int atomicAdd(int* p, int x) { int before=*p; *p+=x; return before; }
static int atomicMax(int* p, int x) { int before=*p; *p=std::max(*p,x); return before; }
'''
        main = r'''
int main() {
  std::mt19937_64 rng(385767);
  const double values[] = {-20,-4,-2,-.1,0,.1,.3,1,2,4,8};
  for (int trial=0; trial<1024; ++trial) {
    const int samples=3, capacity=4096, page_rows=128;
    const bool exclusive=bool(trial%2);
    const double threshold=(trial%3==0 ? 0 : trial%3==1 ? 1 : 4);
    std::vector<std::vector<Node>> storage(capacity/page_rows, std::vector<Node>(page_rows*samples));
    std::vector<Node*> pages; for(auto& page:storage) pages.push_back(page.data());
    std::vector<State> states(samples);
    std::vector<Event> events(32768); int count=0, error=0;
    std::vector<std::unique_ptr<IBD_Segment>> cpu;
    IBD_Pool pool(1024);
    for(int s=0;s<samples;++s) {
      states[s]={}; states[s].top=states[s].start=states[s].end=-1; states[s].free=0;
      Stack stack{pages.data(),page_rows,states[s],s,samples,&error};
      for(int i=0;i<capacity;++i) stack.node(i).next=i+1<capacity?i+1:-1;
      cpu.emplace_back(new IBD_Segment(std::to_string(s),threshold,&pool,exclusive));
      cpu.back()->add_recorder(std::make_shared<CountRecorder>());
    }
    std::ostringstream expected;
    int rows=32+(rng()%512);
    for(int row=0;row<=rows;++row) for(int s=0;s<samples;++s) {
      double lod=row==rows?-INFINITY:values[rng()%11];
      // Include long declining suffixes and plateaus that exercise replay.
      if(trial%7==0 && row<rows) lod=row%40<3?4:-.1;
      if(trial%11==0 && row<rows) lod=0;
      unsigned char bits=row==rows?0:rng()%32;
      U64 pos=row==rows?0:100+row*83;
      Stack stack{pages.data(),page_rows,states[s],s,samples,&error};
      if(!stack.add(pos,lod,bits,row,threshold,exclusive,events.data(),&count,events.size()) || error) return 2;
      cpu[s]->add_lod("1",pos,lod,bits,expected);
    }
    std::sort(events.begin(),events.begin()+count,[](const Event& a,const Event& b) {
      if(a.row!=b.row) return a.row<b.row;
      if(a.sample!=b.sample) return a.sample<b.sample;
      return a.order<b.order;
    });
    std::ostringstream observed;
    for(int i=0;i<count;++i) {
      const auto& e=events[i];
      observed<<e.sample<<"\t1\t"<<e.start<<'\t'<<e.end<<'\t'<<e.lod;
      for(auto v:e.stats) observed<<'\t'<<v;
      observed<<'\n';
    }
    if(expected.str()!=observed.str()) {
      std::cerr<<"Mismatch trial="<<trial<<"\nCPU:\n"<<expected.str()<<"\nCUDA-source host:\n"<<observed.str();
      return 1;
    }
  }
  std::cout<<"1024 trials, 3 samples each; identical intervals, LODs, and nine statistics\n";
}
'''
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source, binary = folder / "parity.cc", folder / "parity"
            source.write_text(prefix + node + stack + main)
            sources = ["IBD_Stack.cc", "IBD_Segment.cc", "Segment_Recorders.cc", "Genotype_Reader.cc", "Sample_Mapper.cc", "Mask_Reader.cc", "lod_calculator.cc"]
            subprocess.run(["g++", "-std=c++17", "-O2", "-ffp-contract=off", "-I", str(upstream / "include"), str(source), *[str(upstream / "src" / p) for p in sources], "-o", str(binary)], check=True, capture_output=True, text=True)
            result = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
            self.assertIn("1024 trials", result.stdout)


if __name__ == "__main__":
    unittest.main()

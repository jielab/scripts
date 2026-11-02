"""Use the actual AF-only schema, real BGZF/tabix and interrupted checkpoints."""
from contextlib import redirect_stdout
import csv
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"f"))
import allofus_shards as a


def row(chrom="chr1", pos=10001, ref="T", alt="A", all_af="0.013301", mid="0.0", **changes):
    values = dict(zip(a.ANNOVAR_COLUMNS, [chrom,str(pos),str(pos),ref,alt,all_af,"0.040816",
        "0.001887","0.014706","0.033333","0.015748",mid,"0.040816","0.0"]))
    values.update(changes)
    return "\t".join(values[k] for k in a.ANNOVAR_COLUMNS)+"\n"


class AllOfUsShards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root/"hg38_allofus.txt.gz"
        self.plan_file = self.root/"plan.json"
        self.prepared = self.root/"prepared"
        self.normalized = self.root/"normalized"
        self.plan = a.make_plan(self.source,"1,X,Y,MT",100000000)

    def tearDown(self): self.temp.cleanup()

    def fixture(self, rows=None):
        rows = rows if rows is not None else [
            row(), row(alt="C",all_af="0.003448",gvs_eas_af="0.0"),
            row(alt="G",all_af="0.001478",mid=""),
            row(pos=10108,ref="C",alt="CA",all_af="0.412381",mid=""),
            row(pos=10110,ref="A",alt="AAC"), row(pos=10110,ref="ACCCT",alt="A"),
            row(pos=100000001,ref="A",alt="C",all_af="0.0",mid=""),
            row(chrom="chrX",pos=20,ref="C",alt="A"),
            row(chrom="chrY",pos=30,ref="G",alt="T"),
            row(chrom="chrM",pos=10,ref="A",alt="G")]
        with gzip.open(self.source,"wt") as handle:
            handle.write("#"+"\t".join(a.ANNOVAR_COLUMNS)+"\n"+"".join(rows))
        a.write_plan(self.plan,self.plan_file)

    def prepare(self):
        with redirect_stdout(io.StringIO()): return a.prepare(self.plan_file,self.prepared,progress_every=0)

    def normalize(self):
        results=[]
        for s in self.plan["shards"]:
            with redirect_stdout(io.StringIO()):
                results.append(a.normalize_shard(self.prepared,s["id"],self.normalized/s["id"]))
        return results

    def table(self,path):
        with gzip.open(path,"rt") if str(path).endswith(".gz") else Path(path).open() as handle:
            reader=csv.DictReader((line[1:] if line.startswith("#") else line for line in handle),delimiter="\t")
            return list(reader)

    def target_file(self, text):
        p=self.root/"targets.tsv";p.write_text("site_id\tchrom\tpos\tref\talts\tbuild\n"+text)
        return p

    def test_plan_is_metadata_only_and_includes_sex_chromosomes(self):
        with patch.object(Path,"open",side_effect=AssertionError("source opened")), patch.object(Path,"stat",side_effect=AssertionError("source stat")):
            plan=a.make_plan("/not/mounted/huge.gz")
        self.assertTrue(300 < len(plan["shards"]) < 400)
        self.assertTrue({"X","Y","MT"} <= set(plan["chromosomes"]))
        self.assertEqual(plan["source_columns"],a.ANNOVAR_COLUMNS)
        self.assertNotIn("MAX",plan["populations"])

    def test_prepare_checkpoints_preserve_blank_fields_and_position_boundaries(self):
        self.fixture(); result=self.prepare()
        self.assertEqual(result["counts"]["source_rows"],10)
        s=self.plan["shards"][0]
        rows=self.table(a.raw_path(self.prepared,s))
        self.assertEqual(rows[2]["gvs_mid_af"],"")
        self.assertEqual({r["Start"] for r in rows[:3]},{"10001"})
        second=self.plan["shards"][1]
        self.assertEqual(self.table(a.raw_path(self.prepared,second))[0]["Start"],"100000001")
        self.assertEqual(len(list((self.prepared/"receipts").glob("*.json"))),len(self.plan["shards"]))

    def test_completed_preparation_reuses_without_reopening_original_gzip(self):
        self.fixture();first=self.prepare()
        times={p.name:p.stat().st_mtime_ns for p in (self.prepared/"raw").glob("*")}
        with patch.object(a,"source_rows",side_effect=AssertionError("reopened")):
            self.assertEqual(self.prepare(),first)
        self.assertEqual(times,{p.name:p.stat().st_mtime_ns for p in (self.prepared/"raw").glob("*")})

    def test_interrupted_preparation_replays_prefix_but_keeps_completed_windows(self):
        self.fixture(); original=a.source_rows
        def interrupted(path):
            for n, values in original(path):
                yield n,values
                if values[1]=="100000001": raise RuntimeError("synthetic interruption")
        with patch.object(a,"source_rows",side_effect=interrupted):
            with self.assertRaisesRegex(RuntimeError,"synthetic interruption"): self.prepare()
        first_raw=a.raw_path(self.prepared,self.plan["shards"][0]); before=first_raw.stat().st_mtime_ns
        self.assertEqual(json.loads((self.prepared/"prepare_manifest.json").read_text())["status"],"failed")
        final=self.prepare()
        self.assertEqual(first_raw.stat().st_mtime_ns,before)
        self.assertEqual(final["counts"]["source_rows"],10)
        self.assertEqual(final["counts"]["reused_nonempty_windows"],1)
        self.assertIn("Replay",final["restart"])

    def test_source_change_cannot_silently_reuse_cache(self):
        self.fixture();self.prepare()
        with self.source.open("ab") as handle: handle.write(b"changed")
        with self.assertRaisesRegex(ValueError,"different source/plan"): self.prepare()

    def test_selected_chromosome_is_checkpointed_before_reading_unselected_remainder(self):
        self.fixture([row(), row(chrom="chr2")]); original=a.source_rows
        def interrupted(path):
            for n, values in original(path):
                yield n,values
                if values[0]=="chr2": raise RuntimeError("unselected-tail interruption")
        with patch.object(a,"source_rows",side_effect=interrupted):
            with self.assertRaisesRegex(RuntimeError,"unselected-tail"): self.prepare()
        first=self.plan["shards"][0]
        receipt=json.loads(a.receipt_path(self.prepared,first).read_text())
        self.assertEqual(receipt["rows"],1)
        self.assertEqual(receipt["status"],"complete")
        self.assertEqual(self.prepare()["counts"]["reused_nonempty_windows"],1)

    def test_order_failure_keeps_existing_checkpoint_and_does_not_claim_complete(self):
        self.fixture([row(pos=10001),row(pos=100000001),row(pos=10000)])
        with self.assertRaisesRegex(ValueError,"sorted"):self.prepare()
        m=json.loads((self.prepared/"prepare_manifest.json").read_text())
        self.assertEqual(m["status"],"failed")
        self.assertTrue(a.receipt_path(self.prepared,self.plan["shards"][0]).is_file())
        with self.assertRaisesRegex(ValueError,"Preparation must complete"):
            a.normalize_shard(self.prepared,self.plan["shards"][0]["id"],self.root/"bad")

    def test_repeated_chromosome_block_is_rejected(self):
        self.fixture([row(),row(chrom="chrX"),row(pos=20000)])
        with self.assertRaisesRegex(ValueError,"contiguous"):self.prepare()

    def test_raw_shard_corruption_fails_inventory_validation(self):
        self.fixture();self.prepare()
        with a.raw_path(self.prepared,self.plan["shards"][0]).open("ab") as handle:handle.write(b"x")
        with self.assertRaisesRegex(ValueError,"has changed"):self.normalize()

    def test_normalize_keeps_eight_populations_wide_and_max_separate(self):
        self.fixture();self.prepare();results=self.normalize()
        first=self.normalized/self.plan["shards"][0]["id"]
        alleles=self.table(first/"wide_alleles.tsv.gz")
        self.assertEqual(len(alleles),3)
        self.assertEqual(set(alleles[0]),set(a.WIDE_FIELDS))
        self.assertEqual(alleles[2]["gvs_mid_af"],"NA")
        self.assertEqual(alleles[0]["gvs_max_af"],"0.040816")
        self.assertNotIn("population",alleles[0])
        self.assertNotIn("AC",alleles[0]);self.assertNotIn("AN",alleles[0])
        site=self.table(first/"wide_loci.tsv.gz")[0]
        self.assertEqual(site["n_positive_ALT_ALL"],"3")
        self.assertEqual(site["n_known_AF_MID"],"2")
        self.assertEqual(results[0]["counts"]["non_snv_rows"],3)
        self.assertEqual(results[0]["populations"]["ALL"]["loci_ge3_positive_ALT"],1)
        self.assertEqual(set(results[0]["populations"]),set(a.POPS))

    def test_actual_tabix_queries_exact_locus_and_zero_is_preserved(self):
        import pysam
        self.fixture();self.prepare();self.normalize()
        first=self.normalized/self.plan["shards"][0]["id"]
        with pysam.TabixFile(str(first/"wide_alleles.tsv.gz")) as t:
            self.assertEqual(len(list(t.fetch("1",10000,10001))),3)
            self.assertEqual(list(t.fetch("1",10001,10002)),[])
        second=self.normalized/self.plan["shards"][1]["id"]
        self.assertEqual(self.table(second/"wide_alleles.tsv.gz")[0]["gvs_all_af"],"0.0")

    def test_normalized_empty_windows_have_no_frequency_zero_rows(self):
        self.fixture();self.prepare();results=self.normalize()
        empty=[r for r in results if r["data_status"]=="no_reported_snv_rows"]
        self.assertTrue(empty)
        for m in empty:
            path=self.normalized/m["shard"]["id"]
            self.assertEqual(self.table(path/"wide_alleles.tsv.gz"),[])
            self.assertEqual(m["index_status"],"no_index_for_empty_table")

    def test_identical_duplicates_deduplicate_and_conflicting_metrics_fail(self):
        self.fixture([row(),row(),row(alt="C")]);self.prepare();r=self.normalize()[0]
        self.assertEqual(r["counts"]["identical_duplicate_rows"],1)
        self.assertEqual(r["counts"]["unique_snv_alleles"],2)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);s=p/"source.gz"
            with gzip.open(s,"wt") as h:h.write("#"+"\t".join(a.ANNOVAR_COLUMNS)+"\n"+row()+row(all_af="0.8"))
            plan=a.make_plan(s,"1",300000000);a.write_plan(plan,p/"plan.json")
            with redirect_stdout(io.StringIO()):a.prepare(p/"plan.json",p/"prep",0)
            with self.assertRaisesRegex(ValueError,"Conflicting duplicate"):
                a.normalize_shard(p/"prep",plan["shards"][0]["id"],p/"norm")

    def test_snv_ref_conflict_fails_but_mixed_indel_refs_do_not(self):
        self.fixture([row(),row(ref="G",alt="A")]);self.prepare()
        with self.assertRaisesRegex(ValueError,"Conflicting canonical SNV REF"):self.normalize()

    def test_summary_requires_every_window_and_keeps_chromosome_metrics(self):
        self.fixture();self.prepare();self.normalize()
        result=a.summarize(self.plan_file,self.normalized,self.root/"summary")
        self.assertEqual(result["counts"]["source_rows"],10)
        self.assertEqual(result["counts"]["unique_snv_alleles"],7)
        rows=self.table(self.root/"summary/population_window_summary.tsv")
        self.assertEqual({r["chrom"] for r in rows},{"1","X","Y","MT"})
        self.assertEqual(len(rows),len(self.plan["shards"])*8)
        missing=self.normalized/self.plan["shards"][-1]["id"]/"manifest.json";missing.unlink()
        with self.assertRaises(FileNotFoundError):a.summarize(self.plan_file,self.normalized,self.root/"summary2")

    def test_indexed_candidate_query_has_additional_alts_missing_and_no_source_scan(self):
        self.fixture();self.prepare();self.normalize()
        targets=self.target_file("s1\t1\t10001\tT\tA,C\t38\ns2\t1\t100000001\tA\tC,G\t38\ns3\tY\t40\tC\tA\t38\n")
        with patch.object(a,"source_rows",side_effect=AssertionError("sequential source scan")):
            result=a.query(self.prepared,self.normalized,targets,self.root/"query")
        rows=self.table(self.root/"query/target_alleles.tsv.gz")
        s1={r["alt"]:r for r in rows if r["site_id"]=="s1"}
        self.assertEqual(set(s1),{"A","C","G"})
        self.assertEqual(s1["G"]["additional_reported_alt"],"1")
        self.assertEqual(s1["G"]["observed_in_target"],"0")
        s2={r["alt"]:r for r in rows if r["site_id"]=="s2"}
        self.assertEqual(s2["C"]["gvs_all_af"],"0.0")
        self.assertEqual(s2["G"]["gvs_all_af"],"NA")
        self.assertEqual(s2["G"]["match_status"],"not_in_summary")
        self.assertEqual(result["queried_windows"],3)

    def test_query_ref_and_build_conflicts_do_not_leak_scores(self):
        self.fixture();self.prepare();self.normalize()
        targets=self.target_file("wrong_ref\t1\t10001\tA\tC,G\t38\n")
        a.query(self.prepared,self.normalized,targets,self.root/"query")
        rows=self.table(self.root/"query/target_alleles.tsv.gz")
        self.assertTrue(all(r["match_status"]=="reference_mismatch" and r["gvs_all_af"]=="NA" for r in rows))
        targets=self.target_file("wrong_build\t1\t10001\tT\tC\t37\n")
        with self.assertRaisesRegex(ValueError,"build conflicts"):a.query(self.prepared,self.normalized,targets,self.root/"bad_build")

    def test_query_all_population_af_missing_is_not_zero_or_frequency_evidence(self):
        self.fixture([row(**{column:"" for column in a.AF_COLUMNS})]);self.prepare();self.normalize()
        targets=self.root/"targets.tsv"
        targets.write_text("site_id\tchrom\tpos\tref\talts\ttarget_build\ns1\t1\t10001\tT\tA\tGRCh38\n")
        a.query(self.prepared,self.normalized,targets,self.root/"query")
        allele=self.table(self.root/"query/target_alleles.tsv.gz")[0]
        site=self.table(self.root/"query/target_sites.tsv")[0]
        self.assertEqual(allele["match_status"],"matched_metrics_missing")
        self.assertEqual(allele["evidence"],"metrics_missing")
        self.assertEqual(allele["gvs_all_af"],"NA")
        self.assertEqual(site["n_known_AF_ALL"],"0")
        self.assertEqual(site["n_positive_ALT_ALL"],"0")
        self.assertEqual(site["evidence"],"metrics_missing")
        targets.write_text(targets.read_text().replace("GRCh38","GRCh37"))
        with self.assertRaisesRegex(ValueError,"build conflicts"):
            a.query(self.prepared,self.normalized,targets,self.root/"bad_target_build")

    def test_plan_cli_only_writes_metadata_without_existing_source(self):
        p=self.root/"cli_plan.json"
        result=subprocess.run([sys.executable,str(ROOT/"f/allofus_shards.py"),"plan","--source",str(self.root/"missing.gz"),
            "--chromosomes","1,X,Y","--out",str(p)],capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(p.read_text())["build"],"GRCh38")


if __name__ == "__main__":unittest.main()

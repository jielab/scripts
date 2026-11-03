"""Public-shell planning and a real small All of Us DAG through batch workers.

The HPC fixture passes a FIFO as a source: opening it for genomic reads on the
submission host would block. Fake bsub records submitted scripts but never runs
their payloads. Local scientific tests call the Python planner directly, leaving
the public shell's site-specific systemd/cgroup launcher to its own checks.
"""
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
ENTRY = ROOT/"allele3.sh"
PLANNER = ROOT/"f/batch_workflows.py"
HEADER = ["#Chr", "Start", "End", "Ref", "Alt", "gvs_all_af", "gvs_max_af",
          "gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af", "gvs_mid_af",
          "gvs_oth_af", "gvs_sas_af"]

FAKE_LSF = r'''import json,os,pathlib,sys
path=pathlib.Path(os.environ["ALLELE3_WORKFLOW_FAKE_LSF"])
state=json.loads(path.read_text());name=pathlib.Path(sys.argv[0]).name
args=sys.argv[1:];state.setdefault("calls",[]).append([name,*args])
if name=="bsub":
    payload=sys.stdin.read();jid=state.get("next",700);state["next"]=jid+1
    state["bsub_n"]=state.get("bsub_n",0)+1
    state.setdefault("jobs",{})[str(jid)]={"payload":payload,"state":"PEND",
        "dependency":args[args.index("-w")+1] if "-w" in args else ""}
    path.write_text(json.dumps(state));print("Job <"+str(jid)+"> is submitted to queue <ser>.")
elif name=="bjobs":
    path.write_text(json.dumps(state))
    for jid in args:
        if jid in state.get("jobs",{}):print(jid,state["jobs"][jid]["state"])
elif name=="bwait":
    path.write_text(json.dumps(state))
else:
    path.write_text(json.dumps(state));sys.exit(2)
'''


def clean_environment():
    exact = {"submit", "SUBMIT", "DATE_SCRATCH", "ANALYSIS_ROOT", "SCRATCH_ANALYSIS_ROOT",
             "SCRIPTS_ROOT", "SHARED_PYTHON", "jobs_N", "SUBMIT_WAIT"}
    prefixes = ("ALLELE_", "PCA3D_", "HERV_", "DIR_", "OS_", "LSB_", "WORKFLOW_")
    env = {k:v for k,v in os.environ.items() if k not in exact and not k.startswith(prefixes)}
    # Preserve a real allocation for tiny local fixtures on a compute node.
    if os.environ.get("LSB_JOBID"):
        env["LSB_JOBID"] = os.environ["LSB_JOBID"]
    env.update(PYTHONDONTWRITEBYTECODE="1", ALLELE_PLANNER_PYTHON=sys.executable)
    return env


def actual_os_config():
    """Use the repository shared script, or the unchanged review-copy sibling."""
    candidates = []
    if os.environ.get("ALLELE_TEST_OS_CONFIG"):
        candidates.append(Path(os.environ["ALLELE_TEST_OS_CONFIG"]))
    candidates.append(ROOT.parent/"0f/0os.sh")
    candidates.extend(p/"hpc_reference/0f/0os.sh" for p in ROOT.parents)
    return next((p for p in candidates if p.is_file()), None)


class BatchWorkflowsIntegration(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="allele3-workflow-qa-")
        self.path = Path(self.temp.name)
        self.env = clean_environment()

    def tearDown(self):
        self.temp.cleanup()

    def run_ok(self, argv, env=None, timeout=30):
        result = subprocess.run([str(a) for a in argv], env=env or self.env,
                                text=True, capture_output=True, timeout=timeout)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        return result

    def shell_env(self):
        config = actual_os_config()
        if config is None:
            self.skipTest("The real shared scripts/0f/0os.sh is required; set ALLELE_TEST_OS_CONFIG")
        return {**self.env, "ALLELE_OS_CONFIG":str(config)}

    def table(self, path):
        with gzip.open(path, "rt") if str(path).endswith(".gz") else path.open() as handle:
            return list(csv.DictReader((line[1:] if line.startswith("#") else line for line in handle), delimiter="\t"))

    def test_atlas_missing_inputs_report_complete_paths_before_submission(self):
        for mode in ('local', 'hpc'):
            with self.subTest(mode=mode):
                data=self.path/mode
                jobs=self.path/(mode+'-jobs')
                expected=data/'annot/allele3/alphagenome/hg38/alphagenome_variant_impact_score_snvs.tsv.gz'
                result=subprocess.run(['bash',str(ENTRY),'batch','atlas','--submit',mode,
                    '--write-cmd','--targets',str(self.path/'targets.tsv'),
                    '--out',str(self.path/'out'),'--job-dir',str(jobs)],
                    env={**self.shell_env(),'DIR_DATA':str(data)},text=True,capture_output=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('Required Atlas file missing',result.stderr)
                for path in (expected, str(expected)+'.tbi', str(expected)+'.csi'):
                    self.assertIn(str(path),result.stderr)
                self.assertFalse(jobs.exists())

    def test_atlas_plan_resolves_config_relative_paths_and_data_override_index(self):
        import importlib.util
        if importlib.util.find_spec('pysam') is None:
            self.skipTest('pysam is needed for the real Atlas worker')
        import pysam
        configdir=self.path/'config';configdir.mkdir()
        plain=configdir/'avi.tsv';plain.write_text('1\t10\tA\tC\t0.1\t2\n')
        data=configdir/'avi.tsv.gz'
        pysam.tabix_compress(str(plain),str(data),force=True)
        pysam.tabix_index(str(data),seq_col=0,start_col=1,end_col=1,zerobased=False)
        with (ROOT/'config/alphagenome_resources.tsv').open() as handle:
            reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
        rows[0].update(path='avi.tsv.gz',index_path='old-index.tbi')
        descriptor=configdir/'resources.tsv'
        with descriptor.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
        targets=self.path/'targets.tsv';targets.write_text('chrom\tpos\tref\talts\n1\t10\tA\tC\n')
        jobs=self.path/'jobs';out=self.path/'out'
        # Explicit data override must clear a configured stale index in both planner and worker.
        self.run_ok([sys.executable,PLANNER,'--entry',ENTRY,'--mode','local','--','batch','atlas',
            '--targets',targets,'--avi-file',data,'--atlas-resources',descriptor,
            '--python',sys.executable,'--out',out],
            {**self.env,'ALLELE_JOB_DIR':str(jobs),'ALLELE_JOB_MEMORY_GIB':'1'},timeout=60)
        manifest=json.loads((out/'manifest.json').read_text())
        self.assertIn(manifest['status'],('complete','complete_with_exclusions'))
        frozen=(jobs/'plan.json').read_text()
        self.assertNotIn('old-index.tbi',frozen)
        self.assertIn(str(data)+'.tbi',frozen)
        self.assertEqual(json.loads((jobs/'submission.json').read_text())['status'],'complete')

    def test_public_hpc_write_cmd_has_hundreds_of_workers_and_never_opens_fifo_source(self):
        source = self.path/"source-that-must-not-be-opened.gz"
        os.mkfifo(source)
        out, jobs = self.path/"cache", self.path/"job directory"
        result = self.run_ok(["bash", ENTRY, "batch", "allofus", "--submit", "hpc",
            "--write-cmd", "--source", source, "--out", out, "--job-dir", jobs,
            "--max-running", "200"], self.shell_env(), timeout=30)
        plan = json.loads((jobs/"plan.json").read_text())
        windows = json.loads((jobs/"allofus_windows.json").read_text())
        self.assertEqual(plan["mode"], "hpc")
        self.assertEqual(plan["max_running"], 200)
        self.assertGreater(len(windows["shards"]), 300)
        self.assertEqual(len(list((jobs/"cmd").glob("*.cmd"))), len(plan["jobs"]))
        self.assertGreater(len(plan["jobs"]), 300)
        byid = {j["id"]:j for j in plan["jobs"]}
        self.assertEqual(byid["allofus_prepare"]["deps"], [])
        normalized = [j for j in plan["jobs"] if j["id"].startswith("aou_chr")]
        self.assertEqual(len(normalized), len(windows["shards"]))
        self.assertTrue(all(j["deps"] == ["allofus_prepare"] for j in normalized))
        self.assertTrue(all(len(j["deps"]) <= 32 for j in plan["jobs"]))
        self.assertTrue({"X", "Y", "MT"} <= {s["chrom"] for s in windows["shards"]})
        self.assertEqual(byid["allofus_prepare"]["outputs"], [str(out/"prepared/prepare_manifest.json")])
        self.assertFalse(out.exists(), "Write-only planning ran a scientific worker")
        self.assertEqual(list((jobs/"runs").glob("*.json")), [])
        self.assertIn("mode=hpc", result.stdout)
        self.assertIn("no login-node payload", (jobs/"cmd/allofus_prepare.cmd").read_text())
        self.assertTrue((jobs/"submit.sh").is_file())

    def test_mode_date_order_and_hpc_default_inputs_follow_the_shared_os_script(self):
        cases = [(["--submit", "hpc", "--date-scratch", "2026-11-02"], "2026-11-02"),
                 (["--date-scratch", "2026-11-03", "--submit", "hpc"], "2026-11-03"),
                 (["--date-scratch=2026-11-04", "--submit=hpc"], "2026-11-04")]
        for i, (flags, date) in enumerate(cases):
            with self.subTest(order=flags):
                out, jobs = self.path/f"out{i}", self.path/f"jobs{i}"
                # Simulate switching from the shared local-mode environment.
                env = {**self.shell_env(), "submit":"local", "SUBMIT":"local",
                       "DIR_DATA":"/mnt/f", "DIR_WORK":"/mnt/d", "DIR_SCRATCH":"/mnt/d",
                       "DIR_SOFTWARE":"/mnt/d/software", "DATE_SCRATCH":"2026-10-31"}
                result = self.run_ok(["bash", ENTRY, *flags, "batch", "allofus", "--write-cmd",
                    "--chromosomes", "1", "--window-bp", "250000000", "--out", out, "--job-dir", jobs], env)
                plan = json.loads((jobs/"plan.json").read_text())
                windows = json.loads((jobs/"allofus_windows.json").read_text())
                roots = plan["provenance"]["roots"]
                self.assertEqual(roots["data"], "/work/sph-huangj/data")
                self.assertEqual(roots["work"], "/work/sph-huangj")
                self.assertEqual(roots["analysis"], f"/scratch/{date}/sph-huangj/analysis/allele3")
                self.assertEqual(windows["source"], "/work/sph-huangj/data/annot/allele3/allofus/hg38_allofus.txt.gz")
                for job in plan["jobs"]:
                    self.assertEqual(job["env"]["DIR_SCRATCH"], f"/scratch/{date}/sph-huangj")
                    self.assertEqual(job["env"]["ALLELE_VCF_DIR"], "/work/sph-huangj/data/annot/gnomad/sites")
                    self.assertEqual(job["env"]["ALLELE_ARCHAIC_ROOT"], "/work/sph-huangj/data/gen")
                    self.assertEqual(job["argv"][0], "/work/sph-huangj/.conda/envs/allele3-evolution/bin/python")
                    resolved_paths = {k:v for k,v in job["env"].items()
                                      if k not in {"PATH", "SCRIPTS_ROOT", "ALLELE_OS_CONFIG"} and v.startswith("/")}
                    # Code binds to the submission checkout, which can be on
                    # /mnt during WSL planning tests; data roots must be HPC.
                    self.assertEqual(job["env"]["ALLELE_OS_CONFIG"], str(actual_os_config()))
                    self.assertNotIn("/mnt/", json.dumps(resolved_paths))
                self.assertIn("date_scratch="+date, result.stdout)
                self.assertFalse(out.exists())

    def test_public_hpc_submission_uses_numeric_dependencies_and_does_not_resubmit(self):
        bindir = self.path/"fake_lsf";bindir.mkdir()
        state = self.path/"scheduler.json";state.write_text("{}")
        for name in ("bsub", "bjobs", "bwait"):
            script = bindir/name
            script.write_text("#!"+sys.executable+"\n"+FAKE_LSF);script.chmod(0o755)
        source = self.path/"source.gz";os.mkfifo(source)
        targets = self.path/"targets.tsv"
        targets.write_text("site_id\tchrom\tpos\tref\talts\tbuild\ns1\t1\t10001\tT\tA,C\t38\n")
        out, jobs = self.path/"out", self.path/"jobs"
        env = {**self.shell_env(), "PATH":str(bindir)+os.pathsep+self.env["PATH"],
               "ALLELE3_WORKFLOW_FAKE_LSF":str(state)}
        command = ["bash", ENTRY, "batch", "allofus", "--submit", "hpc", "--wait", "FALSE",
            "--source", source, "--targets", targets, "--chromosomes", "1", "--window-bp", "250000000",
            "--out", out, "--job-dir", jobs]
        self.run_ok(command, env)
        first = json.loads((jobs/"submission.json").read_text())
        fake = json.loads(state.read_text())
        self.assertEqual(fake["bsub_n"], 4)
        self.assertEqual([fake["jobs"][str(jid)]["dependency"] for jid in range(700,704)],
                         ["", "done(700)", "done(701)", "done(702)"])
        self.assertTrue(all(row["status"] == "submitted" for row in first["jobs"].values()))
        self.assertFalse(out.exists())
        self.run_ok(command, env)
        second = json.loads((jobs/"submission.json").read_text())
        self.assertEqual(json.loads(state.read_text())["bsub_n"], 4)
        self.assertEqual({k:r["job_id"] for k,r in first["jobs"].items()},
                         {k:r["job_id"] for k,r in second["jobs"].items()})

    def test_local_small_allofus_pipeline_queries_extra_alts_na_indels_and_reuses_outputs(self):
        if shutil.which("bsub", path=self.env.get("PATH")):
            self.skipTest("Local payloads are deliberately refused on an LSF submission host")
        source = self.path/"hg38_allofus.txt.gz"
        rows = [
            ["chr1","10001","10001","T","A","0.013301","0.040816","0.001887","0.014706","0.033333","0.015748","0.0","0.040816","0.0"],
            ["chr1","10001","10001","T","C","0.003448","0.005906","0.001887","0.0","0.0","0.005906","0.0","0.0","0.0"],
            ["chr1","10001","10001","T","G","0.001478","0.002953","0.0","0.002953","0.0","0.002953","","0.0","0.0"],
            ["chr1","10108","10108","C","CA","0.412381","0.5","0.378151","0.430769","0.375","0.422034","","0.414286","0.5"],
            ["chr1","10110","10110","A","AAC","0.000165","0.000199","0.000183","0.0","0.0","0.000145","0.0","0.0","0.0"],
            ["chr1","10110","10110","ACCCT","A","0.000237","0.000546","0.000173","0.000304","0.00051","0.000189","0.0","0.000254","0.000546"],
            ["chr1","10200","10200","G","A","0.0","0.0","0.0","0.0","0.0","0.0","","0.0","0.0"],
        ]
        with gzip.open(source,"wt",newline="") as handle:
            writer=csv.writer(handle,delimiter="\t",lineterminator="\n");writer.writerow(HEADER);writer.writerows(rows)
        targets=self.path/"targets.tsv"
        targets.write_text("site_id\tchrom\tpos\tref\talts\tbuild\n"
            "third\t1\t10001\tT\tA,C\t38\nzero\t1\t10200\tG\tA,C\t38\nmissing\t1\t20000\tA\tG\t38\n")
        out, jobs = self.path/"cache", self.path/"jobs"
        env={**self.env,"ALLELE_JOB_DIR":str(jobs),"ALLELE_JOB_MEMORY_GIB":"1","jobs_N":"2"}
        command=[sys.executable,PLANNER,"--entry",ENTRY,"--mode","local","--","batch","allofus",
            "--source",source,"--targets",targets,"--chromosomes","1","--window-bp","250000000",
            "--python",sys.executable,"--out",out]
        result=self.run_ok(command,env,timeout=60)
        frozen=json.loads((jobs/"plan.json").read_text())
        windows=json.loads((jobs/"allofus_windows.json").read_text())
        state=json.loads((jobs/"submission.json").read_text())
        self.assertEqual(state["status"],"complete",result.stdout)
        self.assertEqual(len(frozen["jobs"]),4);self.assertEqual(len(windows["shards"]),1)
        self.assertTrue(all(r["status"]=="succeeded" for r in state["jobs"].values()))
        prepare=json.loads((out/"prepared/prepare_manifest.json").read_text())
        summary=json.loads((out/"summary/manifest.json").read_text())
        self.assertEqual(prepare["status"],"complete")
        self.assertEqual(prepare["counts"]["source_rows"],7)
        self.assertEqual(summary["counts"]["unique_snv_alleles"],4)
        self.assertEqual(summary["counts"]["non_snv_rows"],3)
        shard=out/"normalized"/windows["shards"][0]["id"]
        self.assertTrue((shard/"wide_alleles.tsv.gz.tbi").is_file())
        alleles=self.table(out/"targets/target_alleles.tsv.gz")
        third={r["alt"]:r for r in alleles if r["site_id"]=="third"}
        self.assertEqual(set(third),{"A","C","G"})
        self.assertEqual(third["G"]["additional_reported_alt"],"1")
        self.assertEqual(third["G"]["gvs_mid_af"],"NA")
        self.assertNotIn("AC",third["G"]);self.assertNotIn("AN",third["G"])
        zero={r["alt"]:r for r in alleles if r["site_id"]=="zero"}
        self.assertEqual(zero["A"]["gvs_all_af"],"0.0")
        self.assertEqual(zero["C"]["gvs_all_af"],"NA")
        self.assertEqual(zero["C"]["match_status"],"not_in_summary")
        site=next(r for r in self.table(out/"targets/target_sites.tsv") if r["site_id"]=="third")
        self.assertEqual(site["n_positive_ALT_ALL"],"3")
        self.assertEqual(site["n_known_AF_MID"],"2")
        # A second identical invocation must reuse successful scientific outputs.
        before={str(p.relative_to(out)):p.stat().st_mtime_ns for p in out.rglob("*") if p.is_file()}
        receipts={p.name:p.read_text() for p in (jobs/"runs").glob("*.json")}
        self.run_ok(command,env,timeout=60)
        after={str(p.relative_to(out)):p.stat().st_mtime_ns for p in out.rglob("*") if p.is_file()}
        self.assertEqual(before,after)
        self.assertEqual(receipts,{p.name:p.read_text() for p in (jobs/"runs").glob("*.json")})


if __name__=="__main__":
    unittest.main()

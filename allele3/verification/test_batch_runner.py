"""Real subprocess/DAG fixtures; the fake scheduler never runs scientific jobs."""
from contextlib import redirect_stdout, redirect_stderr
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"f"))
import batch_runner as br

FAKE_SCHEDULER='''import json,os,pathlib,sys,time
path=pathlib.Path(os.environ["ALLELE3_FAKE_SCHEDULER"])
s=json.loads(path.read_text()); name=pathlib.Path(sys.argv[0]).name
args=sys.argv[1:]; s.setdefault("calls",[]).append([name,*args])
def save(): path.write_text(json.dumps(s))
if name=="bsub":
    cmd=sys.stdin.read(); s["bsub_n"]=s.get("bsub_n",0)+1
    if s.get("reject_n")==s["bsub_n"]:
        save(); print("Queue request rejected",file=sys.stderr); sys.exit(1)
    jid=s.get("next",100); s["next"]=jid+1
    s.setdefault("jobs",{})[str(jid)]={"state":"PEND","cmd":cmd,"dependency":args[args.index("-w")+1] if "-w" in args else ""}
    save()
    if s.get("opaque_n")==s["bsub_n"]: print("accepted"); sys.exit(0)
    print("Job <"+str(jid)+"> is submitted to queue <ser>.",flush=True)
    if s.get("slow_n")==s["bsub_n"]: time.sleep(1)
elif name=="bjobs":
    save()
    for jid in [v for v in args if v.isdigit()]:
        if jid in s.get("jobs",{}): print(jid,s["jobs"][jid]["state"])
elif name=="bmod":
    s["bmod_n"]=s.get("bmod_n",0)+1
    if s.get("reject_bmod_n")==s["bmod_n"]:
        save(); print("bmod rejected",file=sys.stderr); sys.exit(1)
    jid=args[-1]; s["jobs"][jid]["dependency"]=args[args.index("-w")+1]
    save(); print("Parameters of job <"+jid+"> changed.")
elif name=="bwait":
    save(); sys.exit(0)
else: sys.exit(2)
'''


class BatchRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix="allele3-batch-qa-")
        self.path=Path(self.temp.name); self.out=self.path/"job directory"
        self.bin=self.path/"scheduler"; self.bin.mkdir()
        for name in ("bsub","bjobs","bmod","bwait"):
            p=self.bin/name; p.write_text("#!"+sys.executable+"\n"+FAKE_SCHEDULER); p.chmod(0o755)
        self.scheduler=self.path/"scheduler.json"; self.scheduler.write_text("{}")
        self.env=patch.dict(os.environ,{"ALLELE3_FAKE_SCHEDULER":str(self.scheduler)})
        self.env.start()

    def tearDown(self):
        self.env.stop(); self.temp.cleanup()

    def job(self,name,deps=(),payload=None,outputs=()):
        if payload is None: payload=[sys.executable,"-c","pass"]
        return dict(id=name,argv=payload,deps=list(deps),cwd=str(self.path),
                    resources=dict(queue="ser",cores=1,memory_gib="0.01",hours=1),outputs=list(map(str,outputs)))

    def plan(self,jobs,**extra):
        return dict(name="fixture",jobs=jobs,driver_python=sys.executable,max_running=2,**extra)

    def run_plan(self,plan,**kw):
        with redirect_stdout(io.StringIO()):
            return br.run_plan(plan,self.out,bsub=str(self.bin/"bsub"),bjobs=str(self.bin/"bjobs"),
                               bmod=str(self.bin/"bmod"),bwait=str(self.bin/"bwait"),**kw)

    def state(self): return json.loads((self.out/"submission.json").read_text())
    def fake(self): return json.loads(self.scheduler.read_text())
    def change_fake(self,**changes):
        d=self.fake();d.update(changes);self.scheduler.write_text(json.dumps(d))

    def test_write_205_cmd_files_does_not_execute_or_submit(self):
        sentinel=self.path/"payload_ran"
        jobs=[self.job("part"+str(i),payload=[sys.executable,"-c",f"open({str(sentinel)!r},'w').write('bad')"]) for i in range(205)]
        plan=self.plan(jobs); plan["max_running"]=200
        result=self.run_plan(plan,write_only=True)
        self.assertEqual(result["status"],"planned")
        self.assertEqual(len(list((self.out/"cmd").glob("*.cmd"))),205)
        self.assertFalse(sentinel.exists());self.assertEqual(self.fake(),{})
        frozen=json.loads((self.out/"plan.json").read_text())
        self.assertEqual(frozen["jobs"][200]["throttle_after"],"part0")
        self.assertTrue((self.out/"submit.sh").exists())

    def test_resources_memory_scope_short_slots_and_worker_allocation(self):
        job=self.job("a");job["resources"].update(queue="short",cores=2,memory_gib=160,hours=3.25,gpu=1)
        self.run_plan(self.plan([job]),write_only=True)
        text=(self.out/"cmd/a.cmd").read_text()
        self.assertIn("#BSUB -n 40",text);self.assertIn("mem=160GB/host",text)
        self.assertIn("#BSUB -W 3:15",text);self.assertIn('num=1:mode=exclusive_process',text)
        self.assertIn("ALLELE_ALLOCATED_MEMORY_GIB=160",text)
        self.assertIn("ALLELE_ALLOCATED_CORES=40",text)
        self.assertIn("ALLELE_BATCH_WORKER=1",text)
        result=subprocess.run(["bash",str(self.out/"cmd/a.cmd")],capture_output=True,text=True,
                              env={k:v for k,v in os.environ.items() if k!="LSB_JOBID"})
        self.assertNotEqual(result.returncode,0);self.assertIn("no login-node payload",result.stderr)

    def test_actual_numeric_dependencies_and_lane_limit(self):
        jobs=[self.job("a"),self.job("b"),self.job("c",["a"]),self.job("d")]
        state=self.run_plan(self.plan(jobs))
        self.assertEqual(state["jobs"]["c"]["dependency"],"done(100)")
        self.assertEqual(state["jobs"]["d"]["dependency"],"ended(101)")
        self.assertEqual(self.fake()["bsub_n"],4)
        self.assertIn("done(100)",self.fake()["jobs"]["102"]["dependency"])
        self.assertTrue(all(r["status"]=="submitted" for r in state["jobs"].values()))

    def test_resubmission_reuses_all_recorded_ids_without_duplicate_bsub(self):
        plan=self.plan([self.job("a"),self.job("b",["a"])])
        first=self.run_plan(plan);second=self.run_plan(plan)
        self.assertEqual(first["jobs"]["a"]["job_id"],second["jobs"]["a"]["job_id"])
        self.assertEqual(self.fake()["bsub_n"],2)

    def test_rejected_submission_preserves_prior_ids_and_resumes(self):
        self.change_fake(reject_n=2)
        plan=self.plan([self.job("a"),self.job("b",["a"]),self.job("c",["b"])])
        with self.assertRaisesRegex(RuntimeError,"submit_failed"):self.run_plan(plan)
        self.assertEqual(self.state()["jobs"]["a"]["job_id"],100)
        self.assertNotIn("c",self.state()["jobs"])
        state=self.run_plan(plan)
        self.assertEqual(state["jobs"]["b"]["job_id"],101)
        self.assertEqual(state["jobs"]["c"]["dependency"],"done(101) && ended(100)")
        self.assertEqual(self.fake()["bsub_n"],4)

    def test_success_without_parseable_id_never_retries_automatically(self):
        self.change_fake(opaque_n=1)
        plan=self.plan([self.job("a")])
        with self.assertRaisesRegex(RuntimeError,"submission_unknown"):self.run_plan(plan)
        with self.assertRaisesRegex(RuntimeError,"outcome unknown"):self.run_plan(plan)
        self.assertEqual(self.fake()["bsub_n"],1)

    def test_timeout_keeps_possible_accepted_id_and_does_not_retry(self):
        self.change_fake(slow_n=1)
        plan=self.plan([self.job("a")])
        with self.assertRaisesRegex(RuntimeError,"submission_unknown"):self.run_plan(plan,submit_timeout=.1)
        self.assertEqual(self.state()["jobs"]["a"]["job_id"],100)
        with self.assertRaisesRegex(RuntimeError,"outcome unknown"):self.run_plan(plan)
        self.assertEqual(self.fake()["bsub_n"],1)

    def test_plan_and_declared_input_change_cannot_reuse_ids(self):
        inp=self.path/"source.tsv";inp.write_text("old\n")
        plan=self.plan([self.job("a")],inputs=[str(inp)])
        self.run_plan(plan)
        inp.write_text("new input\n")
        with self.assertRaisesRegex(ValueError,"Plan/input identity changed"):self.run_plan(plan)
        frozen=json.loads((self.out/"plan.json").read_text())
        with self.assertRaisesRegex(ValueError,"Declared source input changed"):
            br.dispatch_frozen(frozen,self.out)
        self.assertEqual(self.fake()["bsub_n"],1)

    def test_cycles_unknown_dependencies_and_unsafe_headers_fail(self):
        for jobs in ([self.job("a",["b"]),self.job("b",["a"])],[self.job("a",["missing"])],
                     [self.job("a"),self.job("a")]):
            with self.subTest(jobs=jobs),self.assertRaises(ValueError):br.normalize_plan(self.plan(jobs),"hpc")
        job=self.job("a");job["resources"]["queue"]="ser\n#BSUB -q bad"
        with self.assertRaises(ValueError):br.normalize_plan(self.plan([job]),"hpc")
        job=self.job("a");job["argv"]="echo unsafe"
        with self.assertRaises(ValueError):br.normalize_plan(self.plan([job]),"hpc")

    def test_real_local_dag_quotes_arguments_and_resumes_completed_outputs(self):
        first=self.path/"a output.tsv";second=self.path/"b output.tsv"
        literal="literal $(touch SHOULD_NOT_EXIST) `not_a_command` $HOME"
        writer=[sys.executable,"-c","import pathlib,sys; pathlib.Path(sys.argv[1]).write_text(sys.argv[2])",str(first),literal]
        reader=[sys.executable,"-c","import pathlib,sys; pathlib.Path(sys.argv[2]).write_text(pathlib.Path(sys.argv[1]).read_text()+' done')",str(first),str(second)]
        plan=self.plan([self.job("first",payload=writer,outputs=[first]),self.job("second",["first"],reader,[second])])
        self.run_plan(plan,mode="local",jobs_n=2)
        self.assertEqual(second.read_text(),literal+" done");self.assertFalse((self.path/"SHOULD_NOT_EXIST").exists())
        before=second.stat().st_mtime_ns
        self.run_plan(plan,mode="local",jobs_n=2)
        self.assertEqual(second.stat().st_mtime_ns,before)
        second.write_text("tampered")
        with self.assertRaisesRegex(ValueError,"Completed worker output changed"):self.run_plan(plan,mode="local")

    def scientific_fixture(self):
        source=self.path/"scientific_source.tsv";source.write_text("input\n")
        product=self.path/"product.tsv";manifest=self.path/"manifest.json"
        code="""import json,pathlib,sys
source,product,manifest=map(pathlib.Path,sys.argv[1:])
product.write_text('measured result\\n')
def metadata(p):
    s=p.stat();return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
out=metadata(product);out.pop('path')
manifest.write_text(json.dumps(dict(status='complete',inputs={str(source):metadata(source)},output_inventory={product.name:out})))
"""
        writer=self.path/"scientific_fixture.py";writer.write_text(code)
        job=self.job("science",payload=[sys.executable,str(writer),str(source),str(product),str(manifest)],outputs=[manifest])
        plan=self.plan([job]);self.run_plan(plan,mode="local")
        return plan,source,product,manifest

    def test_manifest_inventory_prevents_reuse_after_a_product_is_deleted(self):
        plan,source,product,manifest=self.scientific_fixture()
        before=manifest.stat().st_mtime_ns;product.unlink()
        with self.assertRaisesRegex(ValueError,"Output identity changed"):
            self.run_plan(plan,mode="local")
        self.assertEqual(manifest.stat().st_mtime_ns,before)
        self.assertFalse(product.exists())

    def test_manifest_sources_are_stat_checked_even_when_manifest_is_unchanged(self):
        plan,source,product,manifest=self.scientific_fixture()
        source.write_text("changed scientific input\n")
        with self.assertRaisesRegex(ValueError,"Scientific file identity changed"):
            self.run_plan(plan,mode="local")

    def test_inventory_verification_does_not_read_large_products_or_input_contents(self):
        product=self.path/"large.bin";source=self.path/"large.input"
        for path in (product,source):
            with path.open("wb") as handle: handle.truncate(64*1024*1024)
        expected=br.path_identity(product);expected.pop("path")
        manifest=self.path/"primary_manifest.json"
        data=dict(status="complete",inputs={str(source):br.path_identity(source)},output_inventory={product.name:expected})
        manifest.write_text(json.dumps(data))
        original=Path.open
        def limited_open(path,*args,**kwargs):
            if path in (source,product): raise AssertionError("Large scientific content was read")
            return original(path,*args,**kwargs)
        with patch.object(Path,"open",limited_open):br.validate_output_file(manifest)
        data["output_inventory"][product.name]["sha256"]="0"*64
        manifest.write_text(json.dumps(data))
        with patch.object(Path,"open",limited_open),self.assertRaisesRegex(ValueError,"checksum verification in compute"):
            br.validate_output_file(manifest)

    def test_ordinary_json_marker_is_not_misclassified_as_scientific_manifest(self):
        marker=self.path/"barrier.json";marker.write_text('{"status":"done","count":12}')
        br.validate_output_file(marker)

    def test_malformed_inventory_does_not_verify_an_absent_file(self):
        marker=self.path/"manifest.json"
        marker.write_text(json.dumps(dict(status="complete",output_inventory={"absent.tsv":{}})))
        with self.assertRaisesRegex(ValueError,"lacks verifiable size/mtime"):
            br.validate_output_file(marker)

    def test_empty_outputs_record_execution_but_cannot_claim_verified_cache_reuse(self):
        sentinel=self.path/"undeclared-product"
        job=self.job("unverifiable",payload=[sys.executable,"-c","import pathlib,sys;pathlib.Path(sys.argv[1]).write_text('executed')",str(sentinel)])
        plan=self.plan([job]);self.run_plan(plan,mode="local")
        row=json.loads((self.out/"runs/unverifiable.json").read_text())
        self.assertEqual(row["output_verification"],"unavailable_no_declared_outputs")
        before=sentinel.stat().st_mtime_ns
        with self.assertRaisesRegex(ValueError,"no declared outputs to verify"):
            self.run_plan(plan,mode="local")
        self.assertEqual(sentinel.stat().st_mtime_ns,before)
        self.assertNotIn("reused_verified_outputs",row)

    def test_local_failed_parent_blocks_child(self):
        target=self.path/"should-not-exist"
        plan=self.plan([self.job("bad",payload=[sys.executable,"-c","raise SystemExit(7)"]),
                        self.job("child",["bad"],[sys.executable,"-c",f"open({str(target)!r},'w').write('bad')"])])
        with self.assertRaisesRegex(RuntimeError,"local jobs failed"):self.run_plan(plan,mode="local",jobs_n=2)
        self.assertEqual(self.state()["jobs"]["bad"]["returncode"],7)
        self.assertEqual(self.state()["jobs"]["child"]["status"],"blocked_dependency")
        self.assertFalse(target.exists())

    def test_zero_exit_without_required_output_is_failure(self):
        plan=self.plan([self.job("empty",outputs=[self.path/"absent.json"])])
        with self.assertRaises(RuntimeError):self.run_plan(plan,mode="local")
        row=json.loads((self.out/"runs/empty.json").read_text())
        self.assertEqual(row["status"],"failed");self.assertIn("declared completed output",row["error"])

    def test_local_on_lsf_submission_host_is_refused(self):
        plan=self.plan([self.job("a")])
        with patch.dict(os.environ,{"PATH":str(self.bin)+os.pathsep+os.environ["PATH"]}),self.assertRaisesRegex(RuntimeError,"LSF submission host"):
            self.run_plan(plan,mode="local")
        self.assertNotIn("bsub_n",self.fake())

    def test_site_memory_requires_explicit_confirmation(self):
        plan=self.plan([self.job("a")],lsf_memory_scope="site_host")
        with self.assertRaises(ValueError):br.normalize_plan(plan,"hpc")
        plan["lsf_resource_reserve_per_task"]=False
        n=br.normalize_plan(plan,"hpc");self.assertEqual(n["jobs"][0]["resources"]["lsf_memory"],"1GB")
        plan["lsf_memory_scope"]="site_task";plan["lsf_resource_reserve_per_task"]=True
        plan["jobs"][0]["resources"].update(memory_gib=10,cores=4)
        n=br.normalize_plan(plan,"hpc");self.assertEqual(n["jobs"][0]["resources"]["lsf_memory"],"3GB")

    def test_generated_hpc_cmd_executes_only_in_allocation_with_receipt(self):
        product=self.path/"actual.tsv"
        job=self.job("a",payload=[sys.executable,"-c","import os,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(os.environ['ALLELE_ALLOCATED_CORES'])",str(product)],outputs=[product])
        self.run_plan(self.plan([job]),write_only=True)
        env={**os.environ,"LSB_JOBID":"100","LSB_DJOB_NUMPROC":"1"}
        result=subprocess.run(["bash",str(self.out/"cmd/a.cmd")],capture_output=True,text=True,env=env)
        self.assertEqual(result.returncode,0,result.stderr);self.assertEqual(product.read_text(),"1")
        row=json.loads((self.out/"runs/a.json").read_text())
        self.assertEqual(row["status"],"succeeded");self.assertEqual(row["lsf_job_id"],"100")
        self.assertEqual(row["outputs"][0]["path"],str(product))

    def test_explicit_retry_rewires_pending_dependency_and_retains_ids(self):
        plan=self.plan([self.job("a"),self.job("child",["a"])])
        self.run_plan(plan)
        state=self.fake();state["jobs"]["100"]["state"]="EXIT";self.scheduler.write_text(json.dumps(state))
        result=self.run_plan(plan,retry_failed=True)
        self.assertEqual(result["jobs"]["a"]["job_id"],102)
        self.assertEqual(result["jobs"]["child"]["job_id"],101)
        self.assertEqual(self.fake()["jobs"]["101"]["dependency"],"done(102)")
        mods=[c for c in self.fake()["calls"] if c[0]=="bmod"]
        self.assertEqual(mods,[["bmod","-w","done(100)","101"],["bmod","-w","done(102)","101"]])
        self.assertEqual(result["retry"]["status"],"complete")

    def test_retry_bmod_failure_resumes_without_resubmitting_new_id(self):
        plan=self.plan([self.job("a"),self.job("child",["a"])])
        self.run_plan(plan)
        state=self.fake();state["jobs"]["100"]["state"]="EXIT";state["reject_bmod_n"]=2
        self.scheduler.write_text(json.dumps(state))
        with self.assertRaisesRegex(RuntimeError,"bmod failed"):self.run_plan(plan,retry_failed=True)
        self.assertEqual(self.state()["jobs"]["a"]["job_id"],102)
        with self.assertRaisesRegex(RuntimeError,"explicit retry is incomplete"):self.run_plan(plan)
        self.run_plan(plan,retry_failed=True)
        self.assertEqual(self.fake()["bsub_n"],3)
        self.assertEqual(self.fake()["jobs"]["101"]["dependency"],"done(102)")

    def test_retry_does_not_touch_running_or_unknown_scheduler_states(self):
        plan=self.plan([self.job("a"),self.job("b")]);self.run_plan(plan)
        state=self.fake();state["jobs"]["100"]["state"]="EXIT";state["jobs"]["101"]["state"]="RUN"
        self.scheduler.write_text(json.dumps(state))
        with self.assertRaisesRegex(RuntimeError,"running/unknown jobs are untouched"):self.run_plan(plan,retry_failed=True)
        self.assertEqual(self.fake()["bsub_n"],2);self.assertNotIn("bmod_n",self.fake())

    def test_wait_does_not_infer_success_from_exit_or_empty_scheduler(self):
        plan=self.plan([self.job("a")]);self.run_plan(plan)
        state=self.fake();state["jobs"]["100"]["state"]="EXIT";self.scheduler.write_text(json.dumps(state))
        with self.assertRaisesRegex(RuntimeError,"failed; done dependencies remain"):self.run_plan(plan,wait=True)
        state=self.fake();state["jobs"]={};self.scheduler.write_text(json.dumps(state))
        with self.assertRaisesRegex(RuntimeError,"state unavailable"):self.run_plan(plan,wait=True)

    def test_wait_ignores_receipts_from_a_replaced_lsf_attempt(self):
        plan=self.plan([self.job("a")]);self.run_plan(plan)
        frozen=json.loads((self.out/"plan.json").read_text())
        (self.out/"runs").mkdir(exist_ok=True)
        for status in ("succeeded","failed"):
            (self.out/"runs/a.json").write_text(json.dumps(dict(
                plan_identity=frozen["identity"],status=status,lsf_job_id="99",outputs=[])))
            with self.assertRaisesRegex(RuntimeError,"Wait timeout"):
                self.run_plan(plan,wait=True,wait_timeout=0)
        row=json.loads((self.out/"runs/a.json").read_text())
        row.update(status="succeeded",lsf_job_id="100")
        (self.out/"runs/a.json").write_text(json.dumps(row))
        self.run_plan(plan,wait=True,wait_timeout=0)

    def test_public_cli_frozen_dispatch_does_not_require_payload_files_for_planning(self):
        plan=self.plan([self.job("a",payload=["/not-mounted-yet/python","/not-mounted-yet/work.py"])])
        path=self.path/"input_plan.json";path.write_text(json.dumps(plan))
        command=[sys.executable,str(ROOT/"f/batch_runner.py"),"--plan",str(path),"--out",str(self.out),"--mode","hpc","--write-cmd"]
        result=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        result=subprocess.run(["bash",str(self.out/"submit.sh"),"--bsub",str(self.bin/"bsub")],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr);self.assertEqual(self.fake()["bsub_n"],1)


if __name__=="__main__":unittest.main()

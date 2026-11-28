"""Project-neutral plan/execution lifecycle; definitions only describe the work."""
from __future__ import annotations
from dataclasses import dataclass, field
import json
import os
from pathlib import Path
from workflow_jobs import write_plan, run_plan


@dataclass
class Workflow:
    jobs: list
    root: object
    name: str
    mode: str
    plan_only: bool = False
    parallel: int = 1
    resume: bool = True
    wait: object = None
    summary: dict = field(default_factory=dict)


@dataclass
class ScriptWorkflow:
    path: Path
    lines: list
    resources: dict
    name: str
    logs: Path
    mode: str
    metadata: dict = field(default_factory=dict)
    plan_only: bool = False


def execute_workflow(spec):
    if isinstance(spec, ScriptWorkflow):
        from command_jobs import lsf_arguments, execute_script
        spec.path.write_text('\n'.join(spec.lines) + '\n'); spec.path.chmod(0o755)
        spec.logs.mkdir(parents=True, exist_ok=True)
        argv = lsf_arguments(spec.name, spec.resources, spec.logs/'run.out', spec.logs/'run.err')
        result = execute_script(spec.path, spec.mode, argv=argv, receipt=spec.path.parent/'submission.json',
                                metadata=spec.metadata, plan_only=spec.plan_only)
    else:
        if not isinstance(spec, Workflow): raise TypeError('Expected a workflow definition')
        path = write_plan(spec.jobs, spec.root, mode=spec.mode, name=spec.name)
        if spec.summary: print(json.dumps(spec.summary, ensure_ascii=False), flush=True)
        wait = os.getenv('SUBMIT_WAIT', 'FALSE').upper() == 'TRUE' if spec.wait is None else spec.wait
        result = run_plan(path, submit=spec.mode, jobs=spec.parallel, wait=wait,
                          dry_run=spec.plan_only, resume=spec.resume)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    return result

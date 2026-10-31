#!/usr/bin/env python3
"""Check only the selected stage's actual interpreter and executables."""
import argparse
import importlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['render','trajectory','project','baseline','time','scan','herv','all'],default='trajectory')
    p.add_argument('--plink2',default='plink2')
    a=p.parse_args()
    packages={'render':['numpy','matplotlib','plotly'],'trajectory':['numpy','scipy','matplotlib'],
              'project':['numpy','pgenlib'],'baseline':['numpy','matplotlib'],
              'time':['numpy','scipy','matplotlib','pgenlib','pysam','pyliftover','zstandard'],
              'scan':['numpy','scipy','matplotlib','pysam'],'herv':['numpy','matplotlib','pysam']}
    binaries={'project':[a.plink2],'scan':['bcftools'],'herv':['bcftools','samtools']}
    stages=list(packages) if a.stage=='all' else [a.stage]
    result={'python':str(Path(sys.executable).resolve()),'python_version':sys.version,'stage':a.stage,'packages':{},'executables':{},'errors':[]}
    if 'herv' in stages:
        spec=importlib.util.spec_from_file_location('herv_dependency_check',Path(__file__).with_name('1.herv.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        executable=module.repeatmasker_path(SimpleNamespace(repeatmasker=None,out=Path('/nonexistent-allele3-check-only')))
        binaries['herv'] += [str(executable),str(executable.parent/'rmblastn')]
        runtime=module.pgen_runtime()
        result['herv_pgen_python']=runtime
        if not runtime:result['errors'].append('No HERV PGEN interpreter with phased allele support')
    for name in sorted({x for s in stages for x in packages[s]}):
        try:
            module=importlib.import_module(name)
            result['packages'][name]={'version':str(getattr(module,'__version__','not_exposed')),'path':getattr(module,'__file__',None)}
        except ImportError as exc:result['errors'].append(str(exc))
    for name in sorted({x for s in stages for x in binaries.get(s,[])}):
        found=shutil.which(name)
        if not found:result['errors'].append('Executable not found: '+name);continue
        path=str(Path(found).absolute())
        try:
            option={'RepeatMasker':'-v','rmblastn':'-version'}.get(Path(name).name,'--version')
            env=dict(os.environ,PATH=str(Path(path).parent)+os.pathsep+os.environ.get('PATH',''))
            version=subprocess.run([path,option],capture_output=True,text=True,errors='replace',timeout=20,env=env)
            result['executables'][name]={'path':path,'version':(version.stdout+version.stderr).splitlines()[:3],'exit_code':version.returncode}
            if version.returncode:result['errors'].append('Version check failed: '+path)
        except (OSError,subprocess.SubprocessError) as exc:result['errors'].append(str(exc))
    result['status']='failed' if result['errors'] else 'verified_selected_dependencies'
    print(json.dumps(result,indent=2));return bool(result['errors'])

if __name__=='__main__':raise SystemExit(main())

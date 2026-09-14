"""Versioned provenance: help text is not an analysis change; real changes remain guarded."""
import ast
import hashlib
import json
from pathlib import Path


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


class AnalysisSyntax(ast.NodeTransformer):
    def visit_ClassDef(self,node):
        if node.name=='HelpFormatter':return None
        return self.generic_visit(node)

    def visit_Call(self,node):
        node=self.generic_visit(node)
        if isinstance(node.func,ast.Attribute) and node.func.attr in ('ArgumentParser','add_argument'):
            presentation={'description','epilog','usage','prog','formatter_class','help'}
            node.keywords=[k for k in node.keywords if k.arg not in presentation]
        return node


def semantic_hash(path):
    if path.suffix=='.py':
        tree=AnalysisSyntax().visit(ast.parse(path.read_text()))
        content=ast.dump(tree,include_attributes=False).encode()
    else:content=path.read_bytes()
    return hashlib.sha256(content).hexdigest()


def code_fingerprints(home):
    # Publication rendering has its own source/data/output signature and cannot invalidate model caches.
    paths=[home/'panome.sh']+sorted(p for p in (home/'f').glob('*.py') if p.name not in {'publication.py','final.py'})+sorted((home/'f').glob('*.R'))
    raw={str(p.relative_to(home)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    semantic={str(p.relative_to(home)):semantic_hash(p) for p in paths}
    return raw,semantic


def build_provenance(home,config,inputs,versions,python):
    raw,semantic=code_fingerprints(home)
    payload=dict(cache_schema=2,config=config,inputs=inputs,analysis_code=semantic,versions=versions,python=python)
    return dict(signature=digest(payload),code=raw,**payload)


def mismatch(old,new):
    reasons=[]
    for key in ('inputs','config','versions','python'):
        if old.get(key)!=new.get(key):
            if key=='config':
                changed=sorted(k for k in set(old.get(key,{}))|set(new[key]) if old.get(key,{}).get(k)!=new[key].get(k))
                reasons.append('parameters: '+', '.join(changed))
            else:reasons.append(key)
    if old.get('analysis_code')!=new.get('analysis_code'):reasons.append('analysis code')
    return '; '.join(reasons)


def publication_upgrade_code(before,after,home):
    path=home/'provenance/publication_cache_upgrade.json'
    if not path.exists():return False
    bridge=json.loads(path.read_text())
    approved_previous=[bridge['previous_analysis_code']]+bridge.get('compatible_previous_analysis_codes',[])
    return before in approved_previous and after==bridge['approved_analysis_code']


def known_publication_upgrade(old,new,home):
    if old.get('cache_schema')!=2:return False
    if any(old.get(k)!=new[k] for k in ('config','inputs','versions','python')):return False
    payload={k:old[k] for k in ('cache_schema','config','inputs','analysis_code','versions','python')}
    return old.get('signature')==digest(payload) and publication_upgrade_code(old['analysis_code'],new['analysis_code'],home)


def known_legacy_upgrade(old,new,home):
    # Narrow release bridge: never bless arbitrary code changes in an old manifest.
    bridge_path=home/'provenance/cache_v1_to_v2.json'
    if not bridge_path.exists():return False
    bridge=json.loads(bridge_path.read_text())
    snapshot=home/'provenance/legacy_panome.py.txt'
    if hashlib.sha256(snapshot.read_bytes()).hexdigest()!=bridge['legacy_code']['f/panome.py']:return False
    if old.get('code')!=bridge['legacy_code']:return False
    if new['analysis_code']!=bridge['approved_v2_analysis_code'] and not publication_upgrade_code(bridge['approved_v2_analysis_code'],new['analysis_code'],home):return False
    if any(old.get(k)!=new[k] for k in ('inputs','config','versions','python')):return False
    payload={k:old[k] for k in ('config','inputs','code','versions','python')}
    return old.get('signature')==digest(payload)


def ensure_manifest(path,new,home,replace=False):
    if path.exists() and not replace:
        old=json.loads(path.read_text())
        if old.get('signature')==new['signature']:
            return old  # Preserve original execution provenance, even if current help text differs.
        if (old.get('cache_schema')!=2 and known_legacy_upgrade(old,new,home)) or known_publication_upgrade(old,new,home):
            history=path.parent/'provenance_history';history.mkdir(exist_ok=True)
            backup=history/(old['signature']+'.json')
            if not backup.exists():backup.write_bytes(path.read_bytes())
            new={**new,'compatible_signatures':list(dict.fromkeys([old['signature']]+old.get('compatible_signatures',[]))),
                 'migration':{'reason':'Verified presentation/export and cache bookkeeping upgrade; no analysis refit.',
                              'original_manifest':str(backup.relative_to(path.parent))}}
            print('Cache metadata upgraded: verified existing analysis; original provenance retained.',flush=True)
        else:
            raise ValueError('Cache mismatch ('+mismatch(old,new)+'). Use a new --run-name or --replace only when a new analysis is intended.')
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(new,indent=2,default=str,allow_nan=False));temporary.replace(path)
    return new


def stage_cached(done,manifest):
    if not done.exists():return False
    marker=json.loads(done.read_text())
    accepted=[manifest['signature']]+manifest.get('compatible_signatures',[])
    if marker.get('signature') not in accepted:raise ValueError(f'Stage provenance does not match manifest: {done}')
    return True

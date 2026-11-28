"""Publish durable scientific outputs independently of disposable runtime files."""
from pathlib import Path
import html
import json
import os
import re
import shutil
import tempfile


def publish(layout, status):
    from analysis_integrity import atomic_json
    from workflow_data_check import _layout_destination
    from output_layout import prepare
    from product_manifest import write_inventory
    receipt_path = layout.qc/'layout.json'
    previous = json.loads(receipt_path.read_text()) if receipt_path.is_file() else {}
    if status != 'complete':
        receipt = dict(previous, status=status, public_root=str(layout.output), runtime_root=str(layout.runtime))
        atomic_json(receipt_path, receipt)
        return receipt  # An unsuccessful attempt cannot replace the last good dataset.
    layout.output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.'+layout.output.name+'.publish-', dir=layout.output.parent))
    backup = staging.with_name(staging.name+'.previous')
    public, mappings, support = {}, {}, {}
    excluded = {'shards','variants','chunks','candidate_cache','__pycache__','history','previous','tmp','logs','jobs'}
    try:
        for parent, dirs, names in os.walk(layout.runtime, followlinks=False):
            dirs[:] = sorted(d for d in dirs if d not in excluded and not (Path(parent)/d).is_symlink() and not d.endswith('.partial'))
            for name in sorted(names):
                source = Path(parent)/name
                if source.is_symlink() or not source.is_file() or name.startswith('.') or name.endswith(('.lock','.partial','.tmp','.log','.out','.err','.sh','.cmd')):
                    continue
                relative = source.relative_to(layout.runtime)
                policy = _layout_destination(layout, relative)
                if policy and policy[0] == 'public':
                    target_relative = policy[1]
                elif policy and policy[0] == 'qc':
                    target_relative = Path('metadata')/policy[1]
                elif name.endswith(('.tsv','.tsv.gz','.csv','.bed','.bed.gz','.json','.png','.pdf','.svg','.html','.txt','.txt.gz','.npz')):
                    target_relative = Path('artifacts')/relative
                else:
                    continue
                target = staging/target_relative
                if target.exists(): raise ValueError('Colliding product: '+str(target_relative))
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                mappings[str(relative)] = str(target_relative)
                public[str(target_relative)] = dict(native=str(relative),size=source.stat().st_size)
        if not public: raise ValueError('No scientific products to publish: '+str(layout.runtime))
        # Real native input snapshots support directory readers after tmp is removed.
        for relative, target_relative in mappings.items():
            alias = staging/'artifacts'/relative
            target = staging/target_relative
            if alias == target: continue
            alias.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(target, alias)
        for name, spec in public.items():
            target = staging/name
            if target.suffix != '.html': continue
            native_parent = (layout.runtime/spec['native']).parent
            def link(match):
                attribute, quote, url = match.groups()
                if re.match(r'^(?:[a-zA-Z]+:|#|//)',url):return match.group(0)
                path, sep, fragment = html.unescape(url).partition('#')
                absolute = Path(os.path.normpath(str(native_parent/path)))
                if not absolute.is_relative_to(layout.runtime):return match.group(0)
                mapped = mappings.get(str(absolute.relative_to(layout.runtime)))
                if mapped is None:return match.group(0)
                url = os.path.relpath(staging/mapped,target.parent)+(('#'+fragment) if sep else '')
                return attribute+'='+quote+html.escape(url,quote=True)+quote
            target.write_text(re.sub(r'''(href|src)=(["'])(.*?)\2''',link,target.read_text()))
            spec['size']=target.stat().st_size
        receipt = dict(schema='allele3-module-layout-1',status='complete',public_root=str(layout.output),
            runtime_root=str(layout.runtime),public_files=public,support_files=support,
            native_to_public=mappings,directory_map={'.':'artifacts'},
            provenance='Scientific products and metadata are permanent; artifacts is a portable native input view.')
        atomic_json(staging/'metadata/layout.json',receipt)
        if layout.runtime.is_relative_to(layout.output):
            # Runtime and products share the module root. Never rotate that root:
            # it contains the active worker, lock, logs and scheduler commands.
            backup.mkdir()
            installed = []
            try:
                for child in list(staging.iterdir()):
                    target = layout.output/child.name
                    if target.exists() or target.is_symlink():target.rename(backup/child.name)
                    child.rename(target)
                    installed.append(child.name)
                prepare(layout.output)
                write_inventory(layout.output)
            except BaseException:
                for name in installed:
                    target=layout.output/name
                    if target.is_dir() and not target.is_symlink():shutil.rmtree(target)
                    elif target.exists() or target.is_symlink():target.unlink()
                for child in backup.iterdir():child.rename(layout.output/child.name)
                raise
            atomic_json(receipt_path,receipt)
            shutil.rmtree(backup)
            return receipt
        if layout.output.exists():layout.output.rename(backup)
        try:
            staging.rename(layout.output)
            prepare(layout.output)
            write_inventory(layout.output)
        except BaseException:
            if layout.output.exists():shutil.rmtree(layout.output)
            if backup.exists():backup.rename(layout.output)
            raise
        atomic_json(receipt_path,receipt)
        if backup.exists():shutil.rmtree(backup)
        return receipt
    finally:
        if staging.exists():shutil.rmtree(staging)

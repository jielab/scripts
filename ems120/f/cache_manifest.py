"""Atomic completion receipts for declared scientific inputs and output files.

No source-code fingerprints or historical cache adoption are used. Selected
training columns have a semantic digest; other inputs use paths, sizes and
modification times. Output completeness uses the recorded file set and sizes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


def state(value):
    item = value if isinstance(value, dict) else {'path': value}
    p = Path(item['path']).resolve()
    if not p.exists():
        return dict(path=str(p), missing=True)
    if item.get('columns'):
        import re
        import zipfile
        import xml.etree.ElementTree as ET
        ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(p) as book:
            shared = []
            if 'xl/sharedStrings.xml' in book.namelist():
                shared = [''.join(n.itertext()) for n in ET.fromstring(book.read('xl/sharedStrings.xml'))]
            rows = ET.fromstring(book.read('xl/worksheets/sheet1.xml')).findall('.//s:row', ns)
            values = []
            for row in rows:
                cells = {}
                for cell in row.findall('s:c', ns):
                    col = re.sub(r'\d+', '', cell.attrib['r'])
                    node = cell.find('s:v', ns)
                    value = '' if node is None else node.text or ''
                    if cell.get('t') == 's':
                        value = shared[int(value)]
                    elif cell.get('t') == 'inlineStr':
                        inline = cell.find('s:is', ns)
                        value = '' if inline is None else ''.join(inline.itertext())
                    cells[col] = value
                values.append(cells)
            header = values[0] if values else {}
            chosen = {v: k for k, v in header.items() if v in item['columns']}
            if set(chosen) != set(item['columns']):
                raise ValueError('Training columns missing')
            data = [[row.get(chosen[c], '') for c in item['columns']] for row in values[1:]]
        return dict(path=str(p), columns=item['columns'],
                    training_sha256=hashlib.sha256(json.dumps(data, ensure_ascii=False).encode()).hexdigest())
    if p.is_dir():
        return dict(path=str(p), files=[state(str(f)) for f in sorted(p.rglob('*'))
                                       if f.is_file() and not f.name.endswith('.lock')])
    info = p.stat()
    return dict(path=str(p), size=info.st_size, mtime_ns=info.st_mtime_ns)


def output_state(path):
    p = Path(path).resolve()
    if not p.is_file() or p.stat().st_size <= 0:
        raise ValueError(f'Missing/empty output: {p}')
    if p.suffix == '.json':
        json.loads(p.read_text(encoding='utf-8'))
    return dict(path=str(p), size=p.stat().st_size)


def identity(spec):
    return dict(schema='ems120-completion-v1', stage=spec['stage'],
                parameters=spec.get('parameters', {}), inputs=[state(p) for p in spec.get('inputs', [])])


def check(receipt, spec=None, seen=None, required_outputs=None, validate_inputs=True):
    try:
        path = Path(receipt).resolve()
        seen = set() if seen is None else set(seen)
        if path in seen:
            return False
        seen.add(path)
        saved = json.loads(path.read_text(encoding='utf-8'))
        expected = saved['spec'] if spec is None else spec
        if saved.get('status') != 'complete':
            return False
        if validate_inputs and saved['identity'] != identity(expected):
            return False
        outputs = saved['outputs']
        if required_outputs is not None:
            required = {str(Path(p).resolve()) for p in required_outputs}
            outputs = [x for x in outputs if x['path'] in required]
            if {x['path'] for x in outputs} != required:
                return False
        if not outputs or any(output_state(x['path']) != x for x in outputs):
            return False
        return not validate_inputs or all(check(p, seen=seen) for p in expected.get('dependencies', []))
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return False


def seal(receipt, spec, outputs):
    if not all(check(p) for p in spec.get('dependencies', [])):
        raise ValueError('Upstream cache is incomplete or changed; cannot publish dependent completion')
    current = [output_state(p) for p in outputs]
    if not current:
        raise ValueError('No cache outputs')
    value = dict(status='complete', spec=spec, identity=identity(spec), outputs=current)
    path = Path(receipt)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.receipt-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, sort_keys=True, ensure_ascii=False)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument('action', choices=['check', 'seal'])
    p.add_argument('--receipt', required=True)
    p.add_argument('--spec')
    p.add_argument('--outputs-only', action='store_true', help='validate completed files without requiring training inputs')
    p.add_argument('--output', action='append', default=[])
    a = p.parse_args()
    spec = json.loads(Path(a.spec).read_text()) if a.spec else None
    if a.action == 'check':
        raise SystemExit(0 if check(a.receipt, spec, validate_inputs=not a.outputs_only) else 1)
    if spec is None:
        p.error('seal requires --spec')
    seal(a.receipt, spec, a.output)

"""Bounded, resumable derived FASTA inputs; source references are never modified."""
import gzip
import json
import hashlib
from pathlib import Path
from contextlib import ExitStack

from common import atomic_json
from library import digest, file_state

# MetaCache defaults: k=16, window=127, stride=112. Align cuts to its window grid.
CORE = (1 << 20) // 112 * 112
OVERLAP = 15
SHARD_BASES = 512 << 20
VERSION = 'bounded-fasta-v2'


def pieces(path, core=CORE):
    """Yield (unique synthetic ID, bases) with bounded reads, even for unwrapped DNA."""
    if core < 112 or core % 112:
        raise ValueError('chunk core must be a multiple of the 112-base window stride')
    with ExitStack() as stack:
        raw = stack.enter_context(open(path, 'rb', buffering=8 << 20))
        stream = stack.enter_context(gzip.GzipFile(fileobj=raw)) if str(path).endswith('.gz') else raw
        header = None
        record = 0
        offset = 0
        buf = bytearray()
        line_start = True
        source_id = hashlib.sha256(Path(path).name.encode()).hexdigest()[:24]

        def identifier():
            return f'metagchunk-{header.hex()}-{source_id}-{record}-{offset}'

        while True:
            line = stream.readline(1 << 20)
            if not line or (line_start and line.startswith(b'>')):
                if header is not None and buf and (offset == 0 or len(buf) > OVERLAP):
                    yield identifier(), bytes(buf)
                if not line:
                    break
                if not line.endswith(b'\n') or len(line) > 65536:
                    raise ValueError(f'Invalid/oversized FASTA header: {path}')
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError(f'Empty FASTA header: {path}')
                header = tokens[0]
                record += 1
                offset = 0
                buf = bytearray()
            else:
                if header is None or b'>' in line:
                    raise ValueError(f'Malformed FASTA: {path}')
                buf.extend(line.translate(None, b'\r\n \t'))
                while len(buf) >= core + OVERLAP:
                    yield identifier(), bytes(buf[:core + OVERLAP])
                    offset += core
                    del buf[:core]
            line_start = line.endswith(b'\n')
        if header is None:
            raise ValueError(f'Empty FASTA: {path}')


def original_accession(identifier):
    if not identifier.startswith('metagchunk-'):
        return identifier
    try:
        fields = identifier.split('-')
        if len(fields) == 5:
            _, encoded, source_id, record, offset = fields
            if len(source_id) != 24 or any(c not in '0123456789abcdef' for c in source_id):
                raise ValueError('invalid source id')
        else:
            _, encoded, record, offset = fields
        int(record); int(offset)
        return bytes.fromhex(encoded).decode('utf-8')
    except (ValueError, UnicodeError) as exc:
        raise RuntimeError('Invalid derived reference identifier') from exc


def valid_receipt(path, identity):
    try:
        saved = json.loads(path.read_text())
        return saved if (saved['identity'] == identity and saved['outputs'] and
                         all(file_state(p) == s for p, s in saved['outputs'].items())) else None
    except (OSError, ValueError, KeyError):
        return None


def shards(part, cache, identity, max_bases=SHARD_BASES):
    """Publish each shard before yielding; a whole-input receipt avoids rescanning on reuse.

    In an interrupted conversion the source is streamed again, but matching published
    shards are not rewritten and their completed databases remain reusable.
    """
    if max_bases < CORE + OVERLAP:
        raise ValueError('shard base limit is smaller than a sequence chunk')
    sources = sorted(p for p in part.iterdir() if p.name != 'assembly_summary.txt')
    source_states = {str(p.resolve()): file_state(p) for p in sources}
    key = digest([identity, VERSION, CORE, OVERLAP, max_bases, source_states])
    root = cache / key
    root.mkdir(parents=True, exist_ok=True)
    manifest = root / 'complete.json'
    try:
        saved = json.loads(manifest.read_text())
        if saved['key'] == key and all(valid_receipt(root/name/'receipt.json', key) for name in saved['shards']):
            for name in saved['shards']:
                yield root/name, key
            return
    except (OSError, ValueError, KeyError):
        pass
    names = []
    folder = None
    handle = None
    used = 0
    reuse = False
    outputs = {}
    current_name = None

    def close_file():
        nonlocal handle
        if handle is not None:
            handle.close()
            target = folder/current_name
            target.with_name(target.name+'.partial').replace(target)
            outputs[str(target)] = file_state(target)
            handle = None

    def finish():
        close_file()
        if not reuse:
            atomic_json(folder/'receipt.json', {'identity': key, 'bases': used, 'outputs': outputs})
        print(f'[MetaCache inputs] {part.name}/{folder.name}: {used:,} bases; '+
              ('reused' if reuse else 'published'), flush=True)

    try:
        for source in sources:
            for identifier, sequence in pieces(source):
                if folder is None or used + len(sequence) > max_bases:
                    if folder is not None:
                        finish()
                        yield folder, key
                    folder = root/f'shard{len(names)+1:05d}'
                    names.append(folder.name)
                    folder.mkdir(exist_ok=True)
                    reuse = bool(valid_receipt(folder/'receipt.json', key))
                    used = 0
                    outputs = {}
                    current_name = None
                    if not reuse:
                        (folder/'receipt.json').unlink(missing_ok=True)
                        for stale in folder.iterdir():
                            if stale.is_file():
                                stale.unlink()
                        mapping = folder/'assembly_summary.txt'
                        mapping.write_bytes((part/'assembly_summary.txt').read_bytes())
                        outputs[str(mapping)] = file_state(mapping)
                if not reuse:
                    if current_name != source.name:
                        close_file()
                        current_name = source.name
                        handle = (gzip.open(folder/(current_name+'.partial'), 'wb', compresslevel=1)
                                  if current_name.endswith('.gz') else
                                  open(folder/(current_name+'.partial'), 'wb'))
                    handle.write(b'>'+identifier.encode('ascii')+b'\n'+sequence+b'\n')
                used += len(sequence)
            close_file()
            current_name = None
            if file_state(source) != source_states[str(source.resolve())]:
                raise RuntimeError(f'Reference changed during conversion: {source}')
        if folder is None:
            raise RuntimeError(f'No FASTA bases in {part}')
        finish()
        yield folder, key
        if any(file_state(p) != state for p, state in source_states.items()):
            raise RuntimeError('Reference changed during conversion')
        atomic_json(manifest, {'key': key, 'shards': names})
    finally:
        if handle is not None:
            handle.close()

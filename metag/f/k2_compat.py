#!/usr/bin/env python3
"""Run installed k2 2.17.1 with checked download and library-copy fixes.

Keep the installed package untouched. Refuse unfamiliar source rather than
silently applying a partial patch. Original assembly summaries remain intact.
"""
import shutil
import re
import textwrap
import time
from pathlib import Path


def write_taxid_blocks(in_file, out_file, filepath, manifest_to_taxid,
                       accession_to_taxid, filepath_to_url, *, binary=False,
                       block_size=4 * 1024 * 1024):
    """Rewrite FASTA headers in bounded batches, preserving k2's text semantics.

    Input is still a text stream (including universal-newline handling). Only
    headers require Python callbacks; sequence lines pass through in bulk. The
    same bytes are fed to the unchanged masker, with the same sequence/URL map.
    """
    sequences_added = ch_added = 0
    sequence_to_url = {}

    def header(line):
        nonlocal sequences_added, ch_added
        ch_added -= len(line) - int(line.endswith('\n'))
        taxid = manifest_to_taxid[filepath]
        if not taxid:
            accession = re.search(r'GC[AF]_[0-9]{9}\.\d+', line)
            if not accession or accession.group() not in accession_to_taxid:
                raise RuntimeError(f'Unable to assign taxid to sequence {line!r} in {filepath}')
            taxid = accession_to_taxid[accession.group()]
        line = line.replace('>', '>kraken:taxid|' + taxid + '|', 1)
        sequence_to_url[line] = filepath_to_url[filepath]
        sequences_added += 1
        return line

    while block := in_file.read(block_size):
        # Complete a line before the next block, including headers at boundaries.
        if not block.endswith('\n'):
            block += in_file.readline()
        ch_added += len(block) - block.count('\n')
        # Keep upstream's len(line)-1 accounting even without a final newline.
        if not block.endswith('\n') and not block.rsplit('\n', 1)[-1].startswith('>'):
            ch_added -= 1
        # str.split scans sequence text in C; process only header-containing
        # pieces in Python. Retain '>' characters inside header descriptions.
        parts = block.split('\n>')
        if parts[0].startswith('>'):
            line, newline, sequence = parts[0].partition('\n')
            parts[0] = header(line + newline) + sequence
        for i in range(1, len(parts)):
            line, newline, sequence = parts[i].partition('\n')
            parts[i] = header('>' + line + newline)[1:] + sequence
        block = '\n>'.join(parts)
        out_file.write(block.encode() if binary else block)
    return sequences_added, ch_added, sequence_to_url


def merge_library(filepaths, library_filename, log):
    """Concatenate the processed FASTAs in manifest order with bounded memory."""
    paths = [Path(filepath).with_suffix('') if str(filepath).endswith('.gz')
             else Path(filepath) for filepath in filepaths]
    total = sum(path.stat().st_size for path in paths)
    copied = 0
    started = reported = time.monotonic()
    # These files were already written by assign_taxids/masking with LF newlines.
    # Avoid text decoding/encoding and 64-KiB Python read/write calls on DrvFS.
    with open(library_filename, 'wb') as out_file:
        for path in paths:
            with open(path, 'rb') as in_file:
                while block := in_file.read(8 * 1024 * 1024):
                    out_file.write(block)
                    copied += len(block)
                    now = time.monotonic()
                    if now - reported >= 60:
                        rate = copied / max(now - started, 0.001)
                        log.info(
                            'Merging {}: {:.1f}/{:.1f} GiB ({:.1f}%), '
                            '{:.1f} MiB/s, ETA {:.1f} h\n'.format(
                                library_filename, copied / 2**30, total / 2**30,
                                100 * copied / total if total else 100,
                                rate / 2**20, max(total - copied, 0) / rate / 3600))
                        reported = now


def patched_source(path):
    source = path.read_text()
    if 'SCRIPT_VERSION = "2.17.1"' not in source:
        raise RuntimeError('This download compatibility layer requires k2 2.17.1')
    replacements = [
        ('import shutil\n',
         'import shutil\nfrom reference_layout import source_path as metag_source_path\n'
         'from k2_compat import merge_library as metag_merge_library\n'
         'from k2_compat import write_taxid_blocks as metag_write_taxid_blocks\n'),
        ('''    with open(library_filename, "w") as out_file:
        for filepath in filepaths:
            if filepath.endswith(".gz"):
                filepath = os.path.splitext(filepath)[0]
            with open(filepath, "r") as in_file:
                shutil.copyfileobj(in_file, out_file)''',
         '    metag_merge_library(filepaths, library_filename, LOG)'),
        ('        tmp_local_name = local_name + ".tmp"',
         '        local_name = metag_source_path(local_name)\n        tmp_local_name = local_name + ".tmp"'),
        ('if not os.path.exists(filepath.strip()):',
         'if not os.path.exists(metag_source_path(filepath.strip())):'),
        ('with opener(absolute_filepath, "rt") as in_file:',
         'with opener(metag_source_path(absolute_filepath), "rt") as in_file:'),
        ('with open(args.library, "r") as in_file:',
         'with open(metag_source_path(os.path.abspath(args.library)), "r") as in_file:'),
        ('remote_path = ftp_path + "/" + os.path.basename(ftp_path) + suffix',
         'ftp_path = ftp_path.rstrip("/")\n        remote_path = ftp_path + "/" + os.path.basename(ftp_path) + suffix'),
        # --resume trusts independently checked, completed source downloads.
        # Check the actual storage path, not the logical path under library/.
        ('                    abs_filepath = os.path.abspath(filepath.strip())',
         '                    abs_filepath = metag_source_path(os.path.abspath(filepath.strip()))'),
        ('if not os.path.exists(abs_filepath):',
         'if (not os.path.isfile(abs_filepath) or os.path.getsize(abs_filepath) == 0\n'
         '                            or os.path.exists(abs_filepath + ".tmp")\n'
         '                            or os.path.exists(abs_filepath + ".aria2")):'),
        ('(done, not_done) = concurrent.futures.wait(futures)\n            if len(not_done) != 0:',
         '(done, not_done) = concurrent.futures.wait(futures)\n            for completed in done:\n                completed.result()\n            if len(not_done) != 0:'),
        ('elif response.status == 404:\n                    LOG.warning(',
         'elif response.status in (408, 429, 500, 502, 503, 504):\n                    response.read()\n                    raise http.client.HTTPException("Temporary HTTP status " + str(response.status))\n                elif response.status == 404:\n                    response.read()\n                    raise RuntimeError("NCBI file not found; refusing empty reference: " + url)\n                    LOG.warning('),
        ('''                    for checksum in checksums:
                        (md5sum, remote_filename) = checksum.split()
                        remote_filename = os.path.basename(
                            remote_filename.decode()
                        )
                        md5[remote_filename] = md5sum.decode()''',
         '''                    checked_md5 = {}
                    for checksum in checksums:
                        if not checksum.strip():
                            continue
                        # md5sum filenames may contain spaces; error pages are not checksums.
                        fields = checksum.split(maxsplit=1)
                        if len(fields) != 2 or not re.fullmatch(rb"[0-9a-fA-F]{32}", fields[0]):
                            raise http.client.HTTPException("Malformed NCBI checksum response for " + url)
                        try:
                            remote_filename = fields[1].rstrip(b"\\r\\n").removeprefix(b"*").decode()
                        except UnicodeError as e:
                            raise http.client.HTTPException("Invalid NCBI checksum filename for " + url) from e
                        # Checksums are relative to remote_dirname. Subdirectories
                        # legitimately contain different files with the same basename.
                        # Only direct children can match this download/cache key.
                        while remote_filename.startswith("./"):
                            remote_filename = remote_filename[2:]
                        if "/" in remote_filename:
                            continue
                        if not remote_filename:
                            raise http.client.HTTPException("Empty NCBI checksum filename for " + url)
                        md5sum = fields[0].decode().lower()
                        if remote_filename in checked_md5 and checked_md5[remote_filename] != md5sum:
                            raise http.client.HTTPException("Conflicting NCBI checksums for " + url)
                        checked_md5[remote_filename] = md5sum
                    # Publish only after the whole response has passed validation.
                    md5.update(checked_md5)'''),
        ('            if not skip_md5_check and os.path.exists(local_name)\\',
         '            if filename.endswith("_genomic.fna.gz") and filename not in md5:\n                raise http.client.HTTPException("NCBI genome checksum unavailable: " + url)\n            if not skip_md5_check and os.path.exists(local_name)\\'),
        ('skip_md5_check = False\n    while i < num_urls:',
         'skip_md5_check = False\n    download_failures = {}\n    while i < num_urls:'),
        ('except (http.client.HTTPException, http.client.RemoteDisconnected) as e:',
         'except (http.client.HTTPException, OSError) as e:'),
        ('conn.close()\n            time.sleep(0.1)\n            conn.connect()\n            continue',
         'conn.close()\n            download_failures[url] = download_failures.get(url, 0) + 1\n            if download_failures[url] >= 6:\n                raise RuntimeError("Download failed after six attempts: " + url) from e\n            time.sleep(min(2 ** download_failures[url], 30))\n            continue'),
        ('            masker.wait()\n\n    return (sequences_added, ch_added, sequence_to_url)',
         '            if masker.wait() != 0:\n                raise RuntimeError("Reference masking failed: " + filepath)\n\n    return (sequences_added, ch_added, sequence_to_url)'),
        ('                skip_md5_check = True\n                conn.connect()\n                continue',
         '                skip_md5_check = True\n                download_failures[url] = download_failures.get(url, 0) + 1\n                if download_failures[url] >= 6:\n                    raise RuntimeError("Repeated NCBI checksum mismatch: " + url)\n                conn.close()\n                time.sleep(min(2 ** download_failures[url], 30))\n                continue'),
    ]
    for before, after in replacements:
        if source.count(before) != 1:
            raise RuntimeError('Unexpected k2 source; download compatibility patch was not applied')
        source = source.replace(before, after)
    # Replace only the installed version's checked assign_taxids loop. Downloads,
    # low-complexity masking, its failure check, and merged libraries are unchanged.
    start_marker = '        with opener(metag_source_path(absolute_filepath), "rt") as in_file:\n'
    end_marker = '        if out_filepath.endswith(".tmp"):\n'
    if source.count(start_marker) != 1 or source.count(end_marker) != 1:
        raise RuntimeError('Unexpected k2 source; batched FASTA rewrite was not applied')
    start = source.index(start_marker) + len(start_marker)
    end = source.index(end_marker, start)
    source = (source[:start] +
              '            sequences_added, ch_added, sequence_to_url = metag_write_taxid_blocks(\n'
              '                in_file, masker.stdin if masker else out_file, filepath,\n'
              '                manifest_to_taxid, accession_to_taxid, filepath_to_url,\n'
              '                binary=bool(masker))\n' + source[end:])
    # A frozen selection supplies the assembly summary before the upstream downloader
    # creates its manifest. The same checked genome-download/masking code is reused.
    before = 'def download_genomic_library(args):\n'
    if source.count(before) != 1:
        raise RuntimeError('Unexpected k2 source; selection hook was not applied')
    source = source.replace(before, before + '    from reference_selection import summary_for_k2\n'
                            '    _metag_summary = summary_for_k2(args)\n')
    start = source.index('        remote_dir_name = args.library\n')
    end = source.index('        if args.library == "human":\n            with open(', start)
    download_summary = source[start:end]
    source = (source[:start] + '        if _metag_summary is not None:\n'
              '            shutil.copyfile(_metag_summary, "assembly_summary.txt")\n'
              '        else:\n' + textwrap.indent(download_summary, '    ') + source[end:])
    # Optional resumable genome batches. Legacy callers/tests keep their exact path.
    source = source.replace('    _metag_summary = summary_for_k2(args)\n',
                            '    from reference_batches import frozen_summary, freeze_current_summary\n'
                            '    _metag_summary = frozen_summary(args, summary_for_k2(args))\n')
    before = '        with open("assembly_summary.txt", "r") as f:\n'
    if source.count(before) != 1:
        raise RuntimeError('Unexpected k2 source; assembly freeze hook was not applied')
    source = source.replace(before, '        freeze_current_summary(args)\n' + before)
    before = ('def assign_taxid_to_sequences(args,\n'
              '                              manifest_to_taxid,\n'
              '                              accession_to_taxid={},\n'
              '                              filepath_to_url={}):\n')
    if source.count(before) != 1:
        raise RuntimeError('Unexpected k2 source; genome batch hook was not applied')
    source = source.replace(before, before +
        '    from reference_batches import enabled, prepare_library\n'
        '    if enabled(args):\n'
        '        return prepare_library(args, manifest_to_taxid, accession_to_taxid,\n'
        '                               filepath_to_url, scan_fasta_file, LOG)\n')
    # A batch publishes its own map and manifest, with no merged FASTA to scan.
    start = source.index('def download_genomic_library(args):\n')
    begin = source.index('        with open(library_filename, "r") as in_file:\n', start)
    end = source.index('    elif args.library in ["plasmid", "plastid", "mitochondrion"]:', begin)
    block = source[begin:end]
    source = (source[:begin] + '        if sequence_to_url is not None:\n' +
              textwrap.indent(block, '    ') + source[end:])
    return source


if __name__ == '__main__':
    executable = shutil.which('k2')
    if not executable:
        raise SystemExit('k2 is missing from the metag environment')
    original = Path(executable).resolve()
    # Compiling against the installed path keeps k2's binary lookup working.
    # Use the real __main__ globals so its worker functions remain picklable.
    code = compile(patched_source(original), str(original), 'exec')
    __file__ = str(original)
    exec(code, globals())

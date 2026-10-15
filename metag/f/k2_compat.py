#!/usr/bin/env python3
"""Run installed k2 2.17.1 with checked fixes for NCBI download failures.

Keep the installed package untouched. Refuse unfamiliar source rather than
silently applying a partial patch. Original assembly summaries remain intact.
"""
import shutil
import textwrap
from pathlib import Path


def patched_source(path):
    source = path.read_text()
    if 'SCRIPT_VERSION = "2.17.1"' not in source:
        raise RuntimeError('This download compatibility layer requires k2 2.17.1')
    replacements = [
        ('remote_path = ftp_path + "/" + os.path.basename(ftp_path) + suffix',
         'ftp_path = ftp_path.rstrip("/")\n        remote_path = ftp_path + "/" + os.path.basename(ftp_path) + suffix'),
        ('if not os.path.exists(abs_filepath):',
         'if not os.path.isfile(abs_filepath) or os.path.getsize(abs_filepath) == 0:'),
        ('if resume:\n                nonexistent_filepaths = []',
         'if False:  # Revalidate cached files against NCBI MD5, not existence alone.\n                nonexistent_filepaths = []'),
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
                        remote_filename = os.path.basename(remote_filename)
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

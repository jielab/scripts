#!/usr/bin/env python3
"""Run installed k2 2.17.1 with checked fixes for NCBI download failures.

Keep the installed package untouched. Refuse unfamiliar source rather than
silently applying a partial patch. Original assembly summaries remain intact.
"""
import shutil
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

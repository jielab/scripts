"""Streaming file checksums shared by project integrity checks."""
import hashlib

def sha256_file(path, *, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(chunk_size), b''):
            digest.update(block)
    return digest.hexdigest()

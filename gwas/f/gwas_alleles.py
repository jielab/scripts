"""Indexed FASTA access and sequence allele normalization (1-based positions)."""
from pathlib import Path


class IndexedFasta:
    def __init__(self, path):
        self.handle = open(path, 'rb')
        self.index = {}
        for line in Path(str(path) + '.fai').read_text().splitlines():
            name, length, offset, bases, width = line.split()[:5]
            self.index[name] = tuple(map(int, (length, offset, bases, width)))

    def fetch(self, chrom, pos, length):
        name = {'23': 'X', '24': 'Y', '25': 'MT'}.get(str(chrom), str(chrom))
        if name not in self.index:
            name = 'chr' + ('M' if name == 'MT' else name)
        size, offset, bases, width = self.index[name]
        if pos < 1 or length < 1 or pos + length - 1 > size:
            raise ValueError(f'FASTA interval out of bounds: {chrom}:{pos}+{length}')
        start = pos - 1
        byte_start = offset + start // bases * width + start % bases
        end = start + length - 1
        byte_end = offset + end // bases * width + end % bases
        self.handle.seek(byte_start)
        return self.handle.read(byte_end - byte_start + 1).replace(b'\n', b'').replace(b'\r', b'').decode().upper()

    def close(self):
        self.handle.close()


def normalize(pos, ref, alt, fasta=None, chrom=None):
    """Trim identical sequence and, with FASTA, left-align repeat indels."""
    if ref == alt or not ref or not alt or not set(ref + alt) <= set('ACGT'):
        raise ValueError('Expected two distinct DNA alleles')
    while ref[-1] == alt[-1]:
        if min(len(ref), len(alt)) == 1:
            if fasta is None or pos <= 1 or len(ref) == len(alt):
                break
            previous = fasta.fetch(chrom, pos - 1, 1)
            if previous not in 'ACGT':
                break
            ref, alt, pos = previous + ref, previous + alt, pos - 1
        ref, alt = ref[:-1], alt[:-1]
    while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def id_alleles(snp, pos):
    fields = snp.split(':')
    if len(fields) != 4 or fields[1] != str(pos):
        return None
    ref, alt = fields[2:]
    if not ref or not alt or ref == alt or not set(ref + alt) <= set('ACGT'):
        return None
    return ref, alt

"""Review an exact Kraken species label using the existing per-mate hit lists."""
from array import array
from collections import Counter
from functools import lru_cache
from pathlib import Path
import sys


class MateVotes:
    def __init__(self, line):
        self.line = line
        self.values = array('Q')

    @lru_cache(maxsize=131072)
    def shared_ancestor(self, candidate, hit):
        return self.line.common([candidate, hit])

    def support(self, candidate, text, chain):
        # A run length counts k-mer positions, including positions sharing a minimizer.
        # It is never interpreted as a number of independent minimizers.
        counts = Counter()
        for token in text.split():
            source, length = token.split(b':')
            if source not in (b'A', b'0'):
                counts[self.shared_ancestor(candidate, int(source))] += int(length)
        total = sum(counts.values())
        if not total:
            return None
        below = 0
        for tid in chain:
            below += counts[tid]
            if below >= 3 and 2*below > total:
                return tid
        return 1

    def review(self, ordinal, candidate, row):
        reviewed = candidate
        if self.line.tax.rank.get(candidate) == 'species':
            left, right = row.rsplit(b'\t', 1)[1].strip().split(b'|:|')
            chain = self.line.chain(candidate)
            a, b = (self.support(candidate, end, chain) for end in (left, right))
            # A silent mate permits a provisional higher-rank call on the informative mate.
            if a is None or b is None:
                supported = a or b
                reviewed = self.line.common([self.line.parent_rank(candidate), supported]) if supported else 0
            else:
                reviewed = self.line.common([a, b])
        if ordinal > len(self.values):
            self.values.extend(array('Q', [0]) * max(65536, ordinal-len(self.values)))
        self.values[ordinal-1] = reviewed

    def finish(self, path, count):
        del self.values[count:]
        if len(self.values) != count:
            raise RuntimeError('Mate review count mismatch')
        if sys.byteorder != 'little':
            self.values.byteswap()
        with Path(path).open('wb') as handle:
            self.values.tofile(handle)
        self.values = array('Q')

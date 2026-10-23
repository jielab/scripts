"""Local, reciprocal, single-nucleotide GRCh37/38 harmonization.

No rsID-only joins, network downloads, indel lift-over, or strand guessing.
Reference changes between assemblies are retained as nucleotide identities.
"""
from __future__ import annotations

from gnomad_io import FastaContext, canonical_chrom

BASES = frozenset("ACGT")
AUTOSOMES = frozenset(str(i) for i in range(1, 23))
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def orient(base, strand):
    if base not in BASES or strand not in ("+", "-"):
        raise ValueError("orient requires one canonical base and + or - strand")
    return base if strand == "+" else base.translate(COMPLEMENT)


class MappingError(ValueError):
    def __init__(self, reason, detail=""):
        self.reason = reason
        super().__init__(reason + (": " + detail if detail else ""))


class GenomeMapper:
    """Map a validated source base to both builds; all reported positions are 1-based.

    Both chains are deliberately local files. The caller supplies GRCh38 and
    GRCh37/hs37d5 FASTAs, including indexes. Reject all ambiguous/non-reciprocal
    mappings rather than accepting a highest-scoring chain by default.
    """
    def __init__(self, fasta38, fasta37, chain38to37, chain37to38):
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise ValueError("Ancient integration requires pyliftover; pip install pyliftover==0.4.1") from exc
        if not chain38to37 or not chain37to38:
            raise ValueError("Both local reciprocal chain files are required")
        self.f38 = FastaContext(fasta38)
        try:
            self.f37 = FastaContext(fasta37)
            self.to37 = LiftOver(str(chain38to37))
            self.to38 = LiftOver(str(chain37to38))
        except BaseException:
            self.f38.close()
            if hasattr(self, "f37"):
                self.f37.close()
            raise

    @staticmethod
    def _unique(chain, chrom, pos):
        hits = chain.convert_coordinate("chr" + canonical_chrom(chrom), pos - 1, "+") or []
        # Two chains with the same coordinate still constitute an ambiguous map.
        if not hits:
            raise MappingError("unmapped")
        if len(hits) != 1:
            raise MappingError("ambiguous_mapping")
        c, p, strand, _ = hits[0]
        c = canonical_chrom(c)
        if c not in AUTOSOMES:
            raise MappingError("non_autosomal_target")
        if strand not in ("+", "-"):
            raise MappingError("invalid_mapping_strand")
        return c, p + 1, strand

    def map_coordinates(self, chrom, pos, source_build="GRCh38"):
        """Return reciprocal coordinates and strand without reading either FASTA.

        This is an exact candidate prefilter, not allele harmonization. A result
        must subsequently pass ``map_from_coordinates`` with the modern source
        REF before it can contribute to an intersection. Reversing the direction
        of this reciprocal, unique mapping preserves the accepted coordinate set.
        """
        chrom = canonical_chrom(chrom)
        if chrom not in AUTOSOMES:
            raise MappingError("non_autosomal_source")
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        if not isinstance(pos, int) or pos < 1:
            raise MappingError("invalid_source_position")
        forward, backward = ((self.to37, self.to38) if source_build == "GRCh38"
                             else (self.to38, self.to37))
        target_c, target_p, strand = self._unique(forward, chrom, pos)
        back_c, back_p, back_strand = self._unique(backward, target_c, target_p)
        if (back_c, back_p) != (chrom, pos) or back_strand != strand:
            raise MappingError("non_reciprocal")
        if source_build == "GRCh38":
            c38, p38, c37, p37 = chrom, pos, target_c, target_p
        else:
            c38, p38, c37, p37 = target_c, target_p, chrom, pos
        return dict(chrom38=c38, pos38=p38, chrom37=c37, pos37=p37,
                    strand37to38=strand)

    def _source_context(self, chrom, pos, ref, source_build):
        if chrom not in AUTOSOMES:
            raise MappingError("non_autosomal_source")
        if ref not in BASES:
            raise MappingError("non_snv_source")
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        source_fasta = self.f38 if source_build == "GRCh38" else self.f37
        try:
            context = source_fasta.context(chrom, pos)
        except (KeyError, ValueError) as exc:
            raise MappingError("source_fasta_unavailable", str(exc)) from exc
        if context[1] != ref:
            raise MappingError("source_ref_mismatch", f"{chrom}:{pos} {ref}!={context[1]}")
        return context

    def _with_reference(self, coordinates, ref, source_build, source_context):
        c38, p38 = coordinates["chrom38"], coordinates["pos38"]
        c37, p37 = coordinates["chrom37"], coordinates["pos37"]
        strand = coordinates["strand37to38"]
        to38 = "+" if source_build == "GRCh38" else strand
        try:
            context38 = (source_context if source_build == "GRCh38" else self.f38.context(c38, p38))
            context37 = (source_context if source_build == "GRCh37" else self.f37.context(c37, p37))
        except (KeyError, ValueError) as exc:
            raise MappingError("target_fasta_unavailable", str(exc)) from exc
        ref38, ref37 = context38[1], context37[1]
        if ref38 not in BASES or ref37 not in BASES:
            raise MappingError("noncanonical_reference")
        canonical_context = (context38 if ref38 in "CT" else
                             context38.translate(COMPLEMENT)[::-1])
        return dict(chrom38=c38, pos38=p38, ref38=ref38, chrom37=c37,
                    pos37=p37, ref37=ref37, strand37to38=strand,
                    source_to38=to38, source_ref38=orient(ref, to38),
                    context38=context38, canonical_context=canonical_context,
                    reference_changed=int(orient(ref37, strand) != ref38),
                    site_id=f"{c38}:{p38}:{ref38}")

    def map_from_coordinates(self, coordinates, ref, source_build="GRCh38"):
        """Validate both FASTAs using coordinates returned by this mapper.

        The caller must match the modern coordinate to the appropriate build in
        ``coordinates`` exactly. This method deliberately makes no chain queries;
        it is only for reusing the reciprocal prefilter from ``map_coordinates``.
        """
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        suffix = "38" if source_build == "GRCh38" else "37"
        source_context = self._source_context(coordinates["chrom" + suffix],
                                              coordinates["pos" + suffix], ref, source_build)
        return self._with_reference(coordinates, ref, source_build, source_context)

    def map(self, chrom, pos, ref, source_build="GRCh38"):
        # Keep source-REF validation first for the reference modern-first path.
        chrom = canonical_chrom(chrom)
        source_context = self._source_context(chrom, pos, ref, source_build)
        coordinates = self.map_coordinates(chrom, pos, source_build)
        return self._with_reference(coordinates, ref, source_build, source_context)

    def close(self):
        self.f38.close()
        self.f37.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

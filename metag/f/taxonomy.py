"""NCBI lineage interpretation, independent of any classifier installation."""
from pathlib import Path
import re

GROUPS = ("human", "mammal", "other_vertebrate", "invertebrate", "plant",
          "fungi", "bacteria", "archaea", "virus", "other_eukaryote", "unresolved")


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")[:100] or "unnamed"


class Taxonomy:
    def __init__(self, folder):
        folder = Path(folder)
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        with (folder / "nodes.dmp").open(encoding="utf-8") as handle:
            for line in handle:
                fields = [x.strip() for x in line.split("|")]
                if len(fields) < 3:
                    raise RuntimeError("Malformed NCBI nodes.dmp")
                tid, parent = int(fields[0]), int(fields[1])
                if tid <= 0 or parent <= 0 or tid in self.parent:
                    raise RuntimeError(f"Invalid/duplicate NCBI taxid {tid}")
                self.parent[tid], self.rank[tid] = parent, fields[2]
        if self.parent.get(1) != 1:
            raise RuntimeError("NCBI taxonomy has no valid root")
        with (folder / "names.dmp").open(encoding="utf-8") as handle:
            for line in handle:
                fields = [x.strip() for x in line.split("|")]
                if len(fields) >= 4 and fields[3] == "scientific name":
                    self.name[int(fields[0])] = fields[1]
        merged = folder / "merged.dmp"
        if merged.is_file():
            with merged.open(encoding="utf-8") as handle:
                for line in handle:
                    fields = [x.strip() for x in line.split("|")]
                    if len(fields) >= 2 and fields[0]:
                        self.merged[int(fields[0])] = int(fields[1])
        self._group_anchors = {9606:"human",40674:"mammal",7742:"other_vertebrate",
                               33208:"invertebrate",33090:"plant",4751:"fungi",2:"bacteria",
                               2157:"archaea",10239:"virus",2759:"other_eukaryote"}
        self._anchor_lineages = {tid:set(self.lineage(tid)) for tid in self._group_anchors if tid in self.parent}
        self._cache = {0: {"taxid": 0, "species_taxid": 0, "rank": "none",
                          "group": "unresolved", "name": "unclassified", "status": "unclassified"}}

    def resolve(self, taxid):
        taxid = int(taxid)
        seen = set()
        while taxid in self.merged:
            if taxid in seen:
                raise RuntimeError("Cycle in merged NCBI taxonomy")
            seen.add(taxid)
            taxid = self.merged[taxid]
        if taxid and taxid not in self.parent:
            raise RuntimeError(f"Classification taxid absent from database taxonomy: {taxid}")
        return taxid

    def lineage(self, taxid):
        current, seen = self.resolve(taxid), set()
        while current:
            if current in seen or current not in self.parent:
                raise RuntimeError(f"Broken/cyclic NCBI lineage at {current}")
            seen.add(current)
            yield current
            parent = self.parent[current]
            if parent == current:
                return
            current = parent

    def describe(self, taxid):
        if taxid in self._cache:
            return self._cache[taxid]
        original = taxid
        taxid = self.resolve(taxid)
        lineage = tuple(self.lineage(taxid))
        ancestors = set(lineage)
        species = next((x for x in lineage if self.rank.get(x) == "species"), 0)
        assigned = species or taxid
        if 9606 in ancestors:
            group = "human"
        elif 40674 in ancestors:
            group = "mammal"
        elif 7742 in ancestors:
            group = "other_vertebrate"
        elif 33208 in ancestors:
            group = "invertebrate"
        elif 33090 in ancestors:
            group = "plant"
        elif 4751 in ancestors:
            group = "fungi"
        elif 2 in ancestors:
            group = "bacteria"
        elif 2157 in ancestors:
            group = "archaea"
        elif 10239 in ancestors:
            group = "virus"
        elif 2759 in ancestors:
            group = "other_eukaryote"
        else:
            group = "unresolved"
        # A parent clade can contain several display groups. Lack of a
        # descendant in its ancestor chain cannot establish invertebrate status.
        if not species and any(taxid in chain and self._group_anchors[anchor] != group
                               for anchor, chain in self._anchor_lineages.items()):
            group = "unresolved"
        result = {"human_ancestor": taxid in self._anchor_lineages.get(9606, set()),
                  "vertebrate_lineage": 7742 in ancestors,
                  "taxid": assigned, "species_taxid": species,
                  "rank": self.rank.get(assigned, "none"), "group": group,
                  "name": self.name.get(assigned, str(assigned)),
                  "status": "assigned_species" if species else "higher_rank"}
        self._cache[original] = result
        return result

    def species(self, taxid):
        return self.describe(taxid)["species_taxid"]

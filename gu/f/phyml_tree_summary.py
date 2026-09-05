#!/usr/bin/env python3
"""Collect PhyML Newick/stat outputs into a table consumable by normalize/Shiny."""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

from csv_compat import enable_wide_csv_fields


enable_wide_csv_fields()

def safe_label(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return value or "tip"


def read_label_map(path: Path) -> dict[str, str]:
    with path.open() as handle:
        return {
            row["phy_label"]: safe_label(row["label"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def relabel_newick(text: str, labels: dict[str, str]) -> str:
    for short in sorted(labels, key=len, reverse=True):
        text = re.sub(
            rf"(?<=[(,]){re.escape(short)}(?=[:),;])",
            labels[short],
            text,
        )
    return "".join(text.split())


def bootstrap_values(newick: str) -> list[float]:
    return [float(x) for x in re.findall(r"\)([0-9]+(?:\.[0-9]+)?)(?=[:),;])", newick)]


@dataclass
class NewickNode:
    label: str = ""
    children: list["NewickNode"] = field(default_factory=list)
    support: float | None = None
    node_id: str = ""


def parse_newick(text: str) -> NewickNode:
    """Parse the simple, unquoted Newick emitted by PhyML.

    Keeping this tiny parser local avoids adding a Biopython dependency to the
    normalization path.  Tip labels have already been sanitized by safe_label.
    """
    source = "".join(text.split())
    index = 0
    serial = 0

    def token() -> str:
        nonlocal index
        start = index
        while index < len(source) and source[index] not in ":,();":
            index += 1
        return source[start:index]

    def branch_length() -> None:
        nonlocal index
        if index < len(source) and source[index] == ":":
            index += 1
            while index < len(source) and source[index] not in ",();":
                index += 1

    def subtree() -> NewickNode:
        nonlocal index, serial
        if index >= len(source):
            raise ValueError("unexpected end of Newick")
        if source[index] != "(":
            label = token()
            if not label:
                raise ValueError("empty Newick tip")
            node = NewickNode(label=label)
            branch_length()
            return node
        index += 1
        children = [subtree()]
        while index < len(source) and source[index] == ",":
            index += 1
            children.append(subtree())
        if index >= len(source) or source[index] != ")":
            raise ValueError("unterminated Newick clade")
        index += 1
        label = token()
        support = None
        try:
            support = float(label) if label else None
        except ValueError:
            pass
        serial += 1
        node = NewickNode(label=label, children=children, support=support, node_id=f"N{serial}")
        branch_length()
        return node

    root = subtree()
    if index < len(source) and source[index] == ";":
        index += 1
    if index != len(source):
        raise ValueError(f"unexpected Newick content at offset {index}")
    return root


def archaic_lineage(label: str) -> str:
    if re.search(r"Altai|Chagyr|Vindija|Neander", label, re.I):
        return "Neanderthal"
    if re.search(r"Denis", label, re.I):
        return "Denisovan"
    if re.search(r"archaic", label, re.I):
        return "Archaic/Other"
    return ""


def candidate_clade(newick: str, min_modern_tips: int = 2) -> dict[str, object]:
    """Find a bootstrap-supported archaic/modern edge without match preselection.

    ``direct_match_pass`` is deliberately not used here: with a small number
    of callable SNPs, a permissive pairwise threshold can select most modern
    haplotypes and make the tree statistic circular.  Instead, both sides of
    every internal edge are enumerated.  An eligible side contains modern
    tips plus one phylogenetically coherent archaic lineage.  We select the
    highest bootstrap, then the side with fewer modern tips, then the side
    containing more references from that lineage.

    PhyML emits an unrooted topology unless an outgroup is supplied.  Testing
    both sides keeps the result invariant to the arbitrary serialization root.
    """
    blank = {
        "candidate_lineage": "", "candidate_clade_pass": 0,
        "candidate_clade_rule": "tree_edge_same_archaic_lineage",
        "tree_has_ancestral_outgroup": 0,
        "candidate_clade_bootstrap": "", "candidate_clade_node": "",
        "candidate_clade_side": "", "candidate_clade_n_tips": "",
        "candidate_clade_modern_tips": "", "candidate_clade_archaic_tips": "",
        "candidate_clade_specificity": "", "candidate_clade_tips": "",
    }
    if not newick:
        return blank
    root = parse_newick(newick)
    descendants: dict[str, set[str]] = {}
    internal: list[NewickNode] = []

    def visit(node: NewickNode) -> set[str]:
        if not node.children:
            return {node.label}
        tips = set().union(*(visit(child) for child in node.children))
        descendants[node.node_id] = tips
        internal.append(node)
        return tips

    all_tips = visit(root)
    has_ancestral = "Ancestral" in all_tips
    blank["tree_has_ancestral_outgroup"] = int(has_ancestral)
    candidates = []
    for node in internal:
        if node is root:
            continue
        down = descendants[node.node_id]
        for side_name, side in (("descendant", down), ("complement", all_tips - down)):
            if has_ancestral and "Ancestral" in side:
                continue
            lineages = {archaic_lineage(x) for x in side if archaic_lineage(x)}
            archaic = {x for x in side if archaic_lineage(x)}
            modern = side - archaic
            if len(lineages) != 1 or len(modern) < min_modern_tips or len(side) >= len(all_tips):
                continue
            support = node.support if node.support is not None else -1.0
            candidates.append((-support, len(modern), -len(archaic), side_name, node, side, archaic, modern, next(iter(lineages))))
    if not candidates:
        return blank
    _, _, _, side_name, node, side, archaic, modern, lineage = min(candidates)
    bootstrap = node.support if node.support is not None else None
    return {
        "candidate_lineage": lineage,
        "tree_has_ancestral_outgroup": int(has_ancestral),
        "candidate_clade_pass": int(bootstrap is not None and bootstrap >= 70),
        "candidate_clade_rule": "tree_edge_same_archaic_lineage;bootstrap>=70;modern_tips>=2",
        "candidate_clade_bootstrap": "" if bootstrap is None else bootstrap,
        "candidate_clade_node": node.node_id,
        "candidate_clade_side": side_name,
        "candidate_clade_n_tips": len(side),
        "candidate_clade_modern_tips": len(modern),
        "candidate_clade_archaic_tips": len(archaic),
        "candidate_clade_specificity": "",
        "candidate_clade_tips": ",".join(sorted(side)),
    }


def phy_inputs(loci_root: Path) -> list[Path]:
    """Return PhyML inputs from flat or multi-window layouts."""
    flat = loci_root / "haplotypes.phy"
    if flat.exists():
        return [flat]
    return sorted(loci_root.glob("*/haplotypes.phy"))


def flat_locus_id(out: Path) -> str:
    summary = out / "final" / "loci.tsv"
    if not summary.is_file():
        raise SystemExit(f"ERROR: flat PhyML output is missing locus summary: {summary}")
    with summary.open() as handle:
        locus_ids = {
            row.get("locus_id", "").strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("locus_id", "").strip()
        }
    if len(locus_ids) != 1:
        raise SystemExit(
            "ERROR: flat PhyML output requires exactly one locus_id in "
            f"{summary}; found {len(locus_ids)}"
        )
    return next(iter(locus_ids))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--default-status",
        choices=("not_requested", "not_run"),
        help="write a row with this status when no current tree exists",
    )
    parser.add_argument(
        "--skip-render",
        action="store_true",
        help="summarize existing Newick/stat files without requiring R or creating a PNG",
    )
    parser.add_argument(
        "--reclassify-existing", action="store_true",
        help="recompute candidate-clade fields from final/trees.tsv when raw PHYLIP inputs are absent",
    )
    args = parser.parse_args()
    loci_root = args.out / "loci"
    if not args.default_status and not args.skip_render:
        rscript = shutil.which("Rscript")
        if not rscript:
            raise SystemExit("ERROR: Rscript is required to render the PhyML tree PNG")
        plotter = Path(__file__).with_name("phyml_tree_plot.R")
        subprocess.run([rscript, str(plotter), "--out", str(args.out)], check=True)
    rows = []
    for phy in phy_inputs(loci_root):
        locus_id = flat_locus_id(args.out) if phy.parent == loci_root else phy.parent.name
        tree = Path(f"{phy}_phyml_tree.txt")
        stats = Path(f"{phy}_phyml_stats.txt")
        plot = tree.with_suffix(".png")
        meta = phy.with_name("haplotypes.phy.meta.tsv")
        tree_exists = (
            tree.is_file()
            and tree.stat().st_size
            and meta.is_file()
        )
        if not tree_exists and args.default_status:
            rows.append(
                {
                    "locus_id": locus_id,
                    "tree_status": args.default_status,
                    "tree_file": "",
                    "stats_file": "",
                    "plot_file": "",
                    "n_bootstrap_nodes": "",
                    "bootstrap_min": "",
                    "bootstrap_median": "",
                    "bootstrap_max": "",
                    **candidate_clade(""),
                    "tree_newick": "",
                }
            )
            continue
        if not tree_exists:
            continue
        newick = relabel_newick(tree.read_text(errors="replace"), read_label_map(meta))
        support = bootstrap_values(newick)
        clade = candidate_clade(newick)
        rows.append(
            {
                "locus_id": locus_id,
                "tree_status": "complete",
                "tree_file": str(tree.resolve()),
                "stats_file": str(stats.resolve()) if stats.is_file() else "",
                "plot_file": str(plot.resolve()) if plot.is_file() and plot.stat().st_size else "",
                "n_bootstrap_nodes": len(support),
                "bootstrap_min": min(support) if support else "",
                "bootstrap_median": median(support) if support else "",
                "bootstrap_max": max(support) if support else "",
                **clade,
                "tree_newick": newick,
            }
        )
    # A locus with insufficient LD/sequence information has no PHYLIP input,
    # but it is still a valid negative result and must remain visible to
    # normalize/Shiny.  Record why tree inference was not run.
    represented = {str(row.get("locus_id", "")) for row in rows}
    existing = args.out / "final" / "trees.tsv"
    if args.reclassify_existing and existing.is_file() and existing.stat().st_size:
        with existing.open() as handle:
            for old in csv.DictReader(handle, delimiter="\t"):
                locus_id = old.get("locus_id", "").strip()
                newick = old.get("tree_newick", "").strip()
                if not locus_id or locus_id in represented or not newick:
                    continue
                old.update(candidate_clade(newick))
                rows.append(old)
                represented.add(locus_id)
    skipped = args.out / "final" / "skipped_loci.tsv"
    if skipped.is_file():
        with skipped.open() as handle:
            for skipped_row in csv.DictReader(handle, delimiter="\t"):
                locus_id = skipped_row.get("locus_id", "").strip()
                if not locus_id or locus_id in represented:
                    continue
                status = re.sub(r"[^A-Za-z0-9_]+", "_", skipped_row.get("status", "skipped").strip()).strip("_") or "skipped"
                rows.append(
                    {
                        "locus_id": locus_id,
                        "tree_status": f"not_run_{status}",
                        "tree_file": "",
                        "stats_file": "",
                        "plot_file": "",
                        "n_bootstrap_nodes": "",
                        "bootstrap_min": "",
                        "bootstrap_median": "",
                        "bootstrap_max": "",
                        **candidate_clade(""),
                        "tree_newick": "",
                    }
                )
    final = args.out / "final"
    final.mkdir(parents=True, exist_ok=True)
    dest = final / "trees.tsv"
    header = [
        "locus_id",
        "tree_status",
        "tree_file",
        "stats_file",
        "plot_file",
        "n_bootstrap_nodes",
        "bootstrap_min",
        "bootstrap_median",
        "bootstrap_max",
        "candidate_lineage",
        "tree_has_ancestral_outgroup",
        "candidate_clade_pass",
        "candidate_clade_rule",
        "candidate_clade_bootstrap",
        "candidate_clade_node",
        "candidate_clade_side",
        "candidate_clade_n_tips",
        "candidate_clade_modern_tips",
        "candidate_clade_archaic_tips",
        "candidate_clade_specificity",
        "candidate_clade_tips",
        "tree_newick",
    ]
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"PHYML TREE SUMMARY: loci={len(rows)} output={dest}")


if __name__ == "__main__":
    main()

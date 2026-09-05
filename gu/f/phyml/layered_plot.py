#!/usr/bin/env python3
"""Render GU PhyML raw, region-first, and filtered candidate trees.

The renderer is deliberately a QC plotter, not a replacement for the circular
population-ring figure used in the MPB manuscript.  It uses only matplotlib
(already listed in gu/environment.yml), reads the role-rich layered manifests
when present, and always preserves the Newick/statistical tables as the source
of truth.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from csv_compat import enable_wide_csv_fields


enable_wide_csv_fields()

def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@dataclass
class Node:
    label: str = ""
    length: float = 0.0
    support: float | None = None
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    serial: int = 0
    x: float = 0.0
    y: float = 0.0
    descendant_tips: set[str] = field(default_factory=set)


def parse_newick(text: str) -> Node:
    source = "".join(text.split())
    index = 0
    serial = 0

    def token() -> str:
        nonlocal index
        start = index
        while index < len(source) and source[index] not in ":,();":
            index += 1
        return source[start:index]

    def branch_length() -> float:
        nonlocal index
        if index >= len(source) or source[index] != ":":
            return 0.0
        index += 1
        value = token()
        try:
            parsed = float(value)
            return parsed if math.isfinite(parsed) and parsed >= 0 else 0.0
        except ValueError:
            return 0.0

    def subtree() -> Node:
        nonlocal index, serial
        if index >= len(source):
            raise ValueError("unexpected end of Newick")
        serial += 1
        node_serial = serial
        if source[index] != "(":
            label = token()
            if not label:
                raise ValueError("empty Newick tip")
            return Node(label=label, length=branch_length(), serial=node_serial)
        index += 1
        children = [subtree()]
        while index < len(source) and source[index] == ",":
            index += 1
            children.append(subtree())
        if index >= len(source) or source[index] != ")":
            raise ValueError("unterminated Newick clade")
        index += 1
        node_label = token()
        support = None
        try:
            candidate = float(node_label) if node_label else None
            if candidate is not None and math.isfinite(candidate) and candidate >= 0:
                support = candidate
        except ValueError:
            pass
        node = Node(label=node_label, length=branch_length(), support=support,
                    children=children, serial=node_serial)
        for child in children:
            child.parent = node
        return node

    root = subtree()
    if index < len(source) and source[index] == ";":
        index += 1
    if index != len(source):
        raise ValueError(f"unexpected Newick content at offset {index}")
    root.length = 0.0
    return root


def tips(root: Node) -> list[Node]:
    out: list[Node] = []

    def visit(node: Node) -> None:
        if not node.children:
            out.append(node)
            return
        for child in node.children:
            visit(child)

    visit(root)
    return out


def internal_nodes(root: Node) -> list[Node]:
    out: list[Node] = []

    def visit(node: Node) -> None:
        if node.children:
            out.append(node)
            for child in node.children:
                visit(child)

    visit(root)
    return out


def layout(root: Node) -> list[Node]:
    leaves = tips(root)
    for i, leaf in enumerate(leaves):
        leaf.y = float(i)

    def assign(node: Node, parent_x: float = 0.0) -> set[str]:
        node.x = parent_x + (node.length if node.parent is not None else 0.0)
        if not node.children:
            node.descendant_tips = {node.label}
            return node.descendant_tips
        descendant: set[str] = set()
        for child in node.children:
            descendant.update(assign(child, node.x))
        node.y = sum(child.y for child in node.children) / len(node.children)
        node.descendant_tips = descendant
        return descendant

    assign(root)
    # All-zero/effectively-zero branch lengths are legal in bootstrap trees but
    # make a phylogram unreadable. Fall back to cladogram depth in that case.
    if max((leaf.x for leaf in leaves), default=0.0) <= 1e-12:
        def depth(node: Node, value: float = 0.0) -> None:
            node.x = value
            for child in node.children:
                depth(child, value + 1.0)
        depth(root)
    return leaves


def relabel_tree(root: Node, labels: dict[str, str]) -> None:
    for leaf in tips(root):
        leaf.label = labels.get(leaf.label, leaf.label)


def manifest_for_phy(phy: Path) -> Path | None:
    candidate = phy.with_name(phy.name.removesuffix(".phy") + ".manifest.tsv")
    return candidate if candidate.is_file() else None


def role_info(phy: Path, root: Node) -> dict[str, dict[str, str]]:
    manifest_path = manifest_for_phy(phy)
    rows = read_tsv(manifest_path) if manifest_path else []
    if rows:
        return {row.get("label", ""): row for row in rows if row.get("label", "")}
    out: dict[str, dict[str, str]] = {}
    for leaf in tips(root):
        label = leaf.label
        lower = label.lower()
        if label == "Ancestral":
            role, lineage = "outgroup", "Ancestral"
        elif re.search(r"altai|chagyr|vindija|neander", lower):
            role, lineage = "archaic", "Neanderthal"
        elif "denis" in lower:
            role, lineage = "archaic", "Denisovan"
        else:
            role, lineage = "modern", "Modern"
        out[label] = {"label": label, "role": role, "lineage": lineage, "n": ""}
    return out


def color_for(info: dict[str, str]) -> str:
    role = info.get("role", "")
    lineage = info.get("lineage", "")
    if role == "candidate":
        return "#b2182b"
    if role == "candidate_context":
        return "#ef8a62"
    if role in {"hard_control", "common_control", "control"}:
        return "#777777"
    if role == "archaic" and lineage == "Neanderthal":
        return "#2166ac"
    if role == "archaic" and lineage == "Denisovan":
        return "#8c510a"
    if role == "outgroup" or lineage == "Ancestral":
        return "#1b7837"
    if info.get("family_ids", ""):
        return "#762a83"
    return "#222222"


def tree_scope(phy: Path, info: dict[str, dict[str, str]]) -> str:
    values = [row.get("tree_scope", "") for row in info.values() if row.get("tree_scope", "")]
    if values:
        return values[0]
    name = phy.name
    if ".region_common" in name:
        return "region_common_unfiltered"
    if ".evidence." in name:
        return "candidate_filtered"
    return "all_unique_exploratory"


def draw_tree(phy: Path, dpi: int, max_labels: int, overwrite: bool) -> Path | None:
    tree_file = Path(f"{phy}_phyml_tree.txt")
    if not tree_file.is_file() or tree_file.stat().st_size == 0:
        return None
    output = tree_file.with_suffix(".png")
    if output.is_file() and output.stat().st_size and not overwrite:
        return output

    meta_rows = read_tsv(Path(f"{phy}.meta.tsv"))
    label_map = {
        row.get("phy_label", ""): row.get("label", "")
        for row in meta_rows if row.get("phy_label", "") and row.get("label", "")
    }
    root = parse_newick(tree_file.read_text(errors="replace"))
    relabel_tree(root, label_map)
    leaves = layout(root)
    info = role_info(phy, root)
    n_tip = len(leaves)
    if n_tip < 2:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    height = min(32.0, max(7.0, 0.20 * n_tip + 2.0))
    fig, ax = plt.subplots(figsize=(15.5, height))
    for node in internal_nodes(root):
        if not node.children:
            continue
        ys = [child.y for child in node.children]
        ax.plot([node.x, node.x], [min(ys), max(ys)], linewidth=0.55, color="#777777")
        for child in node.children:
            ax.plot([node.x, child.x], [child.y, child.y], linewidth=0.55, color="#777777")

    show_all = n_tip <= max_labels
    x_max = max((leaf.x for leaf in leaves), default=1.0)
    label_x = x_max + max(x_max * 0.015, 0.005)
    for leaf in leaves:
        leaf_info = info.get(leaf.label, {"label": leaf.label})
        color = color_for(leaf_info)
        role = leaf_info.get("role", "")
        important = role in {
            "candidate", "candidate_context", "archaic", "outgroup"
        } or bool(leaf_info.get("family_ids", ""))
        ax.scatter([leaf.x], [leaf.y], s=15 if important else 5, color=color, zorder=3)
        if show_all or important:
            n = leaf_info.get("n", "")
            suffix = f"  n={n}" if n not in {"", "NA", "1"} else ""
            ax.text(label_x, leaf.y, leaf.label + suffix, va="center", fontsize=6.8 if show_all else 7.5,
                    color=color, fontweight="bold" if important else "normal")

    if n_tip <= 300:
        for node in internal_nodes(root):
            if node.support is not None and node.support >= 70 and node is not root:
                ax.text(node.x, node.y, f"{node.support:g}", fontsize=6.2,
                        ha="right", va="bottom", color="#7f0000")

    scope = tree_scope(phy, info)
    expected = next((row.get("expected_lineage", "") for row in info.values()
                     if row.get("expected_lineage", "")), "")
    title = f"{phy.parent.name}: {scope}"
    if expected:
        title += f" ({expected})"
    ax.set_title(title, fontsize=12, loc="left")
    ax.set_xlabel("substitutions per site (or cladogram depth when branch lengths are zero)")
    ax.set_ylabel("haplotype groups")
    ax.set_ylim(-1, n_tip)
    ax.set_xlim(left=min(0.0, min((node.x for node in internal_nodes(root)), default=0.0)),
                right=label_x + max(x_max * 0.38, 0.12))
    ax.invert_yaxis()
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.grid(False)
    fig.text(
        0.01, 0.004,
        "QC rendering. Regional unfiltered trees show archaic affinity; filtered candidate trees provide the stricter evidence layer.",
        fontsize=7.0,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    tmp = output.with_name(f".{output.name}.tmp.png")
    fig.savefig(tmp, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    tmp.replace(output)
    return output


def phy_inputs(out: Path) -> Iterable[Path]:
    root = out / "loci"
    if not root.is_dir():
        return []
    return sorted(
        path for path in root.glob("**/haplotypes*.phy")
        if path.is_file() and path.stat().st_size and Path(f"{path}_phyml_tree.txt").is_file()
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--max-labels", type=int, default=180)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if not 72 <= args.dpi <= 600:
        fail("--dpi must be between 72 and 600")
    if args.max_labels < 0:
        fail("--max-labels must be non-negative")
    created = 0
    for phy in phy_inputs(args.out):
        try:
            result = draw_tree(phy, args.dpi, args.max_labels, args.overwrite)
            if result:
                created += 1
                print(f"PHYML LAYERED PLOT: {result}", flush=True)
        except Exception as exc:
            print(f"WARNING: failed to render {phy}: {exc}", flush=True)
    print(f"PHYML LAYERED PLOTS: available={created} output={args.out}", flush=True)


if __name__ == "__main__":
    main()

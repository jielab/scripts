#!/usr/bin/env python3
"""Shared Newick parser and bootstrap support for GWAS risk trees."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from comm import enable_wide_csv_fields


enable_wide_csv_fields()


# 🚩 Bootstrap support
def bootstrap_values(newick: str) -> list[float]:
    return [float(x) for x in re.findall(r"\)([0-9]+(?:\.[0-9]+)?)(?=[:),;])", newick)]


# 🚩 Newick tree structure and parsing
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



#!/usr/bin/env python3
"""Single command entry for GU's PhyML analysis stages."""
import importlib
import sys


COMMANDS = {
    "compare": "core",
    "region": "region_scan",
    "anchor": "anchor",
    "evidence": "evidence",
    "tree": "tree_summary",
    "plot": "layered_plot",
}


def main():
    if len(sys.argv) == 1 or sys.argv[1] in {"-h", "--help"}:
        print("Usage: phyml.py {compare,region,anchor,evidence,tree,plot} [options]\n"
              "Use COMMAND --help for stage options; run analyses with gu.sh phyml.")
        return
    command = sys.argv.pop(1)
    if command not in COMMANDS:
        raise SystemExit(f"Unknown PhyML stage: {command}")
    module = importlib.import_module(f"phyml_{COMMANDS[command]}")
    module.main()


if __name__ == "__main__":
    main()

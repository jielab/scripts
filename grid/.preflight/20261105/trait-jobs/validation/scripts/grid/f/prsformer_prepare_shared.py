#!/usr/bin/env python3
"""Compatibility entry for shared preparation; never trains PRSformer.

The infrastructure is public through 2.csx.sh/3.disco.sh/4.prsformer.sh prepare
and grid.sh. No method's fitted outputs are prerequisites.
"""
from shared_cohort import main

if __name__ == "__main__":
    main()

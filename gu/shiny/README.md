# GU Shiny

Run `./gu.sh normalize` and then `./gu.sh shiny`.

The app reads `gu.sqlite` directly below `GU_ANALYSIS_ROOT` and dynamically calculates build-filtered
overview, locus, regional-segment and per-sample summaries. It reconnects after
an atomic database refresh. `normalize/` contains lightweight summaries and
Shiny artifacts: its `as3/`, `trace/`, `phyml/`, and `ibdmix/` subdirectories
hold method-specific exports, while `summary/` holds cross-method tables. SQLite
records artifact paths relative to `GU_ANALYSIS_ROOT`. To run Shiny elsewhere,
copy `gu.sqlite` and `normalize/` together and keep them adjacent.

Published AS3 population callsets are stored separately in
`reference_callsets` with `reference_role=external_reference`. The Segments page
can overlay those intervals and summarize GU-AS3 overlap/reciprocal overlap by
1KG population. They never enter model inputs, `segments`, or individual
predictions. The Individuals page exposes raw-call, non-redundant interval-union,
and consensus-catalog burden definitions.

The Overview is QC-first: PhyML carriers come from the supported tree edge,
not the permissive pairwise identity flag. Native caller statistics remain
separate, and unavailable, unsupported, and exploratory methods are displayed
as distinct states. `consensus_catalog` includes only catalog tracts supported
by at least two distinct methods.

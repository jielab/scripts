# GU Shiny haplotype viewer v6

## What changed

- Corrects archaic-match values stored as fraction, percent, or percent multiplied by 100. Values are normalized to the valid 0–100% range; for example, `0.224`, `22.4`, and `2240` all display as `22.4%`.
- Removes the `[MATCH]` prefix from modern haplotype row labels.
- Appends 10 stable-random, non-matched 1KG phased haplotypes as negative controls. Their `SAMPLE | hapN` labels are red and they are separated at the bottom of the matrix.
- Adds an embedded IGV-Web panel synchronized to the selected locus, plus full-window IGV and UCSC links. The default view adds 250 kb on either side and the flank is adjustable.
- `./gu.sh normalize` now also writes:
  - `Rshiny/loci.browser.bed`
  - `Rshiny/loci.browser_links.tsv`

## Replace and run

Replace the existing `gu` directory with this directory, then run:

```bash
chmod +x gu.sh install.sh
./gu.sh normalize
./gu.sh shiny
```

The embedded browser requires internet access to `igv.org`; the UCSC button requires access to `genome.ucsc.edu`. If iframe embedding is restricted by the browser/network, use the full-window buttons.

## Negative-control definition

A control is sampled from all phased copies in `kg.tsv` after excluding:

1. the displayed matched copies; and
2. every copy present in the locus-level `hap_match.tsv` report.

If that strict pool is empty, v6 falls back to all non-displayed copies. Sampling is deterministic per locus.

# TCGA publication plots

Generated from the uncommitted `hiteshv3` working tree.

- Base commit: `d70d3a5554c2982e6bb5cfb67546bb0c0233a4ae`
- Tracked working-tree diff SHA-256: `b1ef0bd20eeff9927c48ec3844566bf0014cd186cb52ea42589926a41f3ed25a`
- Reference: `hg38`
- Width: `12` inches
- Aspect ratio: `0.5` per graph/cycle panel
- Interval offset: `0.08`
- Minimum coordinate width: `0.08`
- Gene font size: `9`
- PNG resolution: `300` DPI
- Combined graph/cycle output: enabled
- Displayed reconstructions: cycles `1`, `2`, and `3`
- Discordant-edge linewidth: read-count-scaled from `0.30` to `4.00` points
- CN and coverage axes: coupled by the graph-derived coverage-per-copy fit,
  with vertically synchronized ticks; both limits retain the earlier shared
  expansion required by the raw discordant-edge arc heights
- Combined layout: graph and cycle panels share the complete graph interval
  axis; interval guides bridge the compact inter-panel space; graph-panel gene
  and coordinate annotations are cropped so they appear only below the cycles
- Combined titles: `<sample>-Amplicon Graph` and `<sample>-Cycles`
- Cycle plot colors: antiquewhite/dim-gray genomic sequence edges, exact
  reference-image blue `#0000F5` for cycle SV edges, and orange path SV edges;
  the same colors are used in the legend
- Cycle reconstruction linewidth: `1.5` points for sequence outlines and SV
  connectors, including the matching legend swatches
- Cycle labels: larger bold cycle/path heading with a smaller CN value below
- Gene annotations: exact italic-label width measurement, dynamic
  collision-free lanes, and a taller shared graph/cycle gene track

Gene subsets:

- `TCGA-FX-A3NK-01A`: 12 genes in `inputs/TCGA-FX-A3NK-01A/gene_subset.txt`
- `TCGA-3B-A9HO-01A`: 15 genes in `inputs/TCGA-3B-A9HO-01A/gene_subset.txt`
- `TCGA-DX-A1KW-01A` amplicon 2: 15 genes in
  `inputs/TCGA-DX-A1KW-01A/gene_subset.txt`

The TCGA-DX-A1KW-01A amplicon 2 graph and cycles inputs were read directly
from the local `CoRAL graph file examples/TCGA-DX-A1KW-01A/..._AA_results`
folder; no download was performed.

These plots use the current uncommitted `coral/plot_amplicons.py` directly.
The older `run_thin_arc_plot.py` wrapper is retained only as provenance for
the previous render and was not used for this run.

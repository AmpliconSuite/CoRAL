# TCGA-DX-A1KW-01A amplicon 2 plotting matrix

All cases use the hg19 AA graph and cycle files from the local
`CoRAL graph file examples/TCGA-DX-A1KW-01A` directory. Both PNG and PDF
outputs were generated.

| Case | Options exercised | Result |
|---|---|---|
| `01_default` | New defaults: 12-inch graph width, content-aware cycle height, total offset `0.10`, minimum coordinate width `0`, robust coverage, colorblind-safe cycle palette, 300 DPI. `--num-cycles 3` was used because this file contains 83 reconstructions. | Completes. As expected, plotting every gene and coordinate is heavily crowded for this amplicon. |
| `02_typical` | `--width 12 --aspect-ratio 0.5 --offset 0.08 --min-coord-width 0.08 --coverage-scale robust --cycle-list c1,c2,c3 --cycle-color-file ... --gene-fontsize 9 --dpi 300 --combined` | Clean publication-oriented graph/cycle/combined output. |
| `03_edge_compact` | `--width 5 --aspect-ratio 1 --offset 0 --min-coord-width 1 --coverage-scale full --cycle-list c10,c32 --dpi 72 --combined` | Completes on a deliberately undersized square canvas. Full scaling expands coverage to about 1,750; endpoint coordinates are suppressed while the chromosome label remains. Long titles and the complete gene set are intentionally too dense at this size. |
| `04_edge_wide_dense` | `--width 18 --aspect-ratio 0.3 --offset 0.5 --min-coord-width 0 --coverage-scale robust --cycle-list p13,c31,p34 --font-size 1.5 --dpi 150 --combined` | Completes with large total gaps, every coordinate, large fonts, complex paths/cycle, and custom colors. This is a deliberate density stress test, not a recommended layout. |
| `05_edge_no_text` | `--width 8 --aspect-ratio 0.4 --offset 0.2 --min-coord-width 0.2 --coverage-scale robust --cycle-list c1 --font-size 0 --dpi 100 --combined` | Completes with text, ticks, and axes hidden while graph/cycle geometry and the gene track remain. |

Validation was also checked with `--cycle-list c13`. CoRAL exits with status 2
and the concise message `Invalid value: ID 13 is a path, not a cycle`; it no
longer prints an internal traceback.

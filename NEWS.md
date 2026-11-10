# ploure (development version)

* Tighter legends: `ThemePublication()`, `ThemeEmbedding()` and `LegendGrob()` size the legend keys to the legend
  text (1.4 x its size, 0.25 x between rows: about 15 pt per row at 9 pt text instead of ggplot2's 23 pt). Larger
  glyphs still widen their keys. `ZoomEmbedding()` leaves room between the corner axes and the legend titles, and
  aligns the two legends at the top.
* Continuous color bars get a black border (`legend.frame`), matching the axis lines: in `ThemePublication()`,
  `ThemeEmbedding()`, `CohortMatrix()` and `LegendGrob()` (so also `ZoomEmbedding()`), and so in every function's
  color bar under these themes.
* Differential abundance (Milo): `DABeeswarmPlot()` (neighbourhood fold changes per group, quasirandom points under
  an empty boxplot, significant neighbourhoods colored, groups ordered by mean fold change, the means on a right-hand
  axis) and `DAEmbeddingPlot()` (the neighbourhood graph on an embedding, sized by neighbourhood size, with corner
  axes and optional group labels). Both take the DA results as a data frame with any fold change and significance
  columns; the embedding plot takes the neighbourhood positions from a Milo object (`miloR`, suggested) or a data
  frame. Generalized from PLEURA's `DABoxPlot()` and `DAReductionPlot()`, which they reproduce
  (`validation/pleura_milo.R`). Example data: `sim_nhoods`.
* Text sizes follow the theme. The text drawn inside plots by `EulerPlot()`, `VolcanoPlot()`, `DECountPlot()`,
  `StageCountPlot()` and `QCFlowChart()` (counts, labels, set names, box text) is sized from the current theme's tick
  labels with the new `ThemeTextSize()`, so it stays at or below the tick-label size and scales with
  `ThemePublication(base_size = )`; at the default size it matches the previous look. The theme sizes these
  functions, `CohortMatrix()` and `LegendGrob()` set themselves are now relative to the theme. `EulerPlot()` keeps the
  current theme (minus axes) instead of `theme_void()`, so its title is the theme's title size, above the set names.
  `ThemePublication()` and `ThemeEmbedding()` scale the legend key size with `base_size` (unchanged at 11). The
  ported embedding functions keep their fixed sizes.
* Figure font: `FigureFont()`, `SetFigureFont()` and `FigurePDF()` use Arial (embedded through `cairo_pdf()`) when a
  licensed copy is installed, and R's Helvetica through `pdf()` otherwise, with a warning; ggrepel labels included.
  `CheckFigureFonts()` lists a PDF's fonts with `pdffonts` and flags anything else. Ported from PLEURA.
* `ThemePublication()` places strips beside an axis outside it (`strip.placement = 'outside'`), so panels stacked
  with patchwork keep their axes on the panels. `DistributionStack()` gains `color_title`. README: assembling
  figures with `patchwork::free()`.
* `VolcanoPlot()` and `DECountPlot()`: a volcano plot and diverging bars of the significant genes per group, sharing
  one classification (significant or not, split by direction and by |log2(FC)| >= `fc_cutoff`) and its colors. They
  take any results table by column name (`fc_col`, `p_col`, `gene_col`; e.g. DESeq2 `log2FoldChange` / `padj`,
  FindMarkers `avg_log2FC` / `p_val_adj`, edgeR `logFC` / `FDR`), with the genes in a column or as rownames, and
  custom axis titles for other effect sizes. The volcano labels given genes, bold
  genes and the top significant genes. Generalized from PLEURA's chemistry supplementary figure; the volcano is
  reproduced exactly, and the count bars have the same bars with the counts placed clear of the zero line
  (`validation/pleura_de_plots.R`). Example data: `sim_de`.
* `EulerPlot()`: a two-set Euler diagram in plain ggplot (circle areas proportional to the set sizes, overlap to the
  shared elements), with the overlap statistics below: Jaccard index, overlap coefficient and, given a universe,
  Fisher's exact test. `SetOverlapStats()` returns the same statistics as a table. Generalized from PLEURA's v1/v2
  barcode comparison, which it reproduces (`validation/pleura_euler.R`).

# ploure 0.2.0

* Palettes from other packages: `PaletteColors()`, `ColorMap()`, `PaletteFunction()` and `ListPalettes()` take
  palettes from MetBrewer (`'met:Hiroshige'`), NatParksPalettes (`'parks:Yellowstone'`), RColorBrewer, ggsci
  (journal palettes, `'ggsci:npg'`, and `'ggsci:igv'` with 51 colors), Polychrome (many distinguishable colors,
  `'polychrome:glasbey'`, near-white colors dropped), scico (scientific color maps, `'scico:batlow'`), viridisLite
  and base R HCL palettes. Every colors argument (`color_df`, `full_colors`, `colors`, `colfunc`, cohort-matrix
  blocks, `GroupedPalette()` ramps) also takes a palette name.
* `QCFlowChart()`: a left-to-right QC flow chart from a table of steps (a shared trunk that fans out into one row per
  branch) and a table of removals by category. Generalized from PLEURA's v2 QC flow chart, which it reproduces
  (`validation/pleura_qc_flow_chart.R`). Example data: `sim_qc_steps`, `sim_qc_removals`.
* `CohortMatrix()`: one column per sample, grouped into facets, with an optional bar plot on top and blocks of square
  tile rows below, each with its own color scale and legend. Generalized from PLEURA's v2 cohort matrix, which it
  reproduces (`validation/pleura_cohort_matrix.R`). Example data: `sim_samples`.
* QC panels: `StageCountPlot()` (cells kept at each filtering stage, in total stacked by a sample group or per
  sample) and `MetricDensityPlot()` (QC metric densities per library with the QC limits, a transform per metric).
  `DistributionStack()` gains `show_legend`. Generalized from PLEURA's v2 QC summary; the four-panel figure itself is
  project-specific and is assembled from these panels in `validation/pleura_qc_summary.R`.

# ploure 0.1.0

First release.

* Data core: `FetchCellData()` for Seurat objects and plain data frames.
* Embeddings: `PlotEmbedding()` and `FeatureEmbedding()` (ported from PLEURA / COMPACT, arguments and defaults
  unchanged), `ZoomEmbedding()`; layers `EmbeddingAxes()`, `OutlineCells()` / `OutlinePolygon()`, `GroupLabels()`;
  legends `LegendGrob()`, `ExtractLegend()`.
* Themes and palettes: `ThemePublication()`, `ThemeEmbedding()`, `ReadColors()` / `WriteColors()`,
  `GroupedPalette()`, `CheckPaletteCVD()`.
* Distributions: `QuasiBoxPlot()`, `DistributionStack()`. Marker dot plot: `DotPlotData()`, `MarkerDotPlot()`.
* Validated against the PLEURA figures: all 91 compared plot layers are identical (`validation/`).

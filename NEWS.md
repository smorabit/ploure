# ploure (development version)

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

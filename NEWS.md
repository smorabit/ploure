# ploure (development version)

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

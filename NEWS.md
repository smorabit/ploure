# ploure (development version)

* `QCFlowChart()`: a left-to-right QC flow chart from a table of steps (a shared trunk that fans out into one row per
  branch) and a table of removals by category. Generalized from PLEURA's v2 QC flow chart, which it reproduces
  (`validation/pleura_qc_flow_chart.R`). Example data: `sim_qc_steps`, `sim_qc_removals`.

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

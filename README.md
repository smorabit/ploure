# ploure

**PLOtting Utilities in R**, and Catalan for "to rain". A lightweight R package of plotting functions shared across
single-cell projects.

Private for now; it may become public when the PLEURA paper is published.

## Scope

A plot belongs here only if:

- it takes standard inputs: a Seurat object, or a data frame with cell / sample / group columns;
- it would be useful in another single-cell project unchanged.

Project-specific figures stay in their project repository.

The package is R only (no matplotlib counterpart). Module plots stay in hdWGCNA for now.

## Design: data-frame core, thin object adapters

`FetchCellData()` turns any supported input into one tidy per-cell data frame (`cell`, `x`, `y`, and the requested
columns), and every plotting function draws from it:

- **Seurat object:** coordinates from a reduction (or two metadata columns); values from the metadata or, for genes,
  an assay layer (default `data`).
- **Data frame:** a `cell` column (or rownames), two coordinate columns (`reduction = c('umap_1', 'umap_2')`, or
  `reduction = 'umap'` to find `umap_1` / `umap_2`).
- Later: SingleCellExperiment, and tables exported from h5ad.

Cached embedding tables can be plotted directly, without wrapping them in a Seurat object.

## Functions (v0.1)

| Module | Functions |
|---|---|
| Data | `FetchCellData()`, `sim_cells` (simulated example data) |
| Embeddings | `PlotEmbedding()`, `FeatureEmbedding()`, `ZoomEmbedding()` |
| Embedding layers | `EmbeddingAxes()`, `EmbeddingAxisLabels()`, `OutlineCells()`, `OutlinePolygon()`, `GroupLabels()` |
| Legends | `LegendGrob()`, `ExtractLegend()` |
| Themes | `ThemePublication()`, `ThemeEmbedding()` |
| Palettes | `ReadColors()`, `WriteColors()`, `GroupedPalette()`, `CheckPaletteCVD()` |
| Distributions | `QuasiBoxPlot()`, `DistributionStack()` |
| Dot plot | `DotPlotData()`, `MarkerDotPlot()` |
| Overview (v0.2) | `QCFlowChart()`, `CohortMatrix()` |

`PlotEmbedding()` and `FeatureEmbedding()` keep the arguments and defaults of the PLEURA / COMPACT versions, so
existing figures are reproduced.

Conventions for continuous colors (`FeatureEmbedding()`, `ZoomEmbedding()`): limits are percentiles (`'q99'`) or
numbers; values above the upper limit are clipped; cells at or below the lower limit are drawn in light grey
underneath; high values are drawn on top; the default palette is viridis inferno, reversed.

`ThemePublication()` is a complete theme: set it once with `theme_set(ThemePublication())` rather than adding it to
each plot, so the plotting functions' own theme settings (italic gene names, strip placement) are kept.

A worked example that draws the four kinds of zoom figure from simulated data: `inst/examples/zoom_figures.R`.

## Validation against PLEURA

`validation/pleura_comparison.R` redraws PLEURA figures with the original utilities and with ploure, in PLEURA's
own R environment (ploure installed into `validation/lib`), and compares the built plot data layer by layer. Its
figures and library are local only (git-ignored). At v0.1.0 all 91 compared layers are identical: PlotEmbedding,
FeatureEmbedding (table and Seurat inputs), the four ZoomEmbedding figures, the marker dot plot, the boxes and count
bars of the per-cluster QC stack, and the grouped patient palette. The QC stack's point subsample differs by design
(fixed per-group sampling).

## Install into a project

Pin a tag, so published figures keep regenerating when defaults change later.

From a local clone (works with pixi, which does not install R packages from GitHub itself). In the project's
`pixi.toml`:

```toml
[tasks]
install-ploure = "git -C /path/to/ploure checkout v0.1.0 && R CMD INSTALL /path/to/ploure"
```

Or from GitHub with a token (private repository):

```r
remotes::install_github('smorabit/ploure', ref = 'v0.1.0', auth_token = Sys.getenv('GITHUB_PAT'),
                        upgrade = 'never', dependencies = FALSE)
```

## Development

The development environment is in `pixi.toml` (R 4.5, the package dependencies, roxygen2, testthat, vdiffr).

```bash
pixi run document   # roxygen: man/ and NAMESPACE
pixi run test       # testthat + vdiffr snapshots
pixi run check      # R CMD check
pixi run install    # install into the pixi environment
```

Interactively: `devtools::load_all()`.

- Code style follows Sam's `STYLE.md`: 4-space indentation, `<-`, single quotes, one verb per line in ggplot chains,
  PascalCase for exported functions, a leading `.` for internal helpers.
- Tests and examples use simulated data only (`data-raw/sim_cells.R`); never patient data.
- Snapshot tests draw points without rasterization. After an intended visual change, review the new snapshots with
  `testthat::snapshot_review()` and accept them.

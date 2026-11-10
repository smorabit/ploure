# visual snapshots of each exported plot (vdiffr); points are not rasterized, so the SVGs stay stable and readable

skip_if_not_installed('vdiffr')
skip_if_not_installed('ggrepel')
small <- sim_small(400)
# the publication theme is set globally (as in a notebook) and restored at the end of the file
old_theme <- ggplot2::theme_set(ThemePublication())

test_that('PlotEmbedding', {
    set.seed(1)
    p <- PlotEmbedding(small, group.by = 'cluster', reduction = 'umap', raster = FALSE, point_size = 0.5,
                       corner_axes = TRUE, plot_theme = ThemeEmbedding())
    vdiffr::expect_doppelganger('plot-embedding', p)
    set.seed(1)
    p <- PlotEmbedding(small, group.by = 'cluster', reduction = 'umap', raster = FALSE, point_size = 0.5,
                       split.by = 'lineage', plot_under = TRUE, label = FALSE, selected = c('B', 'NK'))
    vdiffr::expect_doppelganger('plot-embedding-split', p$Lymphoid)
})

test_that('FeatureEmbedding', {
    p <- FeatureEmbedding(small, features = 'score', reduction = 'umap', raster = FALSE, plot_min = 'q20',
                          plot_max = 'q95', corner_axes = TRUE)
    vdiffr::expect_doppelganger('feature-embedding', p)
    p <- FeatureEmbedding(small, features = c('CD3E', 'LYZ'), reduction = 'umap', raster = FALSE, ncol = 2,
                          same_range = TRUE)
    vdiffr::expect_doppelganger('feature-embedding-two', p)
})

test_that('ZoomEmbedding', {
    skip_if_not_installed('MASS')
    sub <- small[small$lineage == 'Lymphoid', ]
    p <- ZoomEmbedding(
        full = small, subset = sub, full_reduction = 'umap', subset_reduction = c('sub_umap_1', 'sub_umap_2'),
        full_color_by = 'cluster', subset_color_by = 'patient', subset_label_by = 'patient',
        subset_colors = GroupedPalette(small$patient, small$disease),
        legend_colors = list(full = NULL, subset = attr(GroupedPalette(small$patient, small$disease), 'families')),
        label_map = c(P01 = 'P01 LC'), raster = FALSE, titles = c('All cells', 'Lymphoid')
    )
    vdiffr::expect_doppelganger('zoom-embedding-categorical', p)
    p <- ZoomEmbedding(
        full = small, subset = sub, full_reduction = 'umap', subset_reduction = c('sub_umap_1', 'sub_umap_2'),
        full_color_by = 'CD3E', subset_color_by = 'score', plot_min = c(0, 'q5'), plot_max = 'q99',
        outline = FALSE, raster = FALSE
    )
    vdiffr::expect_doppelganger('zoom-embedding-continuous', p)
})

test_that('outline and labels on any embedding', {
    skip_if_not_installed('MASS')
    p <- ggplot2::ggplot(small, ggplot2::aes(x = umap_1, y = umap_2)) +
        ggplot2::geom_point(size = 0.3, color = 'grey50') +
        OutlineCells(small$umap_1, small$umap_2, small$lineage == 'Myeloid') +
        GroupLabels(small$umap_1, small$umap_2, small$lineage) +
        EmbeddingAxes(small$umap_1, small$umap_2) +
        ggplot2::coord_equal() +
        ThemeEmbedding()
    vdiffr::expect_doppelganger('outline-labels', p)
})

test_that('distributions', {
    skip_if_not_installed('ggbeeswarm')
    p <- QuasiBoxPlot(small, x = 'sample', y = 'n_genes', color_by = 'disease', facet_by = 'disease',
                      raster = FALSE)
    vdiffr::expect_doppelganger('quasi-boxplot', p)
    p <- DistributionStack(small, group_by = 'cluster', metrics = c('UMIs' = 'n_umi', 'Mito (%)' = 'pct_mt'),
                           color_by = 'lineage', facet_by = 'lineage', raster = FALSE, x_label = 'Cluster')
    vdiffr::expect_doppelganger('distribution-stack', p)
})

test_that('marker dot plot', {
    dots <- DotPlotData(small, features = c('CD3E', 'CD8A', 'MS4A1', 'LYZ', 'CD14'), group_by = 'cluster')
    p <- MarkerDotPlot(dots, gene_groups = list(T = c('CD3E', 'CD8A'), B = 'MS4A1', Myeloid = c('LYZ', 'CD14')),
                       x_groups = with(unique(small[c('cluster', 'lineage')]), setNames(lineage, cluster)),
                       x_label = 'Cluster')
    vdiffr::expect_doppelganger('marker-dot-plot', p)
})

test_that('palette check and legend', {
    skip_if_not_installed('colorspace')
    vdiffr::expect_doppelganger('palette-cvd', CheckPaletteCVD(GroupedPalette(small$patient, small$disease)))
    leg <- LegendGrob(colors = c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c'), title = 'Lineage', ncol = 2)
    vdiffr::expect_doppelganger('legend-grob', patchwork::wrap_elements(leg))
})


test_that('QC flow chart', {
    p <- QCFlowChart(sim_qc_steps, sim_qc_removals,
                     branch_colors = c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c', Stromal = '#33a02c'),
                     category_labels = c('Damaged / low quality' = 'Damaged / LQ'), titles = c('All cells', 'Per lineage'))
    vdiffr::expect_doppelganger('qc-flow-chart', p)
})

test_that('cohort matrix', {
    cells <- as.data.frame(table(sample = small$sample, lineage = small$lineage), responseName = 'n')
    p <- CohortMatrix(
        sim_samples, group_col = 'disease', bar_data = cells, bar_fill = 'lineage', bar_title = 'Lineage',
        blocks = list(
            list(vars = c('Chemistry' = 'chemistry', 'Site' = 'site'), title = 'Chemistry / site',
                 colors = c("5' v2" = '#1b9e77', 'GEM-X' = '#d95f02', 'Site A' = '#7570b3', 'Site B' = '#e7298a')),
            list(vars = c('Loaded CD45+' = 'pct_cd45', 'Loaded EPCAM+' = 'pct_epcam'), title = 'Cells loaded (%)',
                 limits = c(0, 100)),
            list(vars = c('Sex' = 'sex', 'Age' = 'age_group'), title = 'Sex / age')
        )
    )
    vdiffr::expect_doppelganger('cohort-matrix', p)
})

test_that('stage counts', {
    stages <- c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final')
    vdiffr::expect_doppelganger('stage-count-total', StageCountPlot(sim_samples, stages, fill_by = 'disease'))
    vdiffr::expect_doppelganger('stage-count-sample',
                                StageCountPlot(sim_samples, stages, by_sample = TRUE, facet_by = 'disease'))
})

test_that('metric densities', {
    p <- MetricDensityPlot(small, metrics = c('UMIs per cell' = 'n_umi', 'Mitochondrial (%)' = 'pct_mt'),
                           group_by = 'sample', color_by = 'disease', thresholds = list(pct_mt = 15),
                           trans = c(n_umi = 'log10', pct_mt = 'sqrt'))
    vdiffr::expect_doppelganger('metric-density', p)
})

test_that('Euler diagram', {
    p <- EulerPlot(list(Atlas = paste0('g', 1:300), Ours = paste0('g', 201:450)), universe = 2000,
                   title = 'Marker genes')
    vdiffr::expect_doppelganger('euler', p)
})

test_that('differential expression', {
    de <- sim_de[sim_de$lineage == 'Myeloid', ]
    vdiffr::expect_doppelganger('volcano', VolcanoPlot(de, label = c('TECH1', 'TECH3'), label_bold = 'MARKD',
                                                       label_top = 3, raster = FALSE))
    vdiffr::expect_doppelganger('de-count', DECountPlot(sim_de, group_by = 'lineage'))
})

test_that('differential abundance', {
    set.seed(1)
    vdiffr::expect_doppelganger('da-beeswarm', DABeeswarmPlot(sim_nhoods, color_limit = 3, raster = FALSE))
    vdiffr::expect_doppelganger('da-embedding', DAEmbeddingPlot(sim_nhoods, x_col = 'umap_1', y_col = 'umap_2',
                                                                color_limit = 3, raster = FALSE))
})

ggplot2::theme_set(old_theme)

test_that('DABeeswarmPlot orders groups by mean fold change and colors only significant neighbourhoods', {
    skip_if_not_installed('ggbeeswarm')
    p <- DABeeswarmPlot(sim_nhoods, raster = FALSE)
    expect_s3_class(p, 'ggplot')
    expect_false('Mixed' %in% p$data$group)
    means <- tapply(p$data$fc, p$data$group, mean)
    expect_equal(names(means), levels(p$data$group))
    expect_true(all(diff(means) <= 0))
    expect_equal(!is.na(p$data$fc_color), p$data$SpatialFDR <= 0.05)
    # mean labels: the right-hand axis, one per group in the group order
    b <- ggplot2::ggplot_build(p)
    right <- b$layout$panel_params[[1]]$y.sec$get_labels()
    expect_equal(as.numeric(right), as.numeric(means), tolerance = 0.01)
    p0 <- DABeeswarmPlot(sim_nhoods, show_mean = FALSE, raster = FALSE)
    expect_null(p0$scales$get_scales('y')$secondary.axis)
})

test_that('DABeeswarmPlot keeps the factor order, takes other columns, and checks its input', {
    skip_if_not_installed('ggbeeswarm')
    p <- DABeeswarmPlot(sim_nhoods, order_groups = FALSE, exclude_groups = NULL, raster = FALSE)
    expect_equal(levels(p$data$group), levels(sim_nhoods$cluster))
    de <- sim_nhoods
    names(de)[names(de) == 'logFC'] <- 'effect'
    de$label <- de$cluster
    p <- DABeeswarmPlot(de, group_col = 'label', fc_col = 'effect', signif_col = 'PValue', signif_thresh = 0.01,
                        raster = FALSE)
    expect_equal(!is.na(p$data$fc_color), p$data$PValue <= 0.01)
    expect_error(DABeeswarmPlot(sim_nhoods, fc_col = 'log2FoldChange'), 'not found.*The table has')
    expect_error(DABeeswarmPlot(sim_nhoods, colors = c('blue', 'red')), 'three colors')
})

test_that('DAEmbeddingPlot draws significant neighbourhoods on top, clips outliers and adds corner axes', {
    p <- DAEmbeddingPlot(sim_nhoods, x_col = 'umap_1', y_col = 'umap_2', raster = FALSE)
    expect_s3_class(p, 'ggplot')
    ns <- p$layers[[1]]$data
    sig <- p$layers[[2]]$data
    expect_true(all(ns$SpatialFDR > 0.05) && all(sig$SpatialFDR <= 0.05))
    expect_equal(nrow(ns) + nrow(sig), nrow(sim_nhoods))
    expect_false(is.unsorted(abs(sig$fc)))
    lims <- stats::quantile(sim_nhoods$logFC, c(0.001, 0.999))
    expect_true(all(p$data$fc >= lims[1] - 1e-9 & p$data$fc <= lims[2] + 1e-9))
    # corner axes, labeled from the coordinate columns
    texts <- Filter(function(l) inherits(l$geom, 'GeomText'), p$layers)
    expect_equal(unname(vapply(texts, function(l) l$aes_params$label, character(1))), c('UMAP 1', 'UMAP 2'))
    expect_length(DAEmbeddingPlot(sim_nhoods, x_col = 'umap_1', y_col = 'umap_2', corner_axes = FALSE)$layers, 2)
})

test_that('DAEmbeddingPlot joins separate DA results and labels groups from the cells', {
    skip_if_not_installed('ggrepel')
    nhoods <- sim_nhoods[, c('Nhood', 'umap_1', 'umap_2', 'nhood_size')]
    da <- sim_nhoods[sample(nrow(sim_nhoods)), c('Nhood', 'logFC', 'SpatialFDR', 'cluster')]
    p <- DAEmbeddingPlot(nhoods, da, x_col = 'umap_1', y_col = 'umap_2', label_by = 'cluster', cells = sim_cells,
                         raster = FALSE)
    m <- merge(p$data[, c('Nhood', 'x', 'logFC')], sim_nhoods[, c('Nhood', 'umap_1', 'logFC')], by = 'Nhood')
    expect_equal(m$x, m$umap_1)
    expect_equal(m$logFC.x, m$logFC.y)
    repel <- Filter(function(l) inherits(l$geom, 'GeomTextRepel'), p$layers)
    expect_length(repel, 1)
    expect_setequal(repel[[1]]$data$text, levels(sim_cells$cluster))
    expect_error(DAEmbeddingPlot(nhoods, da, x_col = 'umap_1', y_col = 'umap_2', label_by = 'cluster'), 'cells')
    expect_error(DAEmbeddingPlot(nhoods, x_col = 'umap_1', y_col = 'umap_2'), 'logFC')
})

test_that('DAEmbeddingPlot reads a Milo object', {
    skip_if_not_installed('miloR')
    skip_if_not_installed('SingleCellExperiment')
    set.seed(1)
    cells <- sim_cells[sample(nrow(sim_cells), 600), ]
    counts <- matrix(stats::rpois(600 * 50, 2), nrow = 50, dimnames = list(paste0('g', 1:50), cells$cell))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts, logcounts = log1p(counts)))
    sce$cluster <- cells$cluster
    emb <- as.matrix(cells[, c('umap_1', 'umap_2')])
    rownames(emb) <- cells$cell
    SingleCellExperiment::reducedDim(sce, 'PCA') <- emb
    SingleCellExperiment::reducedDim(sce, 'UMAP') <- emb
    milo <- miloR::Milo(sce)
    milo <- miloR::buildGraph(milo, k = 15, d = 2, reduced.dim = 'PCA')
    milo <- miloR::makeNhoods(milo, prop = 0.1, k = 15, d = 2, reduced_dims = 'PCA', refined = TRUE)
    milo <- miloR::buildNhoodGraph(milo)
    n <- ncol(miloR::nhoods(milo))
    da <- data.frame(Nhood = seq_len(n), logFC = stats::rnorm(n), SpatialFDR = stats::runif(n))
    expect_error(DAEmbeddingPlot(milo, da, reduction = 'TSNE'), 'not found')
    p <- DAEmbeddingPlot(milo, da, reduction = 'UMAP', label_by = 'cluster', raster = FALSE)
    expect_equal(nrow(p$data), n)
    index <- as.numeric(igraph::vertex_attr(miloR::nhoodGraph(milo))$name)
    expect_equal(p$data$x[order(p$data$Nhood)], unname(cells$umap_1[index]))
    texts <- Filter(function(l) inherits(l$geom, 'GeomText'), p$layers)
    expect_equal(unname(vapply(texts, function(l) l$aes_params$label, character(1))), c('UMAP 1', 'UMAP 2'))
})

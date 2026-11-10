#---------------------------------------------------------#
# Per-group distributions: quasirandom points under empty boxplots
#
# The boxes summarize every cell; the points are a fixed random subsample per group (rasterized), so large groups do
# not swamp the plot and the figure is reproducible.
#---------------------------------------------------------#

#' Quasirandom points under an empty boxplot
#'
#' One metric per group (e.g. per sample or cluster): an empty-fill boxplot of all cells, over quasirandom points of
#' a fixed random subsample of up to `n_points` cells per group.
#'
#' @param object A data frame with one row per cell, or a Seurat object (columns from its metadata).
#' @param x Grouping column (x axis).
#' @param y Metric column.
#' @param color_by Optional column to color the points by (default: all grey).
#' @param colors Named colors, or a palette name ([PaletteColors()]), for `color_by`.
#' @param facet_by Optional column to split the groups into column facets (e.g. disease), with free x space.
#' @param n_points Points drawn per group (and facet).
#' @param point_size,point_alpha,point_width Quasirandom point size, transparency and spread.
#' @param box_width,box_linewidth Boxplot width and line width.
#' @param raster,raster_dpi Rasterize the points, and at which resolution.
#' @param seed Seed for the point subsample.
#'
#' @return A ggplot.
#' @export
#' @examples
#' QuasiBoxPlot(sim_cells, x = 'sample', y = 'n_genes', color_by = 'disease', facet_by = 'disease')
QuasiBoxPlot <- function(
    object,
    x,
    y,
    color_by = NULL,
    colors = NULL,
    facet_by = NULL,
    n_points = 300,
    point_size = 0.12,
    point_alpha = 0.5,
    point_width = 0.35,
    box_width = 0.65,
    box_linewidth = 0.25,
    raster = TRUE,
    raster_dpi = 600,
    seed = 12345
){
    df <- FetchCellData(object, vars = c(x, y, color_by, facet_by))
    df$group <- .AsFactor(df[[x]])
    df$value <- df[[y]]
    df$facet <- if(is.null(facet_by)) factor('all') else .AsFactor(df[[facet_by]])
    df$color <- if(is.null(color_by)) factor('all') else .AsFactor(df[[color_by]])
    points_df <- .SamplePerGroup(df, c('group', 'facet'), n_points, seed)

    p <- ggplot(df, aes(x = group, y = value)) +
        .QuasiLayer(points_df, color_by, point_size, point_alpha, point_width, raster, raster_dpi) +
        geom_boxplot(outlier.shape = NA, fill = NA, width = box_width, linewidth = box_linewidth) +
        scale_y_continuous(labels = scales::label_comma()) +
        labs(x = NULL, y = y, color = color_by)
    p <- p + .PointColorScale(color_by, colors, levels(df$color))
    if(!is.null(facet_by)) p <- p + facet_grid(~ facet, scales = 'free_x', space = 'free_x')
    p
}

#' Stacked per-group distributions of several metrics
#'
#' Several metrics per group (e.g. UMIs, genes and % mitochondrial per cluster), one row per metric with its own y
#' scale and shared group columns, optionally with a row of cell counts per group on top. Each row is a quasirandom
#' point subsample under an empty boxplot of all cells, as in [QuasiBoxPlot()].
#'
#' @inheritParams QuasiBoxPlot
#' @param group_by Grouping column (x axis).
#' @param metrics Metric columns; names, if given, are used as row labels (e.g. `c('Genes' = 'n_genes')`).
#' @param count_bar Add a row with the number of cells per group (filled by `color_by`, if given).
#' @param count_label Row label of the count bar.
#' @param x_label X axis title.
#' @param show_legend Show the legend of `color_by`.
#' @param color_title Legend title (default `color_by`).
#'
#' @return A ggplot.
#' @export
#' @examples
#' DistributionStack(sim_cells, group_by = 'cluster',
#'                   metrics = c('UMIs' = 'n_umi', 'Genes' = 'n_genes', 'Mito (%)' = 'pct_mt'),
#'                   color_by = 'lineage', x_label = 'Cluster')
DistributionStack <- function(
    object,
    group_by,
    metrics,
    color_by = NULL,
    colors = NULL,
    facet_by = NULL,
    count_bar = TRUE,
    count_label = 'Cells',
    n_points = 300,
    point_size = 0.1,
    point_alpha = 0.5,
    point_width = 0.35,
    box_width = 0.7,
    box_linewidth = 0.2,
    raster = TRUE,
    raster_dpi = 600,
    x_label = group_by,
    seed = 12345,
    show_legend = FALSE,
    color_title = color_by
){
    if(is.null(names(metrics))) names(metrics) <- metrics
    names(metrics)[names(metrics) == ''] <- metrics[names(metrics) == '']
    row_levels <- c(if(count_bar) count_label, names(metrics))

    df <- FetchCellData(object, vars = c(group_by, unname(metrics), color_by, facet_by))
    df$group <- .AsFactor(df[[group_by]])
    df$facet <- if(is.null(facet_by)) factor('all') else .AsFactor(df[[facet_by]])
    df$color <- if(is.null(color_by)) factor('all') else .AsFactor(df[[color_by]])

    # long table: one row per cell and metric
    long <- do.call(rbind, lapply(names(metrics), function(m){
        data.frame(group = df$group, facet = df$facet, color = df$color, metric = m, value = df[[metrics[[m]]]])
    }))
    long$metric <- factor(long$metric, levels = row_levels)
    points_df <- .SamplePerGroup(long, c('group', 'facet', 'metric'), n_points, seed)

    p <- ggplot(mapping = aes(x = group, y = value))
    if(count_bar){
        counts <- stats::aggregate(list(value = rep(1, nrow(df))),
                                   by = list(group = df$group, facet = df$facet, color = df$color), FUN = sum)
        counts$metric <- factor(count_label, levels = row_levels)
        p <- p + if(is.null(color_by)){
            geom_col(data = counts, width = 0.8, fill = 'grey60')
        } else {
            geom_col(data = counts, aes(fill = color), width = 0.8)
        }
    }
    p <- p +
        .QuasiLayer(points_df, color_by, point_size, point_alpha, point_width, raster, raster_dpi) +
        geom_boxplot(data = long, outlier.shape = NA, fill = NA, width = box_width, linewidth = box_linewidth) +
        facet_grid(if(is.null(facet_by)) metric ~ . else metric ~ facet, scales = 'free', space = 'free_x',
                   switch = 'y') +
        scale_y_continuous(labels = scales::label_number(big.mark = ',', accuracy = NULL), n.breaks = 4) +
        labs(x = x_label, y = NULL) +
        theme(
            strip.background = element_blank(),
            strip.placement = 'outside',
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            panel.spacing.y = grid::unit(0.4, 'lines')
        )
    if(!is.null(color_by)){
        colors <- .AsNamedColors(colors, levels(df$color))
        p <- p +
            scale_color_manual(values = colors, aesthetics = c('color', 'fill'), name = color_title) +
            guides(color = if(show_legend) guide_legend(override.aes = list(size = 2, alpha = 1)) else 'none',
                   fill = 'none')
    }
    p
}

# up to n rows per combination of the `by` columns, reproducibly
.SamplePerGroup <- function(df, by, n, seed){
    old <- if(exists('.Random.seed', envir = globalenv())) get('.Random.seed', envir = globalenv()) else NULL
    on.exit(if(!is.null(old)) assign('.Random.seed', old, envir = globalenv()))
    set.seed(seed)
    key <- interaction(df[by], drop = TRUE)
    rows <- unlist(lapply(split(seq_len(nrow(df)), key), function(i) if(length(i) > n) sample(i, n) else i))
    df[sort(rows), , drop = FALSE]
}

.QuasiLayer <- function(points_df, color_by, size, alpha, width, raster, raster_dpi){
    if(!requireNamespace('ggbeeswarm', quietly = TRUE)) stop('quasirandom points need the ggbeeswarm package.')
    layer <- if(is.null(color_by)){
        ggbeeswarm::geom_quasirandom(data = points_df, color = 'grey40', size = size, alpha = alpha, width = width)
    } else {
        ggbeeswarm::geom_quasirandom(data = points_df, aes(color = color), size = size, alpha = alpha,
                                     width = width)
    }
    if(raster) ggrastr::rasterise(layer, dpi = raster_dpi) else layer
}

.PointColorScale <- function(color_by, colors, levels){
    if(is.null(color_by)) return(NULL)
    scale_color_manual(values = .AsNamedColors(colors, levels),
                       guide = guide_legend(override.aes = list(size = 2, alpha = 1)))
}

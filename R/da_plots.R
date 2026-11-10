#---------------------------------------------------------#
# Differential abundance (Milo neighbourhoods)
#
# Generalized from PLEURA's utils/milo_functions.R (DABoxPlot, DAReductionPlot), as drawn in
# v2/composition/milo_figures.Rmd. Both take the neighbourhood DA results as a data frame (one row per neighbourhood,
# e.g. miloR::testNhoods() output with a label column from miloR::annotateNhoods()); the embedding plot also takes
# the neighbourhood positions, from a Milo object or a data frame.
#---------------------------------------------------------#

#' Neighbourhood fold changes per group (beeswarm)
#'
#' One row per group (e.g. cell type): the log fold change of every neighbourhood as quasirandom points under an empty
#' boxplot, the significant neighbourhoods colored by their fold change and the others grey, and the mean fold change
#' of each group as the labels of a right-hand axis. Groups are ordered by their mean fold change.
#'
#' @param da Neighbourhood DA results: a data frame with one row per neighbourhood (e.g. `miloR::testNhoods()` output
#'   with a group label column from `miloR::annotateNhoods()`).
#' @param group_col Group column (default `'cluster'`).
#' @param fc_col,signif_col Fold change and significance columns (default `'logFC'` and `'SpatialFDR'`; e.g.
#'   `'PValue'` for nominal p values).
#' @param signif_thresh Neighbourhoods with `signif_col <= signif_thresh` are colored.
#' @param exclude_groups Groups to leave out (default `'Mixed'`, the label `annotateNhoods()` cutoffs often give
#'   mixed neighbourhoods).
#' @param order_groups Order the groups by mean fold change (highest on top); otherwise keep the factor or table
#'   order.
#' @param colors Three colors: lower, not significant (and zero), higher.
#' @param color_limit Optional fold change at which the color scale saturates (symmetric, e.g. 3).
#' @param show_mean Label each group's mean fold change on the right.
#' @param point_size Point size.
#' @param raster,raster_dpi Rasterize the points, and at which resolution.
#' @param x_label X axis title.
#' @param title Plot title.
#'
#' @return A ggplot.
#' @export
#' @examples
#' DABeeswarmPlot(sim_nhoods, color_limit = 3)
DABeeswarmPlot <- function(
    da,
    group_col = 'cluster',
    fc_col = 'logFC',
    signif_col = 'SpatialFDR',
    signif_thresh = 0.05,
    exclude_groups = 'Mixed',
    order_groups = TRUE,
    colors = c('dodgerblue', 'lightgrey', 'darkgoldenrod1'),
    color_limit = NULL,
    show_mean = TRUE,
    point_size = 1,
    raster = TRUE,
    raster_dpi = 300,
    x_label = expression(log[2]*'(fold change)'),
    title = NULL
){
    if(!requireNamespace('ggbeeswarm', quietly = TRUE)) stop('DABeeswarmPlot needs the ggbeeswarm package.')
    df <- .DATable(da, c(group_col, fc_col, signif_col))
    colors <- .DAColors(colors)
    df$group <- as.character(df[[group_col]])
    df <- df[!(df$group %in% exclude_groups), ]
    df$fc <- df[[fc_col]]
    df$signif <- !is.na(df[[signif_col]]) & df[[signif_col]] <= signif_thresh

    # groups in factor (or first-seen) order, or by mean fold change
    levels <- if(is.factor(da[[group_col]])) intersect(levels(da[[group_col]]), df$group) else unique(df$group)
    means <- data.frame(group = levels, mean = vapply(levels, function(g) mean(df$fc[df$group == g]), numeric(1)))
    if(order_groups) means <- means[order(means$mean, decreasing = TRUE), ]
    means$group <- factor(means$group, levels = means$group)
    means$label <- trimws(format(means$mean, digits = 1))
    df$group <- factor(df$group, levels = levels(means$group))

    range <- max(abs(df$fc)) + 0.25
    df$fc_color <- ifelse(df$signif, df$fc, NA)
    # not significant first, so the colored points are drawn on top
    df <- rbind(df[!df$signif, ], df[df$signif, ])

    points <- ggbeeswarm::geom_quasirandom(aes(color = fc_color), method = 'pseudorandom', size = point_size,
                                           orientation = 'y')
    p <- ggplot(df, aes(x = fc, y = group)) +
        (if(raster) ggrastr::rasterise(points, dpi = raster_dpi) else points) +
        geom_boxplot(fill = NA, outlier.shape = NA) +
        scale_x_continuous(limits = c(-range, range)) +
        # the group means as a right-hand axis, clear of the points
        scale_y_discrete(sec.axis = if(show_mean) dup_axis(labels = means$label, name = NULL) else waiver()) +
        .DAColorScale(colors, colors[2], color_limit) +
        guides(color = 'none') +
        labs(x = x_label, y = NULL, title = title) +
        theme(
            axis.line = element_blank(),
            panel.border = element_rect(linewidth = 0.6, color = 'black', fill = NA),
            panel.grid.major.x = element_line(linewidth = 0.4, color = 'grey85'),
            panel.grid.major.y = element_line(linewidth = 0.25, color = 'grey85'),
            axis.ticks.y.right = element_blank()
        )
    p
}

#' Neighbourhood fold changes on an embedding
#'
#' The Milo neighbourhood graph on a 2D embedding: one point per neighbourhood at its index cell, sized by the
#' neighbourhood size, the significant neighbourhoods colored by their fold change (the most extreme on top) and the
#' others pale, with corner axes and optional group labels.
#'
#' @param nhoods The neighbourhood positions: a Milo object (with `miloR::buildNhoodGraph()` run; positions from
#'   `reduction`), or a data frame with a `Nhood` column, the coordinates (`x_col`, `y_col`) and the neighbourhood
#'   size (`size_col`), and the DA results too if `da` is NULL.
#' @param da Neighbourhood DA results with a `Nhood` column (e.g. `miloR::testNhoods()` output), joined to `nhoods`.
#' @param reduction For a Milo object: the reduced dimension to place the neighbourhoods on (e.g. `'UMAP'`).
#' @param x_col,y_col,size_col For a data frame: coordinate and neighbourhood size columns.
#' @inheritParams DABeeswarmPlot
#' @param colors Three colors: lower, not significant (and zero), higher.
#' @param adjust_outliers Clip the fold changes to this quantile range (e.g. 0.999: the 0.1% most extreme on each
#'   side), or NA to keep them.
#' @param point_sizes Point size range for the neighbourhood sizes.
#' @param label_by Optional cell-level column to label groups at their centroid: a `colData` column of a Milo object,
#'   or a column of `cells`.
#' @param cells For `label_by` with a data frame `nhoods`: the cells, as a data frame or Seurat object that
#'   [FetchCellData()] reads with `cells_reduction`.
#' @param cells_reduction Coordinates of `cells` (a reduction or two column names; default `c(x_col, y_col)`).
#' @param label_size Group label size in mm (default the current theme's tick-label size; see [ThemeTextSize()]).
#' @param corner_axes Replace the axes with short corner axes ([EmbeddingAxes()]).
#' @param axes_labels Corner axis labels (default from `reduction` or `x_col`, see [EmbeddingAxisLabels()]).
#' @param axes_frac Corner axis length, as a fraction of the data range.
#' @param plot_ratio Equal axis scales (`coord_equal()`).
#'
#' @return A ggplot.
#' @export
#' @examples
#' DAEmbeddingPlot(sim_nhoods, x_col = 'umap_1', y_col = 'umap_2', color_limit = 3)
#' DAEmbeddingPlot(sim_nhoods, x_col = 'umap_1', y_col = 'umap_2', label_by = 'cluster', cells = sim_cells)
DAEmbeddingPlot <- function(
    nhoods,
    da = NULL,
    reduction = 'UMAP',
    x_col = 'x',
    y_col = 'y',
    size_col = 'nhood_size',
    fc_col = 'logFC',
    signif_col = 'SpatialFDR',
    signif_thresh = 0.05,
    colors = c('dodgerblue', 'whitesmoke', 'darkgoldenrod1'),
    color_limit = NULL,
    adjust_outliers = 0.999,
    point_sizes = c(0.1, 4),
    label_by = NULL,
    cells = NULL,
    cells_reduction = c(x_col, y_col),
    label_size = ThemeTextSize(),
    corner_axes = TRUE,
    axes_labels = NULL,
    axes_frac = 0.15,
    plot_ratio = TRUE,
    raster = TRUE,
    raster_dpi = 300,
    title = NULL
){
    colors <- .DAColors(colors)
    is_milo <- inherits(nhoods, 'Milo')
    if(is_milo){
        layout <- .MiloNhoodLayout(nhoods, reduction)
        if(is.null(axes_labels)) axes_labels <- EmbeddingAxisLabels(reduction)
    } else {
        layout <- .DATable(nhoods, c('Nhood', x_col, y_col, size_col))
        layout$x <- layout[[x_col]]
        layout$y <- layout[[y_col]]
        layout$nhood_size <- layout[[size_col]]
        if(is.null(axes_labels)) axes_labels <- EmbeddingAxisLabels(x_col)
    }
    if(is.null(da)){
        df <- layout
    } else {
        da <- .DATable(da, 'Nhood')
        layout <- layout[, c('Nhood', 'x', 'y', 'nhood_size')]
        df <- merge(layout, da[, setdiff(colnames(da), c('x', 'y', 'nhood_size')), drop = FALSE], by = 'Nhood')
        df <- df[order(df$Nhood), ]
    }
    df <- .DATable(df, c(fc_col, signif_col))
    df$fc <- df[[fc_col]]
    if(!is.na(adjust_outliers)){
        lims <- stats::quantile(df$fc, c(1 - adjust_outliers, adjust_outliers), na.rm = TRUE)
        df$fc <- pmin(pmax(df$fc, lims[1]), lims[2])
    }
    df$signif <- !is.na(df[[signif_col]]) & df[[signif_col]] <= signif_thresh
    # by effect size, not significant first: the strongest significant neighbourhoods are drawn on top
    df <- df[order(abs(df$fc)), ]
    df <- rbind(df[!df$signif, ], df[df$signif, ])

    points <- function(layer) if(raster) ggrastr::rasterise(layer, dpi = raster_dpi) else layer
    p <- ggplot(df, aes(x = x, y = y, color = fc, size = nhood_size)) +
        points(geom_point(data = df[!df$signif, ], color = colors[2])) +
        points(geom_point(data = df[df$signif, ])) +
        .DAColorScale(colors, NULL, color_limit) +
        scale_size(range = point_sizes) +
        labs(color = expression(log[2]*'(FC)'), size = 'Neighborhood\nsize', title = title)

    if(!is.null(label_by)){
        cell_df <- .DALabelCells(nhoods, is_milo, reduction, label_by, cells, cells_reduction)
        p <- p + GroupLabels(cell_df$x, cell_df$y, cell_df[[label_by]], position = 'centroid', size = label_size)
    }
    if(plot_ratio) p <- p + coord_equal()
    p <- p + .NoAxesTheme()
    if(corner_axes) p <- p + EmbeddingAxes(df$x, df$y, frac = axes_frac, labels = axes_labels)
    p
}

# the table as a data frame, with the needed columns (the error lists the table's columns)
.DATable <- function(df, cols){
    df <- as.data.frame(df)
    missing <- setdiff(cols, colnames(df))
    if(length(missing) > 0){
        stop('column(s) not found: ', paste(missing, collapse = ', '), '. The table has: ',
             paste(utils::head(colnames(df), 20), collapse = ', '), if(ncol(df) > 20) ', ...', call. = FALSE)
    }
    df
}

.DAColors <- function(colors){
    if(length(colors) != 3) stop('colors: give three colors (lower, not significant, higher).')
    unname(colors)
}

# diverging fold change scale; NA (not significant) in na_color; optional symmetric limits with squished outliers
.DAColorScale <- function(colors, na_color, color_limit){
    args <- list(low = colors[1], mid = colors[2], high = colors[3], midpoint = 0,
                 limits = if(is.null(color_limit)) NULL else c(-color_limit, color_limit), oob = scales::squish)
    if(!is.null(na_color)) args$na.value <- na_color
    do.call(scale_color_gradient2, args)
}

# neighbourhood positions (index cells on the reduction) and sizes from a Milo object's neighbourhood graph
.MiloNhoodLayout <- function(milo, reduction){
    for(pkg in c('miloR', 'igraph', 'SingleCellExperiment')){
        if(!requireNamespace(pkg, quietly = TRUE)) stop('a Milo object needs the ', pkg, ' package.')
    }
    graph <- miloR::nhoodGraph(milo)
    if(igraph::vcount(graph) == 0) stop('the Milo object has no neighbourhood graph: run miloR::buildNhoodGraph().')
    if(!(reduction %in% SingleCellExperiment::reducedDimNames(milo))){
        stop("reduction '", reduction, "' not found; the object has: ",
             paste(SingleCellExperiment::reducedDimNames(milo), collapse = ', '))
    }
    index <- as.numeric(igraph::vertex_attr(graph)$name)
    emb <- SingleCellExperiment::reducedDim(milo, reduction)[index, 1:2, drop = FALSE]
    data.frame(Nhood = seq_along(index), x = unname(emb[, 1]), y = unname(emb[, 2]),
               nhood_size = as.numeric(igraph::vertex_attr(graph)$size))
}

# cell coordinates and the label column, for the group labels
.DALabelCells <- function(nhoods, is_milo, reduction, label_by, cells, cells_reduction){
    if(is_milo){
        meta <- as.data.frame(SingleCellExperiment::colData(nhoods))
        if(!(label_by %in% colnames(meta))) stop("label_by: column '", label_by, "' not in the Milo object's colData.")
        emb <- SingleCellExperiment::reducedDim(nhoods, reduction)
        out <- data.frame(x = emb[, 1], y = emb[, 2])
        out[[label_by]] <- meta[[label_by]]
        return(out)
    }
    if(is.null(cells)) stop('label_by with a data frame of neighbourhoods needs the cells (cells = ...).')
    FetchCellData(cells, reduction = cells_reduction, vars = label_by)
}

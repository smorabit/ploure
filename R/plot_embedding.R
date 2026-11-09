#---------------------------------------------------------#
# PlotEmbedding / FeatureEmbedding
#
# Ported from PLEURA's utils/plotting_functions.R (originally from the COMPACT plotting functions; PlotEmbedding ~
# Seurat::DimPlot, FeatureEmbedding ~ Seurat::FeaturePlot). Arguments and defaults are unchanged, so existing figures
# are reproduced; both now also take a plain data frame, through FetchCellData().
#---------------------------------------------------------#

#' Embedding colored by a categorical variable
#'
#' A 2D embedding (UMAP, PCA, ...) of cells colored by a categorical variable, with optional split panels, group
#' labels, custom colors and corner axes.
#'
#' @param seurat_obj A Seurat object or a data frame with one row per cell (see [FetchCellData()]).
#' @param group.by Column (metadata) to color by.
#' @param reduction Coordinates: a reduction name, or two coordinate columns (see [FetchCellData()]). For a data
#'   frame, a single name such as `'umap'` uses the columns `umap_1` and `umap_2`.
#' @param split.by Optional column; returns one plot per level.
#' @param plot_under If TRUE with `split.by`, draws the other cells in grey underneath.
#' @param label If TRUE, labels each group (needs ggrepel).
#' @param point_size,legend_point_size,text_size Point / legend key / label sizes.
#' @param raster Rasterize the point layer (ggrastr).
#' @param order_points `'shuffle'` randomizes the draw order to reduce overplotting.
#' @param x,y Optional columns to use as coordinates instead of `reduction`.
#' @param raster_dpi,raster_scale Rasterization resolution and scale.
#' @param selected Optional subset of groups to show; the other cells are drawn in light grey.
#' @param plot_theme Optional ggplot2 theme to add.
#' @param color_df Colors: a data frame with the groups in its first column and the colors in its second (e.g.
#'   `group`, `colour`), a named color vector, or a palette name such as `'met:Hiroshige'` ([PaletteColors()]).
#'   Default: ggplot's hue palette.
#' @param plot_ratio If TRUE, fix the aspect ratio (`coord_equal()`).
#' @param corner_axes If TRUE, replace the axes with short corner axes ([EmbeddingAxes()]).
#' @param axes_labels Corner-axis labels; default from the reduction name ('UMAP 1', 'UMAP 2', ...).
#' @param axes_frac Corner-axis length as a fraction of the data range.
#' @param label_position Where each label sits: `'centroid'` (mean coordinates) or `'medoid'` (a cell on the data).
#' @param label_map Optional named vector mapping group names to label text.
#' @param label_halo Draw a white halo behind the labels.
#'
#' @return A ggplot (no `split.by`) or a named list of ggplots (one per split level).
#' @export
#' @examples
#' PlotEmbedding(sim_cells, group.by = 'cluster', reduction = 'umap', corner_axes = TRUE,
#'               plot_theme = ThemeEmbedding())
PlotEmbedding <- function(
    seurat_obj,
    group.by,
    reduction = 'umap',
    split.by = NULL,
    plot_under = FALSE,
    label = TRUE,
    point_size = 1,
    legend_point_size = 5,
    text_size = 3,
    raster = TRUE,
    order_points = 'shuffle',
    x = NULL,
    y = NULL,
    raster_dpi = 400,
    raster_scale = 1,
    selected = NULL,
    plot_theme = NULL,
    color_df = NULL,
    plot_ratio = TRUE,
    corner_axes = FALSE,
    axes_labels = NULL,
    axes_frac = 0.15,
    label_position = c('centroid', 'medoid'),
    label_map = NULL,
    label_halo = FALSE
){
    label_position <- match.arg(label_position)
    coords <- if(!is.null(x) && !is.null(y)) c(x, y) else reduction
    plot_df <- FetchCellData(seurat_obj, reduction = coords, vars = c(group.by, split.by))
    plot_df$value <- .AsFactor(plot_df[[group.by]])

    # labels over all cells (the same in every split panel)
    label_layer <- NULL
    if(label){
        label_groups <- plot_df$value
        if(!is.null(selected)) label_groups[!(label_groups %in% selected)] <- NA
        # ggrepel's own spacing defaults, as in the original PlotEmbedding
        label_layer <- GroupLabels(plot_df$x, plot_df$y, label_groups, position = label_position,
                                   label_map = label_map, size = text_size, halo = label_halo,
                                   min.segment.length = 0.5, segment.size = 0.5, box.padding = 0.25)
    }

    colors <- .AsNamedColors(color_df, levels(plot_df$value))

    # restrict to selected groups
    if(!is.null(selected)){
        plot_df$value[!(plot_df$value %in% selected)] <- NA
        plot_df$value <- droplevels(plot_df$value)
    }

    if(identical(order_points, 'shuffle')) plot_df <- plot_df[sample(nrow(plot_df)), ]

    if(is.null(axes_labels)) axes_labels <- EmbeddingAxisLabels(coords[1])
    finish <- function(p, title){
        p <- p +
            scale_color_manual(values = colors, na.value = 'grey90', name = group.by) +
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5)) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size)))
        if(!is.null(label_layer)) p <- p + label_layer
        if(!is.null(plot_theme)) p <- p + plot_theme
        if(plot_ratio) p <- p + coord_equal()
        # corner axes over the full data range (the same position in every split panel)
        if(corner_axes) p <- p + .NoAxesTheme() + EmbeddingAxes(plot_df$x, plot_df$y, frac = axes_frac,
                                                                labels = axes_labels)
        p
    }
    points <- function(df, ...){
        layer <- geom_point(data = df, ...)
        if(raster) ggrastr::rasterise(layer, dpi = raster_dpi, scale = raster_scale) else layer
    }

    if(is.null(split.by)){
        p <- ggplot(plot_df, aes(x = x, y = y, color = value)) +
            points(plot_df, size = point_size)
        return(finish(p, group.by))
    }

    split_groups <- unique(plot_df[[split.by]])
    plot_list <- lapply(split_groups, function(cur_split){
        cur_df <- plot_df[plot_df[[split.by]] %in% cur_split, ]
        p <- ggplot(cur_df, aes(x = x, y = y, color = value))
        if(plot_under){
            under_df <- plot_df[!(plot_df[[split.by]] %in% cur_split), ]
            p <- p + points(under_df, aes(x = x, y = y), inherit.aes = FALSE, size = point_size / 2,
                            color = 'lightgrey')
        }
        p <- p + points(cur_df, size = point_size)
        finish(p, as.character(cur_split))
    })
    names(plot_list) <- as.character(split_groups)
    plot_list
}

#' Embedding colored by a continuous feature
#'
#' Overlays a continuous feature (gene expression or a numeric column) on a 2D embedding. The color limits are
#' percentiles (`'q99'`) or numbers; values above the upper limit are clipped, cells at or below the lower limit are
#' drawn in light grey underneath the others, and high values are drawn on top.
#'
#' @param seurat_obj A Seurat object or a data frame with one row per cell (see [FetchCellData()]).
#' @param features Genes and/or numeric columns.
#' @param raster Rasterize the point layers (ggrastr).
#' @param facet_by Optional column to facet a single-feature plot by.
#' @param layer,assay Assay layer and assay for genes (default: `data` of the default assay).
#' @param reduction Coordinates: a reduction name, or two coordinate columns (see [FetchCellData()]).
#' @param dpi,dpi_scale Rasterization resolution and scale.
#' @param ncol Columns when combining several features.
#' @param combine If TRUE, return one patchwork; else a named list of ggplots.
#' @param order_points TRUE draws high values on top; `'shuffle'` randomizes the order.
#' @param point_size Point size.
#' @param plot_max,plot_min Upper and lower color limits: a percentile string such as `'q99'`, or a number.
#' @param same_range If TRUE, all features share one color range.
#' @param colfunc Palette function (integer -> colors) or a palette name ([PaletteColors()]); default
#'   `viridis::inferno`.
#' @param rev_colors Reverse the palette.
#' @param corner_axes If TRUE, add short corner axes ([EmbeddingAxes()]).
#' @param axes_labels Corner-axis labels; default from the reduction name.
#' @param axes_frac Corner-axis length as a fraction of the data range.
#'
#' @return A ggplot, a patchwork, or a named list of ggplots.
#' @export
#' @examples
#' FeatureEmbedding(sim_cells, features = 'score', reduction = 'umap', plot_min = 'q5', plot_max = 'q99',
#'                  corner_axes = TRUE)
FeatureEmbedding <- function(
    seurat_obj,
    features,
    raster = TRUE,
    facet_by = NULL,
    layer = 'data',
    assay = NULL,
    reduction = 'umap',
    dpi = 200,
    dpi_scale = 0.5,
    ncol = 4,
    combine = TRUE,
    order_points = TRUE,
    point_size = 0.5,
    plot_max = 'q100',
    plot_min = 'q0',
    same_range = FALSE,
    colfunc = NULL,
    rev_colors = TRUE,
    corner_axes = FALSE,
    axes_labels = NULL,
    axes_frac = 0.15
){
    data <- FetchCellData(seurat_obj, reduction = reduction, vars = c(features, facet_by), layer = layer,
                          assay = assay)
    for(f in features){
        if(!is.numeric(data[[f]])) stop("feature '", f, "' is not numeric; use PlotEmbedding for categories.")
    }

    # shared color range across features
    if(same_range){
        pooled <- unlist(data[features])
        plot_min <- .ResolveLimit(pooled, plot_min)
        plot_max <- .ResolveLimit(pooled, plot_max)
    }

    pal <- .ContinuousPalette(colfunc, rev_colors)
    if(is.null(axes_labels)) axes_labels <- EmbeddingAxisLabels(reduction[1])
    points <- function(layer) if(raster) ggrastr::rasterise(layer, dpi = dpi, scale = dpi_scale) else layer

    plot_list <- lapply(features, function(feature){
        plot_df <- data
        plot_df$value <- plot_df[[feature]]
        plot_df <- .PrepareValues(plot_df, plot_min, plot_max, order_points = order_points)
        low_df <- plot_df[plot_df$low, ]
        hi_df <- plot_df[!plot_df$low, ]

        p <- ggplot(hi_df, aes(x = x, y = y, color = value)) +
            points(geom_point(data = low_df, aes(x = x, y = y), inherit.aes = FALSE, color = 'lightgrey',
                              size = point_size)) +
            points(geom_point(size = point_size)) +
            labs(color = feature) +
            scale_color_gradientn(colors = pal, limits = attr(plot_df, 'limits')) +
            coord_equal() +
            theme(
                plot.title = element_text(hjust = 0.5, face = 'plain'),
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_blank()
            )
        if(corner_axes) p <- p + EmbeddingAxes(data$x, data$y, frac = axes_frac, labels = axes_labels)
        p
    })
    names(plot_list) <- features

    if(length(plot_list) == 1){
        p <- plot_list[[1]]
        if(!is.null(facet_by)) p <- p + facet_wrap(vars(.data[[facet_by]]), ncol = ncol)
        return(p)
    }
    if(combine) patchwork::wrap_plots(plot_list, ncol = ncol) else plot_list
}

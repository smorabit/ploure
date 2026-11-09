#---------------------------------------------------------#
# ZoomEmbedding: two-panel "zoom-in" embedding plot
#
# Ported from PLEURA's utils/zoom_embedding.R. Both panels are drawn in one ggplot: the subset embedding is rescaled
# into a square placed to the right of the full embedding, so the zoom lines join the two in the same coordinates.
# The left panel uses the color scale and the right panel the fill scale (points with shape 21 and no stroke), so each
# panel keeps its own scale whatever its type.
#---------------------------------------------------------#

#' Two-panel zoom-in embedding
#'
#' Left: an embedding of the full dataset, with a subset of its cells (for example one lineage) outlined by a polygon
#' ([OutlinePolygon()]). Right: the same subset in its own embedding, drawn in a square panel joined to the outline by
#' zoom lines. Each panel is colored by a categorical variable (with a legend of groups) or by a continuous feature
#' (as in [FeatureEmbedding()]: percentile limits, values above the upper limit clipped, cells at or below the lower
#' limit in light grey underneath, high values on top). Each legend sits under its own panel.
#'
#' @param full,subset Seurat objects or data frames (see [FetchCellData()]): the full dataset and the subset shown in
#'   the zoom panel.
#' @param full_reduction,subset_reduction A reduction name, or the names of two coordinate columns.
#' @param full_color_by,subset_color_by Column or gene to color each panel by. Numeric values give a continuous
#'   feature plot; anything else is categorical.
#' @param highlight_cells Cells of the full dataset to outline. Default: the cells of `subset` found in `full`.
#' @param full_colors,subset_colors Named color vectors or palette names ([PaletteColors()]) for categorical panels
#'   (default: ggplot hue palette).
#' @param layer,assay Assay layer and assay for genes in Seurat objects (default: `data` of the default assay).
#' @param plot_min,plot_max Color limits for continuous panels, as a `'qNN'` percentile string or a number; one value
#'   for both panels or two (full, subset).
#' @param colfunc,rev_colors Palette function or palette name for continuous panels (default `viridis::inferno`,
#'   reversed).
#' @param full_label_by,subset_label_by Optional categorical column to label in each panel, at each group's medoid
#'   cell ([GroupLabels()]).
#' @param label_map Optional named vector mapping group names to label text (e.g. `c(P07 = 'P07 LC')`).
#' @param label_size Label text size.
#' @param outline If TRUE, outline the highlighted cells with a density-contour polygon; else a rectangle.
#' @param outline_level,outline_trim,outline_bandwidth,outline_min_frac Passed to [OutlinePolygon()] as `level`,
#'   `trim`, `bandwidth` and `min_frac`.
#' @param corner_axes Draw corner axes ([EmbeddingAxes()]) in the full and subset panels: one value or two.
#' @param axes_labels Optional `list(full = , subset = )` of corner-axis labels; default from each reduction's name.
#' @param axes_frac Corner-axis length as a fraction of each panel's data range.
#' @param titles Panel titles (full, subset).
#' @param legend_titles Legend titles (full, subset); default the color_by names.
#' @param legend_ncol Legend columns (full, subset) for categorical legends.
#' @param legend_colors Optional `list(full = , subset = )` of named color vectors that replace a categorical panel's
#'   legend, e.g. a few color families when the points are colored by patient and the patients are labeled.
#' @param point_size Point sizes (full, subset).
#' @param raster Rasterize the point layers.
#' @param raster_dpi Resolution of the rasterized point layers.
#' @param gap Gap between the panels, as a fraction of the full embedding's x range.
#' @param dash Line style shared by the outline and the zoom lines.
#' @param legend_height Height of the legend row relative to the plots.
#' @param seed Seed for the point shuffle and the medoid sampling (one random stream, so labels are placed as in the
#'   PLEURA version).
#'
#' @return A patchwork object (the plot, with a legend under each panel).
#' @export
#' @examples
#' lymphoid <- sim_cells[sim_cells$lineage == 'Lymphoid', ]
#' ZoomEmbedding(
#'     full = sim_cells, subset = lymphoid,
#'     full_reduction = c('umap_1', 'umap_2'), subset_reduction = c('sub_umap_1', 'sub_umap_2'),
#'     full_color_by = 'cluster', subset_color_by = 'score', plot_min = 'q5', plot_max = 'q99',
#'     titles = c('All cells', 'Lymphoid cells')
#' )
ZoomEmbedding <- function(
    full,
    subset,
    full_reduction = 'umap',
    subset_reduction = 'umap',
    full_color_by,
    subset_color_by,
    highlight_cells = NULL,
    full_colors = NULL,
    subset_colors = NULL,
    layer = 'data',
    assay = NULL,
    plot_min = 'q0',
    plot_max = 'q100',
    colfunc = NULL,
    rev_colors = TRUE,
    full_label_by = NULL,
    subset_label_by = NULL,
    label_map = NULL,
    label_size = 2.3,
    outline = TRUE,
    outline_level = 0.995,
    outline_trim = 0.01,
    outline_bandwidth = 0.04,
    outline_min_frac = 0.01,
    corner_axes = TRUE,
    axes_labels = NULL,
    axes_frac = 0.15,
    titles = c('All cells', 'Subset'),
    legend_titles = NULL,
    legend_ncol = c(3, 2),
    legend_colors = list(full = NULL, subset = NULL),
    point_size = c(0.15, 0.4),
    raster = TRUE,
    raster_dpi = 600,
    gap = 0.25,
    dash = list(color = 'black', linewidth = 0.4, linetype = 'dashed'),
    legend_height = 0.3,
    seed = 12345
){
    if(is.null(legend_titles)) legend_titles <- c(full_color_by, subset_color_by)
    corner_axes <- rep(corner_axes, length.out = 2)
    if(is.null(axes_labels)) axes_labels <- list(full = NULL, subset = NULL)
    if(is.null(axes_labels$full)) axes_labels$full <- EmbeddingAxisLabels(full_reduction[1])
    if(is.null(axes_labels$subset)) axes_labels$subset <- EmbeddingAxisLabels(subset_reduction[1])
    plot_min <- rep(plot_min, length.out = 2)
    plot_max <- rep(plot_max, length.out = 2)
    set.seed(seed)

    # panel data: cell, x, y, value (+ label)
    full_df <- .ZoomData(full, full_reduction, full_color_by, full_label_by, layer, assay)
    sub_df <- .ZoomData(subset, subset_reduction, subset_color_by, subset_label_by, layer, assay)
    full_df <- .PrepareValues(full_df, plot_min[1], plot_max[1])
    sub_df <- .PrepareValues(sub_df, plot_min[2], plot_max[2])

    # cells to outline
    if(is.null(highlight_cells)) highlight_cells <- intersect(sub_df$cell, full_df$cell)
    hl <- full_df$cell %in% highlight_cells
    if(!any(hl)) stop('no highlighted cells found in the full dataset (check cell names or highlight_cells).')

    #---------------------------------------------------------#
    # geometry
    #---------------------------------------------------------#

    x_rng <- range(full_df$x)
    y_rng <- range(full_df$y)
    shape <- OutlinePolygon(full_df$x, full_df$y, hl, shape = if(outline) 'contour' else 'rectangle',
                            level = outline_level, trim = outline_trim, bandwidth = outline_bandwidth,
                            min_frac = outline_min_frac)

    # the subset panel: a square with the height of the full embedding, after a gap
    panel_size <- diff(y_rng)
    gap_size <- gap * diff(x_rng)
    panel <- data.frame(xmin = x_rng[2] + gap_size, xmax = x_rng[2] + gap_size + panel_size,
                        ymin = y_rng[1], ymax = y_rng[2])
    # inner margin of the subset panel; wider with corner axes, so the axes and their labels clear the frame
    inner <- (if(corner_axes[2]) 0.1 else 0.04) * panel_size
    s <- (panel_size - 2 * inner) / max(diff(range(sub_df$x)), diff(range(sub_df$y)))
    sub_df$x <- panel$xmin + inner + (sub_df$x - min(sub_df$x)) * s +
        (panel_size - 2 * inner - diff(range(sub_df$x)) * s) / 2
    sub_df$y <- panel$ymin + inner + (sub_df$y - min(sub_df$y)) * s +
        (panel_size - 2 * inner - diff(range(sub_df$y)) * s) / 2
    sub_rng <- list(x = range(sub_df$x), y = range(sub_df$y))

    # zoom lines from the top and bottom of the outline to the left corners of the panel
    top <- shape[which.max(shape$y), ]
    bottom <- shape[which.min(shape$y), ]
    zoom_lines <- data.frame(x = c(top$x, bottom$x), y = c(top$y, bottom$y),
                             xend = panel$xmin, yend = c(panel$ymax, panel$ymin))

    title_df <- data.frame(x = c(mean(x_rng), mean(c(panel$xmin, panel$xmax))),
                           y = y_rng[2] + 0.06 * diff(y_rng), label = titles)

    #---------------------------------------------------------#
    # plot: left panel on the color scale, right panel on the fill scale
    #---------------------------------------------------------#

    pal <- .ContinuousPalette(colfunc, rev_colors)
    full_colors <- .ZoomColors(full_df, full_colors)
    subset_colors <- .ZoomColors(sub_df, subset_colors)
    points <- function(layer) if(raster) ggrastr::rasterise(layer, dpi = raster_dpi) else layer

    p <- ggplot()
    if(attr(full_df, 'continuous')){
        p <- p + points(geom_point(data = full_df[full_df$low, ], aes(x = x, y = y), color = 'lightgrey',
                                   size = point_size[1], stroke = 0))
        full_df <- .SubsetRows(full_df, !full_df$low)
    }
    p <- p +
        points(geom_point(data = full_df, aes(x = x, y = y, color = value), size = point_size[1], stroke = 0)) +
        geom_path(data = shape, aes(x = x, y = y, group = group), color = dash$color, linewidth = dash$linewidth,
                  linetype = dash$linetype) +
        geom_segment(data = zoom_lines, aes(x = x, y = y, xend = xend, yend = yend), color = dash$color,
                     linewidth = dash$linewidth, linetype = dash$linetype) +
        geom_rect(data = panel, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA,
                  color = 'black', linewidth = 0.4)
    if(attr(sub_df, 'continuous')){
        p <- p + points(geom_point(data = sub_df[sub_df$low, ], aes(x = x, y = y), color = 'lightgrey',
                                   size = point_size[2], stroke = 0))
        sub_df <- .SubsetRows(sub_df, !sub_df$low)
    }
    p <- p +
        points(geom_point(data = sub_df, aes(x = x, y = y, fill = value), shape = 21, size = point_size[2],
                          stroke = 0)) +
        geom_text(data = title_df, aes(x = x, y = y, label = label), size = 4.2) +
        .ZoomScale('color', full_df, full_colors, pal) +
        .ZoomScale('fill', sub_df, subset_colors, pal)

    # corner axes (over each panel's full data range, grey cells included)
    if(corner_axes[1]) p <- p + EmbeddingAxes(x_rng, y_rng, frac = axes_frac, labels = axes_labels$full)
    if(corner_axes[2]) p <- p + EmbeddingAxes(sub_rng$x, sub_rng$y, frac = axes_frac, labels = axes_labels$subset)

    # labels at each group's medoid
    if(!is.null(full_label_by)){
        p <- p + GroupLabels(full_df$x, full_df$y, full_df$label, label_map = label_map, size = label_size,
                             fontface = 'bold', medoid_seed = NULL)
    }
    if(!is.null(subset_label_by)){
        p <- p + GroupLabels(sub_df$x, sub_df$y, sub_df$label, label_map = label_map, size = label_size,
                             fontface = 'bold', medoid_seed = NULL)
    }

    p <- p +
        coord_equal(clip = 'off') +
        theme_void() +
        theme(legend.position = 'none', plot.margin = margin(15, 5, 0, 5))

    #---------------------------------------------------------#
    # legends, each under its own panel
    #---------------------------------------------------------#

    leg_full <- .ZoomLegend(full_df, full_colors, pal, legend_titles[1], legend_ncol[1], legend_colors$full)
    leg_sub <- .ZoomLegend(sub_df, subset_colors, pal, legend_titles[2], legend_ncol[2], legend_colors$subset)
    legends <- (patchwork::wrap_elements(leg_full) | patchwork::wrap_elements(leg_sub)) +
        patchwork::plot_layout(widths = c(diff(x_rng) + gap_size / 2, panel_size + gap_size / 2))
    p / legends + patchwork::plot_layout(heights = c(1, legend_height))
}

# one panel's cell, x, y, value (+ label) from the fetch core
.ZoomData <- function(obj, reduction, color_by, label_by, layer, assay){
    df <- FetchCellData(obj, reduction = reduction, vars = c(color_by, label_by), layer = layer, assay = assay)
    out <- data.frame(cell = df$cell, x = df$x, y = df$y, value = df[[color_by]])
    if(!is.null(label_by)) out$label <- df[[label_by]]
    out
}

.ZoomColors <- function(df, colors){
    if(attr(df, 'continuous')) return(colors)
    .AsNamedColors(colors, levels(df$value))
}

.ZoomScale <- function(aesthetic, df, colors, pal){
    if(attr(df, 'continuous')){
        scale_continuous <- if(aesthetic == 'color') scale_color_gradientn else scale_fill_gradientn
        scale_continuous(colors = pal, limits = attr(df, 'limits'))
    } else {
        scale_manual <- if(aesthetic == 'color') scale_color_manual else scale_fill_manual
        scale_manual(values = colors, drop = FALSE)
    }
}

.ZoomLegend <- function(df, colors, pal, title, ncol, override = NULL){
    if(attr(df, 'continuous')) return(LegendGrob(limits = attr(df, 'limits'), palette = pal, title = title))
    if(!is.null(override)) return(LegendGrob(colors = override, title = title, ncol = ncol))
    LegendGrob(colors = colors[levels(df$value)], title = title, ncol = ncol)
}

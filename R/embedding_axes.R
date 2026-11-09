#---------------------------------------------------------#
# Corner axes for embeddings
#---------------------------------------------------------#

#' Short corner axes for an embedding
#'
#' Two arrows from the lower-left corner of the data, each a fixed fraction of the data range, labeled (e.g. "UMAP
#' 1", "UMAP 2"), as is usual for UMAPs in papers. The arrows are drawn as annotation segments rather than theme axes,
#' so they extend the plot limits as needed and are never clipped. The plot's own axes should be blank (e.g.
#' [ThemeEmbedding()] or `theme_void()`).
#'
#' @param x,y The plotted coordinates (vectors), used for the data range.
#' @param frac Length of each arm as a fraction of the data range (default 0.15).
#' @param labels Axis labels (x, y). See [EmbeddingAxisLabels()] to derive them from a reduction name.
#' @param equal_length If TRUE, both arms take `frac` of the larger of the two ranges, so they look equal under
#'   `coord_equal()`; if FALSE, each takes `frac` of its own range.
#' @param offset Gap between the data and the arms, as a fraction of each range.
#' @param linewidth,color Arm line width and color.
#' @param text_size Label size (ggplot text size, mm).
#' @param arrow A [grid::arrow()] for the arm ends, or NULL for plain lines.
#'
#' @return A list of ggplot layers, to add to a plot.
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(sim_cells, aes(x = umap_1, y = umap_2)) +
#'     geom_point(size = 0.3) +
#'     EmbeddingAxes(sim_cells$umap_1, sim_cells$umap_2) +
#'     coord_equal() +
#'     ThemeEmbedding()
EmbeddingAxes <- function(
    x,
    y,
    frac = 0.15,
    labels = c('UMAP 1', 'UMAP 2'),
    equal_length = TRUE,
    offset = 0.02,
    linewidth = 0.4,
    color = 'black',
    text_size = 3,
    arrow = grid::arrow(length = grid::unit(0.06, 'inches'), type = 'closed')
){
    xr <- range(x, na.rm = TRUE)
    yr <- range(y, na.rm = TRUE)
    len <- frac * c(diff(xr), diff(yr))
    if(equal_length) len <- rep(max(len), 2)
    x0 <- xr[1] - offset * diff(xr)
    y0 <- yr[1] - offset * diff(yr)
    list(
        annotate('segment', x = x0, y = y0, xend = x0 + len[1], yend = y0, linewidth = linewidth, color = color,
                 arrow = arrow, linejoin = 'mitre'),
        annotate('segment', x = x0, y = y0, xend = x0, yend = y0 + len[2], linewidth = linewidth, color = color,
                 arrow = arrow, linejoin = 'mitre'),
        annotate('text', x = x0 + len[1] / 2, y = y0, label = labels[1], vjust = 1.5, size = text_size,
                 color = color),
        annotate('text', x = x0, y = y0 + len[2] / 2, label = labels[2], angle = 90, vjust = -0.6,
                 size = text_size, color = color)
    )
}

#' Corner-axis labels from a reduction name
#'
#' `umap` gives "UMAP 1", "UMAP 2"; `pca` gives "PC 1", "PC 2"; `tsne` gives "tSNE 1", "tSNE 2". Matching ignores case
#' and works on names such as `umap_scvi` or `X_umap`. Anything else is used as the stem.
#'
#' @param reduction A reduction name (or the first of two coordinate column names).
#'
#' @return A character vector of length 2.
#' @export
#' @examples
#' EmbeddingAxisLabels('umap_harmony')
#' EmbeddingAxisLabels('pca')
EmbeddingAxisLabels <- function(reduction){
    reduction <- reduction[1]
    stem <- if(grepl('umap', reduction, ignore.case = TRUE)){
        'UMAP'
    } else if(grepl('tsne', reduction, ignore.case = TRUE)){
        'tSNE'
    } else if(grepl('pca|^pc', reduction, ignore.case = TRUE)){
        'PC'
    } else {
        reduction
    }
    paste(stem, 1:2)
}

# blank axes, for plots that use EmbeddingAxes
.NoAxesTheme <- function(){
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )
}

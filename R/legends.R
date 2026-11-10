#---------------------------------------------------------#
# Legend helpers
#---------------------------------------------------------#

#' Build a standalone legend
#'
#' A legend grob for a categorical color key or a continuous color bar, to place under a panel with patchwork
#' (`p / patchwork::wrap_elements(legend)`). A categorical key can show any named colors, e.g. a few color families
#' instead of the 44 patients a panel is colored by.
#'
#' @param colors Named color vector for a categorical key (names = key labels, in order).
#' @param limits For a color bar: the value range.
#' @param palette For a color bar: the colors along the bar.
#' @param title Legend title.
#' @param ncol Columns of a categorical key.
#' @param key_size Point size of the categorical keys.
#' @param bar_width,bar_height Color-bar size in cm.
#' @param position Legend position (`'bottom'`, `'right'`, ...).
#' @param title_size,text_size Title and key-label text sizes in points (default: the current theme's legend title
#'   and legend text sizes).
#'
#' @return A grob.
#' @export
#' @examples
#' leg <- LegendGrob(colors = c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c'), title = 'Lineage')
#' grid::grid.draw(leg)
LegendGrob <- function(
    colors = NULL,
    limits = NULL,
    palette = NULL,
    title = NULL,
    ncol = 1,
    key_size = 2.5,
    bar_width = 5,
    bar_height = 0.35,
    position = 'bottom',
    title_size = ThemeTextSize(element = 'legend.title', unit = 'pt'),
    text_size = ThemeTextSize(element = 'legend.text', unit = 'pt')
){
    base <- theme_void() +
        theme(
            legend.position = position,
            legend.title = element_text(size = title_size, color = 'black'),
            legend.text = element_text(size = text_size, color = 'black'),
            legend.key.spacing.x = grid::unit(0.5, 'cm'),
            legend.frame = .ColorbarFrame()
        )
    if(!is.null(colors)){
        groups <- names(colors)
        p <- ggplot(data.frame(x = 1, y = 1, value = factor(groups, levels = groups)),
                    aes(x = x, y = y, color = value)) +
            geom_point() +
            scale_color_manual(values = colors, name = title, drop = FALSE) +
            guides(color = guide_legend(override.aes = list(size = key_size), ncol = ncol,
                                        title.position = 'top')) +
            base
    } else if(!is.null(limits)){
        if(is.null(palette)) palette <- .ContinuousPalette()
        p <- ggplot(data.frame(x = 1, y = 1, value = limits), aes(x = x, y = y, color = value)) +
            geom_point() +
            scale_color_gradientn(colors = palette, limits = limits, name = title) +
            guides(color = guide_colorbar(title.position = 'top', theme = theme(
                legend.key.width = grid::unit(bar_width, 'cm'),
                legend.key.height = grid::unit(bar_height, 'cm')
            ))) +
            base
    } else {
        stop('give either colors (categorical key) or limits (color bar).')
    }
    ExtractLegend(p)
}

#' Extract the legend of a ggplot
#'
#' @param p A ggplot.
#'
#' @return The legend as a grob (a zero grob if the plot has no legend).
#' @export
ExtractLegend <- function(p){
    g <- ggplotGrob(p)
    idx <- grep('^guide-box', g$layout$name)
    boxes <- g$grobs[idx]
    boxes <- boxes[!vapply(boxes, inherits, logical(1), 'zeroGrob')]
    if(length(boxes) == 0) return(grid::nullGrob())
    boxes[[1]]
}

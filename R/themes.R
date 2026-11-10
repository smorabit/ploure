#---------------------------------------------------------#
# Themes
#---------------------------------------------------------#

#' Publication theme
#'
#' A clean theme (cowplot-like) that encodes the figure rules: black text, axis lines and ticks; centered plain
#' titles; bold panel tags and nothing else bold; blank strip backgrounds; a black border around color bars. The remaining rules are for the labels you
#' write: sentence case ("Cells per sample", not "Cells Per Sample") and American spelling.
#'
#' It is a complete theme, so adding it to a plot (`p + ThemePublication()`) replaces the theme settings the plotting
#' functions set (italic gene names, strip placement). Set it once instead, with `theme_set(ThemePublication())`.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param line_width Width of axis lines and ticks.
#' @param title_size,tag_size Plot title and panel tag sizes (default `base_size + 1` and `base_size + 2`).
#'
#' @return A ggplot2 theme.
#' @export
#' @examples
#' library(ggplot2)
#' old <- theme_set(ThemePublication())
#' ggplot(sim_cells, aes(x = cluster)) +
#'     geom_bar() +
#'     labs(title = 'Cells per cluster', tag = 'a')
#' theme_set(old)
ThemePublication <- function(
    base_size = 11,
    base_family = '',
    line_width = 0.3,
    title_size = base_size + 1,
    tag_size = base_size + 2
){
    theme_classic(base_size = base_size, base_family = base_family) +
        theme(
            text = element_text(color = 'black'),
            axis.text = element_text(color = 'black', size = base_size - 1),
            axis.line = element_line(color = 'black', linewidth = line_width),
            axis.ticks = element_line(color = 'black', linewidth = line_width),
            plot.title = element_text(face = 'plain', hjust = 0.5, size = title_size),
            plot.subtitle = element_text(hjust = 0.5),
            plot.tag = element_text(face = 'bold', size = tag_size),
            strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = base_size - 1),
            # strips beside an axis go outside it; plots stacked with patchwork then keep their axes on the panels
            strip.placement = 'outside',
            legend.title = element_text(size = base_size - 1),
            legend.text = element_text(size = base_size - 2),
            legend.key.size = .LegendKeySize(base_size),
            legend.frame = .ColorbarFrame(line_width)
        )
}

# black border around continuous color bars (ggplot2 >= 3.5 draws legend.frame on color bars only)
.ColorbarFrame <- function(line_width = 0.3) element_rect(color = 'black', fill = NA, linewidth = line_width)

# ggplot2's legend key size (1.2 lines at grid's 12 pt: 17.28 big points) scaled to the base size; unchanged at 11
.LegendKeySize <- function(base_size) grid::unit(17.28 * base_size / 11, 'bigpts')

#' Embedding theme
#'
#' No axes, a centered plain title, a 1:1 aspect ratio and a black border around color bars; pair with
#' [EmbeddingAxes()] for corner axes.
#'
#' @param title_size Title size.
#' @param aspect_ratio Panel aspect ratio (NULL to leave it free, e.g. when the plot already uses `coord_equal()`).
#' @param base_size Base font size of the legend text.
#'
#' @return A ggplot2 theme.
#' @export
ThemeEmbedding <- function(title_size = 11, aspect_ratio = 1, base_size = 11){
    theme_void(base_size = base_size) +
        theme(
            plot.title = element_text(face = 'plain', hjust = 0.5, size = title_size),
            plot.tag = element_text(face = 'bold', size = base_size + 2),
            legend.title = element_text(size = base_size - 1),
            legend.text = element_text(size = base_size - 2),
            legend.key.size = .LegendKeySize(base_size),
            legend.frame = .ColorbarFrame(),
            aspect.ratio = aspect_ratio
        )
}

#' Text size from the current theme
#'
#' The size of a theme text element (by default the axis tick labels) of the current theme, in mm for `geom_text()`
#' and `annotate()` (`size = `) or in points for `element_text()`. The plotting functions size the text they draw
#' inside a plot with it, so that text follows the theme set with `theme_set()` (call that first) and stays at or
#' below the tick labels: the plot title is the largest text, then the axis and legend titles, then the tick labels,
#' legend keys and annotations.
#'
#' @param rel Multiplier of the element's size.
#' @param element Theme element (e.g. `'axis.text'`, `'axis.title'`, `'legend.text'`).
#' @param unit `'mm'` (geom text sizes) or `'pt'` (theme sizes).
#' @param theme The theme (default the current one, [ggplot2::theme_get()]).
#'
#' @return A number.
#' @export
#' @examples
#' library(ggplot2)
#' old <- theme_set(ThemePublication())
#' ThemeTextSize()                     # tick labels: 10 pt, in mm
#' ThemeTextSize(0.8)                  # smaller annotations
#' ThemeTextSize(element = 'axis.title', unit = 'pt')
#' ggplot(sim_cells, aes(x = cluster)) +
#'     geom_bar() +
#'     annotate('text', x = 1, y = 5, label = 'n', size = ThemeTextSize(0.9))
#' theme_set(old)
ThemeTextSize <- function(rel = 1, element = 'axis.text', unit = c('mm', 'pt'), theme = ggplot2::theme_get()){
    unit <- match.arg(unit)
    size <- ggplot2::calc_element(element, theme)$size
    if(is.null(size)) size <- ggplot2::calc_element('text', theme)$size
    if(is.null(size)) size <- 11
    size <- size * rel
    if(unit == 'mm') size / ggplot2::.pt else size
}

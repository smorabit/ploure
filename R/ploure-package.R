#' ploure: plotting utilities for single-cell data
#'
#' Embeddings (with corner axes, cell outlines and zoom-in panels), publication themes and palettes, per-group metric
#' distributions and marker dot plots. Every plot draws from a tidy per-cell data frame ([FetchCellData()]), taken
#' from a Seurat object or a plain data frame.
#'
#' @import ggplot2
#' @keywords internal
'_PACKAGE'

# columns referenced in aes()
utils::globalVariables(c(
    'x', 'y', 'xend', 'yend', 'xmin', 'xmax', 'ymin', 'ymax', 'value', 'group', 'label', 'text', 'color', 'hex',
    'condition', 'gene', 'frac', 'expr', 'facet', 'metric', 'gene_group', 'x_group', 'w', 'h', 'fill', 'border',
    'linetype', 'title', 'title_y', 'body_y', 'body', 'row', 'stage', 'n', 'facet_value'
))

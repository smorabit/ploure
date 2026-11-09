#---------------------------------------------------------#
# Outline a set of cells on any embedding
#
# A density contour around the highlighted cells: the 1% tails on each axis are trimmed first, so scattered cells do
# not add islands, then the contour that keeps 99.5% of the remaining cells is drawn, with every island holding at
# least 1% of them (a subset can sit in several islands, e.g. B and plasma cells).
#---------------------------------------------------------#

#' Outline polygon around a set of cells
#'
#' @param x,y Coordinates of all cells of the embedding (the full range sets the bandwidth and padding).
#' @param highlight Logical vector (or row indices) of the cells to outline; default all cells.
#' @param shape `'contour'` (density contour) or `'rectangle'` (bounding box of the trimmed cells).
#' @param level Fraction of the (trimmed) highlighted cells kept inside the contour.
#' @param trim Highlighted cells outside these per-axis quantiles (`trim`, `1 - trim`) are left out of the density.
#' @param bandwidth Kernel bandwidth as a fraction of the x range of all cells.
#' @param min_frac Contour islands holding at least this fraction of the highlighted cells are kept (at least the
#'   largest island is always kept).
#' @param n Grid size of the density estimate.
#'
#' @return A data frame of closed paths with columns `x`, `y` and `group` (one group per island).
#' @export
#' @examples
#' poly <- OutlinePolygon(sim_cells$umap_1, sim_cells$umap_2, sim_cells$lineage == 'Lymphoid')
#' head(poly)
OutlinePolygon <- function(
    x,
    y,
    highlight = NULL,
    shape = c('contour', 'rectangle'),
    level = 0.995,
    trim = 0.01,
    bandwidth = 0.04,
    min_frac = 0.01,
    n = 300
){
    shape <- match.arg(shape)
    if(is.null(highlight)) highlight <- rep(TRUE, length(x))
    if(is.logical(highlight)) highlight <- which(highlight)
    if(length(highlight) == 0) stop('no cells to outline.')
    hx <- x[highlight]
    hy <- y[highlight]
    x_rng <- range(x, na.rm = TRUE)
    pad <- 0.03 * diff(x_rng)
    q <- c(trim, 1 - trim)
    core <- hx >= stats::quantile(hx, q[1]) & hx <= stats::quantile(hx, q[2]) &
        hy >= stats::quantile(hy, q[1]) & hy <= stats::quantile(hy, q[2])
    cx <- hx[core]
    cy <- hy[core]

    if(shape == 'rectangle'){
        xr <- range(cx) + c(-1, 1) * pad
        yr <- range(cy) + c(-1, 1) * pad
        return(data.frame(x = xr[c(1, 2, 2, 1, 1)], y = yr[c(1, 1, 2, 2, 1)], group = 1))
    }

    if(!requireNamespace('MASS', quietly = TRUE)) stop("shape = 'contour' needs the MASS package.")
    dens <- MASS::kde2d(cx, cy, n = n, h = rep(bandwidth * diff(x_rng), 2),
                        lims = c(range(cx) + c(-1, 1) * pad, range(cy) + c(-1, 1) * pad))
    cell_dens <- dens$z[cbind(findInterval(cx, dens$x), findInterval(cy, dens$y))]
    contours <- grDevices::contourLines(dens$x, dens$y, dens$z, levels = stats::quantile(cell_dens, 1 - level))
    # keep the islands that hold enough of the highlighted cells (at least the largest one)
    frac <- sapply(contours, function(cc) mean(.InPolygon(cx, cy, cc$x, cc$y)))
    keep <- which(frac >= min_frac)
    if(length(keep) == 0) keep <- which.max(frac)
    do.call(rbind, lapply(seq_along(keep), function(i){
        cc <- contours[[keep[i]]]
        data.frame(x = c(cc$x, cc$x[1]), y = c(cc$y, cc$y[1]), group = i)
    }))
}

#' Outline a set of cells on an embedding plot
#'
#' Adds the [OutlinePolygon()] of the highlighted cells to any embedding plot (`p + OutlineCells(...)`).
#'
#' @inheritParams OutlinePolygon
#' @param color,linewidth,linetype Outline line style.
#' @param ... Further arguments to [OutlinePolygon()] (`shape`, `level`, `trim`, `bandwidth`, `min_frac`).
#'
#' @return A ggplot layer.
#' @export
#' @examples
#' p <- PlotEmbedding(sim_cells, group.by = 'cluster', x = 'umap_1', y = 'umap_2', label = FALSE, raster = FALSE)
#' p + OutlineCells(sim_cells$umap_1, sim_cells$umap_2, sim_cells$lineage == 'Lymphoid')
OutlineCells <- function(
    x,
    y,
    highlight = NULL,
    color = 'black',
    linewidth = 0.4,
    linetype = 'dashed',
    ...
){
    poly <- OutlinePolygon(x, y, highlight, ...)
    geom_path(data = poly, aes(x = x, y = y, group = group), inherit.aes = FALSE, color = color,
              linewidth = linewidth, linetype = linetype)
}

# point-in-polygon by ray casting (vectorized over points)
.InPolygon <- function(px, py, vx, vy){
    inside <- rep(FALSE, length(px))
    j <- length(vx)
    for(i in seq_along(vx)){
        crosses <- ((vy[i] > py) != (vy[j] > py)) &
            (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
        inside <- xor(inside, crosses)
        j <- i
    }
    inside
}

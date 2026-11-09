#---------------------------------------------------------#
# Data fetch core: any supported input -> one tidy per-cell data frame
#
# Every plotting function in the package draws from the data frame returned here, so supporting a new input type
# (SingleCellExperiment, a table exported from h5ad) only means adding a branch to FetchCellData().
#---------------------------------------------------------#

#' Fetch a per-cell data frame for plotting
#'
#' Turns a Seurat object or a plain data frame into one tidy data frame with a `cell` column, optional `x` / `y`
#' coordinates and one column per requested variable. All plotting functions in the package use it, so cached
#' embedding tables can be plotted directly, without wrapping them in a Seurat object.
#'
#' @param object A Seurat object or a data frame (one row per cell).
#' @param reduction Where the coordinates come from, or NULL for none:
#'   * one string: a reduction of a Seurat object (its first two dimensions are used); for a data frame, the
#'     columns `<reduction>_1` and `<reduction>_2` (any case, e.g. `umap_1`, `UMAP_2`);
#'   * two strings: two columns holding the coordinates (metadata columns of a Seurat object, or columns of a data
#'     frame).
#' @param vars Columns (metadata or data-frame columns) or, for a Seurat object, features of the assay to fetch.
#'   Metadata columns take precedence over features with the same name.
#' @param layer,assay Assay layer and assay for features of a Seurat object (default: the `data` layer of the
#'   default assay).
#' @param cell_col For a data frame: the column holding cell names (default `'cell'`; rownames if absent).
#'
#' @return A data frame with `cell`, then `x` and `y` (if `reduction` is given), then the columns in `vars` under their
#'   own names.
#' @export
#' @examples
#' plot_df <- FetchCellData(sim_cells, reduction = c('umap_1', 'umap_2'), vars = c('cluster', 'score'))
#' head(plot_df)
FetchCellData <- function(
    object,
    reduction = NULL,
    vars = NULL,
    layer = 'data',
    assay = NULL,
    cell_col = 'cell'
){
    vars <- unique(vars)
    if(inherits(object, 'Seurat')){
        .FetchSeurat(object, reduction, vars, layer, assay)
    } else if(is.data.frame(object)){
        .FetchDataFrame(object, reduction, vars, cell_col)
    } else {
        stop('unsupported input of class ', class(object)[1], ': use a Seurat object or a data frame.')
    }
}

.FetchSeurat <- function(object, reduction, vars, layer, assay){
    meta <- object@meta.data
    out <- data.frame(cell = rownames(meta))
    if(length(reduction) == 1){
        if(!(reduction %in% names(object@reductions))) stop("reduction '", reduction, "' not found.")
        emb <- object@reductions[[reduction]]@cell.embeddings
        out$x <- emb[out$cell, 1]
        out$y <- emb[out$cell, 2]
    } else if(length(reduction) == 2){
        missing <- setdiff(reduction, colnames(meta))
        if(length(missing) > 0) stop('coordinate columns not found in the metadata: ', paste(missing, collapse = ', '))
        out$x <- meta[[reduction[1]]]
        out$y <- meta[[reduction[2]]]
    }
    in_meta <- vars[vars %in% colnames(meta)]
    features <- setdiff(vars, in_meta)
    for(v in in_meta) out[[v]] <- meta[[v]]
    if(length(features) > 0){
        if(is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
        X <- SeuratObject::LayerData(object, assay = assay, layer = layer)
        missing <- setdiff(features, rownames(X))
        if(length(missing) > 0){
            stop('not found in the metadata or the features of assay ', assay, ': ', paste(missing, collapse = ', '))
        }
        X <- X[features, out$cell, drop = FALSE]
        for(f in features) out[[f]] <- as.numeric(X[f, ])
    }
    out[, c('cell', intersect(c('x', 'y'), colnames(out)), vars), drop = FALSE]
}

.FetchDataFrame <- function(object, reduction, vars, cell_col){
    object <- as.data.frame(object, check.names = FALSE)
    cells <- if(cell_col %in% colnames(object)) as.character(object[[cell_col]]) else rownames(object)
    out <- data.frame(cell = cells)
    if(length(reduction) == 1){
        # a reduction name: look for its <name>_1 / <name>_2 columns (any case)
        cand <- match(tolower(paste0(reduction, c('_1', '_2'))), tolower(colnames(object)))
        if(anyNA(cand)){
            stop("for a data frame, the reduction must name its two coordinate columns (or columns '", reduction,
                 "_1' and '", reduction, "_2' must exist).")
        }
        reduction <- colnames(object)[cand]
    }
    if(length(reduction) == 2){
        missing <- setdiff(reduction, colnames(object))
        if(length(missing) > 0) stop('coordinate columns not found: ', paste(missing, collapse = ', '))
        out$x <- object[[reduction[1]]]
        out$y <- object[[reduction[2]]]
    }
    missing <- setdiff(vars, colnames(object))
    if(length(missing) > 0) stop('columns not found in the data frame: ', paste(missing, collapse = ', '))
    for(v in vars) out[[v]] <- object[[v]]
    out
}

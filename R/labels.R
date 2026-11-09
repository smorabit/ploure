#---------------------------------------------------------#
# Group labels on embeddings
#---------------------------------------------------------#

#' Repelled group labels on an embedding
#'
#' One label per group, placed at the group's medoid cell (a position on the data, also for curved or split groups)
#' or at its centroid, and repelled with ggrepel. Labels carry a white halo by default.
#'
#' @param x,y Cell coordinates.
#' @param group Group of each cell; one label per group (NA cells are ignored).
#' @param position `'medoid'` (the cell closest on average to the others in its group) or `'centroid'` (mean
#'   coordinates).
#' @param label_map Optional named vector mapping group names to label text (e.g. `c(P07 = 'P07 LC')`); unmapped
#'   groups keep their name.
#' @param size Label text size.
#' @param color Label color.
#' @param fontface Label font face.
#' @param halo Draw a white halo behind the text.
#' @param seed Seed for the label repulsion.
#' @param medoid_seed Seed for the medoid sampling (groups of more than 500 cells are sampled); NULL continues the
#'   session's random stream. The global random state is restored afterwards unless NULL.
#' @param ... Further arguments to [ggrepel::geom_text_repel()].
#'
#' @return A ggplot layer.
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(sim_cells, aes(x = umap_1, y = umap_2, color = cluster)) +
#'     geom_point(size = 0.3) +
#'     GroupLabels(sim_cells$umap_1, sim_cells$umap_2, sim_cells$cluster) +
#'     ThemeEmbedding() +
#'     theme(legend.position = 'none')
GroupLabels <- function(
    x,
    y,
    group,
    position = c('medoid', 'centroid'),
    label_map = NULL,
    size = 3,
    color = 'black',
    fontface = 'plain',
    halo = TRUE,
    seed = 1,
    medoid_seed = seed,
    ...
){
    if(!requireNamespace('ggrepel', quietly = TRUE)) stop('labels need the ggrepel package.')
    centers <- .GroupCenters(x, y, group, match.arg(position), medoid_seed)
    centers$text <- .MapLabels(centers$group, label_map)
    args <- list(
        data = centers,
        mapping = aes(x = x, y = y, label = text),
        inherit.aes = FALSE,
        size = size,
        color = color,
        fontface = fontface,
        min.segment.length = 0.2,
        segment.size = 0.2,
        box.padding = 0.15,
        max.overlaps = Inf,
        seed = seed
    )
    if(halo) args <- c(args, list(bg.color = 'white', bg.r = 0.12))
    extra <- list(...)
    args[names(extra)] <- extra
    do.call(ggrepel::geom_text_repel, args)
}

# label positions per group (factor order kept); seed NULL uses the current random stream
.GroupCenters <- function(x, y, group, position = 'medoid', seed = 1){
    keep <- !is.na(group)
    group <- .AsFactor(group[keep])
    x <- x[keep]
    y <- y[keep]
    with_seed <- function(code){
        if(is.null(seed)) return(code)
        old <- if(exists('.Random.seed', envir = globalenv())) get('.Random.seed', envir = globalenv()) else NULL
        on.exit(if(!is.null(old)) assign('.Random.seed', old, envir = globalenv()))
        set.seed(seed)
        code
    }
    idx <- split(seq_along(group), group, drop = TRUE)
    centers <- with_seed(lapply(idx, function(i){
        if(position == 'centroid') return(c(mean(x[i]), mean(y[i])))
        .Medoid(x[i], y[i])
    }))
    data.frame(
        group = factor(names(idx), levels = levels(group)),
        x = sapply(centers, `[`, 1),
        y = sapply(centers, `[`, 2),
        row.names = NULL
    )
}

.MapLabels <- function(groups, label_map){
    groups <- as.character(groups)
    if(is.null(label_map)) return(groups)
    ifelse(groups %in% names(label_map), label_map[groups], groups)
}

# the cell closest on average to the others (from up to n_max sampled cells): a label position on the data
.Medoid <- function(x, y, n_max = 500){
    i <- if(length(x) > n_max) sample(length(x), n_max) else seq_along(x)
    d <- as.matrix(stats::dist(cbind(x[i], y[i])))
    j <- i[which.min(rowSums(d))]
    c(x[j], y[j])
}

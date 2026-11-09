#---------------------------------------------------------#
# Marker dot plot
#---------------------------------------------------------#

#' Per-group expression summary for a dot plot
#'
#' Mean expression and fraction of expressing cells (value > 0) of each feature in each group.
#'
#' @param object A Seurat object, or a data frame with one row per cell and one column per feature.
#' @param features Genes (or numeric columns).
#' @param group_by Grouping column (e.g. clusters).
#' @param layer,assay Assay layer and assay for a Seurat object (default: `data`, i.e. log-normalized expression).
#'
#' @return A data frame with `group`, `gene`, `mean_expr` and `frac_expr`.
#' @export
#' @examples
#' DotPlotData(sim_cells, features = c('CD3E', 'MS4A1', 'LYZ'), group_by = 'cluster')
DotPlotData <- function(object, features, group_by, layer = 'data', assay = NULL){
    df <- FetchCellData(object, vars = c(group_by, features), layer = layer, assay = assay)
    groups <- .AsFactor(df[[group_by]])
    do.call(rbind, lapply(features, function(f){
        v <- df[[f]]
        data.frame(
            group = factor(levels(groups), levels = levels(groups)),
            gene = f,
            mean_expr = as.numeric(tapply(v, groups, mean)),
            frac_expr = as.numeric(tapply(v > 0, groups, mean))
        )
    }))
}

#' Marker dot plot
#'
#' Genes on the y axis (italic), in facets by gene group; groups (e.g. clusters) on the x axis, optionally in facets
#' by a parent group with the strip under the panel. Color is the mean expression scaled to 0-1 per gene; dot size is
#' the fraction of cells expressing.
#'
#' @param data A summary from [DotPlotData()], or any data frame with group, gene, mean and fraction columns.
#' @param gene_groups Optional named list of genes per gene group, in display order (e.g. `list(T = c('CD3E',
#'   'CD8A'), B = 'MS4A1')`). Without it, genes keep their order in `data`, in one facet.
#' @param x_groups Optional named vector mapping each x group to a parent group (e.g. cluster to lineage); its order
#'   sets the facet order.
#' @param x,gene,expr,frac Column names in `data`.
#' @param scale_expr Scale the mean expression to 0-1 per gene (if FALSE, `expr` is plotted as is).
#' @param colors Low and high colors of the expression scale.
#' @param max_size Largest dot size.
#' @param x_label X axis title.
#'
#' @return A ggplot.
#' @export
#' @examples
#' dots <- DotPlotData(sim_cells, features = c('CD3E', 'CD8A', 'MS4A1', 'LYZ', 'CD14'), group_by = 'cluster')
#' MarkerDotPlot(dots,
#'               gene_groups = list(T = c('CD3E', 'CD8A'), B = 'MS4A1', Myeloid = c('LYZ', 'CD14')),
#'               x_groups = with(unique(sim_cells[c('cluster', 'lineage')]), setNames(lineage, cluster)))
MarkerDotPlot <- function(
    data,
    gene_groups = NULL,
    x_groups = NULL,
    x = 'group',
    gene = 'gene',
    expr = 'mean_expr',
    frac = 'frac_expr',
    scale_expr = TRUE,
    colors = c('grey92', '#08306b'),
    max_size = 3.2,
    x_label = NULL
){
    plot_df <- data.frame(x = data[[x]], gene = as.character(data[[gene]]), expr = data[[expr]],
                          frac = data[[frac]])
    if(scale_expr){
        rng <- tapply(plot_df$expr, plot_df$gene, range)
        lo <- sapply(rng, `[`, 1)[plot_df$gene]
        hi <- sapply(rng, `[`, 2)[plot_df$gene]
        plot_df$expr <- ifelse(hi > lo, (plot_df$expr - lo) / (hi - lo), 0)
    }

    if(is.null(gene_groups)) gene_groups <- list(all = unique(plot_df$gene))
    gene_map <- stats::setNames(rep(names(gene_groups), lengths(gene_groups)), unlist(gene_groups))
    plot_df <- plot_df[plot_df$gene %in% names(gene_map), ]
    plot_df$gene_group <- factor(gene_map[plot_df$gene], levels = names(gene_groups))
    plot_df$gene <- factor(plot_df$gene, levels = rev(unique(names(gene_map))))

    plot_df$x <- .AsFactor(plot_df$x)
    if(!is.null(x_groups)){
        plot_df$x_group <- factor(x_groups[as.character(plot_df$x)], levels = unique(x_groups))
        # x groups ordered by parent group, then their own order
        plot_df$x <- factor(plot_df$x, levels = levels(plot_df$x)[order(match(x_groups[levels(plot_df$x)],
                                                                                  unique(x_groups)))])
    }

    p <- ggplot(plot_df, aes(x = x, y = gene)) +
        geom_point(aes(size = frac, color = expr)) +
        scale_color_gradient(low = colors[1], high = colors[2], limits = if(scale_expr) c(0, 1),
                             breaks = if(scale_expr) c(0, 0.5, 1) else waiver(),
                             name = if(scale_expr) 'Scaled\nexpression' else 'Mean\nexpression') +
        scale_size_area(max_size = max_size, name = 'Fraction\nexpressing', breaks = c(0.25, 0.5, 0.75)) +
        labs(x = x_label, y = NULL) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(face = 'italic'),
            strip.background = element_blank(),
            strip.placement = 'outside',
            strip.text.x = element_text(angle = 90, hjust = 1),
            strip.text.y = element_text(angle = 0, hjust = 0),
            panel.spacing = grid::unit(0.15, 'lines'),
            panel.border = element_rect(color = 'grey70', fill = NA, linewidth = 0.3)
        )
    show_rows <- length(gene_groups) > 1
    if(!is.null(x_groups) && show_rows){
        p <- p + facet_grid(gene_group ~ x_group, scales = 'free', space = 'free', switch = 'x')
    } else if(!is.null(x_groups)){
        p <- p + facet_grid(. ~ x_group, scales = 'free', space = 'free', switch = 'x')
    } else if(show_rows){
        p <- p + facet_grid(gene_group ~ ., scales = 'free', space = 'free')
    }
    p
}

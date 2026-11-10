#---------------------------------------------------------#
# Differential expression: volcano and diverging count bars
#
# Generalized from PLEURA's eda/chemistry_supplementary_figure.Rmd. Both plots share one classification of the genes:
# significant (p < p_cutoff) or not, split by direction and by |log2(FC)| >= fc_cutoff, with the same colours (dark
# for strong, light for weak, blue down, red up). They take a plain results table (limma, DESeq2, edgeR, Seurat
# FindMarkers, ...) with a fold-change column and a p value column.
#---------------------------------------------------------#

#' Volcano plot
#'
#' Log2 fold change against -log10 p value, with the genes colored by direction and effect size (significant at
#' `p < p_cutoff`; strong at `|log2(FC)| >= fc_cutoff`), the cutoffs as dashed / dotted lines, and optional gene labels
#' (circled, with repelled text). Uses the same categories and colors as [DECountPlot()].
#'
#' @param de A data frame of differential expression results, one row per gene.
#' @param fc_col,p_col Effect size (usually log2 fold change) and p value (usually adjusted) columns. Defaults are
#'   limma's; e.g. DESeq2 `'log2FoldChange'`, `'padj'`; Seurat FindMarkers `'avg_log2FC'`, `'p_val_adj'`; edgeR
#'   `'logFC'`, `'FDR'`.
#' @param gene_col Gene name column; `NULL` (or a column that is not there, with the rownames set) uses the rownames,
#'   as in FindMarkers and DESeq2 results.
#' @param p_cutoff,fc_cutoff Significance and effect-size cutoffs. `fc_cutoff = NULL` splits by direction only.
#' @param direction_labels Legend text for the two directions (e.g. `c('Lower in GEM-X', 'Higher in GEM-X')`).
#' @param colors Named colors by category (`lower_big`, `lower_small`, `higher_small`, `higher_big`, `ns`; or `lower`,
#'   `higher`, `ns` without `fc_cutoff`), colors in that order, or a palette name (low to high). Default blue / red.
#' @param label Genes to label.
#' @param label_bold Genes to label in bold (e.g. the genes of interest, with `label` for the rest).
#' @param label_top Also label the `label_top` significant genes with the smallest p values.
#' @param symmetric Center the x axis on zero and make the panel square.
#' @param point_size,point_alpha Point size and transparency.
#' @param raster,raster_dpi Rasterize the points, and at which resolution.
#' @param text_size Label text size.
#' @param legend_ncol Legend columns.
#' @param x_label,y_label Axis titles (default log2(fold change), and -log10 of the p value, 'adj. P' when `p_col`
#'   looks adjusted).
#' @param title Plot title.
#'
#' @return A ggplot.
#' @export
#' @examples
#' de <- sim_de[sim_de$lineage == 'Myeloid', ]
#' VolcanoPlot(de, label = c('TECH1', 'TECH2', 'TECH3', 'TECH4'), label_bold = 'MARKD', label_top = 5,
#'             direction_labels = c('Lower in tumour', 'Higher in tumour'), title = 'Myeloid')
VolcanoPlot <- function(
    de,
    fc_col = 'logFC',
    p_col = 'adj.P.Val',
    gene_col = 'gene',
    p_cutoff = 0.05,
    fc_cutoff = 1,
    direction_labels = c('Down', 'Up'),
    colors = NULL,
    label = NULL,
    label_bold = NULL,
    label_top = 0,
    symmetric = FALSE,
    point_size = 0.4,
    point_alpha = 0.5,
    raster = TRUE,
    raster_dpi = 300,
    text_size = 2.8,
    legend_ncol = 2,
    x_label = expression(log[2] * '(fold change)'),
    y_label = .PLabel(p_col),
    title = NULL
){
    df <- .DETable(de, fc_col, p_col, gene_col, p_cutoff, fc_cutoff)
    df$y <- -log10(df$p)
    colors <- .DEColors(colors, levels(df$status))

    # labelled genes: given, bold, and the top significant ones
    if(is.null(df$gene) && (length(c(label, label_bold)) > 0 || label_top > 0)){
        stop('labels need gene names: set gene_col, or give the genes as rownames.')
    }
    sig <- df[df$status != 'ns', ]
    top <- if(label_top > 0) utils::head(sig$gene[order(sig$p, -abs(sig$fc))], label_top) else NULL
    label_df <- df[df$gene %in% c(label, label_bold, top), ]
    label_df$face <- ifelse(label_df$gene %in% label_bold, 'bold', 'plain')
    missing <- setdiff(c(label, label_bold), df$gene)
    if(length(missing) > 0) warning('genes to label not found: ', paste(missing, collapse = ', '))

    points <- geom_point(aes(color = status), size = point_size, alpha = point_alpha)
    p <- ggplot(df, aes(x = fc, y = y)) +
        (if(raster) ggrastr::rasterise(points, dpi = raster_dpi) else points) +
        geom_hline(yintercept = -log10(p_cutoff), linetype = 'dashed', color = 'grey50')
    if(!is.null(fc_cutoff)) p <- p + geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 'dotted',
                                                color = 'grey50')
    if(nrow(label_df) > 0){
        if(!requireNamespace('ggrepel', quietly = TRUE)) stop('gene labels need the ggrepel package.')
        p <- p +
            geom_point(data = label_df, shape = 21, fill = NA, color = 'black', size = 1.6) +
            ggrepel::geom_text_repel(data = label_df, aes(label = gene, fontface = face), size = text_size,
                                     min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4, seed = 1)
    }
    p <- p +
        scale_color_manual(values = colors, breaks = levels(df$status),
                           labels = .DELegendLabels(direction_labels, fc_cutoff)) +
        guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size = 2.5, alpha = 1))) +
        labs(x = x_label, y = y_label, title = title) +
        theme(legend.title = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7),
              plot.title = element_text(hjust = 0.5))
    if(symmetric){
        x_max <- ceiling(max(abs(df$fc)))
        p <- p +
            scale_x_continuous(limits = c(-x_max, x_max)) +
            theme(aspect.ratio = 1)
    }
    p
}

#' Diverging bars of differentially expressed genes per group
#'
#' One horizontal bar per group (e.g. lineage or cell type) with the number of significant genes, decreases to the
#' left of zero and increases to the right, each split into strong (`|log2(FC)| >= fc_cutoff`, dark, outside) and weak
#' (light, next to zero) segments, with the counts above the weak segments and at the ends of the bars. Uses the same categories and colors as
#' [VolcanoPlot()].
#'
#' @inheritParams VolcanoPlot
#' @param de A data frame of differential expression results, one row per gene and group.
#' @param group_by Group column; its factor levels (or sorted values) run from top to bottom.
#' @param bar_width Bar height (as a fraction of the space per group).
#' @param x_label X axis title.
#'
#' @return A ggplot.
#' @export
#' @examples
#' DECountPlot(sim_de, group_by = 'lineage', direction_labels = c('Lower in tumour', 'Higher in tumour'))
DECountPlot <- function(
    de,
    group_by,
    fc_col = 'logFC',
    p_col = 'adj.P.Val',
    p_cutoff = 0.05,
    fc_cutoff = 1,
    direction_labels = c('Down', 'Up'),
    colors = NULL,
    bar_width = 0.6,
    text_size = 2.3,
    legend_ncol = 2,
    x_label = sprintf('Significant genes (%s < %s)', .PName(p_col), p_cutoff)
){
    df <- .DETable(de, fc_col, p_col, NULL, p_cutoff, fc_cutoff)
    df$group <- .AsFactor(de[[group_by]])
    categories <- setdiff(levels(df$status), 'ns')
    colors <- .DEColors(colors, levels(df$status))

    # one row per group and category, with the segment ends: weak next to zero, strong outside
    tab <- table(df$group, factor(df$status, levels = categories))
    groups <- factor(rownames(tab), levels = levels(df$group))
    n_of <- function(cat) as.vector(tab[, cat])
    lo_in <- n_of(if(is.null(fc_cutoff)) 'lower' else 'lower_small')
    hi_in <- n_of(if(is.null(fc_cutoff)) 'higher' else 'higher_small')
    lo_out <- if(is.null(fc_cutoff)) rep(0, length(groups)) else n_of('lower_big')
    hi_out <- if(is.null(fc_cutoff)) rep(0, length(groups)) else n_of('higher_big')
    segs <- data.frame(
        group = rep(groups, 4),
        side = rep(c('lower', 'lower', 'higher', 'higher'), each = length(groups)),
        inner = rep(c(FALSE, TRUE, TRUE, FALSE), each = length(groups)),
        n = c(lo_out, lo_in, hi_in, hi_out),
        xmin = c(-(lo_in + lo_out), -lo_in, rep(0, length(groups)), hi_in),
        xmax = c(-lo_in, rep(0, length(groups)), hi_in, hi_in + hi_out)
    )
    segs$category <- if(is.null(fc_cutoff)) segs$side else
        paste0(segs$side, ifelse(segs$inner, '_small', '_big'))
    if(is.null(fc_cutoff)) segs <- segs[segs$inner, ]
    segs$category <- factor(segs$category, levels = categories)
    segs$y <- nlevels(df$group) + 1 - as.integer(factor(segs$group, levels = levels(df$group)))
    segs <- segs[segs$n > 0, ]

    # counts: above the weak segments, beside the zero line (clear of it); at the bar ends for the strong ones (or for
    # both sides without fc_cutoff)
    half <- bar_width / 2
    at_end <- if(is.null(fc_cutoff)) rep(TRUE, nrow(segs)) else !segs$inner
    above <- segs[!at_end, ]
    above$hjust <- ifelse(above$side == 'lower', 1.15, -0.15)
    ends <- segs[at_end, ]
    ends$x <- ifelse(ends$side == 'lower', ends$xmin, ends$xmax)
    ends$hjust <- ifelse(ends$side == 'lower', 1.25, -0.25)

    p <- ggplot(segs) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - half, ymax = y + half, fill = category)) +
        geom_vline(xintercept = 0, color = 'black') +
        geom_text(data = above, aes(x = 0, y = y + half, label = n, hjust = hjust), size = text_size, vjust = -0.3) +
        geom_text(data = ends, aes(x = x, y = y, label = n, hjust = hjust), size = text_size) +
        scale_fill_manual(values = colors, breaks = categories,
                          labels = .DELegendLabels(direction_labels, fc_cutoff)[seq_along(categories)],
                          drop = FALSE) +
        scale_y_continuous(breaks = seq_len(nlevels(df$group)), labels = rev(levels(df$group)),
                           limits = c(0.5, nlevels(df$group) + 0.5), expand = expansion(add = 0.1)) +
        scale_x_continuous(labels = abs, expand = expansion(mult = 0.1)) +
        labs(x = x_label, y = NULL) +
        guides(fill = guide_legend(ncol = legend_ncol)) +
        theme(legend.title = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7),
              axis.line.y = element_blank(), axis.ticks.y = element_blank())
    p
}

# rownames that are names, not the default 1..n
.hasRownames <- function(df){
    rn <- attr(df, 'row.names')
    is.character(rn) && !identical(rn, as.character(seq_len(nrow(df))))
}

# fold change, p value (zeros set to the smallest nonzero p) and category of each gene
.DETable <- function(de, fc_col, p_col, gene_col, p_cutoff, fc_cutoff){
    de <- as.data.frame(de)
    has_rownames <- .hasRownames(de)
    # genes from the rownames (FindMarkers, DESeq2) when there is no gene column
    use_rownames <- !is.null(gene_col) && !(gene_col %in% colnames(de)) && has_rownames
    missing <- setdiff(c(fc_col, p_col, if(!use_rownames) gene_col), colnames(de))
    if(length(missing) > 0){
        stop('columns not found: ', paste(missing, collapse = ', '), '. Set fc_col / p_col / gene_col to the names ',
             'in your table (columns: ', paste(utils::head(colnames(de), 12), collapse = ', '), ').')
    }
    df <- data.frame(fc = de[[fc_col]], p = de[[p_col]])
    if(use_rownames || (is.null(gene_col) && has_rownames)){
        df$gene <- rownames(de)
    } else if(!is.null(gene_col)){
        df$gene <- as.character(de[[gene_col]])
    }
    if(any(df$p == 0, na.rm = TRUE)) df$p[df$p == 0 & !is.na(df$p)] <- min(df$p[df$p > 0], na.rm = TRUE)
    df$status <- .DECategory(df$fc, df$p, p_cutoff, fc_cutoff)
    df
}

.DECategory <- function(fc, p, p_cutoff, fc_cutoff){
    sig <- !is.na(p) & !is.na(fc) & p < p_cutoff
    if(is.null(fc_cutoff)){
        status <- ifelse(!sig, 'ns', ifelse(fc < 0, 'lower', 'higher'))
        return(factor(status, levels = c('lower', 'higher', 'ns')))
    }
    big <- abs(fc) >= fc_cutoff
    status <- ifelse(!sig, 'ns', paste0(ifelse(fc < 0, 'lower', 'higher'), ifelse(big, '_big', '_small')))
    factor(status, levels = c('lower_big', 'lower_small', 'higher_small', 'higher_big', 'ns'))
}

.DEColors <- function(colors, categories){
    defaults <- c(lower_big = 'dodgerblue4', lower_small = 'lightskyblue', higher_small = 'lightpink',
                  higher_big = 'firebrick3', lower = 'dodgerblue4', higher = 'firebrick3', ns = 'grey85')
    if(is.null(colors)) return(defaults[categories])
    sig_categories <- setdiff(categories, 'ns')
    if(.IsPaletteName(colors)) colors <- PaletteColors(colors, n = length(sig_categories))
    if(is.null(names(colors))){
        if(!(length(colors) %in% c(length(sig_categories), length(categories)))){
            stop('give ', length(sig_categories), ' colors (', paste(sig_categories, collapse = ', '),
                 '), optionally with one for ns.')
        }
        names(colors) <- categories[seq_along(colors)]
    }
    out <- defaults[categories]
    out[intersect(names(colors), categories)] <- colors[intersect(names(colors), categories)]
    out
}

# legend labels (plotmath), e.g. "Down, |log2(FC)| >= 1"
.DELegendLabels <- function(direction_labels, fc_cutoff){
    if(length(direction_labels) != 2) stop('direction_labels needs two labels (down, up).')
    if(is.null(fc_cutoff)) return(c(direction_labels, 'Not significant'))
    lo <- direction_labels[1]
    hi <- direction_labels[2]
    labels <- list(
        bquote(.(lo) * ', |' * log[2] * '(FC)|' >= .(fc_cutoff)),
        bquote(.(lo) * ', |' * log[2] * '(FC)|' < .(fc_cutoff)),
        bquote(.(hi) * ', |' * log[2] * '(FC)|' < .(fc_cutoff)),
        bquote(.(hi) * ', |' * log[2] * '(FC)|' >= .(fc_cutoff)),
        'Not significant'
    )
    do.call(expression, labels)
}

# p value name for the axis titles
.PName <- function(p_col){
    if(grepl('adj|fdr|^q|qval', p_col, ignore.case = TRUE)) 'adj. P' else 'P'
}

.PLabel <- function(p_col){
    if(.PName(p_col) == 'adj. P') expression(-log[10] * '(adj. P)') else expression(-log[10] * '(P)')
}

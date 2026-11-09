#---------------------------------------------------------#
# QC summary
#
# Generalized from PLEURA's v2/overview/qc_summary_figure.Rmd: cells retained at each stage (a), QC metric
# distributions per library with the QC limits (b), cells per sample at each stage (c) and cell quality per sample (d,
# DistributionStack), with the sample columns of c and d aligned.
#
# Stage counts are a table with one row per sample and one column per stage (e.g. 'Cell Ranger', 'CellBender',
# 'Final'), plus any sample columns (label, disease, ...).
#---------------------------------------------------------#

#' Cells retained at each stage
#'
#' Stacked bars of the cells kept at each filtering stage, summed over samples and stacked by a sample group (e.g.
#' disease), with the totals on top (e.g. "357.3k").
#'
#' @param stage_counts A data frame with one row per sample and one column per stage.
#' @param stages The stage columns, in order; names, if given, are the stage labels.
#' @param fill_by Optional column of `stage_counts` to stack the bars by.
#' @param colors Named colors (or a palette name) for `fill_by`.
#' @param title Plot title.
#' @param fill_title Legend title (default `fill_by`).
#'
#' @return A ggplot.
#' @export
#' @examples
#' StageCountPlot(sim_samples, stages = c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final'),
#'                fill_by = 'disease')
StageCountPlot <- function(
    stage_counts,
    stages,
    fill_by = NULL,
    colors = NULL,
    title = 'Cells retained',
    fill_title = fill_by
){
    long <- .StageLong(stage_counts, stages)
    long$fill <- if(is.null(fill_by)) factor('all') else .AsFactor(long[[fill_by]])
    sums <- stats::aggregate(list(n = long$n), by = list(stage = long$stage, fill = long$fill), FUN = sum)
    totals <- stats::aggregate(list(n = long$n), by = list(stage = long$stage), FUN = sum)
    totals$text <- scales::label_number(scale = 1e-3, suffix = 'k', accuracy = 0.1)(totals$n)

    p <- ggplot(sums, aes(x = stage, y = n)) +
        geom_col(aes(fill = fill), width = 0.7, color = 'black', linewidth = 0.2) +
        geom_text(data = totals, aes(label = text), vjust = -0.4, size = 3.2) +
        scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = 'Cells', fill = fill_title, title = title) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom',
              legend.justification = 'center')
    if(is.null(fill_by)){
        p + scale_fill_manual(values = c(all = 'grey60'), guide = 'none')
    } else {
        p +
            scale_fill_manual(values = .AsNamedColors(colors, levels(sums$fill))) +
            guides(fill = guide_legend(nrow = 2, title.position = 'top', title.hjust = 0.5))
    }
}

#' Cells per sample at each stage
#'
#' Dodged bars of the cells kept per sample at each stage (a grey ramp from the first to the last stage), with the
#' samples optionally in facets (e.g. by disease).
#'
#' @inheritParams StageCountPlot
#' @param sample_col,label_col Sample column and optional column of x axis labels.
#' @param facet_by Optional column to facet the samples by, with free x space.
#' @param facet_labels Optional named vector of facet labels.
#' @param stage_colors Optional colors per stage (default a grey ramp, light to black).
#' @param show_x Show the sample labels.
#'
#' @return A ggplot.
#' @export
#' @examples
#' SampleStagePlot(sim_samples, stages = c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final'),
#'                 facet_by = 'disease')
SampleStagePlot <- function(
    stage_counts,
    stages,
    sample_col = 'sample',
    label_col = NULL,
    facet_by = NULL,
    facet_labels = NULL,
    stage_colors = NULL,
    show_x = TRUE,
    title = 'Cells per sample'
){
    long <- .StageLong(stage_counts, stages)
    long$label <- .SampleLabels(stage_counts, sample_col, label_col, facet_by)[long$.row]
    if(!is.null(facet_by)) long$facet_value <- .AsFactor(long[[facet_by]])
    if(is.null(stage_colors)){
        stage_colors <- grDevices::colorRampPalette(c('grey85', 'black'))(nlevels(long$stage))
    }
    stage_colors <- .AsNamedColors(stage_colors, levels(long$stage))
    if(is.null(names(stage_colors))) names(stage_colors) <- levels(long$stage)

    p <- ggplot(long, aes(x = label, y = n, fill = stage)) +
        geom_col(position = position_dodge(width = 0.85), width = 0.85) +
        scale_fill_manual(values = stage_colors) +
        scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.05))) +
        labs(x = NULL, y = 'Cells', fill = NULL, title = title) +
        theme(
            axis.text.x = if(show_x) element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7.5) else element_blank(),
            axis.ticks.x = if(show_x) element_line() else element_blank(),
            strip.background = element_blank(),
            panel.spacing.x = grid::unit(0.4, 'lines')
        )
    if(!is.null(facet_by)) p <- p + .SampleFacet(facet_labels)
    p
}

#' QC metric distributions per library
#'
#' One density curve per library (or sample) for each QC metric, colored by a library variable (e.g. chemistry),
#' with the QC limits as dashed lines. Each metric gets its own panel and x axis transform (e.g. log10 for UMIs and
#' genes, square root for percent mitochondrial).
#'
#' @param cells A data frame with one row per cell (or a Seurat object; columns from its metadata).
#' @param metrics Metric columns; names, if given, are the axis titles (e.g. `c('UMIs per cell' = 'nCount_RNA')`).
#' @param group_by Column with one curve per value (e.g. library).
#' @param color_by Optional column to color the curves by.
#' @param colors Named colors (or a palette name) for `color_by`.
#' @param thresholds Optional QC limits per metric: a named list or vector by metric column, one or more values each.
#' @param trans Optional x axis transform per metric (named by metric column): `'log10'`, `'sqrt'` or `'identity'`.
#'   Values that the transform cannot take (zero or below for log10) are left out.
#' @param breaks Optional x axis breaks per metric (named list by metric column).
#' @param color_title Legend title (default `color_by`).
#'
#' @return A patchwork of one panel per metric, stacked, with one legend.
#' @export
#' @examples
#' MetricDensityPlot(sim_cells, metrics = c('UMIs per cell' = 'n_umi', 'Mitochondrial (%)' = 'pct_mt'),
#'                   group_by = 'sample', color_by = 'disease', thresholds = list(n_umi = 500, pct_mt = 20),
#'                   trans = c(n_umi = 'log10', pct_mt = 'sqrt'))
MetricDensityPlot <- function(
    cells,
    metrics,
    group_by,
    color_by = NULL,
    colors = NULL,
    thresholds = NULL,
    trans = NULL,
    breaks = NULL,
    color_title = color_by
){
    if(is.null(names(metrics))) names(metrics) <- metrics
    names(metrics)[names(metrics) == ''] <- metrics[names(metrics) == '']
    df <- FetchCellData(cells, vars = c(group_by, color_by, unname(metrics)))
    df$group <- as.character(df[[group_by]])
    df$color <- if(is.null(color_by)) factor('all') else .AsFactor(df[[color_by]])
    colors <- if(is.null(color_by)) c(all = 'grey30') else .AsNamedColors(colors, levels(df$color))

    panels <- lapply(seq_along(metrics), function(i){
        m <- metrics[[i]]
        cur_trans <- if(!is.null(trans) && m %in% names(trans)) trans[[m]] else 'identity'
        cur <- data.frame(group = df$group, color = df$color, value = df[[m]])
        cur <- cur[!is.na(cur$value), ]
        if(cur_trans == 'log10') cur <- cur[cur$value > 0, ]
        if(cur_trans == 'sqrt') cur <- cur[cur$value >= 0, ]
        cur_breaks <- if(!is.null(breaks) && m %in% names(breaks)) breaks[[m]] else waiver()
        p <- ggplot(cur, aes(x = value, group = group, color = color)) +
            geom_density(linewidth = 0.3, alpha = 0.8) +
            scale_x_continuous(transform = cur_trans, breaks = cur_breaks,
                               labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
            scale_color_manual(values = colors, guide = if(is.null(color_by) || i > 1) 'none' else 'legend',
                               name = color_title) +
            guides(color = if(is.null(color_by) || i > 1) 'none' else
                       guide_legend(override.aes = list(linewidth = 1))) +
            labs(x = names(metrics)[i], y = 'Density')
        if(!is.null(thresholds) && m %in% names(thresholds)){
            p <- p + geom_vline(xintercept = thresholds[[m]], linetype = 'dashed', color = 'grey30')
        }
        p
    })
    patchwork::wrap_plots(panels, ncol = 1, guides = 'collect')
}

#' QC summary figure
#'
#' Four panels: (a) cells retained at each stage ([StageCountPlot()]), (b) QC metric distributions per library with
#' the QC limits ([MetricDensityPlot()]), (c) cells per sample at each stage ([SampleStagePlot()]) and (d) cell
#' quality per sample (quasirandom points under boxplots, [DistributionStack()]). The sample columns of c and d line
#' up.
#'
#' Write it with `cairo_pdf()` (the default stage labels and limits may use non-ASCII characters); about 15 x 16
#' inches suits 50 samples.
#'
#' @param cells A data frame with one row per cell, with the sample, library, QC metric and color columns.
#' @param stage_counts A data frame with one row per sample and one column per stage, plus `sample_col` and the
#'   optional label and group columns; its row order (within groups) sets the sample order.
#' @param stages The stage columns of `stage_counts`, in order; names, if given, are the stage labels.
#' @param metrics Metric columns of `cells` for panel b; names are the axis titles.
#' @param quality_metrics Metric columns of `cells` for panel d (default `metrics`); names are the row labels (e.g.
#'   `c('UMIs (log10)' = 'log10_umi')`).
#' @param sample_col,label_col Sample column (in both tables) and optional label column of `stage_counts`.
#' @param library_col Column of `cells` with one density curve per value (default `sample_col`).
#' @param group_col Optional column of `stage_counts` grouping the samples into facets (e.g. disease); also stacks
#'   panel a.
#' @param group_colors Named colors (or a palette name) for `group_col` in panel a.
#' @param group_labels Optional named vector of facet labels for c and d (e.g. `c('Lung cancer' = 'LC')`).
#' @param color_by Optional column of `cells` coloring the curves (b) and points (d), e.g. chemistry.
#' @param colors Named colors (or a palette name) for `color_by`.
#' @param thresholds,trans,breaks Passed to [MetricDensityPlot()].
#' @param pass_col Optional logical column of `cells`: panel d shows only the cells passing (b shows all).
#' @param n_points Points per sample and metric in panel d.
#' @param raster_dpi Resolution of the rasterized points in panel d.
#' @param titles Titles of panels a, c and d.
#' @param group_title,color_title Legend titles for `group_col` and `color_by` (default the column names; e.g.
#'   'Disease', 'Chemistry').
#'
#' @return A patchwork.
#' @export
#' @examples
#' QCSummary(
#'     sim_cells, sim_samples,
#'     stages = c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final'),
#'     metrics = c('UMIs per cell' = 'n_umi', 'Genes per cell' = 'n_genes', 'Mitochondrial (%)' = 'pct_mt'),
#'     group_col = 'disease', color_by = 'disease',
#'     thresholds = list(n_umi = 500, n_genes = 250, pct_mt = 20),
#'     trans = c(n_umi = 'log10', n_genes = 'log10', pct_mt = 'sqrt')
#' )
QCSummary <- function(
    cells,
    stage_counts,
    stages,
    metrics,
    quality_metrics = metrics,
    sample_col = 'sample',
    label_col = NULL,
    library_col = sample_col,
    group_col = NULL,
    group_colors = NULL,
    group_labels = NULL,
    color_by = NULL,
    colors = NULL,
    thresholds = NULL,
    trans = NULL,
    breaks = NULL,
    pass_col = NULL,
    n_points = 300,
    raster_dpi = 600,
    titles = c('Cells retained', 'Cells per sample', 'Cell quality per sample'),
    group_title = group_col,
    color_title = color_by
){
    stage_counts <- as.data.frame(stage_counts, check.names = FALSE)
    cells <- as.data.frame(cells, check.names = FALSE)
    labels <- .SampleLabels(stage_counts, sample_col, label_col, group_col)
    # one color map when the samples are grouped and colored by the same variable (e.g. disease)
    if(!is.null(group_col) && identical(group_col, color_by)){
        lv <- levels(.AsFactor(stage_counts[[group_col]]))
        if(is.null(colors)) colors <- if(is.null(group_colors)) .AsNamedColors(NULL, lv) else group_colors
        if(is.null(group_colors)) group_colors <- .AsNamedColors(colors, lv)
        cells[[color_by]] <- factor(as.character(cells[[color_by]]), levels = lv)
    }
    missing <- setdiff(unique(as.character(cells[[sample_col]])), stage_counts[[sample_col]])
    if(length(missing) > 0) stop('samples in cells but not in stage_counts: ', paste(missing, collapse = ', '))

    p_a <- StageCountPlot(stage_counts, stages, fill_by = group_col, colors = group_colors, title = titles[1],
                          fill_title = group_title)
    p_b <- MetricDensityPlot(cells, metrics, group_by = library_col, color_by = color_by, colors = colors,
                             thresholds = thresholds, trans = trans, breaks = breaks, color_title = color_title)
    p_c <- SampleStagePlot(stage_counts, stages, sample_col = sample_col, label_col = label_col,
                           facet_by = group_col, facet_labels = group_labels, show_x = FALSE, title = titles[2])

    # panel d: per-sample distributions, columns as in panel c
    quality <- if(is.null(pass_col)) cells else cells[cells[[pass_col]] %in% TRUE, ]
    quality$.label <- labels[match(as.character(quality[[sample_col]]), stage_counts[[sample_col]])]
    if(!is.null(group_col)){
        group <- .AsFactor(stage_counts[[group_col]])
        quality$.group <- group[match(as.character(quality[[sample_col]]), stage_counts[[sample_col]])]
        quality$.group <- factor(.MapLabels(quality$.group, group_labels),
                                 levels = .MapLabels(levels(group), group_labels))
    }
    p_d <- DistributionStack(quality, group_by = '.label', metrics = quality_metrics, color_by = color_by,
                             colors = colors, facet_by = if(is.null(group_col)) NULL else '.group',
                             count_bar = FALSE, n_points = n_points, point_size = 0.12, box_width = 0.65,
                             box_linewidth = 0.25, raster_dpi = raster_dpi, x_label = NULL,
                             show_legend = !is.null(color_by)) +
        labs(title = titles[3], color = color_title) +
        theme(
            strip.text.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7.5),
            strip.text.y.left = element_text(angle = 90, hjust = 0.5),
            panel.spacing.x = grid::unit(0.4, 'lines')
        )

    top <- (p_a + labs(tag = 'a') | patchwork::wrap_elements(full = p_b) + labs(tag = 'b')) +
        patchwork::plot_layout(widths = c(1, 1.6))
    bottom <- (p_c + labs(tag = 'c')) / (p_d + labs(tag = 'd')) +
        patchwork::plot_layout(heights = c(1.1, 3), guides = 'collect')
    patchwork::wrap_elements(top) / patchwork::wrap_elements(bottom) +
        patchwork::plot_layout(heights = c(1, 2.1))
}

# stage counts in long form: one row per sample and stage (.row = the sample's row)
.StageLong <- function(stage_counts, stages){
    stage_counts <- as.data.frame(stage_counts, check.names = FALSE)
    if(is.null(names(stages))) names(stages) <- stages
    names(stages)[names(stages) == ''] <- stages[names(stages) == '']
    missing <- setdiff(stages, colnames(stage_counts))
    if(length(missing) > 0) stop('stage columns not found: ', paste(missing, collapse = ', '))
    long <- do.call(rbind, lapply(names(stages), function(s){
        cur <- stage_counts
        cur$.row <- seq_len(nrow(cur))
        cur$stage <- s
        cur$n <- stage_counts[[stages[[s]]]]
        cur
    }))
    long$stage <- factor(long$stage, levels = names(stages))
    long
}

# sample labels as a factor in row order within groups
.SampleLabels <- function(stage_counts, sample_col, label_col, group_col){
    labels <- as.character(stage_counts[[if(is.null(label_col)) sample_col else label_col]])
    if(anyDuplicated(labels)) stop('sample labels must be unique.')
    ord <- if(is.null(group_col)) seq_along(labels) else order(.AsFactor(stage_counts[[group_col]]), seq_along(labels))
    factor(labels, levels = labels[ord])
}

.SampleFacet <- function(facet_labels){
    labeller <- if(is.null(facet_labels)) 'label_value' else as_labeller(function(x) .MapLabels(x, facet_labels))
    facet_grid(cols = vars(facet_value), scales = 'free_x', space = 'free_x', labeller = labeller)
}

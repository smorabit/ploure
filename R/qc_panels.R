#---------------------------------------------------------#
# QC panels
#
# Generalized from PLEURA's v2/overview/qc_summary_figure.Rmd: cells kept at each filtering stage, in total or per
# sample (StageCountPlot), and QC metric distributions per library with the QC limits (MetricDensityPlot). Each is
# one panel; a project's QC summary figure combines them (with e.g. DistributionStack) in its own layout, as
# validation/pleura_qc_summary.R does for PLEURA.
#
# Stage counts are a table with one row per sample and one column per stage (e.g. 'Cell Ranger', 'CellBender',
# 'Final'), plus any sample columns (label, disease, ...).
#---------------------------------------------------------#

#' Cells kept at each filtering stage
#'
#' Bar plot of the cells kept at each stage of a QC pipeline, from a table with one row per sample and one column per
#' stage. Two views:
#'
#' * `by_sample = FALSE` (default): one bar per stage, summed over samples, optionally stacked by a sample group
#'   (`fill_by`, e.g. disease), with the totals on top (e.g. "357.3k").
#' * `by_sample = TRUE`: one group of bars per sample, one bar per stage (a grey ramp from the first to the last
#'   stage by default). Samples are ordered by `facet_by`, then by their row order in `stage_counts`.
#'
#' @param stage_counts A data frame with one row per sample and one column per stage.
#' @param stages The stage columns, in order; names, if given, are the stage labels.
#' @param by_sample One group of bars per sample instead of one bar per stage.
#' @param sample_col,label_col Sample column and optional column of x axis labels (`by_sample = TRUE`).
#' @param fill_by Optional column to stack the stage totals by (`by_sample = FALSE`).
#' @param facet_by Optional column to split the bars into column facets, with free x space.
#' @param facet_labels Optional named vector of facet labels (e.g. `c('LC' = 'Lung cancer')`).
#' @param colors Fill colors (named colors or a palette name): of the `fill_by` levels when `by_sample = FALSE`, of
#'   the stages when `by_sample = TRUE`.
#' @param show_totals Print the stage totals on top of the bars (`by_sample = FALSE`).
#' @param show_x Show the x axis labels.
#' @param title Plot title (default 'Cells retained', or 'Cells per sample' when `by_sample = TRUE`).
#' @param fill_title Legend title (default `fill_by`; none for the stages).
#'
#' @return A ggplot.
#' @export
#' @examples
#' stages <- c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final')
#' StageCountPlot(sim_samples, stages, fill_by = 'disease')
#' StageCountPlot(sim_samples, stages, by_sample = TRUE, facet_by = 'disease')
StageCountPlot <- function(
    stage_counts,
    stages,
    by_sample = FALSE,
    sample_col = 'sample',
    label_col = NULL,
    fill_by = NULL,
    facet_by = NULL,
    facet_labels = NULL,
    colors = NULL,
    show_totals = !by_sample,
    show_x = TRUE,
    title = if(by_sample) 'Cells per sample' else 'Cells retained',
    fill_title = fill_by
){
    if(by_sample && !is.null(fill_by)) warning('fill_by is ignored with by_sample = TRUE (bars are filled by stage).')
    long <- .StageLong(stage_counts, stages)
    long$facet_value <- if(is.null(facet_by)) factor('all') else .AsFactor(long[[facet_by]])

    if(by_sample){
        long$x <- .SampleLabels(stage_counts, sample_col, label_col, facet_by)[long$.row]
        long$fill <- long$stage
        if(is.null(colors)) colors <- grDevices::colorRampPalette(c('grey85', 'black'))(nlevels(long$stage))
        colors <- .AsNamedColors(colors, levels(long$stage))
        if(is.null(names(colors))) names(colors) <- levels(long$stage)
        bars <- geom_col(position = position_dodge(width = 0.85), width = 0.85)
        fill_title <- NULL
    } else {
        long$x <- long$stage
        long$fill <- if(is.null(fill_by)) factor('all') else .AsFactor(long[[fill_by]])
        colors <- if(is.null(fill_by)) c(all = 'grey60') else .AsNamedColors(colors, levels(long$fill))
        bars <- geom_col(width = 0.7, color = 'black', linewidth = 0.2)
    }
    sums <- stats::aggregate(list(n = long$n), by = list(x = long$x, fill = long$fill, facet_value = long$facet_value),
                             FUN = sum)

    p <- ggplot(sums, aes(x = x, y = n, fill = fill)) +
        bars +
        scale_y_continuous(labels = scales::label_comma(),
                           expand = expansion(mult = c(0, if(show_totals) 0.12 else 0.05))) +
        labs(x = NULL, y = 'Cells', fill = fill_title, title = title) +
        theme(
            axis.text.x = if(!show_x) element_blank() else if(by_sample)
                element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7.5) else element_text(angle = 45, hjust = 1),
            axis.ticks.x = if(show_x) element_line() else element_blank(),
            strip.background = element_blank(),
            panel.spacing.x = grid::unit(0.4, 'lines')
        )
    if(show_totals){
        totals <- stats::aggregate(list(n = sums$n), by = list(x = sums$x, facet_value = sums$facet_value), FUN = sum)
        totals$text <- scales::label_number(scale = 1e-3, suffix = 'k', accuracy = 0.1)(totals$n)
        p <- p + geom_text(data = totals, aes(x = x, y = n, label = text), inherit.aes = FALSE,
                           vjust = -0.4, size = 3.2)
    }
    if(!by_sample && is.null(fill_by)){
        p <- p + scale_fill_manual(values = colors, guide = 'none')
    } else {
        p <- p + scale_fill_manual(values = colors)
        if(!by_sample) p <- p +
            theme(legend.position = 'bottom', legend.justification = 'center') +
            guides(fill = guide_legend(nrow = 2, title.position = 'top', title.hjust = 0.5))
    }
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
            scale_color_manual(values = colors, name = color_title) +
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

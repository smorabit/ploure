#---------------------------------------------------------#
# Cohort matrix
#
# Generalized from PLEURA's v2/overview/cohort_matrix.Rmd: one column per sample, grouped into facets (e.g. by
# disease), with an optional bar plot on top (e.g. cells per sample) and blocks of tile rows below, each block with
# its own color scale and legend. Panels are sized in absolute units so the tiles stay square (coord_fixed() cannot
# be combined with free-width facets).
#---------------------------------------------------------#

#' Cohort matrix
#'
#' One column per sample, grouped into facets (e.g. by disease): an optional bar plot on top (e.g. cells per sample,
#' stacked by lineage) and blocks of tile rows below (sample design, data available, clinical variables). Each block
#' has its own color scale and legend: categorical (named colors) or, when all its variables are numeric, a
#' continuous gradient. Tiles are square, sized in absolute units, and the legends sit in columns to the right.
#'
#' @param samples A data frame with one row per sample, in display order within each group, holding the tile
#'   variables.
#' @param blocks A list of tile blocks, top to bottom. Each block is a list with:
#'   * `vars`: the columns of `samples` to draw, one tile row each; names, if given, are the row labels
#'     (e.g. `c('Chemistry' = 'chemistry', 'Site' = 'site')`);
#'   * `colors`: for a categorical block, named colors (their order is the legend order; values without a color are
#'     drawn as missing) or a palette name ([PaletteColors()]); for a continuous block, the colors from low to high
#'     (default grey to dark blue) or a palette name;
#'   * `title`: the legend title (default: the first row label);
#'   * `limits`: for a continuous block, the color limits.
#' @param sample_col Column of `samples` with the sample names (matched against `bar_data`).
#' @param group_col Optional column of `samples` to group the columns into facets; its factor levels (or order of
#'   appearance) set the facet order.
#' @param label_col Optional column with the column labels (default `sample_col`).
#' @param group_labels Optional named vector of facet labels.
#' @param group_counts Add the number of samples to the facet labels (e.g. "LC (28)").
#' @param bar_data Optional long data frame for the bar plot: `sample_col`, `bar_value` and optionally `bar_fill`.
#' @param bar_value,bar_fill Columns of `bar_data` with the bar heights and the stacking groups.
#' @param bar_colors Named colors for `bar_fill`.
#' @param bar_title,bar_label Bar legend title and y axis label.
#' @param tile_cm Tile size (cm); the matrix is one tile wide per sample and one tile high per row.
#' @param gap_cm Gap between facets (cm).
#' @param bar_cm Height of the bar plot (cm).
#' @param legend_ncol Number of legend columns.
#' @param legend_cm Width of a legend column (cm).
#' @param na_color Color of missing values.
#'
#' @return A patchwork. The attribute `pdf_size` holds a fitting page size in inches (width, height), e.g.
#'   `cairo_pdf(file, width = attr(p, 'pdf_size')[1], height = attr(p, 'pdf_size')[2])`; the attribute `panels`
#'   holds the individual ggplots (`bars`, `block1`, `block2`, ...), with their legends.
#' @export
#' @examples
#' cells <- as.data.frame(table(sample = sim_cells$sample, lineage = sim_cells$lineage), responseName = 'n')
#' p <- CohortMatrix(
#'     sim_samples, group_col = 'disease', bar_data = cells, bar_fill = 'lineage',
#'     bar_colors = c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c', Stromal = '#33a02c'), bar_title = 'Lineage',
#'     blocks = list(
#'         list(vars = c('Chemistry' = 'chemistry', 'Site' = 'site'), title = 'Chemistry / site',
#'              colors = c("5' v2" = '#1b9e77', 'GEM-X' = '#d95f02', 'Site A' = '#7570b3', 'Site B' = '#e7298a')),
#'         list(vars = c('Loaded CD45+' = 'pct_cd45', 'Loaded EPCAM+' = 'pct_epcam'), title = 'Cells loaded (%)',
#'              limits = c(0, 100)),
#'         list(vars = c('Sex' = 'sex', 'Age' = 'age_group'), title = 'Sex / age',
#'              colors = c(Female = '#c51b7d', Male = '#4d9221', '< 60' = '#fdd49e', '60-70' = '#fc8d59',
#'                         '> 70' = '#b30000')),
#'         list(vars = c('Bulk RNA-seq' = 'bulk_rnaseq'), title = 'Data available',
#'              colors = c(Yes = 'grey35', No = 'grey92'))
#'     )
#' )
#' p
CohortMatrix <- function(
    samples,
    blocks,
    sample_col = 'sample',
    group_col = NULL,
    label_col = NULL,
    group_labels = NULL,
    group_counts = TRUE,
    bar_data = NULL,
    bar_value = 'n',
    bar_fill = NULL,
    bar_colors = NULL,
    bar_title = NULL,
    bar_label = 'Cells',
    tile_cm = 0.45,
    gap_cm = 0.15,
    bar_cm = 4,
    legend_ncol = 2,
    legend_cm = 4.5,
    na_color = 'grey97'
){
    samples <- as.data.frame(samples, check.names = FALSE)
    if(!(sample_col %in% colnames(samples))) stop("samples needs a '", sample_col, "' column.")
    if(anyDuplicated(samples[[sample_col]])) stop('samples must have one row per sample.')
    if(is.null(label_col)) label_col <- sample_col
    blocks <- lapply(blocks, .CohortBlockSpec)
    missing <- setdiff(unlist(lapply(blocks, function(b) unname(b$vars))), colnames(samples))
    if(length(missing) > 0) stop('columns not found in samples: ', paste(missing, collapse = ', '))

    # columns: sample labels in row order, ordered by group
    group <- if(is.null(group_col)) factor(rep('all', nrow(samples))) else
        factor(as.character(samples[[group_col]]), levels = .GroupLevels(samples[[group_col]]))
    samples <- samples[order(group), , drop = FALSE]
    group <- group[order(group)]
    labels <- as.character(samples[[label_col]])
    meta <- data.frame(sample = as.character(samples[[sample_col]]), label = factor(labels, levels = labels),
                       group = group)
    strips <- .MapLabels(levels(group), group_labels)
    if(group_counts) strips <- paste0(strips, ' (', as.numeric(table(group)), ')')
    names(strips) <- levels(group)
    show_strip <- !is.null(group_col)

    #---------------------------------------------------------#
    # panels
    #---------------------------------------------------------#

    facet <- facet_grid(. ~ group, scales = 'free_x', space = 'free_x', labeller = as_labeller(strips))
    # sizes follow the current theme's base size (11 for the default themes)
    base_size <- ThemeTextSize(element = 'text', unit = 'pt')
    base_family <- ggplot2::calc_element('text', theme_get())$family
    base <- theme_classic(base_size = base_size, base_family = if(is.null(base_family)) '' else base_family) +
        theme(
            text = element_text(color = 'black'),
            axis.title.x = element_blank(),
            axis.text = element_text(color = 'black'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.spacing.x = grid::unit(gap_cm, 'cm'),
            legend.title = element_text(size = base_size - 1),
            legend.text = element_text(size = base_size - 2),
            legend.key.size = grid::unit(0.35, 'cm'),
            legend.frame = .ColorbarFrame(),
            legend.justification = c(0, 1),
            plot.margin = margin(1, 5, 1, 5)
        )

    panels <- list()
    if(!is.null(bar_data)){
        bars <- as.data.frame(bar_data)
        bars$sample <- as.character(bars[[sample_col]])
        bars <- bars[bars$sample %in% meta$sample, ]
        bars$label <- meta$label[match(bars$sample, meta$sample)]
        bars$group <- meta$group[match(bars$sample, meta$sample)]
        bars$value <- bars[[bar_value]]
        p <- if(is.null(bar_fill)){
            ggplot(bars, aes(x = label, y = value)) +
                geom_col(width = 0.85, fill = 'grey40')
        } else {
            bars$fill <- .AsFactor(bars[[bar_fill]])
            ggplot(bars, aes(x = label, y = value, fill = fill)) +
                geom_col(width = 0.85) +
                scale_fill_manual(values = .AsNamedColors(bar_colors, levels(bars$fill)),
                                  name = if(is.null(bar_title)) bar_fill else bar_title)
        }
        panels$bars <- p +
            facet +
            scale_y_continuous(labels = scales::label_comma(), expand = c(0, 0)) +
            scale_x_discrete(expand = c(0, 0)) +
            labs(y = bar_label) +
            base
    }
    for(i in seq_along(blocks)){
        panels[[paste0('block', i)]] <- .CohortBlock(samples, meta, blocks[[i]], na_color) +
            facet +
            base +
            theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank())
    }

    # strips on the top panel, sample labels under the bottom one
    if(show_strip) panels[[1]] <- panels[[1]] + theme(strip.text = element_text(size = base_size))
    n_panels <- length(panels)
    panels[[n_panels]] <- panels[[n_panels]] +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size - 3))

    #---------------------------------------------------------#
    # assemble: absolute panel sizes, legends in columns to the right
    #---------------------------------------------------------#

    n_rows <- sapply(blocks, function(b) length(b$vars))
    heights <- c(if(!is.null(bar_data)) bar_cm, n_rows * tile_cm)
    panel_w <- nrow(meta) * tile_cm + (nlevels(group) - 1) * gap_cm

    legends <- lapply(panels, function(p) ExtractLegend(p))
    n_keys <- sapply(panels, .LegendKeys)
    has_legend <- !sapply(legends, inherits, 'null')
    legends <- legends[has_legend]
    n_keys <- n_keys[has_legend]
    # fill the legend columns in order, splitting at about equal numbers of keys
    col_of <- pmin(legend_ncol, floor((cumsum(n_keys + 2) - (n_keys + 2) / 2) / sum(n_keys + 2) * legend_ncol) + 1)
    legend_cols <- lapply(seq_len(legend_ncol), function(k){
        cur <- which(col_of == k)
        if(length(cur) == 0) return(patchwork::plot_spacer())
        patchwork::wrap_plots(c(lapply(legends[cur], patchwork::wrap_elements), list(patchwork::plot_spacer())),
                              ncol = 1, heights = c(n_keys[cur] + 2, max(1, sum(n_keys + 2) / legend_ncol -
                                                                              sum(n_keys[cur] + 2))))
    })

    stripped <- lapply(panels, function(p) p + theme(legend.position = 'none'))
    matrix <- patchwork::wrap_plots(stripped, ncol = 1) +
        patchwork::plot_layout(heights = grid::unit(heights, 'cm'), widths = grid::unit(panel_w, 'cm'))
    label_cm <- 0.25 * max(nchar(labels))
    out <- patchwork::wrap_elements(full = matrix) +
        patchwork::wrap_plots(legend_cols, ncol = legend_ncol) +
        patchwork::plot_layout(widths = c(panel_w + 5, legend_ncol * legend_cm))
    attr(out, 'pdf_size') <- c(
        width = (panel_w + 5 + legend_ncol * legend_cm + 1) / 2.54,
        height = (sum(heights) + label_cm + 2.5) / 2.54
    )
    attr(out, 'panels') <- panels
    out
}

# a block spec with vars (named), colors, title, limits
.CohortBlockSpec <- function(block){
    if(!is.list(block) || is.null(block$vars)) stop('each block needs vars (a list with vars, colors, title).')
    vars <- block$vars
    if(is.null(names(vars))) names(vars) <- vars
    names(vars)[names(vars) == ''] <- vars[names(vars) == '']
    block$vars <- vars
    if(is.null(block$title)) block$title <- names(vars)[1]
    block
}

.GroupLevels <- function(values){
    if(is.factor(values)) return(levels(droplevels(values)))
    unique(as.character(values))
}

# one block of tile rows with its own scale
.CohortBlock <- function(samples, meta, block, na_color){
    vars <- block$vars
    continuous <- all(sapply(samples[unname(vars)], is.numeric))
    plot_df <- do.call(rbind, lapply(names(vars), function(v){
        data.frame(label = meta$label, group = meta$group, row = v,
                   value = if(continuous) samples[[vars[[v]]]] else as.character(samples[[vars[[v]]]]))
    }))
    plot_df$row <- factor(plot_df$row, levels = rev(names(vars)))
    p <- ggplot(plot_df, aes(x = label, y = row, fill = value)) +
        geom_tile(color = 'white', linewidth = 0.4, show.legend = TRUE) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0))
    if(continuous){
        colors <- if(is.null(block$colors)) c('grey97', '#08306b') else block$colors
        if(.IsPaletteName(colors)) colors <- PaletteColors(colors, 100, type = 'continuous')
        p + scale_fill_gradientn(colors = colors, limits = block$limits, name = block$title, na.value = na_color)
    } else {
        values <- sort(unique(stats::na.omit(plot_df$value)))
        colors <- .AsNamedColors(block$colors, values)
        p + scale_fill_manual(values = colors, limits = names(colors), name = block$title, na.value = na_color,
                              na.translate = FALSE)
    }
}

# number of legend keys of a panel (a color bar counts as five)
.LegendKeys <- function(p){
    b <- ggplot_build(p)
    sc <- b$plot$scales$get_scales('fill')
    if(is.null(sc)) return(0)
    if(sc$is_discrete()) return(length(sc$get_breaks()))
    5
}

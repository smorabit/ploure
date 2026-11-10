# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_de_plots.R

#---------------------------------------------------------#
# The PLEURA chemistry volcano (myeloid) and chemistry-sensitive genes per lineage
# (eda/chemistry_supplementary_figure.Rmd) drawn with VolcanoPlot and DECountPlot
#
# Draws both with the notebook's code and with ploure, compares the built layers, and writes
# validation/figures/de_plots_ploure.pdf (notebook panels on the left, ploure on the right).
#---------------------------------------------------------#

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrastr)
library(ggrepel)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

theme_set(theme_cowplot())
theme_update(plot.title = element_text(face = 'plain'))
setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

adj_p <- 0.05
lineages <- c('myeloid', 'CD4', 'CD8_NK', 'B', 'stromal', 'mesothelial', 'malignant')

#---------------------------------------------------------#
# the notebook's code, verbatim
#---------------------------------------------------------#

fc_colors <- c(lower_big = 'dodgerblue4', lower_small = 'lightskyblue', higher_small = 'lightpink',
               higher_big = 'firebrick3', ns = 'grey85')
fc_labels <- c(
    lower_big = expression('Lower in GEM-X, |' * log[2] * '(FC)|' >= 1),
    lower_small = expression('Lower in GEM-X, |' * log[2] * '(FC)|' < 1),
    higher_small = expression('Higher in GEM-X, |' * log[2] * '(FC)|' < 1),
    higher_big = expression('Higher in GEM-X, |' * log[2] * '(FC)|' >= 1),
    ns = expression('Not significant')
)
fc_category <- function(logFC, sig){
    case_when(
        !sig ~ 'ns',
        logFC < 0 & abs(logFC) >= 1 ~ 'lower_big',
        logFC < 0 ~ 'lower_small',
        abs(logFC) >= 1 ~ 'higher_big',
        TRUE ~ 'higher_small'
    )
}

chem <- read.csv('analysis/eda/pseudobulk/data/chemistry_limma_per_lineage.csv') %>%
    mutate(lineage = factor(lineage, levels = lineages), sig = adj.P.Val < adj_p)

label_technical <- c('MIF', 'RPS28', 'MT2A', 'CRIP1', 'CYBA')
label_biology <- c('APOE', 'C1QC', 'ISG15', 'CD74', 'HLA-DRA', 'SPP1', 'CXCL8')

volcano_df <- chem %>%
    filter(lineage == 'myeloid') %>%
    mutate(
        status = factor(fc_category(logFC, sig), levels = names(fc_colors)),
        label_group = case_when(
            gene %in% label_technical ~ 'shared technical',
            gene %in% label_biology ~ 'lineage biology',
            TRUE ~ NA_character_
        )
    )
label_df <- filter(volcano_df, !is.na(label_group))
v_orig <- volcano_df %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
    rasterise(geom_point(aes(color = status), size = 0.4, alpha = 0.5), dpi = 300) +
    geom_hline(yintercept = -log10(adj_p), linetype = 'dashed', color = 'grey50') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dotted', color = 'grey50') +
    geom_point(data = label_df, shape = 21, fill = NA, color = 'black', size = 1.6) +
    geom_text_repel(
        data = label_df,
        aes(label = gene, fontface = ifelse(label_group == 'shared technical', 'plain', 'bold')),
        size = 10 * 0.28, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4, seed = 1
    ) +
    scale_color_manual(values = fc_colors, breaks = names(fc_colors), labels = fc_labels) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.5, alpha = 1))) +
    xlab(expression(log[2] * '(fold change)')) +
    ylab(expression(-log[10] * '(adj. P)')) +
    ggtitle('Myeloid: GEM-X vs 5\' v2') +
    theme(legend.title = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 10))

size_df <- chem %>%
    filter(sig) %>%
    mutate(category = fc_category(logFC, sig)) %>%
    count(lineage, category)
bar_df <- size_df %>%
    pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
    mutate(y = match(lineage, rev(lineages))) %>%
    reframe(
        lineage = rep(lineage, 4),
        y = rep(y, 4),
        category = rep(c('lower_big', 'lower_small', 'higher_small', 'higher_big'), each = n()),
        n = c(lower_big, lower_small, higher_small, higher_big),
        xmin = c(-(lower_small + lower_big), -lower_small, rep(0, n()), higher_small),
        xmax = c(-lower_small, rep(0, n()), higher_small, higher_small + higher_big)
    ) %>%
    mutate(category = factor(category, levels = names(fc_colors)))
bar_half <- 0.3
b_orig <- bar_df %>%
    ggplot() +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - bar_half, ymax = y + bar_half, fill = category)) +
    geom_vline(xintercept = 0, color = 'black') +
    geom_text(data = filter(bar_df, category %in% c('lower_big', 'lower_small')),
              aes(x = (xmin + xmax) / 2, y = y + bar_half, label = n), size = 2.3, vjust = -0.3) +
    geom_text(data = filter(bar_df, category == 'higher_small'),
              aes(x = 0, y = y + bar_half, label = n), size = 2.3, vjust = -0.3, hjust = -0.15) +
    geom_text(data = filter(bar_df, category == 'higher_big'),
              aes(x = xmax, y = y, label = n), size = 2.3, hjust = -0.25) +
    scale_fill_manual(values = fc_colors, breaks = names(fc_colors)[1:4], labels = fc_labels[1:4]) +
    scale_y_continuous(breaks = seq_along(lineages), labels = rev(lineages), expand = expansion(add = 0.5)) +
    scale_x_continuous(labels = abs, expand = expansion(mult = c(0.03, 0.08))) +
    xlab('Chemistry-sensitive genes (adj.P < 0.05)') +
    ylab('') +
    guides(fill = guide_legend(ncol = 2)) +
    theme(legend.title = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())

#---------------------------------------------------------#
# ploure
#---------------------------------------------------------#

directions <- c('Lower in GEM-X', 'Higher in GEM-X')
v_new <- VolcanoPlot(filter(chem, lineage == 'myeloid'), label = label_technical, label_bold = label_biology,
                     direction_labels = directions, text_size = 10 * 0.28, title = 'Myeloid: GEM-X vs 5\' v2') +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
b_new <- DECountPlot(chem, group_by = 'lineage', direction_labels = directions,
                     x_label = 'Chemistry-sensitive genes (adj.P < 0.05)')

#---------------------------------------------------------#
# layer comparison (positions, colors, labels; row order ignored)
#---------------------------------------------------------#

cols <- c('x', 'y', 'xmin', 'xmax', 'ymin', 'ymax', 'xintercept', 'yintercept', 'colour', 'fill', 'label', 'fontface',
          'size', 'alpha', 'shape', 'linetype')
canon <- function(d){
    d <- as.data.frame(d)[, intersect(cols, colnames(d)), drop = FALSE]
    d[] <- lapply(d, function(v) if(is.factor(v)) as.character(v) else v)
    d <- d[do.call(order, unname(as.list(d))), , drop = FALSE]
    rownames(d) <- NULL
    d
}
compare <- function(name, a, b){
    la <- ggplot_build(a)$data
    lb <- ggplot_build(b)$data
    for(i in seq_len(min(length(la), length(lb)))){
        same <- isTRUE(all.equal(canon(la[[i]]), canon(lb[[i]]), check.attributes = FALSE))
        cat(sprintf('%s layer %d (%s): %s\n', name, i, class(a$layers[[i]]$geom)[1],
                    if(same) 'identical' else 'DIFFERENT'))
    }
}
compare('volcano', v_orig, v_new)
compare('counts', b_orig, b_new)

# count labels: placed differently by design (strong segments at the bar ends on both sides), but the same numbers
count_labels <- function(p){
    d <- ggplot_build(p)$data
    sort(as.numeric(unlist(lapply(d[sapply(d, function(x) 'label' %in% colnames(x))], `[[`, 'label'))))
}
cat('count labels, same numbers:', identical(count_labels(b_orig), count_labels(b_new)), '\n')

# free(): the volcano axis titles sit next to their own axes, not pushed out by the long lineage names below
fig <- wrap_plots(free(v_orig, type = 'label', side = 'l'), free(v_new, type = 'label', side = 'l'), b_orig, b_new,
                  ncol = 2)
cairo_pdf(paste0(fig_dir, 'de_plots_ploure.pdf'), width = 11, height = 9.5)
print(fig)
invisible(dev.off())
png(paste0(fig_dir, 'de_plots_ploure.png'), width = 11, height = 9.5, units = 'in', res = 80, type = 'cairo')
print(fig)
invisible(dev.off())

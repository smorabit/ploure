# launch an interactive node
# interactive 4 64G 2:00:00 genD
#
# to launch (PLEURA's R environment, with ploure installed into validation/lib, not into the environment):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     R CMD INSTALL -l validation/lib .
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_comparison.R

#---------------------------------------------------------#
# Redraw PLEURA figures with the original utilities and with ploure, and compare them
#
# The PLEURA utilities are sourced read-only into their own environment (`old`); ploure is loaded from validation/lib.
# Ported functions (PlotEmbedding, FeatureEmbedding, ZoomEmbedding) are compared on the built plot data, layer by
# layer; figures that ploure generalizes from notebook code (marker dot plot, per-cluster QC stack) are compared
# visually. Every figure is written twice to validation/figures/ (*_pleura.pdf, *_ploure.pdf), which is git-ignored.
#---------------------------------------------------------#

# plotting
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrastr)
library(ggrepel)
library(ggbeeswarm)
library(scales)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

theme_set(theme_cowplot())

select <- dplyr::select
filter <- dplyr::filter

setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- 'analysis_v2/overview/data/'

# the original PLEURA utilities, in their own environment
old <- new.env()
for(f in c('plotting_functions.R', 'zoom_embedding.R', 'figure_functions.R')){
    sys.source(file.path('pleura_repo/utils', f), envir = old)
}

#---------------------------------------------------------#
# helpers
#---------------------------------------------------------#

results <- list()

# compare the built data of two plots, layer by layer, on the drawn aesthetics
compare_built <- function(name, p_old, p_new, cols = c('x', 'y', 'colour', 'fill', 'label')){
    b_old <- ggplot_build(p_old)$data
    b_new <- ggplot_build(p_new)$data
    if(length(b_old) != length(b_new)){
        results[[name]] <<- data.frame(figure = name, layer = NA, identical = FALSE,
                                       note = sprintf('%d vs %d layers', length(b_old), length(b_new)))
        return(invisible())
    }
    results[[name]] <<- do.call(rbind, lapply(seq_along(b_old), function(i){
        cur_cols <- intersect(intersect(cols, colnames(b_old[[i]])), colnames(b_new[[i]]))
        # factors as text (a factor label and a character label draw the same)
        d_old <- as.data.frame(lapply(b_old[[i]][cur_cols], function(v) if(is.factor(v)) as.character(v) else v))
        d_new <- as.data.frame(lapply(b_new[[i]][cur_cols], function(v) if(is.factor(v)) as.character(v) else v))
        same <- isTRUE(all.equal(d_old, d_new, check.attributes = FALSE))
        note <- if(same) '' else if(nrow(d_old) != nrow(d_new)) sprintf('%d vs %d rows', nrow(d_old), nrow(d_new)) else
            paste(cur_cols[!sapply(cur_cols, function(cc) isTRUE(all.equal(d_old[[cc]], d_new[[cc]],
                                                                               check.attributes = FALSE)))],
                  collapse = ', ')
        data.frame(figure = name, layer = i, identical = same, note = note)
    }))
}

save_pair <- function(name, p_old, p_new, width, height){
    for(v in c('pleura', 'ploure')){
        cairo_pdf(paste0(fig_dir, name, '_', v, '.pdf'), width = width, height = height)
        print(if(v == 'pleura') p_old else p_new)
        invisible(dev.off())
    }
}

#---------------------------------------------------------#
# palettes and tables (as in v2/overview/report_figures.Rmd and malignant_zoom_umap.Rmd)
#---------------------------------------------------------#

lineage_colors <- read_csv('pleura_repo/v2/lineages/annotations/global_lineage_colors.csv', show_col_types = FALSE) %>%
    deframe()
lv1_colors <- read_csv('pleura_repo/v2/lineages/annotations/global_lv1_colors.csv', show_col_types = FALSE) %>%
    deframe()
names(lv1_colors) <- recode(names(lv1_colors), 'gd T cells' = 'γδ T cells')
cancer_colors <- read_csv('pleura_repo/v2/lineages/annotations/malignant_colors.csv', show_col_types = FALSE) %>%
    deframe()

gl <- read_csv(paste0(data_dir, 'global_umap.csv.gz'), show_col_types = FALSE, col_types = cols(sample = 'c')) %>%
    mutate(lv1 = factor(recode(lv1, 'gd T cells' = 'γδ T cells'), levels = names(lv1_colors)))

disease <- read_csv('analysis/eda/pseudobulk/data/tumor_burden.csv', show_col_types = FALSE) %>%
    distinct(patient, disease)
cancersea <- read_csv('analysis_v2/lineages/malignant/annotation/cancersea_ucell_scores.csv.gz',
                      show_col_types = FALSE)
mal <- read_csv(paste0(data_dir, 'malignant_pca_umap.csv.gz'), show_col_types = FALSE) %>%
    left_join(disease, by = 'patient') %>%
    left_join(cancersea, by = 'cell') %>%
    mutate(lv2 = factor(lv2, levels = names(cancer_colors)))

#---------------------------------------------------------#
# PlotEmbedding: the lv2-UMAP path of the supplementary figures (PlotLv2Umap), on the global table
#---------------------------------------------------------#

umap_df <- gl %>%
    transmute(cell, umap_1, umap_2, state = factor(as.character(lv1), levels = names(lv1_colors)))

set.seed(1)
p_old <- old$PlotEmbedding(old$EmbeddingObject(umap_df), group.by = 'state', x = 'umap_1', y = 'umap_2',
                           label = TRUE, point_size = 0.02, text_size = 2.5, raster_dpi = 300,
                           color_df = old$PaletteDf(lv1_colors), plot_theme = old$umap_theme()) +
    ggtitle('Cell type (lv1)')
set.seed(1)
p_new <- ploure::PlotEmbedding(umap_df, group.by = 'state', x = 'umap_1', y = 'umap_2', label = TRUE,
                               point_size = 0.02, text_size = 2.5, raster_dpi = 300, color_df = lv1_colors,
                               plot_theme = ThemeEmbedding(aspect_ratio = NULL)) +
    ggtitle('Cell type (lv1)')
compare_built('plot_embedding_lv1', p_old, p_new)
save_pair('plot_embedding_lv1', p_old, p_new, 8, 7)

# corner axes and split panels
set.seed(1)
p_old <- old$PlotEmbedding(old$EmbeddingObject(gl), group.by = 'lineage', x = 'umap_1', y = 'umap_2',
                           split.by = 'chemistry', plot_under = TRUE, label = FALSE, point_size = 0.05,
                           corner_axes = TRUE, color_df = old$PaletteDf(lineage_colors),
                           plot_theme = old$umap_theme())
set.seed(1)
p_new <- ploure::PlotEmbedding(gl, group.by = 'lineage', x = 'umap_1', y = 'umap_2', split.by = 'chemistry',
                               plot_under = TRUE, label = FALSE, point_size = 0.05, corner_axes = TRUE,
                               color_df = lineage_colors, plot_theme = ThemeEmbedding(aspect_ratio = NULL))
for(cur in names(p_old)){
    compare_built(paste0('plot_embedding_split_', cur), p_old[[cur]], p_new[[cur]])
}
save_pair('plot_embedding_split', wrap_plots(p_old), wrap_plots(p_new), 14, 7)

#---------------------------------------------------------#
# FeatureEmbedding: a numeric column (cached table) and a gene (Seurat object)
#---------------------------------------------------------#

mal_obj <- old$EmbeddingObject(mal, coord_cols = c('umap_1', 'umap_2'))
p_old <- old$FeatureEmbedding(mal_obj, features = c('proliferation', 'hypoxia'), plot_min = 'q1',
                              plot_max = 'q99', corner_axes = TRUE, combine = FALSE)
p_new <- ploure::FeatureEmbedding(mal, features = c('proliferation', 'hypoxia'), reduction = 'umap',
                                  plot_min = 'q1', plot_max = 'q99', corner_axes = TRUE, combine = FALSE)
for(cur in names(p_old)) compare_built(paste0('feature_embedding_', cur), p_old[[cur]], p_new[[cur]])
save_pair('feature_embedding_cancersea', wrap_plots(p_old), wrap_plots(p_new), 12, 6)

b_plasma <- readRDS('analysis_v2/objects/b_plasma_v2.rds')
p_old <- old$FeatureEmbedding(b_plasma, features = 'JCHAIN', reduction = 'umap', plot_min = 0, plot_max = 'q99')
p_new <- ploure::FeatureEmbedding(b_plasma, features = 'JCHAIN', reduction = 'umap', plot_min = 0, plot_max = 'q99')
compare_built('feature_embedding_jchain_seurat', p_old, p_new)
save_pair('feature_embedding_jchain', p_old, p_new, 7, 6)

#---------------------------------------------------------#
# ZoomEmbedding: the four figures of malignant_zoom_umap.Rmd
#---------------------------------------------------------#

ramps <- list(
    LC = c('#c7e9c0', '#41ab5d', '#00441b', '#08519c', '#6baed6'),
    BC = c('#ffeda0', '#fec44f', '#fd8d3c', '#e31a1c', '#800026'),
    OC = c('#c51b7d', '#f1b6da', '#8c510a')
)
patients <- mal %>%
    distinct(patient, disease) %>%
    mutate(disease = factor(disease, levels = names(ramps))) %>%
    arrange(disease, patient)
patient_colors <- patients %>%
    group_by(disease) %>%
    mutate(color = colorRampPalette(ramps[[as.character(disease[1])]])(n())) %>%
    ungroup() %>%
    select(patient, color) %>%
    deframe()
patient_labels <- setNames(paste(patients$patient, patients$disease), patients$patient)
family_colors <- setNames(sapply(ramps, function(r) colorRampPalette(r)(5)[3]),
                          c('Lung cancer (LC)', 'Breast cancer (BC)', 'Ovarian cancer (OC)'))

# the grouped palette of ploure gives the same patient and family colors
patient_colors_new <- GroupedPalette(patients$patient, as.character(patients$disease), ramps = ramps,
                                     family_labels = c(LC = 'Lung cancer (LC)', BC = 'Breast cancer (BC)',
                                                       OC = 'Ovarian cancer (OC)'))
results$grouped_palette <- data.frame(
    figure = 'grouped_palette', layer = NA,
    identical = identical(unname(patient_colors_new[names(patient_colors)]), unname(patient_colors)) &&
        identical(unname(attr(patient_colors_new, 'families')[names(family_colors)]), unname(family_colors)),
    note = 'GroupedPalette vs the notebook ramps'
)

zoom_args <- list(
    full = gl, full_reduction = c('umap_scvi_1', 'umap_scvi_2'), full_color_by = 'lv1', full_colors = lv1_colors,
    subset = mal, subset_reduction = c('umap_1', 'umap_2'),
    titles = c('All cells (scVI)', sprintf('Malignant cells (PCA), n = %s', scales::comma(nrow(mal))))
)
zoom_variants <- list(
    malignant_zoom_umap = list(subset_color_by = 'lv2', subset_colors = cancer_colors,
                               legend_titles = c('Cell type (lv1)', 'Cancer type')),
    malignant_zoom_umap_patient = list(subset_color_by = 'patient', subset_colors = patient_colors,
                                       subset_label_by = 'patient', label_map = patient_labels,
                                       legend_colors = list(full = NULL, subset = family_colors),
                                       legend_ncol = c(3, 1),
                                       legend_titles = c('Cell type (lv1)', 'Patient, colored by cancer type')),
    malignant_zoom_umap_proliferation = list(subset_color_by = 'proliferation', plot_min = 'q1', plot_max = 'q99',
                                             legend_titles = c('Cell type (lv1)', 'CancerSEA proliferation (UCell)'))
)
for(cur in names(zoom_variants)){
    p_old <- do.call(old$ZoomEmbedding, c(zoom_args, zoom_variants[[cur]]))
    p_new <- do.call(ploure::ZoomEmbedding, c(zoom_args, zoom_variants[[cur]]))
    compare_built(cur, p_old[[1]], p_new[[1]])
    save_pair(cur, p_old, p_new, 13, 7.5)
}

# the Seurat path, with a gene
jchain_args <- list(
    full = gl, full_reduction = c('umap_scvi_1', 'umap_scvi_2'), full_color_by = 'lv1', full_colors = lv1_colors,
    subset = b_plasma, subset_reduction = 'umap', subset_color_by = 'JCHAIN', plot_min = 0, plot_max = 'q99',
    titles = c('All cells (scVI)', sprintf('B/plasma cells (scVI), n = %s', scales::comma(ncol(b_plasma)))),
    legend_titles = c('Cell type (lv1)', 'JCHAIN')
)
p_old <- do.call(old$ZoomEmbedding, jchain_args)
p_new <- do.call(ploure::ZoomEmbedding, jchain_args)
compare_built('b_plasma_zoom_umap_jchain', p_old[[1]], p_new[[1]])
save_pair('b_plasma_zoom_umap_jchain', p_old, p_new, 13, 7.5)
rm(b_plasma)
gc()

#---------------------------------------------------------#
# marker dot plot (report_figures.Rmd, dataset_wide_markers.pdf) vs MarkerDotPlot
#---------------------------------------------------------#

lineage_labels <- c('t_nk' = 'T/NK', 'b_plasma' = 'B/plasma', 'myeloid' = 'Myeloid', 'stromal' = 'Stromal',
                    'malignant' = 'Malignant', 'low_quality' = 'Low quality (removed)',
                    'doublet' = 'Doublets (removed)')
dw_colors <- c(lineage_colors, 'Low quality (removed)' = '#d9d9d9', 'Doublets (removed)' = '#252525')

dw <- read_csv(paste0(data_dir, 'dataset_wide_umap.csv.gz'), show_col_types = FALSE,
               col_types = cols(leiden_1 = 'c', sample = 'c')) %>%
    mutate(lineage = factor(lineage_labels[lineage], levels = names(dw_colors)))
cluster_lineage <- read_csv('analysis_v2/dataset_wide/data/lineage_final_leiden_1.csv', show_col_types = FALSE,
                            col_types = cols(leiden_1 = 'c')) %>%
    select(leiden_1, lineage) %>%
    mutate(lineage = factor(lineage_labels[lineage], levels = names(dw_colors)))
cluster_order <- dw %>%
    count(leiden_1) %>%
    left_join(cluster_lineage, by = 'leiden_1') %>%
    arrange(lineage, desc(n)) %>%
    pull(leiden_1)

panel_raw <- read_csv('analysis_v2/dataset_wide/data/marker_panel_leiden_1.csv', show_col_types = FALSE,
                      col_types = cols(leiden_1 = 'c')) %>%
    mutate(group = recode(group, 'Epithelial/tumour' = 'Epithelial/tumor', 'B/Plasma' = 'B/plasma'))

# the notebook version, verbatim
panel <- panel_raw %>%
    group_by(gene) %>%
    mutate(scaled = (mean_expr - min(mean_expr)) / (max(mean_expr) - min(mean_expr))) %>%
    ungroup() %>%
    left_join(cluster_lineage, by = 'leiden_1') %>%
    mutate(
        leiden_1 = factor(leiden_1, levels = cluster_order),
        lineage = recode(lineage, 'Low quality (removed)' = 'Low quality', 'Doublets (removed)' = 'Doublets'),
        lineage = factor(lineage, levels = c(names(lineage_colors), 'Low quality', 'Doublets')),
        group = factor(group, levels = unique(group)),
        gene = factor(gene, levels = rev(unique(gene)))
    )
p_old <- ggplot(panel, aes(x = leiden_1, y = gene)) +
    geom_point(aes(size = frac_expr, color = scaled)) +
    facet_grid(group ~ lineage, scales = 'free', space = 'free', switch = 'x') +
    scale_color_gradient(low = 'grey92', high = '#08306b', name = 'Scaled\nexpression', breaks = c(0, 0.5, 1)) +
    scale_size_area(max_size = 3.2, name = 'Fraction\nexpressing', breaks = c(0.25, 0.5, 0.75)) +
    labs(x = 'Leiden cluster (resolution 1)', y = NULL) +
    theme(
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8, face = 'italic'),
        strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text.x = element_text(size = 9, angle = 90, hjust = 1),
        strip.text.y = element_text(size = 8, angle = 0, hjust = 0),
        panel.spacing = unit(0.15, 'lines'),
        panel.border = element_rect(color = 'grey70', fill = NA, linewidth = 0.3),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
    )

# ploure: gene groups as a list, clusters mapped to lineages
gene_groups <- panel_raw %>%
    distinct(group, gene) %>%
    mutate(group = factor(group, levels = unique(group))) %>%
    split(.$group) %>%
    lapply(function(d) d$gene)
x_groups <- cluster_lineage %>%
    mutate(lineage = recode(as.character(lineage), 'Low quality (removed)' = 'Low quality',
                            'Doublets (removed)' = 'Doublets')) %>%
    arrange(factor(lineage, levels = c(names(lineage_colors), 'Low quality', 'Doublets'))) %>%
    select(leiden_1, lineage) %>%
    deframe()
dots <- panel_raw %>%
    mutate(leiden_1 = factor(leiden_1, levels = cluster_order))
p_new <- MarkerDotPlot(dots, gene_groups = gene_groups, x_groups = x_groups, x = 'leiden_1',
                       x_label = 'Leiden cluster (resolution 1)') +
    theme(
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
    )
compare_built('marker_dot_plot', p_old, p_new, cols = c('x', 'y', 'colour', 'size', 'PANEL'))
save_pair('dataset_wide_markers', p_old, p_new, 11, 11.5)

#---------------------------------------------------------#
# per-cluster QC stack (report_figures.Rmd, dataset_wide_overview.pdf panel c) vs DistributionStack
#---------------------------------------------------------#

qc_levels <- c('Cells', 'UMIs (log10)', 'Genes', 'Mito (%)', 'Doublet score')
qc_cells <- dw %>%
    mutate(
        leiden_1 = factor(leiden_1, levels = cluster_order),
        lineage = factor(cluster_lineage$lineage[match(leiden_1, cluster_lineage$leiden_1)], levels = names(dw_colors)),
        log10_umi = log10(nCount_cb)
    )
qc_long <- qc_cells %>%
    select(leiden_1, lineage, log10_umi, nFeature_cb, pct_mt_cb, scDblFinder_score) %>%
    pivot_longer(c(log10_umi, nFeature_cb, pct_mt_cb, scDblFinder_score), names_to = 'metric') %>%
    mutate(metric = factor(recode(metric, 'log10_umi' = 'UMIs (log10)', 'nFeature_cb' = 'Genes',
                                  'pct_mt_cb' = 'Mito (%)', 'scDblFinder_score' = 'Doublet score'),
                           levels = qc_levels))
set.seed(12345)
qc_points <- qc_long %>%
    group_by(leiden_1, metric) %>%
    slice_sample(n = 300) %>%
    ungroup()
qc_counts <- dw %>%
    count(leiden_1, name = 'value') %>%
    left_join(cluster_lineage, by = 'leiden_1') %>%
    mutate(leiden_1 = factor(leiden_1, levels = cluster_order), metric = factor('Cells', levels = qc_levels))
p_old <- ggplot(mapping = aes(x = leiden_1, y = value)) +
    geom_col(data = qc_counts, aes(fill = lineage), width = 0.8) +
    rasterise(geom_quasirandom(data = qc_points, aes(color = lineage), size = 0.1, alpha = 0.5, width = 0.35),
              dpi = 300) +
    geom_boxplot(data = qc_long, outlier.shape = NA, fill = NA, width = 0.7, linewidth = 0.2) +
    facet_grid(metric ~ ., scales = 'free_y', switch = 'y') +
    scale_fill_manual(values = dw_colors, guide = 'none') +
    scale_color_manual(values = dw_colors, guide = 'none') +
    scale_y_continuous(labels = label_number(big.mark = ',', accuracy = NULL), n.breaks = 4) +
    labs(x = 'Leiden cluster (resolution 1)', y = NULL) +
    theme(
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9),
        panel.spacing.y = unit(0.4, 'lines')
    )
p_new <- DistributionStack(qc_cells, group_by = 'leiden_1',
                           metrics = c('UMIs (log10)' = 'log10_umi', 'Genes' = 'nFeature_cb',
                                       'Mito (%)' = 'pct_mt_cb', 'Doublet score' = 'scDblFinder_score'),
                           color_by = 'lineage', colors = dw_colors, raster_dpi = 300,
                           x_label = 'Leiden cluster (resolution 1)') +
    theme(
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9)
    )
# the boxes (all cells) and the count bars should match; the point subsamples differ (per-group sampling order)
# rows sorted by panel and position (the stacking order of the bars differs, not their values)
sort_rows <- function(d) as.data.frame(d)[order(d$PANEL, d$x), ]
b_old <- lapply(ggplot_build(p_old)$data, sort_rows)
b_new <- lapply(ggplot_build(p_new)$data, sort_rows)
cols_box <- c('x', 'ymin', 'lower', 'middle', 'upper', 'ymax', 'PANEL')
cols_bar <- c('x', 'y', 'fill', 'PANEL')
results$distribution_stack <- data.frame(
    figure = 'distribution_stack', layer = c(1, 3),
    identical = c(isTRUE(all.equal(as.data.frame(b_old[[1]][cols_bar]), as.data.frame(b_new[[1]][cols_bar]),
                                   check.attributes = FALSE)),
                  isTRUE(all.equal(as.data.frame(b_old[[3]][cols_box]), as.data.frame(b_new[[3]][cols_box]),
                                   check.attributes = FALSE))),
    note = c('count bars', 'boxplots (all cells)')
)
save_pair('dataset_wide_qc_stack', p_old, p_new, 12, 6)

#---------------------------------------------------------#
# summary
#---------------------------------------------------------#

summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL
write_csv(summary_df, paste0(fig_dir, 'comparison.csv'))
print(summary_df, right = FALSE)

# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_milo.R

#---------------------------------------------------------#
# The PLEURA v2 Milo DA figures (v2/composition/milo_figures.Rmd) drawn with DABeeswarmPlot and DAEmbeddingPlot
#
# For each lineage's primary run (seed 12345, the default k) and the first primary contrast, draws the beeswarm and
# the neighbourhood graph on the UMAP with PLEURA's DABoxPlot / DAReductionPlot (utils/milo_functions.R, sourced, not
# modified) and with ploure, compares the built layers, and writes validation/figures/milo_ploure.pdf (PLEURA's panel
# then ploure's, per lineage).
#---------------------------------------------------------#

suppressPackageStartupMessages({
    library(tidyverse)
    library(miloR)
    library(SingleCellExperiment)
    library(patchwork)
})

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

setwd(proj_root)
pleura <- new.env()
sys.source('pleura_repo/utils/milo_functions.R', envir = pleura)

theme_set(ThemePublication())
SetFigureFont()
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml('pleura_repo/v2/composition/milo_config.yml')
signif_thresh <- cfg$params$signif
color_limit <- 3
default_k <- c(t_nk = 100, myeloid = 100, b_plasma = 100, stromal = 35)
contrasts <- c('lung_vs_breast', 'breast_vs_hf')
run_dir <- function(lin) sprintf('analysis_v2/composition/milo/%s/seed%s_k%s/', lin, cfg$seeds$primary, default_k[[lin]])

same_layer <- function(p_old, p_new, i_old, i_new, cols){
    set.seed(1)
    a <- layer_data(p_old, i_old)
    set.seed(1)
    b <- layer_data(p_new, i_new)
    if(nrow(a) == 0 || nrow(b) == 0) return(nrow(a) == nrow(b))  # e.g. no significant neighbourhoods
    isTRUE(all.equal(as.data.frame(a[, cols]), as.data.frame(b[, cols]), check.attributes = FALSE))
}

results <- list()
pages <- list()
for(lin in names(default_k)){
    contrast <- contrasts[file.exists(paste0(run_dir(lin), contrasts, '_milo_DA.csv'))][1]
    da <- read_csv(paste0(run_dir(lin), contrast, '_milo_DA.csv'), show_col_types = FALSE) %>% as.data.frame()
    milo_obj <- readRDS(paste0(run_dir(lin), 'milo.rds'))
    title <- paste0(lin, ': ', contrast, ' (SpatialFDR < ', signif_thresh, ')')

    # beeswarm
    p_old <- pleura$DABoxPlot(da, signif_col = 'SpatialFDR', signif_thresh = signif_thresh,
                              color_limit = color_limit) + ggtitle(title)
    p_new <- DABeeswarmPlot(da, signif_col = 'SpatialFDR', signif_thresh = signif_thresh,
                            color_limit = color_limit, title = title)
    old_means <- trimws(layer_data(p_old, 3)$label)
    new_means <- ggplot_build(p_new)$layout$panel_params[[1]]$y.sec$get_labels()
    old_groups <- levels(p_old$data$cluster)
    results[[paste(lin, 'beeswarm')]] <- c(
        points = same_layer(p_old, p_new, 1, 1, c('x', 'y', 'colour')),
        boxes = same_layer(p_old, p_new, 2, 2, c('xmin', 'xlower', 'xmiddle', 'xupper', 'xmax', 'y')),
        groups = identical(old_groups, levels(p_new$data$group)),
        # PLEURA's labels are in its group order, like the levels of the right-hand axis
        means = identical(old_means, as.character(new_means))
    )

    # neighbourhood graph on the UMAP
    e_old <- pleura$DAReductionPlot(milo_obj, da, reduction = 'UMAP', signif_col = 'SpatialFDR',
                                    signif_thresh = signif_thresh, color_limit = color_limit,
                                    point_sizes = c(0.1, 4)) + ggtitle(title) + theme(aspect.ratio = 1)
    e_new <- DAEmbeddingPlot(milo_obj, da, reduction = 'UMAP', signif_col = 'SpatialFDR',
                             signif_thresh = signif_thresh, color_limit = color_limit, title = title)
    cols <- c('x', 'y', 'colour', 'size')
    results[[paste(lin, 'embedding')]] <- c(
        not_significant = same_layer(e_old, e_new, 1, 1, cols),
        significant = same_layer(e_old, e_new, 2, 2, cols)
    )
    pages[[lin]] <- list((p_old + labs(subtitle = 'PLEURA DABoxPlot')) | (p_new + labs(subtitle = 'ploure')),
                         (e_old + labs(subtitle = 'PLEURA DAReductionPlot')) | (e_new + labs(subtitle = 'ploure')))
    rm(milo_obj); gc()
}

for(name in names(results)) cat(sprintf('%-22s %s\n', name, paste(names(results[[name]]), results[[name]], sep = ': ',
                                                                    collapse = '  ')))
stopifnot(all(unlist(results)))

FigurePDF(paste0(fig_dir, 'milo_ploure.pdf'), width = 15, height = 7.5)
for(page in unlist(pages, recursive = FALSE)) print(page)
invisible(dev.off())

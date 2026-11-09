# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_qc_summary.R

#---------------------------------------------------------#
# The PLEURA v2 QC summary (v2/overview/qc_summary_figure.Rmd) drawn with QCSummary
#
# Reads every CellBender cell (analysis_v2/doublets/*_scdblfinder.csv.gz) and the cells per sample at each stage that
# the notebook wrote (analysis_v2/overview/data/cells_per_sample_stage.csv), and writes
# validation/figures/qc_summary_ploure.pdf, to compare with analysis_v2/overview/figures/qc_summary.pdf.
#---------------------------------------------------------#

library(data.table)
library(tidyverse)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

theme_set(ThemePublication())

setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- 'analysis_v2/overview/data/'

qc <- list(min_umi = 500, min_genes = 250, max_pct_mt = 50)
chemistry_colors <- c("5' v2" = '#1b9e77', 'GEM-X OCM' = '#d95f02')
chemistry_labels <- c('5p_v2' = "5' v2", 'GEM-X_OCM' = 'GEM-X OCM')
disease_labels <- c('LC' = 'Lung cancer', 'BC' = 'Breast cancer', 'OC' = 'Ovarian cancer', 'HF' = 'Heart failure')
disease_colors <- c('Lung cancer' = '#6baed6', 'Breast cancer' = '#e7969c', 'Ovarian cancer' = '#e6c35c',
                    'Heart failure' = '#8fbc8f')
stage_levels <- c('Cell Ranger', 'CellBender', 'Permissive QC', 'Dataset-wide', 'Final')

# per-sample stage counts, in the notebook's order (disease, chemistry, sample ID)
per_sample <- read_csv(paste0(data_dir, 'cells_per_sample_stage.csv'), col_types = cols(sample = 'c'),
                       show_col_types = FALSE) %>%
    mutate(disease = factor(disease, levels = disease_labels))

cb_cells <- rbindlist(lapply(Sys.glob('analysis_v2/doublets/*_scdblfinder.csv.gz'), fread,
                             colClasses = list(character = c('sample', 'library', 'barcode')))) %>%
    as_tibble() %>%
    select(sample, library, chemistry, nCount_cb, nFeature_cb, pct_mt_cb, passes_loose_qc) %>%
    mutate(chemistry = factor(chemistry_labels[chemistry], levels = names(chemistry_colors)),
           log10_umi = log10(nCount_cb))

stopifnot(nrow(cb_cells) == 511539, all(cb_cells$sample %in% per_sample$sample))

p <- QCSummary(
    cb_cells, per_sample,
    stages = setNames(stage_levels, stage_levels),
    metrics = c('UMIs per cell' = 'nCount_cb', 'Genes per cell' = 'nFeature_cb', 'Mitochondrial (%)' = 'pct_mt_cb'),
    quality_metrics = c('UMIs (log10)' = 'log10_umi', 'Genes' = 'nFeature_cb', 'Mitochondrial (%)' = 'pct_mt_cb'),
    sample_col = 'sample', label_col = 'label', library_col = 'library',
    group_col = 'disease', group_colors = disease_colors, group_labels = setNames(names(disease_labels), disease_labels),
    color_by = 'chemistry', colors = chemistry_colors,
    thresholds = list(nCount_cb = qc$min_umi, nFeature_cb = qc$min_genes, pct_mt_cb = qc$max_pct_mt),
    trans = c(nCount_cb = 'log10', nFeature_cb = 'log10', pct_mt_cb = 'sqrt'),
    breaks = list(nCount_cb = c(100, 1000, 10000, 100000), nFeature_cb = c(100, 1000, 10000),
                  pct_mt_cb = c(0, 5, 10, 25, 50, 100)),
    pass_col = 'passes_loose_qc', raster_dpi = 300, group_title = 'Disease', color_title = 'Chemistry'
)

cairo_pdf(paste0(fig_dir, 'qc_summary_ploure.pdf'), width = 15, height = 16)
print(p)
invisible(dev.off())

# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_cohort_matrix.R

#---------------------------------------------------------#
# The PLEURA v2 cohort matrix (v2/overview/cohort_matrix.Rmd) drawn with CohortMatrix
#
# Builds the same sample table as the notebook and writes validation/figures/cohort_matrix_ploure.pdf, to compare with
# analysis_v2/overview/figures/cohort_matrix.pdf. The tile values of each block are also checked against the
# notebook's long tables.
#---------------------------------------------------------#

library(tidyverse)
library(patchwork)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- 'analysis_v2/overview/data/'

#---------------------------------------------------------#
# palettes and the sample table, as in the notebook
#---------------------------------------------------------#

disease_labels <- c('LC' = 'Lung cancer', 'BC' = 'Breast cancer', 'OC' = 'Ovarian cancer', 'HF' = 'Heart failure')
disease_colors <- c('Lung cancer' = '#6baed6', 'Breast cancer' = '#e7969c', 'Ovarian cancer' = '#e6c35c',
                    'Heart failure' = '#8fbc8f')
chemistry_colors <- c("5' v2" = '#1b9e77', 'GEM-X OCM' = '#d95f02')
site_colors <- c('Sant Pau' = '#7570b3', 'Lleida' = '#e7298a')
available_colors <- c('Yes' = 'grey35', 'No' = 'grey92')
histology_colors <- c('LUAD' = '#2171b5', 'LUAD (to confirm)' = '#9ecae1', 'LUSC' = '#6a51a3', 'SCLC' = '#cb181d',
                      'NSCLC-NOS' = '#fd8d3c', 'Unknown' = '#bdbdbd')
sex_colors <- c('Female' = '#c51b7d', 'Male' = '#4d9221')
age_colors <- c('< 60' = '#fdd49e', '60-70' = '#fc8d59', '> 70' = '#b30000')
result_colors <- c('Positive' = '#b2182b', 'Negative' = '#4393c3', 'Unclear' = '#d9d9d9')
ict_colors <- c('Yes' = '#542788', 'No' = '#d8daeb')
lineage_colors <- read_csv('pleura_repo/v2/lineages/annotations/global_lineage_colors.csv', show_col_types = FALSE) %>%
    deframe()

samples <- read_csv('analysis/eda/pseudobulk/data/tumor_burden.csv', col_types = cols(sample = 'c'),
                    show_col_types = FALSE) %>%
    select(sample, sample_id, patient, disease, chemistry)
loading <- read_csv('data/metadata/sample_loading.csv', show_col_types = FALSE) %>%
    select(sample_id, starts_with('load_frac'))
clinical <- read_csv('data/metadata/clinical_basic.csv', show_col_types = FALSE) %>%
    rename(sample_id = patient)
histology <- read_csv('data/metadata/lung_histology.csv', show_col_types = FALSE) %>%
    mutate(histology = case_when(
        lung_histology == 'unknown' ~ 'Unknown',
        str_detect(source, 'to confirm') ~ paste(lung_histology, '(to confirm)'),
        TRUE ~ lung_histology
    )) %>%
    select(patient, histology)

meta <- samples %>%
    left_join(loading, by = 'sample_id') %>%
    left_join(clinical, by = 'sample_id') %>%
    left_join(histology, by = 'patient') %>%
    mutate(
        disease = factor(disease_labels[disease], levels = disease_labels),
        chemistry = recode(chemistry, '5p_v2' = "5' v2", 'GEM-X_OCM' = 'GEM-X OCM'),
        label = if_else(sample_id != patient, paste0(sample_id, '*'), sample_id),
        scrna = 'Yes',
        sctcr = 'Yes',
        across(starts_with('load_frac'), ~ 100 * .x)
    ) %>%
    arrange(disease, chemistry, sample_id)

cells <- read_csv(paste0(data_dir, 'cells_per_sample_lineage.csv'), col_types = cols(sample = 'c'),
                  show_col_types = FALSE) %>%
    select(sample, lineage, n) %>%
    mutate(lineage = factor(lineage, levels = names(lineage_colors)))

stopifnot(nrow(meta) == 55, sum(cells$n) == 282391)

#---------------------------------------------------------#
# the matrix
#---------------------------------------------------------#

blocks <- list(
    list(vars = c('Disease' = 'disease'), colors = disease_colors),
    list(vars = c('Chemistry' = 'chemistry', 'Site' = 'site'), colors = c(chemistry_colors, site_colors),
         title = 'Chemistry / site'),
    list(vars = c('Loaded CD45+' = 'load_frac_cd45pos', 'Loaded EPCAM+' = 'load_frac_epcam',
                  'Loaded DN' = 'load_frac_dn', 'Loaded CD45−' = 'load_frac_cd45neg_unsplit'),
         limits = c(0, 100), title = 'Cells loaded (%)'),
    list(vars = c('scRNA-seq' = 'scrna', 'scTCR-seq' = 'sctcr', 'Bulk RNA-seq' = 'bulk_rnaseq'),
         colors = available_colors, title = 'Data available'),
    list(vars = c('Lung histology' = 'histology'), colors = histology_colors),
    list(vars = c('Sex' = 'sex', 'Age' = 'age_group'), colors = c(sex_colors, age_colors), title = 'Sex / age'),
    list(vars = c('Cytology' = 'cytology', 'Biopsy' = 'biopsy'), colors = result_colors,
         title = 'Tumor cells detected'),
    list(vars = c('ICT after effusion' = 'ict_after_effusion'), colors = ict_colors)
)

p <- CohortMatrix(
    meta, blocks, sample_col = 'sample', group_col = 'disease', label_col = 'label',
    group_labels = setNames(names(disease_labels), disease_labels),
    bar_data = cells, bar_fill = 'lineage', bar_colors = lineage_colors, bar_title = 'Lineage (v2)',
    bar_label = 'Cells (v2)'
)
print(attr(p, 'pdf_size'))

cairo_pdf(paste0(fig_dir, 'cohort_matrix_ploure.pdf'), width = 17, height = 6.8)
print(p)
invisible(dev.off())

#---------------------------------------------------------#
# tile values: each block's fill per sample and row, against the sample table
#---------------------------------------------------------#

panels <- attr(p, 'panels')
n_bad <- 0
for(i in seq_along(blocks)){
    plot_df <- panels[[paste0('block', i)]]$data
    vars <- blocks[[i]]$vars
    for(v in names(vars)){
        drawn <- plot_df %>%
            filter(row == v)
        expected <- meta[[vars[[v]]]][match(as.character(drawn$label), meta$label)]
        same <- identical(as.character(drawn$value), as.character(expected))
        if(!same) n_bad <- n_bad + 1
        cat(sprintf('%-20s %s\n', v, if(same) 'ok' else 'DIFFERENT'))
    }
}
bar_built <- panels$bars$data
cat(sprintf('%-20s %s\n', 'bars (total cells)', if(sum(bar_built$value) == 282391) 'ok' else 'DIFFERENT'))
cat('columns:', paste(levels(panels$block1$data$label)[1:5], collapse = ', '), '...\n')
cat(n_bad, 'rows differ\n')

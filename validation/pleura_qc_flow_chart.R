# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_qc_flow_chart.R

#---------------------------------------------------------#
# The PLEURA v2 QC flow chart (v2/overview/qc_flow_chart.Rmd) drawn with QCFlowChart
#
# Reads the same tables (analysis_v2/overview/data/qc_flow_steps.csv and qc_flow_removals.csv), maps them to the
# QCFlowChart inputs (one row per box, removals keyed by the step they are shown under) and writes
# validation/figures/qc_flow_chart_ploure.pdf, to compare with analysis_v2/overview/figures/qc_flow_chart.pdf.
#---------------------------------------------------------#

library(tidyverse)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- 'analysis_v2/overview/data/'

steps <- read_csv(paste0(data_dir, 'qc_flow_steps.csv'), show_col_types = FALSE)
removals <- read_csv(paste0(data_dir, 'qc_flow_removals.csv'), show_col_types = FALSE)
lineage_colors <- read_csv('pleura_repo/v2/lineages/annotations/global_lineage_colors.csv', show_col_types = FALSE) %>%
    deframe()
lineage_labels <- c('t_nk' = 'T/NK', 'b_plasma' = 'B/plasma', 'myeloid' = 'Myeloid', 'stromal' = 'Stromal',
                    'malignant' = 'Malignant')

get_step <- function(cur_lineage, cur_stage, cur_step){
    steps %>%
        filter(lineage == cur_lineage, stage == cur_stage, step == cur_step)
}

#---------------------------------------------------------#
# trunk: the library steps and the dataset-wide pass
#---------------------------------------------------------#

dbl <- get_step('all', 'library', 'scDblFinder')
qc <- get_step('all', 'library', 'Loose QC')
n_doublets <- as.numeric(str_extract(dbl$criterion, '^[0-9]+'))

trunk <- tribble(
    ~step, ~n_in, ~n_out, ~label, ~excluded, ~body,
    'Cell Ranger', NA, get_step('all', 'library', 'Cell Ranger')$n_out, '55 samples\n37 GEM wells', NA, NA,
    'CellBender', NA, get_step('all', 'library', 'CellBender')$n_out, 'called\nambient RNA removed', NA, NA,
    'Loose QC', qc$n_in, qc$n_out, '≥ 500 UMIs\n≥ 250 genes\n≤ 50% mito', 'failed loose QC', NA,
    'scDblFinder', NA, dbl$n_out, NA, NA,
        paste0(scales::comma(n_doublets), ' doublets\nflagged (', scales::percent(n_doublets / qc$n_out, accuracy = 0.1),
               '), kept\nfor cluster review'),
    'Dataset-wide pass', get_step('all', 'dataset-wide', 'Dataset-wide clusters')$n_in,
        get_step('all', 'dataset-wide', 'Dataset-wide clusters')$n_out, 'PCA, Leiden\n5 broad lineages', NA, NA
) %>%
    mutate(lineage = 'all')

#---------------------------------------------------------#
# branches: mito limit, stage 1 and stage 2 reviews, final (with the annotation removals under it)
#---------------------------------------------------------#

branch_steps <- lapply(names(lineage_labels), function(cur_lineage){
    mito <- get_step(cur_lineage, 'stage1', 'Mito limit')
    s1 <- get_step(cur_lineage, 'stage1', 'Cluster review')
    s2 <- get_step(cur_lineage, 'stage2', 'Cluster review')
    final <- get_step(cur_lineage, 'final', 'Final')
    annot <- get_step(cur_lineage, 'annotation', 'Annotation')
    states <- str_extract(final$criterion, '^[0-9]+ states')
    if(cur_lineage == 'malignant') states <- str_replace(states, 'states', 'cancer types')
    out <- tribble(
        ~step, ~n_in, ~n_out, ~label, ~excluded, ~body,
        'Mito limit', mito$n_in, mito$n_out, str_replace(mito$criterion, '<=', '≤'), 'above mito limit', NA,
        'Stage 1 review', s1$n_in, s1$n_out, 'BigSur + scVI\nLeiden clusters', NA, NA
    )
    if(nrow(s2) == 1){
        out <- add_row(out, step = 'Stage 2 review', n_in = s2$n_in, n_out = s2$n_out,
                       label = 'scVI re-run\nLeiden clusters')
    }
    final_in <- if(nrow(annot) == 1) annot$n_in else final$n_in
    add_row(out, step = 'Final', n_in = final_in, n_out = final$n_out, label = states,
            excluded = if(nrow(annot) == 1) 'unlisted subclusters' else NA) %>%
        mutate(lineage = lineage_labels[[cur_lineage]])
})

flow_steps <- bind_rows(trunk, branch_steps)

# removals keyed by the box they are shown under (annotation removals use the excluded text instead)
flow_removals <- removals %>%
    filter(stage != 'annotation') %>%
    mutate(
        step = recode(stage, 'dataset-wide' = 'Dataset-wide pass', 'stage1' = 'Stage 1 review',
                      'stage2' = 'Stage 2 review'),
        lineage = if_else(lineage == 'all', 'all', lineage_labels[lineage])
    )

p <- QCFlowChart(
    flow_steps, flow_removals, branch_col = 'lineage', branch_colors = lineage_colors,
    step_order = c('Mito limit', 'Stage 1 review', 'Stage 2 review', 'Final'),
    skip_label = 'no stage 2\n(stage 1 sufficient)',
    category_labels = c('Damaged / low quality' = 'Damaged / LQ'),
    titles = c('Dataset level (per GEM well / all cells)', 'Lineage level (one scVI per lineage)')
)

cairo_pdf(paste0(fig_dir, 'qc_flow_chart_ploure.pdf'), width = 19, height = 10)
print(p)
invisible(dev.off())

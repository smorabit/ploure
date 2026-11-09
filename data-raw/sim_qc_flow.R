# simulated QC flow tables for the QCFlowChart examples and tests (no real data)
#
# a trunk of three shared steps that fans out into three lineages; the stromal lineage has a second review stage
#
# run from the package root: Rscript data-raw/sim_qc_flow.R

sim_qc_steps <- data.frame(
    branch = c(rep('all', 3), rep('Lymphoid', 3), rep('Myeloid', 3), rep('Stromal', 4)),
    step = c('Cell calling', 'Cell QC', 'Lineage split',
             'Mito limit', 'Cluster review', 'Final',
             'Mito limit', 'Cluster review', 'Final',
             'Mito limit', 'Cluster review', 'Second review', 'Final'),
    n_in = c(5000, 5000, 4200, 1900, 1820, 1700, 1300, 1250, 1100, 800, 770, 700, 650),
    n_out = c(5000, 4200, 4000, 1820, 1700, 1700, 1250, 1100, 1100, 770, 700, 650, 650),
    label = c('8 samples', '>= 500 UMIs\n>= 250 genes', '3 lineages',
              '<= 10% mito', 'Leiden clusters', '4 states',
              '<= 15% mito', 'Leiden clusters', '3 states',
              '<= 20% mito', 'Leiden clusters', 'Leiden clusters', '2 states'),
    excluded = c(NA, 'failed cell QC', NA,
                 'above mito limit', NA, NA,
                 'above mito limit', NA, NA,
                 'above mito limit', NA, NA, NA)
)

sim_qc_removals <- data.frame(
    branch = c('all', 'all', 'Lymphoid', 'Lymphoid', 'Myeloid', 'Myeloid', 'Myeloid', 'Stromal', 'Stromal'),
    step = c('Lineage split', 'Lineage split', 'Cluster review', 'Cluster review', 'Cluster review',
             'Cluster review', 'Cluster review', 'Cluster review', 'Second review'),
    cluster = c(NA, NA, '7', '3', '2', '5', '9', '4', '1'),
    category = c('Damaged / low quality', 'Doublets', 'Damaged / low quality', 'Doublets', 'Doublets',
                 'Damaged / low quality', 'Dissociation stress', 'Doublets', 'Damaged / low quality'),
    n_removed = c(150, 50, 80, 40, 70, 60, 20, 70, 50)
)

save(sim_qc_steps, file = 'data/sim_qc_steps.rda', compress = 'xz')
save(sim_qc_removals, file = 'data/sim_qc_removals.rda', compress = 'xz')

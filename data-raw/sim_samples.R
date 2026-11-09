# simulated sample table for the CohortMatrix examples and tests (no real data): one row per sample of sim_cells,
# with sample design, data, basic clinical variables and cells per QC stage
#
# run from the package root: Rscript data-raw/sim_samples.R

load('data/sim_cells.rda')
set.seed(20261110)

sim_samples <- unique(sim_cells[c('sample', 'patient', 'disease')])
sim_samples <- sim_samples[order(factor(sim_samples$disease, levels = c('LC', 'BC', 'OC')), sim_samples$sample), ]
n <- nrow(sim_samples)
sim_samples$disease <- factor(sim_samples$disease, levels = c('LC', 'BC', 'OC'))
sim_samples$chemistry <- sample(c("5' v2", 'GEM-X'), n, replace = TRUE)
sim_samples$site <- sample(c('Site A', 'Site B'), n, replace = TRUE, prob = c(0.7, 0.3))
sim_samples$pct_cd45 <- round(runif(n, 20, 90))
sim_samples$pct_epcam <- 100 - sim_samples$pct_cd45
sim_samples$bulk_rnaseq <- sample(c('Yes', 'No'), n, replace = TRUE)
sim_samples$sex <- sample(c('Female', 'Male'), n, replace = TRUE)
sim_samples$age_group <- sample(c('< 60', '60-70', '> 70', NA), n, replace = TRUE, prob = c(0.3, 0.3, 0.3, 0.1))
# cells per sample at three QC stages (final = the cells in sim_cells)
sim_samples$cells_final <- as.integer(table(sim_cells$sample)[sim_samples$sample])
sim_samples$cells_qc <- round(sim_samples$cells_final * runif(n, 1.05, 1.25))
sim_samples$cells_called <- round(sim_samples$cells_qc * runif(n, 1.1, 1.5))
rownames(sim_samples) <- NULL

save(sim_samples, file = 'data/sim_samples.rda', compress = 'xz')

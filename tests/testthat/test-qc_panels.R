stages <- c('Called' = 'cells_called', 'QC' = 'cells_qc', 'Final' = 'cells_final')

test_that('stage counts are reshaped to one row per sample and stage', {
    long <- .StageLong(sim_samples, stages)
    expect_equal(nrow(long), 3 * nrow(sim_samples))
    expect_equal(levels(long$stage), names(stages))
    expect_equal(sum(long$n[long$stage == 'Final']), nrow(sim_cells))
    expect_error(.StageLong(sim_samples, 'nope'), 'nope')
})

test_that('sample labels follow the group order, then the row order', {
    shuffled <- sim_samples[c(9, 1, 5, 2, 8, 3, 4, 6, 7), ]
    lab <- .SampleLabels(shuffled, 'sample', NULL, 'disease')
    # groups in factor order (LC, BC, OC), samples in their row order within each group
    expect_equal(levels(lab), c('P01', 'P02', 'P03', 'P05', 'P04', 'P06', 'P09', 'P08', 'P07'))
    expect_error(.SampleLabels(rbind(sim_samples, sim_samples[1, ]), 'sample', NULL, NULL), 'unique')
})

test_that('stage totals and per-sample bars match the table', {
    p <- StageCountPlot(sim_samples, stages, fill_by = 'disease')
    b <- ggplot2::ggplot_build(p)$data[[1]]
    expect_equal(sum(b$y[b$x == 3] - b$ymin[b$x == 3]), nrow(sim_cells))
    p <- StageCountPlot(sim_samples, stages, by_sample = TRUE, facet_by = 'disease')
    b <- ggplot2::ggplot_build(p)$data[[1]]
    expect_equal(nrow(b), 3 * nrow(sim_samples))
    expect_equal(sum(b$y), sum(sim_samples[unname(stages)]))
    expect_warning(StageCountPlot(sim_samples, stages, by_sample = TRUE, fill_by = 'disease'), 'ignored')
})

test_that('stage totals can be faceted, with one total per facet and stage', {
    p <- StageCountPlot(sim_samples, stages, fill_by = 'chemistry', facet_by = 'disease')
    totals <- ggplot2::ggplot_build(p)$data[[2]]
    expect_equal(nrow(totals), 3 * length(unique(sim_samples$disease)))
})

test_that('density panels leave out values the transform cannot take', {
    cells <- sim_cells
    cells$n_umi[1:10] <- 0
    p <- MetricDensityPlot(cells, c('UMIs' = 'n_umi'), group_by = 'sample', trans = c(n_umi = 'log10'))
    expect_equal(nrow(p[[1]]$data), nrow(cells) - 10)
})

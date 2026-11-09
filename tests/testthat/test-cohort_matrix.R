cohort_blocks <- list(
    list(vars = c('Chemistry' = 'chemistry', 'Site' = 'site'),
         colors = c("5' v2" = '#1b9e77', 'GEM-X' = '#d95f02', 'Site A' = '#7570b3', 'Site B' = '#e7298a')),
    list(vars = c('Loaded CD45+' = 'pct_cd45'), limits = c(0, 100)),
    list(vars = 'sex')
)

test_that('block specs get row labels and titles', {
    b <- .CohortBlockSpec(list(vars = c('Sex' = 'sex', 'age_group')))
    expect_equal(names(b$vars), c('Sex', 'age_group'))
    expect_equal(b$title, 'Sex')
    expect_error(.CohortBlockSpec(list(colors = 'red')), 'vars')
})

test_that('blocks are categorical or continuous, rows top to bottom', {
    meta <- data.frame(label = factor(sim_samples$sample), group = sim_samples$disease)
    p <- .CohortBlock(sim_samples, meta, .CohortBlockSpec(cohort_blocks[[1]]), 'grey97')
    expect_equal(levels(p$data$row), c('Site', 'Chemistry'))
    expect_true(p$scales$get_scales('fill')$is_discrete())
    p <- .CohortBlock(sim_samples, meta, .CohortBlockSpec(cohort_blocks[[2]]), 'grey97')
    expect_false(p$scales$get_scales('fill')$is_discrete())
})

test_that('columns follow the group order, then the row order', {
    shuffled <- sim_samples[c(9, 1, 5, 2, 8, 3, 4, 6, 7), ]
    p <- CohortMatrix(shuffled, cohort_blocks, group_col = 'disease')
    expect_s3_class(p, 'patchwork')
    expect_named(attr(p, 'pdf_size'), c('width', 'height'))
    expect_error(CohortMatrix(sim_samples, list(list(vars = 'nope'))), 'nope')
    expect_error(CohortMatrix(rbind(sim_samples, sim_samples[1, ]), cohort_blocks), 'one row per sample')
})

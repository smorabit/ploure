test_that('dot plot data gives the mean and fraction expressing per group', {
    dots <- DotPlotData(sim_cells, features = c('CD3E', 'LYZ'), group_by = 'cluster')
    expect_equal(nrow(dots), 2 * nlevels(sim_cells$cluster))
    b <- sim_cells$cluster == 'B'
    row <- dots[dots$gene == 'CD3E' & dots$group == 'B', ]
    expect_equal(row$mean_expr, mean(sim_cells$CD3E[b]))
    expect_equal(row$frac_expr, mean(sim_cells$CD3E[b] > 0))
})

test_that('the same summary comes from a Seurat object', {
    skip_if_not_installed('SeuratObject')
    small <- sim_small()
    expect_equal(DotPlotData(sim_seurat(small), 'CD3E', 'cluster'), DotPlotData(small, 'CD3E', 'cluster'))
})

test_that('expression is scaled 0-1 per gene', {
    dots <- DotPlotData(sim_cells, features = c('CD3E', 'LYZ'), group_by = 'cluster')
    p <- MarkerDotPlot(dots)
    expect_equal(tapply(p$data$expr, p$data$gene, range), list(LYZ = c(0, 1), CD3E = c(0, 1)), ignore_attr = TRUE)
})

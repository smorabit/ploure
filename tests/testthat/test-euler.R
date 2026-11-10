a <- paste0('g', 1:300)
b <- paste0('g', 201:450)

test_that('overlap statistics match the sets', {
    ov <- SetOverlapStats(a, b, universe = 2000)
    expect_equal(c(ov$a_only, ov$shared, ov$b_only, ov$union), c(200, 100, 150, 450))
    expect_equal(ov$jaccard, 100 / 450)
    expect_equal(ov$overlap_coef, 100 / 250)
    expect_equal(ov$expected, 300 * 250 / 2000)
    ref <- fisher.test(matrix(c(100, 150, 200, 2000 - 450), nrow = 2), alternative = 'greater')
    expect_equal(ov$p_value, ref$p.value)
    expect_equal(ov$odds_ratio, unname(ref$estimate))
    # a universe vector gives the same as its size; counts give the same as the sets
    expect_equal(SetOverlapStats(a, b, universe = paste0('g', 1:2000)), ov)
    expect_equal(SetOverlapStats(counts = c(200, 100, 150), universe = 2000), ov)
    expect_true(is.na(SetOverlapStats(a, b)$p_value))
})

test_that('overlap statistics check their input', {
    expect_error(SetOverlapStats(a, b, universe = paste0('g', 1:400)), 'not in the universe')
    expect_error(SetOverlapStats(a, b, universe = 100), 'smaller than the union')
    expect_error(SetOverlapStats(counts = c(1, 2)), 'three')
    expect_error(SetOverlapStats(counts = c(0, 0, 5)), 'at least one')
})

test_that('named region counts give the set names', {
    counts <- .RegionCounts(c(v1_only = 10, v2_only = 30, shared = 20))
    expect_equal(as.numeric(counts), c(10, 20, 30))
    expect_equal(attr(counts, 'set_names'), c('v1', 'v2'))
})

test_that('the circle overlap area equals the shared elements', {
    geo <- .EulerGeometry(300, 250, 100)
    expect_equal(.LensArea(geo$d, geo$r1, geo$r2), 100, tolerance = 1e-3)
    expect_equal(pi * geo$r1^2, 300)
    # disjoint sets do not touch; a subset touches the edge of the larger set
    geo <- .EulerGeometry(300, 250, 0)
    expect_gt(geo$d, geo$r1 + geo$r2)
    geo <- .EulerGeometry(300, 100, 100)
    expect_equal(geo$d, geo$r1 - geo$r2)
})

test_that('the Euler plot shows the region counts and the statistics', {
    p <- EulerPlot(list(Atlas = a, Ours = b), universe = 2000)
    expect_equal(p$layers[[2]]$data$label, c('200', '100', '150'))
    expect_match(p$labels$caption, 'Jaccard index 0.22')
    expect_match(p$labels$caption, "Fisher's exact test: OR 5.16, p = 7.8e-26")
    expect_equal(levels(p$layers[[1]]$data$set), c('Atlas', 'Ours'))
    # no universe: Jaccard only; no statistics at all
    expect_no_match(EulerPlot(list(a, b))$labels$caption, 'Fisher')
    expect_null(EulerPlot(list(a, b), stats = NULL)$labels$caption)
    # empty regions are not labelled
    p <- EulerPlot(list(a, a[1:50]))
    expect_equal(p$layers[[2]]$data$label, c('250', '50'))
})

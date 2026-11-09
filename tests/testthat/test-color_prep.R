test_that('percentile strings and numbers resolve to limits', {
    v <- 0:100
    expect_equal(.ResolveLimit(v, 'q90'), 90)
    expect_equal(.ResolveLimit(v, 3), 3)
    expect_equal(.ResolveLimit(v, '3'), 3)
    expect_error(.ResolveLimit(v, 'p90'), 'percentile')
})

test_that('continuous values are clipped, low cells flagged and high values drawn last', {
    df <- data.frame(cell = letters[1:6], value = c(5, 0, 3, 10, 1, 0))
    out <- .PrepareValues(df, plot_min = 0, plot_max = 6)
    expect_true(attr(out, 'continuous'))
    expect_equal(attr(out, 'limits'), c(0, 6))
    expect_equal(max(out$value), 6)
    expect_equal(out$value, sort(out$value))
    expect_equal(sum(out$low), 2)
    kept <- .SubsetRows(out, !out$low)
    expect_equal(attr(kept, 'limits'), c(0, 6))
})

test_that('categorical values become factors with sorted levels', {
    df <- data.frame(cell = 1:4, value = c('b', 'a', 'c', 'a'))
    out <- .PrepareValues(df)
    expect_false(attr(out, 'continuous'))
    expect_equal(levels(out$value), c('a', 'b', 'c'))
})

test_that('colors come from a named vector, a data frame or the hue palette', {
    expect_equal(.AsNamedColors(data.frame(group = 'a', colour = 'red')), c(a = 'red'))
    expect_equal(names(.AsNamedColors(NULL, c('x', 'y'))), c('x', 'y'))
})

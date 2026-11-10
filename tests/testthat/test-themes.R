test_that('ThemeTextSize follows the current theme', {
    old <- ggplot2::theme_set(ThemePublication())
    on.exit(ggplot2::theme_set(old))
    expect_equal(ThemeTextSize(unit = 'pt'), 10)
    expect_equal(ThemeTextSize(), 10 / ggplot2::.pt)
    expect_equal(ThemeTextSize(0.5, unit = 'pt'), 5)
    expect_equal(ThemeTextSize(element = 'axis.title', unit = 'pt'), 11)
    expect_equal(ThemeTextSize(element = 'legend.text', unit = 'pt'), 9)
    ggplot2::theme_set(ThemePublication(base_size = 7))
    expect_equal(ThemeTextSize(unit = 'pt'), 6)
    expect_equal(ThemeTextSize(unit = 'pt', theme = ggplot2::theme_grey()), 8.8)
})

test_that('in-plot text stays at or below the tick labels and scales with the theme', {
    old <- ggplot2::theme_set(ThemePublication(base_size = 7))
    on.exit(ggplot2::theme_set(old))
    p <- EulerPlot(list(A = letters[1:10], B = letters[5:15]))
    sizes <- vapply(p$layers[2:3], function(l) l$aes_params$size, numeric(1)) * ggplot2::.pt
    expect_equal(unname(sizes), c(6, 7))  # counts at tick size, set names at axis-title size, below the 8 pt title
    p <- VolcanoPlot(sim_de[sim_de$lineage == 'Myeloid', ], label_top = 2, raster = FALSE)
    repel <- Filter(function(l) inherits(l$geom, 'GeomTextRepel'), p$layers)[[1]]
    expect_lt(repel$aes_params$size * ggplot2::.pt, 6)
})

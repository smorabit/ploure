test_that('color files round-trip and keep their order', {
    f <- tempfile(fileext = '.csv')
    cols <- c(Zeta = '#000000', Alpha = '#FFFFFF')
    WriteColors(cols, f, name_col = 'lv2')
    expect_equal(readLines(f)[1], 'lv2,hex_color')
    expect_equal(ReadColors(f), cols)
})

test_that('grouped palettes give one shade per member and one family color per parent', {
    pal <- GroupedPalette(sim_cells$patient, sim_cells$disease, family_labels = c(LC = 'Lung cancer'))
    expect_length(pal, 9)
    expect_equal(names(pal)[1:3], c('P04', 'P05', 'P06'))
    expect_equal(names(attr(pal, 'families')), c('BC', 'Lung cancer', 'OC'))
    expect_false(anyDuplicated(pal) > 0)
    ramps <- list(LC = c('white', 'darkgreen'), BC = c('white', 'darkred'), OC = c('white', 'navy'))
    pal2 <- GroupedPalette(sim_cells$patient, sim_cells$disease, ramps = ramps)
    expect_equal(unname(pal2['P04']), '#FFFFFF')
    expect_error(GroupedPalette(c('a', 'a'), c('x', 'y')), 'one parent')
})

test_that('the CVD check reports the closest pair per condition', {
    skip_if_not_installed('colorspace')
    tab <- CheckPaletteCVD(c(red = '#FF0000', green = '#00AA00', blue = '#0000FF'), plot = FALSE)
    expect_equal(tab$condition, c('Normal', 'Deuteranopia', 'Protanopia', 'Tritanopia'))
    # red and green collapse under deuteranopia
    expect_lt(tab$min_delta_e[2], tab$min_delta_e[1])
    expect_equal(tab$closest_pair[2], 'red / green')
})

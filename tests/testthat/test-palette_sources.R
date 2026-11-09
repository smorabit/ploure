test_that('palettes resolve with or without a source prefix', {
    expect_equal(PaletteColors('brewer:Set2', 3), RColorBrewer::brewer.pal(3, 'Set2'))
    expect_equal(PaletteColors('Set2', 3), RColorBrewer::brewer.pal(3, 'Set2'))
    expect_equal(PaletteColors('viridis:mako', 4), viridisLite::mako(4))
    expect_equal(PaletteColors('hcl:Zissou 1', 4), grDevices::hcl.colors(4, 'Zissou 1'))
    expect_equal(PaletteColors('viridis:mako', 4, reverse = TRUE), rev(viridisLite::mako(4)))
    expect_error(PaletteColors('nope:Set2'), 'unknown palette source')
    expect_error(PaletteColors('brewer:Nope'), 'not found')
    expect_error(PaletteColors('NotAPalette'), 'not found')
    expect_error(PaletteColors('viridis:mako'), 'needs n')
})

test_that('MetBrewer and NatParksPalettes match their own functions', {
    skip_if_not_installed('MetBrewer')
    for(n in c(3, 7)){
        expect_equal(PaletteColors('met:Hiroshige', n), as.character(MetBrewer::met.brewer('Hiroshige', n)))
    }
    expect_equal(PaletteColors('met:Hiroshige', 20, type = 'continuous'),
                 as.character(MetBrewer::met.brewer('Hiroshige', 20, type = 'continuous')))
    skip_if_not_installed('NatParksPalettes')
    expect_equal(PaletteColors('parks:Yellowstone', 4),
                 as.character(NatParksPalettes::natparks.pals('Yellowstone', 4)))
})

test_that('discrete palettes are interpolated past their length, with a message', {
    expect_message(cols <- PaletteColors('brewer:Set2', 12), 'interpolating')
    expect_length(cols, 12)
})

test_that('ColorMap names colors by group level', {
    cm <- ColorMap(sim_cells$cluster, 'brewer:Paired')
    expect_equal(names(cm), levels(sim_cells$cluster))
    expect_equal(unname(cm), RColorBrewer::brewer.pal(9, 'Paired'))
})

test_that('palette names work wherever colors are taken', {
    expect_true(.IsPaletteName('viridis:mako'))
    expect_true(.IsPaletteName('Set2'))
    expect_false(.IsPaletteName('red'))
    expect_false(.IsPaletteName(c(a = 'red')))
    expect_equal(.AsNamedColors('brewer:Set2', c('a', 'b', 'c')), c(a = '#66C2A5', b = '#FC8D62', c = '#8DA0CB'))
    expect_equal(.ContinuousPalette('viridis:mako', rev_colors = FALSE),
                 PaletteColors('viridis:mako', 256, type = 'continuous'))
    p <- PlotEmbedding(sim_cells, group.by = 'cluster', reduction = 'umap', color_df = 'brewer:Paired',
                       label = FALSE, raster = FALSE)
    built <- ggplot2::ggplot_build(p)$data[[1]]
    expect_setequal(unique(built$colour), RColorBrewer::brewer.pal(9, 'Paired'))
    pal <- GroupedPalette(sim_cells$patient, sim_cells$disease,
                          ramps = list(LC = 'brewer:Greens', BC = 'brewer:Reds', OC = 'brewer:Blues'))
    expect_length(pal, 9)
})

test_that('palettes can be listed and filtered', {
    tab <- ListPalettes(c('brewer', 'viridis'))
    expect_named(tab, c('palette', 'source', 'name', 'n', 'type', 'colorblind'))
    expect_true('brewer:Set2' %in% tab$palette)
    expect_false('viridis:turbo' %in% ListPalettes('viridis', colorblind = TRUE)$palette)
    # every listed palette resolves
    for(cur in ListPalettes()$palette) expect_type(PaletteColors(cur, 5), 'character')
})

test_that('ggsci palettes match the package, by family or family_variant', {
    skip_if_not_installed('ggsci')
    # ggsci's functions append an alpha channel ('#E64B35FF')
    rgb <- function(x) toupper(substr(unname(x), 1, 7))
    expect_equal(rgb(PaletteColors('ggsci:npg', 5)), rgb(ggsci::pal_npg()(5)))
    expect_equal(rgb(PaletteColors('ggsci:igv', 40)), rgb(ggsci::pal_igv()(40)))
    expect_equal(rgb(PaletteColors('ggsci:d3_category20', 20)), rgb(ggsci::pal_d3('category20')(20)))
    expect_equal(rgb(PaletteColors('ggsci:d3', 3)), rgb(ggsci::pal_d3('category10')(3)))
})

test_that('Polychrome palettes match the package, without near-white colors', {
    skip_if_not_installed('Polychrome')
    glasbey <- unname(Polychrome::glasbey.colors())
    cols <- PaletteColors('polychrome:glasbey')
    expect_false('#FFFFFF' %in% cols)
    expect_equal(cols, glasbey[glasbey %in% cols])
    expect_equal(PaletteColors('polychrome:palette36', 10), unname(Polychrome::palette36.colors())[1:10])
})

test_that('scico palettes match the package', {
    skip_if_not_installed('scico')
    expect_equal(PaletteColors('scico:batlow', 6), scico::scico(6, palette = 'batlow'))
    expect_equal(PaletteColors('vik', 6), scico::scico(6, palette = 'vik'))
    tab <- ListPalettes('scico')
    expect_equal(tab$type[tab$name == 'vik'], 'diverging')
    expect_equal(tab$type[tab$name == 'romaO'], 'cyclic')
})

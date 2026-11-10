test_that('the figure font is Arial or Helvetica', {
    expect_true(FigureFont() %in% c('Arial', 'Helvetica'))
    expect_type(.FcFamilies(), 'character')
})

test_that('SetFigureFont sets the theme, geom and option', {
    old_theme <- ggplot2::theme_get()
    old_opt <- options(figure_font = NULL)
    on.exit({ggplot2::theme_set(old_theme); options(old_opt)})
    ggplot2::theme_set(ThemePublication())
    expect_warning(SetFigureFont('Helvetica'), 'Arial is not installed')
    expect_equal(getOption('figure_font'), 'Helvetica')
    expect_equal(ggplot2::theme_get()$text$family, 'Helvetica')
    if(utils::packageVersion('ggplot2') >= '4.0.0') expect_equal(ggplot2::theme_get()$geom$family, 'Helvetica')
    # FigurePDF opens pdf() with Helvetica
    f <- tempfile(fileext = '.pdf')
    FigurePDF(f, 3, 3)
    print(ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0, label = 'x'))
    grDevices::dev.off()
    expect_true(file.exists(f))
    skip_if(Sys.which('pdffonts') == '', 'pdffonts not installed')
    expect_true(all(CheckFigureFonts(f)$ok))
})

test_that('CheckFigureFonts flags other fonts', {
    skip_if(Sys.which('pdffonts') == '' || !capabilities('cairo'), 'pdffonts or cairo not available')
    f <- tempfile(fileext = '.pdf')
    grDevices::cairo_pdf(f, 3, 3, family = 'DejaVu Sans')
    grid::grid.text('x')
    grDevices::dev.off()
    res <- suppressWarnings(CheckFigureFonts(f))
    skip_if(all(res$ok), 'DejaVu Sans not installed')
    expect_warning(CheckFigureFonts(f), 'other than Arial')
})

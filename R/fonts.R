#---------------------------------------------------------#
# Figure font: Arial (embedded, via cairo_pdf) when installed, R's Helvetica (pdf) as the fallback
#
# Ported from PLEURA's utils/plotting_functions.R. Journals (Nature, Cell Press) accept only Arial or Helvetica;
# look-alikes (Liberation Sans, Arimo, DejaVu Sans) are named as such in the PDF. Arial is licensed: point fontconfig
# at a licensed copy (e.g. a fonts.conf with a <dir> for it, and FONTCONFIG_FILE set in the pixi environment) and never
# commit it. Usage:
#
#     theme_set(ThemePublication())
#     SetFigureFont()
#     FigurePDF('figures/x.pdf', width = 8, height = 6)
#     print(p)
#     dev.off()
#---------------------------------------------------------#

#' The font family for figures
#'
#' `'Arial'` when a real Arial is installed (as seen by fontconfig, through systemfonts), otherwise `'Helvetica'`.
#'
#' @return `'Arial'` or `'Helvetica'`.
#' @export
#' @examples
#' FigureFont()
FigureFont <- function(){
    families <- if(requireNamespace('systemfonts', quietly = TRUE)){
        systemfonts::system_fonts()$family
    } else {
        .FcFamilies()
    }
    if('Arial' %in% families) 'Arial' else 'Helvetica'
}

#' Set the figure font
#'
#' Sets the font for theme text and for text drawn by geoms (`geom_text()`, `annotate()`, ggrepel labels) in the current
#' ggplot theme, and the null device cowplot uses to measure extracted legends. Call it after `theme_set()`. The family
#' is kept in `options(figure_font)` for [FigurePDF()]. Warns when Arial is not installed.
#'
#' @param family Font family (default [FigureFont()]).
#'
#' @return The family, invisibly.
#' @export
#' @examples
#' \dontrun{
#' theme_set(ThemePublication())
#' SetFigureFont()
#' }
SetFigureFont <- function(family = FigureFont()){
    if(family != 'Arial'){
        warning('Arial is not installed: figures fall back to Helvetica (pdf(), not embedded; \u2265, \u2212 and Greek ',
                'letters are substituted). Point fontconfig at a licensed Arial (see ?FigurePDF).', call. = FALSE)
    }
    options(figure_font = family)
    theme_update(text = element_text(family = family))
    if(utils::packageVersion('ggplot2') >= '4.0.0'){
        theme_update(geom = do.call(ggplot2::element_geom, list(family = family)))
    } else {
        for(g in c('text', 'label')) update_geom_defaults(g, list(family = family))
    }
    # ggrepel labels do not read the theme's geom defaults
    if(requireNamespace('ggrepel', quietly = TRUE)){
        for(g in c(ggrepel::GeomTextRepel, ggrepel::GeomLabelRepel)) update_geom_defaults(g, list(family = family))
    }
    if(requireNamespace('cowplot', quietly = TRUE)){
        cowplot::set_null_device(if(family == 'Arial') 'cairo' else 'pdf')
    }
    invisible(family)
}

#' Open a PDF device with the figure font
#'
#' Opens `cairo_pdf()` (embeds Arial) or `pdf()` (names Helvetica), following [SetFigureFont()]. Close it with
#' `dev.off()`. Also usable as a ggsave() device (`ggsave(..., device = FigurePDF)`).
#'
#' Arial is licensed, so it is not shipped: make a licensed copy visible to fontconfig, for example with a `fonts.conf`
#' that adds its folder (`<dir>`) and maps Helvetica and sans to Arial, and `FONTCONFIG_FILE` set in the pixi
#' environment's `[activation.env]`. Check the result with `pdffonts` ([CheckFigureFonts()]).
#'
#' @param filename Output path.
#' @param width,height Size in inches.
#' @param ... Passed to `cairo_pdf()` or `pdf()`.
#'
#' @return Opens the device (invisibly NULL).
#' @export
#' @examples
#' \dontrun{
#' FigurePDF('figure.pdf', width = 6, height = 4)
#' print(p)
#' dev.off()
#' }
FigurePDF <- function(filename, width, height, ...){
    family <- getOption('figure_font', FigureFont())
    if(family == 'Arial'){
        grDevices::cairo_pdf(filename, width = width, height = height, family = 'Arial', ...)
    } else {
        grDevices::pdf(filename, width = width, height = height, family = 'Helvetica', ...)
    }
    invisible(NULL)
}

#' Check the fonts of figure PDFs
#'
#' Lists the fonts in each PDF with `pdffonts` (poppler) and flags any that are not Arial or Helvetica, such as DejaVu
#' Sans from `cairo_pdf()` without a family, or a look-alike like Liberation Sans.
#'
#' @param files PDF files.
#' @param pdffonts The `pdffonts` command (e.g. `'pixi exec --spec poppler pdffonts'` when poppler is not installed).
#'
#' @return A data frame with one row per file and font: `file`, `font`, `ok`; invisibly, after printing a warning for
#'   the files with other fonts.
#' @export
CheckFigureFonts <- function(files, pdffonts = 'pdffonts'){
    rows <- lapply(files, function(f){
        out <- suppressWarnings(system(paste(pdffonts, shQuote(f)), intern = TRUE))
        if(length(out) < 2 || !grepl('^name', out[1])) stop('could not run pdffonts on ', f, ': ', paste(out, collapse = ' '))
        fonts <- sub('^(\\S+).*', '\\1', out[-(1:2)])
        fonts <- sub('^[A-Z]{6}\\+', '', fonts)  # subset prefix, e.g. ABCDEF+ArialMT
        if(length(fonts) == 0) fonts <- NA_character_
        data.frame(file = f, font = fonts, ok = is.na(fonts) | grepl('^(Arial|Helvetica)', fonts))
    })
    res <- do.call(rbind, rows)
    bad <- res[!res$ok, ]
    if(nrow(bad) > 0){
        warning('fonts other than Arial / Helvetica: ', paste(unique(paste0(basename(bad$file), ' (', bad$font, ')')),
                                                              collapse = ', '), call. = FALSE)
    }
    invisible(res)
}

# font families known to fontconfig, without systemfonts
.FcFamilies <- function(){
    out <- tryCatch(suppressWarnings(system('fc-list : family', intern = TRUE)), error = function(e) character(0))
    unique(trimws(unlist(strsplit(out, ','))))
}

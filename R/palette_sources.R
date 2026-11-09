#---------------------------------------------------------#
# Palettes from other packages
#
# One adapter over several palette sources, so any ploure argument that takes colors also takes a palette name:
# - met: MetBrewer (palettes from the Metropolitan Museum of Art)
# - parks: NatParksPalettes (US national parks)
# - brewer: RColorBrewer
# - ggsci: journal and other palettes (npg, aaas, nejm, lancet, jama, igv, d3_category20, ...)
# - polychrome: Polychrome palettes of many distinguishable colors (glasbey, palette36, alphabet, kelly, ...);
#   near-white colors are dropped, since they vanish on a white background
# - scico: Crameri's perceptually uniform scientific color maps (batlow, roma, vik, ...)
# - viridis: viridisLite (magma, inferno, plasma, viridis, cividis, rocket, mako, turbo)
# - hcl: grDevices::hcl.colors (base R)
# A palette is named 'source:name' (e.g. 'met:Hiroshige'), or by its name alone, which is looked up in that order.
#---------------------------------------------------------#

.palette_sources <- c(met = 'MetBrewer', parks = 'NatParksPalettes', brewer = 'RColorBrewer', ggsci = 'ggsci',
                      polychrome = 'Polychrome', scico = 'scico', viridis = 'viridisLite', hcl = 'grDevices')

.polychrome_palettes <- c('glasbey', 'palette36', 'alphabet', 'green.armytage', 'kelly', 'dark', 'light', 'sky')
.scico_diverging <- c('bam', 'berlin', 'broc', 'cork', 'lisbon', 'managua', 'roma', 'tofino', 'vanimo', 'vik')

#' Colors from a named palette
#'
#' Takes colors from MetBrewer, NatParksPalettes, RColorBrewer, ggsci, Polychrome, scico, viridisLite or base R HCL
#' palettes through one interface. Polychrome palettes lose their near-white colors (e.g. the white that starts
#' `glasbey`), which vanish on a white background. Discrete colors follow each palette's own order (as `MetBrewer::met.brewer()` does); when more colors
#' are asked for than a palette has, or for `type = 'continuous'`, the palette is interpolated.
#'
#' @param palette `'source:name'` (e.g. `'met:Hiroshige'`, `'parks:Yellowstone'`, `'brewer:Set2'`, `'ggsci:npg'`,
#'   `'ggsci:igv'`, `'polychrome:glasbey'`, `'scico:batlow'`, `'viridis:mako'`, `'hcl:Zissou 1'`) or a name alone, looked up in the sources in that order. See
#'   [ListPalettes()].
#' @param n Number of colors; default the palette's own number (required for scico, viridis and hcl).
#' @param type `'discrete'` (the palette's colors, in its preferred order) or `'continuous'` (interpolated along the
#'   palette).
#' @param reverse Reverse the colors.
#'
#' @return A character vector of hex colors.
#' @export
#' @examples
#' PaletteColors('brewer:Set2', 4)
#' PaletteColors('viridis:mako', 5)
#' if(requireNamespace('MetBrewer', quietly = TRUE)) PaletteColors('met:Hiroshige', 6)
PaletteColors <- function(palette, n = NULL, type = c('discrete', 'continuous'), reverse = FALSE){
    type <- match.arg(type)
    pal <- .ResolvePalette(palette)
    base <- pal$colors
    if(is.null(base)){
        # generated palettes (viridis, hcl) are drawn at n directly
        if(is.null(n)) stop("palette '", palette, "' needs n.")
        out <- pal$fun(n)
    } else if(type == 'continuous'){
        if(is.null(n)) n <- length(base)
        out <- grDevices::colorRampPalette(base)(n)
    } else {
        if(is.null(n)) n <- length(base)
        if(n <= length(base)){
            out <- base[pal$order <= n]
        } else {
            message("palette '", palette, "' has ", length(base), ' colors; interpolating ', n, '.')
            out <- grDevices::colorRampPalette(base)(n)
        }
    }
    if(reverse) out <- rev(out)
    out
}

#' Named colors for a set of groups from a palette
#'
#' @param values Group values (e.g. a cluster column) or the groups themselves; one color per level (factor levels,
#'   else sorted unique values).
#' @param palette A palette name, as in [PaletteColors()].
#' @param reverse Reverse the colors.
#'
#' @return A named color vector, for any colors argument (e.g. `color_df` in [PlotEmbedding()]).
#' @export
#' @examples
#' ColorMap(sim_cells$cluster, 'brewer:Paired')
ColorMap <- function(values, palette, reverse = FALSE){
    lv <- levels(.AsFactor(stats::na.omit(values)))
    stats::setNames(PaletteColors(palette, length(lv), reverse = reverse), lv)
}

#' A palette as a continuous color function
#'
#' @inheritParams PaletteColors
#'
#' @return A function of `n` returning `n` colors along the palette, e.g. for `colfunc` in [FeatureEmbedding()] (whose
#'   `rev_colors = TRUE` reverses it again).
#' @export
#' @examples
#' f <- PaletteFunction('viridis:rocket')
#' f(5)
PaletteFunction <- function(palette, reverse = FALSE){
    .ResolvePalette(palette)
    function(n) PaletteColors(palette, n, type = 'continuous', reverse = reverse)
}

#' List the available palettes
#'
#' @param source Optional sources to list (`'met'`, `'parks'`, `'brewer'`, `'ggsci'`, `'polychrome'`, `'scico'`,
#'   `'viridis'`, `'hcl'`); sources whose package is not installed are left out.
#' @param colorblind If TRUE, only palettes flagged as colorblind-friendly by their package (MetBrewer,
#'   NatParksPalettes, RColorBrewer) or designed to be (scico, viridis except turbo).
#'
#' @return A data frame with `palette` (the full name to use), `source`, `name`, `n` (colors; NA for generated
#'   palettes), `type` (qualitative, sequential, diverging, cyclic, where known) and `colorblind`.
#' @export
#' @examples
#' head(ListPalettes('brewer'))
ListPalettes <- function(source = NULL, colorblind = NULL){
    if(is.null(source)) source <- names(.palette_sources)
    out <- do.call(rbind, lapply(source, .SourceTable))
    if(is.null(out)) return(data.frame())
    out$palette <- paste0(out$source, ':', out$name)
    out <- out[c('palette', 'source', 'name', 'n', 'type', 'colorblind')]
    if(isTRUE(colorblind)) out <- out[out$colorblind %in% TRUE, ]
    rownames(out) <- NULL
    out
}

# is a single string a palette name (rather than a color)?
.IsPaletteName <- function(x){
    if(!is.character(x) || length(x) != 1 || !is.null(names(x))) return(FALSE)
    if(grepl(':', x, fixed = TRUE)) return(TRUE)
    !is.null(tryCatch(.ResolvePalette(x), error = function(e) NULL))
}

# a palette: colors + order (fixed palettes), or fun (generated palettes)
.ResolvePalette <- function(palette){
    if(grepl(':', palette, fixed = TRUE)){
        source <- sub(':.*', '', palette)
        name <- sub('^[^:]*:', '', palette)
        if(!(source %in% names(.palette_sources))){
            stop("unknown palette source '", source, "'; use one of: ", paste(names(.palette_sources), collapse = ', '))
        }
        pal <- .SourcePalette(source, name)
        if(is.null(pal)){
            pkg <- .palette_sources[[source]]
            if(!requireNamespace(pkg, quietly = TRUE)) stop("palette source '", source, "' needs the ", pkg, ' package.')
            stop("palette '", name, "' not found in ", pkg, '.')
        }
        return(pal)
    }
    for(source in names(.palette_sources)){
        pal <- .SourcePalette(source, palette)
        if(!is.null(pal)) return(pal)
    }
    stop("palette '", palette, "' not found (see ListPalettes()).")
}

.SourcePalette <- function(source, name){
    pkg <- .palette_sources[[source]]
    if(source != 'hcl' && !requireNamespace(pkg, quietly = TRUE)) return(NULL)
    if(source %in% c('met', 'parks')){
        pals <- if(source == 'met') MetBrewer::MetPalettes else NatParksPalettes::NatParksPalettes
        if(!(name %in% names(pals))) return(NULL)
        return(list(colors = pals[[name]][[1]], order = pals[[name]][[2]]))
    }
    if(source == 'brewer'){
        info <- RColorBrewer::brewer.pal.info
        if(!(name %in% rownames(info))) return(NULL)
        n_max <- info[name, 'maxcolors']
        return(list(colors = RColorBrewer::brewer.pal(n_max, name), order = seq_len(n_max)))
    }
    if(source == 'ggsci'){
        tab <- .GgsciTable()
        i <- match(name, tab$name)
        # a family name alone takes its first variant (e.g. 'd3' -> 'd3_category10')
        if(is.na(i)) i <- match(name, tab$family)
        if(is.na(i)) return(NULL)
        cols <- tab$colors[[i]]
        return(list(colors = cols, order = seq_along(cols)))
    }
    if(source == 'polychrome'){
        if(!(name %in% .polychrome_palettes)) return(NULL)
        cols <- .DropNearWhite(unname(getExportedValue('Polychrome', paste0(name, '.colors'))()))
        return(list(colors = cols, order = seq_along(cols)))
    }
    if(source == 'scico'){
        if(!(name %in% scico::scico_palette_names())) return(NULL)
        return(list(fun = function(n) scico::scico(n, palette = name)))
    }
    if(source == 'viridis'){
        opts <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo')
        if(!(name %in% opts)) return(NULL)
        return(list(fun = function(n) viridisLite::viridis(n, option = name)))
    }
    if(!(name %in% grDevices::hcl.pals())) return(NULL)
    list(fun = function(n) grDevices::hcl.colors(n, name))
}

.SourceTable <- function(source){
    pkg <- .palette_sources[[source]]
    if(source != 'hcl' && !requireNamespace(pkg, quietly = TRUE)) return(NULL)
    if(source %in% c('met', 'parks')){
        pals <- if(source == 'met') MetBrewer::MetPalettes else NatParksPalettes::NatParksPalettes
        return(data.frame(source = source, name = names(pals), n = sapply(pals, function(p) length(p[[1]])),
                          type = NA_character_, colorblind = sapply(pals, function(p) isTRUE(p$colorblind))))
    }
    if(source == 'brewer'){
        info <- RColorBrewer::brewer.pal.info
        types <- c(qual = 'qualitative', seq = 'sequential', div = 'diverging')
        return(data.frame(source = source, name = rownames(info), n = info$maxcolors,
                          type = unname(types[as.character(info$category)]), colorblind = info$colorblind))
    }
    if(source == 'ggsci'){
        tab <- .GgsciTable()
        return(data.frame(source = source, name = tab$name, n = lengths(tab$colors),
                          type = ifelse(tab$family == 'gsea', 'diverging', NA_character_), colorblind = NA))
    }
    if(source == 'polychrome'){
        return(data.frame(source = source, name = .polychrome_palettes,
                          n = sapply(.polychrome_palettes, function(x) length(.SourcePalette('polychrome', x)$colors)),
                          type = 'qualitative', colorblind = FALSE))
    }
    if(source == 'scico'){
        opts <- scico::scico_palette_names()
        type <- ifelse(grepl('O$', opts), 'cyclic', ifelse(opts %in% .scico_diverging, 'diverging', 'sequential'))
        return(data.frame(source = source, name = opts, n = NA_integer_, type = type, colorblind = TRUE))
    }
    if(source == 'viridis'){
        opts <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo')
        return(data.frame(source = source, name = opts, n = NA_integer_, type = 'sequential',
                          colorblind = opts != 'turbo'))
    }
    hcl_types <- c(qualitative = 'qualitative', sequential = 'sequential', diverging = 'diverging',
                   divergingx = 'diverging')
    types <- unlist(lapply(names(hcl_types), function(t) stats::setNames(rep(hcl_types[[t]],
                                                                             length(grDevices::hcl.pals(t))),
                                                                         grDevices::hcl.pals(t))))
    data.frame(source = source, name = grDevices::hcl.pals(), n = NA_integer_,
               type = unname(types[grDevices::hcl.pals()]), colorblind = NA)
}

# the fixed ggsci palettes: one row per family and variant, named 'family' (one variant or 'default') or 'family_variant'
.GgsciTable <- function(){
    db <- utils::getFromNamespace('ggsci_db', 'ggsci')
    do.call(rbind, lapply(names(db), function(family){
        variants <- names(db[[family]])
        name <- ifelse(variants == 'default' | length(variants) == 1, family, paste0(family, '_', variants))
        out <- data.frame(family = family, name = name)
        out$colors <- lapply(variants, function(v) unname(db[[family]][[v]]))
        # fixed palettes only (gephi stores generator settings, not colors)
        out[sapply(out$colors, is.character), ]
    }))
}

# drop near-white colors (CIELAB lightness above 93)
.DropNearWhite <- function(colors){
    lab <- grDevices::convertColor(t(grDevices::col2rgb(colors)) / 255, from = 'sRGB', to = 'Lab')
    colors[lab[, 1] <= 93]
}

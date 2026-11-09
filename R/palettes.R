#---------------------------------------------------------#
# Palettes: named color files, grouped (hue family) palettes and a color-vision-deficiency check
#---------------------------------------------------------#

#' Read and write named color files
#'
#' Color files are two-column CSVs: the group names, then the hex colors (e.g. `lv2,hex_color`). The file order is
#' kept, so it can also set the factor levels.
#'
#' @param file Path to the CSV.
#' @param colors Named color vector.
#' @param name_col,color_col Column names to write.
#'
#' @return `ReadColors()`: a named color vector. `WriteColors()`: the file path, invisibly.
#' @export
#' @examples
#' f <- tempfile(fileext = '.csv')
#' WriteColors(c(Lymphoid = '#1f78b4', Myeloid = '#e31a1c'), f, name_col = 'lineage')
#' ReadColors(f)
ReadColors <- function(file){
    df <- utils::read.csv(file, colClasses = 'character', check.names = FALSE)
    if(ncol(df) < 2) stop('a color file needs two columns: names, then colors.')
    stats::setNames(df[[2]], df[[1]])
}

#' @rdname ReadColors
#' @export
WriteColors <- function(colors, file, name_col = 'group', color_col = 'hex_color'){
    colors <- .AsNamedColors(colors)
    df <- stats::setNames(data.frame(names(colors), unname(colors)), c(name_col, color_col))
    utils::write.csv(df, file, row.names = FALSE, quote = FALSE)
    invisible(file)
}

#' Grouped palette: one hue family per parent group
#'
#' Colors members (e.g. patients) so that their parent group (e.g. cancer type) reads from the hue and the member from
#' the shade: each parent gets a color ramp, and its members take evenly spaced shades along it.
#'
#' @param members Member names (e.g. patients), one per member, or one per cell (duplicates are dropped).
#' @param parents Parent group of each entry of `members`, same length.
#' @param ramps Optional named list of color ramps, one per parent (e.g. `list(LC = c('#c7e9c0', '#00441b'))`).
#'   Default: one light-to-dark HCL ramp per parent, hues evenly spaced.
#' @param family_labels Optional named vector of legend labels for the parents (e.g. `c(LC = 'Lung cancer')`).
#'
#' @return A named color vector over the members, ordered by parent then member. The attribute `families` holds one
#'   representative color per parent (the middle of its ramp), for a compact legend (e.g. `legend_colors` in
#'   [ZoomEmbedding()]).
#' @export
#' @examples
#' pal <- GroupedPalette(sim_cells$patient, sim_cells$disease)
#' pal
#' attr(pal, 'families')
GroupedPalette <- function(members, parents, ramps = NULL, family_labels = NULL){
    if(length(members) != length(parents)) stop('members and parents must have the same length.')
    pairs <- unique(data.frame(member = as.character(members), parent = parents))
    if(anyDuplicated(pairs$member)) stop('each member must belong to exactly one parent.')
    pairs$parent <- .AsFactor(pairs$parent)
    pairs <- pairs[order(pairs$parent, pairs$member), ]
    lv <- levels(droplevels(pairs$parent))
    if(is.null(ramps)){
        hues <- seq(15, 375, length.out = length(lv) + 1)[seq_along(lv)]
        ramps <- stats::setNames(lapply(hues, function(h) grDevices::hcl(h, c = c(45, 70, 50), l = c(82, 55, 30))),
                                 lv)
    }
    missing <- setdiff(lv, names(ramps))
    if(length(missing) > 0) stop('no ramp for: ', paste(missing, collapse = ', '))
    colors <- unlist(lapply(lv, function(cur_parent){
        cur <- pairs$member[pairs$parent == cur_parent]
        stats::setNames(grDevices::colorRampPalette(ramps[[cur_parent]])(length(cur)), cur)
    }))
    families <- sapply(lv, function(cur_parent) grDevices::colorRampPalette(ramps[[cur_parent]])(5)[3])
    names(families) <- .MapLabels(lv, family_labels)
    attr(colors, 'families') <- families
    colors
}

#' Check a palette for color-vision deficiencies
#'
#' Simulates deuteranopia, protanopia and tritanopia (colorspace) and measures how distinguishable the colors remain:
#' the smallest pairwise CIELAB distance (Delta E 1976) under each condition. Pairs below about 10 are hard to tell
#' apart.
#'
#' @param colors Named (or unnamed) color vector.
#' @param plot If TRUE, return swatches of the palette under each condition; if FALSE, a table.
#' @param severity Simulation severity (0-1).
#'
#' @return A ggplot of swatches (with the closest pair per condition in the row labels), or a data frame with one row
#'   per condition: `condition`, `min_delta_e`, `closest_pair`.
#' @export
#' @examples
#' CheckPaletteCVD(c(a = '#1b9e77', b = '#d95f02', c = '#7570b3'), plot = FALSE)
CheckPaletteCVD <- function(colors, plot = TRUE, severity = 1){
    if(!requireNamespace('colorspace', quietly = TRUE)) stop('CheckPaletteCVD needs the colorspace package.')
    colors <- .AsNamedColors(colors)
    if(is.null(names(colors))) names(colors) <- colors
    sims <- list(
        Normal = colors,
        Deuteranopia = colorspace::deutan(colors, severity),
        Protanopia = colorspace::protan(colors, severity),
        Tritanopia = colorspace::tritan(colors, severity)
    )
    summary_df <- do.call(rbind, lapply(names(sims), function(cond){
        lab <- grDevices::convertColor(t(grDevices::col2rgb(sims[[cond]])) / 255, from = 'sRGB', to = 'Lab')
        d <- as.matrix(stats::dist(lab))
        diag(d) <- Inf
        ij <- sort(which(d == min(d), arr.ind = TRUE)[1, ])
        data.frame(condition = cond, min_delta_e = round(min(d), 1),
                   closest_pair = paste(names(colors)[ij], collapse = ' / '))
    }))
    if(!plot) return(summary_df)

    row_labels <- sprintf('%s\nmin. \u0394E %s (%s)', summary_df$condition, summary_df$min_delta_e,
                          summary_df$closest_pair)
    swatch_df <- do.call(rbind, lapply(seq_along(sims), function(i){
        data.frame(condition = row_labels[i], color = names(colors), hex = unname(sims[[i]]))
    }))
    swatch_df$condition <- factor(swatch_df$condition, levels = rev(row_labels))
    swatch_df$color <- factor(swatch_df$color, levels = names(colors))
    ggplot(swatch_df, aes(x = color, y = condition, fill = hex)) +
        geom_tile(color = 'white', linewidth = 0.5) +
        scale_fill_identity() +
        labs(x = NULL, y = NULL) +
        ThemePublication() +
        theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        )
}

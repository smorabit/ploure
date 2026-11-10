#---------------------------------------------------------#
# Two-set Euler diagram with overlap statistics
#
# Generalized from PLEURA's v2/comparison/compare_v1_v2_figure.Rmd (panel b): each circle's area is proportional to
# the size of its set, and the distance between the centres is solved numerically so that the lens area equals the
# shared elements. Plain ggplot (no eulerr / VennDiagram). Below the diagram, the overlap statistics: the Jaccard
# index and, given a universe, Fisher's exact test.
#---------------------------------------------------------#

#' Overlap statistics of two sets
#'
#' Sizes, Jaccard index, overlap coefficient and, given the universe (all elements that could be in either set),
#' Fisher's exact test of the overlap.
#'
#' @param a,b Two vectors (the sets; duplicates are dropped).
#' @param counts Instead of `a` and `b`: the region counts, `c(a_only, shared, b_only)`, or a vector named with
#'   `'shared'` and two `'<set>_only'` names (e.g. `c(v1_only = 46541, v2_only = 81540, shared = 200851)`).
#' @param universe Optional universe: a vector of all elements (with `a` and `b`) or its size. Needed for Fisher's
#'   test.
#' @param alternative Alternative hypothesis of Fisher's test: `'greater'` (more overlap than expected by chance,
#'   the default), `'two.sided'` or `'less'`.
#'
#' @return A one-row data frame: `n_a`, `n_b`, `a_only`, `b_only`, `shared`, `union`, `universe`, `jaccard`
#'   (shared / union), `overlap_coef` (shared / smaller set), `expected` (shared elements expected by chance),
#'   `odds_ratio` and `p_value` (Fisher's test; NA without a universe).
#' @export
#' @examples
#' SetOverlapStats(paste0('gene', 1:300), paste0('gene', 201:450), universe = 2000)
SetOverlapStats <- function(a = NULL, b = NULL, counts = NULL, universe = NULL, alternative = 'greater'){
    alternative <- match.arg(alternative, c('greater', 'two.sided', 'less'))
    if(is.null(counts)){
        if(is.null(a) || is.null(b)) stop('give two sets (a, b) or the region counts.')
        a <- unique(a)
        b <- unique(b)
        if(length(universe) > 1){
            universe <- unique(universe)
            outside <- setdiff(union(a, b), universe)
            if(length(outside) > 0) stop(length(outside), ' set elements are not in the universe.')
        }
        counts <- c(length(setdiff(a, b)), length(intersect(a, b)), length(setdiff(b, a)))
    } else {
        counts <- .RegionCounts(counts)
    }
    if(length(universe) > 1) universe <- length(universe)
    a_only <- counts[[1]]
    shared <- counts[[2]]
    b_only <- counts[[3]]
    n_a <- a_only + shared
    n_b <- b_only + shared
    union <- a_only + shared + b_only
    if(n_a == 0 || n_b == 0) stop('both sets need at least one element.')

    out <- data.frame(n_a = n_a, n_b = n_b, a_only = a_only, b_only = b_only, shared = shared, union = union,
                      universe = if(is.null(universe)) NA_real_ else universe, jaccard = shared / union,
                      overlap_coef = shared / min(n_a, n_b), expected = NA_real_, odds_ratio = NA_real_,
                      p_value = NA_real_)
    if(!is.null(universe)){
        if(universe < union) stop('the universe (', universe, ') is smaller than the union of the sets (', union, ').')
        tab <- matrix(c(shared, b_only, a_only, universe - union), nrow = 2)
        test <- stats::fisher.test(tab, alternative = alternative)
        out$expected <- n_a * n_b / universe
        out$odds_ratio <- unname(test$estimate)
        out$p_value <- test$p.value
    }
    out
}

#' Two-set Euler diagram
#'
#' Two circles with areas proportional to the set sizes and an overlap proportional to the shared elements, the
#' region counts inside, the set names and sizes above, and the overlap statistics below (Jaccard index and, given a
#' universe, Fisher's exact test; see [SetOverlapStats()]).
#'
#' @inheritParams SetOverlapStats
#' @param sets A list of two vectors; its names are the set names.
#' @param set_names Set names (default: the names of `sets` or `counts`, else 'A' and 'B').
#' @param colors Two colors (fill and outline), named by set or in set order, or a palette name.
#' @param alpha,linewidth Circle fill transparency and outline width.
#' @param stats Statistics to show below the diagram: any of `'jaccard'`, `'overlap'` (overlap coefficient) and
#'   `'fisher'` (odds ratio and p value; only with a `universe`). `NULL` shows none.
#' @param count_size,label_size,stats_size Text sizes (mm) of the region counts, set labels and statistics (default:
#'   the current theme's tick-label size, its axis-title size and 0.9 of the tick-label size; see [ThemeTextSize()]).
#' @param title Plot title.
#'
#' @return A ggplot, with the statistics as the attribute `stats` ([SetOverlapStats()]).
#' @export
#' @examples
#' sets <- list(Atlas = paste0('gene', 1:300), Ours = paste0('gene', 201:450))
#' EulerPlot(sets, universe = 2000, title = 'Marker genes')
#' EulerPlot(counts = c(v1_only = 46541, v2_only = 81540, shared = 200851), stats = 'jaccard')
EulerPlot <- function(
    sets = NULL,
    counts = NULL,
    set_names = NULL,
    universe = NULL,
    alternative = 'greater',
    colors = c('grey65', '#2c7fb8'),
    alpha = 0.45,
    linewidth = 0.5,
    stats = c('jaccard', 'fisher'),
    count_size = ThemeTextSize(),
    label_size = ThemeTextSize(element = 'axis.title'),
    stats_size = ThemeTextSize(0.9),
    title = NULL
){
    if(!is.null(sets)){
        if(length(sets) != 2) stop('EulerPlot draws two sets.')
        if(is.null(set_names)) set_names <- names(sets)
        ov <- SetOverlapStats(sets[[1]], sets[[2]], universe = universe, alternative = alternative)
    } else {
        if(is.null(set_names)) set_names <- attr(.RegionCounts(counts), 'set_names')
        ov <- SetOverlapStats(counts = counts, universe = universe, alternative = alternative)
    }
    if(is.null(set_names) || length(set_names) != 2) set_names <- c('A', 'B')
    stats <- if(length(stats) == 0) character(0) else match.arg(stats, c('jaccard', 'overlap', 'fisher'),
                                                               several.ok = TRUE)

    geo <- .EulerGeometry(ov$n_a, ov$n_b, ov$shared)
    r1 <- geo$r1
    r2 <- geo$r2
    d <- geo$d
    circles <- rbind(.CircleDf(0, r1), .CircleDf(d, r2))
    circles$set <- factor(rep(set_names, each = nrow(circles) / 2), levels = set_names)
    regions <- data.frame(x = c((-r1 + d - r2) / 2, (d - r2 + r1) / 2, (r1 + d + r2) / 2), y = 0,
                          n = c(ov$a_only, ov$shared, ov$b_only))
    regions <- regions[regions$n > 0, ]
    regions$label <- scales::label_comma()(regions$n)
    # set names above the outer upper edge of each circle, clear of the region counts
    set_labels <- data.frame(x = c(-r1 * 0.75, d + r2 * 0.75), y = c(r1, r2) * 0.85, hjust = c(1, 0),
                             label = sprintf('%s (%s)', set_names, scales::label_comma()(c(ov$n_a, ov$n_b))))

    colors <- .AsNamedColors(colors, set_names)
    if(is.null(names(colors))) names(colors) <- set_names
    p <- ggplot() +
        geom_polygon(data = circles, aes(x = x, y = y, fill = set, color = set), alpha = alpha, linewidth = linewidth) +
        geom_text(data = regions, aes(x = x, y = y, label = label), size = count_size) +
        geom_text(data = set_labels, aes(x = x, y = y, label = label, hjust = hjust), size = label_size,
                  lineheight = 0.9) +
        scale_fill_manual(values = colors, guide = 'none') +
        scale_color_manual(values = colors, guide = 'none') +
        coord_equal(xlim = c(-r1 * 1.8, d + r2 * 1.8), clip = 'off') +
        labs(title = title) +
        # the current theme without axes or panel, so the title follows it
        theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
              panel.grid = element_blank())

    stats_text <- .OverlapText(ov, stats)
    if(!is.null(stats_text)){
        p <- p +
            labs(caption = stats_text) +
            theme(plot.caption = element_text(size = stats_size * ggplot2::.pt, hjust = 0.5, lineheight = 1.1,
                                              margin = margin(t = 6)))
    }
    attr(p, 'stats') <- ov
    p
}

# region counts c(a_only, shared, b_only) from an ordered or named vector; set names (from '<set>_only') as attribute
.RegionCounts <- function(counts){
    counts <- unlist(counts)
    if(length(counts) != 3 || any(counts < 0)) stop('counts must be three non-negative numbers: a_only, shared, b_only.')
    set_names <- NULL
    if(!is.null(names(counts)) && 'shared' %in% names(counts)){
        only <- names(counts)[names(counts) != 'shared']
        counts <- counts[c(only[1], 'shared', only[2])]
        set_names <- sub('_only$', '', only)
    }
    counts <- as.numeric(counts)
    attr(counts, 'set_names') <- set_names
    counts
}

# area of the intersection of two circles with radii r1, r2 and centre distance d
.LensArea <- function(d, r1, r2){
    if(d >= r1 + r2) return(0)
    if(d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
    a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
    a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
    a3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    a1 + a2 - a3
}

# radii and centre distance; disjoint sets get a small gap, a set inside the other touches its edge
.EulerGeometry <- function(n_a, n_b, shared){
    r1 <- sqrt(n_a / pi)
    r2 <- sqrt(n_b / pi)
    d <- if(shared == 0){
        r1 + r2 + 0.1 * min(r1, r2)
    } else if(shared >= min(n_a, n_b)){
        abs(r1 - r2)
    } else {
        stats::uniroot(function(d) .LensArea(d, r1, r2) - shared,
                       interval = c(abs(r1 - r2) + 1e-9, r1 + r2 - 1e-9))$root
    }
    list(r1 = r1, r2 = r2, d = d)
}

.CircleDf <- function(x0, r, n = 360){
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = x0 + r * cos(theta), y = r * sin(theta))
}

# one line of overlap statistics, e.g. "Jaccard index 0.61   Fisher's exact test: OR 4.2, p = 1.3e-25"
.OverlapText <- function(ov, stats){
    parts <- c(
        if('jaccard' %in% stats) sprintf('Jaccard index %.2f', ov$jaccard),
        if('overlap' %in% stats) sprintf('Overlap coefficient %.2f', ov$overlap_coef),
        if('fisher' %in% stats && !is.na(ov$p_value)) sprintf("Fisher's exact test: OR %s, %s",
                                                               .FormatOddsRatio(ov$odds_ratio), .FormatP(ov$p_value))
    )
    if(length(parts) == 0) return(NULL)
    paste(parts, collapse = '\n')
}

.FormatOddsRatio <- function(x){
    if(is.infinite(x)) return('Inf')
    formatC(x, digits = 2, format = if(x >= 0.01 && x < 1000) 'f' else 'g')
}

.FormatP <- function(p){
    if(p < 1e-300) return('p < 1e-300')
    if(p >= 0.001) return(sprintf('p = %s', formatC(p, digits = 2, format = 'fg')))
    sprintf('p = %s', formatC(p, digits = 1, format = 'e'))
}

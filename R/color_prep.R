#---------------------------------------------------------#
# Shared color preparation for embedding plots
#
# Type detection (numeric = continuous, anything else = categorical) and the conventions of FeatureEmbedding for
# continuous values: 'qNN' percentile or numeric limits, values above the upper limit clipped, cells at or below the
# lower limit flagged to be drawn in light grey underneath, high values drawn last. Categorical values become factors
# (sorted levels unless already a factor) in a shuffled draw order.
#---------------------------------------------------------#

# resolve a 'qNN' percentile string (e.g. 'q99') against the values, or pass a number through
.ResolveLimit <- function(values, limit){
    if(is.character(limit)){
        # a number that became a string in a mixed vector, e.g. c(0, 'q99')
        if(!is.na(suppressWarnings(as.numeric(limit)))) return(as.numeric(limit))
        if(!grepl('^q[0-9.]+$', limit)) stop("limits must be numbers or percentile strings such as 'q99'.")
        return(as.numeric(stats::quantile(values, as.numeric(sub('q', '', limit)) / 100, na.rm = TRUE)))
    }
    limit
}

.IsContinuous <- function(values){
    is.numeric(values)
}

.AsFactor <- function(values){
    if(is.factor(values)) return(values)
    factor(as.character(values), levels = sort(unique(as.character(values))))
}

# prepare the `value` column of a plot data frame for drawing; adds `low` (continuous) and sets the attributes
# `continuous` and `limits`
.PrepareValues <- function(df, plot_min = 'q0', plot_max = 'q100', order_points = TRUE, shuffle = TRUE){
    continuous <- .IsContinuous(df$value)
    if(continuous){
        lo <- .ResolveLimit(df$value, plot_min)
        hi <- .ResolveLimit(df$value, plot_max)
        df$value <- pmin(df$value, hi)
        df$low <- is.na(df$value) | df$value <= lo
        if(isTRUE(order_points)){
            df <- df[order(df$value), ]
        } else if(identical(order_points, 'shuffle')){
            df <- df[sample(nrow(df)), ]
        }
        attr(df, 'limits') <- c(lo, hi)
    } else {
        df$value <- .AsFactor(df$value)
        if(shuffle) df <- df[sample(nrow(df)), ]
    }
    attr(df, 'continuous') <- continuous
    df
}

# row subset that keeps the `continuous` / `limits` attributes
.SubsetRows <- function(df, rows){
    keep <- attributes(df)[intersect(c('continuous', 'limits'), names(attributes(df)))]
    df <- df[rows, , drop = FALSE]
    attributes(df)[names(keep)] <- keep
    df
}

# continuous palette: colfunc(256), reversed by default (inferno: light = low, dark = high)
.ContinuousPalette <- function(colfunc = NULL, rev_colors = TRUE){
    if(is.null(colfunc)) colfunc <- .DefaultColfunc()
    if(is.character(colfunc)) colfunc <- PaletteFunction(colfunc)
    pal <- colfunc(256)
    if(rev_colors) pal <- rev(pal)
    pal
}

.DefaultColfunc <- function(){
    if(requireNamespace('viridis', quietly = TRUE)) return(viridis::inferno)
    function(n) grDevices::hcl.colors(n, 'Inferno')
}

# named color vector from a named vector, a data frame with the groups in the first column and the colors in the
# second (e.g. group/colour, as in PlotEmbedding's color_df), or a palette name (e.g. 'met:Hiroshige', see
# PaletteColors); default: ggplot's hue palette over the levels
.AsNamedColors <- function(colors, levels = NULL){
    if(.IsPaletteName(colors)){
        if(is.null(levels)) stop('a palette name needs the groups to color.')
        return(ColorMap(levels, colors))
    }
    if(is.null(colors)){
        if(is.null(levels)) return(NULL)
        return(stats::setNames(scales::hue_pal()(length(levels)), levels))
    }
    if(is.data.frame(colors)) return(stats::setNames(as.character(colors[[2]]), as.character(colors[[1]])))
    if(is.list(colors)) colors <- unlist(colors)
    colors
}

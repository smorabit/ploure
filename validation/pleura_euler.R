# to launch (PLEURA's R environment, with ploure installed into validation/lib):
# pixi run --manifest-path /scratch_isilon/groups/singlecell/smorabito/pleura/pleura_repo/envs/r/pixi.toml \
#     Rscript validation/pleura_euler.R

#---------------------------------------------------------#
# The PLEURA v1 vs v2 barcode Euler diagram (v2/comparison/compare_v1_v2_figure.Rmd, panel b) drawn with EulerPlot
#
# Draws the panel with the notebook's code and with EulerPlot (statistics off), compares the built layers, and writes
# validation/figures/euler_ploure.pdf (with the Jaccard index below) next to the notebook's panel
# (analysis_v2/comparison_v1/figures/panels/b_barcode_euler.pdf).
#---------------------------------------------------------#

library(tidyverse)
library(scales)

pkg_root <- '/scratch_isilon/groups/singlecell/smorabito/bin/R-plotting-utils'
proj_root <- '/scratch_isilon/groups/singlecell/smorabito/pleura'
library(ploure, lib.loc = file.path(pkg_root, 'validation/lib'))

setwd(proj_root)
fig_dir <- file.path(pkg_root, 'validation/figures/')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
version_colors <- c(v1 = 'grey65', v2 = '#2c7fb8')

overlap <- read.csv('analysis_v2/comparison_v1/barcode_overlap.csv') %>%
    with(setNames(cells, set))

#---------------------------------------------------------#
# the notebook's panel b, verbatim
#---------------------------------------------------------#

lens_area <- function(d, r1, r2){
    if (d >= r1 + r2) return(0)
    if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
    a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
    a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
    a3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    a1 + a2 - a3
}

circle_df <- function(x0, r, n = 360){
    theta <- seq(0, 2 * pi, length.out = n)
    tibble(x = x0 + r * cos(theta), y = r * sin(theta))
}

n_v1 <- overlap[['v1_only']] + overlap[['shared']]
n_v2 <- overlap[['v2_only']] + overlap[['shared']]
r1 <- sqrt(n_v1 / pi)
r2 <- sqrt(n_v2 / pi)
d <- uniroot(function(d) lens_area(d, r1, r2) - overlap[['shared']],
             interval = c(abs(r1 - r2) + 1e-9, r1 + r2 - 1e-9))$root

circles <- bind_rows(
    circle_df(0, r1) %>% mutate(version = 'v1'),
    circle_df(d, r2) %>% mutate(version = 'v2')
)
region_labels <- tibble(
    x = c((-r1 + d - r2) / 2, (d - r2 + r1) / 2, (r1 + d + r2) / 2),
    y = 0,
    label = comma(c(overlap[['v1_only']], overlap[['shared']], overlap[['v2_only']]))
)
set_labels <- tibble(
    x = c(-r1 * 0.75, d + r2 * 0.75), y = c(r1, r2) * 0.85, hjust = c(1, 0),
    label = sprintf('%s (%s)', c('v1', 'v2'), comma(c(n_v1, n_v2)))
)

p_orig <- ggplot() +
    geom_polygon(data = circles, aes(x = x, y = y, fill = version, color = version), alpha = 0.45, linewidth = 0.5) +
    geom_text(data = region_labels, aes(x = x, y = y, label = label), size = 3.5) +
    geom_text(data = set_labels, aes(x = x, y = y, label = label, hjust = hjust), size = 3.8, lineheight = 0.9) +
    scale_fill_manual(values = version_colors, guide = 'none') +
    scale_color_manual(values = version_colors, guide = 'none') +
    coord_equal(xlim = c(-r1 * 1.8, d + r2 * 1.8), clip = 'off') +
    labs(title = 'Final barcodes, v1 vs v2') +
    theme_void() +
    theme(plot.title = element_text(size = 11, hjust = 0.5))

#---------------------------------------------------------#
# ploure, and the layer comparison
#---------------------------------------------------------#

p_new <- EulerPlot(counts = overlap, colors = version_colors, stats = NULL, title = 'Final barcodes, v1 vs v2')

layer_data <- function(p){
    lapply(ggplot_build(p)$data, function(x) x[, intersect(c('x', 'y', 'fill', 'colour', 'label', 'hjust', 'size',
                                                              'alpha', 'linewidth'), colnames(x))])
}
same <- mapply(function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE)), layer_data(p_orig), layer_data(p_new))
cat(sprintf('identical layers: %d / %d\n', sum(same), length(same)))
stopifnot(all(same))

p_stats <- EulerPlot(counts = overlap, colors = version_colors, stats = c('jaccard', 'overlap'),
                     title = 'Final barcodes, v1 vs v2')
print(attr(p_stats, 'stats'))
cairo_pdf(paste0(fig_dir, 'euler_ploure.pdf'), width = 5, height = 4.4)
print(p_stats)
invisible(dev.off())

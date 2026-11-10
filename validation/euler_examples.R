# Simulated Euler examples (no project data): pixi run Rscript validation/euler_examples.R
pkgload::load_all('.', quiet = TRUE)
library(patchwork)
set.seed(1)
universe <- paste0('gene', 1:15000)
sig <- sample(universe, 900)
plots <- list(
    EulerPlot(list('Cohort A' = sig, 'Cohort B' = c(sample(sig, 350), sample(setdiff(universe, sig), 450))),
              universe = universe, title = 'Strong overlap'),
    EulerPlot(list('Cohort A' = sig, 'Cohort B' = sample(universe, 800)), universe = universe,
              colors = 'met:Hiroshige', title = 'Chance overlap'),
    EulerPlot(list('All markers' = sig, 'Top markers' = sig[1:200]), universe = universe,
              stats = c('jaccard', 'overlap', 'fisher'), colors = c('#b2abd2', '#5e3c99'), title = 'Subset')
)
dir.create('validation/figures', showWarnings = FALSE)
cairo_pdf('validation/figures/euler_examples.pdf', width = 14, height = 4.6)
print(wrap_plots(plots, nrow = 1))
invisible(dev.off())

# simulated Milo neighbourhood DA results for the DABeeswarmPlot / DAEmbeddingPlot examples and tests (no real data)
#
# 300 neighbourhoods with index cells drawn from sim_cells: their embedding position, a size, the index cell's
# cluster as the label (a few 'Mixed'), and a fold change shifted per cluster (higher in CD8 T and NK, lower in the
# monocytes and macrophages), with edgeR-like p values and a BH-adjusted stand-in for Milo's SpatialFDR
#
# run from the package root: Rscript data-raw/sim_nhoods.R

load('data/sim_cells.rda')
set.seed(7)
n <- 300
index <- sort(sample(nrow(sim_cells), n))
effects <- c('CD4 T' = 0, 'CD8 T' = 1.5, 'NK' = 1.2, 'B' = -0.3, 'Monocyte' = -1.4, 'Macrophage' = -0.8,
             'DC' = 0, 'Fibroblast' = 0.4, 'Endothelial' = 0)
cluster <- as.character(sim_cells$cluster[index])
stopifnot(all(cluster %in% names(effects)))
se <- 0.6
logfc <- effects[cluster] + stats::rnorm(n, sd = 0.8)
p <- 2 * stats::pnorm(-abs(logfc / se))
label <- cluster
label[sample(n, 12)] <- 'Mixed'

sim_nhoods <- data.frame(
    Nhood = seq_len(n),
    umap_1 = round(sim_cells$umap_1[index], 3),
    umap_2 = round(sim_cells$umap_2[index], 3),
    nhood_size = round(stats::rgamma(n, shape = 6, rate = 0.15)),
    logFC = round(unname(logfc), 4),
    PValue = signif(p, 4),
    SpatialFDR = signif(stats::p.adjust(p, 'BH'), 4),
    cluster = factor(label, levels = c(levels(sim_cells$cluster), 'Mixed'))
)

save(sim_nhoods, file = 'data/sim_nhoods.rda', compress = 'xz')

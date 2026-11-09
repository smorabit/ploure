# simulated single-cell table for the package examples and tests (no real data)
#
# 3,000 cells in 9 clusters of 3 lineages on a 2D "UMAP" of Gaussian blobs, with a separate embedding of the lymphoid
# cells (sub_umap_*), 9 patients in 3 disease groups, QC metrics, a continuous score and 8 marker genes
# (log-normalized-like values with dropout).
#
# run from the package root: Rscript data-raw/sim_cells.R

set.seed(20261109)

clusters <- data.frame(
    cluster = c('CD4 T', 'CD8 T', 'NK', 'B', 'Monocyte', 'Macrophage', 'DC', 'Fibroblast', 'Endothelial'),
    lineage = rep(c('Lymphoid', 'Myeloid', 'Stromal'), c(4, 3, 2)),
    cx = c(-6, -4, -2.5, -7, 4, 6, 5, 1, 2),
    cy = c(2, 4, 1.5, -2, 4, 2.5, 6, -5, -7.5),
    n = c(550, 450, 250, 400, 500, 350, 100, 250, 150)
)
markers <- list(
    'CD4 T' = c('CD3E'), 'CD8 T' = c('CD3E', 'CD8A'), 'NK' = c('NKG7'), 'B' = c('MS4A1'),
    'Monocyte' = c('LYZ', 'CD14'), 'Macrophage' = c('LYZ'), 'DC' = c('LYZ'),
    'Fibroblast' = c('COL1A1'), 'Endothelial' = c('PECAM1')
)
genes <- c('CD3E', 'CD8A', 'NKG7', 'MS4A1', 'LYZ', 'CD14', 'COL1A1', 'PECAM1')

sim_cells <- do.call(rbind, lapply(seq_len(nrow(clusters)), function(i){
    cur <- clusters[i, ]
    data.frame(
        cluster = cur$cluster,
        lineage = cur$lineage,
        umap_1 = rnorm(cur$n, cur$cx, 0.8),
        umap_2 = rnorm(cur$n, cur$cy, 0.8)
    )
}))
n <- nrow(sim_cells)
sim_cells$cell <- sprintf('cell_%04d', seq_len(n))

# patients nested in diseases
patients <- data.frame(patient = sprintf('P%02d', 1:9), disease = rep(c('LC', 'BC', 'OC'), each = 3))
sim_cells$patient <- sample(patients$patient, n, replace = TRUE)
sim_cells$disease <- patients$disease[match(sim_cells$patient, patients$patient)]
sim_cells$sample <- sim_cells$patient

# qc metrics
sim_cells$n_umi <- round(exp(rnorm(n, 8, 0.6)))
sim_cells$n_genes <- round(sim_cells$n_umi^0.85 * runif(n, 0.6, 0.9))
sim_cells$pct_mt <- round(pmin(100, rgamma(n, 2, 0.4)), 2)

# continuous score: high in CD8 T and NK cells
sim_cells$score <- round(rnorm(n, ifelse(sim_cells$cluster %in% c('CD8 T', 'NK'), 0.6, 0.2), 0.12), 3)

# marker genes with dropout
for(g in genes){
    is_marker <- sapply(sim_cells$cluster, function(cl) g %in% markers[[cl]])
    expr <- rnorm(n, ifelse(is_marker, 2.5, 0.4), 0.5)
    detected <- runif(n) < ifelse(is_marker, 0.85, 0.1)
    sim_cells[[g]] <- round(ifelse(detected, pmax(expr, 0.1), 0), 3)
}

# separate embedding of the lymphoid cells
lym <- sim_cells$lineage == 'Lymphoid'
sub_centers <- data.frame(cluster = c('CD4 T', 'CD8 T', 'NK', 'B'), cx = c(-3, 3, 4, -2), cy = c(3, 3, -2, -4))
k <- match(sim_cells$cluster[lym], sub_centers$cluster)
sim_cells$sub_umap_1 <- NA_real_
sim_cells$sub_umap_2 <- NA_real_
sim_cells$sub_umap_1[lym] <- rnorm(sum(lym), sub_centers$cx[k], 1)
sim_cells$sub_umap_2[lym] <- rnorm(sum(lym), sub_centers$cy[k], 1)

sim_cells$cluster <- factor(sim_cells$cluster, levels = clusters$cluster)
sim_cells <- sim_cells[, c('cell', 'umap_1', 'umap_2', 'sub_umap_1', 'sub_umap_2', 'cluster', 'lineage',
                           'patient', 'disease', 'sample', 'n_umi', 'n_genes', 'pct_mt', 'score', genes)]
rownames(sim_cells) <- NULL

save(sim_cells, file = 'data/sim_cells.rda', compress = 'xz')

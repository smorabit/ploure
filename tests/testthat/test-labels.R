test_that('axis labels follow the reduction name', {
    expect_equal(EmbeddingAxisLabels('umap_scvi'), c('UMAP 1', 'UMAP 2'))
    expect_equal(EmbeddingAxisLabels('X_tsne'), c('tSNE 1', 'tSNE 2'))
    expect_equal(EmbeddingAxisLabels('pca'), c('PC 1', 'PC 2'))
    expect_equal(EmbeddingAxisLabels('diffmap'), c('diffmap 1', 'diffmap 2'))
})

test_that('group centers are cells of the group (medoid) or means (centroid)', {
    df <- sim_small()
    med <- .GroupCenters(df$umap_1, df$umap_2, df$cluster, 'medoid')
    expect_equal(as.character(med$group), levels(droplevels(df$cluster)))
    b <- df$cluster == 'B'
    expect_true(any(df$umap_1[b] == med$x[med$group == 'B'] & df$umap_2[b] == med$y[med$group == 'B']))
    cen <- .GroupCenters(df$umap_1, df$umap_2, df$cluster, 'centroid')
    expect_equal(cen$x[cen$group == 'B'], mean(df$umap_1[b]))
})

test_that('label maps replace only the mapped groups', {
    expect_equal(.MapLabels(c('P01', 'P02'), c(P01 = 'P01 LC')), c('P01 LC', 'P02'))
})

test_that('the medoid sampling does not change the global random state', {
    set.seed(3)
    before <- .Random.seed
    .GroupCenters(sim_cells$umap_1, sim_cells$umap_2, sim_cells$cluster)
    expect_identical(.Random.seed, before)
})

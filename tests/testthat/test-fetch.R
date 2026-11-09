test_that('a data frame gives cell, x, y and the requested columns', {
    df <- FetchCellData(sim_cells, reduction = c('umap_1', 'umap_2'), vars = c('cluster', 'score'))
    expect_named(df, c('cell', 'x', 'y', 'cluster', 'score'))
    expect_equal(df$x, sim_cells$umap_1)
    expect_equal(df$cell, sim_cells$cell)
})

test_that('a single reduction name finds <name>_1 / <name>_2 columns in a data frame', {
    df <- FetchCellData(sim_cells, reduction = 'UMAP')
    expect_equal(df$y, sim_cells$umap_2)
    expect_error(FetchCellData(sim_cells, reduction = 'tsne'), 'coordinate columns')
})

test_that('rownames are used when there is no cell column', {
    df <- sim_cells[1:5, c('umap_1', 'umap_2')]
    rownames(df) <- letters[1:5]
    expect_equal(FetchCellData(df, reduction = c('umap_1', 'umap_2'))$cell, letters[1:5])
})

test_that('missing columns and unsupported inputs fail clearly', {
    expect_error(FetchCellData(sim_cells, vars = 'nope'), 'nope')
    expect_error(FetchCellData(matrix(1:4, 2)), 'unsupported')
})

test_that('a Seurat object gives metadata, genes and reduction coordinates', {
    skip_if_not_installed('SeuratObject')
    small <- sim_small()
    obj <- sim_seurat(small)
    df <- FetchCellData(obj, reduction = 'umap', vars = c('cluster', 'CD3E'))
    expect_named(df, c('cell', 'x', 'y', 'cluster', 'CD3E'))
    expect_equal(df$x, small$umap_1)
    expect_equal(df$CD3E, small$CD3E)
    # metadata columns as coordinates
    df2 <- FetchCellData(obj, reduction = c('umap_1', 'umap_2'))
    expect_equal(df2$y, small$umap_2)
    expect_error(FetchCellData(obj, reduction = 'pca'), 'reduction')
    expect_error(FetchCellData(obj, vars = 'NOTAGENE'), 'NOTAGENE')
})

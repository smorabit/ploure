test_that('the outline polygon is closed and holds most highlighted cells', {
    hl <- sim_cells$lineage == 'Lymphoid'
    poly <- OutlinePolygon(sim_cells$umap_1, sim_cells$umap_2, hl)
    expect_named(poly, c('x', 'y', 'group'))
    inside <- rep(FALSE, sum(hl))
    for(g in unique(poly$group)){
        cur <- poly[poly$group == g, ]
        expect_equal(cur[1, c('x', 'y')], cur[nrow(cur), c('x', 'y')], ignore_attr = TRUE)
        inside <- inside | .InPolygon(sim_cells$umap_1[hl], sim_cells$umap_2[hl], cur$x, cur$y)
    }
    expect_gt(mean(inside), 0.95)
    # the stromal cells sit outside
    st <- sim_cells$lineage == 'Stromal'
    expect_false(any(.InPolygon(sim_cells$umap_1[st], sim_cells$umap_2[st], poly$x, poly$y)))
})

test_that('separate islands are kept', {
    hl <- sim_cells$cluster %in% c('B', 'Fibroblast')
    poly <- OutlinePolygon(sim_cells$umap_1, sim_cells$umap_2, hl)
    expect_gte(length(unique(poly$group)), 2)
})

test_that('a rectangle outline is the trimmed bounding box', {
    poly <- OutlinePolygon(sim_cells$umap_1, sim_cells$umap_2, sim_cells$cluster == 'B', shape = 'rectangle')
    expect_equal(nrow(poly), 5)
})

test_that('the point-in-polygon test works on a square', {
    expect_equal(.InPolygon(c(0.5, 2), c(0.5, 0.5), c(0, 1, 1, 0), c(0, 0, 1, 1)), c(TRUE, FALSE))
})

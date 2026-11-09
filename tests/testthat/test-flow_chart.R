test_that('step orders of the branches are merged keeping each order', {
    expect_equal(.MergeOrder(list(c('a', 'b', 'd'), c('a', 'b', 'c', 'd'))), c('a', 'b', 'c', 'd'))
    expect_equal(.MergeOrder(list(c('a', 'd'), c('b', 'c'))), c('a', 'd', 'b', 'c'))
    expect_error(.MergeOrder(list(c('a', 'b'), c('b', 'a'))), 'conflict')
})

test_that('exclusion text lists categories by size, or falls back to the excluded text', {
    st <- .FlowSteps(sim_qc_steps, 'branch')
    rm_df <- .FlowRemovals(sim_qc_removals, 'branch')
    cur <- st[st$branch == 'Myeloid' & st$step == 'Cluster review', ]
    expect_equal(.ExcludedText(cur, rm_df, c('Damaged / low quality' = 'Damaged / LQ')),
                 '−150 cells\nDoublets: 70\nDamaged / LQ: 60\nDissociation stress: 20')
    cur <- st[st$branch == 'all' & st$step == 'Cell QC', ]
    expect_equal(.ExcludedText(cur, rm_df, NULL), '−800 cells\nfailed cell QC')
})

test_that('n_removed defaults to n_in - n_out, and steps are unique within a branch', {
    st <- .FlowSteps(sim_qc_steps, 'branch')
    expect_equal(st$n_removed, sim_qc_steps$n_in - sim_qc_steps$n_out)
    expect_error(.FlowSteps(rbind(sim_qc_steps, sim_qc_steps[1, ]), 'branch'), 'unique')
    expect_error(.FlowSteps(sim_qc_steps[, -1], 'branch'), 'branch')
})

test_that('the chart has one box per step, entry, placeholder and exclusion', {
    p <- QCFlowChart(sim_qc_steps, sim_qc_removals)
    tiles <- ggplot2::ggplot_build(p)$data[[4]]
    n_excl <- sum(sim_qc_steps$n_in > sim_qc_steps$n_out)
    # 13 steps + 3 entry boxes + 2 skipped second reviews + exclusions
    expect_equal(nrow(tiles), 13 + 3 + 2 + n_excl)
    # a trunk-only chart works too
    p <- QCFlowChart(sim_qc_steps[sim_qc_steps$branch == 'all', ])
    expect_s3_class(p, 'ggplot')
})

de <- sim_de[sim_de$lineage == 'Myeloid', ]

test_that('genes are classified by significance, direction and effect size', {
    status <- .DECategory(c(-2, -0.5, 0.5, 2, 3, NA), c(0.01, 0.01, 0.01, 0.01, 0.2, 0.01), 0.05, 1)
    expect_equal(as.character(status), c('lower_big', 'lower_small', 'higher_small', 'higher_big', 'ns', 'ns'))
    expect_equal(as.character(.DECategory(c(-2, 0.5, 1), c(0.01, 0.01, NA), 0.05, NULL)), c('lower', 'higher', 'ns'))
})

test_that('zero p values become the smallest nonzero p value', {
    df <- .DETable(data.frame(gene = c('a', 'b', 'c'), logFC = c(1, 2, 3), padj = c(0, 1e-10, 0.5)), 'logFC', 'padj',
                   'gene', 0.05, 1)
    expect_equal(df$p, c(1e-10, 1e-10, 0.5))
    expect_error(.DETable(de, 'log2FoldChange', 'padj', 'gene', 0.05, 1), 'log2FoldChange')
})

test_that('colors fill in the categories', {
    expect_equal(names(.DEColors(NULL, levels(.DECategory(0, 1, 0.05, 1)))),
                 c('lower_big', 'lower_small', 'higher_small', 'higher_big', 'ns'))
    cols <- .DEColors(c('#111111', '#222222'), c('lower', 'higher', 'ns'))
    expect_equal(unname(cols), c('#111111', '#222222', 'grey85'))
    expect_equal(unname(.DEColors(c(higher_big = 'black'), levels(.DECategory(0, 1, 0.05, 1)))[4]), 'black')
    expect_error(.DEColors(c('red', 'blue', 'green', 'pink'), c('lower', 'higher', 'ns')), 'give 2 colors')
})

test_that('the volcano labels the given, bold and top genes', {
    p <- VolcanoPlot(de, label = 'TECH1', label_bold = 'MARKD', label_top = 3)
    labs <- p$layers[[5]]$data
    top <- utils::head(de$gene[de$adj.P.Val < 0.05][order(de$adj.P.Val[de$adj.P.Val < 0.05])], 3)
    expect_setequal(labs$gene, unique(c('TECH1', 'MARKD', top)))
    expect_equal(labs$face[labs$gene == 'MARKD'], 'bold')
    expect_warning(VolcanoPlot(de, label = 'NOPE'), 'NOPE')
    b <- ggplot2::ggplot_build(VolcanoPlot(de, raster = FALSE))$data[[1]]
    expect_equal(b$y, -log10(de$adj.P.Val))
})

test_that('the count bars match the significant genes per group', {
    p <- DECountPlot(sim_de, group_by = 'lineage')
    rects <- ggplot2::ggplot_build(p)$data[[1]]
    sig <- sim_de[sim_de$adj.P.Val < 0.05, ]
    for(g in levels(sim_de$lineage)){
        y <- nlevels(sim_de$lineage) + 1 - match(g, levels(sim_de$lineage))
        cur <- rects[abs((rects$ymin + rects$ymax) / 2 - y) < 1e-9, ]
        expect_equal(-min(cur$xmin), sum(sig$lineage == g & sig$logFC < 0))
        expect_equal(max(cur$xmax), sum(sig$lineage == g & sig$logFC > 0))
    }
    p <- DECountPlot(sim_de, group_by = 'lineage', fc_cutoff = NULL)
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2 * nlevels(sim_de$lineage))
})

test_that('other DE tools: column names and genes as rownames', {
    # DESeq2-style: log2FoldChange / padj, genes as rownames
    deseq <- data.frame(baseMean = 10, log2FoldChange = de$logFC, pvalue = de$P.Value, padj = de$adj.P.Val,
                        row.names = de$gene)
    p <- VolcanoPlot(deseq, fc_col = 'log2FoldChange', p_col = 'padj', label = 'TECH1')
    expect_equal(p$layers[[5]]$data$gene, 'TECH1')
    ref <- VolcanoPlot(de, label = 'TECH1')
    expect_equal(ggplot2::ggplot_build(p)$data[[1]][, c('x', 'y', 'colour')],
                 ggplot2::ggplot_build(ref)$data[[1]][, c('x', 'y', 'colour')])
    # Seurat FindMarkers-style: avg_log2FC / p_val_adj, genes as rownames, gene_col = NULL
    fm <- data.frame(p_val = de$P.Value, avg_log2FC = de$logFC, p_val_adj = de$adj.P.Val, row.names = de$gene)
    p <- VolcanoPlot(fm, fc_col = 'avg_log2FC', p_col = 'p_val_adj', gene_col = NULL, label_top = 2)
    expect_length(p$layers[[5]]$data$gene, 2)
    # an unadjusted p value column gets a plain P axis title; custom axis titles
    expect_equal(as.character(VolcanoPlot(fm, 'avg_log2FC', 'p_val', NULL)$labels$y), as.character(.PLabel('p_val')))
    expect_equal(VolcanoPlot(fm, 'avg_log2FC', 'p_val_adj', NULL, x_label = 'AUC')$labels$x, 'AUC')
    # edgeR-style counts per group
    edger <- data.frame(genes = de$gene, logFC = sim_de$logFC, FDR = sim_de$adj.P.Val, cluster = sim_de$lineage)
    p <- DECountPlot(edger, group_by = 'cluster', p_col = 'FDR')
    expect_equal(ggplot2::ggplot_build(p)$data[[1]],
                 ggplot2::ggplot_build(DECountPlot(sim_de, group_by = 'lineage'))$data[[1]])
    # helpful errors
    expect_error(VolcanoPlot(deseq), 'Set fc_col')
    plain <- data.frame(log2FoldChange = de$logFC, padj = de$adj.P.Val)
    expect_error(VolcanoPlot(plain, 'log2FoldChange', 'padj', NULL, label = 'TECH1'), 'need gene names')
})

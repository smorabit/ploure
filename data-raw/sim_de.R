# simulated differential expression results for the VolcanoPlot / DECountPlot examples and tests (no real data)
#
# one limma-style table per lineage (3 lineages x 2000 genes): most genes null, some shifted down or up, a few
# strongly; gene names are made up (Gene0001, ...), with a few marker-like names for the label examples
#
# run from the package root: Rscript data-raw/sim_de.R

set.seed(42)
n_genes <- 2000
genes <- sprintf('Gene%04d', seq_len(n_genes))
genes[1:8] <- c('MARKA', 'MARKB', 'MARKC', 'MARKD', 'TECH1', 'TECH2', 'TECH3', 'TECH4')
lineages <- c('Lymphoid' = 0.10, 'Myeloid' = 0.20, 'Stromal' = 0.05)

sim_de <- do.call(rbind, lapply(names(lineages), function(cur_lineage){
    frac <- lineages[[cur_lineage]]
    effect <- rep(0, n_genes)
    shifted <- sample(n_genes, round(frac * n_genes))
    effect[shifted] <- sample(c(-1, 1), length(shifted), replace = TRUE, prob = c(0.6, 0.4)) *
        stats::rgamma(length(shifted), shape = 2, rate = 2.5)
    # the TECH genes shift strongly in every lineage
    effect[5:8] <- c(-2.5, -1.8, 2.2, 1.5)
    ave_expr <- stats::rnorm(n_genes, mean = 4, sd = 2)
    se <- 0.15 + 0.3 / (1 + exp(ave_expr - 2))
    logfc <- effect + stats::rnorm(n_genes, sd = se)
    t <- logfc / se
    p <- 2 * stats::pt(-abs(t), df = 20)
    data.frame(gene = genes, lineage = cur_lineage, logFC = round(logfc, 4), AveExpr = round(ave_expr, 3),
               t = round(t, 3), P.Value = signif(p, 4), adj.P.Val = signif(stats::p.adjust(p, 'BH'), 4))
}))
sim_de$lineage <- factor(sim_de$lineage, levels = names(lineages))

save(sim_de, file = 'data/sim_de.rda', compress = 'xz')
